#struct BaseAnnotation <: AnnotationStyle
#    type::String
#    name::String
#    ref::Vector{Int}
#    a::Vector{Int}
#    t::Vector{Int}
#    g::Vector{Int}
#    c::Vector{Int}
#    gap::Vector{Int}
#    ins::Vector{Int}
#end
#
#function BaseAnnotation(feature::Interval{Annotation}, base_coverage::BaseCoverage)
#    left = leftposition(feature)
#    right = rightposition(feature)
#    seq = base_coverage.ref_seq[left:right]
#    ispositivestrand(feature) || reverse_complement!(seq)
#    count = ispositivestrand(feature) ? errorcov.fcount : errorcov.rcount
#    r = ispositivestrand(feature) ? (left:right) : (right:-1:left)
#    ref = Int[(seq[i] in (DNA_A, DNA_T, DNA_G, DNA_C)) ? count[seq[i]][ii] : 0 for (i, ii) in enumerate(r)]
#    BaseAnnotation(type(feature), name(feature), ref, count[ref][:A][r], count[:T][r], count[:G][r], count[:C][r], count[:Gap][r], count[:Ins][r])
#end
#
#function coverage(feature::Interval{BaseAnnotation})
#   return feature.metadata.a .+ feature.metadata.t .+ feature.metadata.g .+ feature.metadata.c .+ feature.metadata.gap
#end
#
#refcount(feature::Interval{BaseAnnotation}) = feature.metadata.ref

struct BaseCoverage
    genome::Genome
    fcount::Dict{String, Dict{Symbol, Vector{Int}}}
    rcount::Dict{String, Dict{Symbol, Vector{Int}}}
end

function BaseCoverage(bam_file::String, genome::Genome; include_secondary_alignments=true, is_reverse_complement=false,
    only_positive_strand=false, quality_cut = 0x01)
    fcount = Dict(chr=>zeros(Int, 7, length(genome.chroms[chr])) for chr in keys(genome.chroms))
    rcount = Dict(chr=>zeros(Int, 7, length(genome.chroms[chr])) for chr in keys(genome.chroms))
    record = BAM.Record()
    reader = BAM.Reader(open(bam_file))
    seq = LongDNA{4}(0)
    while !eof(reader)
        read!(reader, record)
        BAM.ismapped(record) || continue
        hasxatag(record) && continue
        !isprimary(record) && !include_secondary_alignments && continue
        is_reverse = is_reverse_complement #(isread2(record) != is_reverse_complement)
        ref_name::String = BAM.refname(record)
        seq = BAM.sequence(record)
        qual = BAM.quality(record)
        length(seq) > 0 || continue
        is_positive = (BAM.ispositivestrand(record) != is_reverse) || only_positive_strand
        is_reverse && reverse_complement!(seq)
        is_positive || complement!(seq)
        count = is_positive ? fcount : rcount
        offset, nops = BAM.cigar_position(record)
        current_ref::Int = BAM.leftposition(record)
        current_seq::Int = 0
        for data_index in offset:4:offset + (nops - 1) * 4
            x = unsafe_load(Ptr{UInt32}(pointer(record.data, data_index)))
            op = BioAlignments.Operation(x & 0x0F)
            n = x >> 4
            r::UnitRange{Int} = current_ref:current_ref+n-1
            if op === OP_INSERT
                current_seq += n
                count[ref_name][7, current_ref] += 1
            elseif op === OP_SOFT_CLIP
                current_seq += n
            elseif BioAlignments.isdeleteop(op)
                for ref_pos in r
                    count[ref_name][5, ref_pos] += 1
                end
                current_ref += n
            elseif BioAlignments.ismatchop(op)
                for (ii, ref_pos) in enumerate(r)
                    if qual[current_seq + ii] >= quality_cut
                        base = seq[current_seq + ii]
                        base === DNA_A ? count[ref_name][1, ref_pos] += 1 :
                        base === DNA_T ? count[ref_name][2, ref_pos] += 1 :
                        base === DNA_G ? count[ref_name][3, ref_pos] += 1 :
                        base === DNA_C ? count[ref_name][4, ref_pos] += 1 :
                        count[ref_name][6, ref_pos] += 1
                    else
                        count[ref_name][6, ref_pos] += 1
                    end
                end
                current_ref += n
                current_seq += n
            end
        end
    end
    fdict = Dict(chr=>Dict(dna=>fcount[chr][i, :] for (i, dna) in enumerate((:A, :T, :G, :C, :Gap, :N, :Ins))) for chr in keys(genome.chroms))
    rdict = Dict(chr=>Dict(dna=>rcount[chr][i, :] for (i, dna) in enumerate((:A, :T, :G, :C, :Gap, :N, :Ins))) for chr in keys(genome.chroms))
    return BaseCoverage(genome, fdict, rdict)
end

function mismatchfractions(base_coverage::BaseCoverage, from::Symbol, to::Symbol; direction=:both)
    from in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for from::Symbol not supported!"))
    to in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for to::Symbol not supported!"))
    d = Dict(:A=>DNA_A, :T=>DNA_T, :G=>DNA_G, :C=>DNA_C, :Gap=>DNA_Gap, :N=>DNA_N)
    base_from = d[from]
    direction in (:forward, :reverse, :both) || throw(AssertionError("direction has to be :forward, :reverse or :both"))
    fstats = Dict(chr=>Float64[] for chr in keys(base_coverage.genome.chroms))
    rstats = Dict(chr=>Float64[] for chr in keys(base_coverage.genome.chroms))
    for chr in keys(base_coverage.genome.chroms)
        findex = base_coverage.genome[chr] .=== base_from
        rindex = BioSequences.complement(base_coverage.genome[chr]) .=== base_from
        f = base_coverage.fcount[chr][to][findex] ./ sum(base_coverage.fcount[chr][s] for s in (:A, :T, :G, :C, :Gap, :N))[findex]
        r = base_coverage.rcount[chr][to][rindex] ./ sum(base_coverage.rcount[chr][s] for s in (:A, :T, :G, :C, :Gap, :N))[rindex]
        append!(fstats[chr], filter!(x->(!isnan(x) && !iszero(x)), f))
        append!(rstats[chr], filter!(x->(!isnan(x) && !iszero(x)), r))
    end
    return fstats, rstats
end

function mismatchpositions(base_coverage::BaseCoverage, from::Symbol, to::Symbol; ratio_cut=0.9)
    from in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for from::Symbol not supported!"))
    to in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for to::Symbol not supported!"))
    d = Dict(:A=>DNA_A, :T=>DNA_T, :G=>DNA_G, :C=>DNA_C, :Gap=>DNA_Gap, :N=>DNA_N)
    base_from = d[from]
    pos = Interval{Tuple{Int, Float64}}[]
    for chr in keys(base_coverage.genome.chroms)
        index = base_coverage.genome[chr] .=== base_from
        ratios = base_coverage.fcount[chr][to] ./ sum(base_coverage.fcount[chr][s] for s in (:A, :T, :G, :C, :Gap, :N))
        ratio_index = (ratios .> ratio_cut) .& index
        for (r, p) in zip(ratios[ratio_index], findall(ratio_index .=== true))
            mc = base_coverage.fcount[chr][to][p]
            push!(pos, Interval(chr, p:p, STRAND_POS, (mc, round(r; digits=3))))
        end
        index = BioSequences.complement(base_coverage.genome[chr]) .=== base_from
        ratios = base_coverage.rcount[chr][to] ./ sum(base_coverage.rcount[chr][s] for s in (:A, :T, :G, :C, :Gap, :N))
        ratio_index = (ratios .> ratio_cut) .& index
        for (r, p) in zip(ratios[ratio_index], findall(ratio_index .=== true))
            mc = base_coverage.rcount[chr][to][p]
            push!(pos, Interval(chr, p:p, STRAND_NEG, (mc, round(r; digits=3))))
        end
    end
    return pos
end

function mismatchpositions(base_coverage::BaseCoverage; check=(:A, :T, :G, :C), ratio_cut=0.9)
    pos = Interval{Annotation}[]
    for from in check, to in check
        from === to && continue
        mpos = mismatchpositions(base_coverage, from, to; ratio_cut=ratio_cut)
        append!(pos, [Interval(refname(i), leftposition(i), leftposition(i), strand(i),
        Annotation("SNP", "$from->$to", Dict("mutation"=>"$from->$to",
                "count"=>"$(i.metadata[1])",
                "frequency"=>"$(i.metadata[2])"))) for i in mpos])
    end
    return Features(pos)
end

function deletionpositions(base_coverage::BaseCoverage; check=(:A, :T, :G, :C), ratio_cut=0.9)
    pos = Interval{Annotation}[]
    for from in check
        mpos = mismatchpositions(base_coverage, from, :Gap; ratio_cut=ratio_cut)
        append!(pos, [Interval(refname(i), leftposition(i), leftposition(i), strand(i),
        Annotation("DEL", "$from", Dict("base"=>"$from",
            "count"=>"$(i.metadata[1])",
            "frequency"=>"$(i.metadata[2])"))) for i in mpos])
    end
    merged_pos = Interval{Annotation}[]
    merger = Interval{Annotation}[]
    for check in (true, false)
        for pos_interval in pos
            ispositivestrand(pos_interval) == check || continue
            show(pos_interval)
            if !isempty(merger) && (leftposition(pos_interval) != (rightposition(merger[end])+1))
                push!(merged_pos, length(merger) == 1 ? merger[1] : Interval(refname(merger[end]),
                                            leftposition(merger[1]),
                                            rightposition(merger[end]),
                                            strand(merger[end]),
                                            Annotation("DEL", join(name(m) for m in merger),
                                                        Dict("base"=>join(name(m) for m in merger),
                                                            "count"=>"$(mean(param(m, "count", Int) for m in merger))",
                                                            "frequency"=>"$(mean(param(m, "frequency", Float64) for m in merger))"))))
                empty!(merger)
            end
            push!(merger, pos_interval)
        end
    end
    return Features(merged_pos)
end

all_perm(xs, n) = vec(map(collect, Iterators.product(ntuple(_ -> xs, n)...)))
function mismatchcontexthist(base_coverage::BaseCoverage, from::Symbol, to::Symbol; pm=2, bins=20, ratio_cut=0.0)
    from in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for from::Symbol not supported!"))
    to in (:A, :T, :G, :C, :Gap, :N, :Ins) || raise(AssertionError("Value for to::Symbol not supported!"))
    kmer_histograms = Dict(LongDNA{4}(c)=>zeros(Int, bins) for c in all_perm([DNA_A, DNA_T, DNA_G, DNA_C], 2*pm+1))
    mpos = mismatchpositions(base_coverage, from, to; ratio_cut)
    for mismatch_interval in mpos
        p = leftposition(mismatch_interval)
        freq_index = Int(floor(mismatch_interval.metadata[2] * (bins-1))) + 1
        ref = refname(mismatch_interval)
        kmer = base_coverage.genome[ref][p-pm:p+pm]
        kmer_histograms[kmer][freq_index] += 1
    end
    return kmer_histograms
end