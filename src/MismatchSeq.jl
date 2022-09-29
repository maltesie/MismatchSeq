module MismatchSeq

using BioSequences, BioAlignments, GenomicFeatures, BioGenerics, FASTX.FASTA, XAM.BAM

export BaseCoverage, mismatchfractions, mismatchcontexthist

include("sequence.jl")
include("basecoverage.jl")

end # module MismatchSeq
