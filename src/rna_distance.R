library(msa)
library(DECIPHER)

rna <- RNAStringSet(c("ACUGG", "AACUG","ACUAG" ))
 
# msa_seqinr <- msaConvert(msa(dna, order="input"), type = "seqinr::alignment")
# seqinr::dist.alignment(msa_seqinr)
# dist.hamming(msa_seqinr)

DistanceMatrix(RNAStringSet(msa(rna, order="input"))) 
