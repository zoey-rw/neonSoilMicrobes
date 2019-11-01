# helper functions


# from the dada2 tutorial
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}


# check for primers
checkPrimers <- function(fwd, rev, FWD.orients, REV.orients){
  print(
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fwd), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = rev), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fwd), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = rev))
)
}