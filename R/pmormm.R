# originally called 'getnature' in affyR

pmormm <- function(cdf) {
  a.i <- which(cdf@pbase.levels == "A")
  t.i <- which(cdf@pbase.levels == "T")
  g.i <- which(cdf@pbase.levels == "G")
  c.i <- which(cdf@pbase.levels == "C")
  md <- dim(cdf@name)

  # init to 'NA'
  nature <- matrix(NA,md[1],md[2])
  i <- which(
             ((cdf@pbase == a.i) & (cdf@tbase == t.i)) | 
             ((cdf@pbase == t.i) & (cdf@tbase == a.i)) |
             ((cdf@pbase == g.i) & (cdf@tbase == c.i)) |
             ((cdf@pbase == c.i) & (cdf@tbase == g.i))
             )
  # set the PM to TRUE
  nature[i] <- TRUE
  i <- which((! is.na(cdf@atom)) & is.na(nature))
  # set the MM to FALSE
  nature[i] <- FALSE
  return(nature)   
}


















