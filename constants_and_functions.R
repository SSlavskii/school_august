haplotype_constants = c(11,12,21,22)
names(haplotype_constants) <- c("AB", 'Ab', 'aB', 'ab')

get_haplotypes_freq <- function(p_a, p_b, r){
  p_ab = r * sqrt(p_a * (1 - p_a) * p_b * (1 - p_b)) + p_a*p_b
  p_aB = p_a - p_ab
  p_Ab = p_b - p_ab
  
  if ((p_Ab<0) | (p_aB<0)) {
    print("Wrong allele frequences for given r")
    return -1
  }
  p_AB = 1 - p_ab - p_aB - p_Ab
  return(c(p_AB, p_Ab, p_aB, p_ab))
} 


get_vector_of_haplotypes <- function(number_of_individuals, p_a, p_b, r){
  sample(haplotype_constants, size=number_of_individuals*2, replace=TRUE, 
         prob=get_haplotypes_freq(p_a = p_a, p_b = p_b, r = r))
}

haplotypes_to_genotypes <- function(vector_of_haplotypes){
  N = length(vector_of_haplotypes)/2
  genotypes_a <- vector_of_haplotypes[1:N] %/% 10 + vector_of_haplotypes[(N+1):(2*N)] %/% 10 - 2
  genotypes_b <- vector_of_haplotypes[1:N] %% 10 + vector_of_haplotypes[(N+1):(2*N)] %% 10 - 2
  return(matrix(c(genotypes_a, genotypes_b),nrow = 2, byrow = T))
}

get_r_from_haplotypes <- function(vector_of_haplotypes){
  actual_hapl_freq = sapply(haplotype_constants, function(x) 
    { sum(vector_of_haplotypes==x)/length(vector_of_haplotypes) })
  P = actual_hapl_freq["AB"]
  Q = actual_hapl_freq["Ab"]
  R = actual_hapl_freq["aB"]
  S = actual_hapl_freq["ab"]
  r = (P*S - Q*R) / sqrt((P+Q)*(R+S)*(P+R)*(Q+S))
  return(unname(r))
}


get_r_from_genotypes <- function(vector_of_haplotypes){
  matrix_of_genotypes <- haplotypes_to_genotypes(vector_of_haplotypes)
  return(cor(matrix_of_genotypes[1,], matrix_of_genotypes[2,]))
}


library(boot)

LD <- function(data, indices) {
  d <- data[indices] # allows boot to select sample
  return(get_r_from_haplotypes(d))
}