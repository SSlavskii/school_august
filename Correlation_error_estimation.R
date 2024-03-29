library(bindata)
source("constants_and_functions.R")
real_r = 0.8
sample_size = 200
p_a = 0.1
p_b = 0.1
number_of_simulations = 100 # snp pairs
'
get_r_rmvbin <- function(p_a, p_b, r, sample_size){
  "sample  = mvrnorm(n=sample_size, mu=c(p_a,p_b), 
                    Sigma = matrix(c(1,real_r,real_r,1), nrow = 2)*0.25)"
  sample = rmvbin(sample_size, margprob = c(p_a, p_b), 
                  bincorr = (1-r)*diag(2)+r)
  cor(sample[,1], sample[,2])
}

start = Sys.time()
corr_coeff <- replicate(number_of_simulations, get_r_rmvbin(p_a=p_a, p_b=p_b, r=real_r, sample_size=sample_size))
Sys.time() - start
'


start = Sys.time()
corr_coeff <- replicate(number_of_simulations, get_r_from_haplotypes(get_vector_of_haplotypes(
  number_of_individuals = sample_size, p_a=p_a, p_b=p_b, r=real_r)))
Sys.time() - start


hist(corr_coeff)

var_real = var(corr_coeff)

var_in_r = (1-real_r**2)/(sample_size-2)
var_in_r/var_real

z_s = atanh(corr_coeff)
var_after_z_transform = var(z_s[is.finite(z_s)])
hist(z_s)
var_in_z = 1/(sample_size - 3)
var_in_z/var_after_z_transform

'
var_from_r <- sapply(seq(0,1,0.05), function(real_r) {
  
  corr_coeff = replicate(number_of_simulations, get_r_rmvbin(p_a=p_a, p_b=p_b, r=real_r, sample_size=sample_size))
  corr_coeffs <- replicate(number_of_simulations, get_r(p_a=p_a, p_b=p_b, r=real_r, sample_size=sample_size))
  
  return( c(var(corr_coeff), var(corr_coeffs)) )
})
plot(x=var_from_r[1,], y=var_from_r[2,], xlim = c(0,0.02), ylim = c(0,0.02))
abline(a=0, b=1)
plot(x=seq(0,1,0.05), y=var_from_r[1,])
plot(x=seq(0,1,0.05), y=var_from_r[2,])
'



number_of_simulations = 100
sample_size = 1000
var_from_p <- sapply(seq(0.1,0.9,0.1), function(p) {
  p_a=p
  p_b=p
  real_r = 0.5
  corr_coeffs = replicate(number_of_simulations, 
                          get_r_from_haplotypes(get_vector_of_haplotypes(
                            sample_size, p_a, p_b, r = real_r)))
  var(corr_coeffs)
})
plot(x=seq(0.1,0.9,0.1), y=var_from_p, xlab='p')
abline(a=1/sample_size, b=0)

## bootstraping

real_r = 0.5
sample_size = 100
p_a = 0.1
p_b = 0.1
number_of_simulations = 500
number_of_bootstrap_iterations = 200
sample = get_vector_of_haplotypes(number_of_individuals = sample_size, p_a, p_b, r = real_r)


results <- boot(data=sample, statistic=LD, R=number_of_bootstrap_iterations)
results
boot.ci(results, type="bca")


CIs  = sapply(c(1:number_of_simulations), function(x){
  sample = get_vector_of_haplotypes(number_of_individuals = sample_size, p_a, p_b, r = real_r)
  results <- boot(data=sample, statistic=LD, R=number_of_bootstrap_iterations)
  boot.ci(results, type="bca")$bca[1,c(4,5)]
})
sum(CIs[1,]>real_r | CIs[2,]<real_r)/number_of_simulations*100 # type I error rate



# bootstraping of delta r
number_of_bootstrap_iterations = 200
number_of_simulations = 100
r1 = 0.2
p_a1 = 0.2
p_b1 = 0.2
r2 = 0.8
p_a2 = 0.4
p_b2 = 0.4
sample_size = 100

get_CI_for_delta_r <- function(x){
  sample1 = get_vector_of_haplotypes(number_of_individuals = sample_size, p_a, p_b, r = r1)
  sample2 = get_vector_of_haplotypes(number_of_individuals = sample_size, p_a, p_b, r = r2)
  
  delta_r <- function(data, indices) {
    d <- data[indices,] # allows boot to select sample 
    return(get_r_from_haplotypes(d[,2]) - get_r_from_haplotypes(d[,1]))
  } 
  
  results <- boot(data=cbind(sample1,sample2), statistic=delta_r, R=number_of_bootstrap_iterations)
  results
  boot.ci(results, type="bca")$bca[1,c(4,5)]
  } 

CIs  = sapply(c(1:number_of_simulations), get_CI_for_delta_r)
sum(CIs[1,]>(r2-r1) | CIs[2,]<(r2-r1))/number_of_simulations*100 
