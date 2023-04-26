
#' Effective sample size (ESS) theta
#'
#' Compute the effective sample size of the prior distribution of theta.
#'
#' @param R number of subjects randomised to MMF for every one on control
#'   (E:C allocation ratio R:1) in a hypothetical study
#' @param mu,sigma2 mean and variance of the prior distribution of \deqn{\theta~N(mu, sigma2)}
#' @param a,b parameters of the prior distribution of \deqn{pC~Beta(a,b)}
#'
#' @return  prior effective sample size of theta prior (expressed in terms of the effective
#'   total number of observations accrued across MMF and CYC)
#' @export
#'
ess_theta <- function(R, mu, sigma2, a, b){
  
  gride = seq(0.001, 0.999, by=0.001)
  midp1 = (0.00001 + 0.001)*0.5
  midp2 = (0.99999 + 0.999)*0.5
  gride = append(c(0.00001, midp1), gride)
  gride = append(gride, c(midp2, 0.99999)) 
  
  le = length(gride)
  if(floor(le/2.0) == (le/2.0)){
    stop("Stop: mistake calculating mesh for prior ESSs: mesh contains even number of elements when odd number are expected.")
  }
  we  = vector(mode="numeric", length=le)
  we[1] = (gride[3]-gride[1])/6.0
  we[le] = (gride[le] - gride[le-2])/6.0
  
  for(i in seq(2, (le-1), by=2)){
    we[i] = 4*(gride[i+1] - gride[i-1])/6.0
  }
  for(i in seq(3, (le-2), by=2)){
    we[i] = (gride[i+2] - gride[i-2])/6.0
  }
  
  ## Set up a mesh for integrating over pC
  gridc = seq(0.001, 0.999, by=0.001)
  midp1 = (0.00001 + 0.001)*0.5
  midp2 = (0.99999 + 0.999)*0.5
  gridc = append(c(0.00001, midp1), gridc)
  gridc = append(gridc, c(midp2, 0.99999)) 
  lc = length(gridc)
  if(floor(lc/2.0) == (lc/2.0)){
    stop("Stop: mistake calculating mesh for prior ESSs: mesh contains even number of elements when odd number are expected.")
  }
  wc  = vector(mode="numeric", length=lc)
  wc[1] = (gridc[3]-gridc[1])/6.0
  wc[lc] = (gridc[lc] - gridc[lc-2])/6.0
  for(i in seq(2, (lc-1), by=2)){
    wc[i] = 4*(gridc[i+1] - gridc[i-1])/6.0
  }
  for(i in seq(3, (lc-2), by=2)){
    wc[i] = (gridc[i+2] - gridc[i-2])/6.0
  }
  
  int = vector(mode="numeric", length = lc)
  dens =  vector(mode="numeric", length=le)
  dens1 = vector(mode="numeric", length=le)
  pbar = vector(mode="numeric", length=le)
  check = vector(mode="numeric", length =lc)
  
  ## Integrating over the joint prior density for (pE, pC)
  for(i in 1:lc){
    pbar = (R*gride + gridc[i])/(R+1)
    dens = (gridc[i]^(a-1))*((1-gridc[i])^(b-1))/(gride*(1-gride))
    dens1 = (-0.5/sigma2)*((log(gride*(1-gridc[i])/(gridc[i]*(1-gride))) - mu)^2)
    check[i] = sum(we*dens*exp(dens1))
    dens = pbar*(1-pbar)*dens*exp(dens1)
    int[i] = sum(we*dens)
  }
  check_dens = sum(wc*check)/(beta(a,b)*sqrt(2*pi*sigma2))
  if(check_dens < 0.99 | check_dens > 1.01){
    stop("Stop evaluating prior effective sample size of log-odds ratio: Insufficient precision attained integrating over joint prior density for (pC, pE) in ess_theta.c")
  }else{
    expect = sum(wc*int)/(beta(a,b)*sqrt(2*pi*sigma2))
    neff = ((R+1)^2)/(sigma2*expect*R)	
    return(neff)
  }
}


#' Effective sample size (ESS) pc
#' 
#' Calculate prior effective sample size of \deqn{log[pc/(1-pC)]}.
#'
#' @param a,b parameters of the prior distribution of \deqn{pC~Beta(a,b)}
#'
#' @return Effective sample size of prior for \deqn{\log[pC/(1-pc)]}
#'   (expressed in terms of effective observations on CYC)
#' @export
#'
ess_pc <- function(a, b){
  ## Calculating prior variance of omega = log[pC/(1-pc)], where pC~B(a,b)
  ## Create a grid for omega assuming it is approximately normally distributed.
  vart = logoddspc(a,b)
  r = as.integer(32)
  mesh = as.integer(6*r -1)
  mesh1 = as.integer(2*mesh-1)
  grid1 = vector(mode="numeric", length= mesh)
  gridt = vector(mode="numeric", length= mesh1)
  ## centre mesh at log[E(pc)/(1-E(pc))], where E(pc) is prior mean of pC~B(a,b)
  mu = a/(a+b)
  mu = log(mu/(1-mu))  	
  for(i in 1:mesh){
    if(i <= (r-1)){
      grid1[i] = mu + sqrt(vart)*(-3-4*log(r/i))
    }else if((i >= r) & (i<= 5*r)){
      grid1[i] = mu + sqrt(vart)*(-3 + 3*(i-r)/(2*r))
    }else{
      grid1[i] = mu + sqrt(vart)*(3+ 4*log(r/(6*r -i)))
    }
  }   
  ## calculating mesh mid-points 
  for(i in seq(1, mesh, by=1)){
    gridt[2*i-1] = grid1[i]	
  }
  for(i in seq(2, mesh1-1, by=2)){
    gridt[i] = (gridt[i+1] + gridt[i-1])/2.0
  }    
  ## calculating Simpson's integration weights
  wtheta  = vector(mode="numeric", length=mesh1)
  wtheta[1] = (gridt[3]-gridt[1])/6.0
  wtheta[mesh1] = (gridt[mesh1] - gridt[mesh1-2])/6.0
  for(i in seq(2, (mesh1-1), by=2)){
    wtheta[i] = 4*(gridt[i+1] - gridt[i-1])/6.0
  }
  for(i in seq(3, (mesh1-2), by=2)){
    wtheta[i] = (gridt[i+2] - gridt[i-2])/6.0
  }
  ## Integrating Fisher's expected information for log[pC/(1-pC)] (divided by neff) across
  ## prior density of log[pC/(1-pC)] (which we assume is approximately normally distributed)
  dens = vector(mode = "numeric", length = mesh1)
  dens1 = vector(mode = "numeric", length = mesh1)
  dens = exp((a+1)*gridt)
  dens1 = (1 + exp(gridt))^(a + b + 2)
  dens = dens/(dens1*beta(a,b))
  
  expect = sum(wtheta*dens)
  neff = 1/(vart*expect)
  
  neff
}
