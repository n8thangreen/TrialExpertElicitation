
#' Posterior pc
#'
#' plot the posterior distribution of pc
#' 
#' @param a,b parameters of prior distribution of pC
#' @param mu,sigma2 parameters of prior distribution of theta
#' @param sc,fc,se,fe number of successes and failures observed on CYC and MMF, respectively 
#'
#' @return dataframe containing posterior expectation, SD, mode,
#'   limits of the 90% credibility interval of pC and the normalisation
#'   constant of the posterior joint distribution of (pC, pE)
#' @export
#'
post_pc <- function(a,b, mu, sigma2, sc, fc, se, fe){
  
  ## Use a mesh for integrating pC which should be adequate if P(0.001 <= pC <= 0.999|data) >= 0.998 
  ## and we check this condition is satisfied later.
  gridc = seq(0.001, 0.999, by=0.001)
  midp1 = (0.00001 + 0.001)*0.5
  midp2 = (0.99999 + 0.999)*0.5
  gridc = append(c(0.00001, midp1), gridc)
  gridc = append(gridc, c(midp2, 0.99999))
  lc = length(gridc)
  
  wc  = vector(mode="numeric", length=lc)
  wc[1] = (gridc[3]-gridc[1])/6
  wc[lc] = (gridc[lc] - gridc[lc-2])/6
  for(i in seq(2, (lc-1), by=2)){
    wc[i] = 4*(gridc[i+1] - gridc[i-1])/6
  }
  for(i in seq(3, (lc-2), by=2)){
    wc[i] = (gridc[i+2] - gridc[i-2])/6
  }
  
  ## Need to integrate out pe from the posterior joint distribution
  gride = seq(0.0001, 0.9999, by=0.0001)
  midp1 = (0.00001 + 0.0001)*0.5
  midp2 = (0.99999 + 0.9999)*0.5
  gride = append(c(0.00001, midp1), gride)
  gride = append(gride, c(midp2, 0.99999))
  le = length(gride)
  if(floor(le/2.0) == (le/2.0)){
    stop("Mistake calculating mesh for control arm remission rate in post_pc.R: mesh contains even number of elements when odd number are expected.")
  }
  we  = vector(mode="numeric", length=le)
  we[1] = (gride[3]-gride[1])/6
  we[le] = (gride[le] - gride[le-2])/6
  for(i in seq(2,(le-1), by=2)){
    we[i] = 4*(gride[i+1] - gride[i-1])/6
  }
  for(i in seq(3, (le-2), by=2)){
    we[i] = (gride[i+2] - gride[i-2])/6
  }
  
  dens =  vector(mode="numeric", length=le)
  dens1 = vector(mode="numeric", length=le)
  int =  vector(mode="numeric", length=lc)
  ## integrating out pe from joint posterior distribution of (pe, pc)
  for(i in 1:lc){
    dens =  (gride^(se-1))*((1-gride)^(fe-1))
    dens1 = (-0.5/sigma2)*((log(gride*(1-gridc[i])/(gridc[i]*(1-gride))) - mu)^2)
    dens = 10000*10000*(gridc[i]^(sc + a-1))*((1-gridc[i])^(fc + b-1))*dens*exp(dens1)
    int[i] = sum(we*dens)
  }
  ## calculate the normalisation constant s.t. the joint posterior density of pC and pE integrates to 1
  norm = 1/sum(wc*int)
  int = norm*int
  
  ## Checking to see whether posterior marginal probability that pC is in [0.001, 0.999] is less
  ## than would be the case under a flat Beta(1,1) density.
  lim1 = as.double(0.001)
  lim2 = as.double(0.999)
  
  gridc1 = gridc[which(gridc >= lim1 & gridc <= lim2)]
  int1 = int[which(gridc >= lim1 & gridc <= lim2)]
  lc1 = length(gridc1)
  
  wc1  = vector(mode="numeric", length=lc1)
  wc1[1] = (gridc1[3]-gridc1[1])/6
  wc1[lc1] = (gridc1[lc1] - gridc1[lc1-2])/6
  for(i in seq(2, (lc1-1), by=2)){
    wc1[i] = 4*(gridc1[i+1] - gridc1[i-1])/6
  }
  for(i in seq(3, (lc1-2), by=2)){
    wc1[i] = (gridc1[i+2] - gridc1[i-2])/6
  }
  istop = sum(wc1*int1)
  if(istop < (lim2-lim1)){
    stop("Error calculating posterior distribution of control arm remission rate: Cannot guarantee accuracy of numerical integration routines when integrating U shaped function")
  }else{
    ## calculate the posterior mean and variance of pc
    expect = sum(wc*gridc*int)
    sd1 = sqrt(sum(wc*gridc*gridc*int) - (expect^2))
    modecal = data.frame(gridc, int)
    modelcal_sort = modecal[order(int, decreasing = TRUE), ]
    mode1 = modelcal_sort$gridc[1]
    ## Call a search routine to find the 5th and 95th percentiles of the distribution 
    fval = vector(mode="numeric", length=2)
    fval[1] =  prob_pC(0.005, 0.05, norm, se, sc, fe, fc, a, b, mu, sigma2)
    fval[2] =  prob_pC(0.9, 0.05, norm, se, sc, fe, fc,  a, b, mu, sigma2)
    z = uniroot(prob_pC, interval=c(0.005, 0.9), 0.05, norm, se, sc, fe, fc,  a, b, mu, sigma2, lower = 0.005, upper=0.9, f.lower=fval[1], f.upper=fval[2])
    ci_low = z$root
    
    fval[1] =  prob_pC(0.1, 0.95, norm, se, sc, fe, fc,  a, b, mu, sigma2)
    fval[2] =  prob_pC(0.995, 0.95, norm, se, sc, fe, fc,  a, b, mu, sigma2)
    z = uniroot(prob_pC, interval=c(0.1, 0.995), 0.95, norm,se, sc, fe, fc,  a, b, mu, sigma2, lower = 0.1, upper=0.995, f.lower=fval[1], f.upper=fval[2])
    ci_upp = z$root   
  }
  
  data.frame(expect, sd1, mode1, ci_low, ci_upp, norm)
}


#' plot the posterior distribution of pe
#'
#' @param a,b parameters of prior distribution of pC
#' @param mu,sigma2 parameters of prior distribution of theta
#' @param sc,fc,se,fe number of successes and failures observed on C and E
#'
#' @return dataframe containing posterior expectation, SD, mode,
#'    limits of the 90% credibility interval of pE and the normalisation
#'    constant of the posterior joint distribution of (pC, pE)
#' @export
#'
post_pe <- function(a,b, mu, sigma2, sc, fc, se, fe){
  
  ## Use a mesh for integrating pE which should be adequate if P(0.001 <= pE <= 0.999|data) >= 0.998 
  ## and we check this condition is satisfied later.
  gride = seq(0.001, 0.999, by=0.001)
  midp1 = (0.00001 + 0.001)*0.5
  midp2 = (0.99999 + 0.999)*0.5
  gride = append(c(0.00001, midp1), gride)
  gride = append(gride, c(midp2, 0.99999))
  le = length(gride)
  
  we  = vector(mode="numeric", length=le)
  we[1] = (gride[3]-gride[1])/6
  we[le] = (gride[le] - gride[le-2])/6
  for(i in seq(2, (le-1), by=2)){
    we[i] = 4*(gride[i+1] - gride[i-1])/6
  }
  for(i in seq(3, (le-2), by=2)){
    we[i] = (gride[i+2] - gride[i-2])/6
  }
  
  gridc = seq(0.0001, 0.9999, by=0.0001)
  midp1 = (0.00001 + 0.0001)*0.5
  midp2 = (0.99999 + 0.9999)*0.5
  gridc = append(c(0.00001, midp1), gridc)
  gridc = append(gridc, c(midp2, 0.99999))
  lc = length(gridc)
  
  lc = length(gridc)
  wc  = vector(mode="numeric", length=lc)
  wc[1] = (gridc[3]-gridc[1])/6
  wc[lc] = (gridc[lc] - gridc[lc-2])/6
  for(i in seq(2, (lc-1), by=2)){
    wc[i] = 4*(gridc[i+1] - gridc[i-1])/6
  }
  for(i in seq(3, (lc-2), by=2)){
    wc[i] = (gridc[i+2] - gridc[i-2])/6
  }
  
  ## Need to integrate out pc from the posterior joint distribution of (pE, pC)
  dens =  vector(mode="numeric", length=lc)
  dens1 = vector(mode="numeric", length=lc)
  int =  vector(mode="numeric", length=le)
  for(i in 1:le){
    dens = (gridc^(sc + a -1))*((1-gridc)^(fc + b-1))
    dens1 = (-0.5/sigma2)*((log(gride[i]*(1-gridc)/(gridc*(1-gride[i]))) - mu)^2)
    dens = 10000*10000*(gride[i]^(se-1))*((1-gride[i])^(fe-1))*dens*exp(dens1)
    int[i] = sum(wc*dens)
  } 	
  
  ## calculate the normalisation constant s.t. the joint posterior density of (pC, pE) integrates to 1
  norm = 1/sum(we*int)
  int = norm*int
  
  ## Checking to see whether posterior marginal probability that pE is in [0.001, 0.999] is less
  ## than would be the case under a flat Beta(1,1) density.
  lim1 = as.double(0.001)
  lim2 = as.double(0.999)
  
  gride1 = gride[which(gride >= lim1 & gride <= lim2)]
  int1 = int[which(gride >= lim1 & gride <= lim2)]
  le1 = length(gride1)
  
  we1  = vector(mode="numeric", length=le1)
  we1[1] = (gride1[3]-gride1[1])/6
  we1[le1] = (gride1[le1] - gride1[le1-2])/6
  for(i in seq(2, (le1-1), by=2)){
    we1[i] = 4*(gride1[i+1] - gride1[i-1])/6
  }
  for(i in seq(3, (le1-2), by=2)){
    we1[i] = (gride1[i+2] - gride1[i-2])/6
  }
  istop = sum(we1*int1)
  if(istop < (lim2 - lim1)){
    stop("Error calculating posterior distribution of Experimental remission rate: posterior density p(Experimental) is U (or L) shaped function")
  }else{
    
    ## Run this part of the code to evaluate mean and variance of pE
    expect = sum(we*gride*int)
    sd1 = sqrt(sum(we*gride*gride*int) - (expect^2))
    modecal = data.frame(gride, int)
    modecal_sort = modecal[order(int, decreasing = TRUE), ]
    mode1 = modecal_sort$gride[1]	   
    ## Call a search routine to find the 5th and 95th percentiles of the distribution 
    fval = vector(mode="numeric", length=2)
    fval[1] =  prob_pE(0.005, 0.05, norm,  a, b, mu, sigma2, sc, fc, se, fe, 1)
    fval[2] =  prob_pE(0.9, 0.05, norm,  a, b, mu, sigma2, sc, fc, se, fe, 1)
    z = uniroot(prob_pE, interval=c(0.005, 0.9), 0.05, norm,  a, b, mu, sigma2, sc, fc, se, fe, 1, lower = 0.005, upper=0.9, f.lower=fval[1], f.upper=fval[2])
    ci_low = z$root
    
    fval[1] =  prob_pE(0.1, 0.95, norm,  a, b, mu, sigma2, sc, fc, se, fe, 1)
    fval[2] =  prob_pE(0.995, 0.95, norm,  a, b, mu, sigma2, sc, fc, se, fe, 1)
    z = uniroot(prob_pE, interval=c(0.1, 0.995), 0.95, norm,  a, b, mu, sigma2, sc, fc, se, fe, 1, lower = 0.1, upper=0.995, f.lower=fval[1], f.upper=fval[2])
    ci_upp = z$root 
    
    data.frame(expect, sd1, mode1, ci_low, ci_upp, norm)
  }
}



#' plot the posterior distribution of theta
#'
#' @param a,b parameters of prior distribution of pC
#' @param sc,se,fc,fe number of successes and failures observed on CYC and MMF
#' @param mu,sigma2 parameters of prior distribution of theta
#' @param c2 non-inferiority margin for the trial
#'
#' @return dataframe containing posterior expectation, SD, mode,
#'    limits of the 90% credibility interval of theta and the normalisation
#'    constant of the posterior joint distribution of (pC, theta)
#' @export
#'
post_theta <- function(a,b, sc, se, fc, fe, mu, sigma2, c2){
  
  gridc = seq(0.0001, 0.9999, by=0.0001)
  midp1 = (0.00001 + 0.0001)*0.5
  midp2 = (0.99999 + 0.9999)*0.5
  gridc = append(c(0.00001, midp1), gridc)
  gridc = append(gridc, c(midp2, 0.99999))
  lc = length(gridc)
  
  wc  = vector(mode="numeric", length=lc)
  wc[1] = (gridc[3]-gridc[1])/6
  wc[lc] = (gridc[lc] - gridc[lc-2])/6
  for(i in seq(2, (lc-1), by=2)){
    wc[i] = 4*(gridc[i+1] - gridc[i-1])/6
  }
  for(i in seq(3, (lc-2), by=2)){
    wc[i] = (gridc[i+2] - gridc[i-2])/6
  }
  
  ## Assume posterior distribution of theta is approximately normal
  
  r = as.integer(16)
  mesh = as.integer(6*r -1)
  mesh1 = as.integer(2*mesh-1)
  grid1 = vector(mode="numeric", length= mesh)
  gridt = vector(mode="numeric", length= mesh1)
  pvar = thetavar(a,b, sc, se, fc, fe, mu, sigma2)
  for(i in 1:mesh){
    if(i <= (r-1)){
      grid1[i] = mu + sqrt(pvar)*(-3-4*log(r/i))
    }else if((i >= r) & (i<= 5*r)){
      grid1[i] = mu + sqrt(pvar)*(-3 + 3*(i-r)/(2*r))
    }else{
      grid1[i] = mu + sqrt(pvar)*(3+ 4*log(r/(6*r -i)))
    }
  }   
  ## calculating mesh mid-points 
  for(i in seq(1, mesh, by=1)){
    gridt[2*i-1] = grid1[i]	
  }
  for(i in seq(2, mesh1-1, by=2)){
    gridt[i] = (gridt[i+1] + gridt[i-1])/2
  }    
  ## calculating Simpson's integration weights
  wtheta  = vector(mode="numeric", length=mesh1)
  wtheta[1] = (gridt[3]-gridt[1])/6
  wtheta[mesh1] = (gridt[mesh1] - gridt[mesh1-2])/6
  for(i in seq(2, (mesh1-1), by=2)){
    wtheta[i] = 4*(gridt[i+1] - gridt[i-1])/6
  }
  for(i in seq(3, (mesh1-2), by=2)){
    wtheta[i] = (gridt[i+2] - gridt[i-2])/6
  }   
  dens = vector(mode = "numeric", length=lc)
  dens1 = vector(mode = "numeric", length=lc)
  int =  vector(mode = "numeric", length=mesh1)
  ## integrating over pc for the joint density of (pc, theta)
  for(i in 1:mesh1){
    dens1 = (gridc^(sc-fe+a))*((1-gridc)^(fc+b-2+fe))
    v = (1 + exp(-gridt[i])*(1-gridc)/gridc)*((1 + exp(gridt[i])*gridc/(1-gridc))^(1/(se+fe-1)))
    v = v^(se+fe-1)
    dens = 10000*10000*dens1/v
    int[i] = sum(wc*dens)*exp((1-fe)*gridt[i])*exp((-0.5/sigma2)*((gridt[i] - mu)^2)) 		
  }
  ## calculating normalising constant of the joint posterior density of (pc, theta)
  norm = 1/sum(wtheta*int)
  int = norm*int
  ## calculate the posterior mean and variance of theta
  
  expect = sum(wtheta*gridt*int)
  sd1 = sqrt(sum(wtheta*gridt*gridt*int) - (expect^2))
  modecal = data.frame(gridt, int)
  modecal_sort = modecal[order(int, decreasing = TRUE), ]
  mode1  = modecal_sort$gridt[1]
  
  ## Call a search routine to find the 5th and 95th percentiles of the distribution 
  fval = vector(mode="numeric", length=2)
  fval[1] =  prob_theta(-6.0, 0.05, norm, se, sc, fe, fc,  a, b, mu, sigma2, expect, sd1)
  fval[2] =  prob_theta(6.0, 0.05, norm, se, sc, fe, fc,  a, b, mu, sigma2, expect, sd1)
  z = uniroot(prob_theta, interval=c(-6.0, 6.0), 0.05, norm, se, sc, fe, fc,  a, b, mu, sigma2, expect, sd1, lower = -6.0, upper=6.0, f.lower=fval[1], f.upper=fval[2])
  ci_low = z$root
  
  fval[1] =  prob_theta(-6.0, 0.95, norm, se, sc, fe, fc,  a, b, mu, sigma2, expect, sd1)
  fval[2] =  prob_theta(6.0, 0.95, norm, se, sc, fe, fc,  a, b, mu, sigma2, expect, sd1)
  z = uniroot(prob_theta, interval=c(-6.0, 6.0), 0.95, norm, se, sc, fe, fc,  a, b, mu, sigma2, expect, sd1, lower = -6.0, upper=6.0, f.lower=fval[1], f.upper=fval[2])
  ci_upp = z$root
  
  data.frame(expect, sd1, mode1, ci_low, ci_upp, norm)
}

#' produce a course estimate of posterior variance of theta
#'
#' @param a,b parameters of prior distribution of pC
#' @param sc,se,fc,fe number of successes and failures observed on CYC and MMF
#' @param mu,sigma2 parameters of prior distribution of theta
#'
#' @return estimate of variance of the posterior distribution of theta.
#' @export
#'
thetavar <- function(a,b, sc, se, fc, fe, mu, sigma2){
  
  gridc = seq(0.001, 0.999, by=0.001)
  midp1 = (0.00001 + 0.001)*0.5
  midp2 = (0.99999 + 0.999)*0.5
  gridc = append(c(0.00001, midp1), gridc)
  gridc = append(gridc, c(midp2, 0.99999))
  lc = length(gridc)
  if(floor(lc/2.0) == (lc/2.0)){
    stop("Error in thetavar.R: integration mesh contains even number of elements when odd number are expected")
  }
  wc  = vector(mode="numeric", length=lc)
  wc[1] = (gridc[3]-gridc[1])/6
  wc[lc] = (gridc[lc] - gridc[lc-2])/6
  for(i in seq(2, (lc-1), by=2)){
    wc[i] = 4*(gridc[i+1] - gridc[i-1])/6
  }
  for(i in seq(3, (lc-2), by=2)){
    wc[i] = (gridc[i+2] - gridc[i-2])/6
  }
  
  r = as.integer(16)
  mesh = as.integer(6*r -1)
  mesh1 = as.integer(2*mesh-1)
  grid1 = vector(mode="numeric", length= mesh)
  gridt = vector(mode="numeric", length= mesh1)
  for(i in 1:mesh){
    if(i <= (r-1)){
      grid1[i] = mu + sqrt(sigma2)*(-3-4*log(r/i))
    }else if((i >= r) & (i<= 5*r)){
      grid1[i] = mu + sqrt(sigma2)*(-3 + 3*(i-r)/(2*r))
    }else{
      grid1[i] = mu + sqrt(sigma2)*(3+ 4*log(r/(6*r -i)))
    }
  }   
  ## calculating mesh mid-points 
  for(i in seq(1, mesh, by=1)){
    gridt[2*i-1] = grid1[i]	
  }
  for(i in seq(2, mesh1-1, by=2)){
    gridt[i] = (gridt[i+1] + gridt[i-1])/2
  }    
  ## calculating Simpson's integration weights
  wtheta  = vector(mode="numeric", length=mesh1)
  wtheta[1] = (gridt[3]-gridt[1])/6
  wtheta[mesh1] = (gridt[mesh1] - gridt[mesh1-2])/6
  for(i in seq(2, (mesh1-1), by=2)){
    wtheta[i] = 4*(gridt[i+1] - gridt[i-1])/6
  }
  for(i in seq(3, (mesh1-2), by=2)){
    wtheta[i] = (gridt[i+2] - gridt[i-2])/6
  }   
  dens = vector(mode = "numeric", length=lc)
  dens1 = vector(mode = "numeric", length=lc)
  int =  vector(mode = "numeric", length=mesh1)
  ## integrating over pc for the joint density of (pc, theta)
  ## We only want a course estimate of posterior variance of theta so don't worry about inaccuracies of numerical integration
  for(i in 1:mesh1){
    dens1 = (gridc^(sc-fe+a))*((1-gridc)^(fc+b-2+fe))
    v = (1 + exp(-gridt[i])*(1-gridc)/gridc)*((1 + exp(gridt[i])*gridc/(1-gridc))^(1/(se+fe-1)))
    v = v^(se+fe-1)
    dens = 10000*10000*dens1/v
    int[i] = sum(wc*dens)*exp((1-fe)*gridt[i])*exp((-0.5/sigma2)*((gridt[i] - mu)^2)) 			
  }
  ## calculating normalising constant of the joint posterior density of (pc, theta)
  norm = 1/sum(wtheta*int)
  int = norm*int
  ## calculate the posterior mean and variance of theta
  expect = sum(wtheta*gridt*int)
  
  sum(wtheta*gridt*gridt*int) - (expect^2)
}

