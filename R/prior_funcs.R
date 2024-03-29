
#' Find the parameters of prior Beta distribution of pc~Beta(a,b)
#'
#' @param mode prior mode of pc as identified by Day 1 elicitation Q1
#' @param percentile25 25th percentile of prior distribution for pC as identified by Day 1 elicitation Q2.
#'
#' @return vector (a,b) (parameters of Beta prior distribution for pC).
#' @export
#'
prior_beta <- function(mode, percentile25) {
  
  no_root_in_range <- !check_interval_valid_for_a(mode, percentile25)
  
  if (no_root_in_range) {
    shinyalert::shinyalert("This combination of Q1 and Q2 does not permit the construction of a prior distribution", "Try different values", type = "error")
    stop("\nGiven answers to elicitation questions Q1 and Q2, ",
         "we cannot determine a Beta prior distribution for control arm/steroid remission rate.\n",
         "Please revise either the answer to elicitation Q1 or Q2.\n",
         "Error in answers to elication questions Q1 and Q2: ",
         "cannot determine Beta prior distribution for control arm/steroid remission rate.")
  }
  
  centred_beta_percentile <-
    function(a, mode) beta_percentile25(a, mode) - percentile25
  
  z <- uniroot(f = centred_beta_percentile, interval = c(0.5, 50), mode)
  
  a_root <- z$root
  b_root <- (a_root - 1)/mode - (a_root - 2) 
  
  check_accuracy_of_uniroot(a_root, b_root) 
  
  c(a = a_root,
    b = b_root)
}

#
check_accuracy_of_uniroot <- function(a, b) {
  # quantile limits
  q_low <- 0.001
  q_high <- 0.999
  
  diff_pbeta <-
    pbeta(q_high, a, b, lower.tail = TRUE) - pbeta(q_low, a, b, lower.tail = TRUE)
  
  if (diff_pbeta < (q_high - q_low)){
    shinyalert::shinyalert("This combination of Q1 and Q2 does not permit the construction of a prior distribution", "Try different values", type = "error")
    stop("Error identifying control arm prior distribution: ",
         "Stop because we cannot guarantee the accuracy of the numerical integration")
  }
}

#
check_interval_valid_for_a <- function(mode, percentile25) {
  # end-points of the interval to be searched for the root
  q0.5 <- beta_percentile25(a = 0.5, mode) - percentile25
  q50 <- beta_percentile25(a = 50, mode) - percentile25
  
  !identical(sign(q0.5), sign(q50))
}

#
beta_percentile25 <- function(a, mode) {
  b <- (a - 1)/mode - (a - 2)
  
  if (b <= 0) {
    shinyalert::shinyalert("This combination of Q1 and Q2 does not permit the construction of a prior distribution", "Try different values", type = "error")
    stop("Beta distribution b parameter is less than zero. Try different a or mode.")
  }
  
  qbeta(0.25, shape1 = a, shape2 = b, lower.tail = TRUE) 
}


#' Identify parameters of a prior normal distribution for theta~N(mu, sigma2)
#'
#' @param pi1 P(pE > pC) as identified by Day 1 Elicitation Q3
#' @param gamma P(pE - pC > -margin) as identified by 1 - Q4
#' @param a,b parameters of prior distribution for pC
#' @param margin non-inferiority margin for the trial
#'
#' @return vector with param[1] = mu (prior mean of theta), param[2] = sigma2 (prior variance of theta).
#' @export
#'
prior_theta <- function(pi1, gamma, a, b, margin){
  mu_sigma = qnorm(pi1, mean=0, sd=1, lower.tail=TRUE)
  fval = vector(mode="numeric", length=2)
  
  low <- 0.01
  upp <- 5
  
  fval[1] = calc_thetavar(sigmainv = low, mu_sigma, a, b, gamma, margin)
  fval[2] = calc_thetavar(sigmainv = upp, mu_sigma, a, b, gamma, margin)
  
  if(identical(sign(fval[1]), sign(fval[2]))){
    shinyalert::shinyalert("This combination of Q3 and Q4 does not permit the construction of a prior distribution", "Try different values", type = "error")
    stop("Given Beta prior for pC and answers to elicitation questions Q3 and Q4, cannot determine a Normal prior distribution for log-odds ratio",
         "Please revise either the answer to elicitation Q3 or Q4.")
  }
  
  ## Search variable here is 1/sigma, which is 1 over the prior sd of theta.
  z = uniroot(f = calc_thetavar, interval=c(low, upp), mu_sigma, a, b, gamma, margin,
              lower = low, upper=upp,
              f.lower=fval[1], f.upper=fval[2])
  param = vector(mode="numeric", length=2)
  
  ## theta ~ normal(mu, sigma^2)
  mu <- mu_sigma/z$root
  sigma2 <- (1/z$root)^2
  
  c(mu = mu,
    sigma2 = sigma2)
}


#' Search for value of 1/sigma such that P{pE - pC > -margin} = gamma as elicited by 1 - Q4
#'
#' Using Simpsons integration to integrate over the joint prior distribution of (pE, pC).
#' 
#' @param sigmainv 1/sigma where sigma is the hypothesized prior sd of theta
#' @param mu_sigma mu/sigma, ratio of the prior mean and sd of theta
#' @param a,b parameters of prior distribution of pC 
#' @param gamma value of P{pE - pC > -margin} elicited form the expert (=1-q4)
#' @param margin non-inferiority margin for the trial
#'
#' @return P(pE - pC > -margin; 1/sigma) - gamma.
#' @export
#'
calc_thetavar <- function(sigmainv, mu_sigma, a, b, gamma, margin){
  
  ## Set up grid of equally spaced points over [0,1] for pC and over [max{0, pc - margin}, 1] for pE
  gridc = seq(0.001, 0.999, by=0.001)
  midp1 = (0.00001 + 0.001)*0.5
  midp2 = (0.99999 + 0.999)*0.5
  gridc = append(c(0.00001, midp1), gridc)
  gridc = append(gridc, c(midp2, 0.99999))  
  lc = length(gridc)
  
  wc  = vector(mode="numeric", length=lc)
  wc[1] = (gridc[3] - gridc[1])/6
  wc[lc] = (gridc[lc] - gridc[lc-2])/6
  
  for(i in seq(2, lc-1, by=2)){
    wc[i] = 4*(gridc[i+1] - gridc[i-1])/6
  }
  for(i in seq(3, lc-2, by=2)){
    wc[i] = (gridc[i+2] - gridc[i-2])/6
  }
  
  int = vector(mode = "numeric", length = lc)
  int1 = 0
  
  for(i in seq_len(lc)){
    ## for each mesh point for pc, integrate joint density over the interval [max{0, pc - margin}, 1]
    upp = max(0, gridc[i] - margin)
    
    if(upp <= 0){
      int[i] = dbeta(gridc[i], shape1=a, shape2=b, ncp=0, log=FALSE)
    }else{		
      ## create a mesh for pE over the interval [upp, 1] (which always has >= 3 points)
      m = floor((0.999 - upp)/0.002)
      u1 = upp + m*0.002
      midp1 = (u1 + 0.99999)/2
      gride = as.double(append(seq(upp, u1, by = 0.001),
                               append(c(midp1), 0.99999) ))
      le = length(gride)
      
      if(floor(le/2.0) == (le/2.0)){
        stop("mesh for Experimental remission rate contains even number of elements when odd number are expected.")	
      }
      
      we = vector(mode="numeric", length=le)
      we[1] = (gride[3] - gride[1])/6
      we[le] = (gride[le] - gride[le-2])/6
      
      for(j in seq(2, le-1, by=2)){
        we[j] = 4*(gride[j+1] - gride[j-1])/6
      }
      
      for(j in seq(3, le-2, by=2)){
        we[j] = (gride[j+2] - gride[j-2])/6
      }
      
      ## Evaluating joint density for (pc, pe) at vector of values of pE for each given gridc[i]
      dens = vector(mode="numeric", length=le)
      dens1 = vector(mode="numeric", length=le)   
      dens = (gridc[i]^(a-1))*((1-gridc[i])^(b-1))/(gride*(1-gride))
      dens1 = -0.5*(( sigmainv*log(gride*(1-gridc[i])/(gridc[i]*(1-gride))) - mu_sigma)^2)
      dens = dens*exp(dens1)
      int[i] = sum(we*dens)
    }
  }
  
  int1 = sum(wc*int)*sigmainv/(beta(a,b)*sqrt(2*pi))
  
  return(int1 - gamma)
}


#' Evaluate the marginal prior distribution of pE
#'
#' @param a,b parameters of prior distribution of pC 
#' @param mu,sigma2 parameters of prior distribution of theta 
#'
#' @return dataframe containing prior expectation, mode, sd,
#'    limits of 90% credibility interval and 25th percentile
#' @export
#'
prior_e <- function(a, b, mu, sigma2){
  
  ## Use a mesh for integrating pE which should be adequate if P(0.001 <= pE <= 0.999) >= 0.998 
  ## and we check this condition is satisfied later.
  gride = seq(0.001, 0.999, by=0.001)
  midp1 = (0.00001 + 0.001)*0.5
  midp2 = (0.99999 + 0.999)*0.5
  gride = append(c(0.00001, midp1), gride)
  gride = append(gride, c(midp2, 0.99999))
  le = length(gride)
  we = vector(mode="numeric", length=le)
  we[1] = (gride[3] - gride[1])/6
  we[le] = (gride[le] - gride[le-2])/6
  
  for(i in seq(2, le-1, by=2)){
    we[i] = 4*(gride[i+1] - gride[i-1])/6
  }
  for(i in seq(3, le-2, by=2)){
    we[i] = (gride[i+2] - gride[i-2])/6
  }
  
  ## Use a fine mesh to integrate the joint prior density (pE, pC) over pC
  gridc = seq(0.001, 0.999, by=0.001)
  midp1 = (0.00001 + 0.001)*0.5
  midp2 = (0.99999 + 0.999)*0.5
  gridc = append(c(0.00001, midp1), gridc)
  gridc = append(gridc, c(midp2, 0.99999))
  lc = length(gridc)
  wc = vector(mode="numeric", length=lc)
  wc[1] = (gridc[3]-gridc[1])/6
  wc[lc] = (gridc[lc] - gridc[lc-2])/6
  
  for(i in seq(2, lc-1, by=2)){
    wc[i] = 4*(gridc[i+1] - gridc[i-1])/6
  }
  for(i in seq(3, lc-2, by=2)){
    wc[i] = (gridc[i+2] - gridc[i-2])/6
  }
  
  dens = vector(mode="numeric", length=lc)
  dens1 = vector(mode="numeric", length=lc)
  int = vector(mode="numeric", length=le)
  
  #?
  for(i in 1:le){
    dens = (gridc^(a-1))*((1-gridc)^(b-1))/(gride[i]*(1-gride[i]))
    dens1 = (-0.5/sigma2)*((log(gride[i]*(1-gridc)/(gridc*(1-gride[i]))) - mu)^2)
    dens = dens*exp(dens1)
    int[i] = sum(wc*dens)    	
  }
  int = int/(beta(a,b)*sqrt(2*pi*sigma2))  
  
  if (check_shape_pE_U_or_L(int, gride)) {
    shinyalert::shinyalert("This combination of Q1 and Q2 does not permit the construction of a prior distribution", "Try different values", type = "error")
    stop("Prior density Experimental remission rate is U (or L) shaped function of pE. Can't guarantee accuracy of numerical integration routines")
  } else {
    ## calculate the prior mean and variance of pe
    expect = sum(we*gride*int)
    sd1 = sqrt(sum(we*gride*gride*int) - (expect^2))
    modecal = data.frame(gride, int)
    modecal_sort = modecal[order(int, decreasing = TRUE), ]
    mode1 = modecal_sort$gride[1]
    
    ## Call a search routine to find the 5th and 95th percentiles of the distribution 
    fval = vector(mode="numeric", length=2)
    fval[1] = prob_pE(0.005, 0.05, 1, a, b, mu, sigma2, 1, 1, 1, 1, 0)
    fval[2] = prob_pE(0.9, 0.05, 1, a, b, mu, sigma2, 1, 1, 1, 1, 0)
    z = uniroot(prob_pE, interval=c(0.005, 0.9), 0.05, 1, a, b, mu, sigma2, 1, 1, 1, 1, 0,
                lower = 0.005, upper=0.9, f.lower=fval[1], f.upper=fval[2])
    ci_low = z$root
    
    fval[1] = prob_pE(0.1, 0.95, 1, a, b, mu, sigma2, 1, 1, 1, 1, 0)
    fval[2] = prob_pE(0.995, 0.95, 1, a, b, mu, sigma2, 1, 1, 1, 1, 0)
    z = uniroot(prob_pE, interval=c(0.1, 0.995), 0.95, 1, a, b, mu, sigma2, 1, 1, 1, 1, 0,
                lower = 0.1, upper=0.995, f.lower=fval[1], f.upper=fval[2])
    ci_upp = z$root
    
    ## Call a search routine to find the 25th percentile of the distribution
    fval[1] = prob_pE(0.005, 0.25, 1, a, b, mu, sigma2, 1, 1, 1, 1, 0)
    fval[2] = prob_pE(0.9, 0.25, 1, a, b, mu, sigma2, 1, 1, 1, 1, 0)
    z = uniroot(prob_pE, interval=c(0.005, 0.9), 0.25, 1, a, b, mu, sigma2, 1, 1, 1, 1, 0,
                lower = 0.005, upper=0.9, f.lower=fval[1], f.upper=fval[2])
    percent25 = z$root
    
    data.frame(expect, mode1, sd1, ci_low, ci_upp, percent25)
  }
}

#' check_shape_pE_U_or_L
#' 
#' Check to see whether marginal prior density for pE is a U or L shaped function of pE
#' Integrate prior pE density over the interval [lim1, lim2].
#' If probability in this interval is less than would be the case under a flat density, we conclude
#' the density could be U or L shaped. 
#' 
#' @param int 
#' @param gride 
#' 
check_shape_pE_U_or_L <- function(int, gride) {
  
  lim1 = 0.001
  lim2 = 0.999
  
  ## Create a grid for pE covering only the interval (lim1, lim2)
  gride2 = seq(0.001, 0.999, by=0.001)
  le2 = length(gride2)
  we2 = vector(mode="numeric", length=le2)
  we2[1] = (gride2[3] - gride2[1])/6
  we2[le2] = (gride2[le2] - gride2[le2-2])/6
  
  for(i in seq(2, le2-1, by=2)){
    we2[i] = 4*(gride2[i+1] - gride2[i-1])/6
  }
  for(i in seq(3, le2-2, by=2)){
    we2[i] = (gride2[i+2] - gride2[i-2])/6
  }
  
  int2 = int[which(gride >= lim1 & gride <= lim2)] 	
  istop = sum(we2*int2)
  
  istop < (lim2 - lim1)
}



#' Calculate variance of prior distribution of log[pc/(1-pc)]
#'
#' @param a,b Parameters of prior distribution of pC
#'
#' @return prior variance of log[pC/(1-pC)] (used for evaluating prior ESS of log-odds)
#' @export
#'
logoddspc <- function(a, b){
  ## set up a grid for theta1 = log[pc/(1-pc)]
  r = as.integer(64)
  mesh = as.integer(6*r - 1)
  mesh1 = as.integer(2*mesh-1)
  grid1 = vector(mode="numeric", length= mesh)
  gridt = vector(mode="numeric", length= mesh1)
  
  ## centre mesh at log(E(pc)/(1-E(pc))), where E(pc) is prior mean of pC~Beta(a,b)
  mu = a/(a+b)
  mu = log(mu/(1-mu)) 
  
  for(i in 1:mesh){
    if(i <= (r-1)){
      grid1[i] = mu + sqrt(3.0)*(-3-4*log(r/i))
    }else if((i >= r) & (i<= 5*r)){
      grid1[i] = mu + sqrt(3.0)*(-3 + 3*(i-r)/(2*r))
    }else{
      grid1[i] = mu + sqrt(3.0)*(3 + 4*log(r/(6*r - i)))
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
  
  dens = vector(mode="numeric", length = mesh1)
  dens1 = vector(mode="numeric", length = mesh1)
  dens = exp(gridt)/(1+exp(gridt))
  dens1 = (dens^a)/((1+exp(gridt))^b)
  dens1 = dens1/beta(a,b)
  expect = sum(gridt*wtheta*dens1)
  
  sum(wtheta*gridt*gridt*dens1) - (expect^2)
}
