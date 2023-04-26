
#' Probability pE
#' 
#' Functions for calculating
#' \deqn{P{pE <= q}, P{pE <= q|data}, P{pC <= q|data}, P{theta <= q|data}}.
#'
#' @param q hypothesized 100*prob percentile of the distribution.	
#' @param prob Probability
#' @param norm normalising constant of joint posterior distribution of (pE, pC) given data
#' @param a,b,mu,sigma2 parameters of priors pC~Beta(a,b), theta~N(mu, sigma2) 
#' @param sc,fc,s,fe number of successes and failures on treatments CYC and MMF 
#' @param posterior indicator (posterior =1 implies inferences concern posterior distribution)
#'
#' @return \deqn{P{pE <= q} - prob}.
#' @export
#'
prob_pE <- function(q, prob, norm, a, b, mu, sigma2, sc, fc, se, fe, posterior){
	
	## Use a fine mesh to integrate the joint prior density (pE, pC) over pC.
	gridc = seq(0.001, 0.999, by=0.001)
    midp1 = (0.00001 + 0.001)*0.5
    midp2 = (0.99999 + 0.999)*0.5
    gridc = append(c(0.00001, midp1), gridc)
    gridc = append(gridc, c(midp2, 0.99999))
  	lc = length(gridc)
  	if(floor(lc/2.0) == (lc/2.0)){
  		stop("Error calculating quantiles of MMF distribution: integration mesh contains even number of elements when odd number are expected.")
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

	## Set up a grid to integrate over [0, q] with midpoints
	x1 = floor((q - 0.001)/0.004)
	uppbd = 0.001 + x1*0.004
	if(x1 == ((q-0.001)/0.004)){
		gride = seq(0.001, q, by = 0.002)
	}else{
		midpt = (uppbd + q)/2.0
		gride = append(seq(0.001, uppbd, by=0.002), c(midpt, q))	
	}
	
	le = length(gride)
	if(floor(le/2) == (le/2.0)){
		stop("Error calculating quantiles of MMF distribution: integration mesh contains even number of elements when odd number are expected.")
	}
	## Calculate Simpson's integration weights
	we  = vector(mode="numeric", length=le)
  	we[1] = (gride[3]-gride[1])/6.0
  	we[le] = (gride[le] - gride[le-2])/6.0
  	if(le-1 >= 2){
 		for(i in seq(2, (le-1), by=2)){
    		we[i] = 4*(gride[i+1] - gride[i-1])/6.0
 		}
 	}
 	if(le-2 >= 3){
 		for(i in seq(3, (le-2), by=2)){
    		we[i] = (gride[i+2] - gride[i-2])/6.0
 		}
 	}

	dens =  vector(mode="numeric", length=lc)
  	dens1 = vector(mode="numeric", length=lc)
  	int = vector(mode="numeric", length=le)
  	
	if(posterior == 0){
		## Evaluate P{pE <= q}
		
 		## Evaluate the prior density of pE at the mesh points. 
  		int =  vector(mode="numeric", length=le)
   		for(i in 1:le){
      		dens = (gridc^(a-1))*((1-gridc)^(b-1))/(gride[i]*(1-gride[i]))
      		dens1 = (-0.5/sigma2)*((log(gride[i]*(1-gridc)/(gridc*(1-gride[i]))) - mu)^2)
      		dens = dens*exp(dens1)
      		int[i] = sum(wc*dens) 	
   		}
   		int = int/(beta(a,b)*sqrt(2*pi*sigma2)) 
	}else{
		## Evaluate P{pE <= q|data}
		for(i in 1:le){
      		dens = (gridc^(sc + a -1))*((1-gridc)^(fc + b-1))
      		dens1 = (-0.5/sigma2)*((log(gride[i]*(1-gridc)/(gridc*(1-gride[i]))) - mu)^2)
      		dens = 10000*10000*(gride[i]^(se-1))*((1-gride[i])^(fe-1))*dens*exp(dens1)
      		int[i] = sum(wc*dens) 	
   		}
   		int = norm*int
	}
	
	return(sum(we*int) - prob)
}

## Function input:	q = hypothesized 100*prob percentile of the distribution.	
##					norm = normalisation constant of joint posterior distribution of (pC, pE)
##					a, b, mu, sigma2 = parameters of priors pC~Beta(a,b), theta~N(mu, sigma2)
##					sc, fc, se, fe = number of successes and failures on treatments CYC and MMF
##					posterior = indicator (posterior =1 implies inferences concern posterior distribution)
## Function returns P{pC <= q|data} - prob.

#' Title
#'
#' @param q 
#' @param prob 
#' @param norm 
#' @param se,sc,fe,fc 
#' @param a,b 
#' @param mu,sigma2 
#'
#' @return
#' @export
#'
prob_pC <-function(q, prob, norm, se, sc, fe, fc, a, b, mu, sigma2){
	
	gride = seq(0.001, 0.999, by=0.001)
    midp1 = (0.00001 + 0.001)*0.5
    midp2 = (0.99999 + 0.999)*0.5
    gride = append(c(0.00001, midp1), gride)
    gride = append(gride, c(midp2, 0.99999))
  	le = length(gride)
  	
  	we  = vector(mode="numeric", length=le)
  	we[1] = (gride[3]-gride[1])/6.0
  	we[le] = (gride[le] - gride[le-2])/6.0
  	for(i in seq(2, (le-1), by=2)){
    	we[i] = 4*(gride[i+1] - gride[i-1])/6.0
 	}
	for(i in seq(3, (le-2), by=2)){
    	we[i] = (gride[i+2] - gride[i-2])/6.0
	}
	dens =  vector(mode="numeric", length=le)
  	dens1 = vector(mode="numeric", length=le)

	## Set up a grid to integrate pC over [0, q] with midpoints
	
	x1 = floor((q - 0.001)/0.004)
	uppbd = 0.001 + x1*0.004
	if(x1 == ((q-0.001)/0.004)){
		gridc = seq(0.001, q, by = 0.002)
	}else{
		midpt = (uppbd + q)/2.0
		gridc = append(seq(0.001, uppbd, by=0.002), c(midpt, q))	
	}
	lc = length(gridc)
	## Calculate Simpson's integration weights
	wc  = vector(mode="numeric", length=lc)
  	wc[1] = (gridc[3]-gridc[1])/6.0
  	wc[lc] = (gridc[lc] - gridc[lc-2])/6.0
  	if(lc-1 >= 2){
 		for(i in seq(2, (lc-1), by=2)){
    		wc[i] = 4*(gridc[i+1] - gridc[i-1])/6.0
 		}
 	}
 	if(lc -2 >= 3){
 		for(i in seq(3, (lc-2), by=2)){
    		wc[i] = (gridc[i+2] - gridc[i-2])/6.0
 		}
 	}
 	
 	int = vector(mode="numeric", length=lc)
 	## integrating out pe from joint posterior distribution of (pe, pc)
   	for(i in 1:lc){
      	dens =  (gride^(se-1))*((1-gride)^(fe-1))
      	dens1 = (-0.5/sigma2)*((log(gride*(1-gridc[i])/(gridc[i]*(1-gride))) - mu)^2)
      	dens = 10000*10000*(gridc[i]^(sc + a-1))*((1-gridc[i])^(fc + b-1))*dens*exp(dens1)
      	int[i] = sum(we*dens)
   	}
  	## multiply by the normalisation constant
   	int = norm*int
	return(sum(wc*int) - prob)
}

## Function input:	q = hypothesized 100*prob percentile of the distribution.	
##					norm = normalisation constant of joint posterior distribution of (pC, theta)
##					sc, fc, se, fe = number of successes and failures on treatments CYC and MMF
##					a, b, mu, sigma2 = parameters of priors pC~Beta(a,b), theta~N(mu, sigma2)
##					postexp, postsd = expectation and standard deviation of posterior distribution of theta. 
## Function returns P{theta <= q|data} - prob.

#' Title
#'
#' @param q 
#' @param prob 
#' @param norm 
#' @param se,sc,fe,fc 
#' @param a,b 
#' @param mu,sigma2 
#' @param postexp,postsd 
#'
#' @return
#' @export
#'
prob_theta <- function(q, prob, norm, se, sc, fe, fc, a, b, mu, sigma2, postexp, postsd){
	
	gridc = seq(0.001, 0.999, by=0.001)
    midp1 = (0.00001 + 0.001)*0.5
    midp2 = (0.99999 + 0.999)*0.5
    gridc = append(c(0.00001, midp1), gridc)
    gridc = append(gridc, c(midp2, 0.99999))
  	lc = length(gridc)
  	
  	wc  = vector(mode="numeric", length=lc)
  	wc[1] = (gridc[3]-gridc[1])/6.0
  	wc[lc] = (gridc[lc] - gridc[lc-2])/6.0
 	  for(i in seq(2, (lc-1), by=2)){
 		   wc[i] = 4*(gridc[i+1] - gridc[i-1])/6.0
   	}
   	for(i in seq(3, (lc-2), by=2)){
    	wc[i] = (gridc[i+2] - gridc[i-2])/6.0
   	}

	r = as.integer(16)
   	mesh = as.integer(6*r -1)
   	index = as.integer(mesh)
   	grid1 = vector(mode="numeric", length= mesh)
   	## Construct a mesh for integrating over theta which is efficient for integrating 
   	## over the posterior density (assuming it is approximately normal)
   	for(i in 1:mesh){
   		if(i <= (r-1)){
   			grid1[i] = postexp + postsd*(-3-4*log(r/i))
   		}else if((i >= r) & (i<= 5*r)){
   			grid1[i] = postexp + postsd*(-3 + 3*(i-r)/(2*r))
   		}else{
   			grid1[i] = postexp + postsd*(3+ 4*log(r/(6*r -i)))
   		}
   	}  
   	
   	for(i in mesh:1){
   		if(grid1[i] > q){
   			grid1[i] = q
   			index = i
   		}
   	}
   	
   	if(index == 1){
   		grid1[1] = q - 0.1
   		grid1[2] = q
   		index = 2
   	}
   	
   	mesh1 = 2*index -1
   	gridt = vector(mode="numeric", length= mesh1)
   	## Keep elements 1:index of grid1 
   	for(i in seq(1, index, by=1)){
   		gridt[2*i-1] = grid1[i]	
   	}
   	for(i in seq(2, mesh1-1, by=2)){
   		gridt[i] = (gridt[i+1] + gridt[i-1])/2.0
   	}
   	   
	## calculating Simpson's integration weights
   	wtheta  = vector(mode="numeric", length=mesh1)
    wtheta[1] = (gridt[3]-gridt[1])/6.0
  	wtheta[mesh1] = (gridt[mesh1] - gridt[mesh1-2])/6.0
  	if(mesh1 -1 >= 2){
  		for(i in seq(2, (mesh1-1), by=2)){
    		wtheta[i] = 4*(gridt[i+1] - gridt[i-1])/6.0
   		}
   	}
   	if(mesh1-2 >= 3){
   		for(i in seq(3, (mesh1-2), by=2)){
			wtheta[i] = (gridt[i+2] - gridt[i-2])/6.0
   		}
   	} 
   	
	dens = vector(mode = "numeric", length=lc)
	dens1 = vector(mode = "numeric", length=lc)
	int =  vector(mode = "numeric", length=mesh1)
	## integrating over pc for the joint posterior density of (pc, theta)
	for(i in 1:mesh1){
		dens1 = (gridc^(sc-fe+a))*((1-gridc)^(fc+b-2+fe))
		v = (1 + exp(-gridt[i])*(1-gridc)/gridc)*((1 + exp(gridt[i])*gridc/(1-gridc))^(1/(se+fe-1)))
		v = v^(se+fe-1)
		dens = 10000*10000*dens1/v
		int[i] = sum(wc*dens)*exp((1-fe)*gridt[i])*exp((-0.5/sigma2)*((gridt[i] - mu)^2)) 		
	}
	
	int = norm*int
	return(sum(wtheta*int) - prob)
}
