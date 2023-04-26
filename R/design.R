
#' Calculate pi
#' 
#' @param se  number of successes observed on MMF in the dataset
#' @param mu,sigma2 parameters of prior distribution of theta
#' @param pmean,pvar posterior mean and variance of theta.
#' @param a,b parameters of prior distribution of pC
#' @param ne,sc,fc number randomised to MMF; number of successes and failures on CYC.
#'
#' @return \deqn{P{pE > pC|data}}
#' @export
#'
calc_pi <- function(se, mu, sigma2, pmean, pvar, a, b, ne, sc, fc) {
  
  fe = double(1)
  fe = ne - se
  norm = as.double(0)
  
  r = as.integer(16)
  mesh = as.integer(6*r -1)	
  
  for(k in 1:2){	
    grid1 = vector(mode="numeric", length= mesh)
    if(k==1){
      calcnorm = as.integer(1)
    }else{
      calcnorm = as.integer(0)
    }
    
    ## set up mesh for theta
    for(i in 1:mesh){
      if(i <= (r-1)){
        grid1[i] = pmean + sqrt(pvar)*(-3-4*log(r/i))
      }else if((i >= r) & (i<= 5*r)){
        grid1[i] = pmean + sqrt(pvar)*(-3 + 3*(i-r)/(2*r))
      }else{
        grid1[i] = pmean + sqrt(pvar)*(3+ 4*log(r/(6*r -i)))
      }
    }
    
    if(calcnorm == 1){
      ## want to calculate normalising constant of joint posterior distribution of (pC, theta)
      crit = as.double(-1000)
    }else{
      crit = as.double(0)
    }
    index = as.integer(1)
    for(i in 1:mesh){
      if(grid1[i] <= crit){
        grid1[i] = crit
        index = i
      }
    }
    if(grid1[mesh] ==crit){
      index = mesh-1
      grid1[mesh] = crit + 0.5
    }
    
    ## put reduced mesh in a new vector retaining only those grid points in the interval [0, infty)
    mesh1 = as.integer(mesh - index+1)
    grid2 = vector(mode="numeric", length= mesh1)
    for(i in 1:mesh1){
      grid2[i] = grid1[i + index-1]
    }
    mesh2 = as.integer(2*mesh1-1)
    gridt = vector(mode="numeric", length= mesh2)
    ## calculating the midpoints of grid2
    if(mesh1 >= 1){
      for(i in seq(1, mesh1, by=1)){
        gridt[2*i-1] = grid2[i]	
      }
    }
    if(mesh2-1 >= 2){
      for(i in seq(2, mesh2-1, by=2)){
        gridt[i] = (gridt[i+1] + gridt[i-1])/2.0
      } 
    }
    
    ## Calculate the Simpsons weights for integrating over this region
    wtheta  = vector(mode="numeric", length=mesh2)
    wtheta[1] = (gridt[3]-gridt[1])/6.0
    wtheta[mesh2] = (gridt[mesh2] - gridt[mesh2-2])/6.0
    for(i in seq(2,(mesh2-1), by=2)){
      wtheta[i] = 4*(gridt[i+1] - gridt[i-1])/6.0
    }
    for(i in seq(3, (mesh2-2), by=2)){
      wtheta[i] = (gridt[i+2] - gridt[i-2])/6.0
    }
    
    ## Evaluate the marginal posterior distribution for theta over [0, infty). Setting up grid for pc and 
    ## integrate the joint posterior distribution for (pc, theta) with respect to pC
    gridc = seq(0.0001, 0.9999, by=0.0001)
    midp1 = (0.00001 + 0.0001)*0.5
    midp2 = (0.99999 + 0.9999)*0.5
    gridc = append(c(0.00001, midp1), gridc)
    gridc = append(gridc, c(midp2, 0.99999))
    lc = length(gridc)
    if(floor(lc/2.0) == (lc/2.0)){
      stop("Error in design.R: integration mesh contains even number of elements when odd number are expected.")
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
    
    dens = vector(mode = "numeric", length=lc)
    dens1 = vector(mode = "numeric", length=lc)
    int =  vector(mode = "numeric", length=mesh2)
    ## integrating over pc for the joint posterior density of (pc, theta)
    for(i in 1:mesh2){
      dens1 = (gridc^(sc-fe+a))*((1-gridc)^(fc+b-2+fe))
      v = 1 + exp(-gridt[i])*(1-gridc)/gridc
      v1 = 1 + exp(gridt[i])*gridc/(1-gridc)
      dens = 10000*10000*dens1/((v^(ne-1))*v1)
      int[i] = sum(wc*dens)*exp((1-fe)*gridt[i])*exp((-0.5/sigma2)*((gridt[i] - mu)^2))	
    }	
    if(calcnorm ==1){
      ## calculating normalising constant of the joint posterior density of (pc, theta)
      norm = sum(wtheta*int)
    }else{
      int = (1/norm)*int
    }
  }
  return(sum(wtheta*int))	
}


#' Calculate gamma
#' 
#' @param se,fe number of successes and failures on CYC
#' @param sc,fc number of successes and failures on MMF
#' @param a,b parameters of prior distribution of pC
#' @param mu,sigma2 parameters of prior distribution of theta
#' @param norm normalising constant of joint posterior distribution of (pC, pE)
#' @param c2 non-inferiority margin for the trial
#'
#' @return P{pE > pC - c2|data}
#' @export
#'
calc_gamma <- function(se, fe, sc, fc, a, b, mu, sigma2, norm, c2) {
  
  gridc = seq(0.001, 0.999, by=0.001)
  midp1 = (0.00001 + 0.001)*0.5
  midp2 = (0.99999 + 0.999)*0.5
  gridc = append(c(0.00001, midp1), gridc)
  gridc = append(gridc, c(midp2, 0.99999))  
  lc = length(gridc)
  
  if(floor(lc/2.0) == (lc/2.0)){
    stop("Error in calc_gamma.R: integration mesh contains even number of elements when odd number are expected.")
  }
  
  wc  = vector(mode="numeric", length=lc)
  wc[1] = (gridc[3]-gridc[1])/6.0
  wc[lc] = (gridc[lc] - gridc[lc-2])/6.0
  
  for(i in seq(2, lc-1, by=2)){
    wc[i] = 4*(gridc[i+1] - gridc[i-1])/6.0
  }
  for(i in seq(3, lc-2, by=2)){
    wc[i] = (gridc[i+2] - gridc[i-2])/6.0
  }
  
  int =  vector(mode="numeric", length=lc)
  
  ## for each mesh point for pc, integrate joint posterior density of (pC, pE)
  ## over the interval [max{0, pc - c2}, 1] with respect to pE
  for(i in 1:lc){
    upp = max(0, gridc[i] - c2)
    if(upp <=0){
      gride = seq(0.0001, 0.9999, by=0.0001)
      midp1 = (0.00001 + 0.0001)*0.5
      midp2 = (0.99999 + 0.9999)*0.5
      gride = append(c(0.00001, midp1), gride)
      gride = append(gride, c(midp2, 0.99999)) 
      le = length(gride)
    }else{
      ## create a mesh for pe over the interval [upp, 1] (which always has >= 3 points)
      m = floor((0.999 - upp)/0.002)
      u1 = upp + m*0.002
      midp1 = (u1 + 0.99999)/2.0 	
      gride = as.double(append(seq(upp, u1, by = 0.001), append(c(midp1), 0.99999) ))
      le = length(gride)
    }
    
    if(floor(le/2.0) == (le/2.0)){
      stop("Error in calc_gamma.R: integration mesh contains even number of elements when odd number are expected.")	
    } 
    we  = vector(mode="numeric", length=le)
    we[1] = (gride[3]-gride[1])/6.0
    we[le] = (gride[le] - gride[le-2])/6.0
    
    for(j in seq(2, (le-1), by=2)){
      we[j] = 4*(gride[j+1] - gride[j-1])/6.0
    }
    for(j in seq(3, (le-2), by=2)){
      we[j] = (gride[j+2] - gride[j-2])/6.0
    }	  
    dens =  vector(mode="numeric", length=le)
    dens1 = vector(mode="numeric", length=le)
    
    ## integrating over the joint posterior density of (pc, pe)  
    dens = (gridc[i]^(sc + a-1))*((1-gridc[i])^(fc+b-1))*(gride^(se-1))*((1-gride)^(fe-1))
    dens1 = (-0.5/sigma2)*((log(gride*(1-gridc[i])/(gridc[i]*(1-gride))) - mu)^2)
    dens = 10000*10000*dens*exp(dens1)
    int[i] = sum(we*dens)
  }
  
  int = int*norm
  
  as.double(sum(wc*int))
}
