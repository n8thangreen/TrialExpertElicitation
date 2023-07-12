## Function to calculate summaries of the prior distributions of pC, pE and theta
## Function inputs: q1, ..., q4 = expert's answers to the Day 1 elicitation questions. 
##					expert =  a character string giving the expert's initials

## 					q1 is the answer to Day 1 elicitation question (i) mode(pc)
## 					q2 is the answer to Day 1 elicitation question (ii) eliciting q s.t. P(pC > q) = 0.75
## 					q3 is the answer to Day 1 elicitation question (iii) elicting P(theta >0)
## 					q4 is the answer to Day 1 elicitation question (iv) eliciting P(pE - pC < -c2)

## Function outputs: a vector summarising properties of the prior distributions of pC, pE and theta
## x[1], x[2] parameters of beta prior distribution for pC, x[3] = ESS of log(pC/(1-pC)) given elicited prior for pC
## x[4] = E(pE), x[5] = mode(pE), x[6] = SD(pE), (x[7], x[8]) = 90% credibility interval, x[14] = quantile such that P(pE <= x[14]) = 0.25
## x[9] = E(theta), x[10] = var(theta), x[11] = P(pE > pC), x[12] = P(pE < pC - 0.1), x[13] = ESS of theta

## The function also outputs a file saved as "expert-priorplot.pdf" in the current working directory containing plots of the prior densities
## and a file "expert-D1answer.txt" storing the expert's answers to the elicitation questions.  

#' priorcall
#'
#' @param q1 
#' @param q2 
#' @param q3 
#' @param q4 
#' @param expert 
#'
#' @return
#' @export
#'
priorcall <- function(q1, q2, q3, q4, expert = "", out_dir = "plots"){
  
  c2 <- 0.1 	# non-inferiority margin cited in Day 1 elicitation question (iv)
  
  # Numerical searches to find the parameters of the prior distributions
  # for pC and theta assuming our statistical model holds 
  bparam = prior_beta(q1, q2)
  tparam = prior_theta(q3, 1-q4, bparam[1], bparam[2], c2)
  
  ## Calculating summaries of prior distributions for pC, pE and theta
  x = vector(mode="numeric", length=14)
  pri_e = prior_e(bparam[1], bparam[2], tparam[1], tparam[2])
  
  ## Summaries of pC prior distribution
  x[1] = bparam[1]
  x[2] = bparam[2]
  x[3] = ess_pc(bparam[1], bparam[2])
  
  ## Summaries of pE prior distribution
  x[4] = pri_e$expect
  x[5] = pri_e$mode1
  x[6] = pri_e$sd1
  x[7] = pri_e$ci_low
  x[8] = pri_e$ci_upp
  x[14] = pri_e$percent25
  
  ## Summaries of theta prior distribution	
  x[9] = tparam[1]
  x[10] = tparam[2]
  x[11] = pnorm(0, mean = tparam[1], sd = sqrt(tparam[2]), lower.tail = FALSE)
  x[12] = q4
  x[13] = ess_theta(1, tparam[1], tparam[2], bparam[1], bparam[2])
  
  ## generate a plot of the prior densities and output to one pdf file
  z = vector(mode="numeric", length =20)
  y1 = distPlot(0,0,0,0, z, x, as.character("pC"), 1)
  y2 = distPlot(0,0,0,0, z, x, as.character("pE"), 1)
  y3 = distPlot(0,0,0,0, z, x, as.character("theta"), 1)
  
  outputfile <- system.file(paste0(out_dir, "/", expert, "-priorplot.pdf"), package = "TrialExpertElicitation", mustWork = TRUE)
  
  title1 <- paste0(expert, "'s prior density of control arm remission rate")
  title2 <- paste0(expert, "'s prior density of MMF remission rate")
  title3 <- paste0(expert, "'s prior density of log-odds ratio")
  title4 <- as.character("Prior density of control arm & MMF remission rate")
  
  pdf(outputfile)
  
  par(mfrow = c(2,2), pty="s")
  
  plot(y1$gridc, y1$dens, type="l", lty=1, lwd=3, col="red", main = title1, xlab = "CYC 6-month remission rate", ylab="Density", xlim =c(0,1), cex.lab = 1.1, cex.axis=1.1, cex.main = 1) 
  plot(y2$gride, y2$dens, type="l", lty=1, lwd=3, col="green", main = title2, xlab = "MMF 6-month remission rate", ylab="Density", xlim =c(0,1), cex.lab = 1.1, cex.axis=1.1, cex.main = 1)
  plot(y3$gridt, y3$dens, type="l", lty=1, lwd=3, col="blue", main = title3, xlab = "log-odds ratio", ylab="Density", xlim =c(-2, 2), cex.lab = 1.1, cex.axis=1.1, cex.main = 1)
  plot(y1$gridc, y1$dens, type="l", lty=1, lwd=3, col="red", main = title4, xlab = "6-month remission rate", ylab="Density", ylim = range(c(y1$dens, y2$dens)), xlim =c(0,1), cex.lab = 1.1, cex.axis=1.1, cex.main = 1)
  lines(y2$gride, y2$dens, type="l", lty=2, lwd=3, col="green")
  dev.off()
  
  ## Write answers to elicitation questions to file 
  outputfile <- system.file(paste0(out_dir, "/", expert, "-D1answer.txt"), package = "TrialExpertElicitation", mustWork = TRUE)
  
  cat(expert, "'s answers to Day 1 prior elicitation questions were: \n", file = outputfile, append=FALSE)
  cat(q1, file = outputfile, sep="\n", append=TRUE)
  cat(q2, file = outputfile, sep="\n", append=TRUE)
  cat(q3, file = outputfile, sep="\n", append=TRUE)
  cat(q4, file = outputfile, sep="\n", append=TRUE)
  
  x
}

## Function to calculate summaries of the posterior distributions of pC, pE and theta
## Function inputs: 	n_mmf, mmf_succ = (nE, SE): number of patients randomised to MMF and number of observed successes on MMF
##						n_cyc, cyc_succ = (nC, SC): number of patients randomised to CYC and number of observed success on CYC
##						priorParm = vector outputted by priorcall() containing summaries of prior distributions
##						posterior40 = logical variable indicating whether wish to calculate posterior distribution assuming nE+nC = 40
##										(if false assume wish to calculate posterior distribution assuming nE + nC = 20)

##	x[1] = E(pC|data), x[2] = mode(pC|data), x[3] = SD(pC|data), (x[4], x[5]) = 90% posterior credibility interval for pC, x[18] = normalising constant of g(pC, pE|data)
##	x[6] = E(pE|data), x[7] = mode(pE|data), x[8] = SD(pE|data), (x[9], x[10]) = 90% posterior credibility interval for pE, x[19] = normalising constant of g(pC, pE|data)
##	x[11] = E(theta|data), x[12] = mode(theta|data), x[13] = SD(theta|data), (x[14], x[15]) = 90% posterior credibility interval for theta, 
## 	x[16] = P{pE > pC|data}, x[17] = P{pE - pC > -c2|data},  x[20] = normalising constant of joint posterior distribution f(theta, pC|data)

#' Summaries of posterior distributions given hypothetical dataset
#'
#' @param n_mmf 
#' @param mmf_succ 
#' @param n_cyc 
#' @param cyc_succ 
#' @param priorParm 
#' @param posterior40 is posterior 40 selected? Logical
#'
#' @return
#' @export
#'
postSumry <- function(n_mmf, mmf_succ, n_cyc, cyc_succ, priorParm, posterior40){
  
  c2 = as.double(0.10) 		## non-inferiority margin for the trial
  x = vector(mode="numeric", length =20)
  
  ## Catching possible input errors for mmf_succ, n_cyc and cyc_succ
  if(mmf_succ >= n_mmf){
    stop("Data on MMF: number of successes exceeds number randomised to MMF")
  }
  
  if(posterior40 & !identical(as.integer(n_mmf + n_cyc), as.integer(40))){
    stop("Total number randomized to MMF and control arm does not sum to 40")
  }
  if(!posterior40 & !identical(as.integer(n_mmf + n_cyc), as.integer(20))){
    stop("Total number randomized to MMF and control arm does not sum to 20")
  }
  
  pc_distn = post_pc(priorParm[1], priorParm[2], priorParm[9], priorParm[10], cyc_succ, n_cyc-cyc_succ, mmf_succ, n_mmf-mmf_succ)	
  pe_distn = post_pe(priorParm[1], priorParm[2], priorParm[9], priorParm[10], cyc_succ, n_cyc-cyc_succ, mmf_succ, n_mmf-mmf_succ)
  theta_distn = post_theta(priorParm[1], priorParm[2], cyc_succ, mmf_succ, n_cyc-cyc_succ, n_mmf-mmf_succ, priorParm[9], priorParm[10], c2)
  
  x[1] = pc_distn$expect
  x[2] = pc_distn$mode1
  x[3] = pc_distn$sd1
  x[4] = pc_distn$ci_low
  x[5] = pc_distn$ci_upp
  x[18] = pc_distn$norm
  
  ## Summaries of the posterior distribution of pE
  x[6] = pe_distn$expect
  x[7] = pe_distn$mode1
  x[8] = pe_distn$sd1
  x[9] = pe_distn$ci_low
  x[10] = pe_distn$ci_upp
  x[19] = pe_distn$norm
  
  ## Summaries of the posterior distribution of theta
  x[11] = theta_distn$expect
  x[12] = theta_distn$mode1
  x[13] = theta_distn$sd1
  x[14] = theta_distn$ci_low
  x[15] = theta_distn$ci_upp
  x[16] = calc_pi(mmf_succ, priorParm[9], priorParm[10], x[11], (x[13])^2, priorParm[1], priorParm[2], n_mmf, cyc_succ, n_cyc-cyc_succ)
  x[17] = calc_gamma(mmf_succ, n_mmf - mmf_succ, cyc_succ, n_cyc-cyc_succ, priorParm[1], priorParm[2], priorParm[9], priorParm[10], x[19], c2)
  x[20] = theta_distn$norm
  
  x
}


## Function inputs: 	n_mmf = number of patients randomised to E (nE)
## 						mmf_succ = number of successes on E (SE)
## 						n_cyc = number randomised to C (nC)
## 						cyc_succ = number of successes on C (SC)
## 						postParm = vector returned by postSumry() summarising posterior distributions of pC, pE, theta
## 						priorParm = vector of outputs returned by priorcall() summarising elicited prior distributions
## 						parmInd = character string containing the parameter whose density we wish to plot 
## 						postind = integer determining whether we wish to plot the prior or posterior density of the stated parameter.
## Function outputs: a dataframe containing a grid of values of the parameter of interest and the marginal density evaluated at those values.

#' Prior and posterior densities of pC, pE and theta ready for plotting
#'
#' @param n_mmf 
#' @param mmf_succ 
#' @param n_cyc 
#' @param cyc_succ 
#' @param postParm 
#' @param priorParm 
#' @param parmInd 
#' @param postind 
#'
#' @return
#' @export
#'
distPlot <- function(n_mmf, mmf_succ, n_cyc, cyc_succ, postParm, priorParm, parmInd, postind){
  
  if(postind){
    ## Wish to generate prior density for the stated parameter
    if(parmInd == "pC"){
      ## We only need evaluate the prior density of pC
      gridc = seq(0.01, 0.99, by=0.01)
      gridc = append(c(0.00001, 0.001), gridc)
      gridc = append(gridc, c(0.999, 0.99999))
      dens = dbeta(gridc, shape1=priorParm[1], shape2=priorParm[2], ncp=0, log=FALSE)
      return(data.frame(gridc, dens))
    }else if(parmInd=="pE"){
      gride = seq(0.01, 0.99, by=0.01)
      gride = append(c(0.00001, 0.001), gride)
      gride = append(gride, c(0.999, 0.99999))
      le = length(gride)
      
      gridc = seq(0.0001, 0.9999, by=0.0001)
      midp1 = (0.00001 + 0.0001)*0.5
      midp2 = (0.99999 + 0.9999)*0.5
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
      
      dens1 =  vector(mode="numeric", length=lc)
      dens2 = vector(mode="numeric", length=lc)
      dens =  vector(mode="numeric", length=le)
      for(i in 1:le){
        dens1 = (gridc^(priorParm[1]-1))*((1-gridc)^(priorParm[2]-1))/(gride[i]*(1-gride[i]))
        dens2 = (-0.5/priorParm[10])*((log(gride[i]*(1-gridc)/(gridc*(1-gride[i]))) - priorParm[9])^2)
        dens1 = dens1*exp(dens2)
        dens[i] = sum(wc*dens1)    		
      }
      dens = dens/(beta(priorParm[1],priorParm[2])*sqrt(2*pi*priorParm[10])) 
      return(data.frame(gride, dens))		
    }else{	
      gridt = seq(-2, 2, by = 0.05)	
      dens = dnorm(gridt, mean=priorParm[9], sd=sqrt(priorParm[10]), log=FALSE)
      return(data.frame(gridt, dens))	
    }	
  }else{
    ## Wish to generate posterior density for the stated parameter
    if(parmInd=="pC"){
      gridc = seq(0.01, 0.99, by=0.01)
      gridc = append(c(0.00001, 0.001), gridc)
      gridc = append(gridc, c(0.999, 0.99999))	
      postd = pc_dens(gridc, n_mmf - mmf_succ, mmf_succ, n_cyc - cyc_succ, cyc_succ, priorParm, postParm[18])
      return(data.frame(gridc, postd))
    }else if(parmInd == "pE"){
      gride = seq(0.01, 0.99, by=0.01)
      gride = append(c(0.00001, 0.001), gride)
      gride = append(gride, c(0.999, 0.99999))
      postd = pe_dens(gride, n_mmf - mmf_succ, mmf_succ, n_cyc -cyc_succ, cyc_succ, priorParm, postParm[19])
      return(data.frame(gride, postd))
    }else{
      gridt = seq(-2, 2, by = 0.05)
      postd = theta_dens(gridt, n_mmf - mmf_succ, mmf_succ, n_cyc - cyc_succ, cyc_succ, priorParm, postParm[20])
      return(data.frame(gridt, postd))
    }
  }		
}


## Function inputs: 	gridc = vector of values of pC
##						fe, se, fc, sc = number of successes and failures on MMF and CYC
##						priorParm = output of priorcall() summarising elicited prior distributions
##						norm = normalising constant of joint posterior distribution g(pC, pE|data)
## Function output: posterior density of pC eveluated at gridc	

#' Posterior density of pC
#'
#' @param gridc 
#' @param fe 
#' @param se 
#' @param fc 
#' @param sc 
#' @param priorParm 
#' @param norm 
#'
#' @return
#' @export
#'
pc_dens <- function(gridc, fe, se, fc, sc, priorParm, norm){
  
  lc = length(gridc)
  gride = seq(0.0001, 0.9999, by=0.0001)
  midp1 = (0.00001 + 0.0001)*0.5
  midp2 = (0.99999 + 0.9999)*0.5
  gride = append(c(0.00001, midp1), gride)
  gride = append(gride, c(midp2, 0.99999))
  le = length(gride)
  we  = vector(mode="numeric", length=le)
  we[1] = (gride[3]-gride[1])/6.0
  we[le] = (gride[le] - gride[le-2])/6.0
  for(i in seq(2,(le-1), by=2)){
    we[i] = 4*(gride[i+1] - gride[i-1])/6.0
  }
  for(i in seq(3, (le-2), by=2)){
    we[i] = (gride[i+2] - gride[i-2])/6.0
  }
  
  dens =  vector(mode="numeric", length=le)
  dens1 = vector(mode="numeric", length=le)
  postd =  vector(mode="numeric", length=lc)
  ## integrating out pe from joint posterior distribution of (pe, pc)
  for(i in 1:lc){
    dens =  (gride^(se-1))*((1-gride)^(fe-1))
    dens1 = (-0.5/priorParm[10])*((log(gride*(1-gridc[i])/(gridc[i]*(1-gride))) - priorParm[9])^2)
    dens = 10000*10000*(gridc[i]^(sc + priorParm[1]-1))*((1-gridc[i])^(fc + priorParm[2]-1))*dens*exp(dens1)
    postd[i] = sum(we*dens)
  }
  postd = norm*postd
  
  postd
}


## Function inputs: 	gride = vector of values of pE
##						priorParm = output of priorcall() summarising elicited prior distributions
##						norm = normalising constant of joint posterior distribution g(pC, pE|data)
## Function output: posterior density of pE eveluated at gride	

#' Posterior density of pE
#'
#' @param gride 
#' @param fe,se,fc,sc number of successes and failures on MMF and CYC
#' @param priorParm 
#' @param norm 
#'
#' @return
#' @export
#'
pe_dens <- function(gride, fe, se, fc, sc, priorParm, norm){
  
  le = length(gride)
  gridc = seq(0.0001, 0.9999, by=0.0001)
  midp1 = (0.00001 + 0.0001)*0.5
  midp2 = (0.99999 + 0.9999)*0.5
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
  
  
  pdens =  vector(mode="numeric", length=lc)
  pdens1 = vector(mode="numeric", length=lc)
  postd =  vector(mode="numeric", length=le)
  for(i in 1:le){
    pdens = (gridc^(sc + priorParm[1] -1))*((1-gridc)^(fc + priorParm[2] - 1))
    pdens1 = (-0.5/priorParm[10])*((log(gride[i]*(1-gridc)/(gridc*(1-gride[i]))) - priorParm[9])^2)
    pdens = 10000*10000*(gride[i]^(se-1))*((1-gride[i])^(fe-1))*pdens*exp(pdens1)
    postd[i] = sum(wc*pdens) 	
  }
  postd = norm*postd
  
  postd
}

#' Evaluate posterior density of theta
#'  
#' @param gridt vector of values of theta
#' @param fe,se,fc,sc number of successes and failures on MMF and CYC 
#' @param priorParm output of priorcall() summarising elicited prior distributions
#' @param norm normalising constant of joint posterior distribution g(pC, theta|data)
#'
#' @return posterior density of theta evaluated at gridt	
#' @export
#'
theta_dens <- function(gridt, fe, se, fc, sc, priorParm, norm){
  
  lt = length(gridt)
  gridc = seq(0.0001, 0.9999, by=0.0001)
  midp1 = (0.00001 + 0.0001)*0.5
  midp2 = (0.99999 + 0.9999)*0.5
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
  
  dens = vector(mode = "numeric", length=lc)
  dens1 = vector(mode = "numeric", length=lc)
  postd =  vector(mode = "numeric", length=lt)
  for(i in 1:lt){
    dens1 = (gridc^(sc-fe+priorParm[1]))*((1-gridc)^(fc+priorParm[2]-2+fe))
    v = (1 + exp(-gridt[i])*(1-gridc)/gridc)*((1 + exp(gridt[i])*gridc/(1-gridc))^(1/(se+fe-1)))
    v = v^(se+fe-1)
    dens = 10000*10000*dens1/v
    postd[i] = sum(wc*dens)*exp((1-fe)*gridt[i])*exp((-0.5/priorParm[10])*((gridt[i] - priorParm[9])^2)) 				
  }
  postd = norm*postd
  
  postd
}

