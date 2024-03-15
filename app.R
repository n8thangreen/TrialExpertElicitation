
#install.packages("shiny")
library(shiny)

#install.packages("rriskDistributions")
library(rriskDistributions)

#install.packages("zipfR")
library(zipfR)

#install.packages("cubature")
library(cubature)


############
# functions

#
to_prob_scale <- function(x)
  x/sum(x)

#
integrand_theta <- function(x, a, b, theta, n, ...) {
  (x^(a-1)*(1-x)^(b-1))/((1-x) + x*exp(theta))^n
}

#
integrand_or <- function(x, a, b, or, n, ...) {
  (x^(a-1)*(1-x)^(b-1))/((1-x) + x*or)^n
}

#
set_integral <- function(integrand, lower = 0, upper = 1) {
  force(integrand)
  force(lower)
  force(upper)
  
  function(...) {
    args <- c(f = integrand,
              lower = lower,
              upper = upper,
              list(...))
    do.call(integrate, args)
  }
}

#
integrate_logor <- set_integral(integrand = integrand_theta)
integrate_or <- set_integral(integrand = integrand_or)

#
post_logOR <- function(a, b, theta, n, sE, mu, sigma) {
  args <- tibble::lst(a, b, theta, n)
  num_integral <- do.call(integrate_logor, args)$value
  exp(theta*sE) * exp(-((theta - mu)^2)/(2*sigma^2)) * num_integral
}

#
post_OR <- function(a, b, or, n, sE, mu, sigma) {
  args <- tibble::lst(a, b, or, n)
  num_integral <- do.call(integrate_or, args)$value
  or^(sE-1) * exp(-((log(or) - mu)^2)/(2*sigma^2)) * num_integral
}

#
make_one_over_K_times_ESS <- function(a, b, mu, stdev) {
  force(a)
  force(b)
  force(mu)
  force(stdev)
  
  one_over_K_integrand <- make_one_over_K_integrand(a, b, mu, stdev)
  
  function(x) {
    pE <- x[1]
    pC <- x[2]
    
    one_over_K_integrand(c(pE, pC))*(stdev^2)*(1/16)*(pE + pC)*(2-pE-pC)
  }
}

#
calc_ESS <- function(a, b, mu, stdev) {
  K <- calc_K(a, b, mu, stdev)
  one_over_K_ESS_integrand <- make_one_over_K_times_ESS(a, b, mu, stdev)
  
  res <- cubature::adaptIntegrate(one_over_K_ESS_integrand,
                                  lowerLimit = c(0,0), upperLimit = c(1,1))
  1/(res$integral*K)
}

#
make_one_over_K_integrand <- function(a, b, mu, stdev) {
  force(a)
  force(b)
  force(mu)
  force(stdev)
  
  function(x) {
    pE <- x[1]
    pC <- x[2]
    
    ((pC^(a-1))*((1-pC)^(b-1))/(pE*(1-pE))*
        exp((-1/(2*(stdev^2))) * ((log((pE*(1-pC))/(pC*(1-pE))) - mu)^2)))
    
  }
}

#
calc_K <- function(a, b, mu, stdev) {
  K_integrand <- make_one_over_K_integrand(a, b, mu, stdev)
  res_K <- cubature::adaptIntegrate(K_integrand, lowerLimit = c(0,0), upperLimit = c(1,1))
  
  1/res_K$integral
}


######
# app

ui <- fluidPage(
  
  # App title ----
  titlePanel("BAR-JDM workshop"),
  
  # Sidebar layout with inputs and outputs ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      tabsetPanel( 
        
        tabPanel("Priors",
                 br(),
                 sliderInput(inputId = "Q1",
                             label = "Q1: I think the response rate on methotrexate will be:",
                             min = 0,
                             max = 1,
                             value = 0.3),
                 
                 sliderInput(inputId = "Q2",
                             label = "Q2: I am 75% sure that the response rate on methotrexate will be less than:",
                             min = 0,
                             max = 1,
                             value = 0.5),
                 
                 br(),
                 
                 sliderInput(inputId = "Q3",
                             label = "Q3: I think that the odds ratio on response for baricitinib relative to methotrexate will exceed 1.0 with probability:",
                             min = 0,
                             max = 100,
                             step = 1,
                             post = " %",
                             value = 5),
                 
                 sliderInput(inputId = "Q4",
                             label = "Q4: I think that the odds ratio on response for baricitinib relative to methotrexate will exceed 2.0 with probability:",
                             min = 0,
                             max = 100,
                             step = 1,
                             post = " %",
                             value = 0)
                 
        ) ,
        
        tabPanel("Hypothetical trial",
                 br(),
                 sliderInput(inputId = "Q5_mtx",
                             label = "Qa: Sample size of methotrexate arm:",
                             min = 0,
                             max = 100,
                             value = 10), 
                 
                 sliderInput(inputId = "Q5_bar",
                             label = "Qb: Sample size of baricitinib arm:",
                             min = 0,
                             max = 100,
                             value = 10), 
                 
                 sliderInput(inputId = "Q6",
                             label = "Qc: Number of successes in methotrexate arm:",
                             min = 0,
                             max = 100,
                             value = 5),
                 
                 sliderInput(inputId = "Q7",
                             label = "Qd: Number of successes in baricitinib arm:",
                             min = 0,
                             max = 100,
                             value = 5)
        ) 
      )
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      tabsetPanel(
        
        tabPanel("Priors",
                 
                 tabsetPanel(
                   
                   # panel 1
                   tabPanel("Control arm prior", plotOutput(outputId = "control_prior")
                            , textOutput("mode")
                            , textOutput("proba_lower_p75")
                            , textOutput("alpha")
                            , textOutput("beta")
                   )
                   ,
                   
                   # panel 2
                   tabPanel("Relative effect prior",
                            plotOutput(outputId = "logOR_prior")
                            , plotOutput(outputId = "OR_prior")
                            , textOutput("ESS")
                            , textOutput("logOR_mean_sd")
                   )
                   ,
                   
                   # panel 3
                   tabPanel("Experimental arm prior",
                            plotOutput(outputId = "marginal_prior_plot")
                            , textOutput("experimental_mode")
                            , textOutput("experimental_proba_lower_p75")
                            , textOutput("experimental_alpha")
                            , textOutput("experimental_beta")
                            )
                   ,         
                   # panel 4
                   tabPanel("All priors"
                            , plotOutput(outputId = "all_priors")
                            #, plotOutput(outputId = "joint_posterior_plot")
                   )
                   
                   
                 ) 
                 
        ), # end of tabset panel for Priors
        
        tabPanel("Posteriors", 
                 
                 tabsetPanel(
                   
                   tabPanel("Control arm posterior",
                            plotOutput(outputId = "marginal_posterior_plot_pc")
                            #, plotOutput(outputId = "joint_posterior_plot")
                            , textOutput("mode_posterior")
                   ),
                   
                   # panel 2
                   tabPanel("Relative effect posterior",
                            plotOutput(outputId = "logOR_posterior")
                            , plotOutput(outputId = "OR_posterior")
                   )
                   ,
                   tabPanel("Experimental arm posterior",
                            plotOutput(outputId = "marginal_posterior_plot_pe"),
                            textOutput("experimental_mode_posterior")
                            # plotOutput(outputId = "joint_posterior_plot")),
                   ),
                   
                   tabPanel("All posteriors", plotOutput(outputId = "all_posteriors")
                   )
                 )    
        )
      ) 
    ) # end of main panel
  ) # end of sidebarLayout
) # end of UI
#)


# Define server logic ----
server <- function(input, output) {
  
  observe({
    updateSliderInput(inputId = "Q6", max = input$Q5_mtx)
  })

  observe({
    updateSliderInput(inputId = "Q7", max = input$Q5_bar)
  })
  
  # a and b parameters of the beta distribution for control arm as reactives
  
  control_beta_a <- reactive({
    mode <- input$Q1
    p75 <- input$Q2
    
    min.SS <- function(params) {
      alpha=params[1]
      beta=params[2]
      
      ((alpha-1)/(alpha+beta-2) - mode)^2 + (Rbeta(p75,a=alpha,b=beta)-0.75)^2
      
    }
    
    optim_result <- optim(par = c(1.5, 2), fn = min.SS, lower=c(1,1), method="L-BFGS-B")
    optim_result$par[1]
  })
  
  control_beta_b <- reactive({
    mode <- input$Q1
    p75 <- input$Q2
    
    min.SS <- function(params) {
      alpha <- params[1]
      beta <- params[2]
      
      ((alpha-1)/(alpha+beta-2) - mode)^2 + (Rbeta(p75,a=alpha,b=beta)-0.75)^2
    }
    
    optim_result <- optim(par = c(1.5, 2), fn = min.SS, lower=c(1,1), method="L-BFGS-B")
    optim_result$par[2]
  })
  
  
  ## Control arm prior plot
  output$control_prior <- renderPlot({
    
    # define range
    p <- seq(0,1, length=100)
    prior_val <- dbeta(p, control_beta_a(), control_beta_b())
    prior_c <- to_prob_scale(prior_val)
    cdf <- cumsum(prior_c)
    mode_c <- round(p[prior_c == max(prior_c)], 2)
    min75 <- round(min(p[cdf >= 0.75]), 2)
    
    # create plot of corresponding Beta distribution
    plot(p, prior_val,
         ylab='density', xlab="Control arm response rate",
         main='Control arm prior density for response rate',
         type='l', lwd=1.5, col='red', xaxt='n')
    polygon(c(p[cdf >= 0.75], 1, min(p[cdf >= 0.75])), c(prior_val[cdf >= 0.75], 0, 0), col="lightblue", border=NA)
    abline(v = mode_c, lty = 2, col = "black")
    
    axis(1, at = sort(c(min75, mode_c, seq(0, 1, 0.2))), label=TRUE)
  })
  
  ## Probability greater than mode as reactive
  
  proba_greater_mode <- reactive({
    
    proba <- pbeta(input$Q1, shape1=control_beta_a(), shape2=control_beta_b(), lower.tail = FALSE)
    round(proba, 2)*100
  })
  
  # output of the probability greater than mode
  output$proba_greater_mode <- renderText({ 
    paste("The parameters for Q1 and Q2 that you have selected correspond to a probability beyond the mode of", proba_greater_mode(),"%")
  })
  
  ## Probability lower p75 as reactive
  
  proba_lower_p75 <- reactive({
    
    proba_p75 <- pbeta(input$Q2, shape1=control_beta_a(), shape2=control_beta_b(), lower.tail = TRUE)
    round(proba_p75, 2)*100
  })
  
  output$mode <- renderText({ 
    glue::glue("The mode, which is the value that you think is most likely for the response rate in the Control arm, is {input$Q1}. ",
               "It is indicated by a vertical dashed line on the plot.")
  })
  
  # output of the probability lower p75
  output$proba_lower_p75 <- renderText({ 
    glue::glue("The blue area to the right of {input$Q2} is 25%. It means that you think that there is one chance out of 4 ",
               "that the response rate will be greater than {input$Q2}.") 
  })

  output$ESS <- renderText({
    ess <- calc_ESS(a = control_beta_a(), b = control_beta_b(),
                    mu = mean_normal(), stdev = sd_normal())
    glue::glue("The Effective Sample Size (ESS) is {round(ess,0)}") 
  })
  
  output$mode_posterior <- renderText({
    dens <- to_prob_scale(marginal_posterior_density_pc())
    mode_c <- round(pe[dens == max(dens)], 2)
    
    glue::glue("The mode, which is the value that you think is most likely for the response rate in the Control arm, is {mode_c}.")
  })
  
  output$experimental_mode_posterior <- renderText({
    dens <- to_prob_scale(marginal_posterior_density_pe())
    mode_c <- round(pe[dens == max(dens)], 2)
    
    glue::glue("The mode, which is the value that you think is most likely for the response rate in the Experimental arm, is {mode_c}.")
  })
  
  # output alpha and beta parameters
  
  output$alpha <- renderText({ 
    # paste("The mode and p75 that you have selected correspond to alpha of", round(control_beta_a(),2))
  })
  
  output$beta <- renderText({ 
    # paste("The mode and p75 that you have selected correspond to beta of", round(control_beta_b(),2))
  })

  output$experimental_mode <- renderText({
    dens <- to_prob_scale(marginal_prior_density())
    mode_c <- round(pe[dens == max(dens)], 2)
    
    glue::glue("The mode, which is the value that you think is most likely for the response rate in the Experimental arm, is {mode_c}. ",
               "It is indicated by a vertical dashed line on the plot.")
  })
  
  output$experimental_proba_lower_p75 <- renderText({
    dens <- to_prob_scale(marginal_prior_density())
    cdf <- cumsum(dens)
    min75 <- round(min(pe[cdf >= 0.75]), 2)
    
    glue::glue("The blue area to the right of {min75} is 25%. It means that you think that there is one chance out of 4 ",
               "that the response rate will be greater than {min75}.") 
  })
  
  output$experimental_alpha <- renderText({ 
    # paste("The mode and p75 that you have selected correspond to alpha of", round(control_beta_a(),2))
  })
  
  output$experimental_beta <- renderText({ 
    # paste("The mode and p75 that you have selected correspond to beta of", round(control_beta_b(),2))
  })
  
  output$logOR_mean_sd <- renderText({
    
    glue::glue("The mean of the log odds ratio is {round(mean_normal(),2)}. ",
               "The standard deviation of the log odds ratio is {round(sd_normal(),2)}")
  })
  
  
  logOR_normal_params <- reactive({
    p1 <- 1 - odd_ratio_gt_1()
    p2 <- 1 - odd_ratio_gt_1.x()
    
    # get.norm.par(p = c(p1, p2), q = c(0, log(1.25)), plot=FALSE, show.output=FALSE)
    get.norm.par(p = c(p1, p2), q = c(0, log(2)), plot=FALSE, show.output=FALSE)
  })

  mean_normal <- reactive({
    as.numeric(logOR_normal_params()[1])
  })
  
  sd_normal <- reactive({
    as.numeric(logOR_normal_params()[2])
  })
  
  # Relative effect prior plots (log OR and OR)
  output$logOR_prior <- renderPlot({
    
    # range for log normal
    p <- seq(-2.5, 2.5, length=100)
    
    plogOR <- to_prob_scale(dnorm(p, mean=mean_normal(), sd=sd_normal()))
    
    # corresponding Normal distribution for log OR
    plot(p, plogOR, ylab='density',
         type ='l', lwd=1.5, xlab="log OR",col='red', main='Prior density for the log odds ratio')
    abline(v = 0, lty = 2)
    abline(v = log(2), lty = 2)
  })
  
  output$OR_prior <- renderPlot({
    
    # range for OR
    x <- seq(-10, 1.5, length=1000)
    
    pOR <- to_prob_scale(exp(dnorm(x, mean=mean_normal(), sd=sd_normal())))
    
    # corresponding distribution for OR
    plot(exp(x), pOR, ylab='density',
         type ='l', lwd=1.5, xlab="OR", col='red', main='Prior density for the odds ratio')
    abline(v = 1, lty = 2)
    abline(v = 2, lty = 2)
  })
  
  # Relative effect prior plots (log OR and OR)
  output$logOR_posterior <- renderPlot({
    theta <- marginal_posterior_density_theta()
    
    plot(theta$val, theta$dens, ylab='density',
         type ='l', lwd=1.5, xlab="Log-odds ratio", 
         col='black', main='Log-odds ratio posterior density')
    abline(v = 0, lty = 2)
  })
  
  output$OR_posterior <- renderPlot({
    etheta <- marginal_posterior_density_OR()
    
    plot(etheta$val, etheta$dens, ylab='density',
         type ='l', lwd=1.5, xlab="Odds ratio", 
         col='black', main='Odds ratio posterior density')
    abline(v = 1, lty = 2)
    abline(v = 2, lty = 2)
  })
  
  #### Experimental arm prior 
  pc <- seq(0.01, 0.99, by = 0.01)
  pe <- seq(0.01, 0.99, by = 0.01)
  
  joint_prior_density <- reactive({ 
    
    z <- matrix(data=NA, nrow=length(pc), ncol=length(pe))

    for(i in 1:length(pe))
    {
      for(j in 1:length(pc))
      {
        term1 <- ((pc[j]^(control_beta_a()-1))*((1-pc[j])^(control_beta_b()-1)))/(pe[i]*(1-pe[i]))
        
        odds_ratio <- (pe[i]*(1-pc[j])) / (pc[j]*(1-pe[i]))
        
        term2 <- exp((-1/(2*sd_normal()^2))*((log(odds_ratio) - mean_normal())^2))
        
        z[i,j] <- term1*term2  
      }
    }
    z
  })
  
  marginal_prior_density <- reactive({
    
    marginal_density <- 100 * to_prob_scale(rowSums(joint_prior_density()))
    #marginal_density <- rowSums(joint_prior_density())
  })
  
  output$marginal_prior_plot <- renderPlot({
    
    dens_val <- marginal_prior_density()
    dens <- to_prob_scale(dens_val)
    
    cdf <- cumsum(dens)
    mode_c <- round(pe[dens == max(dens)], 2)
    min75 <- round(min(pe[cdf >= 0.75]), 2)
    
    plot(pe, dens_val,
         ylab='density', xlab="Experimental arm response rate", 
         main='Experimental arm prior density for response rate',
         type ='l', lwd=1.5, col='purple')
    polygon(c(pe[cdf >= 0.75], 1, min(pe[cdf >= 0.75])), c(dens_val[cdf >= 0.75], 0, 0), col="lightblue", border=NA)
    abline(v = mode_c, lty = 2, col = "black")
  })
  
  
  #https://www.tutorialspoint.com/how-to-create-a-function-in-r-with-two-inputs#:~:text=To%20create%20a%20function%20with,x%20and%20y%20inside%20function.
  #http://www.countbio.com/web_pages/left_object/R_for_biology/R_fundamentals/3D_surface_plot_R.html
  
  output$all_priors <- renderPlot({
    
    # define range
    # p = seq(-2.5,2.5, length=100)
    x <- seq(-10, 1.5, length=1000)
    p <- seq(0, 1, length=100)
    
    par(mfrow=c(2,2))
    
    #create plot of corresponding Beta distribution
    plot(p, to_prob_scale(dbeta(p, control_beta_a(), control_beta_b())), ylab='density',
         type ='l', lwd=1.5, xlab="Control arm response rate", col='red', main='Control arm prior density for response rate')
    
    plot(pe, to_prob_scale(marginal_prior_density()), ylab='density',
         type ='l', lwd=1.5, xlab="Experimental arm response rate", 
         col='purple', main='Experimental arm prior density for response rate')
    
    plot(exp(x), to_prob_scale(exp(dnorm(x, mean=mean_normal(), sd=sd_normal()))), ylab='density',
         type ='l', lwd=1.5, xlab="OR", col='red', main='Prior density for the odds ratio')
    
    plot(p, to_prob_scale(dbeta(p, control_beta_a(), control_beta_b())), ylab='density',
         type ='l', lwd=1.5, xlab="Response rate", col='red', main='Prior densities for response rate')
    lines(pe, to_prob_scale(marginal_prior_density()), lwd=1.5,col='purple')

    ## save figure to file    
    # png("figures/prior_plot.png")
    # plot(p, to_prob_scale(dbeta(p, control_beta_a(), control_beta_b())), ylab='Density',
    #      type ='l', lwd=2, xlab="Response rate", col='green')
    # lines(pe, to_prob_scale(marginal_prior_density()), lty=2, lwd=2, col='red')
    # legend("topleft", legend = c("MTX", "BARI"), col = c("green", "red"), lty = c(1,2))
    # dev.off()
    
    #logOR_prior <- plot(p, dnorm(p, mean=mean_normal(), sd=sd_normal()), ylab='density',
    #                    type ='l', lwd=1.5, xlab="log OR",col='red', main='Prior density for the log odds ratio')
  })
  
  ############################################################ 
  ## hypothetical trial and posteriors
  
  odd_ratio_gt_1 <- reactive({ 
    input$Q3/100
  })
  
  odd_ratio_gt_1.x <- reactive({ 
    input$Q4/100
  })
  
  N_sample_control <- reactive({ 
    input$Q5_mtx
  })
  
  N_sample_exp <- reactive({ 
    input$Q5_bar
  })
  
  success_control <- reactive({ 
    input$Q6
  })
  
  success_exp <- reactive({ 
    input$Q7
  })
  
  failures_control <- reactive({ 
    N_sample_control() - success_control()
  })
  
  failures_exp <- reactive({ 
    N_sample_exp() - success_exp()
  })
  
  
  joint_posterior_density <- reactive({ 
    
    # empty matrix
    x = matrix(data=NA, nrow=length(pc), ncol=length(pe))
    
    for(i in 1:length(pe))
    {
      for(j in 1:length(pc))
      {
        term1=((pc[j]^(input$Q6+control_beta_a()-1))*((1-pc[j])^(failures_control() + control_beta_b()-1)))*((pe[i]^(input$Q7-1))*((1-pe[i])^(failures_exp()-1) ))
        
        term2=exp((-1/(2*sd_normal()^2))*((log((pe[i]*(1-pc[j]))/(pc[j]*(1-pe[i]))) - mean_normal())^2))
        
        x[i,j] = term1*term2  
      }
    }
    x
  })
  
  # output$joint_posterior_plot <- renderPlot({
  #   
  #   persp(pe, pc, joint_posterior_density(), #shade=0.5, 
  #         theta=30, phi=25,
  #         axes=TRUE, scale=TRUE , box=TRUE, expand = 0.5,
  #         nticks=3, #ticktype="detailed", 
  #         xlim = c(0,1), ylim=c(0,1), col="lightgreen",
  #         xlab="Control \narm", ylab="Experimental \narm", zlab="density", main="Joint posterior density")
  # })
  
  marginal_posterior_density_pe <- reactive({
    to_prob_scale(rowSums(joint_posterior_density()))
  })
  
  marginal_posterior_density_pc <- reactive({
    to_prob_scale(colSums(joint_posterior_density()))
  })
  
  # log OR
  marginal_posterior_density_theta <- reactive({
    val <- seq(-2.5, 2.5, length=1000)
    success <- success_exp() + success_control()
    fail <- failures_exp() + failures_control()
    
    dens <- purrr::map_dbl(val, \(x) post_logOR(a = success + control_beta_a(), b = fail + control_beta_b(),
                                                theta = x, n = N_sample_exp(), sE = success_exp(),
                                                mu = mean_normal(), sigma = sd_normal()))
    tibble::lst(val, dens)
  })
  
  # OR
  marginal_posterior_density_OR <- reactive({
    val <- exp(seq(-10, 1.5, length=1000))
    success <- success_exp() + success_control()
    fail <- failures_exp() + failures_control()
    
    dens <- purrr::map_dbl(val, \(x) post_OR(a = success + control_beta_a(), b = fail + control_beta_b(),
                                             or = x, n = N_sample_exp(), sE = success_exp(),
                                             mu = mean_normal(), sigma = sd_normal()))
    tibble::lst(val, dens)
  })
  
  output$marginal_posterior_plot_pe <- renderPlot({
    plot(pe, marginal_posterior_density_pe(),
         xlab="Experimental arm response rate", ylab='density',
         type ='l', lwd=1.5, col='purple', main='Experimental arm posterior density for response rate')
  })
  
  output$marginal_posterior_plot_pc <- renderPlot({
    plot(pc, marginal_posterior_density_pc(),
         xlab="Control arm response rate", ylab='density',
         type ='l', lwd=1.5, col='red', main='Control arm posterior density for response rate')
  })
  
  
  output$all_posteriors <- renderPlot({
    p <- seq(0,1, length=100)
    
    par(mfrow=c(2,3)) 
    
    prior_c <- dbeta(p, control_beta_a(), control_beta_b())
    plot(pc, marginal_posterior_density_pc(), ylab='density',
         type ='l', lwd=1.5, xlab="Control arm response rate", 
         col='red', main='Control arm posterior density for response rate')
    lines(p, to_prob_scale(prior_c), lwd=1.5, lty = 2, col='red')
    
    plot(pe, marginal_posterior_density_pe(), ylab='density',
         type ='l', lwd=1.5, xlab="Experimental arm response rate", 
         col='purple', main='Experimental arm posterior density for response rate')
    lines(pe, to_prob_scale(marginal_prior_density()), lwd=1.5, lty = 2, col='purple')
    
    plot(pc, marginal_posterior_density_pc(), ylab='density',
         type ='l', lwd=1.5, xlab="Control arm response rate", 
         col='red', main='Control arm posterior density for response rate')
    lines(pe, marginal_posterior_density_pe(), type ='l', lwd=1.5, col='purple')
    
    # log OR
    theta <- marginal_posterior_density_theta()
    plot(theta$val, theta$dens, ylab='density',
         type ='l', lwd=1.5, xlab="Log-odds ratio", 
         col='black', main='Log-odds ratio posterior density')
    
    # OR
    etheta <- marginal_posterior_density_OR()
    plot(etheta$val, etheta$dens, ylab='density',
         type ='l', lwd=1.5, xlab="Odds ratio", 
         col='black', main='Odds ratio posterior density')
  })

}

shinyApp(ui = ui, server = server)

#rsconnect::deployApp()
