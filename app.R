
#install.packages("shiny")
library(shiny)

#install.packages("rriskDistributions")
library(rriskDistributions)

#install.packages("zipfR")
library(zipfR)

# Define UI for app ----
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
                 # Input Q1  ----
                 sliderInput(inputId = "Q1",
                             label = "Q1: What's the mode for the Control arm prior:",
                             min = 0,
                             max = 1,
                             value = 0.3),
                 
                 # Input Q2  ----
                 sliderInput(inputId = "Q2",
                             label = "Q2: What's p75:",
                             min = 0,
                             max = 1,
                             value = 0.5),
                 
                 br(),
                 
                 # Input Q3  ----
                 sliderInput(inputId = "Q3",
                             label = "Q3: Chance that Experimental arm will do better:",
                             min = 0,
                             max = 1,
                             value = 0.5),
                 
                 # Input Q4  ----
                 sliderInput(inputId = "Q4",
                             label = "Q4: Chance that Experimental arm will do at least 25% better:",
                             min = 0,
                             max = 1,
                             value = 0.25)
                 
        ) ,
        
        tabPanel("Hypothetical trial",
                 br(),
                 sliderInput(inputId = "Q5",
                             label = "Q5: Sample size per arm of hypothetical trial:",
                             min = 4,
                             max = 40,
                             value = 20), 
                 
                 sliderInput(inputId = "Q6",
                             label = "Q6: # successes in Control arm:",
                             min = 4,
                             max = 40,
                             value = 5),
                 
                 sliderInput(inputId = "Q7",
                             label = "Q7: # successes in Experimental arm:",
                             min = 4,
                             max = 40,
                             value = 5),
                 
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
                            
                            , textOutput("proba_lower_p75")
                            , textOutput("alpha")
                            , textOutput("beta")
                            , textOutput("mode")
                   )
                   ,
                   
                   # panel 2
                   tabPanel("Relative effect prior", plotOutput(outputId = "logOR_prior")
                            , plotOutput(outputId = "OR_prior")
                   )
                   ,
                   
                   # panel 3
                   tabPanel("Experimental arm prior", plotOutput(outputId = "marginal_prior_plot")
                            , plotOutput(outputId = "joint_prior_plot"))
                   ,         
                   # panel 4
                   tabPanel("All priors"
                            ,plotOutput(outputId = "all_priors")
                            #, plotOutput(outputId = "joint_posterior_plot")
                   )
                   
                   
                 ) 
                 
        ), # end of tabset panel for Priors
        
        tabPanel("Posteriors", 
                 
                 tabsetPanel(
                   
                   tabPanel("Control arm posterior",plotOutput(outputId = "marginal_posterior_plot_pc")
                            #, plotOutput(outputId = "joint_posterior_plot")
                   ),
                   
                   tabPanel("Experimental arm posterior",plotOutput(outputId = "marginal_posterior_plot_pe")
                            , plotOutput(outputId = "joint_posterior_plot")),
                   
                   tabPanel("All posteriors",plotOutput(outputId = "all_posteriors")
                   )
                   
                 )    
        ),
        tabPanel("Summary",
                 
                 tabsetPanel(
                   tabPanel(
                     "Priors",
                     br(),
                     h4("Control arm characteristics given your answers to Q1 and Q2:"),
                     textOutput("mode_prior_control"),
                     textOutput("proba_greater_mode"),
                     br(),
                     h4("Relative effect of experimental treatment given your answers to Q3 and Q4:"),
                     textOutput("positive_trt_effect"),
                     textOutput("positive_trt_effect_25pct")
                   ),
                   
                   tabPanel(
                     "Posteriors", 
                     br(),
                     h4("Hypothetical trial:"),
                     textOutput("trial_size"),
                     textOutput("responders_ctrl"),
                     textOutput("responders_exp"),
                     
                     br(),
                     h4("Posterior distributions given prior distributions and trial results:"),
                     #textOutput("proba_post_exp_greater_20pct"),
                     textOutput("proba_post_exp_greater_20pct_text"),
                     textOutput("proba_post_exp_greater_30pct_text"),
                     textOutput("proba_post_exp_greater_40pct_text"),
                     textOutput("proba_post_exp_greater_50pct_text")
                     
                   )
                 ) # end of tabsetPanel within Summary 
                 
        ) # end of summary tabPanel
        
        
      ) 
      
    ) # end of main panel
  ) # end of sidebarLayout
) # end of UI
#)

# Define server logic ----
server <- function(input, output) {
  
  # a and b parameters of the beta distribution for control arm as reactives
  
  control_beta_a <- reactive({
    
    mode=input$Q1
    p75=input$Q2
    
    min.SS <- function(params) {
      
      alpha=params[1]
      beta=params[2]
      
      ((alpha-1)/(alpha+beta-2) - mode)^2 + (Rbeta(p75,a=alpha,b=beta)-0.75)^2
      
    }
    
    optim_result <- optim(par = c(1.5, 2), fn = min.SS, lower=c(1,1), method="L-BFGS-B")
    optim_result$par[1]
    
  })
  
  control_beta_b <- reactive({
    
    mode=input$Q1
    p75=input$Q2
    
    min.SS <- function(params) {
      
      alpha=params[1]
      beta=params[2]
      
      ((alpha-1)/(alpha+beta-2) - mode)^2 + (Rbeta(p75,a=alpha,b=beta)-0.75)^2
      
    }
    
    optim_result <- optim(par = c(1.5, 2), fn = min.SS, lower=c(1,1), method="L-BFGS-B")
    optim_result$par[2]
    
  })
  
  
  ## Control arm prior plot
  output$control_prior <- renderPlot({
    
    #define range
    p = seq(0,1, length=100)
    
    #create plot of corresponding Beta distribution
    control_prior <- plot(p, dbeta(p, control_beta_a(), control_beta_b()), ylab='density',
                          type ='l', lwd=1.5, xlab="Control arm response rate", col='red', main='Control arm prior density for response rate')
    control_prior
  })
  
  ## Proba greater than mode as reactive
  
  proba_greater_mode <- reactive({
    
    proba <- pbeta(input$Q1, shape1=control_beta_a(), shape2=control_beta_b(),lower.tail = FALSE)
    round(proba,2)*100
  })
  
  # output of the proba greater than mode
  output$proba_greater_mode <- renderText({ 
    paste("The parameters for Q1 and Q2 that you have selected correspond to a probability beyond the mode of", proba_greater_mode(),"%")
  })
  
  ## Proba lower p75 as reactive
  
  proba_lower_p75 <- reactive({
    
    proba_p75 <- pbeta(input$Q2, shape1=control_beta_a(), shape2=control_beta_b(),lower.tail = TRUE)
    round(proba_p75,2)*100
  })
  
  # output of the proba lower p75
  output$proba_lower_p75 <- renderText({ 
    paste("The mode and p75 that you have selected correspond to a probability to the left of p75 of", proba_lower_p75(),"%")
  })
  
  # output alpha and beta parameters
  
  output$alpha <- renderText({ 
    paste("The mode and p75 that you have selected correspond to alpha of", round(control_beta_a(),2))
  })
  
  output$beta <- renderText({ 
    paste("The mode and p75 that you have selected correspond to beta of", round(control_beta_b(),2))
  })
  
  # output the mode
  
  output$mode <- renderText({ 
    paste("The current mode is", round((control_beta_a()-1)/(control_beta_a()+control_beta_b()-2),2))
  })
  
  ## parameters of the normal distribution of log OR as reactives
  
  mean_normal <- reactive({
    p1=1-input$Q3
    p2=1-input$Q4
    
    params <- get.norm.par(p = c(p1, p2), q=c(0,log(1.25)),plot=FALSE, show.output=FALSE)
    mean = as.numeric(params[1])
    #sd = as.numeric(params[2])
    mean
  })
  
  sd_normal <- reactive({
    p1=1-input$Q3
    p2=1-input$Q4
    
    params <- get.norm.par(p = c(p1, p2), q=c(0,log(1.25)),plot=FALSE, show.output=FALSE)
    #mean = as.numeric(params[1])
    sd = as.numeric(params[2])
    sd
  })
  
  # Relative effect prior plots (log OR and OR)
  output$logOR_prior <- renderPlot({
    
    #define range for log normal
    p = seq(-2.5,2.5, length=100)
    
    #create plot of corresponding Normal distribution for log OR
    
    logOR_prior <- plot(p, dnorm(p, mean=mean_normal(), sd=sd_normal()), ylab='density',
                        type ='l', lwd=1.5, xlab="log OR",col='red', main='Prior density for the log odds ratio')
    
    logOR_prior
    
  })
  
  output$OR_prior <- renderPlot({
    
    #define range for OR
    x = seq(-10,1, length=1000)
    
    #create plot of corresponding distribution for OR
    
    OR_prior <- plot(exp(x), exp(dnorm(x, mean=mean_normal(), sd=sd_normal())), ylab='density',
                     type ='l', lwd=1.5, xlab="OR", col='red', main='Prior density for the odds ratio')
    
  })
  
  
  #### Experimental arm prior 
  pc =  seq(0.01,0.99,0.01) # seq(0,1, length=100) #
  pe = seq(0.01,0.99,0.01)
  
  joint_prior_density <- reactive({ 
    
    # An empty matrix z
    z = matrix(data=NA, nrow=length(pc), ncol=length(pe))
    z
    
    for(i in 1:length(pe))
    {
      for(j in 1:length(pc))
      {
        term1=((pc[j]^(control_beta_a()-1))*((1-pc[j])^(control_beta_b()-1)))/(pe[i]*(1-pe[i]))
        
        term2=exp((-1/(2*sd_normal()^2))*((log((pe[i]*(1-pc[j]))/(pc[j]*(1-pe[i]))) - mean_normal())^2))
        
        z[i,j] = term1*term2  
      }
    }
    z
  })
  
  
  output$joint_prior_plot <- renderPlot({
    
    persp(pe,pc,joint_prior_density(), #shade=0.5, 
          theta=30, phi=25,
          axes=TRUE,scale=TRUE , box=TRUE, expand = 0.5,
          nticks=3, #ticktype="detailed", 
          xlim = c(0,1),ylim=c(0,1), col="lightblue", xlab="Control \narm", 
          ylab="Experimental \narm", zlab="density", main="Joint prior density")
  })
  
  marginal_prior_density <- reactive({
    
    marginal_density <- 100 * (rowSums(joint_prior_density())/sum(rowSums(joint_prior_density())) )
    #marginal_density <- rowSums(joint_prior_density())
    
  })
  
  output$marginal_prior_plot <- renderPlot({
    
    marginal_prior <- plot(pe,marginal_prior_density(), ylab='density',
                           type ='l', lwd=1.5, xlab="Experimental arm response rate", 
                           col='purple', main='Experimental arm prior density for response rate')
    marginal_prior
  })
  
  
  #https://www.tutorialspoint.com/how-to-create-a-function-in-r-with-two-inputs#:~:text=To%20create%20a%20function%20with,x%20and%20y%20inside%20function.
  #http://www.countbio.com/web_pages/left_object/R_for_biology/R_fundamentals/3D_surface_plot_R.html
  
  output$all_priors <- renderPlot({
    par(mfrow=c(2,2))
    
    # p = seq(-2.5,2.5, length=100)
    
    #create plot of corresponding Normal distribution for log OR
    
    x = seq(-10,1, length=1000)
    
    #create plot of corresponding distribution for OR
    
    #define range
    p = seq(0,1, length=100)
    
    #create plot of corresponding Beta distribution
    control_prior <- plot(p, dbeta(p, control_beta_a(), control_beta_b()), ylab='density',
                          type ='l', lwd=1.5, xlab="Control arm response rate", col='red', main='Control arm prior density for response rate')
    
    marginal_prior <- plot(pe,marginal_prior_density(), ylab='density',
                           type ='l', lwd=1.5, xlab="Experimental arm response rate", 
                           col='purple', main='Experimental arm prior density for response rate')
    
    OR_prior <- plot(exp(x), exp(dnorm(x, mean=mean_normal(), sd=sd_normal())), ylab='density',
                     type ='l', lwd=1.5, xlab="OR", col='red', main='Prior density for the odds ratio')
    
    plot(p, dbeta(p, control_beta_a(), control_beta_b()), ylab='density',
         type ='l', lwd=1.5, xlab="Response rate", col='red', main='Prior densities for response rate')
    lines(pe,marginal_prior_density(), 
          type ='l', lwd=1.5,col='purple')
    
    #logOR_prior <- plot(p, dnorm(p, mean=mean_normal(), sd=sd_normal()), ylab='density',
    #                    type ='l', lwd=1.5, xlab="log OR",col='red', main='Prior density for the log odds ratio')
    
  })
  
  
  ############################################################################## 
  
  ## hypothetical trial and posteriors
  
  ##############################################################################
  
  failures_control <- reactive({ 
    n_failures_control = input$Q5 - input$Q6
  })
  
  failures_exp <- reactive({ 
    n_failures_exp = input$Q5 - input$Q7
  })
  
  
  
  joint_posterior_density <- reactive({ 
    
    # An empty matrix x
    x = matrix(data=NA, nrow=length(pc), ncol=length(pe))
    x
    
    for(i in 1:length(pe))
    {
      for(j in 1:length(pc))
      {
        term1=((pc[j]^(input$Q6+control_beta_a()-1))*((1-pc[j])^(failures_control()+control_beta_b()-1)))*((pe[i]^(input$Q7-1))*((1-pe[i])^(failures_exp()-1) ))
        
        term2=exp((-1/(2*sd_normal()^2))*((log((pe[i]*(1-pc[j]))/(pc[j]*(1-pe[i]))) - mean_normal())^2))
        
        x[i,j] = term1*term2  
      }
    }
    x
  })
  
  
  output$joint_posterior_plot <- renderPlot({
    
    persp(pe,pc,joint_posterior_density(), #shade=0.5, 
          theta=30, phi=25,
          axes=TRUE,scale=TRUE , box=TRUE, expand = 0.5,
          nticks=3, #ticktype="detailed", 
          xlim = c(0,1),ylim=c(0,1), col="lightgreen", xlab="Control \narm", 
          ylab="Experimental \narm", zlab="density", main="Joint posterior density")
  })
  
  marginal_posterior_density_pe <- reactive({
    
    marginal_density_post_pe <- rowSums(joint_posterior_density())/sum(rowSums(joint_posterior_density()))
    
  })
  
  
  marginal_posterior_density_pc <- reactive({
    
    marginal_density_post_pc <- colSums(joint_posterior_density())/sum(colSums(joint_posterior_density()))
    
  })
  
  output$marginal_posterior_plot_pe <- renderPlot({
    
    marginal_posterior_plot_pe <- plot(pe,marginal_posterior_density_pe(), ylab='density',
                                       type ='l', lwd=1.5, xlab="Experimental arm response rate", 
                                       col='purple', main='Experimental arm posterior density for response rate')
    marginal_posterior_plot_pe
  })
  
  output$marginal_posterior_plot_pc <- renderPlot({
    
    marginal_posterior_plot_pc <- plot(pc,marginal_posterior_density_pc(), ylab='density',
                                       type ='l', lwd=1.5, xlab="Control arm response rate", 
                                       col='red', main='Control arm posterior density for response rate')
    marginal_posterior_plot_pc
  })
  
  
  output$all_posteriors <- renderPlot({
    par(mfrow=c(2,2)) 
    
    plot(pc,marginal_posterior_density_pc(), ylab='density',
         type ='l', lwd=1.5, xlab="Control arm response rate", 
         col='red', main='Control arm posterior density for response rate')
    
    plot(pe,marginal_posterior_density_pe(), ylab='density',
         type ='l', lwd=1.5, xlab="Experimental arm response rate", 
         col='purple', main='Experimental arm posterior density for response rate')
    
    plot(pc,marginal_posterior_density_pc(), ylab='density',
         type ='l', lwd=1.5, xlab="Control arm response rate", 
         col='red', main='Posterior densities for response rate')
    lines(pe,marginal_posterior_density_pe(),type ='l', lwd=1.5,col='purple')
    
  })
  
  
  ############################################################################## 
  
  ## summary stats page
  
  ##############################################################################  
  
  output$mode_prior_control <- renderText({ 
    paste("1/ The current mode, which is the value for the response rate in the Control arm that you think is most probable, is", 100*round((control_beta_a()-1)/(control_beta_a()+control_beta_b()-2),2),"%")
  })
  
  # output of the proba greater than mode
  output$proba_greater_mode <- renderText({ 
    paste("2 / The probability that a response rate greater than the mode will be observed is estimated at", proba_greater_mode(),"%")
  })
  
  # output of chance of an positive treatment effect
  output$positive_trt_effect <- renderText({ 
    paste("1 / The probability that there will be an improvement thanks to the experimental treatment strategy is estimated at", 100*input$Q3,"%")
  })
  
  # output of chance of an positive treatment effect greater than 25 pct
  output$positive_trt_effect_25pct <- renderText({ 
    paste("2 / The chance that this improvement will be greater than 25% is estimated at", 100*input$Q4,"%")
  })
  
  # size of hypothetical trial 
  output$trial_size <- renderText({ 
    paste("1 / The trial has a sample size per arm of", input$Q5,"patients")
  })
  
  output$responders_ctrl <- renderText({ 
    paste("2 / The number of responders in the Control arm is", input$Q6,"patients")
  })
  
  output$responders_exp <- renderText({ 
    paste("3 / The number of responders in the Experimental arm is", input$Q7,"patients")
  })
  
  # posterior proba that response rate is greater than 20 pct, 30 pct etc
  
  #posterior_exp_dataset = reactive({data.frame(pe=pe, dens=marginal_posterior_density_pe() )
  # }) 
  
  # posterior_exp_dataset = data.frame(pe=pe, dens=marginal_posterior_density_pe() )
  
  proba_post_exp_greater_20pct <- reactive({
    posterior_exp_dataset = data.frame(pe=pe, dens=marginal_posterior_density_pe() )
    proba_post_exp_greater_20pct = sum(posterior_exp_dataset[which(posterior_exp_dataset$pe>0.2), 2])/sum(posterior_exp_dataset[,2])
    proba_post_exp_greater_20pct
  })        
  
  proba_post_exp_greater_30pct <- reactive({
    posterior_exp_dataset = data.frame(pe=pe, dens=marginal_posterior_density_pe() )
    proba_post_exp_greater_30pct = sum(posterior_exp_dataset[which(posterior_exp_dataset$pe>0.3), 2])/sum(posterior_exp_dataset[,2])
    proba_post_exp_greater_30pct
  })
  
  proba_post_exp_greater_40pct <- reactive({
    posterior_exp_dataset = data.frame(pe=pe, dens=marginal_posterior_density_pe() )
    proba_post_exp_greater_40pct = sum(posterior_exp_dataset[which(posterior_exp_dataset$pe>0.4), 2])/sum(posterior_exp_dataset[,2])
    proba_post_exp_greater_40pct
  })
  
  proba_post_exp_greater_50pct <- reactive({
    posterior_exp_dataset = data.frame(pe=pe, dens=marginal_posterior_density_pe() )
    proba_post_exp_greater_50pct = sum(posterior_exp_dataset[which(posterior_exp_dataset$pe>0.5), 2])/sum(posterior_exp_dataset[,2])
    proba_post_exp_greater_50pct
  })
  
  output$proba_post_exp_greater_20pct_text <- renderText({ 
    paste("1 / The posterior probability that the response rate is greater than 20% is", round(100*proba_post_exp_greater_20pct(),0),"%")
  })
  
  output$proba_post_exp_greater_30pct_text <- renderText({ 
    paste("2 / The posterior probability that the response rate is greater than 30% is", round(100*proba_post_exp_greater_30pct(),0),"%")
  })
  
  output$proba_post_exp_greater_40pct_text <- renderText({ 
    paste("3 / The posterior probability that the response rate is greater than 40% is", round(100*proba_post_exp_greater_40pct(),0),"%")
  })
  
  output$proba_post_exp_greater_50pct_text <- renderText({ 
    paste("4 / The posterior probability that the response rate is greater than 50% is", round(100*proba_post_exp_greater_50pct(),0),"%")
  })
  
}

shinyApp(ui = ui, server = server)

#setwd("S:/SLMS_CTU_Analysis/Hakim files/BAR-JDM/BAR JDM workshop")
#rsconnect::deployApp()
