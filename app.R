
#install.packages("shiny")
library(shiny)

#install.packages("rriskDistributions")
library(rriskDistributions)

#install.packages("zipfR")
library(zipfR)

#
to_prob_scale <- function(x)
  x/sum(x)

# Define UI for app ----
ui <- fluidPage(
  
  titlePanel("Blocks Game"),
  
  # Sidebar layout with inputs and outputs ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      tabsetPanel( 
        
        tabPanel("Prior",
                 br(),
                 # Input Q1  ----
                 sliderInput(inputId = "Q1",
                             label = "Q1: I think the proportion of black blocks is:",
                             min = 0,
                             max = 1,
                             value = 0.3),
                 
                 # Input Q2  ----
                 sliderInput(inputId = "Q2",
                             label = "Q2: I am 75% sure that the proportion of black blocks is less than:",
                             min = 0,
                             max = 1,
                             value = 0.5),
        ) ,
        
        tabPanel("Trial",
                 br(),
                 sliderInput(inputId = "Q5",
                             label = "Sample size:",
                             min = 5,
                             max = 5,
                             value = 5), 
                 
                 sliderInput(inputId = "Q6",
                             label = "Number of black blocks is:",
                             min = 0,
                             max = 5,
                             value = 0),

                 actionButton("button_add_data", "Add repeated sampling"),
                 checkboxInput("chk_box_all_data", "Use all data", FALSE),
                 
                 sliderInput(inputId = "subset_samples",
                             label = "Number of samples plotted",
                             min = 1,
                             max = 1,
                             step = 1,
                             value = 1),
                 
                 tableOutput("combinedTable"),
                 textOutput("totals_hypo_trial"),
                 actionButton("button_new_trial", "Start new trial")
        ) 
        
      )
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      tabsetPanel(
        
        tabPanel("Prior", plotOutput(outputId = "control_prior")
                 , textOutput("mode")
                 , textOutput("r_and_s")
                 , textOutput("beta_credble_intervals")
                 , textOutput("pseudo_obs")
        ),
        
        tabPanel("Posteriors",  plotOutput(outputId = "posterior_plot")
                 , textOutput("posterior_mode")
                 , textOutput("posterior_r_and_s")
                 , textOutput("posterior_beta_credble_intervals")
                 , textOutput("posterior_pseudo_obs")
        )
      ) 
      
    ) # end of main panel
  ) # end of sidebarLayout
) # end of UI

#
server <- function(input, output) {
  
  combined_df <- reactiveVal(NULL)
  
  observeEvent(input$button_add_data, {
    current <- combined_df()
    new_df <- rbind(current, c(`Sample size` = input$Q5,
                               `# black blocks` = input$Q6))
    combined_df(new_df)
  })
  
  output$combinedTable <- renderTable({
    combined_df()
  })
  
  observe({
    val <- nrow(combined_df())
    updateSliderInput(inputId = "subset_samples", max = val)
  })
  
  output$totals_hypo_trial <- renderText({
    if (is.null(combined_df())) {
      " "
    } else {
      glue::glue("Sample size used is {N_sample()} and number of black blocks used is {success_control()}")
    }
  })
  
  observeEvent(input$button_new_trial, {
    combined_df(NULL)
    
    updateSliderInput(inputId = "subset_samples", max = 1)
  })
  
  control_beta_a <- reactive({
    
    mode <- input$Q1
    p75 <- input$Q2
    
    min.SS <- function(params) {
      alpha <- params[1]
      beta <- params[2]
      
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
  
  ## Prior plot
  output$control_prior <- renderPlot({
    
    # define range
    p = seq(0,1, length=100)
    prior_c <- to_prob_scale(dbeta(p, control_beta_a(), control_beta_b()))
    cdf <- cumsum(prior_c)
    mode_c <- round(p[prior_c == max(prior_c)], 2)
    min75 <- round(min(p[cdf >= 0.75]), 2)
    
    # create plot of corresponding Beta distribution
    plot(p, prior_c,
         ylab='density', xlab="Proportion of black blocks",
         main='',
         type='l', lwd=1.5, col='red', xaxt='n')
    polygon(c(p[cdf >= 0.75], 1, min(p[cdf >= 0.75])), c(prior_c[cdf >= 0.75], 0, 0), col="lightblue", border=NA)
    abline(v = mode_c, lty = 2, col = "black")
    
    axis(1, at = sort(c(min75, mode_c, seq(0, 1, 0.2))), label=TRUE)
  })
  
  ## Probability greater than mode as reactive
  proba_greater_mode <- reactive({
    proba <- pbeta(input$Q1, shape1=control_beta_a(), shape2=control_beta_b(), lower.tail = FALSE)
    round(proba,2)*100
  })
  
  # output of the probability greater than mode
  output$proba_greater_mode <- renderText({ 
    paste("The parameters for Q1 and Q2 that you have selected correspond to a probability beyond the mode of",
          proba_greater_mode(), "%")
  })
  
  ## Probability lower p75 as reactive
  proba_lower_p75 <- reactive({
    proba_p75 <- pbeta(input$Q2, shape1=control_beta_a(), shape2=control_beta_b(), lower.tail = TRUE)
    round(proba_p75,2)*100
  })
  
  # output of the probability lower p75
  output$proba_lower_p75 <- renderText({ 
    glue::glue("- The blue area to the right of {input$Q2} is 25%. It means that you think that there is one chance out of 4 that the response rate will be greater than {input$Q2}.") 
  })
  
  output$alpha <- renderText({ 
    # paste("The mode and p75 that you have selected correspond to alpha of", round(control_beta_a(),2))
  })
  
  output$beta <- renderText({ 
    # paste("The mode and p75 that you have selected correspond to beta of", round(control_beta_b(),2))
  })
  
  output$all_priors <- renderPlot({
    x = seq(-10,1, length=1000)
    p = seq(0,1, length=100)
    
    par(mfrow=c(2,2))
    plot(p, to_prob_scale(dbeta(p, control_beta_a(), control_beta_b())), ylab='density',
         type ='l', lwd=1.5, xlab="Proportion of black blocks", col='red', main='')
    
    plot(p, to_prob_scale(dbeta(p, control_beta_a(), control_beta_b())), ylab='density',
         type ='l', lwd=1.5, xlab="Response rate", col='red', main='')
    lines(pe, to_prob_scale(marginal_prior_density()), lwd=1.5,col='purple')
  })
  
  ############################################################################## 
  ## trial and posteriors

  N_sample <- reactive({ 
    if (input$chk_box_all_data) {
      sum(combined_df()[1:input$subset_samples, "Sample size"])
    } else {
      input$Q5
    }
  })
  
  success_control <- reactive({ 
    if (input$chk_box_all_data) {
      sum(combined_df()[1:input$subset_samples, "# black blocks"])
    } else {
      input$Q6
    }
  })

  failures_control <- reactive({
    n_failures_control = N_sample() - success_control()
  })
  
  beta_posterior <- reactive({
    probs <- seq(0, 1, length=100)
    list(p = probs,
         vals = dbeta(probs, shape1 = control_beta_a() + success_control(),
                      shape2 = control_beta_b() + failures_control()))
  })
  
  posterior_a <- reactive({
    control_beta_a() + success_control()
  })
  
  posterior_b <- reactive({
    control_beta_b() + failures_control()
  })
  
  #
  output$posterior_plot <- renderPlot({
    plot(beta_posterior()$p, beta_posterior()$vals,
         ylab='density', xlab="Proportion of black blocks",
         type ='l', lwd=1.5, 
         col='red', main='')
  })
  
  ############################################################################## 
  ## summary stats page
  
  output$mode <- renderText({
    para_a <- control_beta_a()
    para_b <- control_beta_b()
    
    paste0("(i) The mode is ", round((para_a - 1)/(para_a + para_b - 2), 2), ". ",
    "The mean is ", round(para_a/(para_a + para_b), 2), ". ",
    "The standard deviation is ", round(sqrt((para_a*para_b)/((para_a + para_b)^2*(para_a+para_b+1))), 2))
  })
  
  output$r_and_s <- renderText({
    para_a <- control_beta_a()
    para_b <- control_beta_b()
    
    paste0("(ii) Values of r and s are ", round(para_a,2), " and ", round(para_b,2))
  })

  output$beta_credble_intervals <- renderText({
    para_a <- control_beta_a()
    para_b <- control_beta_b()
    
    paste0("(iii) The 50% credible interval is (",
           round(qbeta(0.25, para_a, para_b), 2), ", ",
           round(qbeta(0.75, para_a, para_b), 2), "). ",
           "The 95% credible interval is (",
           round(qbeta(0.025, para_a, para_b), 2), ", ",
           round(qbeta(0.975, para_a, para_b), 2), ")")
  })
  
  output$pseudo_obs <- renderText({
    para_a <- control_beta_a()
    para_b <- control_beta_b()
    
    paste0("(iv) The number of “pseudo-observations” to which the prior is equivalent is ",
           round(para_a + para_b, 2), " and the number of them which are black is ", round(para_a, 2))
  })

  output$proba_greater_mode <- renderText({ 
    paste("2 / The probability that a response rate greater than the mode will be observed is estimated at", proba_greater_mode(),"%")
  })
  
  output$positive_trt_effect <- renderText({ 
    paste("1 / The probability that there will be an improvement thanks to the experimental treatment strategy is estimated at", 100*input$Q3,"%")
  })
  
  # output of chance of an positive treatment effect greater than 25 pct
  output$positive_trt_effect_25pct <- renderText({ 
    paste("2 / The chance that this improvement will be greater than 25% is estimated at", 100*input$Q4,"%")
  })
  
  output$trial_size <- renderText({ 
    paste("1 / The trial has a sample size per arm of", input$Q5,"patients")
  })
  
  output$responders_ctrl <- renderText({ 
    paste("2 / The number of responders is", input$Q6,"patients")
  })
  
  output$posterior_mode <- renderText({
    para_a <- posterior_a()
    para_b <- posterior_b()
    
    paste0("(i) The posterior mode is ", round((para_a - 1)/(para_a + para_b - 2), 2), ". ",
           "The mean is ", round(para_a/(para_a + para_b), 2), ". ",
           "The standard deviation is ", round(sqrt((para_a*para_b)/((para_a + para_b)^2*(para_a+para_b+1))), 2))
  })

  output$posterior_r_and_s <- renderText({
    para_a <- posterior_a()
    para_b <- posterior_b()
    
    paste0("(ii) Values of r and s are ", round(para_a,2), " and ", round(para_b,2))
  })
  
  output$posterior_beta_credble_intervals <- renderText({
    para_a <- posterior_a()
    para_b <- posterior_b()
    
    paste0("(iii) The 50% credible interval is (",
           round(qbeta(0.25, para_a, para_b), 2), ", ",
           round(qbeta(0.75, para_a, para_b), 2), "). ",
           "The 95% credible interval is (",
           round(qbeta(0.025, para_a, para_b), 2), ", ",
           round(qbeta(0.975, para_a, para_b), 2), ")")
  })
  
  output$posterior_pseudo_obs <- renderText({
    para_a <- posterior_a()
    para_b <- posterior_b()
    
    paste0("(iv) The number of “pseudo-observations” to which the prior is equivalent is ",
           round(para_a + para_b, 2), " and the number of them which are black is ", round(para_a, 2))
  })
  
  # control arm
  ##TODO: update for control arm...
  
  proba_post_control_greater_20pct <- reactive({
    posterior_exp_dataset = data.frame(pe=pe, dens=marginal_posterior_density_pe() )
    sum(posterior_exp_dataset[which(posterior_exp_dataset$pe>0.2), 2])/sum(posterior_exp_dataset[,2])
  })        
  
  proba_post_control_greater_30pct <- reactive({
    posterior_exp_dataset = data.frame(pe=pe, dens=marginal_posterior_density_pe() )
    sum(posterior_exp_dataset[which(posterior_exp_dataset$pe>0.3), 2])/sum(posterior_exp_dataset[,2])
  })
  
  proba_post_control_greater_40pct <- reactive({
    posterior_exp_dataset = data.frame(pe=pe, dens=marginal_posterior_density_pe() )
    sum(posterior_exp_dataset[which(posterior_exp_dataset$pe>0.4), 2])/sum(posterior_exp_dataset[,2])
  })
  
  proba_post_control_greater_50pct <- reactive({
    posterior_exp_dataset = data.frame(pe=pe, dens=marginal_posterior_density_pe() )
    sum(posterior_exp_dataset[which(posterior_exp_dataset$pe>0.5), 2])/sum(posterior_exp_dataset[,2])
  })
  
  output$proba_post_control_greater_20pct_text <- renderText({ 
    paste("1 / The posterior probability that the response rate is greater than 20% is", round(100*proba_post_control_greater_20pct(),0),"%")
  })
  
  output$proba_post_control_greater_30pct_text <- renderText({ 
    paste("2 / The posterior probability that the response rate is greater than 30% is", round(100*proba_post_control_greater_30pct(),0),"%")
  })
  
  output$proba_post_control_greater_40pct_text <- renderText({ 
    paste("3 / The posterior probability that the response rate is greater than 40% is", round(100*proba_post_control_greater_40pct(),0),"%")
  })
  
  output$proba_post_control_greater_50pct_text <- renderText({ 
    paste("4 / The posterior probability that the response rate is greater than 50% is", round(100*proba_post_control_greater_50pct(),0),"%")
  })
}

shinyApp(ui = ui, server = server)

#rsconnect::deployApp()
