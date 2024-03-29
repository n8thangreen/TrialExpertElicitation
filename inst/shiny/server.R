# require(TrialExpertElicitation)

shinyServer(function(input, output){
  
  posterior20 <- reactive({
    length(input$hypo_data_size) > 0 && input$hypo_data_size == 20
  })
  
  posterior40 <- reactive({
    length(input$hypo_data_size) > 0 && input$hypo_data_size == 40
  })
  
  output$valid_range_q1 <- #reactive({
    renderText({
      paste("Min:", input$pc_q2 + 0.05, ";",
            "Max:", 0.95)
    })
  # })
  
  output$valid_range_q2 <- #reactive({
    renderText({
      paste("Min:", round(beta_percentile25(a = 0.5, input$pc_q1),2), ";",
            "Max:", round(beta_percentile25(a = 50, input$pc_q1),2))
      })
  # })
  
  ## Need to determine parameters of the prior distribution
  ## Reactive expression to generate the requested distribution. This is called
  ## whenever the inputs change. The renderers defined below then all use the value 
  ## computed from this expression
  priorParam <- reactive({

    ##TODO: reset sliders if values not valid
    # no_root_in_range <- !check_interval_valid_for_a(input$pc_q1, input$pc_q2)
    # 
    # if (no_root_in_range) {
    #   q1 <- 0.7
    #   q2 <- 0.5
    # 
    #   updateSliderInput(inputId = 'pc_q1', value = q1)
    #   updateSliderInput(inputId = 'pc_q2', value = q2)
    # } else {
      q1 <- input$pc_q1
      q2 <- input$pc_q2
    # }
    
    prior_summaries(q1, q2, input$theta_q1/100, input$theta_q2/100, input$expert)	
  })
  
  postParam <- reactive({
    if (!posterior40() && !posterior20())
      return(NULL)
    
    posterior_summaries(as.numeric(input$hypo_data_size) - input$n_cyc, input$mmf_succ,
                        input$n_cyc, input$cyc_succ, priorParam(), posterior40())
  })
  
  postSumary <- reactive({
    if((posterior40() | posterior20()) & input$postsum){
      n <- as.numeric(input$hypo_data_size)
      dtacase <- scan(glue::glue("data/data_scenario{n}.txt"))
      scen1 = posterior_summaries(n - dtacase[1], dtacase[3], dtacase[1], dtacase[2], priorParam(), posterior40())	
      scen2 = posterior_summaries(n - dtacase[4], dtacase[6], dtacase[4], dtacase[5], priorParam(), posterior40())
      scen3 = posterior_summaries(n - dtacase[7], dtacase[9], dtacase[7], dtacase[8], priorParam(), posterior40())
      
      return(data.frame(scen1, scen2, scen3))
    }else{
      NULL
    }
  })
  
  ## generating data to plot the prior and posterior densities of pC, pE and theta. 
  pC_priorDens <- reactive({
    ## Return the data to plot the prior and posterior densities of theta.
    x = vector(mode="numeric", length = 20)
    dist_plot_data(0,0,0,0, x, priorParam(), "pC", TRUE)
  })
  
  pC_postDens <- reactive({
    z = priorParam()
    x = postParam()
    n <- as.numeric(input$hypo_data_size)
    dist_plot_data(n_mmf = n - input$n_cyc,
             mmf_succ = input$mmf_succ,
             n_cyc = input$n_cyc,
             cyc_succ = input$cyc_succ,
             postParm = x,
             priorParm = z,
             parmInd = as.character("pC"),
             postind = FALSE)
  })
  
  pE_priorDens <- reactive({
    x = vector(mode="numeric", length =20)
    dist_plot_data(0,0,0,0, x, priorParam(), "pE", TRUE)
  })
  
  pE_postDens <- reactive({
    z = priorParam()
    x = postParam()
    n <- as.numeric(input$hypo_data_size)
    dist_plot_data(n - input$n_cyc, input$mmf_succ, input$n_cyc, input$cyc_succ, x, z, "pE", FALSE)
  })
  
  theta_priorDens <- reactive({
    x = vector(mode="numeric", length =20)
    dist_plot_data(0,0,0,0,x, priorParam(), "theta", TRUE)
  })
  
  theta_postDens <- reactive({
    z = priorParam()
    x = postParam() 
    n <- as.numeric(input$hypo_data_size)
    dist_plot_data(n - input$n_cyc, input$mmf_succ, input$n_cyc, input$cyc_succ, x, z, "theta", FALSE)
  })
  
  ###########
  ## plots
  
  ## Outputting prior and posterior densities of pC, pE and theta.
  output$pc_density <- renderPlot({
    ## Plot the prior and posterior distributions for pC
    if((posterior40() | posterior20()) & !input$postsum){
      y = pC_priorDens()
      x = pC_postDens()
      
      ns <- as.numeric(input$hypo_data_size)
      leg1 = c(bquote("Posterior" ~ n[C] ==  .(input$n_cyc) ~ S[C]  == .(input$cyc_succ) ~ n[M] == .(ns - input$n_cyc) ~ S[M] ==  .(input$mmf_succ)), "Prior")
      
      postPlot = plot(x$gridc, x$postd, type="l", lty=2, lwd=3, col="dark red",
                      main = "Prior and posterior densities of 6-month remission rate on control arm",
                      xlab = "control arm 6-month remission rate", ylab="Density",
                      xlim = c(0,1), ylim = range(c(x$postd, y$dens)), cex.lab = 1.5, cex.axis=1.5, cex.main = 2) 
      lines(y$gridc, y$dens, lty=1, lwd=3, col="red")
      legend("topleft", as.expression(leg1), col = c("dark red","red"), lty = c(2,1), lwd = c(3,3), cex=1.2, bty="n")	
      
    }else if((posterior40() | posterior20()) & input$postsum){
      y = pC_priorDens()
      z = priorParam()
      w = postSumary()
      n <- as.numeric(input$hypo_data_size)
      dtacase = scan(glue::glue("data/data_scenario{n}.txt"))
      u1 = dist_plot_data(n - dtacase[1], dtacase[3], dtacase[1], dtacase[2], w$scen1, z, "pC", FALSE)
      u2 = dist_plot_data(n - dtacase[4], dtacase[6], dtacase[4], dtacase[5], w$scen2, z, "pC", FALSE)
      u3 = dist_plot_data(n - dtacase[7], dtacase[9], dtacase[7], dtacase[8], w$scen3, z, "pC", FALSE)
      
      postPlot = plot(u1$gridc, u1$postd, type="l", lty=2, lwd=3, col="dark red",
                      main="Prior and posterior densities of 6-month remission rate on control arm/steroids",
                      xlab = "control arm 6-month remission rate", ylab="Density",
                      xlim =c(0,1), ylim = range(c(u1$postd, u2$postd, u3$postd, y$dens)), cex.lab = 1.5, cex.axis=1.5, cex.main = 2)  
      lines(u2$gridc, u2$postd, lty=3, lwd=3, col="dark red")
      lines(u3$gridc, u3$postd, lty=4, lwd=3, col="dark red")
      lines(y$gridc, y$dens, lty=1, lwd=3, col="red")
      
      leg1 = paste("Control arm Succ = ", dtacase[2], "; Experimental arm Succ = ", dtacase[3], sep="")
      leg2 = paste("Control arm Succ = ", dtacase[5], "; Experimental arm Succ = ", dtacase[6], sep="")
      leg3 = paste("Control arm Succ = ", dtacase[8], "; Experimental arm Succ = ", dtacase[9], sep="")
      
      legend("topleft", c(leg1,leg2, leg3,"Prior"),
             col = c("dark red","dark red","dark red","red"), lty = c(2,3,4,1), lwd = c(3,3,3,3), cex=1.2, bty="n")
      
    }else{
      x = pC_priorDens()
      ## Only want to plot the prior density in this scenario
      plot(x$gridc, x$dens, type="l", lty=1, lwd=3, col="red",
           main = "Prior density of 6-month remission rate on control arm/steroids",
           xlab = "Control arm 6-month remission rate", ylab="Density", xlim=c(0,1), cex.lab = 1.5, cex.axis=1.5, cex.main = 2)	
    }
  })
  
  output$pe_density <- renderPlot({
    ## Plot the prior and posterior distributions for pE
    if((posterior40() | posterior20()) & !input$postsum){
      y = pE_priorDens()
      x = pE_postDens()
      ns <- as.numeric(input$hypo_data_size)
      leg1 = c(bquote("Posterior" ~ n[C] ==  .(input$n_cyc) ~ S[C]  == .(input$cyc_succ) ~ n[M] == .(ns - input$n_cyc) ~ S[M] ==  .(input$mmf_succ)), "Prior")
      
      postPlot = plot(x$gride, x$postd, type="l", lty=2, lwd=3, col="dark green",
                      main="Prior and posterior densities of 6-month remission rate on Experimental",
                      xlab = "Experimental arm 6-month remission rate", ylab="Density",
                      xlim =c(0,1), ylim = range(c(x$postd, y$dens)), cex.lab = 1.5, cex.axis=1.5, cex.main = 2) 
      lines(y$gride, y$dens, lty=1, lwd=3, col="green")
      legend("topleft", as.expression(leg1), col = c("dark green","green"), lty = c(2,1), lwd = c(3,3), cex=1.2, bty="n")	
    }else if((posterior40() | posterior20()) & input$postsum){
      y = pE_priorDens()
      z = priorParam()
      w = postSumary()
      n <- as.numeric(input$hypo_data_size)
      dtacase = scan(glue::glue("data/data_scenario{n}.txt"))
      u1 = dist_plot_data(n - dtacase[1], dtacase[3], dtacase[1], dtacase[2], w$scen1, z, "pE", FALSE)
      u2 = dist_plot_data(n - dtacase[4], dtacase[6], dtacase[4], dtacase[5], w$scen2, z, "pE", FALSE)
      u3 = dist_plot_data(n - dtacase[7], dtacase[9], dtacase[7], dtacase[8], w$scen3, z, "pE", FALSE)
      
      postPlot = plot(u1$gride, u1$postd, type="l", lty=2, lwd=3, col="dark green",
                      main="Prior and posterior densities of 6-month remission rate on Experimental",
                      xlab = "Experimental 6-month remission rate", ylab="Density",
                      xlim =c(0,1), ylim = range(c(u1$postd, u2$postd, u3$postd, y$dens)),
                      cex.lab = 1.5, cex.axis=1.5, cex.main = 2)  
      lines(u2$gride, u2$postd, lty=3, lwd=3, col="dark green")
      lines(u3$gride, u3$postd, lty=4, lwd=3, col="dark green")
      lines(y$gride, y$dens, lty=1, lwd=3, col="green")
      
      leg1 = paste("Control arm Succ = ", dtacase[2], "; Experimental arm Succ = ", dtacase[3], sep="")
      leg2 = paste("Control arm Succ = ", dtacase[5], "; Experimental arm Succ = ", dtacase[6], sep="")
      leg3 = paste("Control arm Succ = ", dtacase[8], "; Experimental arm Succ = ", dtacase[9], sep="")
      
      legend("topleft", c(leg1, leg2, leg3, "Prior"),
             col = c("dark green","dark green","dark green","green"), lty = c(2,3,4,1), lwd = c(3,3,3,3), cex=1.2, bty="n")
    }else{
      x = pE_priorDens()
      ## Only want to plot the prior density in this scenario
      priorPlot = plot(x$gride, x$dens, type="l", lty=1, lwd=3, col="green",
                       main = "Prior density of 6-month remission rate on Experimental arm",
                       xlab = "Experimental 6-month remission rate", ylab="Density",
                       xlim =c(0,1), cex.lab = 1.5, cex.axis=1.5, cex.main = 2)	
    }
  })
  
  output$pepc_density <- renderPlot({		
    ## Plot the prior and posterior distributions for pE
    if((posterior40() | posterior20()) & !input$postsum){
      y = pC_postDens()
      x = pE_postDens()
      postPlot = plot(y$gridc, y$postd, type="l", lty=2, lwd=3, col="dark red",
                      main="Posterior densities of 6-month remission rates on Experimental and control arm",
                      xlab = "6-month remission rate", ylab="Density",
                      xlim =c(0,1), ylim = range(c(y$postd, x$postd)), cex.lab = 1.5, cex.axis=1.5, cex.main = 2) 
      lines(x$gride, x$postd, lty=2, lwd=3, col="dark green")
      legend("topleft", c("Control arm/steroids", "Experimental arm/steroids"), col = c("dark red","dark green"), lty = c(2,2), lwd = c(3,3), cex=1.2, bty="n")	
    }else if((posterior40() | posterior20()) & input$postsum){
      z = priorParam()
      w = postSumary()
      ns <- as.numeric(input$hypo_data_size)
      dtacase = scan(glue::glue("data/data_scenario{ns}.txt"))
      u1pe = dist_plot_data(ns - dtacase[1], dtacase[3], dtacase[1], dtacase[2], w$scen1, z, "pE", FALSE)
      u2pe = dist_plot_data(ns - dtacase[4], dtacase[6], dtacase[4], dtacase[5], w$scen2, z, "pE", FALSE)
      u3pe = dist_plot_data(ns - dtacase[7], dtacase[9], dtacase[7], dtacase[8], w$scen3, z, "pE", FALSE)
      
      u1pc = dist_plot_data(ns - dtacase[1], dtacase[3], dtacase[1], dtacase[2], w$scen1, z, "pC", FALSE)
      u2pc = dist_plot_data(ns - dtacase[4], dtacase[6], dtacase[4], dtacase[5], w$scen2, z, "pC", FALSE)
      u3pc = dist_plot_data(ns - dtacase[7], dtacase[9], dtacase[7], dtacase[8], w$scen3, z, "pC", FALSE)
      
      leg1 = bquote("Control arm:" ~ n[C] == .(dtacase[1]) ~ S[C]  == .(dtacase[2]) ~ n[M] == .(ns - dtacase[1]) ~ S[M] ==  .(dtacase[3]))
      leg2 = bquote("Control arm:" ~ n[C] == .(dtacase[4]) ~ S[C]  == .(dtacase[5]) ~ n[M] == .(ns - dtacase[4]) ~ S[M] ==  .(dtacase[6]))
      leg3 = bquote("Control arm:" ~ n[C] == .(dtacase[7]) ~ S[C]  == .(dtacase[8]) ~ n[M] == .(ns - dtacase[7]) ~ S[M] ==  .(dtacase[9]))
      leg4 = bquote("Experimental arm:" ~ S[C]  == .(dtacase[2]) ~ S[M] ==  .(dtacase[3]))
      leg5 = bquote("Experimental arm:" ~ S[C]  == .(dtacase[5]) ~ S[M] ==  .(dtacase[6]))
      leg6 = bquote("Experimental arm:" ~ S[C]  == .(dtacase[8]) ~ S[M] ==  .(dtacase[9]))
      
      postPlot = plot(u1pc$gridc, u1pc$postd, type="l", lty=2, lwd=3, col="dark red",
                      main="Posterior densities of 6-month remission rates on Experimental and control arm",
                      xlab = "6-month remission rate", ylab="Density",
                      xlim =c(0,1), ylim = range(c(u1pc$postd, u2pc$postd, u3pc$postd, u1pe$postd, u2pe$postd, u3pe$postd)),
                      cex.lab = 1.5, cex.axis=1.5, cex.main = 2)  
      lines(u2pc$gridc, u2pc$postd, lty=3, lwd=3, col="dark red")
      lines(u3pc$gridc, u3pc$postd, lty=4, lwd=3, col="dark red")
      lines(u1pe$gride, u1pe$postd, lty=2, lwd=3, col="dark green")
      lines(u2pe$gride, u2pe$postd, lty=3, lwd=3, col="dark green")
      lines(u3pe$gride, u3pe$postd, lty=4, lwd=3, col="dark green")
      
      legend("topleft", as.expression(c(leg1, leg2, leg3, leg4, leg5, leg6)),
             col = c("dark red","dark red","dark red", "dark green","dark green","dark green"),
             lty = c(2,3,4,2,3,4), lwd = c(3,3,3,3,3,3), cex=1.2, bty="n")
    }else{
      y = pC_priorDens()
      x = pE_priorDens()
      ## Only want to plot the prior density in this scenario
      priorPlot = plot(y$gridc, y$dens, type="l", lty=1, lwd=3, col="red",
                       main = "Prior density of 6-month remission rates on Experimental and control arm",
                       xlab = "6-month remission rate", ylab="Density",
                       xlim =c(0,1), ylim = range(c(y$dens, x$dens)), cex.lab = 1.5, cex.axis=1.5, cex.main = 2)	
      lines(x$gride, x$dens, lty=1, lwd=3, col="green")
      legend("topleft", c("Control arm/steroids", "Experimental arm/steroids"), col = c("red","green"), lty = c(1,1), lwd = c(3,3), cex=1.2, bty="n")	
    }
  })
  
  output$theta_density <- renderPlot({	
    ## Plot the prior and posterior distributions for theta
    
    pc_to_pe_transformation <- function(pc, lor)
      1.0/(1 + exp(-(lor + log(pc/(1 - pc)))))
    
    pe_lim1 = pc_to_pe_transformation(input$pc_q1, -2)
    pe_lim2 = pc_to_pe_transformation(input$pc_q1, -1)
    pe_lim3 = pc_to_pe_transformation(input$pc_q1, 0)
    pe_lim4 = pc_to_pe_transformation(input$pc_q1, 1)
    pe_lim5 = pc_to_pe_transformation(input$pc_q1, 2)
    
    label1 <- paste("6-month remission rate on Experimental arm assuming that on control arm is ", round(input$pc_q1, 2), sep="")
    
    if((posterior40() | posterior20()) & !input$postsum){
      y = theta_priorDens()
      x = theta_postDens()
      
      ns <- as.numeric(input$hypo_data_size)
      leg1 = c(bquote("Posterior" ~ n[C] ==  .(input$n_cyc) ~ S[C]  == .(input$cyc_succ) ~ n[M] == .(ns - input$n_cyc) ~ S[M] ==  .(input$mmf_succ)), "Prior")
      
      par(pty="m", plt=c(0.1, 1, 0, 1), omd=c(0.1,1.0,0,0.9))
      par(mfrow = c(2, 1))
      
      postPlot = plot(x$gridt, x$postd, type="l", lty=2, lwd=3, col="dark blue",
                      main="", ylab="Density", xaxt = "n", cex.lab = 1.5, cex.axis=1.5, ylim = range(c(x$postd, y$dens))) 
      lines(y$gridt, y$dens, lty=1, lwd=3, col="blue")
      legend("topleft", as.expression(leg1), col = c("dark blue","blue"), lty = c(2,1), lwd = c(3,3), cex=1.2, bty="n")	
      mtext(side=3, "Prior and posterior densities of the log-odds ratio", line=1.5, cex = 2)
      
      ## Setting up axis 1
      axis(side=1, at=c(-2, -1, 0, 1, 2), labels=c(-2, -1, 0, 1, 2), cex.axis=1.5)
      linloc <- par()$usr[3]
      abline(h=linloc, col="black")
      mtext(side=1, "log-odds ratio", line=2.5, cex=1.5, col="black")
      
      ## axis 2
      xaxis2 = seq(-2, 2, by=1)
      par(plt=c(0.1,1,0.6,1))
      plot(xaxis2, type="n",xaxt="n", xlab="", yaxt="n", ylab="", xlim = c(-2, 2), bty="n")
      axis(side=1, at = c(-2, -1, 0, 1, 2), labels = round(c(pe_lim1, pe_lim2, pe_lim3, pe_lim4, pe_lim5), 2), col="red", col.axis="red", cex.axis=1.5)
      linloc <- par()$usr[3]
      abline(h=linloc, col="red")
      mtext(side=1, label1, line=3.0, cex=1.5, col="red")	
    }else if((posterior40() | posterior20()) & input$postsum){
      y = theta_priorDens()
      z = priorParam()
      w = postSumary()
      
      ns <- as.numeric(input$hypo_data_size)
      dtacase = scan(glue::glue("data/data_scenario{n}.txt"))
      u1 = dist_plot_data(n - dtacase[1], dtacase[3], dtacase[1], dtacase[2], w$scen1, z, "theta", FALSE)
      u2 = dist_plot_data(n - dtacase[4], dtacase[6], dtacase[4], dtacase[5], w$scen2, z, "theta", FALSE)
      u3 = dist_plot_data(n - dtacase[7], dtacase[9], dtacase[7], dtacase[8], w$scen3, z, "theta", FALSE)
      
      par(pty="m", plt=c(0.1, 1, 0, 1), omd=c(0.1,1.0,0,0.9))
      par(mfrow = c(2, 1))
      
      postPlot = plot(u1$gridt, u1$postd, type="l", lty=2, lwd=3, col="dark blue",
                      main="", ylab="Density", xaxt="n", cex.lab = 1.5, cex.axis=1.5, ylim = range(c(u1$postd, u2$postd, u3$postd, y$dens)))
      lines(u2$gridt, u2$postd, lty=3, lwd=3, col="dark blue")
      lines(u3$gridt, u3$postd, lty=4, lwd=3, col="dark blue")
      lines(y$gridt, y$dens, lty=1, lwd=3, col="blue")
      
      leg1 = paste("Control arm Succ = ", dtacase[2], "; Experimental arm Succ = ", dtacase[3], sep="")
      leg2 = paste("Control arm Succ = ", dtacase[5], "; Experimental arm Succ = ", dtacase[6], sep="")
      leg3 = paste("Control arm Succ = ", dtacase[8], "; Experimental arm Succ = ", dtacase[9], sep="")
      
      legend("topright", c(leg1, leg2, leg3, "Prior"),
             col = c("dark blue","dark blue","dark blue","blue"), lty = c(2,3,4,1), lwd = c(3,3,3,3), cex=1.2, bty="n")	
      mtext(side=3, "Prior and posterior densities of the log-odds ratio", line=1.5, cex = 2)
      
      ## Setting up axis 1
      axis(side=1, at=c(-2, -1, 0, 1, 2), labels=c(-2, -1, 0, 1, 2), cex.axis=1.5)
      linloc <- par()$usr[3]
      abline(h=linloc, col="black")
      mtext(side=1, "log-odds ratio", line=2.5, cex=1.5, col="black")
      
      ## axis 2
      xaxis2 = seq(-2, 2, by=1)
      par(plt=c(0.1,1,0.6,1))
      plot(xaxis2, type="n",xaxt="n", xlab="", yaxt="n", ylab="", xlim = c(-2, 2), bty="n")
      axis(side=1, at = c(-2, -1, 0, 1, 2), labels = round(c(pe_lim1, pe_lim2, pe_lim3, pe_lim4, pe_lim5), 2), col="red", col.axis="red", cex.axis=1.5)
      linloc <- par()$usr[3]
      abline(h=linloc, col="red")
      mtext(side=1, label1, line=3.0, cex=1.5, col="red")			
    }else{
      x = theta_priorDens()
      ## Only want to plot the prior density in this scenario
      par(pty="m", plt=c(0.1, 1, 0, 1), omd=c(0.1,1.0,0,0.9))
      par(mfrow = c(3, 1))
      
      priorPlot = plot(x$gridt, x$dens, type="l", lty=1, lwd=3, col="blue", main="", ylab="Density", xaxt="n", cex.lab = 1.5, cex.axis=1.5)
      mtext(side=3, "Prior density of the log-odds ratio", line=1.5, cex = 2)
      
      ## Setting up axis 1
      axis(side=1, at=c(-2, -1, 0, 1, 2), labels=c(-2, -1, 0, 1, 2), cex.axis=1.5)
      linloc <- par()$usr[3]
      abline(h=linloc, col="black")
      mtext(side=1, "log-odds ratio", line=2.5, cex=1.5, col="black")
      
      ## axis 2
      xaxis2 = seq(-2, 2, by=1)
      par(plt=c(0.1,1,0.6,1))
      plot(xaxis2, type="n",xaxt="n", xlab="", yaxt="n", ylab="", xlim = c(-2, 2), bty="n")
      axis(side=1, at = c(-2, -1, 0, 1, 2),
           labels = round(c(pe_lim1, pe_lim2, pe_lim3, pe_lim4, pe_lim5), 2),
           col="red", col.axis="red", cex.axis=1.5)
      linloc <- par()$usr[3]
      abline(h=linloc, col="red")
      mtext(side=1, label1, line=3.0, cex=1.5, col="red")	
      
      ## axis 3
      xaxis2 = seq(-2, 2, by=1)
      par(plt=c(0.1,1,0.6,1))
      plot(xaxis2, type="n",xaxt="n", xlab="", yaxt="n", ylab="", xlim = c(-2, 2), bty="n")
      axis(side=1, at = c(-2, -1, 0, 1, 2),
           labels = round(c(pe_lim1, pe_lim2, pe_lim3, pe_lim4, pe_lim5), 2),
           col="red", col.axis="red", cex.axis=1.5)
      linloc <- par()$usr[3]
      abline(h=linloc, col="red")
      mtext(side=1, label1, line=3.0, cex=1.5, col="red")	
    }
  })
  
  ## Generate a summary of the prior or posterior densities as necessary
  output$summary = renderUI({
    if (posterior40() | posterior20()){
      
      x = postParam()
      ns <- as.numeric(input$hypo_data_size)
      
      summary_text <- paste(
        glue::glue("Posterior distributions update the priors implied by the opinions:<br/> Q1 = {input$pc_q1}, Q2 = {input$pc_q2}, Q3 = {input$theta_q1}%, Q4 = {input$theta_q2}%."),
         glue::glue("<br/><br/>Posterior distributions incorporate the following data:"),
      glue::glue("<ul><li>{input$n_cyc} patients randomised to control arm of whom {input$cyc_succ} were in remission within six months</li>"),
      glue::glue("<li>{ns - input$n_cyc} patients randomised to Experimental of whom {input$mmf_succ} were in remission within six months</li></ul>"),
      glue::glue("<br/><br/>Summary of the posterior distribution of the 6-month remission rate on control arm/steroids:"),
      glue::glue("<ul><li>90% credibility interval for control arm remission rate is ({round(x[4], 2)}, {round(x[5], 2)})</li>"),
      glue::glue("<li>Your distribution has mode {round(x[2], 2)} and mean {round(x[1], 2)} and standard deviation {round(x[3], 2)}</li></ul>"),
      glue::glue("<br/><br/>Summary of the posterior distribution of the 6-month remission rate on Experimental/steroids:"),
      glue::glue("<ul><li>90% credibility interval for Experimental remission rate is ({round(x[9], 2)}, {round(x[10], 2)})</li>"),
      glue::glue("<li>Your distribution has mode {round(x[7], 2)} and mean {round(x[6], 2)} and standard deviation {round(x[8], 2)}</li></ul>"),
      glue::glue("<br/><br/>Summary of the posterior distribution of the log-odds ratio:"),
      glue::glue("<ul><li>90% credibility interval for the log-odds ratio is ({round(x[14], 2)}, {round(x[15], 2)})</li>"),
      glue::glue("<li>Your distribution for the log-odds ratio has mode {round(x[12], 2)} and mean {round(x[11], 2)} and standard deviation {round(x[13], 2)}</li>"),
      glue::glue("<li>Posterior probability that the 6-month remission rate on Experimental exceeds that on control arm is {round(x[16], 2)}</li>"),
      glue::glue("<li>Posterior probability that the 6-month remission rate on Experimental is at worst 10% less than that on control arm is {round(x[17], 2)}</li>"),
      glue::glue("<li>Posterior probability that the 6-month remission rate on control arm exceeds that on Experimental by more than 10% is {round(1 - x[17], 2)}</li></ul>"),
      sep = '')
      
      HTML(summary_text)
    }else{
      x <- priorParam()
      
      pcquantile0_05 = qbeta(0.05, shape1 = x[1], shape2 = x[2], lower.tail=TRUE)
      pcquantile0_95 = qbeta(0.95, shape1 = x[1], shape2 = x[2], lower.tail=TRUE)
      pc_mode = (x[1] - 1)/(x[1] + x[2] - 2)
      pc_mean = x[1]/(x[1] + x[2])
      pc_var = (x[1]*x[2])/(((x[1] + x[2])^2)*(x[1] + x[2] + 1))
      
      tquantile0_05 = qnorm(0.05, mean=x[9], sd = sqrt(x[10]), lower.tail=TRUE)
      tquantile0_95 = qnorm(0.95, mean=x[9], sd = sqrt(x[10]), lower.tail=TRUE)
      
      summary_text <- paste(
        glue::glue("Prior distributions are determined by the opinions: <br/> Q1 = {input$pc_q1} <br/> Q2 = {input$pc_q2} <br/> Q3 = {input$theta_q1}% <br/> Q4 = {input$theta_q2}%"),
        "<br/><br/>Summary of the prior distribution of the 6-month remission rate on control arm/steroids:",
        glue::glue("<ul><li>90% credibility interval for control arm remission rate is ({round(pcquantile0_05, 2)}, {round(pcquantile0_95, 2)})"),
        glue::glue("<li>Your distribution is Beta with parameters {round(x[1], 2)} and {round(x[2], 2)}</li>"),
        glue::glue("<li>with mode {round(pc_mode, 3)} and mean {round(pc_mean, 2)} and standard deviation {round(sqrt(pc_var), 2)}</li>"),
        glue::glue("<li>This information is equivalent to {round(x[3], 0)} observations on control arm</li></ul>"),
        "<br/><br/>Summary of the prior distribution of Experimental/steroid 6-month remission rate:",
        glue::glue("<ul><li>90% credibility interval for Experimental remission rate is ({round(x[7], 2)}, {round(x[8], 2)})</li>"),
        glue::glue("<li>Your distribution has mode {round(x[5], 2)} and mean {round(x[4], 2)} and standard deviation {round(x[6], 2)}</li>"),
        glue::glue("<li>Proportion which we are 75% confident the true 6-month remission rate on Experimental exceeds is {round(x[14], 2)}</li></ul>"),
        "<br/><br/>Summary of the prior distribution of log-odds ratio:",
        glue::glue("<ul><li>90% credibility interval for log-odds ratio is ({round(tquantile0_05, 2)}, {round(tquantile0_95, 2)})</li>"),
        glue::glue("<li>Your distribution is Gaussian with mode {round(x[9], 2)} and mean {round(x[9], 2)} and standard deviation {round(sqrt(x[10]), 2)}</li>"),
        glue::glue("<li>This information is equivalent to what would be generated by a RCT comparing control arm with Experimental recruiting {round(0.5*x[13], 0)} patients per treatment arm</li>"),
        glue::glue("<li>Prior probability that the 6-month remission rate on Experimental exceeds that on control arm is {round(x[11], 2)}</li>"),
        glue::glue("<li>Prior probability that the 6-month remission rate on Experimental is at worst 10% less than that on control arm is {round(1 - x[12], 2)}</li>"),
        glue::glue("<li>Prior probability that the 6-month remission rate on control arm exceeds that on Experimental arm by more than 10% is {round(x[12], 2)}</li></ul>"),
        sep = '')
      
      HTML(summary_text)
    }		
  })
})
