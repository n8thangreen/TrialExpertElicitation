
shinyServer(function(input, output){
  
  posterior20 <- reactive({
    length(input$hypo_data_size) > 0 && input$hypo_data_size == 20
  })
  
  posterior40 <- reactive({
    length(input$hypo_data_size) > 0 && input$hypo_data_size == 40
  })
  
  ## Need to determine parameters of the prior distribution
  ## Reactive expression to generate the requested distribution. This is called
  ## whenever the inputs change. The renderers defined below then all use the value 
  ## computed from this expression
  priorParam <- reactive({
    ## prior returns a vector containing the summaries of the posterior distributions of pC, pE and theta
    priorcall(input$pc_q1, input$pc_q2, input$theta_q1, input$theta_q2, input$expert)	
  })
  
  postParam <- reactive({
    if (!posterior40() && !posterior20())
      return(NULL)
    
    postSumry(as.numeric(input$hypo_data_size) - input$n_cyc, input$mmf_succ, input$n_cyc, input$cyc_succ, priorParam(), posterior40())
  })
  
  postSumary <- reactive({
    if((posterior40() | posterior20()) & input$postsum){
      n <- as.numeric(input$hypo_data_size)
      dtacase <- scan(glue::glue("data_scenario{n}.txt"))
      scen1 = postSumry(n - dtacase[1], dtacase[3], dtacase[1], dtacase[2], priorParam(), posterior40())	
      scen2 = postSumry(n - dtacase[4], dtacase[6], dtacase[4], dtacase[5], priorParam(), posterior40())
      scen3 = postSumry(n - dtacase[7], dtacase[9], dtacase[7], dtacase[8], priorParam(), posterior40())
      
      return(data.frame(scen1, scen2, scen3))	
    }else{
      NULL
    }
  })
  
  ## generating data to plot the prior and posterior densities of pC, pE and theta. 
  pC_priorDens <- reactive({
    ## Return the data to plot the prior and posterior densities of theta.
    z = priorParam()
    x = vector(mode="numeric", length =20)
    distPlot(0,0,0,0, x, z, as.character("pC"), TRUE)
  })
  
  pC_postDens <- reactive({
    z = priorParam()
    x = postParam()
    n <- as.numeric(input$hypo_data_size)
    distPlot(n - input$n_cyc, input$mmf_succ, input$n_cyc, input$cyc_succ, x, z, as.character("pC"), FALSE)
  })
  
  pE_priorDens <- reactive({
    z = priorParam()
    x = vector(mode="numeric", length =20)
    distPlot(0,0,0,0, x, z, as.character("pE"), TRUE)
  })
  pE_postDens <- reactive({
    z = priorParam()
    x = postParam()
    n <- as.numeric(input$hypo_data_size)
    distPlot(n - input$n_cyc, input$mmf_succ, input$n_cyc, input$cyc_succ, x, z, as.character("pE"), FALSE)
  })
  
  theta_priorDens <- reactive({
    z = priorParam()
    x = vector(mode="numeric", length =20)
    
    distPlot(0,0,0,0,x, z, as.character("theta"), TRUE)
  })
  
  theta_postDens <- reactive({
    z = priorParam()
    x = postParam() 
    n <- as.numeric(input$hypo_data_size)
    distPlot(n - input$n_cyc, input$mmf_succ, input$n_cyc, input$cyc_succ, x, z, as.character("theta"), FALSE)
  })
  
  ###########
  ## plots
  
  ## Outputting prior and posterior densities of pC, pE and theta.
  output$pc_density <- renderPlot({
    ## Plot the prior and posterior distributions for pC
    if((posterior40() | posterior20()) & !input$postsum){
      y = pC_priorDens()
      x = pC_postDens()
      
      if(posterior40()){
        leg1 = c(bquote("Posterior" ~ n[C] ==  .(input$n_cyc40) ~ S[C]  == .(input$cyc_succ40) ~ n[M] == .(40 - input$n_cyc40) ~ S[M] ==  .(input$mmf_succ40)), "Prior")
      }else{
        leg1 = c(bquote("Posterior" ~ n[C] ==  .(input$n_cyc20) ~ S[C]  == .(input$cyc_succ20) ~ n[M] == .(20 - input$n_cyc20) ~ S[M] ==  .(input$mmf_succ20)), "Prior")
      }
      
      postPlot = plot(x$gridc, x$postd, type="l", lty=2, lwd=3, col="dark red",
                      main = "Prior and posterior densities of 6-month remission rate on CYC",
                      xlab = "CYC 6-month remission rate", ylab="Density",
                      xlim = c(0,1), ylim = range(c(x$postd, y$dens)), cex.lab = 1.5, cex.axis=1.5, cex.main = 2) 
      lines(y$gridc, y$dens, lty=1, lwd=3, col="red")
      legend("topleft", as.expression(leg1), col = c("dark red","red"), lty = c(2,1), lwd = c(3,3), cex=1.2, bty="n")	
      
    }else if((posterior40() | posterior20()) & input$postsum){
      y = pC_priorDens()
      z = priorParam()
      w = postSumary()
      n <- as.numeric(input$hypo_data_size)
      dtacase = scan(glue::glue("data_scenario{n}.txt"))
      u1 = distPlot(n - dtacase[1], dtacase[3], dtacase[1], dtacase[2], w$scen1, z, as.character("pC"), FALSE)
      u2 = distPlot(n - dtacase[4], dtacase[6], dtacase[4], dtacase[5], w$scen2, z, as.character("pC"), FALSE)
      u3 = distPlot(n - dtacase[7], dtacase[9], dtacase[7], dtacase[8], w$scen3, z, as.character("pC"), FALSE)
      
      postPlot = plot(u1$gridc, u1$postd, type="l", lty=2, lwd=3, col="dark red",
                      main="Prior and posterior densities of 6-month remission rate on CYC/steroids",
                      xlab = "CYC 6-month remission rate", ylab="Density",
                      xlim =c(0,1), ylim = range(c(u1$postd, u2$postd, u3$postd, y$dens)), cex.lab = 1.5, cex.axis=1.5, cex.main = 2)  
      lines(u2$gridc, u2$postd, lty=3, lwd=3, col="dark red")
      lines(u3$gridc, u3$postd, lty=4, lwd=3, col="dark red")
      lines(y$gridc, y$dens, lty=1, lwd=3, col="red")
      
      leg1 = paste("CYC Succ = ", dtacase[2], "; MMF Succ = ", dtacase[3], sep="")
      leg2 = paste("CYC Succ = ", dtacase[5], "; MMF Succ = ", dtacase[6], sep="")
      leg3 = paste("CYC Succ = ", dtacase[8], "; MMF Succ = ", dtacase[9], sep="")
      
      legend("topleft", c(leg1,leg2, leg3,"Prior"), col = c("dark red","dark red","dark red","red"), lty = c(2,3,4,1), lwd = c(3,3,3,3), cex=1.2, bty="n")
      
    }else{
      x = pC_priorDens()
      ## Only want to plot the prior density in this scenario
      plot(x$gridc, x$dens, type="l", lty=1, lwd=3, col="red",
           main = "Prior density of 6-month remission rate on CYC/steroids",
           xlab = "CYC 6-month remission rate", ylab="Density", xlim =c(0,1), cex.lab = 1.5, cex.axis=1.5, cex.main = 2)	
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
                      main="Prior and posterior densities of 6-month remission rate on MMF",
                      xlab = "MMF 6-month remission rate", ylab="Density",
                      xlim =c(0,1), ylim = range(c(x$postd, y$dens)), cex.lab = 1.5, cex.axis=1.5, cex.main = 2) 
      lines(y$gride, y$dens, lty=1, lwd=3, col="green")
      legend("topleft", as.expression(leg1), col = c("dark green","green"), lty = c(2,1), lwd = c(3,3), cex=1.2, bty="n")	
    }else if((posterior40() | posterior20()) & input$postsum){
      y = pE_priorDens()
      z = priorParam()
      w = postSumary()
      n <- as.numeric(input$hypo_data_size)
      dtacase = scan(glue::glue("data_scenario{n}.txt"))
      u1 = distPlot(n - dtacase[1], dtacase[3], dtacase[1], dtacase[2], w$scen1, z, as.character("pE"), FALSE)
      u2 = distPlot(n - dtacase[4], dtacase[6], dtacase[4], dtacase[5], w$scen2, z, as.character("pE"), FALSE)
      u3 = distPlot(n - dtacase[7], dtacase[9], dtacase[7], dtacase[8], w$scen3, z, as.character("pE"), FALSE)
      
      postPlot = plot(u1$gride, u1$postd, type="l", lty=2, lwd=3, col="dark green",
                      main="Prior and posterior densities of 6-month remission rate on MMF",
                      xlab = "MMF 6-month remission rate", ylab="Density",
                      xlim =c(0,1), ylim = range(c(u1$postd, u2$postd, u3$postd, y$dens)),
                      cex.lab = 1.5, cex.axis=1.5, cex.main = 2)  
      lines(u2$gride, u2$postd, lty=3, lwd=3, col="dark green")
      lines(u3$gride, u3$postd, lty=4, lwd=3, col="dark green")
      lines(y$gride, y$dens, lty=1, lwd=3, col="green")
      
      leg1 = paste("CYC Succ = ", dtacase[2], "; MMF Succ = ", dtacase[3], sep="")
      leg2 = paste("CYC Succ = ", dtacase[5], "; MMF Succ = ", dtacase[6], sep="")
      leg3 = paste("CYC Succ = ", dtacase[8], "; MMF Succ = ", dtacase[9], sep="")
      
      legend("topleft", c(leg1, leg2, leg3, "Prior"), col = c("dark green","dark green","dark green","green"), lty = c(2,3,4,1), lwd = c(3,3,3,3), cex=1.2, bty="n")
    }else{
      x = pE_priorDens()
      ## Only want to plot the prior density in this scenario
      priorPlot = plot(x$gride, x$dens, type="l", lty=1, lwd=3, col="green",
                       main = "Prior density of 6-month remission rate on MMF",
                       xlab = "MMF 6-month remission rate", ylab="Density",
                       xlim =c(0,1), cex.lab = 1.5, cex.axis=1.5, cex.main = 2)	
    }
  })
  
  output$pepc_density <- renderPlot({		
    ## Plot the prior and posterior distributions for pE
    if((posterior40() | posterior20()) & !input$postsum){
      y = pC_postDens()
      x = pE_postDens()
      postPlot = plot(y$gridc, y$postd, type="l", lty=2, lwd=3, col="dark red",
                      main="Posterior densities of 6-month remission rates on MMF and CYC",
                      xlab = "6-month remission rate", ylab="Density",
                      xlim =c(0,1), ylim = range(c(y$postd, x$postd)), cex.lab = 1.5, cex.axis=1.5, cex.main = 2) 
      lines(x$gride, x$postd, lty=2, lwd=3, col="dark green")
      legend("topleft", c("CYC/steroids", "MMF/steroids"), col = c("dark red","dark green"), lty = c(2,2), lwd = c(3,3), cex=1.2, bty="n")	
    }else if((posterior40() | posterior20()) & input$postsum){
      z = priorParam()
      w = postSumary()
      ns <- as.numeric(input$hypo_data_size)
      dtacase = scan(glue::glue("data_scenario{ns}.txt"))
      u1pe = distPlot(ns - dtacase[1], dtacase[3], dtacase[1], dtacase[2], w$scen1, z, as.character("pE"), FALSE)
      u2pe = distPlot(ns - dtacase[4], dtacase[6], dtacase[4], dtacase[5], w$scen2, z, as.character("pE"), FALSE)
      u3pe = distPlot(ns - dtacase[7], dtacase[9], dtacase[7], dtacase[8], w$scen3, z, as.character("pE"), FALSE)
      
      u1pc = distPlot(ns - dtacase[1], dtacase[3], dtacase[1], dtacase[2], w$scen1, z, as.character("pC"), FALSE)
      u2pc = distPlot(ns - dtacase[4], dtacase[6], dtacase[4], dtacase[5], w$scen2, z, as.character("pC"), FALSE)
      u3pc = distPlot(ns - dtacase[7], dtacase[9], dtacase[7], dtacase[8], w$scen3, z, as.character("pC"), FALSE)
      
      leg1 = bquote("CYC:" ~ n[C] == .(dtacase[1]) ~ S[C]  == .(dtacase[2]) ~ n[M] == .(ns - dtacase[1]) ~ S[M] ==  .(dtacase[3]))
      leg2 = bquote("CYC:" ~ n[C] == .(dtacase[4]) ~ S[C]  == .(dtacase[5]) ~ n[M] == .(ns - dtacase[4]) ~ S[M] ==  .(dtacase[6]))
      leg3 = bquote("CYC:" ~ n[C] == .(dtacase[7]) ~ S[C]  == .(dtacase[8]) ~ n[M] == .(ns - dtacase[7]) ~ S[M] ==  .(dtacase[9]))
      leg4 = bquote("MMF:" ~ S[C]  == .(dtacase[2]) ~ S[M] ==  .(dtacase[3]))
      leg5 = bquote("MMF:" ~ S[C]  == .(dtacase[5]) ~ S[M] ==  .(dtacase[6]))
      leg6 = bquote("MMF:" ~ S[C]  == .(dtacase[8]) ~ S[M] ==  .(dtacase[9]))
      
      postPlot = plot(u1pc$gridc, u1pc$postd, type="l", lty=2, lwd=3, col="dark red",
                      main="Posterior densities of 6-month remission rates on MMF and CYC",
                      xlab = "6-month remission rate", ylab="Density",
                      xlim =c(0,1), ylim = range(c(u1pc$postd, u2pc$postd, u3pc$postd, u1pe$postd, u2pe$postd, u3pe$postd)),
                      cex.lab = 1.5, cex.axis=1.5, cex.main = 2)  
      lines(u2pc$gridc, u2pc$postd, lty=3, lwd=3, col="dark red")
      lines(u3pc$gridc, u3pc$postd, lty=4, lwd=3, col="dark red")
      lines(u1pe$gride, u1pe$postd, lty=2, lwd=3, col="dark green")
      lines(u2pe$gride, u2pe$postd, lty=3, lwd=3, col="dark green")
      lines(u3pe$gride, u3pe$postd, lty=4, lwd=3, col="dark green")
      
      legend("topleft", as.expression(c(leg1, leg2, leg3, leg4, leg5, leg6)), col = c("dark red","dark red","dark red", "dark green","dark green","dark green"), lty = c(2,3,4,2,3,4), lwd = c(3,3,3,3,3,3), cex=1.2, bty="n")
    }else{
      y = pC_priorDens()
      x = pE_priorDens()
      ## Only want to plot the prior density in this scenario
      priorPlot = plot(y$gridc, y$dens, type="l", lty=1, lwd=3, col="red",
                       main = "Prior density of 6-month remission rates on MMF and CYC",
                       xlab = "6-month remission rate", ylab="Density",
                       xlim =c(0,1), ylim = range(c(y$dens, x$dens)), cex.lab = 1.5, cex.axis=1.5, cex.main = 2)	
      lines(x$gride, x$dens, lty=1, lwd=3, col="green")
      legend("topleft", c("CYC/steroids", "MMF/steroids"), col = c("red","green"), lty = c(1,1), lwd = c(3,3), cex=1.2, bty="n")	
    }
  })
  
  output$theta_density <- renderPlot({	
    ## Plot the prior and posterior distributions for theta
    
    pc_to_pe_transformation <- function(pc, a)
      1.0/(1 + exp(-(a + log(pc/(1 - pc)))))
    
    pe_lim1 = pc_to_pe_transformation(input$pc_q1, -2)
    pe_lim2 = pc_to_pe_transformation(input$pc_q1, -1)
    pe_lim3 = pc_to_pe_transformation(input$pc_q1, 0)
    pe_lim4 = pc_to_pe_transformation(input$pc_q1, 1)
    pe_lim5 = pc_to_pe_transformation(input$pc_q1, 2)
    
    label1 <- paste("6-month remission rate on MMF assuming that on CYC is ", round(input$pc_q1, 2), sep="")
    
    if((posterior40() | posterior20()) & !input$postsum){
      y = theta_priorDens()
      x = theta_postDens()
      
      if(posterior40()){
        leg1 = c(bquote("Posterior" ~ n[C] ==  .(input$n_cyc40) ~ S[C]  == .(input$cyc_succ40) ~ n[M] == .(40 - input$n_cyc40) ~ S[M] ==  .(input$mmf_succ40)), "Prior")
      }else{
        leg1 = c(bquote("Posterior" ~ n[C] ==  .(input$n_cyc20) ~ S[C]  == .(input$cyc_succ20) ~ n[M] == .(20 - input$n_cyc20) ~ S[M] ==  .(input$mmf_succ20)), "Prior")
      }
      
      par(pty="m", plt=c(0.1, 1, 0, 1), omd=c(0.1,1.0,0,0.9))
      par(mfrow = c(2, 1))
      
      postPlot = plot(x$gridt, x$postd, type="l", lty=2, lwd=3, col="dark blue", main="", ylab="Density", xaxt = "n", cex.lab = 1.5, cex.axis=1.5, ylim = range(c(x$postd, y$dens))) 
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
      
      if(posterior40()){
        dtacase = scan("data_scenario40.txt")
        u1 = distPlot(40 - dtacase[1], dtacase[3], dtacase[1], dtacase[2], w$scen1, z, as.character("theta"), FALSE)
        u2 = distPlot(40 - dtacase[4], dtacase[6], dtacase[4], dtacase[5], w$scen2, z, as.character("theta"), FALSE)
        u3 = distPlot(40 - dtacase[7], dtacase[9], dtacase[7], dtacase[8], w$scen3, z, as.character("theta"), FALSE)
      }else{
        dtacase = scan("data_scenario20.txt")
        u1 = distPlot(20 - dtacase[1], dtacase[3], dtacase[1], dtacase[2], w$scen1, z, as.character("theta"), FALSE)
        u2 = distPlot(20 - dtacase[4], dtacase[6], dtacase[4], dtacase[5], w$scen2, z, as.character("theta"), FALSE)
        u3 = distPlot(20 - dtacase[7], dtacase[9], dtacase[7], dtacase[8], w$scen3, z, as.character("theta"), FALSE)
      }
      
      par(pty="m", plt=c(0.1, 1, 0, 1), omd=c(0.1,1.0,0,0.9))
      par(mfrow = c(2, 1))
      
      postPlot = plot(u1$gridt, u1$postd, type="l", lty=2, lwd=3, col="dark blue", main="", ylab="Density", xaxt="n", cex.lab = 1.5, cex.axis=1.5, ylim = range(c(u1$postd, u2$postd, u3$postd, y$dens)))
      lines(u2$gridt, u2$postd, lty=3, lwd=3, col="dark blue")
      lines(u3$gridt, u3$postd, lty=4, lwd=3, col="dark blue")
      lines(y$gridt, y$dens, lty=1, lwd=3, col="blue")
      
      leg1 = paste("CYC Succ = ", dtacase[2], "; MMF Succ = ", dtacase[3], sep="")
      leg2 = paste("CYC Succ = ", dtacase[5], "; MMF Succ = ", dtacase[6], sep="")
      leg3 = paste("CYC Succ = ", dtacase[8], "; MMF Succ = ", dtacase[9], sep="")
      
      legend("topright", c(leg1, leg2, leg3, "Prior"), col = c("dark blue","dark blue","dark blue","blue"), lty = c(2,3,4,1), lwd = c(3,3,3,3), cex=1.2, bty="n")	
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
      par(mfrow = c(2, 1))
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
      axis(side=1, at = c(-2, -1, 0, 1, 2), labels = round(c(pe_lim1, pe_lim2, pe_lim3, pe_lim4, pe_lim5), 2), col="red", col.axis="red", cex.axis=1.5)
      linloc <- par()$usr[3]
      abline(h=linloc, col="red")
      mtext(side=1, label1, line=3.0, cex=1.5, col="red")	
    }
  })
  
  ## Generate a summary of the prior or posterior densities as necessary
  output$summary = renderPrint({
    if(posterior40() | posterior20()){
      
      x = postParam()
      
      cat("Posterior distributions update the priors implied by the opinions: Q1 = ", input$pc_q1, ", Q2 = ", input$pc_q2, ", Q3 = ", input$theta_q1, ", Q4 = ", input$theta_q2, "\n")
      cat("\n")
      cat("Posterior distributions incorporate the following data: \n")
      if(posterior40()){
        cat(input$n_cyc40, "patients randomised to CYC of whom ", input$cyc_succ40, "were in remission within six months \n")
        cat(40 - input$n_cyc40, "patients randomised to MMF of whom ", input$mmf_succ40, "were in remission within six months \n")
      }else{
        cat(input$n_cyc20, "patients randomised to CYC of whom ", input$cyc_succ20, "were in remission within six months \n")
        cat(20 - input$n_cyc20, "patients randomised to MMF of whom ", input$mmf_succ20, "were in remission within six months \n")
      }
      cat("\n")
      cat("\n")
      cat("Summary of the posterior distribution of the 6-month remission rate on CYC/steroids: \n")
      cat("\n")
      cat("90% credibility interval for CYC remission rate is (", round(x[4], 2), ",", round(x[5], 2), ") \n") 
      cat("Your distribution has mode ", round(x[2], 2), " and mean ", round(x[1], 2), "and standard deviation ", round(x[3], 2))	
      cat("\n")
      cat("\n")
      cat("Summary of the posterior distribution of the 6-month remission rate on MMF/steroids: \n")
      cat("\n")
      cat("90% credibility interval for MMF remission rate is (", round(x[9], 2), ",", round(x[10], 2), ") \n")
      cat("Your distribution has mode ", round(x[7], 2), " and mean ", round(x[6], 2), "and standard deviation", round(x[8], 2), "\n")
      cat("\n")
      cat("\n")
      cat("Summary of the posterior distribution of the log-odds ratio: \n")
      cat("\n")
      cat("90% credibility interval for the log-odds ratio is (", round(x[14], 2), ",", round(x[15], 2), ") \n") 
      cat("Your distribution for the log-odds ratio has mode ", round(x[12], 2), " and mean ", round(x[11], 2), "and standard deviation ", round(x[13], 2), "\n")			
      cat("\n")
      cat("\n")
      cat("Posterior probability that the 6-month remission rate on MMF exceeds that on CYC is ", round(x[16], 2), "\n")
      cat("Posterior probability that the 6-month remission rate on MMF is at worst 10% less than that on CYC is", round(x[17], 2), "\n")
      cat("Posterior probability that the 6-month remission rate on CYC exceeds that on MMF by more than 10% is ", round(1 - x[17], 2), "\n")			
    }else{
      
      x = priorParam()
      cat("Prior distributions are determined by the opinions: Q1 = ", input$pc_q1, ", Q2 = ", input$pc_q2, ", Q3 = ", input$theta_q1, ", Q4 = ", input$theta_q2, "\n")
      cat("\n")
      cat("\n")
      cat("Summary of the prior distribution of the 6-month remission rate on CYC/steroids: \n")
      cat("\n")
      
      pcquantile0_05 = qbeta(0.05, shape1 = x[1], shape2 = x[2], lower.tail=TRUE)
      pcquantile0_95 = qbeta(0.95, shape1 = x[1], shape2 = x[2], lower.tail=TRUE)
      pc_mode = (x[1] - 1)/(x[1] + x[2] -2)
      pc_mean = x[1]/(x[1] + x[2])
      pc_var = (x[1]*x[2])/(((x[1] + x[2])^2)*(x[1] + x[2] +1))
      
      tquantile0_05 = qnorm(0.05, mean=x[9], sd = sqrt(x[10]), lower.tail=TRUE)
      tquantile0_95 = qnorm(0.95, mean=x[9], sd = sqrt(x[10]), lower.tail=TRUE)
      
      cat("90% credibility interval for CYC remission rate is (", round(pcquantile0_05, 2), " , ", round(pcquantile0_95, 2), ") \n")
      cat("Your distribution is Beta with parameters ", round(x[1], 2), " and ", round(x[2], 2), "\n")
      cat("with mode ", round(pc_mode, 3), " and mean ",round(pc_mean, 2), "and standard deviation ", round(sqrt(pc_var), 2), "\n")		
      cat("This information is equivalent to ", round(x[3], 0), "observations on CYC. \n")
      cat("\n")
      cat("\n")
      cat("Summary of the prior distribution of MMF/steroid 6-month remission rate: \n")
      cat("\n")
      cat("90% credibility interval for MMF remission rate is (", round(x[7], 2), ",", round(x[8], 2), ") \n")
      cat("Your distribution has mode ", round(x[5], 2), " and mean ", round(x[4], 2), " and standard deviation ", round(x[6], 2), "\n")
      cat("Proportion which we are 75% confident the true 6-month remission rate on MMF exceeds is ", round(x[14], 2), "\n")
      cat("\n")
      cat("\n")
      cat("Summary of the prior distribution of log-odds ratio : \n")
      cat("\n")
      cat("90% credibility interval for log-odds ratio is (", round(tquantile0_05, 2), ",", round(tquantile0_95, 2), ") \n")
      cat("Your distribution is Gaussian with mode ", round(x[9], 2), " and mean", round(x[9], 2), " and standard deviation ", round(sqrt(x[10]), 2), "\n")
      cat("This information is equivalent to what would be generated by a RCT comparing CYC with MMF recruiting ", round(0.5*x[13], 0), "patients per treatment arm \n")
      cat("\n")
      cat("\n")
      cat("Prior probability that the 6-month remission rate on MMF exceeds that on CYC is ", round(x[11], 2),  "\n")
      cat("Prior probability that the 6-month remission rate on MMF is at worst 10% less than that on CYC is", round(1 - x[12], 2), "\n")
      cat("Prior probability that the 6-month remission rate on CYC exceeds that on MMF by more than 10% is ", round(x[12], 2), "\n")
    }		
  })
})
