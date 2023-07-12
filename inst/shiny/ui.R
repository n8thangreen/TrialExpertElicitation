
shinyUI(pageWithSidebar(

	headerPanel("Prior Specification"),
	
	sidebarPanel(
	
		helpText("Before any data are observed, please answer the following questions to specify your prior distributions: "),
		br(),
		
		textInput("expert", "Expert's name: ", ""),
		br(),
		br(),
						
		sliderInput(inputId = "pc_q1",
				label = "Q1: What do you think the 6-month remission rate for children with PAN treated with control arm in combination with steroids is? ",
				min = 0, max = 1, value = 0.7, step=0.05),
		br(),
		br(),		
		sliderInput(inputId = "pc_q2",
				label = "Q2: Provide a proportion such that you are 75% sure that the true 6-month remission rate on control arm/steroids exceeds this value",
				min = 0, max=1, value = 0.5, step=0.05),
		br(),
		br(),
		br(),
		br(),
		
		helpText("Because of the unpleasant side-effects of control arm, Experimental might be considered the preferable treatment even if it is associated with a somewhat lower 6-month remission rate: "),
		
		br(),
		br(),
		sliderInput(inputId = "theta_q1",
				label = "Q3: What is the chance that the 6-month remission rate on Experimental/steroids is higher than that on control arm/steroids?",
				min = 0, max=1, value = 0.3, step=0.05),
		br(),
		br(),

		sliderInput(inputId = "theta_q2",
				label = "Q4: What is the chance that the 6-month remission rate on control arm/steroids exceeds that on Experimental/steroids by more than 10%?",
				min = 0, max=1, value = 0.3, step=0.05),
		br(),
		br(),	

		radioButtons(inputId="hypo_data_size", label="Update MYPAN prior distributions with a hypothetical dataset of patient size ", 
		             choices=c(20, 40), selected = character(0)),

	## Output two tabs: one to plot the density of the (prior/posterior) distribution
	## and one to summarise the (prior/posterior) distribution
	
	br(),		
  
  ## Display this only if data are available to update the prior distribution						
	conditionalPanel(condition = "input.posterior40 == true & input.posterior20 == false",
			numericInput(inputId = "n_cyc40", label = "Number randomized to CYC (out of 40 patients)", value=20, min = 1, max= 40,  step=1),
			br(),
			numericInput(inputId = "cyc_succ40", label = "Number of successes on CYC:", value=14, min = 0, max = 40, step =1),
			br(),
			br(),
			helpText("For the remaining patients randomized to MMF: "),
			numericInput(inputId = "mmf_succ40", label = "Number of successes on MMF:", value=14, min = 0, step=1)
	),
	
	conditionalPanel(condition = "input.posterior20 == true & input.posterior40 == false",
			numericInput(inputId = "n_cyc20", label = "Number randomized to CYC (out of 20 patients)", value=10, min = 1, max= 20,  step=1),
			br(),
			numericInput(inputId = "cyc_succ20", label = "Number of successes on CYC:", value=7, min = 0, max = 40, step =1),
			br(),
			br(),
			helpText("For the remaining patients randomized to MMF: "),
			numericInput(inputId = "mmf_succ20", label = "Number of successes on MMF:", value=7, min = 0, step=1)
	),
	
	conditionalPanel(condition = "input.hypo_data_size == 20 | input.hypo_data_size == 40",
	                 numericInput(inputId = "n_cyc", label = "Number randomized to CYC", value=10, min = 1, max= 20,  step=1),
	                 br(),
	                 numericInput(inputId = "cyc_succ", label = "Number of successes on CYC:", value=7, min = 0, max = 40, step =1),
	                 br(),
	                 helpText("For the remaining patients randomized to MMF: "),
	                 numericInput(inputId = "mmf_succ", label = "Number of successes on MMF:", value=7, min = 0, step=1)
	),
	
	# read in (hypothetical) trial data
	# rows grouped
	conditionalPanel(condition = "input.posterior40 == true | input.posterior20 == true",
			checkboxInput(inputId = "postsum", 
					label = strong("Summarise posterior distributions for a variety of data scenarios"),
					value = FALSE)
	),
	
	submitButton("Update Bayesian Distributions")

	),
	
	mainPanel(
		tabsetPanel(
			tabPanel(title = "Density: CYC Remission Rate", plotOutput(outputId = "pc_density", height="700px")),
			tabPanel(title = "Density: MMF Remission Rate", plotOutput(outputId = "pe_density", height="700px")),
			tabPanel(title = "Density: CYC & MMF Remission", plotOutput(outputId = "pepc_density", height="700px")),
			tabPanel(title = "Density: Log-odds Ratio", plotOutput(outputId = "theta_density", height="700px")),
			tabPanel(title = "Summary", htmlOutput("summary"),
			         tags$head(tags$style("#summary{color: black;
                                 font-size: 15px;
                                 }")
			         )
			)
		)
	)	
))
	