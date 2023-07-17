
shinyUI(pageWithSidebar(

	headerPanel("Prior Specification"),
	
	sidebarPanel(
	
		helpText("Before any data are observed, please answer the following questions to specify your prior distributions: "),
		br(),
		
		textInput("expert", "Expert's name: ", ""),
		br(),
		br(),
						
		sliderInput(inputId = "pc_q1",
				label = "Q1: What do you think the 12-month response rate (clinically inactive disease off steroids) will be for the methotrexate group?",
				min = 0, max = 1, value = 0.7, step=0.05),
		br(),
		br(),		
		sliderInput(inputId = "pc_q2",
				label = "Q2: What is the 12-month response rate (clinically inactive disease off steroids) for the methotrexate group that you are at least 75% sure that will be achieved?",
				min = 0, max=1, value = 0.5, step=0.05),
		br(),
		br(),
		br(),
		br(),
		
		helpText("Based on the case reports we have presented, we hypothesise that baricitinib is superior to methotrexate for the treatment of JDM."),
		
		br(),
		br(),
		sliderInput(inputId = "theta_q1",
				label = "Q3: What is the chance that the response rate in the baricitinib arm is higher than that in methotrexate group?",
				min = 0, max = 100, value = 30, step=5),
		br(),
		br(),

		sliderInput(inputId = "theta_q2",
		            ##TODO: what is the %?
				label = "Q4: What is the chance that the response rate on baricitinib exceeds that on methotrexate by more than %?",
				# label = "Q4: What is the improvement you will be 75% confident it will achieved?",
				# this is the AUC on the left
				# should we use the point value instead?
				# how to modify the code?
				min = 0, max = 100, value = 30, step=5),
		br(),
		br(),	

		radioButtons(inputId="hypo_data_size", label="Update prior distributions with a hypothetical dataset of patient size ", 
		             choices=c(20, 40), selected = character(0)),

	## Output two tabs: one to plot the density of the (prior/posterior) distribution
	## and one to summarise the (prior/posterior) distribution
	
	br(),		
  
  ## Display this only if data are available to update the prior distribution						
	conditionalPanel(condition = "input.posterior40 == true & input.posterior20 == false",
			numericInput(inputId = "n_cyc40", label = "Number randomized to Control arm (out of 40 patients)", value=20, min = 1, max= 40,  step=1),
			br(),
			numericInput(inputId = "cyc_succ40", label = "Number of successes on Control arm:", value=14, min = 0, max = 40, step =1),
			br(),
			br(),
			helpText("For the remaining patients randomized to Experimental arm: "),
			numericInput(inputId = "mmf_succ40", label = "Number of successes on Experimental arm:", value=14, min = 0, step=1)
	),
	
	conditionalPanel(condition = "input.posterior20 == true & input.posterior40 == false",
			numericInput(inputId = "n_cyc20", label = "Number randomized to Control arm (out of 20 patients)", value=10, min = 1, max= 20,  step=1),
			br(),
			numericInput(inputId = "cyc_succ20", label = "Number of successes on Control arm:", value=7, min = 0, max = 40, step =1),
			br(),
			br(),
			helpText("For the remaining patients randomized to Experimental arm: "),
			numericInput(inputId = "mmf_succ20", label = "Number of successes on Experimental arm:", value=7, min = 0, step=1)
	),
	
	conditionalPanel(condition = "input.hypo_data_size == 20 | input.hypo_data_size == 40",
	                 numericInput(inputId = "n_cyc", label = "Number randomized to Control arm", value=10, min = 1, max= 20,  step=1),
	                 br(),
	                 numericInput(inputId = "cyc_succ", label = "Number of successes on Control arm:", value=7, min = 0, max = 40, step =1),
	                 br(),
	                 helpText("For the remaining patients randomized to Experimental arm: "),
	                 numericInput(inputId = "mmf_succ", label = "Number of successes on Experimental arm:", value=7, min = 0, step=1)
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
			tabPanel(title = "Density: Control Arm Remission Rate", plotOutput(outputId = "pc_density", height="700px")),
			tabPanel(title = "Density: Experimental Arm Remission Rate", plotOutput(outputId = "pe_density", height="700px")),
			tabPanel(title = "Density: Control Arm & Experimental Arm Remission", plotOutput(outputId = "pepc_density", height="700px")),
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
	