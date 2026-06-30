#==============================================================================#
#   GENE FLOW SHINY APP                                                        #
#==============================================================================#
#______________________________________________________________________________#
#                                                                              #

# Copyright L Stearman 2022


#==============================================================================#
#   ACQUISITION                                                                #
#==============================================================================#


	require(igraph)
	require(plotrix)
	require(shinybusy)

# Source codes
	source("GeneFlow.R")



#==============================================================================#
#   USER INTERFACE DEFINITION                                                  #
#==============================================================================#


# User Interface Definition
	ui <- navbarPage(title = "Gene Flow in Riverine Systems",

		tags$p("Copyright L. Stearman, 2022"),


		# tabPanel(title = "Introduction",
		# # Adds the main title
		# 	tags$h1("Sexy Fishes: Gene Flow Simulations in Complex Riverine Systems"),

		# # Adds a horizontal bar
		# 	tags$hr(),

		# 	tags$p("Explanatory Text")
		# ),

		tabPanel(title = "Simulation Settings",
			tags$hr(),
			tags$h2("Riverscape Settings"),


			fluidRow(
				column(width = 4,

					fileInput("map.file", label = "Please select an input map file",
						accept = ".csv", buttonLabel = "Browse..."),

					sliderInput("slider.md",
						label = "Maximum Disersal Distance (Patches)",
						min = 1, max = 5, value = 1, step = 1
					),


					tags$br(),
					tags$br(),
					tags$br(),

					tags$h4("Random Disturbance Events"),
					tags$p("These are random extinction events due to a density-independent disturbance which cause a local extinction. For no disturbances, leave patches unchecked. Extinctions can and do occur without disturbances due to population stochasticity."),
					sliderInput("slider.df",
						label = "Disturbance Frequency (Years)", min = 1,
						max = 100, val = 50, step = 1
					),

					tags$p("Enter patches with disturbances separated by a comma (example: 1, 2, 7), leave blank for none, or enter ALL for all patches."),
					textInput("text.disturb",
						label = "Patches with Disturbances",
						value = ""),

				),

				column(width = 4,

					tags$h4("Habitat Degradation Degree"),
					tags$p("This allows users to 'degrade' patches by consistently reducing the percent of K which a population can attain. Individuals can still disperse through these patches but will not reproduce in them."),

					sliderInput("slider.kr",
						label = "Percent Reduction in K", min = 0, max = 100,
						value = 0, step = 1
					),

					tags$p("Enter degraded patches separated by a comma (example: 1, 2, 7), leave blank for none, or enter ALL for all patches."),
					textInput("text.degrade",
						label = "Degraded Patches",
						value = ""),

					tags$br(),
					tags$h4("Habitat Patch Fragmentation Type"),
					tags$p("Users can specify units with a barrier at their downstream end and whether fragmentation is complete or upstream-only. Individuals cannot disperse past these barriers but can reproduce in the patch."),

					radioButtons("radio.fr", label = "Fragmentation Type",
						choiceNames = c("None", "Upstream-Only", "Total"),
						choiceValues = c(1, 2, 3), inline = TRUE
					),


					tags$p("Enter patches fragmented at their downstream end separated by a comma (example: 1, 2, 7), leave blank for none."),
					textInput("text.fragment",
						label = "Fragmented Patches (at downstream end)",
						value = "")
				),

				column(width = 4,

					tags$p("Fig 1. Unrooted dendrogram of habitat patches. Patches are scaled by link magnitude. This diagram reprojects and changes orientation with many changes to riverscape parameters but patch relationships are preserved."),

				# Adds the tree image
					plotOutput("map.plot", width = 400, height = 400),


					tags$p("Fig 2. Patch connectivity plot. Brighter colors indicate higher probabilities of dispersal."),

					plotOutput("move.plot", width = 400, height = 400)

				)
			),

			tags$hr(),
			tags$h2("Species Settings"),

			fluidRow(
				column(width = 4,

					tags$h4("Population and Metapopulation Dynamics"),

					sliderInput("slider.r",
						label = "Instantaneous Rate of Increase (r)",
						min = 1, max = 2.5, value = 1.5, step = 0.05
					),
					
					sliderInput("slider.mp",
						label = "Proportion of Dispersing Individuals",
						min = 0, max = 0.3, value = 0.05, step = 0.01
					),

					tags$br(),

					tags$h4("Genome Settings"),

					sliderInput("slider.lo",
						label = "Number of Loci", min = 5, max = 100, step = 1,
						value = 25
					),

					sliderInput("slider.mu",
						label = "Mutation Rate (10^X)", min = -10, max = -5,
						step = 0.1, value = -7.5
					)

				),

				column(width = 4,

					tags$h4("Habitat Specificity"),

					tags$p("Users can specify habitat specificity relative to stream link magnitude."),

					numericInput("num.idealk",
					    label = "Ideal Link Magnitude (Highest K)", value = 2,
					    min = 1
					),

					numericInput("num.ksd",
						label = "Variance of Habitat Quality (sd of K)",
						value = 3, min = 0, step = 0.25
					),

					numericInput("num.mink",
						label = "Minimum K in any plot",
						value = 25, min = 0
					),

					numericInput("num.maxk",
						label = "Maximum K in any plot",
						value = 125, min = 5
					),

					tags$br(),

					tags$p("Total watershed K"),

					verbatimTextOutput("ktotals"),

				),

				column(width = 4,
					plotOutput("habspecial")
				)
			),


			tags$hr(),
			tags$h2("Simulation Settings"),

			fluidRow(
				column(width = 6,

					radioButtons("radio.start",
						label = "Select starting configuration type",
						choiceNames = c("Start from scratch",
						"Start from a prior run"), choiceValues = c(1, 2),
						inline = TRUE),

					tags$br(),

					tags$h4("Start from scratch"),
					tags$p("Use this input to start from scratch. Colonization can be simulated by only stocking individuals in specific plots (e.g., headwaters or stream outlet. Make sure to check 'Start from Scratch' above or this will not run correctly."),
					
					numericInput("num.stock",
						label = "Initial number of individuals per plot",
						value = 50
					),

					tags$p("Enter patches to stock separated by a comma (example: 1, 2, 7). Enter ALL for all patches."),
					textInput("text.stock",
						label = "Patches to stock",
						value = "ALL"),

					tags$br(),

					tags$h4("Start from a prior run"),

					tags$p("Use this input to start from the end of a prior run. This can be used to do common burn-ins prior to variable treatments. Make sure to check 'Input a prior run' above or this will not run correctly."),
					fileInput("start.config", label = "Load previous run?",
						accept = ".rds", buttonLabel = "Browse..."

					),


				),

				column(width = 6, 

					tags$h4("Population Parameters"),
					tags$p("Simulations with more generations take much longer to run. Consider <= 250 generations for initial trial models. Full burn-in often occurs at around 1,000 generations."),

					numericInput("num.gen",
						label = "Number of Generations", value = 250
					),

					tags$h4("Genetic Data Sampling"),
					tags$p("Users can subsample genetic data to reduce output file size. Subsampling is always inclusive of the last generation."),
				
					radioButtons("radio.step",
						label = "Interval of Genetic Data Sampling and Storage",
						choiceNames = list("Every Generation", "5 Generations",
						"10 Generations", "25 Generations", "100 Generations",
						"Last Generation Only"), choiceValues = list(1, 2, 3,
						4, 5, 6), selected = 1
					),

				),

			),

			fluidRow(
				tags$h4("Prior Run Parameters"),
				tags$p("Make sure that any parameters you wish to hold constant are set to the same values above. Upoading an old run only affects starting numbers of individuals in each plot and starting genotypes."),
				tags$br(),
				tags$p("Number of Generations"),
				verbatimTextOutput("oldGen"),
				tags$br(),


				column(width = 6,
					tags$h4("Riverscape Parameters"),
					tableOutput("riverTable"),

				),

				column(width = 6,
					tags$h4("Species Parameters"),
					tableOutput("speciesTable")
				),

			),

		),

		tabPanel(title = "Laboratory",

		# Sets the record button
			actionButton(inputId = "popsim", label = "Simulate Metapopulations",
				style="color: #fff; background-color: #006633;
				    border-color: #2e6da4; font-size:200%"),

			add_busy_spinner(spin = "orbit", position = "top-left",
				margins = c(500, 100)),

			numericInput("num.pop",
			    label = "Please select a population to plot",
			    value = 1, min = 1
			   ),

			fluidRow(
				column(width = 6, 
						plotOutput("simhist", height = "600px"),
				),

				column(width = 6, 
					plotOutput("patchstats", height = "600px")

				),

			),


			tags$hr(),

			actionButton(inputId = "genesim", label = "Simulate Gene Flow",
				style="color: #fff; background-color: #006633;
				    border-color: #2e6da4; font-size:200%"),

			tags$br(),
			numericInput("num.iter",
				label = "Please select an output generation",
				value = 250, min = 1),

			tags$br(),


			fluidRow(
				column(width = 4, 
					plotOutput("FSTimage", width = 500, height = 500),
		        ),
		        column(width = 4,
		        	plotOutput("FSTmds", width = 500, height = 500),
		        ),
		        column(width = 4,
		        	plotOutput("isobydist", width = 500, height = 500),
		        ),



		    ),


			tags$hr(),

		# Data download button
			downloadButton("downloadData", "Download Simulation")

		)
	)






#==============================================================================#
#   SERVER DEFINITION                                                          #
#==============================================================================#


	server <- function(input, output) {


		options(shiny.maxRequestSize=30*1024^2)


#==============================================================================#
#   INPUT REACTIVES                                                            #



	# Parameters needed for map creation
		# Obtains the maximum movement distance
			max.disp <- reactive(input$slider.md)

		# Obtains rate of dispersal
			disp.rate <- reactive(input$slider.mp)

		# Obtains the disturbed patches
			disturbed <- reactive(input$text.disturb)
			dist.freq <- reactive(input$slider.df)

		# Obtains the degraded patches
			degraded <- reactive(input$text.degrade)
			deg.per <- reactive(input$slider.kr)

		# Obtains the fragmented patches
			fragmented <- reactive(input$text.fragment)
			frag.type <- reactive(input$radio.fr)



	# Parameters needed for metapopulation simulations
		# Stores the r value
			r.val <- reactive(input$slider.r)

		# Stores the ideal link magnitude
			idealk <- reactive(input$num.idealk)

		# Stores standard deviation of K around ideal link magnitude
			ksd <- reactive(input$num.ksd)

		# Stores minimum K value
			mink <- reactive(input$num.mink)

		# Stores maximum K value
			maxk <- reactive(input$num.maxk)

		# Obtains the number of generations
			gens <- reactive(input$num.gen)

		# Obtains stocking values
			num.stock <- reactive(input$num.stock)
			stock.plots <- reactive(input$text.stock)


    # Parameters needed for gene flow simulation

	# Stores number of loci and mutation rates
		loci <- reactive(input$slider.lo)
		muta <- reactive(input$slider.mu)


		radio.step <- reactive(input$radio.step)

		pop.2.plot <- reactive(input$num.pop)

		plot.iter <- reactive(input$num.iter)

		start.config <- reactive(input$start.config)

		start.type <- reactive(input$radio.start)


#==============================================================================#
#   COMPLEX REACTIVE OBJECTS                                                   #


	# Reactive map acquisition
		map.data <- reactive({
    		inFile <- input$map.file
    		if (is.null(inFile)) return(NULL)
    		raw.map <- read.csv(inFile$datapath, header = 1)
    		raw.map <- as.matrix(raw.map)
    		raw.map
  		})
  
	# Reactive old run acquisition
		old.run <- reactive({
			inFile <- input$start.config
			if (is.null(inFile)) return(NULL)
			readRDS(inFile$datapath)
		})


	# Map construction
		Map <- reactive({

			if (!disturbed() %in% c("", "Currently none...")) {
				if (disturbed() %in% c("all", "All", "ALL")) {
					Disturbed <- sort(unique(map.data()))
				} else {
				 	dist.raw <- disturbed()
				 	Disturbed <- as.numeric(strsplit(gsub(" ", "",
				 		dist.raw), ",")[[1]])
				}
			} else {
			 	Disturbed <- NULL
			}

			 if (!degraded() %in% c("", "Currently none...")) {
				if (degraded() %in% c("all", "All", "ALL")) {
					Degraded <- sort(unique(map.data()))
				} else {
				 	deg.raw <- degraded()
				 	Degraded <- as.numeric(strsplit(gsub(" ", "",
				 		deg.raw), ",")[[1]])
				}
			 } else {
			 	Degraded <- NULL
			 }

			 if (!fragmented() %in% c("", "Currently none...")) {
			 	frag.raw <- fragmented()
			 	Fragmented <- as.numeric(strsplit(gsub(" ", "",
			 		frag.raw), ",")[[1]])
			 } else {
			 	Fragmented <- NULL
			 }

			MaxDisp <- max.disp()

			map.out <- mapMaker(FromTo = map.data(), maxdist = MaxDisp,
				frag.type = frag.type(), fragments = Fragmented,
				disturbed = Disturbed,
				degraded = Degraded)

     		map.out

		})

	# Reactive event to create the carrying capacity matrix
		k.values <- reactive({
			link.mags <- V(Map()$network)$link.mag
			unit.ID <- as.numeric(V(Map()$network))
			Mean <- idealk()
			SD = ksd()
			Min = mink()
			Max = maxk()
			raw.dist <- dnorm(link.mags, Mean, SD)
			scaled.dist <- rescale(raw.dist, c(Min, Max))
			cbind(ID = unit.ID, Link = link.mags, K = scaled.dist)
		})

	# Reactive event to create the simulation K table
	# Rows are generations, columns are plots
		k.table <- reactive({

		# Stores the base K values
			k.base <- k.values()[, 3]

		# Stores the degraded units
			if (!degraded() %in% c("", "Currently none...")) {
				if (degraded() %in% c("all", "All", "ALL")) {
					Degraded <- sort(unique(map.data()))
				} else {
				 	deg.raw <- degraded()
				 	Degraded <- as.numeric(strsplit(gsub(" ", "",
				 		deg.raw), ",")[[1]])
				}
			} else {
			 	Degraded <- NULL
			}

		# Degrades habitats if necessary
			if (length(Degraded) > 0) {
				deg.val <- 1 - (deg.per()/100)
				degradeds <- which(k.values[, 1] %in% Degraded)
				k.base[degradeds] <- k.base[degradeds] * deg.val
			} else {

			}

		# Sets up the basic K matrix
			k.mat <- t(replicate(gens(), k.base))

		# Stores the units with disturbances
			if (!disturbed() %in% c("", "Currently none...")) {
				if (disturbed() %in% c("all", "All", "ALL")) {
					Disturbed <- sort(unique(map.data()))
				} else {
				 	dist.raw <- disturbed()
				 	Disturbed <- as.numeric(strsplit(gsub(" ", "",
				 		dist.raw), ",")[[1]])
				}
			} else {
			 	Disturbed <- NULL
			}

		# Sets K values to 0 for disturbance events
			if (length(Disturbed) > 0) {
				for (i in 1:length(Disturbed)) {
					focal.plot <- k.mat[, Disturbed[i]]
					focal.plot <- focal.plot * rbinom(length(focal.plot),
						1, (1 - (1/dist.freq())))
					k.mat[, Disturbed[i]] <- focal.plot
				}
			}

	# Generates poisson distributed K values
		pois.mat <- apply(k.mat, 2, function(x) sapply(x, rpois, n = 1))

	# Output
		pois.mat


	})

    # Reactive creation of starting population size
    	pop.start <- reactive({

    	# Stores the plots
    		start.pops <- vector("numeric", length(as.numeric(V(Map()$network))))
    		names(start.pops) <- as.numeric(V(Map()$network))

    	# Identifies the stocking plots
    		if (stock.plots() %in% c("", "ALL", "All", "all")) {
    			Stocks <- as.numeric(names(start.pops))
    		} else {
    			stock.raw <- stock.plots()
    			Stocks <- as.numeric(strsplit(gsub(" ", "",
    				stock.raw), ",")[[1]])
			}

		# Stocks the stock plots
			start.pops[names(start.pops) %in% Stocks] <- num.stock()

		# Output
			start.pops

    	})



    # Reactive creation of a steps list object
    	step.list <- reactive({
    		step.length <- c(1, 5, 10, 25, 100, gens())[as.numeric(radio.step())]
    		part.seq <- seq(from = 0, to = gens(), by = step.length,)
    		part.seq <- unique(sort(c(part.seq, gens())))
    		part.seq <- part.seq[part.seq > 0]
    		step.df <- data.frame(Input = 1:length(part.seq),
    		    Iteration = part.seq)
    		list(step.length = step.length, step.df = step.df)
    	})



	# Reactive calculation of pairwise FST
		pairwise.fst <- reactive({
			iters <- step.list()$step.df
			G <- iters[iters[, 2] %in% plot.iter(), 1]
			if (is.null(G)) {
				NULL
			} else {

				mod <- gene.mod()
				pops <- mod[[G]]
				OmniFSTPair(pops)
			}
		})



#==============================================================================#
#   MODELS                                                                     #



	# Metapopulation Simulation
		pop.mod <- eventReactive(input$popsim,{

			if (is.null(old.run()) | start.type() == 1) {
				N.start <- pop.start()
			} else {
				N.start <- tail(old.run()$pop.hist$N.pre, 1)
			}

			out.mod <- metaPop(N = N.start,
				r = r.val(),
				K.mat = k.table(),
				G = gens(),
				move.rate = disp.rate(),
				move.map = Map()$move.prob)

			out.mod

		})


		gene.mod <- eventReactive(input$genesim, {


			if (is.null(old.run()) | start.type() == 1) {
				pop.p <- matrix(runif(length(pop.start()) * loci()),
					loci(), length(pop.start()))
				Genotypes <- genoType(N = pop.start(), p.freq = pop.p)

			} else {
				all.geno <- old.run()[[4]]
				gene.arrays <- all.geno[[length(all.geno)]]
				Genotypes <- lapply(gene.arrays, function(x) {
					lapply(seq(dim(x)[3]), function(y) x[, , y])
				})
			}
####
			progress <- shiny::Progress$new()
			progress$set(message = "Computing data", value = 0)
			on.exit(progress$close())

			updateProgress <- function(value = NULL, detail = NULL) {
				if (is.null(value)) {
					value <- progress$getValue()
				}
				progress$set(value = value, detail = detail)
			}

			gene.sim <- geneFlow(pop.hist = pop.mod(), genotypes = Genotypes,
				Step = step.list()$step.length, Mutation = 10^muta(),
				updateProgress)

#####
			gene.sim

		})






#==============================================================================#
#   OUTPUT OBJECTS                                                             #


		output$ktotals <- renderPrint({
			if (is.null(map.data())) {
				"Please input a map file..."
			} else {
				round(sum(k.values()[, 3]))
			}
		})




	# Creates the habitat specialization plot
		output$habspecial <- renderPlot({
			if (is.null(map.data())) {
				plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE,
					ann = FALSE)
				text(0, 0, "Please input a map file...")

			} else {

				kvals <- k.values()
				kvals <- unique(data.frame(kvals))
				kvals <- kvals[order(kvals$Link), ]

				plot(kvals$Link, kvals$K, type = "o", ann = FALSE,
					ylim = c(0, maxk()))
				mtext("Stream Link Mangitude", side = 1, line = 3)
				mtext("Carrying Capacity", side = 2, line = 3)
			}

		})


	# Map plotting
		output$map.plot <- renderPlot({
			if (is.null(map.data())) {
				plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE,
					ann = FALSE)
				text(0, 0, "Please input a map file...")

			} else {

			par(mar = rep(0, 4))
			net <- Map()$network
			plot(net, vertex.size = rescale(V(net)$link.mag, c(8, 15)),
			    edge.color = E(net)$frag.col,
				edge.lty = E(net)$frag.lty,
				vertex.shape = V(net)$disturb.shape,
				vertex.color = V(net)$degraded.col,
				vertex.label.color = "black")
			}
		})

		output$move.plot <- renderPlot({

			if (is.null(map.data())) {
				plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE,
					ann = FALSE)
				text(0, 0, "Please input a map file...")
				
			} else {

					colRamp <- colorRampPalette(c("#00FF0080", "#00FF00FF"),
					    alpha = TRUE)
					par(mar = c(4.5, 4.5, 3, 3))

					moves <- Map()$move.prob
					moves <- ifelse(moves == 0, NA, moves)

					plot.seq <- ncol(moves)

					ax.vals <- seq(0, 1, length.out = plot.seq)
					ax1 <- ax.vals[seq(1, plot.seq, by = 2)]
					ax2 <- ax.vals[seq(2, plot.seq, by = 2)]

					image(t(moves), axes = FALSE, ann = FALSE)
					rect(par("usr")[1], par("usr")[3], par("usr")[2],
					    par("usr")[4], col = "black")
					image(t(moves), col = colRamp(10), add = TRUE)

				grid.vals <- seq(par("usr")[1], par("usr")[2],
				    length.out = plot.seq + 1)

				abline(v = grid.vals, h = grid.vals, lwd = 0.5, col = "grey")

				axis(side = 1, at = ax1, line = -0.5, cex.axis = 0.8,
					labels = as.character(1:plot.seq)[seq(1, 25, by = 2)],
					tick = FALSE)
				axis(side = 1, at = ax2, line = 0.5, cex.axis = 0.8,
					labels = as.character(1:plot.seq)[seq(2, 25, by = 2)],
					tick = FALSE)

				axis(side = 2, at = ax1, line = -0.5, cex.axis = 0.8,
					labels = as.character(1:plot.seq)[seq(1, 25, by = 2)],
					tick = FALSE)
				axis(side = 2, at = ax2, line = 0.5, cex.axis = 0.8,
					labels = as.character(1:plot.seq)[seq(2, 25, by = 2)],
					tick = FALSE)

				axis(side = 3, at = ax1, line = -0.5, cex.axis = 0.8,
					labels = as.character(1:plot.seq)[seq(1, 25, by = 2)],
					tick = FALSE)
				axis(side = 3, at = ax2, line = 0.5, cex.axis = 0.8,
					labels = as.character(1:plot.seq)[seq(2, 25, by = 2)],
					tick = FALSE)

				axis(side = 4, at = ax1, line = -0.5, cex.axis = 0.8,
					labels = as.character(1:plot.seq)[seq(1, 25, by = 2)],
					tick = FALSE)
				axis(side = 4, at = ax2, line = 0.5, cex.axis = 0.8,
					labels = as.character(1:plot.seq)[seq(2, 25, by = 2)],
					tick = FALSE)

				mtext("To Patch", side = 1, line = 3)
				mtext("From Patch", side = 2, line = 3)
			}
		})


		output$simhist <- renderPlot({
			if (is.null(map.data())) {

				plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE,
					ann = FALSE)
				text(0, 0,
				    "Please input a map file...")
				box()
	
			} else {

				pop.id <- as.numeric(pop.2.plot())

				par(mfrow = c(3, 1))
				par(mar = c(2, 4.5, 0.5, 0.5))
				par(oma = c(2.5, 0, 0, 0))

				plot(apply(pop.mod()$N.pre, 1, function(x) sum(x > 0)),
					type = "l", ann = FALSE)
				mtext("Occupied Patches", side = 2, line = 3)

				plot(apply(pop.mod()$N.pre, 1, sum), type = "l", ann = FALSE)
				mtext("Individuals in Simulation", side = 2, line = 3)

				if (pop.id <= ncol(pop.mod()$N.pre)) {
					plot(pop.mod()$N.pre[, pop.id], type = "l", ann = FALSE)
					mtext(paste("Individuals in Plot ", pop.id, sep = ""),
						side = 2, line = 3)
					mtext("Generation", side = 1, line = 3)
				} else {
					plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE,
						ann = FALSE)
				text(0, 0,
				    "Please select an existing population.", cex = 3)

				}
			}
		})



		output$patchstats <- renderPlot({

			if (is.null(map.data())) {

				plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE,
					ann = FALSE)
				text(0, 0,
				    "Please input a map file...")
				box()
	
			} else {

				production <- pop.mod()$N.pre[-nrow(pop.mod()$N.pre), ]

				production <- apply(production, 2, function(x) {
					xcur <- 0
					xi <- 1
					while(xcur == 0 & xi <= length(x)) {
						xcur <- x[xi]
						if (xcur == 0) {
							x[xi] <- NA
						}
						xi <- xi + 1
					}
					x
				})


				emigrants <- tapply(pop.mod()$move.list[, 3],
				    list(factor(pop.mod()$move.list[, 1],
				    levels = 1:nrow(production)),
					factor(pop.mod()$move.list[, 3], levels = 1:ncol(production))),
					length)

				emigrants <- ifelse(is.na(emigrants), 0, emigrants)

				immigrants <-tapply(pop.mod()$move.list[, 4],
				    list(factor(pop.mod()$move.list[, 1],
				    levels = 1:nrow(production)),
					factor(pop.mod()$move.list[, 4], levels = 1:ncol(production))),
					length)

				immigrants <- ifelse(is.na(immigrants), 0, immigrants)


				prop.occ <- apply(ifelse(production > 0, 1, production), 2,
					mean, na.rm = TRUE)
				prop.occ <- ifelse(is.na(prop.occ), 0, prop.occ)

				num.extinct <- apply(production, 2, function(x) {
					sum(diff(ifelse(x > 1, 1, x)) %in% -1, na.rm = TRUE)
				})

				emi.imma <- log10(((apply(emigrants, 2, mean) + 0.1) / 
					(apply(immigrants, 2, mean) + 0.1)))

				par(mfrow = c(3, 1))
				par(mar = c(2, 5.5, 0.5, 0.5))
				par(oma = c(2.5, 0, 0, 0))

				# barplot(prod.sum, ylim = c(0, max(prod.sum) * 1.05), ann = FALSE)
				# box()
				# mtext("Proportion of Generations\nwhere Immigration > Births", side = 2, line = 3.5)


				barplot(prop.occ, ylim = c(0, 1.05), ann = FALSE)
				box()
				mtext("Proportion of Time Occupied\nPost 1st Colonization",
					side = 2, line = 2.5)

				barplot(num.extinct, ylim = c(0, max(num.extinct) * 1.05),
					ann = FALSE)
				box()
				mtext("Number of Extinctions", side = 2, line = 3)

				barplot(emi.imma, ylim = c(min(emi.imma,na.rm = TRUE) * 1.05, max(emi.imma, na.rm = TRUE) * 1.05),
					ann = FALSE)
				box()
				mtext(expression(paste(Log[10], " Emigration:Immigration", sep = "")),
					side = 2, line = 3.5)
				mtext("Source", side = 2, line = 2.5, adj = 1)
				mtext("Sink", side = 2, line = 2.5, adj = 0)
				mtext("Patch", side = 1, line = 3)
			}


		})

		output$oldGen <- renderPrint({
			if (!is.null(old.run())) {
				nrow(old.run()[[3]]$N.pre)
			}
		})

		output$riverTable <- renderTable({
			if (!is.null(old.run())) {
				river.list <- old.run()[[1]]
				river.df <- do.call("rbind", river.list)
				river.df <- data.frame("Parameter" = rownames(river.df),
					"Value" = river.df)
				river.df[3, 2] <- ifelse(river.df[3, 2] == "", "None",
					river.df[3, 2])
				river.df[5, 2] <- ifelse(river.df[5, 2] == "", "None",
					river.df[5, 2])
				river.df[7, 2] <- ifelse(river.df[7, 2] == "", "None",
					river.df[7, 2])

				river.df
			} else {
				
			}
		})


		output$speciesTable <- renderTable({
			if (!is.null(old.run())) {
				species.list <- old.run()[[2]]
				species.df <- do.call("rbind", species.list)
				species.df <- data.frame("Parameter" = rownames(species.df),
					"Value" = species.df)
				species.df
			} else {

			}
		})




		output$FSTimage <- renderPlot({

			maxG <- nrow(pop.mod()$N.pre)

			iters <- step.list()$step.df

			if (sum(iters[, 2] %in% plot.iter()) == 0) {

				plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE,
			    	ann = FALSE)
		    	text(0, 0, "Please select a valid generation.\nRemember that save intervals may not be every generation.")
		   		box()

			} else {

			G <- iters[iters[, 2] %in% plot.iter(), 1]
			cur.pops <- pop.mod()$N.pre[plot.iter(), ]

				if (sum(cur.pops > 0) < 1) {
					plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE,
					    ann = FALSE)
				    text(0, 0, "There were no survivors...")
		 		    box()
				} else {
					if (sum(cur.pops > 0) < 2) {
						plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1),
						    axes = FALSE, ann = FALSE)
			   			text(0, 0,
			   			    expression(paste("Only one population, no pairwise ",
			    			F[ST], " calculated.")))
				    	box()
			

					} else {
				# Extracts the pairwise FST matrix
					fst.mat <- as.matrix(pairwise.fst())
			
				# Removes the upper triangle
					fst.mat[upper.tri(fst.mat)] <- NA

				# Plots the plot
					image(fst.mat, axes = FALSE, ann = FALSE)
					rect(par("usr")[1], par("usr")[3], par("usr")[2],
					    par("usr")[4], col = "black")
					image(fst.mat, add = TRUE,
						col = heat.colors(100), rev = TRUE)
					box()

				# Sets the axis labels
					ax.labs <- rownames(fst.mat)
					x.vals <- seq(0, 1, length.out = length(ax.labs))
					y.vals <- rev(x.vals)

					ind1 <- seq(1, length(x.vals), by = 2)
					ind2 <- seq(2, length(x.vals), by = 2)

					axis(side = 1, at = x.vals[ind1], labels = ax.labs[ind1],
						tick = FALSE, line = -0.5)
					axis(side = 1, at = x.vals[ind2], labels = ax.labs[ind2],
						tick = FALSE, line = 0.5)

					axis(side = 2, at = x.vals[ind1], labels = ax.labs[ind1],
						tick = FALSE, line = -0.5)
					axis(side = 2, at = x.vals[ind2], labels = ax.labs[ind2],
						tick = FALSE, line = 0.5)

					mtext("Patch", side = 1, line = 3)
					mtext("Patch", side = 2, line = 3)
					mtext(expression(paste("Pairwise ", F[ST], sep = "")),
						line = 0.5, side = 3, cex = 1.5)

					par.usr <- par("usr")
					xdist <- par.usr[2] - par.usr[1]
					ydist <- par.usr[4] - par.usr[3]

     				gradient.rect(par.usr[1] + (xdist * 0.1),
					    par.usr[3] + (ydist * 0.7), par.usr[1] + (xdist * 0.2),
					    par.usr[3] + (ydist * 0.95), col = heat.colors(30),
						gradient = "y")

     				text(par.usr[1] + (xdist * 0.15),
     				    par.usr[3] + (ydist * 0.98), round(max(fst.mat,
     				    	na.rm = TRUE), 3), col = "white")
     				text(par.usr[1] + (xdist * 0.15),
     				    par.usr[3] + (ydist * 0.67), round(min(fst.mat,
     				    	na.rm = TRUE), 3), col = "white")



				}

			}
		}
	})


		output$FSTmds <- renderPlot({

			iters <- step.list()$step.df
			G <- iters[iters[, 2] %in% plot.iter(), 1]

			cur.pop <- pop.mod()$N.pre[plot.iter(), ]

			if (sum(iters[, 2] %in% plot.iter()) == 0) {
				plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE,
			    	ann = FALSE)
		    	text(0, 0, "Please select a valid generation.\nRemember that save intervals may not be every generation.")
		   		box()
	
			} else {

				if (sum(cur.pop > 0) < 3) {
					plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE,
				    	ann = FALSE)
			    	text(0, 0, "Fewer than three populations,\nNo PCOA calculated")
			   		box()

				} else {

					sol <- cmdscale(pairwise.fst())
					plot(sol, cex = 0, ann = FALSE)
					text(sol, rownames(sol))
					mtext("PCOA 1", side = 1, line = 3)
					mtext("PCOA 2", side = 2, line = 3)
					mtext(expression(atop("Principle Coordinates Analysis of",
						paste(" Pairwise ", F[ST], sep = ""))), side = 3,
						line = -0.25, cex = 1.5)
				}
			}
		})


		output$isobydist <- renderPlot({


			iters <- step.list()$step.df
			G <- iters[iters[, 2] %in% plot.iter(), 1]

			cur.pop <- pop.mod()$N.pre[plot.iter(), ]

			if (sum(iters[, 2] %in% plot.iter()) == 0) {
				plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE,
			    	ann = FALSE)
		    	text(0, 0, "Please select a valid generation.\nRemember that save intervals may not be every generation.")
		   		box()
	
			} else {

				if (sum(cur.pop > 0) < 3) {
					plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE,
				    	ann = FALSE)
			    	text(0, 0, "Fewer than three populations,\nNo PCOA calculated")
			   		box()
			   	} else {

					map.dists <- as.matrix(Map()$dists)
					fst.dists <- as.matrix(pairwise.fst())

					keeps <- rownames(map.dists) %in% rownames(fst.dists)
					map.dists <- map.dists[keeps, keeps]

					map.dists <- as.dist(map.dists)
					fst.dists <- as.dist(fst.dists)

					plot(c(map.dists), c(fst.dists), ann = FALSE, pch = 21,
						col = rgb(0, 0, 0, 0.5), bg = rgb(0, 0, 0, 0.25), cex = 1.5)
					mtext("Distance Between Plots", side = 1, line = 2.5)
					mtext(expression(paste("Pairwise ", F[ST], sep = "")),
						side = 2, line = 2.5)
					mtext("Isolation by Distance", side = 3, line = 1,
					    cex = 1.5)
				}
			}

		})



		output$allelerich <- renderPlot({




			})










	# Defines the data download event
		output$downloadData <- downloadHandler(
			filename = paste("Gene Flow Simulation ",
				Sys.time(), ".rds", sep = ""),
			content = function(file) {

				riverscape.list <- list("Riverscape File" = input$map.file$name,
					"Max Dispersal Distance" = max.disp(),
					"Disturbed Patches" = disturbed(),
					"Disturbance Frequency" = dist.freq(),
					"Degraded Patches" = degraded(),
					"Percent Habitat Degradation" = deg.per(),
					"Fragmented Patches" = fragmented(),
					"Fragmentation Type" = c("None", "Upstream-Only",
						"Total")[as.numeric(frag.type())])

				species.list <- list("Instantaneous Rate of Increase" = r.val(),
					"Proportion of Dispersers" = disp.rate(),
					"Optimal Stream Order (Mean)" = idealk(),
					"Niche Width (SD)" = ksd(),
					"Minimum K" = mink(),
					"Maximum K" = maxk(),
					"Number of Loci" = loci(),
					"Mutation Rate (10^X)" = muta())

				saveRDS(list(riverscape = riverscape.list,
					species = species.list, pop.hist = pop.mod(),
					GeneSim = gene.mod()), file)
			},
			contentType = "rds"
    	)






	}














#==============================================================================#
#   SERVER EXECUTION                                                           #
#==============================================================================#


	shinyApp(ui = ui, server = server)



#==============================================================================#
#   END OF SCRIPT                                                              #
#==============================================================================#