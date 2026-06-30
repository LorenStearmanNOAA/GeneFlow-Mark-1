#==============================================================================#
#   GENE FLOW SIMULATOR AND ANALYSIS TOOLS                                     #
#==============================================================================#
#______________________________________________________________________________#
#                                                                              #

	

kMat <- function(k.mean, k.adj, Gens, Map, degraded, dper, disturbed, dfreq) {
# A function to calculate a matrix of K values
#
# Args:
# 
# Returns:
#
# Stores base K for each plot adjusted for habitat specificity
	base.vals <-rep(k.mean, ncol(Map)) * k.adj

# Accounts for baseline habitat degradation
	if (length(degraded) > 0) {
		deg.vec <- rep(1, ncol(Map))
		deg.vec[degraded] <- dper
		base.vals <- base.vals * deg.vec
	}

# Expands the values into a matrix
	base.mat <- t(replicate(Gens, base.vals))

# Accounts for disturbance events
	if (length(disturbed) > 0) {
		base.mat <- apply(base.mat, 2, function(x){
		    disturbs <- rbinom(length(x), 1, (1 - (1/dfreq)))
		    x * disturbs
		})
	}

# Generates poisson distributed values
	pois.mat <- apply(base.mat, 2, function(x) sapply(x, rpois, n = 1))

# Output
	pois.mat
}





popGrow <- function(r, n, k, mfc) {
# A function to calculate population growth given basic demographic parameters
#
# Args:
#	r: The instantaneous rate of increase
#	n: The number of individuals
#	k: Carrying capacity
#	mfc: male female count. This prevents populations of only one sex from
# 		reproducing
#
# Returns a new population size

	if (k > 0 & min(mfc) > 0) {
		new <- (r * n) * ((k - n)/k)
		max(c(new + n, 0))
	} else {
		0
	}
}





mapMaker <- function(FromTo, maxdist = 1, frag.type = 1, fragments = NULL,
    disturbed = NULL, degraded = NULL) {
# A function to generate a transitional probability matrix, representing
# dispersal probabilities among units
#
# Args:
#	FromTo: A mx2 matrix specifying hierarchical upstream/downstream
#           relationships between adjacent units
#   maxdist: The maximum distance dispersal permissible
#	frag.type: The fragmentation type. 1 = none, 2 = upstream only, 3 = total
#	fragments: Units where fragmentation occurs at their downstream boundary
#	disturbed: Units with periodic disturbances
#	degraded: Units with "degraded" habitat conditions (lower K)
#
# Returns a matrix of dispersal probabilities

# Identifies the unique units
	units <- sort(unique(c(FromTo[, 1], FromTo[, 2])))

# Determines the direct downstreams
	downstreams <- sapply(units, function(x) {
		curID <- x
		id.vec <- x
		while(!is.na(curID)) {
			curID <- FromTo[FromTo[, 1] == curID, 2]
			id.vec <- c(id.vec, curID)
		}
		id.vec[!is.na(id.vec)]
	})

# Creates an expanded grid of units for pairwise distance calculation
	full.grid <- expand.grid(units, units)
	lower.tri.mat <- matrix(1:(length(units)^2), length(units), length(units))
	lower.tris <- which(lower.tri(lower.tri.mat))
	lower.grid <- full.grid[lower.tris, ]
	lower.grid$Dist = NA

# For loop to assign distance values to a distance matrix
	move.dists <- matrix(0, length(units), length(units))
	for (i in 1:nrow(lower.grid)) {
		s1 <- downstreams[[lower.grid[i, 1]]]
		s2 <- downstreams[[lower.grid[i, 2]]]
		counts <- tapply(c(s1, s2), c(s1, s2), length)
		singles <- as.numeric(names(counts)[which(counts == 1)])
		singles <- unique(c(singles, lower.grid[i, 2]))
		singles <- singles[!singles == lower.grid[i, 1]]
		lower.grid[i, 3] <- length(singles)
	}

	move.dists[lower.tri(move.dists)] <- lower.grid[, 3]
	move.dists <- as.matrix(as.dist(move.dists))

# Determines the upstreams
	upstreams <- lapply(units, function(x) {
		units[sapply(downstreams, function(y) x %in% y)]
	})

# Identifies the downstream most unit
	base.unit <- FromTo[FromTo[, 2] == 0 | is.na(FromTo[, 2]), 1]
	cur.unit <- base.unit

# Identifies first order units
	firsts <- which(sapply(upstreams, length) == 1)

# Calculates link magnitude of each unit
	link.mag <- sapply(upstreams, function(x) sum(x %in% firsts))

# Determines the "real" downstream units
# These are units which an individual could reach by going downstream first from
# a focal unit. Some of these require subsequent upstream movements. 
	downstreams <- lapply(upstreams, function(x) {
		units[!units %in% x]
	})

# Stores the baseline move map
	move.map <- as.matrix(move.dists)

# Applies upstream fragmentation
	if (frag.type > 1 & length(fragments) > 0) {
		for (i in 1:length(fragments)) {
			cols <- upstreams[[fragments[i]]]
			rows <- downstreams[[fragments[i]]]
			rowcol <- expand.grid(rows, cols)

			move.map[rowcol[, 1], rowcol[, 2]] <- 0
		}
	}

# Applies downstream fragmentation
	if (frag.type == 3 & length(fragments) > 0) {
		for (i in 1:length(fragments)) {
			cols <- upstreams[[fragments[i]]]
			rows <- downstreams[[fragments[i]]]
			rowcol <- expand.grid(rows, cols)

			move.map[rowcol[, 2], rowcol[, 1]] <- 0
		}
	}

# Removes units above the maximum dispersal distance
	move.map <- ifelse(move.map > maxdist, 0, move.map)

# Converts the map to movement probabilities (these will be relativized during
# analyses)
	move.prob <- ifelse(move.map > 0, 1/(move.map^2), 0)
	move.prob <- t(apply(move.prob, 1, function(x) x / sum(x)))

# Removes the NA values from the FromTo object
	f2 <- as.matrix(na.omit(FromTo))
	f2 <- f2[order(f2[, 1], f2[, 2]), ]

# Creates a vector of edges for network construction
	f2vec <- vector("numeric", 0)
		for (i in 1:nrow(f2)) {
		vals <- unlist(f2[i, 1:2])
		f2vec <- c(f2vec, vals)
	}

# Constructs a network graph object
	net <- graph(edges = f2vec, n = 24, directed = F)

 # Identifies fragments and types in the network
 	frag.col <- rep("black", nrow(f2))
 	frag.lty <- rep(1, nrow(f2))
 	if (frag.type > 1 & length(fragments) > 0) {
 		frag.edges <- which(f2[, 1] %in% fragments)
 		frag.col[fragments] <- ifelse(frag.type == 2, "grey", "white")
		frag.lty[fragments] <- ifelse(frag.type == 2, 3, 0)
	}

# Identifies disturbed and/or degraded patches
 	disturb.shape <- rep("square", length(V(net)))
 	degraded.col <- rep("white", length(V(net)))

 	if (length(disturbed) > 0) {
 		vertex.disturb <- which(V(net) %in% disturbed)
		disturb.shape[vertex.disturb] <- "circle"
 	}

 	if (length(degraded) > 0) {
 		vertex.degraded <- which(V(net) %in% degraded)
 		degraded.col[vertex.degraded] <- "gray70"
 	}

# Adds attributes to the network
 	V(net)$link.mag <- link.mag
 	V(net)$disturb.shape <- disturb.shape
 	V(net)$degraded.col <- degraded.col
 	E(net)$frag.col <- frag.col 
 	E(net)$frag.lty <- frag.lty

# Stores a layout to provide stable plotting
	lay <- layout_with_fr(net)

# Outputs the result
	out <- list(move.prob = move.prob, network = net,
	    dists = as.dist(move.dists))
	out

}





metaPop <- function(N, r, K.mat, G, move.rate, move.map) {

# Creates empty N matrices for pre-migration and post-migration
	N.pre <- matrix(0, G, ncol(move.map))
	N.post <- N.pre
	N.pre[1, ] <- N

# Creates an empty list to store individuals who moved and their destinations
	move.list <- vector("list", G)

# Creates an empty list to store genealogies
	genealogy <- move.list

# Creates a matrix to store counts of females
	females <- matrix(0, G, ncol(N.pre))

# Nested for-loop to conduct metapopulation dynamics
	for (i in 1:(G - 1)) {
		cur.N <- N.pre[i, ]

		move.hist <- matrix(0, 0, 4)
		colnames(move.hist) <- c("Generation", "Row", "Origin", "Destination")

		move.subs <- vector("list", ncol(move.map))

	# Creates the movement history matrix for Generation i
		for (j in 1:ncol(N.pre)) {

		# Stores the initial starting N
			Nij <- N.pre[i, j]

		# Stores movers (emigrants)
			if (Nij > 0) {
				movers <- (1:Nij)[rbinom(Nij, 1, move.rate) == 1]
			} else {
				movers <- vector("numeric", 0)
			}

		# Handles the case of moves
			if (length(movers) > 1) {

			# Determines destination
				destinations <- sample(25, size = length(movers),
				    replace = TRUE, prob = move.map[j, ])
			
			# Stores a mini matrix of move histories for Nij
				move.hist <- cbind("Generation" = i, "Row" = movers,
				"Origin" = j, "Destination" = destinations)

			} else {

			# Handles the case of no moves and makes an empty matrix
				move.hist <- matrix(0, 0, 4)
				colnames(move.hist) <- c("Generation", "Row", "Origin",
				    "Destination")
			}

		# Stores the moves in the sub list
			move.subs[[j]] <- move.hist
		}

	# Combines the movement history
		move.subs <- do.call("rbind", move.subs)

	# Stores the movement history matrix for posterity
		move.list[[i]] <- move.subs

	# Calculates subtractions and additions from each subpopulation
	# Number of emigrants
		emigrants <- -tapply(move.subs[, 3], factor(move.subs[, 3],
			levels = 1:ncol(move.map)), length)
		emigrants <- ifelse(is.na(emigrants), 0, emigrants)

	# Number of immigrants
		immigrants <- tapply(move.subs[, 4], factor(move.subs[, 4],
			levels = 1:ncol(move.map)), length)
		immigrants <- ifelse(is.na(immigrants), 0, immigrants)

	# Net change
		nets <- emigrants + immigrants

	# Stores the post-migration population sizes
		demes <- N.pre[i, ] + nets
		N.post[i, ] <- demes

	# Stores male and female individuals
		mf.list <- lapply(demes, rbinom, size = 1, prob = 0.5)

	# Calculates the number of males and females in each subpopulation
	# If this ends up being 0 for either sex, no reproduction occurs in that
	# subpopulation
		mf.count <- t(sapply(mf.list, function(x) tapply(x, factor(x,
		    levels = 0:1), length)))
		mf.count <- ifelse(is.na(mf.count), 0, mf.count)

	# Stores the number of females
		females[i, ] <- mf.count[, 2]

	# Calculates population sizes in the next generation
		for (j in 1:ncol(N.pre)) {
			N.pre[(i + 1), j] <- floor(popGrow(r, n = N.post[i, j],
			    k = K.mat[i, j], mfc = mf.count[j, ]))
		}

	# Stores a genealogy list for this generation
		genealogy.sub <- vector("list", ncol(N.pre))

	# Assigns parentage history
		for (j in 1:ncol(N.pre)) {
			sexes <- mf.list[[j]]
			parent.N <- N.post[i, j]
			offspring.N <- N.pre[(i + 1), j]
			if (offspring.N > 0) {
				mom <- sample(parent.N, size = offspring.N, prob = sexes == 1,
				    replace = TRUE)
				dad <- sample(parent.N, size = offspring.N, prob = sexes == 0,
				    replace = TRUE)
				genealogy.sub[[j]] <- cbind(Generation = i + 1,
				    Subpopulation = j, Row = 1:offspring.N, mom, dad)
				} else {

				}

			}

	# Combines the genealogy list into a single matrix and stores it
		genealogy[[i]] <- do.call("rbind", genealogy.sub)

	}

# Combines the move list and genealogy list into single objects
	move.list <- do.call("rbind", move.list)
	genealogy <- do.call("rbind", genealogy)

# Stores the output
	output <- list("N.pre" = N.pre, "N.post" = N.post, move.list = move.list,
		genealogy = genealogy, females = females)

# Output
	output

}





genoType <- function(N, p.freq) {
# A function to create genotypes to start simulations

# Args
#	N: A vector of population sizes
#	p.freq: A matrix of frequency of the p allele (assumes biallelic loci)
#			where columns are populations and rows are loci

# Returns a list of arrays where:
#	row: individual
#	column: locus
#	slice: chromosome

# Stores a blank list object for each subpopulation
	geno.list <- vector("list", length(N))

# For loop to populate the genotype list
	for (i in 1:length(N)) {

	# Stores the p frequencies for the population
		pqi <- p.freq[, i]

	# Stores the population size
		Ni <- N[i]

	# Stops the process if Ni = 0

		if (Ni > 0) {
		# Creates the array for the subpopulation
			geno.i <- asplit(sapply(pqi, function(x) {
				replicate(Ni, sample(2, size = 2, replace = TRUE,
				    prob = c(x, 1 - x)))
				}, simplify = "array"), 2)

		# Stores the array
			geno.list[[i]] <- geno.i
		} else {

		}
	}

# Handles the output
	geno.list

}


	






geneFlow <- function(pop.hist, genotypes, Step = 10, Mutation = 2.5E-8,
	updateProgress = NULL) {
# A function to simulate reproduction and gene flow on top of an existing
# metapolulation simulation.
#
# Args:
#	pop.hist: The population history from function metaPop
#	genotypes: An input list of starting genotypes. These can be randomly
# 		created using function genoType or user specified.
#	

# Stores the step vector
	if (length(Step) == 1) {
		step.vec <- sort(unique(c(seq(0, nrow(pop.hist$N.pre), by = Step),
			nrow(pop.hist$N.pre))))
		step.vec <- step.vec[step.vec > 0]
	} else {
		step.vec <- Step
	}

# Stores a list of genotype histories through time
	geno.list <- vector("list", length(step.vec))
	names(geno.list) <- step.vec

# Extracts the movement history
	move.hist <- pop.hist$move.list
	genealogy <- pop.hist$genealogy

# Sets an initial genotype from the input genotype
	genotypes.i <- genotypes
	geno.list[[1]] <- genotypes.i

# For loop to conduct simulation
	for (i in 1:(nrow(pop.hist$N.pre) - 1)) {

	# Extracts pertinent move histories and genealogies
    	move.i <- move.hist[move.hist[, 1] == i, ]
		genealogy.i <- genealogy[genealogy[, 1] == (i + 1), ]

	# Moves individuals from the subpopulation of their birth to their final
	# destination
		if (nrow(move.i) > 0) {

	# First extracts the migrants
			migrant.list <- vector("list", nrow(move.i))
			for (j in 1:nrow(move.i)) {
				migrant.list[[j]] <- genotypes.i[[move.i[j, 3]]][[move.i[j, 2]]]
			}
			
		# Next removes these individuals from their source populations 
		# (bulk process)
			for (j in 1:length(genotypes.i)) {
				move.ij <- move.i[move.i[, 3] == j, ]
				if (nrow(move.ij) > 0) {
					genotypes.i[[j]] <- genotypes.i[[j]][-move.ij[, 2]]
				}
			}

		# Finally puts the migrants into their end destinations
			for (j in 1:length(migrant.list)) {
				genotypes.i[[move.i[j, 4]]] <- c(genotypes.i[[move.i[j, 4]]],
					migrant.list[j])
			}

		}

	# Combines lists into arrays
		genotypes.i <- lapply(genotypes.i, function(x) {
			if (!is.null(x)) {
				nx <- simplify2array(x)
				aperm(nx, c(3, 2, 1))
			} 
		})

		genotypes.i2 <- vector("list", length(genotypes.i))

	# Gene passage from parent to offspring
		for (j in 1:length(genotypes.i)) {
			gene.ij <- genealogy.i[genealogy.i[, 2] == j, , drop = FALSE]
			if (nrow(gene.ij) > 0) {
				new.genotypes <- apply(gene.ij, 1, function(x) {
					apply(genotypes.i[[j]][x[c(4, 5)], , ], c(1, 2), sample,
						size = 1)
					}, simplify = FALSE)

			# Converts the new gene list to an array to simplify mutation
			# calculation
				new.gene.array <- array(unlist(new.genotypes), dim = c(2,
					max(sapply(new.genotypes, ncol)), length(new.genotypes)))

			# Rescales the scores for different alleles
				new.gene.array <- 2 * (new.gene.array - 1.5)

			# Determines the mutations
				mutations <- 1 - (2 * rbinom(length(new.gene.array), 1,
				    Mutation))

			# Mutates the alleles
				new.gene.array <- new.gene.array * mutations

			# Re-rescales the allele values
				new.gene.array <- (new.gene.array/2) + 1.5

			# Relists the new genotypes
				new.genotypes <- apply(new.gene.array, 3, as.matrix,
				    simplify = "list")

			# Stores the new genotypes
				genotypes.i2[[j]] <- new.genotypes
			}
		}

	
	# Advances the genotypes
		genotypes.i <- genotypes.i2


#####

	# Inserts the new genotypes into the genotype history list
		if ((i + 1) %in% step.vec) {
			geno.list[[which(step.vec == (i + 1))]] <- genotypes.i2
		}
		

	# Progress update
		if (is.function(updateProgress)) {
			text <- paste(round((i/nrow(pop.hist$N.pre)) * 100, 0), "%",
				sep = "")
			updateProgress(detail = text)
		}

#####

	}

# Converts the deep list structure to arrays
	genotype.history <- lapply(geno.list, function(x) lapply(x,
	    simplify2array))

# Handles the output
	genotype.history

}








popCrunch <- function(pop) {
# A function to crunch several summary statistics from a population from a
# simulation iteration
	p.count <- apply(pop, c(3, 2), function(x) sum(x == 1))
	p <- apply(p.count, 2, sum)
	p <- p/(nrow(p.count) * 2)
	q <- 1 - p
	eh <- 1 - (p^2 + q^2)
	pp <- apply(p.count, 2, function(x) sum (x == 2))
	pq <- apply(p.count, 2, function(x) sum(x == 1))
	qq <- apply(p.count, 2, function(x) sum(x == 0))
	oh <- pq/(pp + pq + qq)

	cbind(p, q, eh, oh, pp, pq, qq)

}






FST <- function(p) {
# A function to calculate FST given a vector of p values (allele 1) from
# biallelic genetic data
	q <- 1 - p
	pq <- cbind(p, q)
	Hs <- mean(1 - apply(pq^2, 1, sum))
	Ht <- 1 - sum(apply(pq, 2, mean)^2)
	fst <- (Ht - Hs)/Ht
	ifelse(is.na(fst), 0, fst)
}




FSTPair <- function(p) {
# A function to calculate pairwise F[ST]. This is a wrapper for FST
# Creates an expanded grid for 
	test.grid <- expand.grid(1:length(p), 1:length(p))
	pair.fst.raw <- apply(test.grid, 1, function(x) {
		FST(p[unlist(x)])
		})
	pair.fst <- matrix(pair.fst.raw, length(p), length(p))
	pair.fst
}




OmniFSTPair <- function(sim, dist.out = TRUE, keep = NULL) {
# A function to calculate mean pairwise F[ST] from an iteration of a simulation
#
# Args:
#	sim: The simulation iteration. This should be a list of population genotypes
#	dist.out: Logical. Should the output be returned as a distance object?
#	keep: Either NULL (keep all) or a numeric vector of the population numbers
# 		to retain. These should match the numbers from the movement matrix.
#
# Returns either:
#	A pairwise F[ST] dissimilarity object
#	A pairwise F[ST] matrix

# Determines which populations ended with individuals
	sim.pops <- c(1:length(sim))[!sapply(sim, is.null)]

# Calculates summary genetic statistics for each population
	sim.sum <- lapply(sim, function(x) {
		if (!is.null(dim(x))) {
			popCrunch(x)
		}
	})

# Simplifies the summary statistics to an array. NOTE that this drops empty
# populations
	sim.array <- simplify2array(sim.sum[sim.pops])

# Creates a list of pairwise F[ST] for each locus
	sim.fst.list <- apply(sim.array[, 1, ], 1, FSTPair, simplify = FALSE)

# Calculates average F[ST] for each population pair across all loci
	sim.pair.fst <- Reduce("+", sim.fst.list)/length(sim.fst.list)

# Sets row and column names to existing populations (handles empty pops)
	rownames(sim.pair.fst) <- as.character(sim.pops)
	colnames(sim.pair.fst) <- as.character(sim.pops)

# Assesses for whether a subset of the results should be returned
	if (!is.null(keep)) {
		sim.pair.fst <- sim.pair.fst[sim.pops %in% keep, sim.pops %in% keep]
	}

# Handles the conditional output: if dist.out == TRUE, this will return a
# dissimilarity object. Otherwise, the pairwise F[ST] is returned as a matrix
	if (dist.out == TRUE) {
		as.dist(sim.pair.fst)
	} else {
		sim.pair.fst
	}

}


















#______________________________________________________________________________#
#                                                                              #
#==============================================================================#
#   END OF SCRIPT                                                              #
#==============================================================================#