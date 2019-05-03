#vertices are mutants and genes

#edge between mutant and gene if MODERATE or HIGH priority mutation observed

arguments <- commandArgs(trailingOnly=TRUE)
fileLoc <- arguments[1]
siblingThreshold <- arguments[2]
linesToRemove <- arguments[3:length(arguments)]
if (length(linesToRemove) == 1) {
	if (linesToRemove == "X") {
		linesToRemove <- NULL
	}
}

suppressPackageStartupMessages(library(igraph,quietly=TRUE))

# Get filenames
files <- dir(fileLoc)
filenames <- files[grep("genes.txt",files)]
annFilenames <- files[grep(".ann",files)]

# Set up output objects
output <- NULL
vcfOutput <- NULL

# Object to store SNP information per-line
perLineSnps <- data.frame(
	Line=character(),
	SNPs=integer()
)

# Iterate over files, get per-locus mutation info
for (eachNum in 1:length(filenames)){
	currentMut <- read.delim(paste0(fileLoc,"/",filenames[eachNum]),header=TRUE, sep="\t")
	currentName <- strsplit(filenames[eachNum],"\\.")[[1]][2]
	currentMut$mutant <- currentName

	if(sum(grepl("impact_HIGH",names(currentMut)))==0) {
		currentMut$variants_impact_HIGH <- 0
	}
	if (sum(grepl("impact_MODERATE",names(currentMut)))==0) {
		currentMut$variants_impact_MODERATE <- 0
	}
	if (sum(grepl("impact_LOW",names(currentMut)))==0) {
		currentMut$variants_impact_LOW <- 0
	}
	
	currentMut <- currentMut[c("mutant","GeneName","GeneId","TranscriptId","variants_impact_HIGH","variants_impact_MODERATE","variants_impact_LOW")]
	if (is.null(output)){
		output <- currentMut
	} else {
		output <- rbind(output,currentMut)
	}
	
	# Get annotation file for current set
	currentAnn <- read.delim(paste0(fileLoc,"/",annFilenames[eachNum]),header=TRUE,sep="\t",comment.char="#")
	currentAnn$mutant <- currentName
	
	if (is.null(vcfOutput)) {
		vcfOutput <- currentAnn
	} else {
		vcfOutput <- rbind(vcfOutput,currentAnn)
	}
	perLineSnps <- rbind(perLineSnps,data.frame(Line=currentName,SNPs=dim(currentAnn)[1]))
}


# Remove duplicate rows
output <- output[!duplicated(output[c(1,2)]),]

# Find which loci have mutations
booleanHigh <- output$variants_impact_HIGH>0
booleanModerate <- output$variants_impact_HIGH>0 | output$variants_impact_MODERATE>0

# Make object to store sibling information

# Make SNP graph
vcfGraphInput <- vcfOutput
vcfGraphInput$Identifier <- paste(vcfGraphInput$CHROM, vcfGraphInput$X0, vcfGraphInput$REF, vcfGraphInput$ALT,sep=".")
vcfGraphObj <- graph_from_data_frame(vcfGraphInput[c(11,12,1,2,3,4,5,6,7,8,9,10)])
multDegree <- degree(vcfGraphObj)
multDegreeComp <- multDegree > 1
vsToCut <- V(vcfGraphObj)[multDegreeComp]
vcfGraphCut <- induced_subgraph(vcfGraphObj,vsToCut)

vcfGraphFrame <- as_data_frame(vcfGraphCut)
vcfGraphName <- paste0(fileLoc,"/snpGraph.txt")
write.table(vcfGraphFrame,file=vcfGraphName,sep="\t",quote=FALSE)

# This gives me a graph of mutant-to-mutant connections.
# There's far more random overlap than I expected but honestly that's dumb of me, that's the whole POINT of this project.
# I need a metric for distinguishing sibs from non-sibs.
# Carly's idea: compare number of overlaps to number of SNPs
# I can think of a few ways to quantify these comparisons
# 1: Do some up-front work pre-generating some thresholds based on Arabidopsis
# 2: Do some on-the-fly simulation work for each sib pair

# Notes: SNP density should be similar between sibs.
# How similar?
# What's the simulated 

# multiple criteria?


# Get the vertex IDs for vertices that represent mutant lines, not genes
mutantLineVertices <- V(vcfGraphCut)[V(vcfGraphCut)$name %in% perLineSnps$Line]

# Look up shortest path from vertex to vertex
mutantLineDistances <- distances(vcfGraphCut,v=mutantLineVertices,to=mutantLineVertices)

# Only consider comparisons where shortest path is two edges
mut2mut <- data.frame(from=character(),to=character())

# Iterate over each set of shared SNPs, add data to output object
for (eachRow in 1:nrow(mutantLineDistances)) {
	for (eachCol in 1:ncol(mutantLineDistances)) {
		if (mutantLineDistances[eachRow,eachCol] == 2) {
			# store the value
			mut2mut <- rbind(mut2mut,data.frame(from=rownames(mutantLineDistances)[eachRow],
				to=colnames(mutantLineDistances)[eachCol]))
		}
	}
}

# Create the duplicate column
mut2mut$duplicate <- FALSE

# Need to remove duplicated rows from this table
for (eachRow in 1:nrow(mut2mut)) {
	currentFrom <- mut2mut$from[eachRow]
	currentTo <- mut2mut$to[eachRow]
	
	# Need to find reciprocal entries
	# I check from the current position to the end
	# If I see a duplicate in that range, I mark the current one false
	# If I don't see a duplicate in that range, the current one must be a duplicate of an earlier one
	# Leave it false
	fromBool <- mut2mut$from[eachRow:nrow(mut2mut)] %in% currentTo
	toBool <- mut2mut$to[eachRow:nrow(mut2mut)] %in% currentFrom
	if ( any(fromBool & toBool)) {
		mut2mut$duplicate[eachRow] <- TRUE
	}
}

# Subset to remove duplicate rows
mut2mut <- mut2mut[mut2mut$duplicate,c(1,2)]

# Set up columns for counting SNPs, shared and not
mut2mut$sharedSnps <- 0
mut2mut$fromSnps <- 0
mut2mut$toSnps <- 0

# Count how many SNPs are shared for each comparison
for (eachRow in 1:nrow(mut2mut) ){
	# Figure out the current set of mutant lines
	partner1 <- mut2mut$from[eachRow]
	partner2 <- mut2mut$to[eachRow]
	
	# Find neighboring SNPs for each mutant line
	p1Neighbors <- neighbors(vcfGraphCut,as.character(partner1),mode="all")
	p2Neighbors <- neighbors(vcfGraphCut,as.character(partner2),mode="all")
	
	# Count shared SNPs between mutant lines
	linkCounter <- sum(p1Neighbors %in% p2Neighbors)

	# Store data
	mut2mut$sharedSnps[eachRow] <- linkCounter
	mut2mut$fromSnps[eachRow] <- perLineSnps$SNPs[perLineSnps$Line %in% partner1]
	mut2mut$toSnps[eachRow] <- perLineSnps$SNPs[perLineSnps$Line %in% partner2]
}

# Generate ratio of shared to all SNPs
mut2mut$fromRatio <- mut2mut$sharedSnps / mut2mut$fromSnps
mut2mut$toRatio <- mut2mut$sharedSnps / mut2mut$toSnps
putativeSibs <- mut2mut[mut2mut$fromRatio > siblingThreshold | mut2mut$toRatio > siblingThreshold,]
# I struggled with whether to make the above line an AND or an OR condition.
# My thinking with AND initially was I didn't want spurious sibling inferences so I wanted to make sure it was a good call
# But then I thought, maybe one of the sibling partners wasn't sequenced well.
# So I switched to OR

# Identify sibling sets
sibGraph <- graph_from_data_frame(putativeSibs,directed=FALSE)
sibComponents <- components(sibGraph)

# Object for component to sibling mappings
sibMap <- data.frame(Line=names(sibComponents$membership),Group=paste0("group.",sibComponents$membership))

# Generate subset tables
for (eachNum in 1:2){
	currentSet <- list(booleanModerate,booleanHigh)[[eachNum]]
	currentName <- c("Moderate+High","High")[eachNum]
	
	# Set up output complete graph
	currentGraph <- output[currentSet,]
	currentGraphName <- paste0(fileLoc,"/",currentName,".full.txt")
	
	# Copy graph so we can change names
	currentGraphCollapsed <- currentGraph

	# Iterate over vertices
	for (eachRow in 1:nrow(currentGraphCollapsed)) {
		# Check if this mutant is in a sibling group
		if (currentGraphCollapsed$mutant[eachRow] %in% sibMap$Line) {
			# Replace mutant name with sibling group ID
			currentGraphCollapsed$mutant[eachRow] <- sibMap$Group[sibMap$Line==currentGraph$mutant[eachRow]]
		}
	}	
		
	# Set up graph objects to manipulate
	graphObj <- graph_from_data_frame(currentGraph[c(1,3,2,4,5,6,7)])
	graphObjCollapsed <- graph_from_data_frame(currentGraphCollapsed[c(1,3,2,4,5,6,7)])
	
	# Get genes shared between lines from complete graph	
	for (tempNum in 1:2){
		currentObj <- list(graphObj,graphObjCollapsed)[[tempNum]]
		multDegree <- degree(currentObj)
		multDegreeComp <- multDegree>1
		vsToCut <- V(currentObj)[multDegreeComp]
		
		# Use eachNum to decide which variable to write to
		if (tempNum == 1) {
			currentGraphShared <- induced_subgraph(currentObj,vsToCut)
		} else {
			currentGraphSharedCollapsed <- induced_subgraph(currentObj,vsToCut)
		}
	}
	
	# Make data frames and write to file
	currentGraphSharedFrame <- as_data_frame(currentGraphShared)
	currentGraphSharedName <- paste0(fileLoc,"/",currentName,".shared.txt")
	currentGraphSharedCollapsedFrame <- as_data_frame(currentGraphSharedCollapsed)
	currentGraphCollapsedName <- paste0(fileLoc,"/",currentName,".collapsed.txt")
	currentGraphSharedCollapsedName <- paste0(fileLoc,"/",currentName,".shared.collapsed.txt")
	
	# Remove bad lines from complete graph
	# TODO: integrate collapsed graphs into this section
	if (is.null(linesToRemove) == FALSE) {
		cutObj <- graphObj

		currentGraphCut <- as_data_frame(cutObj,what="both")
		currentGraphCutName <- paste0(fileLoc,"/",currentName,".linesRemoved.full.txt")
		
		for (eachLine in linesToRemove) {
			if (eachLine %in% vertex_attr(cutObj)$name) {
				cutObj <- cutObj - neighbors(cutObj,eachLine)
			}
		}
		
		# Get genes shared between lines from graph with bad lines removed
		cutDegree <- degree(cutObj)
		cutDegreeComp <- cutDegree>1
		vsToCut <- V(cutObj)[cutDegreeComp]
		currentGraphCutShared <- induced_subgraph(cutObj,vsToCut)
		currentGraphCutSharedFrame <- as_data_frame(currentGraphCutShared)
		currentGraphCutSharedName <- paste0(fileLoc,"/",currentName,".linesRemoved.shared.txt")		
		
		write.table(currentGraphCut$edges,currentGraphCutName,sep="\t",quote=FALSE,row.names=FALSE)
		
		write.table(currentGraphCutSharedFrame,currentGraphCutSharedName,sep="\t",quote=FALSE,row.names=FALSE)
	
	}
	currentGraphFrame <- as_data_frame(graphObj)
	currentGraphCollapsedFrame <- as_data_frame(graphObjCollapsed)
	
	# Write collapsed output
	write.table(currentGraphCollapsedFrame,currentGraphCollapsedName,quote=FALSE,sep="\t",row.names=FALSE)
	write.table(currentGraphSharedCollapsedFrame,currentGraphSharedCollapsedName,sep="\t",quote=FALSE,row.names=FALSE)
	
	# Write output
	write.table(currentGraphFrame,currentGraphName,sep="\t",quote=FALSE,row.names=FALSE)
	write.table(currentGraphSharedFrame,currentGraphSharedName,sep="\t",quote=FALSE,row.names=FALSE)
}
