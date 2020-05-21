#vertices are mutants and genes

#edge between mutant and gene if MODERATE or HIGH priority mutation observed
options(stringsAsFactors=FALSE)
arguments <- commandArgs(trailingOnly=TRUE)
fileLoc <- arguments[1]
siblingThreshold <- as.numeric(arguments[2])
sibDetectThreshold <- as.numeric(arguments[3])
linesToRemove <- arguments[4:length(arguments)]
if (length(linesToRemove) == 1) {
	if (linesToRemove == "X") {
		linesToRemove <- NULL
	}
}

suppressPackageStartupMessages(library(igraph,quietly=TRUE))


removeSingles <- function(currentObj,vsToKeep=NULL) {
	# currentObj = igraph graph object
	# vsToKeep = character list of vertex names you want to keep
	multDegree <- degree(currentObj)
	multDegreeComp <- multDegree>1
	multDegreeComp[names(multDegreeComp) %in% as.character(vsToKeep)] <- TRUE
	vsToCut <- V(currentObj)[multDegreeComp]
	reducedGraph <- induced_subgraph(currentObj,vsToCut)
	return(reducedGraph)
}

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
vcfGraphCut <- removeSingles(vcfGraphObj)

vcfGraphFrame <- as_data_frame(vcfGraphCut)
vcfGraphName <- paste0(fileLoc,"/snpGraph.txt")
write.table(vcfGraphFrame,file=vcfGraphName,sep="\t",quote=FALSE)

# This gives me a graph of mutant-to-mutant connections.

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
putativeSibs <- mut2mut[(mut2mut$fromRatio > siblingThreshold & mut2mut$fromSnps > sibDetectThreshold) | (mut2mut$toRatio > siblingThreshold & mut2mut$toSnps > sibDetectThreshold),]
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

	# if not, check currentGraphCollapsed for a line that goes from that gene to that sib group, and remove it

	# Set up graph objects to manipulate
	graphObj <- graph_from_data_frame(currentGraph[c(1,3,2,4,5,6,7)])
	graphObjCollapsed <- graph_from_data_frame(currentGraphCollapsed[c(1,3,2,4,5,6,7)])
	
	# For collapsed graph, remove edges from genes to sib group that aren't present in all of the sibs in that group
	groupNames <- names(V(graphObjCollapsed))
	
	presentGroups <- sibMap$Group[sibMap$Group %in% groupNames]
	
	for (eachGroup in presentGroups) {
		# Get vertex ID for the current group
		currentGroup <- V(graphObjCollapsed)[groupNames==eachGroup]
		
		# Get all the genes in that sib group
		groupNeighbors <- neighbors(graphObjCollapsed,as.character(eachGroup),mode="all")
		# For each gene, check that all of the sibs are in its neighbors
		listOfSibs <- sibMap$Line[sibMap$Group==eachGroup]
		for (eachSnp in groupNeighbors) {
			geneNeighbors <- neighbors(graphObj,eachSnp,mode="all")
			if (sum(names(geneNeighbors) %in% listOfSibs) < length(listOfSibs)) {
				graphObjCollapsed <- delete_edges(graphObjCollapsed,E(graphObjCollapsed,c(eachSnp,currentGroup),directed=FALSE))
			}
		}
	}

		

	# Get genes shared between lines from complete graph	
	currentGraphShared <- removeSingles(graphObj)
	currentGraphSharedCollapsed <- removeSingles(graphObjCollapsed)
	
	# Make data frames and write to file
	currentGraphSharedFrame <- as_data_frame(currentGraphShared)
	currentGraphSharedName <- paste0(fileLoc,"/",currentName,".shared.txt")
	currentGraphSharedCollapsedFrame <- as_data_frame(currentGraphSharedCollapsed)
	currentGraphCollapsedName <- paste0(fileLoc,"/",currentName,".collapsed.txt")
	currentGraphSharedCollapsedName <- paste0(fileLoc,"/",currentName,".shared.collapsed.txt")
	
	# Generate sib group networks and write to file
	for (eachComponent in 1:sibComponents$no) {
		currentComponentNames <- names(sibComponents$membership[sibComponents$membership==eachComponent])
		
		# It's possible a putative sibling has no High impact SNPs
		# Need to add them back in
		missingLines <- as.character(perLineSnps$Line[! perLineSnps$Line %in% names(V(graphObj))])
		graphObjWithMissing <- add_vertices(graphObj,length(missingLines),name=missingLines)
		
		# subset graph for just these vertices and their first neighbors
		egoSets <- make_ego_graph(graphObjWithMissing,order=1,nodes=currentComponentNames,mode="all",mindist=0)
		egoVertices <- NULL
		for (eachEgo in 1:length(egoSets)) {
			egoVertices <- c(egoVertices,names(V(egoSets[[eachEgo]])))
		}
		egoVertices <- unique(egoVertices)
		egoGraph <- induced_subgraph(graphObjWithMissing,egoVertices)
		egoGraphCut <- removeSingles(egoGraph,as.character(perLineSnps$Line))
		
		egoGraphName <- paste0(fileLoc,"/","siblingset.",eachComponent,".",currentName,".txt")
		egoGraphCutName <- paste0(fileLoc,"/","siblingset.",eachComponent,".",currentName,".shared.txt")
		
		egoGraphFrame <- as_data_frame(egoGraph)
		egoGraphCutFrame <- as_data_frame(egoGraphCut)
		
		write.table(egoGraphFrame,file=egoGraphName,sep="\t",quote=FALSE,row.names=FALSE)
		write.table(egoGraphCutFrame,file=egoGraphCutName,sep="\t",quote=FALSE,row.names=FALSE)
	}
	
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
		currentGraphCutShared <- removeSingles(cutObj)
		
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
