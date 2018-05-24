#vertices are mutants and genes

#edge between mutant and gene if HIGH priority mutation observed

#read in gene sets

#filter down to just priority columns
#add column with mutant name
arguments <- commandArgs(trailingOnly=TRUE)
fileLoc <- arguments[1]
linesToRemove <- arguments[2:length(arguments)]
if (length(linesToRemove) == 1) {
	if (linesToRemove == "X") {
		linesToRemove <- NULL
	}
}

suppressPackageStartupMessages(library(igraph,quietly=TRUE))

files <- dir(fileLoc)
filenames <- files[grep("genes.txt",files)]

output <- NULL

# Iterate over files, get per-locus mutation info

for (eachNum in 1:length(filenames)){
	currentMut <- read.delim(paste0(fileLoc,"/",filenames[eachNum]),header=TRUE)
	currentName <- strsplit(filenames[eachNum],"\\.")[[1]][2]
	currentMut$mutant <- currentName

	if(sum(grepl("impact_HIGH",names(currentMut)))==0) {
		currentMut$variants_impact_HIGH <- 0
	} else if (sum(grepl("impact_MODERATE",names(currentMut)))==0) {
		currentMut$variants_impact_moderate <- 0
	} else if (sum(grepl("impact_LOW",names(currentMut)))==0) {
		currentMut$variants_impact_LOW <- 0
	}
	
	currentMut <- currentMut[c("mutant","GeneName","GeneId","TranscriptId","variants_impact_HIGH","variants_impact_MODERATE","variants_impact_LOW")]
	if (is.null(output)){
		output <- currentMut
	} else {
		output <- rbind(output,currentMut)
	}
}


# Remove duplicate rows
output <- output[!duplicated(output[c(1,2)]),]

# Find which loci have mutations
booleanHigh <- output$variants_impact_HIGH>0
booleanModerate <- output$variants_impact_HIGH>0 | output$variants_impact_MODERATE>0


# Generate subset tables
for (eachNum in 1:2){
	currentSet <- list(booleanModerate,booleanHigh)[[eachNum]]
	currentName <- c("Moderate","High")[eachNum]
	
	# Set up output complete graph
	currentGraph <- output[currentSet,]
	currentGraphName <- paste0(fileLoc,"/",currentName,".full.txt")
	
	# Set up graph object to manipulate
	graphObj <- graph_from_data_frame(currentGraph[c(1,3,2,4,5,6,7)])
	# Get genes shared between lines from complete graph
	multDegree <- degree(graphObj)
	multDegreeComp <- multDegree>1
	vsToCut <- V(graphObj)[multDegreeComp]
	currentGraphShared <- induced_subgraph(graphObj,vsToCut)
	currentGraphSharedFrame <- as_data_frame(currentGraphShared)
	currentGraphSharedName <- paste0(fileLoc,"/",currentName,".shared.txt")
	
	if (is.null(linesToRemove) == FALSE) {
		# Remove bad lines from complete graph
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
	currentGraph <- as_data_frame(graphObj)
	
	# Write output
	write.table(currentGraph,currentGraphName,sep="\t",quote=FALSE,row.names=FALSE)
	write.table(currentGraphSharedFrame,currentGraphSharedName,sep="\t",quote=FALSE,row.names=FALSE)
}

#matrixHigh <- as.matrix(outputHigh[c(1,2,3)])
#matrixModerate <- as.matrix(outputModerate[c(1,2,3)])
#write.table(matrixHigh,paste0(fileLoc,"/matrixHigh.txt"),quote=FALSE,row.names=FALSE,sep="\t")
#write.table(matrixModerate,paste0(fileLoc,"/matrixModerate.txt"),quote=FALSE,row.names=FALSE,sep="\t")
