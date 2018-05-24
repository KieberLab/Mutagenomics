options(stringsAsFactors=FALSE)

arguments <- commandArgs(trailingOnly=TRUE)
hetFolder <- arguments[1]
hetPartial <- arguments[2]
hetFull <- arguments[3]
homFolder <- arguments[4]
homPartial <- arguments[5]
homFull <- arguments[6]

# Read in metadata and external data
metadata <- read.delim("gene_aliases_20140331.txt")
atgenData <- read.delim("atgenDataMean.txt")
atgenMetadata <- read.delim("atgen_metadata.txt")
atgenData$AGI.code <- toupper(atgenData$AGI.code)
chip <- read.delim("ARR10.targets.Schaller.txt")
KevinData <- read.delim("Kevin Gene Lists.txt")

# Construct filenames
hetPartialName <- paste0(hetFolder,"/",hetPartial)
hetFullName <- paste0(hetFolder,"/",hetFull)
homPartialName <- paste0(homFolder,"/",homPartial)
homFullName <- paste0(homFolder,"/",homFull)

# Read in data
hetPartialFile <- read.delim(hetPartialName)
hetFullFile <- read.delim(hetFullName)
homPartialFile <- read.delim(homPartialName)
homFullFile <- read.delim(homFullName)

# Make single gene lists
hetPartialGenes <- unique(c(hetPartialFile$from,hetPartialFile$to))
hetFullGenes <- unique(c(hetFullFile$from,hetFullFile$to))
homPartialGenes <- unique(c(homPartialFile$from,homPartialFile$to))
homFullGenes <- unique(c(homFullFile$from,homFullFile$to))

# Make output table
outputData <- data.frame(hetAll=hetFullGenes,homAll="N/A",hetShared="N/A",homShared="N/A")

# Populate output sheet
# homozygous unique overlap
for (each in 1:dim(outputData)[1]) {
	if ( outputData$hetAll[each] %in% homFullGenes) {
		outputData$homAll[each] <- outputData$hetAll[each]
	}
}

# heterozygous shared genes
for (each in 1:dim(outputData)[1]) {
	if ( outputData$hetAll[each] %in% hetPartialGenes) {
		outputData$hetShared[each] <- outputData$hetAll[each]
	}
}

# homozygous shared genes
for (each in 1:dim(outputData)[1]) {
	if ( outputData$hetAll[each] %in% homPartialGenes) {
		outputData$homShared[each] <- outputData$hetAll[each]
	}
}

# Match annotation data
annotationComp <- match(outputData$hetAll,metadata$locus_name)
annotationSymbol <- metadata$symbol[annotationComp]
annotationName <- metadata$full_name[annotationComp]
naSearch <- is.na(annotationSymbol)
annotationSymbol[naSearch] <- ""
annotationName[naSearch] <- ""

outputData$locus_symbol <- annotationSymbol
outputData$full_name <- annotationName

outputData$chip <- ""
outputData$Kevin <- ""
outputData$tissue <- ""
atgenNames <- names(atgenData)

for (each in 1:dim(outputData)[1]) {
	currentLocus <- outputData$hetAll[each]
	if (currentLocus %in% atgenData$AGI.code) {
		locusData <- atgenData[atgenData$AGI.code %in% currentLocus,]
		locusDataSorted <- locusData[c(-1,-2)][order(locusData[c(-1,-2)],decreasing=TRUE)]
		locusDatasets <- names(locusDataSorted)[1:3]
		locusTissue <- ""
		for (eachSet in locusDatasets) {
			currentTissue <- atgenMetadata$Tissue[atgenMetadata$Sample_ID==eachSet]
			if (locusTissue == "") {
				locusTissue <- currentTissue
			} else {
				locusTissue <- paste(locusTissue,currentTissue,sep="; ")
			}
		}
		#locusDataset <- atgenNames[2+which(locusData[c(-1,-2)] == max(locusData[c(-1,-2)]))]
		#outputData$tissue[each] <- atgenMetadata$Tissue[atgenMetadata$Sample_ID==locusDataset]
		outputData$tissue[each] <- locusTissue
	}
	if (currentLocus %in% chip$geneName) {
		outputData$chip[each] <- "BA-specific binding"
	}
	if (currentLocus %in% KevinData$Gene) {
		KevinCut <- KevinData[KevinData$Gene %in% currentLocus,]
		KevinName <- paste(KevinCut$Experiment,collapse=", ")
		outputData$Kevin[each] <- KevinName
	}
}

# Reorder data
outputData <- outputData[order(outputData$homShared),]
# Write to output
write.table(outputData,"outputData.txt",sep="\t",quote=FALSE,row.names=FALSE)