# This script removes background SNPs and also does a little metadata generation and filtering
# Expects: Rscript backgroundSNPscreen.R WTfile alleleThresholdHet alleleThresholdHom pipeline
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly=TRUE)
WTfile <- args[1]
alleleThresholdHet <- as.numeric(args[2])
alleleThresholdHom <- as.numeric(args[3])
pipeline <- args[4]

# Read in WT data
wt <- read.delim(WTfile,comment.char="#")
names(wt) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","otherinfo")
wtSnps <- wt[c("CHROM","POS","ID","REF","ALT","QUAL")]

# Get filenames for mutant files
filenames <- dir("filtered/")

mutFileSearch <- grep("mut",filenames)

mutFiles <- filenames[mutFileSearch]

# Filter characteristics


# Iterate over mutant files, remove background snps
for (eachName in mutFiles) {
	currentMut <- read.delim(paste0("filtered/",eachName),comment.char="#")
	outputNameHet <- paste0("mutsnps.het/het.",eachName,".txt")
	outputNameHom <- paste0("mutsnps.hom/hom.",eachName,".txt")
	outputTableHet <- NULL
	outputTableHom <- NULL
	names(currentMut) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","otherinfo")
	
	for (eachLine in 1:dim(currentMut)[1]) {
		currentLine <- currentMut[eachLine,]
		
		# Check for presence in WT background
		wtCut <- wt[currentLine$CHROM %in% wt$CHROM & currentLine$POS %in% wt$POS,]
		if (dim(wtCut)[1] == 0) {
			
			# check EMS mutation
			if ((currentLine$REF == "G" & currentLine$ALT == "A") | (currentLine$REF == "C" & currentLine$ALT == "T") ) {
				# check quality score
				if (currentLine$FILTER == "PASS") {
					# check allelic depth
					
					if (pipeline == "samtools") {
						infoSet <- strsplit(currentLine$otherinfo,":")[[1]]
						ADset <- as.numeric(strsplit(infoSet[5],",")[[1]])
						altRatio <- ADset[2] / sum(ADset)
						ADtest <- grepl("AD",currentLine$FORMAT)
						if (altRatio >= alleleThresholdHet & ADtest) {
							# write line to table, it's not a background SNP
							if (is.null(outputTableHet)) {
								outputTableHet <- currentLine
							} else {
								outputTableHet <- rbind(outputTableHet,currentLine)
							}
						}
						if (altRatio >= alleleThresholdHom & ADtest) {
							# write line to table, it's not a background SNP
							if (is.null(outputTableHom)) {
								outputTableHom <- currentLine
							} else {
								outputTableHom <- rbind(outputTableHom,currentLine)
							}
						}
					} else if (pipeline == "gatk") {
						ADtest <- grepl("AD",currentLine$FORMAT)
						infoSet <- strsplit(currentLine$otherinfo,":")[[1]]
						ADset <- as.numeric(strsplit(infoSet[2],",")[[1]])
						altRatio <- ADset[2] / sum(ADset)

						if (sum(ADset)==0 |  any(is.na(ADset))) {
							altRatio <- 0
						}
						
						if (altRatio >= alleleThresholdHet & ADtest) {
							# write line to table, it's not a background SNP
							if (is.null(outputTableHet)) {
								outputTableHet <- currentLine
							} else {
								outputTableHet <- rbind(outputTableHet,currentLine)
							}
						}
						if (altRatio >= alleleThresholdHom & ADtest) {
							# write line to table, it's not a background SNP
							if (is.null(outputTableHom)) {
								outputTableHom <- currentLine
							} else {
								outputTableHom <- rbind(outputTableHom,currentLine)
							}
						}
					}
				}
			}
		}
	}
	write.table(outputTableHom,outputNameHom,quote=FALSE,sep="\t",row.names=FALSE)
	write.table(outputTableHet,outputNameHet,quote=FALSE,sep="\t",row.names=FALSE)
}
