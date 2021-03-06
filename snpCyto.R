options(stringsAsFactors=FALSE)
library(RCy3)
ck <- read.delim("output.hom/Moderate.shared.txt")
ck$from <- as.character(ck$from) #important for making networks, source and target need to be character, not integers
ckCut <- ck[c(1,3)]
names(ckCut) <- c("source","target")
ckNet <- createNetworkFromDataFrames(edges=ckCut[c(1,2)])

nodeList <- getAllNodes(ckNet)
mutList <- nodeList[nodeList %in% ckCut$source]


setNodeColorBypass(mutList,new.colors="#ffff00",network=ckNet)
setNodeShapeDefault(new.shape="ELLIPSE")
setNodeShapeBypass(mutList,new.shapes="ROUND_RECTANGLE",network=ckNet)

