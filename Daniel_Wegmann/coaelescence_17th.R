plotTrees <- function(trees, popCol=c("black", "orange2", "purple")){
  maxHeight <- max(unlist(lapply(trees, function(x){ max(branching.times(x)) })))
  nCols <- ceiling(sqrt(length(trees)))
  nRows <- ceiling(length(trees) / nCols)
  par(mfrow = c(nRows, nCols), oma=c(0,4,0,0), las=1, xpd=NA)
  
  for(tr in 1:length(trees)){
    plotTree(trees[[tr]], direction="downwards", ylim=c(0, maxHeight),
             lwd=0.6, ftype="off", mar=c(0.1,0.7,0.1,0.7))
    
    if(tr %% nCols == 1){ axis(side = 2)}
    
    #add tips and color by population
    pop <- as.numeric(unlist(lapply(strsplit(trees[[1]]$tip.label, "[.]"), '[', 2) ))
    
    nTips <- length(trees[[tr]]$tip.label);
    symbols(1:nTips, rep(0, nTips), circles=rep(0.4, nTips), add=TRUE, inches=0.03, fg=NA, bg=popCol[pop])
  }
}

library(phytools)
trees <- read.nexus("C:\\Users\\kundu\\OneDrive\\Desktop\\LAB\\EMBO\\Course\\Daniel\\constsize_1_true_trees.trees")
plotTrees(trees)
