\name{AOD.cat}
\alias{AOD.cat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
AOD.cat
}
\description{
Main function for constructing categorical AoD biplots}
\usage{
AOD.cat(X = Ocotea.cat[, 3:8], group.vec = 1:37, exp.factor = 1.2, scaffolding = 
c(1, 2), parplotmar = c(3, 3, 3, 3), dist.metric = c("EMC", "ChiSq"), 
plot.CLPs = 1, plot.CatCentroids = 1, CatCentr.group = c("All", 1, 2, 3), 
CatCentroids.col = rep("purple", 5), CatCentroids.pch = rep(24, 5), 
CatCentroids.cex = 1, CatCentroids.labels.col = 1:ncol(X), 
CatCentroids.labels.cex = 0.8, predregions = NULL, x.grid = 0.015, 
y.grid = 0.015, pred.region.symbol = 15, pred.region.symbol.size = 0.2, 
pch.sample = rep(16, nrow(X)), pch.samples.col = rep("green", 37), pch.group = 
rep(15, 37), pch.group.col = rep("green", 37), pch.group.cex = 1, 
plot.CLPs.pch = rep(24, ncol(X)), plot.CLPs.col = rep("red", 5), 
plot.CLPs.cex = 1, plot.CLPs.labels.cex = 0.8, legend.plot = FALSE, 
CLP.PointsInPredRegs = TRUE, CLP.LabelsInPredRegs = TRUE, region.colours = 
c("pink", "orange", "yellow", "cyan", "grey"), sample.labels = TRUE, 
sample.labels.cex = 0.7, text.size = c(1, 1), colours = UBcolours, 
cij.null = FALSE, predict.samples = NULL, labels = NULL, CLP.multiplier = 1, 
jitter = FALSE, jitter.factor = 5)
}
