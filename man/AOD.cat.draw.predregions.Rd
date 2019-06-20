\name{AOD.cat.draw.predregions}
\alias{AOD.cat.draw.predregions}
\title{
AOD.cat.draw.predregions
}
\description{
Draws prediction regions for categorical AoD.
} 
\usage{
AOD.cat.draw.predregions(X = Ocotea.cat[, 3:8], group.vec = Ocotea.cat[, 2], 
exp.factor = 1.2, parplotmar = c(3, 3, 3, 3), dist.metric = c("EMC", "ChiSq"), 
scaffolding = c(1, 2), plot.CLPs = NULL, predregions = 1, x.grid = 0.015, 
y.grid = 0.015, pred.region.symbol = 15, pred.region.symbol.size = 0.2, 
pch.sample = c(22, 24, 25), pch.group = c(22, 24, 25), pch.sample.lwd = 2, 
legend.plot = TRUE, CLP.LabelsInPredRegs = FALSE, region.colours = 
UBcolours[1:5], sample.labels = TRUE, sample.labels.cex = 0.7, text.size = 
c(1, 1), colours = UBcolours, cij.null = FALSE, predict.samples = NULL, 
labels = NULL, CLP.multiplier = 1, jitter.factor = 0)
}
