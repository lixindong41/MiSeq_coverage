#!/usr/bin/env Rscript
#--- Read command line args
args <- commandArgs(trailingOnly=TRUE)
file_coverage = as.character(args[1])
file_corr = as.character(args[2])


library(grid)
library(gridBase)
library(gridExtra)

#--- Get coverage and correlation data

coverage <- read.table(file_coverage,header=F,sep="\t")
colnames(coverage) <- c("start","count")
corr <- read.table(file_corr,header=T,sep="\t")

#--- Print coverage across the genome
pdf(gsub(".tsv", ".pdf", file_coverage))
plot.ts(coverage$start,coverage$count, 
	main = "Coverage across the genome", 
	xlab = "Position across reference", ylab = "count",
	type = "l", col = "red")
dev.off()


#--- Report correlation between the nucleotide content of the reference and the coverage
pdf(gsub(".txt", ".pdf", file_corr))
layout(matrix(c(1:8), 4, 2, byrow = TRUE))

#Define a function for ploting statistics table
plotstatable <- function(regmodel) {	
df <- t(data.frame("R2" = signif(summary(regmodel)$r.squared, 5),
					"Adj R2" = signif(summary(regmodel)$adj.r.squared, 5),
					"Intercept" =signif(regmodel$coef[[1]],5 ),
					"Slope" =signif(regmodel$coef[[2]], 5),
					"P" =signif(summary(regmodel)$coef[2,4], 5)))
frame()
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot) 
grid.table(df)
popViewport(3)
}

#Plot results for A content
plot(corr$num_A, corr$coverage, pch=19, col=rgb(0.5, 0.5, 0.5, 0.5), 
	cex=0.8, xlab = "A-content", ylab = "coverage")
regmodel.A <- lm(corr$coverage ~ corr$num_A)
abline(regmodel.A, col="red")       
plotstatable(regmodel.A)

#Plot results for C content
plot(corr$num_C, corr$coverage, pch=19, col=rgb(0.5, 0.5, 0.5, 0.5), 
	cex=0.8, xlab = "C-content", ylab = "coverage")
regmodel.C <- lm(corr$coverage ~ corr$num_C)
abline(regmodel.C, col="red")       
plotstatable(regmodel.C)

#Plot results for G content
plot(corr$num_G, corr$coverage, pch=19, col=rgb(0.5, 0.5, 0.5, 0.5), 
	cex=0.8, xlab = "G-content", ylab = "coverage")
regmodel.G <- lm(corr$coverage ~ corr$num_G)
abline(regmodel.G, col="red")       
plotstatable(regmodel.G)

#Plot results for T content
plot(corr$num_T, corr$coverage, pch=19, col=rgb(0.5, 0.5, 0.5, 0.5), 
	cex=0.8, xlab = "T-content", ylab = "coverage")
regmodel.T <- lm(corr$coverage ~ corr$num_T)
abline(regmodel.T, col="red")       
plotstatable(regmodel.T)

dev.off()
