#!/usr/bin/env Rscript

# The data should be in the format
# 1st line: (space separated) true labels
# 2d  line: (space separated) predicted labels from system 1
# 3d  line: (space separated) predicted labels from system 2

library('stats')
library('BSDA')

StatTest <- function(dt) {
  x <- as.numeric(dt[1,] == dt[2,])
  y <- as.numeric(dt[1,] == dt[3,])

  print(x)
  print(y)

  print(mean(x))
  print(mean(y))

  #wilcox.test(x, y, paired=T)  
  SIGN.test(x,y)
}

args <- commandArgs(trailingOnly = TRUE)
inputFile <- args[1] 
print(paste('SIGN.test, input file: ', inputFile))
dt <- read.csv(inputFile, sep=' ',header=F)

StatTest(dt)


