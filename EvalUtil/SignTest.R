#!/usr/bin/env Rscript

# The data should be in the format
# 1st line: (space separated) true labels
# 2d  line: (space separated) predicted labels from system 1
# 3d  line: (space separated) predicted labels from system 2
# and so on...
#
# The methods can be selected using arguments:
# enumeration of methods starts from one.
# Method #1 means the method  whose performance numbers are given in the second row.
# Recall the the first row contains true lables (gold standard).


library('stats')
library('BSDA')

StatTest <- function(dt, row1, row2) {
  x <- as.numeric(dt[1,] == dt[row1,])
  y <- as.numeric(dt[1,] == dt[row2,])

  print(x)
  print(y)

  print(mean(x))
  print(mean(y))

  #wilcox.test(x, y, paired=T)  
  SIGN.test(x,y)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage <file name> <mehod #1 num> <method #2 num>");
}
inputFile <- args[1] 
print(paste('SIGN.test, input file: ', inputFile))
dt <- read.csv(inputFile, sep=' ',header=F)


StatTest(dt, as.numeric(args[2])+1, as.numeric(args[3])+1)


