library(tidyverse)

# Directory that has raw data
setwd("C:\\Users\\miar\\Desktop\\data")

# Read in data and store as object
fragmentTable <- read_delim("fragmentInfoTable.tsv")

view(fragmentTable)
