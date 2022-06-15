library(tidyverse)

# Directory that has raw data
setwd("C:\\Users\\miar\\Desktop\\data")

# Read in data and store as object
# FYI, the table exported from python is quite whack (should fix)
fragmentTable <- read_delim("fragmentsInfoTable.tsv", delim = '|')

view(fragmentTable)
