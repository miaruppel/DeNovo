library(tidyverse)

# Directory that has raw data
setwd("C:\\Users\\miar\\Desktop\\data")

# Read in data and store as object
fragmentTable <- read_delim("fragmentInfoTable.tsv", delim=" ", col_types = "iciicic")
#view(fragmentTable)

# Create histogram for distribution of missing fragments per peptide
#ggplot(data = fragmentTable, aes(x = Number_Missing)) + geom_histogram(binwidth = 1, fill='blue')

x <- fragmentTable %>% 
  separate_rows(Missing_Fragment_Locations, sep=",", convert=T) %>%
  mutate(cTerminusEnd = Sequence_Length - Missing_Fragment_Locations) %>%
  rename(nTerminusEnd = Missing_Fragment_Locations)


# Create histogram for distribution of missing fragments on c-terminus side
#ggplot(data = x, aes(x = cTerminusEnd)) + geom_histogram(binwidth = 1, fill='red')

# Create histogram for distribution of missing fragments on n-terminus side
#ggplot(data = x, aes(x = nTerminusEnd)) + geom_histogram(binwidth = 1, fill='orange')

# Save altered table 
write.table(x, file='fragmentInfoTable_altered.tsv', sep='/t', row.names = FALSE)
