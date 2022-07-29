library(tidyverse)

data <- read_csv("C:/Users/miar/Desktop/data/predicted.csv")



p <- (ggplot(data, aes(x=mz, y=Intensity, label=mz, color='red')) +
        geom_point() +
        geom_segment( aes(x=mz, xend=mz, y=0, yend=Intensity), size=2, color='red') #+ 
        #geom_text(size=3)
      )


print(p)
