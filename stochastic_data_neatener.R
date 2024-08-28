#This code should take the output of the stochastic model as a csv file, remove all the trailing junk and generally neaten the data, and then return that as a csv file

library(tidyverse)
library(prodlim)

files <- c(1,3,5,6,9,10)

if(!dir.exists("./../Desktop/SAT_models/trimmed")){
  dir.create("./../Desktop/SAT_models/trimmed")
}


for(i in files){
  path <- paste0("./../Desktop/SAT_models/",i,"/1/data.csv")
  data <- read.csv("./../Desktop/SAT_models/1/1/data.csv") %>% select(-X)
  match <- c(rep(0,11))
  n <- row.match(match,data)
  data <- data[1:(n-1),]
  path_out <- paste0("./../Desktop/SAT_models/trimmed/data_",i,".csv")
  write.csv(data,path_out,row.names = FALSE,)
}




