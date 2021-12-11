# library
library(dplyr)

# initialize
fn <- "D:/Michal_Jan2021/R-pivot-table/test2.txt" # CSV filename (I used a modified version of a sample dataset)
data <- read.csv( fn, header=FALSE, sep="\t" )

#check data
head(data)

#remove NA rows (to avoid errors later)
data <- na.omit(data)

## basic code to pivot the data
data %>%
    group_by(variety) %>%
    summarise(median=median(sepal.length), sd=round(sd(sepal.length),2))

## group_by multiple columns
data %>%
    group_by(.dots=c("variety","year") ) %>%
    summarise_all(mean)

## select specific columns for summary
data %>%
    group_by(.dots=c("variety","year") ) %>%
    select(petal.length, petal.width) %>%    
    summarise_all(median)

# save pivot table data matrix
save_pivot <- data %>%
    group_by(.dots=c("variety","year") ) %>%
    select(petal.length, petal.width) %>%    
    summarise_all(median)

# export as a tab-delimited text file 
write.table(save_pivot,"saved.txt")
