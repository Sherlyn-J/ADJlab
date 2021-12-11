library(dplyr)
library(ggplot2)


a <- read.csv("D:/Ito-data-plots-Michal-paper-Feb21/A_diff.csv")
b <- read.csv("D:/Ito-data-plots-Michal-paper-Feb21/B_diff.csv")
c <- read.csv("D:/Ito-data-plots-Michal-paper-Feb21/C_diff.csv")

a <- na.omit(a)
b <- na.omit(b)
c <- na.omit(c)

a <- a %>%
  dplyr::filter( log2FC < -1.5 | log2FC > 1.5)
write.table(a,"D:/Ito-data-plots-Michal-paper-Feb21/A_lfcsig.csv")

b <- b %>%
  dplyr::filter( log2FC < -1.5 | log2FC > 1.5)
write.table(b,"D:/Ito-data-plots-Michal-paper-Feb21/B_lfcsig.csv")

c <- c %>%
  dplyr::filter( log2FC < -1.5 | log2FC > 1.5)
write.table(c,"D:/Ito-data-plots-Michal-paper-Feb21/C_lfcsig.csv")

# x <- Reduce(intersect, list(a$NAME,b$NAME,c$NAME))
