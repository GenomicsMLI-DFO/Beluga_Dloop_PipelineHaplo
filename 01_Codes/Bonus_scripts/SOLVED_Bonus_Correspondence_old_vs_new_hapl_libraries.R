library(dplyr)
"%nin%" <- Negate("%in%")

## Compare haplo libraries

o234 <- read.csv("../Downloads/librairie_51_haplotypes234.csv", stringsAsFactors = F)
o615 <- read.csv("../Downloads/librairie_137_haplotypes615.csv", stringsAsFactors = F)

n234 <- read.csv("./GitHub/Beluga_Dloop_PipelineHaplo/02_Results/00_libraries/librairie_53_haplotypes234.csv", stringsAsFactors = F)
h615 <- read.csv("./GitHub/Beluga_Dloop_PipelineHaplo/02_Results/00_libraries/librairie_144_haplotypes615.csv", stringsAsFactors = F)

n615 <- left_join(h615, o615, by = "hapl")
colnames(n615)[2:3] <- c("seq.new","seq.old")
head(n615)
n615$check <- n615$seq.new == n615$seq.old
table(n615$check)
which(n615$seq.new %nin% n615$seq.old)  # 136 137 138 139 140 141 142 143 144
which(n615$seq.old %nin% n615$seq.new)  # 104, 130 (138:144 are NAs)

for(i in c(104:135)){
  print(which(n615$seq.old %in% n615$seq.new[i]))
}  # corresponding sequence were 1 haplo after (104 doesn't exist anymore)


