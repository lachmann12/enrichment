setwd("~/OneDrive/eclipse/FastEnrichment/")

ranks = read.table("output/rank_human.txt", sep="\t", quote="", stringsAsFactors=F)
rownames(ranks) = ranks[,1]
ranks = data.matrix(ranks[,-1])
