setwd("~/OneDrive/eclipse/FastEnrichment/")

ranks = read.table("output/rank_human.txt", sep="\t", quote="", stringsAsFactors=F)
rownames(ranks) = ranks[,1]
si = ranks[,2]
ranks = data.matrix(ranks[,-c(1,2)])


pval = read.table("output/pval_human.txt", sep="\t", quote="", stringsAsFactors=F)
rownames(pval) = pval[,1]
pval = data.matrix(pval[,-c(1,2)])

rr = rowMeans(ranks)

ll = log(pval)


coco = cor(ll)

png("output/heat_user1.png", width = 1480, height = 1480)
heatmap(coco[1:4000,1:4000], scale="none")
dev.off()

png("output/heat_user2.png", width = 1480, height = 1480)
heatmap(coco[1:100,1:100], scale="none")
dev.off()


coco = cor(t(ll)[,1:1000])

png("output/heat1.png", width = 1480, height = 1480)
heatmap(coco[1:1000,1:1000], scale="none")
dev.off()

png("output/heat2.png", width = 1480, height = 1480)
heatmap(coco[1:100,1:100], scale="none")
dev.off()

