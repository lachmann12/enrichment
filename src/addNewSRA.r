tt = read.table("pvals_test.txt", sep="\t", stringsAsFactors=F)
tt = data.matrix(tt)

ltt = log(tt)
cs = colSums(ltt)
rs = rowSums(ltt)

plot(density((rs), na.rm=T), lwd=3, col="blue")
plot(density((cs), na.rm=T), lwd=3)

for(i in 1:10){
    plot(density(ltt[i,], na.rm=T), lwd=3, col="blue")
}


for(i in 1:10){
    plot(density(ltt[,i], na.rm=T), lwd=3, col="red")
}




plot(density(scale(rs), na.rm=T), lwd=3, col="blue")
lines(density(scale(cs), na.rm=T), lwd=3)








