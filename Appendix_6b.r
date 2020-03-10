#APPENDIX 6B
#The following script was written for analyses in R of GC and codon usage data - Spearman's and plotting



#Optimal codon validation in E. coli - comparison to pre-programmed optimal codons

filelist <- list.files (path = ".", pattern = "AP010953_gc.csv")

for (i in filelist) {

	data1 <- read.csv(i, header = TRUE)

	Egc3 <- unlist(data1[2])
	Efgc3 <- unlist(data1[3])
	Etgc3 <- unlist(data1[4])
}


filelist <- list.files (path = ".", pattern = "AP010953_f.csv")

for (i in filelist) {

	data2 <- read.csv(i, header = TRUE)

	Efop <- unlist(data2[2])
}


filelist <- list.files (path = ".", pattern = "AP010953_fd.csv")

for (i in filelist) {

	data3 <- read.csv(i, header = TRUE)

	Efopd <- unlist(data3[2])
}


corr.test <- cor.test(Efopd, Efop, method = "spearman")
Rho.corr.test <- round(corr.test$estimate, digits = 3)


png("Ec_proof.png", width = 400, height = 400)
par(oma = c(1,1,0,1) + 0.1, mar = c(2,5,1,2) + 0.1)
par(pty="s")


plot(Efopd, Efop, xlab= "Pre-programmed codons", ylab = "Calculated codons", cex.axis=1.1, mgp = c(3, 1, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 1.7)

	abline(lm(Efop~Efopd), col = "black")
	
	
mtext(paste0("r = ", Rho.corr.test), side = 3, adj = 0.95, line = -1.2, cex = 1, col = "dodgerblue")

dev.off()



Ecor.fgc3.fop <- cor.test(Efgc3, Efop, method = "spearman")
Erho.fgc3.fop <- round(Ecor.fgc3.fop$estimate, digits = 3)
EPval.fgc3.fop  <- round(Ecor.fgc3.fop$p.value, digits = 3)

Ecor.tgc3.fop <- cor.test(Etgc3, Efop, method = "spearman")
Erho.tgc3.fop <- round(Ecor.tgc3.fop$estimate, digits = 3)
EPval.tgc3.fop  <- round(Ecor.tgc3.fop$p.value, digits = 3)






png("Ec_termini.png", width = 800, height = 450)
par(mfrow=c(1,2), oma = c(0,0,0,0) + 0.1, mar = c(4,4.5,4,2) + 0.1)
par(pty="s")


plot(Efgc3, Efop, xlab=expression("5' GC"[3]*"%"), ylab = "Gene Expression (Fop)", cex.axis=1.4, mgp = c(3, 1, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 1.7)


	abline(lm(Efop~Efgc3), col = "black")
	
	title("A", adj = 0, cex.main = 2.5)
	
	mtext(paste0("r = ", Erho.fgc3.fop), side = 3, adj = 0.05, line = -1.5, cex = 1.2)



plot(Etgc3, Efop, xlab=expression("3' GC"[3]*"%"), ylab = "Gene Expression (Fop)", cex.axis=1.4, mgp = c(3, 1, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 1.7)


	abline(lm(Efop~Etgc3), col = "black")
	
	title("B", adj = 0, cex.main = 2.5)
	
	mtext(paste0("r = ", Erho.tgc3.fop), side = 3, adj = 0.05, line = -1.5, cex = 1.2)



dev.off()



#Analyses with sample of six genomes

filelist <- list.files (path = ".", pattern = "CP006681_gc.csv")

for (i in filelist) {

	data1 <- read.csv(i, header = TRUE)

	gc3.L <- unlist(data1[2])
	fgc3.L <- unlist(data1[3])
	tgc3.L <- unlist(data1[4])
}

filelist <- list.files (path = ".", pattern = "CP006681_fop.csv")

for (i in filelist) {

	data2 <- read.csv(i, header = TRUE)

	fop.L <- unlist(data2[2])
}



filelist <- list.files (path = ".", pattern = "AP014521_gc.csv")

for (i in filelist) {

	data1a <- read.csv(i, header = TRUE)

	gc3.La <- unlist(data1a[2])
	fgc3.La <- unlist(data1a[3])
	tgc3.La <- unlist(data1a[4])
}

filelist <- list.files (path = ".", pattern = "AP014521_fop.csv")

for (i in filelist) {

	data2a <- read.csv(i, header = TRUE)

	fop.La <- unlist(data2a[2])
}




filelist <- list.files (path = ".", pattern = "CP000038_gc.csv")

for (i in filelist) {

	data3 <- read.csv(i, header = TRUE)

	gc3.M <- unlist(data3[2])
	fgc3.M <- unlist(data3[3])
	tgc3.M <- unlist(data3[4])
}

filelist <- list.files (path = ".", pattern = "CP000038_fop.csv")

for (i in filelist) {

	data4 <- read.csv(i, header = TRUE)

	fop.M <- unlist(data4[2])
}




filelist <- list.files (path = ".", pattern = "CP003415_gc.csv")

for (i in filelist) {

	data3a <- read.csv(i, header = TRUE)

	gc3.Ma <- unlist(data3a[2])
	fgc3.Ma <- unlist(data3a[3])
	tgc3.Ma <- unlist(data3a[4])
}

filelist <- list.files (path = ".", pattern = "CP003415_fop.csv")

for (i in filelist) {

	data4a <- read.csv(i, header = TRUE)

	fop.Ma <- unlist(data4a[2])
}





filelist <- list.files (path = ".", pattern = "CP002810_gc.csv")

for (i in filelist) {

	data5 <- read.csv(i, header = TRUE)

	gc3.H <- unlist(data5[2])
	fgc3.H <- unlist(data5[3])
	tgc3.H <- unlist(data5[4])
}

filelist <- list.files (path = ".", pattern = "CP002810_fop.csv")

for (i in filelist) {

	data6 <- read.csv(i, header = TRUE)

	fop.H <- unlist(data6[2])
}


filelist <- list.files (path = ".", pattern = "CP000769_gc.csv")

for (i in filelist) {

	data5a <- read.csv(i, header = TRUE)

	gc3.Ha <- unlist(data5a[2])
	fgc3.Ha <- unlist(data5a[3])
	tgc3.Ha <- unlist(data5a[4])
}

filelist <- list.files (path = ".", pattern = "CP000769_fop.csv")

for (i in filelist) {

	data6a <- read.csv(i, header = TRUE)

	fop.Ha <- unlist(data6a[2])
}





cor.gc3.fop.L <- cor.test(gc3.L, fop.L, method = "spearman")
rho.gc3.fop.L <- round(cor.gc3.fop.L$estimate, digits = 4)
Pval.gc3.fop.L  <- round(cor.gc3.fop.L$p.value, digits = 4)

cor.fgc3.fop.L <- cor.test(fgc3.L, fop.L, method = "spearman")
rho.fgc3.fop.L <- round(cor.fgc3.fop.L$estimate, digits = 4)
Pval.fgc3.fop.L  <- round(cor.fgc3.fop.L$p.value, digits = 4)

cor.tgc3.fop.L <- cor.test(tgc3.L, fop.L, method = "spearman")
rho.tgc3.fop.L <- round(cor.tgc3.fop.L$estimate, digits = 4)
Pval.tgc3.fop.L  <- round(cor.tgc3.fop.L$p.value, digits = 4)


cor.gc3.fop.La <- cor.test(gc3.La, fop.La, method = "spearman")
rho.gc3.fop.La <- round(cor.gc3.fop.La$estimate, digits = 4)
Pval.gc3.fop.La  <- round(cor.gc3.fop.La$p.value, digits = 4)

cor.fgc3.fop.La <- cor.test(fgc3.La, fop.La, method = "spearman")
rho.fgc3.fop.La <- round(cor.fgc3.fop.La$estimate, digits = 4)
Pval.fgc3.fop.La  <- round(cor.fgc3.fop.La$p.value, digits = 4)

cor.tgc3.fop.La <- cor.test(tgc3.La, fop.La, method = "spearman")
rho.tgc3.fop.La <- round(cor.tgc3.fop.La$estimate, digits = 4)
Pval.tgc3.fop.La  <- round(cor.tgc3.fop.La$p.value, digits = 4)

cor.gc3.fop.M <- cor.test(gc3.M, fop.M, method = "spearman")
rho.gc3.fop.M <- round(cor.gc3.fop.M$estimate, digits = 4)
Pval.gc3.fop.M  <- round(cor.gc3.fop.M$p.value, digits = 4)

cor.fgc3.fop.M <- cor.test(fgc3.M, fop.M, method = "spearman")
rho.fgc3.fop.M <- round(cor.fgc3.fop.M$estimate, digits = 4)
Pval.fgc3.fop.M  <- round(cor.fgc3.fop.M$p.value, digits = 4)

cor.tgc3.fop.M <- cor.test(tgc3.M, fop.M, method = "spearman")
rho.tgc3.fop.M <- round(cor.tgc3.fop.M$estimate, digits = 4)
Pval.tgc3.fop.M  <- round(cor.tgc3.fop.M$p.value, digits = 4)



cor.gc3.fop.Ma <- cor.test(gc3.Ma, fop.Ma, method = "spearman")
rho.gc3.fop.Ma <- round(cor.gc3.fop.Ma$estimate, digits = 4)
Pval.gc3.fop.Ma  <- round(cor.gc3.fop.Ma$p.value, digits = 4)

cor.fgc3.fop.Ma <- cor.test(fgc3.Ma, fop.Ma, method = "spearman")
rho.fgc3.fop.Ma <- round(cor.fgc3.fop.Ma$estimate, digits = 4)
Pval.fgc3.fop.Ma  <- round(cor.fgc3.fop.Ma$p.value, digits = 4)

cor.tgc3.fop.Ma <- cor.test(tgc3.Ma, fop.Ma, method = "spearman")
rho.tgc3.fop.Ma <- round(cor.tgc3.fop.Ma$estimate, digits = 4)
Pval.tgc3.fop.Ma  <- round(cor.tgc3.fop.Ma$p.value, digits = 4)

cor.gc3.fop.H <- cor.test(gc3.H, fop.H, method = "spearman")
rho.gc3.fop.H <- round(cor.gc3.fop.H$estimate, digits = 4)
Pval.gc3.fop.H  <- round(cor.gc3.fop.H$p.value, digits = 4)

cor.fgc3.fop.H <- cor.test(fgc3.H, fop.H, method = "spearman")
rho.fgc3.fop.H <- round(cor.fgc3.fop.H$estimate, digits = 4)
Pval.fgc3.fop.H  <- round(cor.fgc3.fop.H$p.value, digits = 4)

cor.tgc3.fop.H <- cor.test(tgc3.H, fop.H, method = "spearman")
rho.tgc3.fop.H<- round(cor.tgc3.fop.H$estimate, digits = 4)
Pval.tgc3.fop.H  <- round(cor.tgc3.fop.H$p.value, digits = 4)

cor.gc3.fop.Ha <- cor.test(gc3.Ha, fop.Ha, method = "spearman")
rho.gc3.fop.Ha <- round(cor.gc3.fop.Ha$estimate, digits = 4)
Pval.gc3.fop.Ha  <- round(cor.gc3.fop.Ha$p.value, digits = 4)

cor.fgc3.fop.Ha <- cor.test(fgc3.Ha, fop.Ha, method = "spearman")
rho.fgc3.fop.Ha <- round(cor.fgc3.fop.Ha$estimate, digits = 4)
Pval.fgc3.fop.Ha  <- round(cor.fgc3.fop.Ha$p.value, digits = 4)

cor.tgc3.fop.Ha <- cor.test(tgc3.Ha, fop.Ha, method = "spearman")
rho.tgc3.fop.Ha<- round(cor.tgc3.fop.Ha$estimate, digits = 4)
Pval.tgc3.fop.Ha  <- round(cor.tgc3.fop.Ha$p.value, digits = 4)





png("Figure8.png", width = 1100, height = 1550)
par(mfrow=c(3,2), oma = c(2,1,1,1) + 0.1, mar = c(7.5,7.5,5.5,3) + 0.1)
par(pty="s")



plot(gc3.L, fop.L, xlab=expression("CDS GC"[3]*"%"), ylab = "Gene Expression (Fop)", cex.axis=2.5, mgp = c(5, 2, 0), col = "dodgerblue", cex=1.4, pch=20, cex.lab = 3)

	title("CP006681", adj = 0.5, cex.main = 3.5)	
	abline(lm(fop.L~gc3.L), col = "black")


mtext(paste0("r = ", rho.gc3.fop.L), side = 3, adj = 0.05, line = -2.5, cex = 1.5)
mtext(paste0("p = ", Pval.gc3.fop.L), side = 3, adj = 0.05, line = -4.5, cex = 1.5)



plot(gc3.La, fop.La, xlab=expression("CDS GC"[3]*"%"), ylab = "Gene Expression (Fop)", cex.axis=2.5, mgp = c(5, 2, 0), col = "dodgerblue", cex=1.4, pch=20, cex.lab = 3)

	title("AP014521", adj = 0.5, cex.main = 3.5)	
	abline(lm(fop.La~gc3.La), col = "black")


mtext(paste0("r = ", rho.gc3.fop.La), side = 3, adj = 0.05, line = -2.5, cex = 1.5)
mtext(paste0("p = ", Pval.gc3.fop.La), side = 3, adj = 0.05, line = -4.5, cex = 1.5)


plot(gc3.Ma, fop.Ma, xlab=expression("CDS GC"[3]*"%"), ylab = "Gene Expression (Fop)", cex.axis=2.5, mgp = c(5, 2, 0), col = "dodgerblue3", cex=1.4, pch=20, cex.lab = 3)

	title("CP003415", adj = 0.5, cex.main = 3.5)	
	abline(lm(fop.Ma~gc3.Ma), col = "black")


mtext(paste0("r = ", rho.gc3.fop.Ma), side = 3, adj = 0.05, line = -2.5, cex = 1.5)
mtext(paste0("p = ", Pval.gc3.fop.Ma), side = 3, adj = 0.05, line = -4.5, cex = 1.5)


plot(gc3.M, fop.M, xlab=expression("CDS GC"[3]*"%"), ylab = "Gene Expression (Fop)", cex.axis=2.5, mgp = c(5, 2, 0), col = "dodgerblue3", cex=1.4, pch=20, cex.lab = 3)

	title("CP000038", adj = 0.5, cex.main = 3.5)	
	abline(lm(fop.M~gc3.M), col = "black")


mtext(paste0("r = ", rho.gc3.fop.M), side = 3, adj = 0.05, line = -2.5, cex = 1.5)
mtext(paste0("p = ", Pval.gc3.fop.M), side = 3, adj = 0.05, line = -4.5, cex = 1.5)



plot(gc3.Ha, fop.Ha, xlab=expression("CDS GC"[3]*"%"), ylab = "Gene Expression (Fop)", cex.axis=2.5, mgp = c(5, 2, 0), col = "dodgerblue4", cex=1.4, pch=20, cex.lab = 3)

	title("CP000769", adj = 0.5, cex.main = 3.5)	
	abline(lm(fop.Ha~gc3.Ha), col = "black")


mtext(paste0("r = ", rho.gc3.fop.Ha), side = 3, adj = 0.05, line = -2.5, cex = 1.5)
mtext(paste0("p = ", Pval.gc3.fop.Ha), side = 3, adj = 0.05, line = -4.5, cex = 1.5)

plot(gc3.H, fop.H, xlab=expression("CDS GC"[3]*"%"), ylab = "Gene Expression (Fop)", cex.axis=2.5, mgp = c(5, 2, 0), col = "dodgerblue4", cex=1.4, pch=20, cex.lab = 3)

	title("CP002810", adj = 0.5, cex.main = 3.5)	
	abline(lm(fop.H~gc3.H), col = "black")


mtext(paste0("r = ", rho.gc3.fop.H), side = 3, adj = 0.05, line = -2.5, cex = 1.5)
mtext(paste0("p = ", Pval.gc3.fop.H), side = 3, adj = 0.05, line = -4.5, cex = 1.5)


dev.off()

