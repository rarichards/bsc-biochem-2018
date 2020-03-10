#APPENDIX 6A

#The following script was written for analyses in R for GC and temperature data - Spearman's tests and plotting

filelist <- list.files (path = ".", pattern = "arc_temps.csv")

for (i in filelist) {

	data1 <- read.csv(i, header = TRUE)

#	first I want to see if there is a correlation between GC and GC3

	accn.arc <- unlist(data1[1])
	genGC.arc <- unlist(data1[4])
	cdsGC3.arc <- unlist(data1[5])
	fveGC3.arc <- unlist(data1[6])
	thrGC3.arc <- unlist(data1[7])
	srnaGC.arc <- unlist(data1[8])
	res.arc <- unlist(data1[9])
	temp.arc <- unlist(data1[10])
}




filelist <- list.files (path = ".", pattern = "bac_temps.csv")

for (i in filelist) {

	data2 <- read.csv(i, header = TRUE)

#	first I want to see if there is a correlation between GC and GC3

	accn.bac <- unlist(data2[1])
	genGC.bac <- unlist(data2[4])
	cdsGC3.bac <- unlist(data2[5])
	fveGC3.bac <- unlist(data2[6])
	thrGC3.bac <- unlist(data2[7])
	srnaGC.bac <- unlist(data2[8])
	res.bac <- unlist(data2[9])
	temp.bac <- unlist(data2[10])
}


###	FIGURE 1 ###

bac.cor.gc <- cor.test(cdsGC3.bac, genGC.bac, method = "spearman")
bac.rho.gc <- round(bac.cor.gc$estimate, digits = 3)
bac.Pval.gc  <- round(bac.cor.gc$p.value, digits = 3)


arc.cor.gc <- cor.test(cdsGC3.arc, genGC.arc, method = "spearman")
arc.rho.gc <- round(arc.cor.gc$estimate, digits = 3)
arc.Pval.gc  <- round(arc.cor.gc$p.value, digits = 3)



png("Figure1.png", width = 800, height = 450)
par(mfrow=c(1,2), oma = c(0,0,0,0) + 0.1, mar = c(4,4.5,4,2) + 0.1)
par(pty="s")



plot(cdsGC3.bac, genGC.bac, xlab=expression('Genomic CDS GC'[3]*'%'), ylab = "Genomic GC%", cex.axis=1.4, mgp = c(3, 1, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 1.7)

	title("A", adj = 0, cex.main = 2.5)
	abline(lm(genGC.bac~cdsGC3.bac), col = "black")

	
	mtext(paste0("r = ", bac.rho.gc), side = 3, adj = 0.05, line = -1.5, cex = 1.2)


plot(cdsGC3.arc, genGC.arc, xlab=expression('Genomic CDS GC'[3]*'%'), ylab = "Genomic GC%", cex.axis=1.4, mgp = c(3, 1, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 1.7)

	title("B", adj = 0, cex.main = 2.5)
	abline(lm(genGC.arc~cdsGC3.arc), col = "black")

	
	mtext(paste0("r = ", arc.rho.gc), side = 3, adj = 0.05, line = -1.5, cex = 1.2)

dev.off()


###	FIGURE 2 ###

png("Figure2.png", width = 800, height = 450)
par(mfrow=c(1,2), oma = c(0,0,0,0) + 0.1, mar = c(4,4.5,4,2) + 0.1)
par(pty="s")

plot(cdsGC3.bac, cdsGC3.bac, pch=20, col = "black", cex=0.8, xlab = expression('Genomic CDS GC'[3]*'%'), ylab =expression('GC'[3]*'%'), cex.axis=1.4, mgp = c(3,1,0), cex.lab = 1.7)
points(cdsGC3.bac, fveGC3.bac, pch=20, col = "dodgerblue", cex=0.8)
points(cdsGC3.bac, thrGC3.bac, pch=20, cex=0.8, col = "red")

legend("topleft", c("CDS", "5'", "3'"), col=c("black", "dodgerblue", "red"), pch=20, cex = 1.2)
title("A", adj = 0, cex.main = 2.3)

plot(cdsGC3.arc, cdsGC3.arc, pch=20, col = "black", cex=0.8, xlab = expression('Genomic CDS GC'[3]*'%'), ylab =expression('GC'[3]*'%'), cex.axis=1.4, mgp = c(3,1,0), cex.lab = 1.7)
points(cdsGC3.arc, fveGC3.arc, pch=20, col = "dodgerblue", cex=0.8)
points(cdsGC3.arc, thrGC3.arc, pch=20, cex=0.8, col = "red")

legend("topleft", c("CDS", "5'", "3'"), col=c("black", "dodgerblue", "red"), pch=20, cex = 1.2)
title("B", adj = 0, cex.main = 2.3)


dev.off()


###	FIGURE 3 ###



bac.cor.rna <- cor.test(cdsGC3.bac, srnaGC.bac, method = "spearman")

bac.rho.rna <- round(bac.cor.rna$estimate, digits = 4)
bac.Pval.rna  <- round(bac.cor.rna$p.value, digits = 4)




arc.cor.rna <- cor.test(cdsGC3.arc, srnaGC.arc, method = "spearman")

arc.rho.rna <- round(arc.cor.rna$estimate, digits = 4)
arc.Pval.rna  <- round(arc.cor.rna$p.value, digits = 4)





png("Figure3.png", width = 800, height = 450)
par(mfrow=c(1,2), oma = c(0,0,0,0) + 0.1, mar = c(4,4.5,4,2) + 0.1)
par(pty="s")


plot(cdsGC3.bac, srnaGC.bac, xlab=expression('Genomic CDS GC'[3]*'%'), ylab = "Structural RNA GC%", cex.axis=1.4, xlim=c(0, 100), ylim=c(35, 70), mgp = c(3, 1, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 1.7)
 
 	abline(lm(srnaGC.bac~cdsGC3.bac), col = "black")
 	title("A", adj = 0, cex.main = 2.5)
 
plot(cdsGC3.arc, srnaGC.arc, xlab=expression('Genomic CDS GC'[3]*'%'), ylab = "Structural RNA GC%", cex.axis=1.4, mgp = c(3, 1, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 1.7)
 
 	abline(lm(srnaGC.arc~cdsGC3.arc), col = "black")
 	title("B", adj = 0, cex.main = 2.5)
 	
 	
 
 dev.off()








###	FIGURE 4 ###

filelist <- list.files (path = ".", pattern = "arct_temps.csv")

for (i in filelist) {

	data3 <- read.csv(i, header = TRUE)

#	first I want to see if there is a correlation between GC and GC3

	accn.arct <- unlist(data3[1])
	genGC.arct <- unlist(data3[4])
	cdsGC3.arct <- unlist(data3[5])
	fveGC3.arct <- unlist(data3[6])
	thrGC3.arct <- unlist(data3[7])
	srnaGC.arct <- unlist(data3[8])
	res.arct <- unlist(data3[9])
	temp.arct <- unlist(data3[10])
}




filelist <- list.files (path = ".", pattern = "bact_temps.csv")

for (i in filelist) {

	data4 <- read.csv(i, header = TRUE)

#	first I want to see if there is a correlation between GC and GC3

	accn.bact <- unlist(data4[1])
	genGC.bact <- unlist(data4[4])
	cdsGC3.bact <- unlist(data4[5])
	fveGC3.bact <- unlist(data4[6])
	thrGC3.bact <- unlist(data4[7])
	srnaGC.bact <- unlist(data4[8])
	res.bact <- unlist(data4[9])
	temp.bact <- unlist(data4[10])
}

#FOR FILTERED BACTERIA


#residuals vs temperature
cor.temp.resb <- cor.test(temp.bact, res.bact, method = "spearman")
rho.temp.resb <- round(cor.temp.resb$estimate, digits = 3)
Pval.temp.resb  <- round(cor.temp.resb$p.value, digits = 3)


#FOR FILTERED ARCHAEA
#residuals vs temperature
cor.temp.resa <- cor.test(temp.arct, res.arct, method = "spearman")
rho.temp.resa <- round(cor.temp.resa$estimate, digits = 3)
Pval.temp.resa  <- round(cor.temp.resa$p.value, digits = 3)




png("Figure4.png", width = 800, height = 450)
par(mfrow=c(1,2), oma = c(0,0,0,0) + 0.1, mar = c(4,4.5,4,2) + 0.1)
par(pty="s")


plot(temp.bact, res.bact, xlab="Optimal Growth Temperature (°C)", ylab = "Residuals (sRNA GC% ~ CDS GC%)", cex.axis=1.4, mgp = c(3, 1, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 1.55)

	abline(lm(res.bact~temp.bact), col = "black")
	title("A", adj = 0, cex.main = 2.5)
	

	mtext(paste0("r = ", rho.temp.resb), side = 3, adj = 0.05, line = -1.5, cex = 1.2)


plot(temp.arct, res.arct, xlab="Optimal Growth Temperature (°C)", ylab = "Residuals (sRNA GC% ~ CDS GC%)", cex.axis=1.4, mgp = c(3, 1, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 1.55)

	abline(lm(res.arct~temp.arct), col = "black")
	title("B", adj = 0, cex.main = 2.5)
	

	mtext(paste0("r = ", rho.temp.resa), side = 3, adj = 0.05, line = -1.5, cex = 1.2)


dev.off()




###	FIGURE 5 ###


bac.cor.cds <- cor.test(res.bac, cdsGC3.bac, method = "spearman")
bac.rho.cds <- round(bac.cor.cds$estimate, digits = 3)
bac.Pval.cds  <- round(bac.cor.cds$p.value, digits = 3)


bac.cor5.res <- cor.test(res.bac, fveGC3.bac, method = "spearman")

bac.rho5.res <- round(bac.cor5.res$estimate, digits = 3)
bac.Pval5.res  <- round(bac.cor5.res$p.value, digits = 3)



bac.cor3.res <- cor.test(res.bac, thrGC3.bac, method = "spearman")

bac.rho3.res <- round(bac.cor3.res$estimate, digits = 3)
bac.Pval3.res  <- round(bac.cor3.res$p.value, digits = 3)





png("Figure5c_bac.png", width = 1250, height = 450)
par(mfrow=c(1,3), oma = c(0,0,0,0) + 0.1, mar = c(6,7.5,6,3) + 0.1)
par(pty="s")


 plot(res.bac, cdsGC3.bac, xlab="Growth temperature (res)", ylab = expression("CDS GC"[3]*"%"), xlim=c(-10, 10), ylim=c(0, 80), cex.axis=2.5, mgp = c(5, 2, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 3)

 
 	title("A", adj = 0, cex.main = 5)
 	abline(lm(cdsGC3.bac~res.bac), col = "black")
  
 	mtext(paste0("r = ", bac.rho.cds), side = 3, adj = 0.05, line = -2.5, cex = 1.5)
 
 
 plot(res.bac, fveGC3.bac, xlab="Growth temperature (res)", ylab = expression("5' GC"[3]*"%"), xlim=c(-10, 10), ylim=c(0, 80), cex.axis=2.5, mgp = c(5, 2, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 3)

 
 	title("B", adj = 0, cex.main = 5)
 	abline(lm(fveGC3.bac~res.bac), col = "black")
  
 	mtext(paste0("r = ", bac.rho5.res), side = 3, adj = 0.05, line = -2.5, cex = 1.5)
 
 plot(res.bac, thrGC3.bac, xlab="Growth temperature (res)", ylab = expression("3' GC"[3]*"%"), xlim=c(-10, 10), ylim=c(0, 80), cex.axis=2.5, mgp = c(5, 2, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 3)

 
 	title("C", adj = 0, cex.main = 5)
 	abline(lm(thrGC3.bac~res.bac), col = "black")
  
 	mtext(paste0("r = ", bac.rho3.res), side = 3, adj = 0.05, line = -2.5, cex = 1.5)
 
dev.off()

#For archaeal genomes...

arc.cor.cds <- cor.test(res.arc, cdsGC3.arc, method = "spearman")
arc.rho.cds <- round(arc.cor.cds$estimate, digits = 3)
arc.Pval.cds  <- round(arc.cor.cds$p.value, digits = 3)


arc.cor5.res <- cor.test(res.arc, fveGC3.arc, method = "spearman")

arc.rho5.res <- round(arc.cor5.res$estimate, digits = 3)
arc.Pval5.res  <- round(arc.cor5.res$p.value, digits = 3)



arc.cor3.res <- cor.test(res.arc, thrGC3.arc, method = "spearman")

arc.rho3.res <- round(arc.cor3.res$estimate, digits = 3)
arc.Pval3.res  <- round(arc.cor3.res$p.value, digits = 3)





png("Figure5c_arc.png", width = 1250, height = 450)
par(mfrow=c(1,3), oma = c(0,0,0,0) + 0.1, mar = c(6,7.5,6,3) + 0.1)
par(pty="s")


 plot(res.arc, cdsGC3.arc, xlab="Growth temperature (res)", ylab = expression("CDS GC"[3]*"%"), xlim=c(-10, 10), ylim=c(0, 80), cex.axis=2.5, mgp = c(5, 2, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 3)

 
 	title("A", adj = 0, cex.main = 5)
 	abline(lm(cdsGC3.arc~res.arc), col = "black")
  
 	mtext(paste0("r = ", arc.rho.cds), side = 3, adj = 0.05, line = -2.5, cex = 1.5)
 
 
 plot(res.arc, fveGC3.arc, xlab="Growth temperature (res)", ylab = expression("5' GC"[3]*"%"), xlim=c(-10, 10), ylim=c(0, 80), cex.axis=2.5, mgp = c(5, 2, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 3)

 
 	title("B", adj = 0, cex.main = 5)
 	abline(lm(fveGC3.arc~res.arc), col = "black")
  
 	mtext(paste0("r = ", arc.rho5.res), side = 3, adj = 0.05, line = -2.5, cex = 1.5)
 
 plot(res.arc, thrGC3.arc, xlab="Growth temperature (res)", ylab = expression("3' GC"[3]*"%"), xlim=c(-10, 10), ylim=c(0, 80), cex.axis=2.5, mgp = c(5, 2, 0), col = "dodgerblue", cex=0.8, pch=20, cex.lab = 3)

 
 	title("C", adj = 0, cex.main = 5)
 	abline(lm(thrGC3.arc~res.arc), col = "black")
  
 	mtext(paste0("r = ", arc.rho3.res), side = 3, adj = 0.05, line = -2.5, cex = 1.5)
 
dev.off()


