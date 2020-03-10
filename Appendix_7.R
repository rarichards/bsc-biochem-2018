#APPENDIX 7
#The following script was used to carry out nonparametric Monte Carlo simulation on Topt data. 

filelist1 <- list.files (path = ".", pattern = "bac_genomes_res.csv")

for (i in filelist1) {

	data1 <- read.csv(i, header = TRUE)

	cdsGC3.bac <- unlist(data1[3])
	fveGC3.bac <- unlist(data1[4])
	thrGC3.bac <- unlist(data1[5])
	res.bac <- unlist (data1[7])
}


bac.cor.cds <- cor.test(res.bac, cdsGC3.bac, method = "spearman")
bac.rho.cds <- round(bac.cor.cds$estimate, digits = 3)
bac.Pval.cds  <- round(bac.cor.cds$p.value, digits = 3)


bac.cor5.res <- cor.test(res.bac, fveGC3.bac, method = "spearman")
bac.rho5.res <- round(bac.cor5.res$estimate, digits = 3)
bac.Pval5.res  <- round(bac.cor5.res$p.value, digits = 3)


bac.cor3.res <- cor.test(res.bac, thrGC3.bac, method = "spearman")
bac.rho3.res <- round(bac.cor3.res$estimate, digits = 3)
bac.Pval3.res  <- round(bac.cor3.res$p.value, digits = 3)

od.rho5 <- bac.rho5.res - bac.rho.cds
od.rho3 <- bac.rho3.res - bac.rho.cds

fve.all <- c(cdsGC3.bac, fveGC3.bac)
thr.all <- c(cdsGC3.bac, thrGC3.bac)
temp.all <- c(res.bac, res.bac)

nums <- c(1:length(temp.all))

lst5 <- c()
lst3 <- c()

 for (i in c(1:1000)) {
 	
 x <- length(temp.all)	
 rV <- sample(x, (x/2), replace = FALSE)
 
 
 r5GC.a <- fve.all[rV]
 r5GC.b <- fve.all[-rV]
 
 r3GC.a <- thr.all[rV]
 r3GC.b <- thr.all[-rV]
 
 rTemp.a <- temp.all[rV]
 rTemp.b <- temp.all[-rV]
 
 
 all5.test.a <- cor.test(r5GC.a, rTemp.a)
 all5.rho.a <- round(all5.test.a$estimate, digits = 3)
 
 all5.test.b <- cor.test(r5GC.b, rTemp.b)
 all5.rho.b <- round(all5.test.b$estimate, digits = 3)

 rd.rho5 <- all5.rho.a - all5.rho.b

 lst5 <- c(lst5, rd.rho5)
 


 all3.test.a <- cor.test(r3GC.a, rTemp.a)
 all3.rho.a <- round(all3.test.a$estimate, digits = 3)
  
 all3.test.b <- cor.test(r3GC.b, rTemp.b)
 all3.rho.b <- round(all3.test.b$estimate, digits = 3)
 
 rd.rho3 <- all3.rho.a - all5.rho.b
 
 lst3 <- c(lst3, rd.rho3)
 
 }
 

mrthn5 <- lst5[lst5 >= od.rho5]
pval5 <- ((length(mrthn5) + 1)/(1000 + 1))


mrthn3 <- lst3[rd.rho3 >= od.rho3]
pval3 <- ((length(mrthn3) + 1)/(1000 + 1))

