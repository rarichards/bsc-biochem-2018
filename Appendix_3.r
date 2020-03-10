#APPENDIX 3

#The following script was used to read in the GC data csv file and calculate the residual values in the regression of structural RNA GC3 against CDS GC3, to use as a predictor of optimal growth temperature

#For BACTERIAL GENOMES

#Read file in

filelist <- list.files (path = ".", pattern = "bac_genomes.csv")

for (i in filelist) {

	data1 <- read.csv(i, header = TRUE)

	accn.bac <- unlist(data1[1])
	genGC.bac <- unlist(data1[2])
	cdsGC3.bac <- unlist(data1[3])
	fveGC3.bac <- unlist(data1[4])
	thrGC3.bac <- unlist(data1[5])
	srnaGC.bac <- unlist(data1[6])
}

#Find residuals in the regression

fit <- lm(srnaGC.bac~cdsGC3.bac)
res <- residuals(fit)

#add a column, titled 'RNAres', with value res, to data1
data1$RNAres=res


# write data1 to a new file
write.csv(data1, file = "bac_genomes_res.csv", row.names = FALSE)




#Repeat for ARCHAEAL GENOMES

filelist <- list.files (path = ".", pattern = "arc_genomes.csv")

for (i in filelist) {

	data2 <- read.csv(i, header = TRUE)

	accn.arc <- unlist(data1[1])
	genGC.arc <- unlist(data1[2])
	cdsGC3.arc <- unlist(data1[3])
	fveGC3.arc <- unlist(data1[4])
	thrGC3.arc <- unlist(data1[5])
	srnaGC.arc <- unlist(data1[6])
}


fit <- lm(srnaGC.arc~cdsGC3.arc)
res <- residuals(fit)

data1$RNAres=res

write.csv(data2, file = "arc_genomes_res.csv", row.names = FALSE)

