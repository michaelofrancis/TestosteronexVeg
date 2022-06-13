suppressMessages(library(qqman))
suppressMessages(library(tidyverse))

pheno<-c("SHBG", "logSHBG", "RankTest", "RankbioavailableTest")

exp<-c("Veg1", "Veg2")
m<-c(1,2)

#m=1

dir<-c("/scratch/mf91122/TestxVeg/GEM")

mancolors<-list(c("#11DFF0", "#0696A3"),
		c("#F01D99", "#A30A63"),
		c("#F0741D", "#A34A0A"),
		c("#05F07E", "#02A355"))


outdir<-"/scratch/mf91122/TestxVeg/GEM/plot"

for (m in 1:length(m)){ #loop over models

for (i in 1:length(pheno)){ #loop over phenotypes

for (e in 1:length(exp)){ #loop over exposure

#m=1;i=1; e=1

file<-list()
for (c in 1:22){
	file[[c]]<-as_tibble(read.table(paste(dir,"/", 
		pheno[i], "/", exp[e], "/M", m,  "/chr", c, sep=""), header=T))
}

file2<-do.call(rbind, file)
Nplot<-file2$N_Samples[1]

file2<-file2%>%select(CHR, RSID, POS,P_Value_Interaction)
colnames(file2)<-c("CHR", "SNP" , "BP", "P")

suppressWarnings(dir.create(outdir))

plotoutputpath=paste(outdir, "/", pheno[i], "-", exp[e], "-M", m, ".png", sep="")

#Make Manhattan plot P_BOLT_LMM
png(filename=plotoutputpath, type="cairo", width=1000,height=500)

print(
manhattan(file2,
	main = paste("UKB (N=", Nplot, ") ", 
		pheno[i]," ", exp[e], " M", m, " P-interaction", sep=""),
	ylim = c(0, 12), cex = 0.6,
	annotatePval = FALSE, annotateTop = FALSE,
	cex.axis = 0.9, 
	col = mancolors[[i]],
	suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
	chrlabs = as.character(1:22)
	)#end manhattan
)

dev.off() 


qqoutputpath=paste(outdir, "/", pheno[i], "-", exp[e], "-M", m, ".qq.png", sep="")
png(filename=qqoutputpath, type="cairo", width=500,height=500)

qq(file2$P,
	main=paste("UKB ", pheno[i]," ", exp[e], " M", m, " P-interaction", sep="")
)

dev.off()

}#end loop over models
} #end loop over phenotypes
}
