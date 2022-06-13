suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
suppressMessages(library(RNOmni))

setwd("/work/kylab/mike/TestxVeg/pheno")

source("manyColsToDummy.R")
source("/work/kylab/mike/TestxVeg/pheno/cbat.R")

hormonemeds<-read.csv("/work/kylab/mike/TestxVeg/pheno/sexHormoneMeds.csv")
source("/work/kylab/mike/TestxVeg/pheno/SSRI_Meds.R")


withdrawn<-read.csv("w48818_20210809.csv", header = FALSE)

QCids<-read.table("/scratch/mf91122/LipidsxVeg/pheno/bd_QC-keep.txt",header=TRUE) #same initial participant QC as this project

###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Load data=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#Load UK Biobank datasets-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
source('/work/kylab/mike/PUFA-GWAS/pheno/load_UKBphenotables.b.R') #20 min

#Phenotypes  ------------------------------------------------------------------------------
#Covariates 
#Model 1: sex, age, age squared, genotyping array, and assessment center indicators (sites of recruitment); 
#lipid medication, socioeconomic status measured by Townsend deprivation index;  

pheno<-bd%>%select(f.eid, f.21003.0.0, f.31.0.0, 
			f.189.0.0, f.30800.0.0,
                   f.30890.0.0, f.30600.0.0, 
                   f.54.0.0, f.22000.0.0,
                f.30850.0.0, f.30830.0.0)

colnames(pheno)<-c("IID", "Age", "Sex",  
			"Townsend", "Oestradiol",
                   "VitaminD_Blood", "Albumin", 
                   "Assessment_center", "Geno_batch",
			"Test", "SHBG")

pheno2<-bd_join4%>%select(f.eid, f.74.0.0, f.30050.0.0,
                          f.21001.0.0, f.20116.0.0, f.20160.0.0,
                          f.6177.0.0,f.6153.0.0,
				        f.3166.0.0,f.2724.0.0
                          )
colnames(pheno2)<-c("IID","Fasting_time", "Mean_corpuscular_haemoglobin",
                    "BMI", "SmokeStatus","Ever_smoked",                    
                    "lipid_med", "lipid_med_plushormones",
			        "blood_draw_time", "Menopause"
                    )

new<-left_join(pheno, pheno2, by="IID")
new<-as_tibble(new)

#Remove withdrawn participants------------------------------------
new<-new[!(new$IID %in% withdrawn$V1), ]

#QC participants via output of UKB_participantQC.R----------------
new<-new[(new$IID %in% QCids$IID),]

#Age squared----------------------------
new$Age2<-new$Age^2

#Make dummy 0/1 cols for each assessment center----------------------
#table(pheno$Assessment_center)
centers<-unique(new$Assessment_center)
centercols<-paste("center", 1:22, sep="")
new[centercols]<-0

for (i in 1:length(centers)){
    new[new$Assessment_center==centers[i],][centercols[i]]<-1
}

new<-new%>%select(-Assessment_center)
new

#Genotype batch------------------------------------
new$Geno_batch1<-0
new$Geno_batch1[new$Geno_batch>0]<-1
#sum(pheno$Geno_batch1) #[1] 438313
new$Geno_batch<-new$Geno_batch1
new<-new%>%select(-Geno_batch1)
#table(new$Geno_batch) #it worked


#Switch sex values to numeric---------------------------
new$Sex<-mapvalues(as.character(new$Sex), 
                     c("Male", "Female"), c(0,1))
new$Sex<-as.numeric(new$Sex)


#ADD VEGETARIAN COLUMNS---------------
VEG<-read.table("/scratch/mf91122/LipidsxVeg/pheno/Vegetarian/vegQC2_04032021.txt", 
		header=T, stringsAsFactors=F)

new<-inner_join(new, VEG, by="IID")


#ADD PCs-------------------------------
pcfile<-"/scratch/mf91122/LipidsxVeg/PLINKPCA/outpc2-10/fullpca-10.eigenvec"
pc<-as_tibble(read.table(pcfile, header=T))
new<-left_join(new, pc)


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#EXTRA PROCESSING FOR T COLUMNS---------------------------------------
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# SHBG Transformations---------------------------------------------
new$SHBG2<-new$SHBG^2
new$logSHBG<-log(new$SHBG)


#Calculate bioavailable T
cbatoutput<-cbat(new$Test, SHBG = new$SHBG, 
                 ALB = new$Albumin)
new$freeTest<-cbatoutput$freeTest
new$bioavailableTest<-cbatoutput$bioavailableTest

#Menopause 
table(new$Menopause, useNA = "always")
new$Menopause[new$Menopause==-3]<-NA
new$Menopause[new$Menopause==2]<-NA
new$Menopause[new$Menopause==3]<-NA


# Taking hormone medications  ----------------------------------

medcols<-c("f.eid", sprintf("f.20003.0.%s", 0:47))
medcodes<-hormonemeds$code #sourced from sexHormoneMeds.csv
names(medcodes)<-hormonemeds$name
manyColsToDummy(medcodes, bd_join4[,medcols], "medoutput") 
colnames(medoutput)<-names(medcodes)
medoutput$tookMed<-apply(medoutput, 1, sum, na.rm=TRUE)
nrow(medoutput[medoutput$tookMed>0,]) #12754
medoutput$tookMed[medoutput$tookMed>1]<-1
medoutput$IID<-bd_join4$f.eid
medoutput2<-medoutput%>%select(IID, tookMed)
new<-left_join(new, medoutput2, by= "IID")

# Taking SSRIs -------------------------------------------------

manyColsToDummy(ssricodes, bd_join4[,medcols], "ssrioutput")
colnames(ssrioutput)<-names(ssricodes)
ssrioutput$tookSSRI<-apply(ssrioutput, 1, sum, na.rm=TRUE)
nrow(ssrioutput[ssrioutput$tookSSRI>0,])
ssrioutput$tookSSRI[ssrioutput$tookSSRI>1]<-1
ssrioutput$IID<-bd_join4$f.eid
ssrioutput2<-ssrioutput%>%select(IID, tookSSRI)
new<-left_join(new, ssrioutput2, by= "IID")


#Z-transformation of blood draw times-----------------------------
new$zblood<-str_split(new$blood_draw_time, " ", 
                        simplify = TRUE)[,2]
new$zblood<-str_split(new$zblood, ":", simplify=TRUE)[,c(1,2)]
new$zblood[,2]<-as.numeric(new$zblood[,2])/60
new$zblood<-as.numeric(new$zblood[,1])+as.numeric(new$zblood[,2])
#hist(new$zblood, na.rm=TRUE)
new$zblood<-(new$zblood - mean(new$zblood, na.rm=TRUE)) / 
    sd(new$zblood, na.rm=TRUE)

#Save before transform
outdir="/scratch/mf91122/TestxVeg/pheno"
new$FID<-new$IID
new<-new%>%select(FID, everything())
colnames(new)[44:45]<-c("Veg1", "Veg2")
participants<-new%>%select(FID,IID)


#Full
write.csv(new,
	paste(outdir, "/TestxVeg-pretransform.csv", sep=""),
        row.names=FALSE, quote=FALSE)


#---------------------------------
#PART 2 starts here===============
#---------------------------------


#REMOVE participants with maf 0.01, mind 0.05, geno 0.02, hwe 1e-06, info >=0.5

outdir="/scratch/mf91122/TestxVeg/pheno"
new<-as_tibble(read.csv(paste(outdir, "/TestxVeg-pretransform.csv", sep="")))

#sapply(new, summary)

#center21 is empty, center 22 for dummy var removal
new<-new%>%select(-center21, -center22)

#Remove 4,373 people on hormone meds
new<-new%>%filter(tookMed!=1)

#Split up into SHBG and Testosterone cohorts
test<-new[complete.cases(new$bioavailableTest),]
shbg<-new[complete.cases(new$SHBG),]

sapply(shbg, summary)


shbg1<-shbg%>%select(FID, IID, SHBG, logSHBG, SHBG2, 
    Age, Sex, Geno_batch,
    paste("center", 1:20, sep=""),
    Veg1, Veg2, paste("PC", 1:10, sep=""))

nrow(shbg1) #132,972

#Write participant IDs
write.table(shbg1, 
	paste(outdir, "/SHBGxVeg_phenoQC_IDS.M1.txt",sep=""), 
	row.names=FALSE, quote=FALSE)

#Full
write.csv(shbg1,
	paste(outdir, "/SHBGxVeg_pheno.M1.csv", sep=""),
        row.names=FALSE, quote=FALSE)

#Males only
write.csv(shbg1[shbg1$Sex==0,],
        paste(outdir, "/SHBGxVeg_pheno.M.M1.csv", sep=""),
        row.names=FALSE, quote=FALSE)


#Female only
write.csv(shbg1[shbg1$Sex==1,],
        paste(outdir, "/SHBGxVeg_pheno.F.M1.csv", sep=""),
        row.names=FALSE, quote=FALSE)


shbg2<-shbg%>%select(FID, IID, SHBG, logSHBG, SHBG2,
    Age, Sex, Geno_batch,
    paste("center", 1:20, sep=""),
    Veg1, Veg2, paste("PC", 1:10, sep=""),
	Mean_corpuscular_haemoglobin, zblood,
	Menopause, BMI)



#Write participant IDs
write.table(shbg2, 
	paste(outdir, "/SHBGxVeg_phenoQC_IDS.M2.txt",sep=""), 
	row.names=FALSE, quote=FALSE)

#Full
write.csv(shbg2,
	paste(outdir, "/SHBGxVeg_pheno.M2.csv", sep=""),
        row.names=FALSE, quote=FALSE)

#Males only
write.csv(shbg2[shbg2$Sex==0,],
        paste(outdir, "/SHBGxVeg_pheno.M.M2.csv", sep=""),
        row.names=FALSE, quote=FALSE)

#Female only
write.csv(shbg2[shbg2$Sex==1,],
        paste(outdir, "/SHBGxVeg_pheno.F.M2.csv", sep=""),
        row.names=FALSE, quote=FALSE)


#=====================================================
#RINT transform=======================================
#=====================================================

sapply(test, summary)


test<-test%>%select(-Oestradiol, Menopause)


###FULL COHORT
test1<-test%>%select(FID, IID, Test, bioavailableTest, 
    Age, Sex, Geno_batch,
    paste("center", 1:20, sep=""),
    Veg1, Veg2, paste("PC", 1:10, sep="")
    )

test1$RankTest<-RankNorm(test1$Test)
test1$RankbioavailableTest<-RankNorm(test1$bioavailableTest)

test2<-test%>%select(FID, IID, Test, bioavailableTest, 
    Age, Sex, Geno_batch,
    paste("center", 1:20, sep=""),
    Veg1, Veg2, paste("PC", 1:10, sep=""),
    Fasting_time, 
    Mean_corpuscular_haemoglobin, zblood,
    BMI
    )

test2<-test2[complete.cases(test2),]
test2$RankTest<-RankNorm(test2$Test)
test2$RankbioavailableTest<-RankNorm(test2$bioavailableTest)

###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###WRITE OUTPUT=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

outdir="/scratch/mf91122/TestxVeg/pheno"


#Participant list
write.table(test1, 
	paste(outdir, "/TestxVeg_phenoQC_IDS.M1.txt",sep=""), 
	row.names=FALSE, quote=FALSE)

#Full
write.csv(test1,
	paste(outdir, "/TestxVeg_pheno.M1.csv", sep=""),
        row.names=FALSE, quote=FALSE)

#Participant list
write.table(test2, 
	paste(outdir, "/TestxVeg_phenoQC_IDS.M2.txt",sep=""), 
	row.names=FALSE, quote=FALSE)

#Full
write.csv(test2,
	paste(outdir, "/TestxVeg_pheno.M2.csv", sep=""),
        row.names=FALSE, quote=FALSE)




#*Sex Stratify TEST pheno tables-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Re-rank transform
#Male=1, female=2
#Split into sex-stratified groups
#dir.create("phenotables")
models<-list(model0, model1, model2, model3, model4, model5)
modelsM<-list()
modelsF<-list()
for (i in 1:length(models)){
    write.csv(models[[i]], paste("phenotables/T-model", i-1, 
                                 "full-05282021.csv", sep=""),
              quote=FALSE, row.names=FALSE)
    modelsM[[i]]<-models[[i]]%>%filter(Sex==0)
    modelsF[[i]]<-models[[i]]%>%filter(Sex==1)
}

#INT Male tables
for (i in 1:length(modelsM)){
    modelsM[[i]]$RankTest<-RankNorm(modelsM[[i]]$Test)
    modelsM[[i]]$RankbioavailableTest<-RankNorm(modelsM[[i]]$bioavailableTest)
    modelsM[[i]]$RankfreeTest<-RankNorm(modelsM[[i]]$freeTest)
    write.csv(modelsM[[i]], paste("phenotables/T-model", i-1, 
                                  "M-05282021.csv", sep=""),
              quote=FALSE, row.names=FALSE)
}



#Female models 2-5: add menopause column
for (i in 3:6){
    modelsF[[i]]<-phenoQCCauc2%>%select(IID, Menopause)%>% 
        right_join(modelsF[[i]], by="IID")
}

for (i in 1:length(modelsF)){
    modelsF[[i]]<-modelsF[[i]][complete.cases(modelsF[[i]]),]
    modelsF[[i]]$RankTest<-RankNorm(modelsF[[i]]$Test)
    modelsF[[i]]$RankbioavailableTest<-RankNorm(modelsF[[i]]$bioavailableTest)
    modelsF[[i]]$RankfreeTest<-RankNorm(modelsF[[i]]$freeTest)
    write.csv(modelsF[[i]], paste("phenotables/T-model", i-1, 
                                  "F-05282021.csv", sep=""),
              quote=FALSE, row.names=FALSE)
}

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Summarize variables-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

outdir="/scratch/mf91122/TestxVeg/pheno"
new<-as_tibble(read.csv(paste(outdir, "/TestxVeg_pheno.csv", sep=""), header=T))

paste(colnames(new[3:50]), " ",	
        sapply(new[3:50], mean, na.rm=T), "(",
        sapply(new[3:50], sd, na.rm=T), ")", sep="")


SHBG2<-as_tibble(read.csv(paste(outdir, "/SHBGxVeg_pheno.M2.csv", sep=""), header=T)) #Same participants as M2
Test1<-as_tibble(read.csv(paste(outdir, "/TestxVeg_pheno.M1.csv", sep=""), header=T))
Test2<-as_tibble(read.csv(paste(outdir, "/TestxVeg_pheno.M2.csv", sep=""), header=T)) #not same

cols<-c(3,6,7,29,30,41:44)

paste(colnames(SHBG2[cols]), " ",
        sapply(SHBG2[cols], mean, na.rm=T), "(",
        sapply(SHBG2[cols], sd, na.rm=T), ")", sep="")

cols<-c(3:6,28,29)

paste(colnames(Test1[cols]), " ",
        sapply(Test1[cols], mean, na.rm=T), "(",
        sapply(Test1[cols], sd, na.rm=T), ")", sep="")

cols<-c(3:6,28,29,41,43)

paste(colnames(Test2[cols]), " ",
        sapply(Test2[cols], mean, na.rm=T), "(",
        sapply(Test2[cols], sd, na.rm=T), ")", sep="")
