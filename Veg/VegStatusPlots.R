library(plyr)
library(dplyr)
library(tidyverse)

#number1
vegqc1<-read.delim("vegQC1_06102022.txt", header=TRUE, sep="\t")
vegqc1<-as_tibble(vegqc1)

vegqc2<-read.delim("vegQC2_04032021.txt", header=TRUE, sep="\t")
vegqc2<-as_tibble(vegqc2)

vegqc2<-inner_join(vegqc1,vegqc2)

colnames(vegqc1)<-c("IID", "Veg1", "Mixed")
vegqc1$Veg1[vegqc1$Veg1==0]<-"No"
vegqc1$Veg1[vegqc1$Veg1==1]<-"Yes"

vegqc1$Mixed[vegqc1$Mixed==0]<-"Consistent"
vegqc1$Mixed[vegqc1$Mixed==1]<-"Inconsistent"

vegqc1$Veg1<-as.factor(vegqc1$Veg1)
vegqc1$Mixed<-as.factor(vegqc1$Mixed)

#number2--------------
vegqc2<-vegqc2%>%select(-Consistent_Self_Reported_Vegetarian_across_all_24hr)
vegqc2
colnames(vegqc2)<-c("IID", "Veg1", "Mixed", "Veg2")
vegqc2$AteMeat<-0
vegqc2$AteMeat[vegqc2$Mixed==0 & vegqc2$Veg1==1 & vegqc2$Veg2==0]<-1

vegqc2$ReasonNoVeg<-""
vegqc2$ReasonNoVeg[vegqc2$AteMeat==1]<-"Ate Meat"
vegqc2$ReasonNoVeg[vegqc2$Mixed==1]<-"Inconsistent Veg self-identification"

vegqc2$Veg2<-as.factor(vegqc2$Veg2)
vegqc2$ReasonNoVeg<-as.factor(vegqc2$ReasonNoVeg)
vegqc2

# for full set of 24HR participants -----------------------------------------

ggplot(vegqc1, aes(Veg1)) + 
    geom_bar(aes(fill=forcats::fct_rev(Mixed)), 
             position="stack")+
    xlab("Special diet indicated as vegetarian and/or vegan
         across 24HR instances (AKA 'Veg1')")+
    ylab("Number of participants")+
    ggtitle("UK Biobank all 24HR participants Veg1 Status")+
    labs(fill = "Consistent special 
    diet response")+
    theme_bw()+
    stat_count(geom = "text", colour = "black", size = 3,
               aes(label = ..count..),
               position=position_stack(vjust=1.05))


ggplot(vegqc2, aes(Veg2)) + 
    geom_bar(aes(fill=forcats::fct_rev(ReasonNoVeg)) ,position="stack")+
    xlab("Self-reported vegetarian and/or vegan
    who also didn't indicate eating meat AKA 'Veg2'")+
    ylab("Number of participants")+
    ggtitle("UK Biobank all 24HR participants Veg2 Status")+
    labs(fill = "Reason for veg disqualification")+
    theme_bw()+
    stat_count(geom = "text", colour = "black", size = 3,
               aes(label = ..count..),
               position=position_stack(vjust=1.05))


# for participants in our analysis only -------------------------
#COPY all the code from above and just change to the list of participants for each analysis

#Lipids list of participants

lip<-read.table("../LipidsxVeg/pheno/LipidsxVeg_phenoQC_IDS.txt", header=T)
nrow(lip) #[1] 155346
#Looking for 147903..?
#Make participant characteristics tables first to get the numbers clear

#Testosterone list of participants


#number1
vegqc1<-read.delim("vegQC1_06102022.txt", header=TRUE, sep="\t")
vegqc1<-as_tibble(vegqc1)

vegqc2<-read.delim("vegQC2_04032021.txt", header=TRUE, sep="\t")
vegqc2<-as_tibble(vegqc2)

vegqc2<-inner_join(vegqc1,vegqc2)

colnames(vegqc1)<-c("IID", "Veg1", "Mixed")
vegqc1$Veg1[vegqc1$Veg1==0]<-"No"
vegqc1$Veg1[vegqc1$Veg1==1]<-"Yes"

vegqc1$Mixed[vegqc1$Mixed==0]<-"Consistent"
vegqc1$Mixed[vegqc1$Mixed==1]<-"Inconsistent"

vegqc1$Veg1<-as.factor(vegqc1$Veg1)
vegqc1$Mixed<-as.factor(vegqc1$Mixed)

#number2--------------
vegqc2<-vegqc2%>%select(-Consistent_Self_Reported_Vegetarian_across_all_24hr)
vegqc2
colnames(vegqc2)<-c("IID", "Veg1", "Mixed", "Veg2")
vegqc2$AteMeat<-0
vegqc2$AteMeat[vegqc2$Mixed==0 & vegqc2$Veg1==1 & vegqc2$Veg2==0]<-1

vegqc2$ReasonNoVeg<-""
vegqc2$ReasonNoVeg[vegqc2$AteMeat==1]<-"Ate Meat"
vegqc2$ReasonNoVeg[vegqc2$Mixed==1]<-"Inconsistent Veg self-identification"

vegqc2$Veg2<-as.factor(vegqc2$Veg2)
vegqc2$ReasonNoVeg<-as.factor(vegqc2$ReasonNoVeg)
vegqc2

# for full set of 24HR participants -----------------------------------------

ggplot(vegqc1, aes(Veg1)) + 
    geom_bar(aes(fill=forcats::fct_rev(Mixed)), 
             position="stack")+
    xlab("Special diet indicated as vegetarian and/or vegan
         across 24HR instances (AKA 'Veg1')")+
    ylab("Number of participants")+
    ggtitle("UK Biobank all 24HR participants Veg1 Status")+
    labs(fill = "Consistent special 
    diet response")+
    theme_bw()+
    stat_count(geom = "text", colour = "black", size = 3,
               aes(label = ..count..),
               position=position_stack(vjust=1.05))


ggplot(vegqc2, aes(Veg2)) + 
    geom_bar(aes(fill=forcats::fct_rev(ReasonNoVeg)) ,position="stack")+
    xlab("Self-reported vegetarian and/or vegan
    who also didn't indicate eating meat AKA 'Veg2'")+
    ylab("Number of participants")+
    ggtitle("UK Biobank all 24HR participants Veg2 Status")+
    labs(fill = "Reason for veg disqualification")+
    theme_bw()+
    stat_count(geom = "text", colour = "black", size = 3,
               aes(label = ..count..),
               position=position_stack(vjust=1.05))
