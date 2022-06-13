library(tidyverse)

tab<-as_tibble(read.csv("../TestxVeg/pheno/06112022/old/TestxVeg_pheno.csv"))
tab

pheno<-c("SHBG", "Test", "bioavailableTest")
phe2<-tab%>%select(all_of(pheno))
phe2<-phe2[complete.cases(phe2),]

A<-cov(phe2)
ev <- eigen(A)
values <- ev$values
values

top<-(sum(values))^2
bottom<-sum(values^2)

top/bottom
#[1] 1.097479

5e-08 / 1.097479
#[1] 4.555896e-08
