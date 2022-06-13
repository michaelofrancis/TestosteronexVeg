library(gridExtra)
library(png)
library(grobblR)
library(qpdf)
library(pdftools)





#Order files---------

order_files_nopop<-function(pngfiles){
    #Use:
    #fileorder<-order_files(pngfiles)
    fileorder<-NULL
    for (i in 1:length(pheno)){
        found<-pngfiles[grep(pheno[i], pngfiles)]
        fileorder<-c(fileorder, found)
    }
    return(fileorder)
}


# Testosterone Manhattan qq plots --------------------------------------
dir<-"/Users/mike/Documents/Research/TestxVeg/plot"
pngfiles<-list.files(dir,full.names =  T)[grep(".png", list.files(dir))]
pheno<-c("logSHBG", "SHBG", "RankTest", "RankbioavailableTest")
#order files
fileorder<-order_files_nopop(pngfiles)
fileorder<-unique(fileorder)
qqpng<-fileorder[grep(".qq.",fileorder)]
qqpng
manpng<-fileorder[!fileorder %in% qqpng]
manpng
page<-list()
rows=4
bord=F
numpages<-ceiling(length(manpng)/rows)
#%% = remainder
remainder<-length(manpng) %% rows

#at 
floor(length(manpng)/rows)

for (i in 1:(floor(length(manpng)/rows))){
    page[[i]]<-grob_layout(
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-1)]), 
            grob_col(
                qqpng[i*rows - (rows-1)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-2)]), 
            grob_col(
                qqpng[i*rows - (rows-2)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-3)]), 
            grob_col(
                qqpng[i*rows - (rows-3)], 
                border=bord)),
        grob_row(
            grob_col(
                border=bord, 
                manpng[i*rows - (rows-4)]), 
            grob_col(
                qqpng[i*rows - (rows-4)], 
                border=bord))
        #height = 100,
        #width = 100,
        #padding = 0
    ) 
}
# #Do the last page
# i=ceiling(length(manpng)/rows)
# page[[i]]<-grob_layout(
#     grob_row( #1 row
#         grob_col(
#             border=bord, 
#             manpng[i*rows - (rows-1)]), 
#         grob_col(
#             qqpng[i*rows - (rows-1)], 
#             border=bord)),
#     grob_row( #2 rows
#         grob_col(
#             border=bord, 
#             manpng[i*rows - (rows-2)]), 
#         grob_col(
#             qqpng[i*rows - (rows-2)], 
#             border=bord)),
#     grob_row( #3 rows
#         grob_col(
#             border=bord, 
#             NA)),
#     grob_row( #4 rows
#         grob_col(
#             border=bord, 
#             NA))
#     
# )


dir.create(paste(dir, "/pdf", sep=""))
outdir<-paste(dir, "/pdf/", sep="")

grob_to_pdf(
    page,
    file_name = paste(outdir, "Testosterone.ManhattanQQ.pdf", sep=""),
    add_page_numbers = F,
    meta_data_title = "PUFA-GWAS supplementary fig PDF"
)
