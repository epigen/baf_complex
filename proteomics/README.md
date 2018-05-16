# R code to analyze the proteomics data

* the data to analyze is taken from Proteome Discoverer export (as tab separated files with R friendly headers)
* myCode.R contains the code to analyze and visualize the data. Examine the log of the R session history in the file Rhistory to see how the data were analyzed. The data exported from PD are loaded into dataframes as follows:

    BRG1.rr2 <- read.table("CC50_Set1-3_BRG1-IP_Arid-Norm_NonMerged-TMT_Proteins.txt",header=T,stringsAsFactors=F)
    ARID1A.rr2 <- read.table("CC50_Set1-3_A1A-IP_Arid-Norm_NonMerged-TMT_Proteins.txt",header=T, stringsAsFactors=F)

    ARID1A.pep2 <- read.table("CC50_Set1-3_A1A-IP_Arid-Norm_NonMerged-TMT_PeptideGroups.txt",header=T, stringsAsFactors=F)
    BRG.pep2 <- read.table("CC50_Set1-3_BRG1-IP_Arid-Norm_NonMerged-TMT_PeptideGroups.txt",header=T, stringsAsFactors=F)

    BRG1.ab2 <- read.table("CC50_Set1-3_BRG1_Precursor_NoNorm_Proteins.txt",header=T,stringsAsFactors=F)
    ARID1A.ab2 <- read.table("CC50_Set1-3_A1A_Precursor_NoNorm_Proteins.txt",header=T,stringsAsFactors=F)
