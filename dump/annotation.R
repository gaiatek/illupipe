
## Annotation packages generation 
#---------------------------------------------
# To annotate a probeset specific to a microArray, you need first to load AnnotationDbi
library(AnnotationDbi)
# You then have to make sure you have the DB schema for the species you create the annotation for 
source("http://bioconductor.org/biocLite.R")
biocLite("mouse.db0")
biocLite("human.db0")
biocLite("rhesus.db0")

# The input used by AnnotationDbi is a tab-delimited text where you have two columns :
# First column for the probeID (or BeadID)
# Second column for the RefSeq accession number
# Sometimes it works best if you leave out the decimal after NM_xxxxx.x and only keep NM_xxxxx
# The decimal part stand for the revision number used for a transcript and only the latest one is kept
# in the database. Keeping it may result in <NA> for a transcript with a newer revision.

# This input file is the 'filename' argument of the fucntion makeDBPackage. 
# It is good practice to keep them in you annotation directory (param$dirCommon), if you want to 
# update them later.



###################################################
### chunk number 7: SQLForge  Human
###################################################
library(AnnotationDbi)

makeDBPackage("HUMANCHIP_DB",
              affy = FALSE,
              prefix = "MyIlluminaHumanRef8V2",
              fileName = paste(param$dirCommon,"/","HumanRef-8_V2_probe2RefSeq.csv",sep=""),
              baseMapType = "refseq",
              outputDir = param$dirCommon,
              version="3.0.0",
              manufacturer = "Illumina",
              chipName = "Illumina Human Ref8V2 Array",
              manufacturerUrl = "http://www.illumina.com")


makeDBPackage("HUMANCHIP_DB",
              affy = FALSE,
              prefix = "MyIlluminaHumanRef8V3",
              fileName = paste(param$dirCommon,"/","HumanRef-8_V3_probe2RefSeq.csv",sep=""),
              baseMapType = "refseq",
              outputDir = param$dirCommon,
              version="3.0.0",
              manufacturer = "Illumina",
              chipName = "Illumina Human Ref8V3 Array",
              manufacturerUrl = "http://www.illumina.com")


# Human HT12 V4
library(AnnotationDbi)
temp = read.delim("HumanHT-12_V4_0_R1_15002873_B.txt")
#temp$Accession  = gsub('\\..$','',temp$Accession,perl=TRUE)
writeLines(paste(temp$ProbeID,temp$Accession,sep="\t"),"temp.txt")
makeDBPackage("HUMANCHIP_DB",
              affy = FALSE,
              prefix = "MyIlluminaHumanHT12V4",
              fileName = "temp.txt", # Uses Accession column from HumanHT-12_V4_0_R1_15002873_B.txt
              baseMapType = "gbNRef",
              outputDir = "/projects/annotations",
              version="2.0.0",
              manufacturer = "Illumina",
              chipName = "Illumina Human HT-12 V4",
              manufacturerUrl = "http://www.illumina.com")



###################################################
### chunk number 7: SQLForge Mouse
###################################################
dirCommon="/home/jp/R/annotations"
library(AnnotationDbi)
makeDBPackage("MOUSECHIP_DB",
    affy=FALSE,
    prefix="MyIlluminaMouseRef8V2",
    fileName=paste(param$dirCommon,"/","MouseRef-8_V2_probe2RefSeq.csv",sep=""),
    baseMapType="gb",
    outputDir = param$dirCommon,
    version="4.0.0",
    manufacturer = "Illumina",
    chipName = "Illumina Mouse Ref8V2 Array",
    manufacturerUrl = "http://www.illumina.com")



###############################################################################
# Affymetrix HG-U133A to Human
###############################################################################

# What is packages AnnotationDbi, matchprobes, altcdf....
library(AnnotationDbi)
library(altcdfenvs)
dirCommon="/home/lefebvrf/R_libraries/myannotations"

# 1) match probes to NM
# 2) build your own cdf with some fancy new approaches
# 3) create anoation package with other fancy annotations other than the gene




###############################################################################
# Library Building, with Lumi data as a starting point
###############################################################################


## INTRO
#
# Although building BioC packages provides up-to-date annotations, it is still dependent on
# Illumina's annotation files as a seed. The situation of ILMN_1798181 (IRF7) provides
# an example where the RefSeqID provided by Illumina has been retired and an important
# probe is thus left un-annotated by AnnotationDBI.
# Here we take advantage of the the lumi data which involves re-blasting the probe
# sequences on RefSeq to get more accurate RefSeqID. Here who we proceeed:
# 1- Load Illumina annotation file. This file is downloaded from Illumina's website. Rows
#   other than probes must be manually deleted prior to loading.
# 2- Load the lumi annations.
# 3- Cross those two annotation to finally get a table with : the probe address and
#   lumi's refseq
# 4- Write the resulting table to temp.txt
# 5- Re-build package as usual.
illu.annot  = read.delim("/projects/annotations/fromIllumina/HumanRef-8_V2_0_R4_11223162_A.txt",stringsAsFactors=FALSE)[,c("RefSeq_ID","Probe_Id","Array_Address_Id","Probe_Sequence")] # Extracted the zip from Illumina WS and left only probe rows
lumi.annot = read.delim("/projects/annotations/fromLumi/all_human_probes_mapped.txt",header=FALSE,stringsAsFactors=FALSE)[,c("V8","V16")]
annot = merge(illu.annot,lumi.annot,by.x="Probe_Sequence",by.y="V8",all.x=TRUE,all.y=FALSE)[,c("Array_Address_Id","V16")]
write.table(annot,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t",file="temp.txt")
makeDBPackage("HUMANCHIP_DB",
              affy = FALSE,
              prefix = "MyIlluminaHumanRef8V2",
              fileName = "/projects/annotations/ourPackages/temp.txt",
              baseMapType = "refseq",
              outputDir = "/projects/annotations/ourPackages/",
              version="1.0.0",
              manufacturer = "Illumina",
              chipName = "Illumina Human Ref8V2 Array",
              manufacturerUrl = "http://www.illumina.com")

