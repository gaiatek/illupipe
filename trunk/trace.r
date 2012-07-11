###############################################################################
# illutrace.R  v.3.0 beta                                                     #
#-----------------------------------------------------------------------------#
# Illumina MicroArray Pipeline History                                        #
#                                                                             #
# This is the trace of function calls to the pipeline for a given project     #
# for illumina beadArrays transcriptomics.                                    #
#                                                                             #
# This file is an history of the data analysis done in this working directory.#
# N.B. : It is recommended to save workspace before interrupting the pipeline.#
#-----------------------------------------------------------------------------#
# @author : Jean-Philippe Goulet                                              #
# @date   : November 23rd 2010                                                #
###############################################################################
# source 'illupipe.r' and 'illupipe_report.r' first
# Source the pipeline functions :
source("/projects/GeneExpressionAnalysis/trunk/illupipe.r") 
#source("../GeneExpressionAnalysis/illupipe_ComBat.R")

# load libraries and parameters 
init() # use if you just reloaded your workspace.

param = getParams("params.txt")


## Define the data directory for the present version
## Make sure you copied the trunk/data template folder into your working data directory
## AND copied the appropriate 'probes.txt' file for your chip version!
dirData = "data" 
# Where the data is stored. Could be different for different versions
# previously : dirData = getwd()

###################################
## annotation package
###################################
# eg. From V3
annDB = "MyIlluminaHumanRef8V3.db"
library(annDB, character.only=T)

# import raw data
eset.raw <- readBeads(chipsDir='chips')

################ Load the data #################
## Load the raw data (eventually will require a function to read it)
eset.raw <- load.raw(dirData,replace.zeroes=NA)
#eset.raw = load.raw(third.id="DonorID",replace.zeroes=NA) # Should there be  a third table in worksheet.xls
#

#################### QC  #################################################################

## Diagnostics plots. (eventually ggplot2?)
qcDiagnose(eset.raw)

## Supporting figures for outlier removal (might have to pre-process first...)
temp = eset.raw
exprs(temp) = ppNormalize(exprs(temp),bg=80, method="quantile") 
temp = ppFilter(temp,"bg"=80, "bgCount"=2,"iqr"=0.3,report=FALSE)
unac = unaCluster(temp,labels=pData(temp)$Tissue,report=FALSE, width=1500, height=800, cex=1,ret=TRUE)
pdf(file=file.path(dirData,"outliers_supporting_figures.pdf"),width=40,height=20)
plot(unac$cor,labels=paste(pData(temp)$ArrayID,pData(temp)$ChipType,pData(temp)$Tissue,pData(temp)$Vaccine,sep="|"))
dev.off()


#PostPP
#qcPostPP(temp,dest.dir="data/diagnostic_plot/postprocessing/")



## Remove outliers
eset.raw <- qcRemoveOutliers(eset.raw,verbose=TRUE)


#################### Pre-processing (Could be merged in one wrapper function)  #############
eset <- eset.raw

# Look at batch profile
#table(pData(eset)$ChipType[pData(eset)$Tissue=="LNa"],pData(eset)$Vaccine[pData(eset)$Tissue=="LNa"])

# Edit parameters for filtering.
# Please put your initials in the "author" slot and set the project name in slot "project"!
param <- list("project"="Project Name", "author"="NA","bg"=80, "bgCount"=2,"normalizationMethod"="quantile","iqr"=0.3,"author"="FL"
,"TimeStamp"=Sys.time(),"batch"="h.Set","batch.covariates"=c()) #
attr(eset,"param") <- param   


# Normalization
exprs(eset) <- ppNormalize(exprs(eset))

# Batch correction
eset <- ppComBat(eset,batch=param$batch,covariates=param$batch.covariates)

# Technical Replicates
#eset = ppAverageTechnicalReplicates(eset)


# Non-specific filtering
save(eset,file=file.path(dirData,"eset_unfiltered.Rdata"))
eset <- ppFilter(eset)

# Finally, write the pre-processed eset to the dirData
save(eset,file=file.path(dirData,"eset.Rdata"))



#retain only immunological genes
#immune.genes = unique(read.delim("http://amigo.geneontology.org/cgi-bin/amigo/term-assoc.cgi?term=GO:0006955&format=go_assoc&session_id=7709amigo1287769728",header=FALSE)$V3)
#eset = eset[fData(eset)$SYMBOL %in% immune.genes,]
#################### Pre-processing (Could be merged in one wrapper function)  #############

## Exploratory Analysis figures
mds <- as.data.frame(cmdscale(euc(eset)))
mds <- cbind(mds,pData(eset))
p   <- qplot(x=V1,y=V2,data=mds,color=Tissue,shape=FullVaccine,facets=~Vaccine)
pdf(file=file.path(dirData,"exploratory_analysis.pdf"),height=10,width=15)
print(p)
dev.off()




######################### Linear Modeling for differential Gene Expression ################################


# Initialize the fits object
fits <- list()
## Linear model 1 : Stim*Vaccine*Time : GIVE AN example of two contrast sets, and two different models. ALSO continuous variable dummy
fits2 <- list()
variables   <- c("Tissue","Vaccine")
model.name  <- paste(variables,collapse=", ")
fit         <- dgLimmaCategoricalFit(eset,variables)    # fit the model

# Comparisons to Liver
contrasts   <- dgANOVAContrasts(eset,variables,"Vaccine","unchallenged","Tissue","LV")[["Vaccine:Tissue"]]# contrasts should be renamed here
fits2[["LV, unchall as baselines, LAV only,no LNm"]]  <- dgLimmaContrastsFit(fit,contrasts[!grepl("LNm",contrasts)&!grepl("CMV",contrasts)]) 
fits2[["LV, unchall as baselines, LAV only"]]         <- dgLimmaContrastsFit(fit,contrasts[!grepl("CMV",contrasts)])
fits2[["LV, unchall as baselines, no LNm"]]           <- dgLimmaContrastsFit(fit,contrasts[!grepl("LNm",contrasts)])  
fits2[["LV, unchall as baselines"]]                   <- dgLimmaContrastsFit(fit,contrasts) 



# Comparisons to unchallenged
#contrasts   = unique(c(dgANOVAContrasts(eset,variables,"Vaccine","unchallenged","Tissue","LNi")$Vaccine
#		,dgANOVAContrasts(eset,variables,"Vaccine","LAV","Tissue","LV")$Vaccine ))# contrasts should be renamed here
#contrasts = sort(contrasts[!grepl('unchallenged-',contrasts)])
#contrasts.by.tissue = sapply(unique(pData(eset)$Tissue),function(tis)contrasts[grepl(tis,contrasts)],simplify=FALSE)
#fits2[["Comparing to unchallenged, LAV only"]]       = dgLimmaContrastsFit(fit,contrasts[!grepl('CMV',contrasts,perl=TRUE) & grepl('unchallenged$',contrasts,perl=TRUE)])
#fits2[["Comparing to unchallenged, all tissues"]]       = dgLimmaContrastsFit(fit,contrasts[grepl('unchallenged$',contrasts,perl=TRUE)])  
#fits2 = c(fits2,sapply(contrasts.by.tissue,function(cts)dgLimmaContrastsFit(fit,cts),simplify=FALSE))

 
#Finally, save fits to an object. List names will be used as aliases for the model name and contrast sets
fits[[model.name]] <- list("fit"=fit,"fits2"=fits2)
save(fits,file=file.path(dirData,"fits.Rdata"))






################ PART 2 : Reporting #################
dirReport  <- "Report" # Where the report will be created. Could be different.
dirData    <- "data" # Where the source data for the report is stored

create.report(project.name=param$project,report.author=param$author,third.id="DonorID",report.dir = dirReport, data.dir=dirData
			,annots.items = list("raw"=TRUE,"samples"=TRUE,"arrays"=TRUE,"probes"=TRUE,"gene.sets"=TRUE)
			,removal.items= list("removed.arrays"=TRUE,"missing.samples"=TRUE)
			,preprocessing.items = TRUE
			,exploratory.items   = TRUE 
			 ,dgea.items=list("Tissue, Vaccine"=names(fits2))                                   
			#,dgea.items=list("Tissue, Vaccine"=c("LV, unchall as baselines","Comparing to Liver, LAV only","Comparing to Liver, CMV only","Comparing to Liver, unchallenged only","Comparing to Liver, all vaccines","Comparing to unchallenged, LAV only","Comparing to unchallenged, all tissues")) 
			,dgea.params = list(n.genes.toplists=500,n.genes.heatmaps=50,heatmaps.alt.text=c("DonorID","FullVaccine","Tissue"),n.genesets.toplists=100,p.thresh=0.05,q.thresh=0.05,heatmap.bl.subtract=FALSE) #
)

finalize.report(version="preliminary",desc="Preliminary report ...", 

stop("End of analysis")


### Heatmap of genes from selected pathways
# Load unfiltered eset
# Load the gene and pathway memberships
# subset the eset, plot the columns colored accoring to class




