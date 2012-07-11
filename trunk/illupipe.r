###############################################################################
# illupipe.R  v.2.0.0.3                                                       #
#-----------------------------------------------------------------------------#
# Illumina MicroArray Pipeline Toolkit                                        #
# Functions for data import, quality control, preprocessing,                  #
# unsupervised analysis, differential gene expression analysis and            #
# pathway analyses of illumina beadArrays transcriptomics                     #
#-----------------------------------------------------------------------------#
# @author : Jean-Philippe Goulet                                              #
# @date   : January 25th  2011                                                #
###############################################################################
# must include file for report functions
# source("illupipe_report.r")
# please include script for batch effect correction if needed :
# source("illupipe_ComBat.R")


# init : load all needed libraries for any task in the pipeline
# N.B. : you still need to load the gene annotation library 
# eg. library("myIlluminaHumanRef8V3.db") value stored in annDB
init <- function() {
	
	library(limma)
	library(gplots)
	library(genefilter)
	library(bioDist)
	library(annaffy)
	library(annotate)
	library(R2HTML)
	library(beadarray)
	library(RColorBrewer)
	library(latticeExtra)
	library(affyQCReport)
	library(statmod)
	library(ggplot2)
	library(Vennerable)
	options(stringsAsFactors=FALSE)
	library(gdata)
	library(hwriter)
	library(WriteXLS)
	library(SortableHTMLTables)
	library(multtest)
	library(impute)
}
init()


###################
## Data Import
###################
# readBeads : extracts raw expression data from Illumina beadArrays contained
#             in 'chips' directory. At least the .txt files for every arrays 
#             imported is needed, plus the Metrics.txt file  
#             for every chip. Some Quality Control diagnostics plots are
#             generated. They are place in a folder named "diagnostic_plot".
#             The chips must be placed under the "chips" folder in the working
#             directory to be read.
#------------------------------------------------------------------------------
# Imports a set of chips from the same beadArray version only (cf.chipsVersion) 
#------------------------------------------------------------------------------
# vgti         : logical. If true, will use "_perBeadFile.txt" instead of ".txt"
# impute       : logical. Whether missing values need to be imputed using KNN.
# chipsDir     : character string. Where are located the chips raw data.
# chipsVersion : character string. What version of beadArrays is imported.
#                (i.e. either "Humanv4","Humanv3","Humanv2" or "Mousev2")
# chipsSet     : optional. A string describing the hybridisation set. Useful
#                when you already have a 'eset_raw.Rdata' and want to avoid to
#                overwrite it. (cf. output) (eg. "batch2") 
# data.dir     : character string. Where to put the diagnostic plots.(eg.'data')
# report       : logical. Deprecated. 
# output       : an ExpressionSet of raw data. Output object is saved as 
#                "eset_raw_<chipsSet>.Rdata" in the data directory.
readBeads <- function(chipsDir='chips', chipsVersion="Humanv4", vgti=T, 
impute=TRUE, report=F, chipsSet=NULL,data.dir=dirData) {
	
	# Get chip list from dirData/chips
	chip = dir(chipsDir) 
	
	# Retrieve control probes information
	data(ExpressionControlData)
	chipControl = ExpressionControlData[[chipsVersion]]
	qclist = chipControl$Array_Address
	
	# read gene annotation file before reading beads (in case there\'s an error)
	gene.ann    <- getProbesAnn(data.dir)
	
	
	##########################################################################
	## Read Data                                                             #
	##########################################################################
	# - reads the beadArrays .txt files in the "chips/<barcode>" directories #
	# - appends expression values to dat                                     #
	# - generates phenotypic annotation to pheno                             #
	# - produces QC plots and reports for each chip                          #
	#------------------------------------------------------------------------#
	if(!file.exists(file.path(data.dir, "diagnostic_plot"))) {
		dir.create(file.path(data.dir, "diagnostic_plot"))
	}
	if(vgti) {
		textType <- "_perBeadFile.txt"
	}
	else {
		textType <- ".txt"
	}
	dat <- NULL
	
	# Private functions to control the bead data summarization
	.myMedian <- function(x) median(x, na.rm = TRUE) 
	#.mySe <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))
	.mySd <- function(x) sd(x, na.rm = TRUE)
	.myGreenChannelTransform <- function(BLData, array){
		x = getBeadData(BLData, array=array,what="Grn")
	}
	
	greenChannel = new("illuminaChannel", .myGreenChannelTransform,
					   illuminaOutlierMethod, .myMedian, .mySd, "G")
	
	
	for(i in 1:length(chip)) {
		# Check and remove unwanted text files 
		chipFiles     <- dir(file.path(chipsDir,chip[i]))
		filesToRemove <- chipFiles[grep('beadType',chipFiles)]
		if(length(filesToRemove) > 0) {
			for(j in 1:length(filesToRemove)) {
				system(paste('rm ', file.path(chipsDir,chip[i],filesToRemove[j])))
			}
		}
		filesToRemove <- chipFiles[grep('_qc.txt',chipFiles)]
		if(length(filesToRemove) > 0) {system(paste('rm ', file.path(chipsDir,chip[i],filesToRemove[1])))}
		
		# To switch probeIDs from ArrayAddressIDs to ILMN_IDs : 
		# put illuminaAnnotation=chipsVersion in readIllumina()
		BLData <- readIllumina(dir = file.path(chipsDir,chip[i]),
							   useImages=F, metrics=T, illuminaAnnotation=chipsVersion)
		
		an     <- sectionNames(BLData)
		
		BSData <- beadarray::summarize(BLData, channelList=list(greenChannel))
		featureNames(BSData) = fData(BSData)$ArrayAddressID
		
		if (i == 1) { 
			dat = exprs(BSData) 
			pheno = pData(BSData)
		}
		else {
			dat = cbind(dat,exprs(BSData))
			pheno = rbind(pheno,pData(BSData))
		}
		## Quality Control
		##---------------------------------------------------------------------
		# MA plot
		png(file=file.path(data.dir, "diagnostic_plot",paste("maXY_",chip[i],".png",sep="")),
			width = 800, height = 800)
		plotMAXY(exprs(BSData), arrays = 1:dim(BSData)[[2]],pch = 16,
				 labels=substr(an,12,12), main=chip[i])
		dev.off()
		
		# Quality Controls
		# STEP 1 : generate qc_info
		#----------------------------------------------------------------------
		# select existing control probes from chip
		control = chipControl[qclist %in% fData(BSData)$ArrayAddressID,1:3]
		qcLabel = as.character(control$Array_Address)
		qcGroup = control$Reporter_Group_Name
		# quality control probes summary for chip[i]
		qc = cbind(exprs(BSData)[qcLabel,],qcGroup)
		
		qcInfo <- matrix(ncol=ncol(qc)-1,nrow=nlevels(qcGroup))
		rownames(qcInfo) = levels(qcGroup)
		colnames(qcInfo) = colnames(exprs(BSData))
		
		for(j in 1:nlevels(qcGroup)) {
			for(k in 1:ncol(qc)-1){
				qcInfo[j,k] <- mean(qc[grep(j,qc[,ncol(qc)]),k])
			}
		}
		qc = qc[,1:ncol(qc)-1,drop=FALSE]		
		write.table(qcInfo, file=file.path(data.dir, "diagnostic_plot",paste(
																			 chip[i],"_gene_qcinfo.txt",sep="")), sep="\t")
		
		# STEP 2 : generate hybridization control plots
		#----------------------------------------------------------------------
		color      <- rainbow(ncol(qc))
		
		# Probes intensities Distribution accross arrays.
		# this section is a test
		#------------------------------------------------
		qcGroup = control$Reporter_Group_id
		qcGroup = sub('-.....-..','', qcGroup)
		qcGroup = as.factor(qcGroup)
		png(file=file.path(data.dir, "diagnostic_plot", paste("qcProbes_", 
															  chip[i],".png",sep="")), width = 1024, height = 1024)	
		par(mfrow=c(3,3))
		# mouse arrays have 12 controls instead of 9
		if(chipsVersion == "Mousev2" || chipsVersion == "Humanv4") par(mfrow=c(4,3))
		
		for(j in 1:nlevels(qcGroup)) {
			boxplot(log2(qc[which(qcGroup == levels(qcGroup)[j]),])~col(qc)
					[which(qcGroup == levels(qcGroup)[j]),], main=levels(qcGroup)[j],
					ylim=c(4,16),col=color)
			
		}
		dev.off()
	# end of test section.
		
	}
	
	##########################################
	## Create ExpressionSet out of read data
	##########################################
	dat = dat[rownames(dat)%in%rownames(gene.ann),] 
	# dat[rownames(gene.ann),]
	pheno = pheno[colnames(dat),]
	
	#pheno$Class <- factor(pheno$Class)
	pheno <- new("AnnotatedDataFrame", data = pheno,
				 varMetadata = data.frame(labelDescription = colnames(pheno)))
	feature <- new("AnnotatedDataFrame", data = gene.ann[rownames(dat),],
				   varMetadata = data.frame(labelDescription = colnames(gene.ann)))
	eset.raw  <- new("ExpressionSet", exprs= dat, phenoData=pheno, 
				 featureData=feature, annotation=chipsVersion)
	fData(eset.raw)$ProbeID = rownames(exprs(eset.raw)) #To avoid NAs
	
	# KNN imputation
	if(impute){
		exprs(eset.raw) = impute.knn(exprs(eset.raw))$data
	}
	if(is.null(chipsSet)) {
		esetName = "eset_raw.Rdata"
	}
	else {
		esetName = paste("eset_raw_", chipsSet, ".Rdata")
	}
	save(eset.raw,file=file.path(data.dir,esetName))
	# output the expressionSet
	eset.raw
}


# load.raw : Loads the raw data object. The resulting eset will contain probes 
#            that are also in the probe annotation and reannotate the raw eset.
# data.dir       : Directory where is stored the eset_raw.Rdata
# strictHybs     : cf. annotate.eset
# third.id       : cf. annotate.eset
# replace.zeroes : A value to use in order to impute expression value instead of keeping 0.
# impute         : whether NAs should be imputed by KNN (bioCpackage impute, default arguments)
load.raw <- function(data.dir=dirData,third.id=NULL,strictHybs=TRUE,replace.zeroes=NA,impute=TRUE) #,trim.controls=TRUE
{
	load(file.path(data.dir,"eset_raw.Rdata")) # This loads an object names eset.raw
	eset.raw = annotate.eset(eset.raw,data.dir=data.dir,third.id=third.id,strictHybs=strictHybs)
	#if(trim.controls)
	#{
    #	#attr(eset.raw,"Controls") = eset.raw[fData(eset.raw)$Control,]
    #	eset.raw = eset.raw[!fData(eset.raw)$Control ,]
    #	fData(eset.raw)$Control = NULL
    #}
	
    # Correct for zeroes (which makes no sense). Replace them by something more meaningful
	if(!is.null(replace.zeroes)){
		exprs(eset.raw)[which( exprs(eset.raw)<=0  ,arr.ind=TRUE)] = replace.zeroes
	}
	
    # KNN imputation
	if(impute & any(is.na(exprs(eset.raw)))){
		exprs(eset.raw) = impute.knn(exprs(eset.raw))$data
	}
	
	return(eset.raw)	
}

# This function loads the eset object, annotating it anew along the way.
load.eset <-function(data.dir=dirData,strictHybs=FALSE,third.id=NULL,unfiltered=FALSE)
{
	if(unfiltered){
		load(file.path(data.dir,"eset_unfiltered.Rdata"))
	}else{
		load(file.path(data.dir,"eset.Rdata"))
	}
	eset  = annotate.eset(eset,data.dir,third.id=third.id,strictHybs=FALSE)
    # fData(eset)$Control = NULL
	return(eset)
}


# This function loads the fits object, annotating its fit$genes anew along the way.
load.fits <-function(data.dir=dirData)
{
	load(file.path(data.dir,"fits.Rdata"))
	probesAnn = getProbesAnn(data.dir)
	for(model in names(fits))
	{
		fit               = fits[[model]]$fit
		fit$genes         = probesAnn[rownames(fit$genes),]
		fits[[model]]$fit = fit
		for(contset in names(fits[[model]]$fits2))
		{
			fit2                           = fits[[model]]$fits2[[contset]]
			fit2$genes                     = probesAnn[rownames(fit2$genes),]
			fits[[model]]$fits2[[contset]] = fit2
		}
	}
	return(fits)
}


# This function takes an eset object and populates its pData and fData with 
# db annotations
#-----------------------------------------------------------------------------#
# We assume consistency between eset and Array annotation. i.e. the arrayAnnot
# is exactly the same as eset (this one has hidden arrays).
# In other words the annotation sets everything
#-----------------------------------------------------------------------------#
# x : The ExpressionSet (ie. eset)
# data.dir       : Directory where is stored the eset_raw.Rdata
# strictHybs     : cf. annotate.eset *TODO describe strictHybs *
# third.id       : Used if you have a third table in your worksheet that describes treatment for patients.
#x = eset;data.dir=dirData;third.id="DonorID"# TEMP
annotate.eset <- function(x,data.dir=dirData,third.id=NULL,strictHybs=TRUE, param=param)
{
	
	param = attr(x,"param") # params should be retained!
	
    # clears the eset of all previous annotations
	x = new("ExpressionSet",exprs=exprs(x))
	
    # Append array annotations
	print("Populating pData")
	arraysAnn =  getArraysAnn(data.dir)  
	if( any( ! sampleNames(x) %in% rownames(arraysAnn)       ) ){
		print("Warning: the following arrays in the eset are not present in the  visible array annotation. Make sure they were hidden since they will be discarded from the resulting eset")
		print(sampleNames(x)[! sampleNames(x) %in% rownames(arraysAnn)  ])
		x = x[,sampleNames(x) %in% rownames(arraysAnn)]
	}
	if( !setequal( rownames(arraysAnn),sampleNames(x) )  & strictHybs ) # Strict hybs check
	{
		print(setdiff(rownames(arraysAnn),sampleNames(x)))
		stop("Warning: The above non-hidden arrays are in your annotation but not in the ExpressionSet object.")
	}
	pData(x) = arraysAnn[sampleNames(x),]
	
	
    # Append sample annotation
	samplesAnn =  getSamplesAnn(data.dir)
	if( any( ! pData(x)$SampleID %in%  rownames(samplesAnn)      ) ){stop("Some SampleIDs in the array annot absent from sample annotation")}
	samplesAnn = samplesAnn[pData(x)$SampleID,] # preparer merger
	samplesAnn$SampleID = NULL
	rownames(samplesAnn) = NULL
	pData(x) = data.frame(samplesAnn,pData(x),row.names="ArrayID",stringsAsFactors=FALSE)	
	pData(x)$ArrayID = sampleNames(x)
	print("done.")
	
    # Populate the fData
	print("Populating fData")
	probesAnn = getProbesAnn(data.dir) #,controls=TRUE
	common   = intersect(featureNames(x),rownames(probesAnn))
	if( any(!  featureNames(x) %in% rownames(probesAnn)       ) ){warning(paste("Some probes in the eset are absent from the visible probe annotation, they are probably control probes that were hidden. They will be dropped from the eset! ", length(common)," features left out of ",nrow(x),".",sep=""))}
	x        = x[common,]
	fData(x) = probesAnn[common,]
	print("done.")
	
    # Append third annotation if necessary
	if(!is.null(third.id))	
	{	
		thirdAnn = getThirdAnn(data.dir,third.id)
		if( any( ! pData(x)[[third.id]] %in%  rownames(thirdAnn)      ) ){stop("Some ids in the samples annot absent from 3rd annotation")}
		thirdAnn = thirdAnn[pData(x)[[third.id]],] # preparer merger
		thirdAnn[[third.id]] = NULL
		rownames(thirdAnn) = NULL   
		pData(x) = data.frame(thirdAnn,pData(x),row.names="ArrayID",stringsAsFactors=FALSE)
	}	
	pData(x)$ArrayID = sampleNames(x)
	
    # Rehabilitate params
	attr(x,"param")  = param
	
	return(x)
}

###################
## Quality Control
###################

# Remove Outliers from ExpressionSet.
#-----------------------------------------------------------------------------#
# You need to have a column named "outlierFlag" within your sample annotation
# sheet. In order to remove a sample, write a short sentence that justifies
# the removal of this sample. Any row with a value stored in its outlierFlag
# will be systematically removed.
#-----------------------------------------------------------------------------#
qcRemoveOutliers <- function(x,data.dir=dirData,verbose=FALSE) 
{
	pa = getArraysAnn(data.dir=data.dir)
	outliers = pa$ArrayID[pa$outlierFlag!=""]
	x = x[,!sampleNames(x) %in% outliers]
	if(verbose){print(pa[outliers,])}
  	return(x)
}

# qcHeatmap : generate plot for QC Heatmap
# eset   : an expressionSet
# output : a levelplot to draw 
qcHeatmap <- function(eset) {

	dd = dist2(log2(exprs(eset)))
	diag(dd) = 0
	dd.row = as.dendrogram(hclust(as.dist(dd)))
	row.ord = order.dendrogram(dd.row)
	legend = 	list(top=list(fun=dendrogramGrob,args=list(x=dd.row,side="top")))
	tr = levelplot(dd[row.ord,row.ord],scales=list(x=list(rot=90)),
				   xlab="",ylab="",legend=legend)
	return(tr)
}


# qcPostPP : generates various diagnostic plots for an object
#            of class ExpressionSet.
# x    : an ExpressionSet. Its phenoData slot must have rownames in the 
#        chipID_strip format and contain the column "Class" 
# log2 : Logical. Defaults to TRUE. Whether or not exprs(x) is log scaled.
# maplots  Logical: whether maplots should be drawn (because this is usually very slow)
# destdir : The directory to which the resulting plots should be written.
#           Defaults to "diagnostic_plot/postprocessing/"
#-----------------------------------------------------------------------------
# TODO : prcomp analysis followed by LMM model fitting to report on the effect
#        of other covariates. 
qcPostPP <- function(x,log2=TRUE,dest.dir=file.path(dirData,"diagnostic_plot"),maplots=TRUE)
{
	#create directories
	dest.dir = file.path(dest.dir,"postprocessing")
	if(!file.exists(dest.dir)) {
		dir.create(dest.dir,recursive=TRUE)
	}
  	if(log2) {
		exprs(x) = 2**exprs(x)
  	}
	
	
	# Principal component analysis
  	pca  = prcomp(log2(exprs(x)),center=TRUE,scale=TRUE)
  	for(covr in colnames(pData(x))) 
	{
		co = pData(x)[[covr]]
		if( sum(!is.na(co)) > 0 ) # to avoid blank columns
        {
			if(class(co) %in% c("numeric","integer")){ 
				mycols = heat.colors(100)[  round(  100*(co-min(co)) / (max(co)-min(co)) )] # a crude color code for intensity
			}else{ # categorial covariates
				f = as.factor(co)
				mycols = rainbow(nlevels(f))
				mycols = sapply(f,function(l){mycols[which(l==levels(f))]})
			}
			
			
			png(file=file.path(dest.dir,paste("pca_",covr,".png",sep="")),width = 2000, height = 2000)
			pairs(pca$rotation[,1:5],pch=19,col=mycols,labels=paste(colnames(pca$rotation[,1:5]),"\n",100*round(pca$sdev**2/sum(pca$sdev**2),3),'%',sep=""))
			par(xpd=TRUE)
			
			if(class(co) %in% c("numeric","integer")){ 
				legend(0,1, c(min(co),max(co)) ,fill= c(heat.colors(100)[1],heat.colors(100)[100]) )
			}else{ # categorial covariates
				legend(0,1,levels(f),fill=rainbow(nlevels(f))) 
			}
			
			dev.off()
        }
	}
	
  	if(maplots)
	{
		## Within-chip MAXY plots. 
		pData(x)$chip.id  = as.factor(sapply( strsplit(rownames(pData(x)),split="_"),function(y) y[1]   )) #first retrieve the chip numbers and strip from the ArrayID.
		pData(x)$strip    = as.factor(sapply( strsplit(rownames(pData(x)),split="_"),function(y) y[2]   ))
		for(chip.id in unique(pData(x)$chip.id))
        {
			x.chip =  x[,which(pData(x)$chip.id==chip.id)] # this splits the original eset.
			png(file=file.path(dest.dir,paste("maXY_",chip.id,".png",sep="")),width = 800, height = 800)
			plotMAXY(exprs(x.chip), arrays = 1:ncol(x.chip),pch = 16, labels=pData(x.chip)$strip, main=chip.id)
			dev.off()
        }
		
		## Within-class MAXY plots : Allows for the appreciation of the difference between biological replicates, especially across different batches 
		for(clss in unique(pData(x)$Class))
        {
			x.clss =  x[,which(pData(x)$Class==clss)] # this splits the original eset.
			batch.labels = paste("\nBatch:", pData(x.clss)$batch )
			if(length( unique(pData(x.clss)$batch)) < 2 ){batch.labels=""} # batch id not shown if comes only from a single batch
			
			if(ncol(x.clss)>1) # because some classes may have only one chip
            {   
				png(file=file.path(dest.dir,paste("within_class_maXY_",clss,".png",sep="")),width = 800, height = 800)
				plotMAXY(exprs(x.clss), arrays = 1:ncol(x.clss),pch = 16, labels=  paste( rownames(pData(x.clss)), batch.labels  ,sep=" ") ,main=clss)
				dev.off()
            }
        }
	}
	
	if(ncol(x) > 100) {
		width = 12 * ncol(x)					   
	} else {
		width = 1100
	}	
	
	
	# Boxplots
  	png(file=file.path(dest.dir,paste("boxplot.png",sep="")),width = width, height = 1100)
	boxplot(log2(exprs(x))~col(exprs(x)), names = rownames(pData(x)),
			las=2,cex.axis = 0.8,pch=".")
  	dev.off()
	
	# density plots
  	png(file=file.path(dest.dir,paste("densityplot.png",sep="")),width = 800, height = 600)
  	plotDensitiesMat(exprs(x)) # log done inside this function by default
  	dev.off()
	
	
	# QCheatmap
  	png(file=file.path(dest.dir,paste("QC_Heatmap.png",sep="")),
		width = width + 100, height = width + 100)
  	print(qcHeatmap(x)) # log2 performed inside qcHeatmap
  	dev.off()
  	
}



###################
## Preprocessing
###################


# ppComBat : A wrapper for the ComBat function that implements
#            the ComBat batch effect adjustment method described in
#            Biostatistics (2007), 8, 1, pp. 118-127. It returns
#            the object x with its exprs() slot substituted with
#            the adjusted values.
# x        : an object of class ExpressionSet. Its pData has to have columns
#             "Class" and "batch".
# batch    : the name of the column of pData(x) to use as a batch. Defaults to "batch"
#             for backkward compatibility.
#
# covariates : a character vector listing the covariates to include in the linear model.
#              Defaults to "Class" for backkward compatibility. covariates=c() for no covariates
# output   : an ExpressionSet with batch-adjusted expression values
ppComBat <- function(x=eset,batch="batch",covariates=c("Class"),data.dir=dirData)
{
	
	#create directory
	if(!file.exists( file.path(data.dir,"diagnostic_plot")  )) {
		dir.create(file.path(data.dir,"diagnostic_plot"))
	}
	
	# Checks
	if(length(unique(pData(x)[[batch]]))<2){print("Just one batch. Nothing to correct here!");return(x)}
	# TODO : other checks, e.g. check if model can be fitted at all
	
	print("Adjusting batch effects using ComBat.")
	write.table(
				rbind( c("gene",colnames(exprs(x)))   ,cbind( rownames(exprs(x)) ,exprs(x) ) )
				, quote     = FALSE
				, row.names = FALSE
				, col.names = FALSE
				, file="ComBatExprs.txt"
				, sep="\t") # write the expressions to input in Combat
	
	
	df = data.frame(
					"Array name"  = colnames(exprs(x))
					,"Sample name"= colnames(exprs(x))
					,"Batch"      = pData(x)[[batch]] )
	
	if(length(covariates)>0)
	{
		for (i in 1:length(covariates))
		{
			df[[paste('Covariate',i)]] = pData(x)[[ covariates[i] ]]
		}
		#df[["Covariate 1"]] = apply(pData(x)[,covariates,drop=FALSE],MARGIN=1,FUN=function(ro)paste(ro,collapse="_"))
	}
	
	write.table(
				df
				, quote = FALSE
				, row.names = FALSE
				, col.names = TRUE
				, file="ComBatSamplesAnnots.txt"
				, sep="\t")# Write the sample description file for combat
	
	
	png(file=file.path(data.dir,"diagnostic_plot","comBat_priors_qqplots.png"))
	combat.exprs = ComBat("ComBatExprs.txt","ComBatSamplesAnnots.txt"
						  , type='txt', write=F, covariates='all', par.prior=T, filter=F, skip=1, prior.plots=T)
	dev.off()
	
	exprs(x)[1:nrow(x),1:ncol(x)] = as.matrix(combat.exprs[,-1])# replace the values in exprs(x)
	
	file.remove(c("ComBatSamplesAnnots.txt","ComBatExprs.txt")) # clean up temp files
	return(x)
}


# preprocessing : Filter genes
# eset   : A (normalized) expressionSet
# iqr    : the interquartile range
# bg     : background treshold for intensity filter.
# bgCount: number of samples required to be over background
# report : Logical. Deprecated in version 2, but left there for backwards compatibility.
# output : a filtered expressionSet
ppFilter <- function(eset, iqr=param$iqr, bg=param$bg*1.3, bgCount=param$bgCount, 
	report=F) {
	
	x = exprs(eset)
	
	f2       <- function(y)(IQR(y)>iqr)
	f1       <- kOverA(bgCount,log2(bg))
	
	ff       <- filterfun(f1,f2)
	selected <- genefilter(x,ff)
	print(paste(length(which(selected)),"genes left"))
	
	eset[selected,]
}

# preprocessing : Normalize - returns log2 and (background)surrogated data 
# x      : expression data from an expressionSet (eg. exprs(eset))
# bg     : background value for surrogate replacement policy
# method : normalization method (eg. "quantile")
# output : a normalized expressionSet 
ppNormalize <- function(x, bg=param$bg, method=param$normalizationMethod) {
	# Normalize
	dataNormalized <- normalizeBetweenArrays(x, method=tolower(method))
	#bg <- 2^ceiling(log2(bg))
	datNormFilteredMin <- replace(dataNormalized,which(dataNormalized<bg),bg)
	datNormFilteredMin <- replace(dataNormalized,which(is.na(dataNormalized)),bg)
	log2(datNormFilteredMin)
	
}


# qcDiagnose : generate png diagnostic plots for exprs(exprSet)
# x : an expressionSet
# title : nametag to append to each plot (eg. "raw" or "norm")
# width, height : dimensions for the boxplot.
# scaling: Logical. Whether you want the dimensions scaled to accomodate large eset.
# report : Logical. Deprecated in version 2, but left there for backwards compatibility.
qcDiagnose <- function(eset, title="raw", width=1100, height=850, data.dir=dirData,
						scaling=T, report=F) 
{
	#create directory
	if(!file.exists( file.path(data.dir,"diagnostic_plot")  )) {
		dir.create(file.path(data.dir,"diagnostic_plot"))
	}
	
	
	##########################################
	## Density plot
	##########################################
	png(file= file.path(data.dir,"diagnostic_plot",paste("densityplot_",title,".png",sep="")) ,
		width = 800, height = 600)
	plotDensitiesMat(exprs(eset))
	dev.off()
	
	if(scaling) {
		if(ncol(eset) > 100) {
			width = 12 * ncol(eset)					   
		}
	}
	##########################################
	## Boxplots
	##########################################
	png(file= file.path(data.dir,"diagnostic_plot",paste("boxplot_",title,"1.png",sep=""))   ,
		width = width, height = height)
	if(title=="raw") { 
		boxplot(log2(exprs(eset))~col(exprs(eset)), names = colnames(exprs(eset)),
				las=2,cex.axis = 0.8,pch=".", 
				main="Gene expression intensities for all arrays")}
	else { boxplot(exprs(eset)~col(exprs(eset)), 
				   names = colnames(exprs(eset)),las=2,cex.axis = 0.8,pch=".", 
				   main="Gene expression intensities for all arrays")}
	dev.off()
	
	png(file=file.path(data.dir,"diagnostic_plot",paste("boxplot_",title,"2.png",sep="")),
		width = width, height = height)
	if(title=="raw") {
		boxplot(log2(exprs(eset))~col(exprs(eset)), names = pData(eset)$Class,
				las=2,cex.axis = 0.8,pch=".")}
	else {boxplot(exprs(eset)~col(exprs(eset)), names = pData(eset)$Class, 
				  las=2,cex.axis = 0.8,pch=".",
				  main="Gene expression intensities for all arrays")}
	dev.off()
	
	## Heatmap QC plot
	pg = qcHeatmap(eset)
	png(file=file.path(data.dir,"diagnostic_plot","QC_Heatmap.png") ,
		width = width+100, height = width+100)
	print(pg) 
	dev.off()
	
	
}



## TODO : implement.
pp <- function(x=eset.raw)
{
    # This could be a wrapper for one-shot preprocessing given the raw eset (with parameters as attribute of course)
}

ppAverageTechnicalReplicates <-function(x,by="SampleID")
{
	replicated = names(which(table(pData(x)[[by]])>1))
	for(sid in replicated)
	{
		aids =  sampleNames(x)[pData(x)[[by]]==sid]
		exprs(x)[,sampleNames(x)==aids[1] ] = apply(exprs(x[,sampleNames(x) %in% aids]),MARGIN=1,mean)
		x = x[, ! sampleNames(x) %in% aids[2:length(aids)] ] # remove others from the eset
	}
	return(x)
}


#########################
## Unsupervised Analysis
#########################


# unaCluster : unsupervised analysis - samples clustering
#              Generates 4 clusterings : - 2 with Euclidian Distance
#                                        - 2 with correlation distance
#          For each distance metrics, generates two different labels,
#          one with ArrayID, one other with "labels" parameter.
# x      : an expressionSet
# labels : what should appear as the phenotypes. By default it will show 
#          pData(x)$Class. 
# report : logical. Deprecated.
# ret    : Set to TRUE if you want to return a list of clustering objects
# width, 
# height : dimensions for the .png files. Sometimes for large data sets
#          it is needed to have a larger diagram.
# cex    : labels font size.
unaCluster <- function(x,labels=NULL,report=F, width=1600, height=800, 
cex=1, ret=FALSE, data.dir=dirData)
{
	if(!file.exists(   file.path(data.dir,"diagnostic_plot") )) {
		dir.create(file.path(data.dir,"diagnostic_plot"))
	}
	#setwd("diagnostic_plot")
	ddCor <- cor.dist(t(exprs(x)), abs = FALSE, diag = FALSE, upper = FALSE)
	ddEuc <- euc(t(exprs(x)), diag = FALSE, upper = FALSE)
	clustCor <- hclust(ddCor, method="average")
	clustEuc <- hclust(ddEuc, method="average")
	
	clustFile = c("EucClust.png","EucClustpheno.png",
				  "CorClust.png", "CorClustpheno.png")
	
	#set default labels if none is specified
	if(is.null(labels)) {
		labels = as.character(pData(x)$Class) 
	}
	
	for(i in 1:length(clustFile)) {
		png(file= file.path(data.dir,"diagnostic_plot",clustFile[i])       ,
			width = width, height = height)
		par(cex=0.8)
		
		if(i == 1) {
			plot(clustEuc,label=colnames(exprs(x)), 
				 main="hierarchical Euclidan clustering dendogram",
				 sub=paste("average linkage, correlation, number of genes : ",
						   nrow(exprs(x)),sep=""), cex=cex)	
		}
		if(i == 2) {
			plot(clustEuc,label=labels,
				 main="hierarchical clustering dendogram",
				 sub=paste("average linkage, correlation, number of genes : ",
						   nrow(exprs(x)),sep=""), cex=cex)
		}
		if( i == 3) {
			plot(clustCor,label=colnames(exprs(x)),
				 main="hierarchical clustering dendogram",
				 sub=paste("average linkage, correlation, number of genes : ", 
						   nrow(exprs(x)),sep=""), cex=cex)
			
		}
		if(i == 4) {
			plot(clustCor,label=labels,
				 main="hierarchical clustering dendogram",
				 sub=paste("average linkage, correlation, number of genes : ",
						   nrow(exprs(x)),sep=""), cex=cex)
		}
		
		dev.off()
	}
	
	if(ret){return(list("euc"=clustEuc,"cor"=clustCor))}
}

# Description : MDS plot from an eset 
#
# eset : Object of class ExpressionSet
# fn : file name. Will print a multi-page pdf
# color.by.column,shape.by.column : column names of pData(eset) to be used for coloring and shaping the dots
# alt.color.column.name,alt.shape.column.name : Prettier names for the previous.
unaMDS <- function(eset,fn=file.path(dirData,"diagnostic_plot","mds.pdf")
,color.by.column,shape.by.column,alt.color.column.name=color.by.column, alt.shape.column.name=shape.by.column,ret=FALSE)
{
    mds  = cbind( as.data.frame(cmdscale(euc(eset),k=3)),  pData(eset)      )
    colnames(mds)[1:3] = gsub('V','Dimension',colnames(mds)[1:3],perl=TRUE)
    colnames(mds)[grepl(color.by.column,colnames(mds))] = alt.color.column.name
    colnames(mds)[grepl(shape.by.column,colnames(mds))] = alt.shape.column.name
#print(head(mds))
    pdf(fn,width=7,height=7)
    
	( p1 <- qplot(x=Dimension1,y=Dimension2,data=mds ,main="front",color=mds[[alt.color.column.name]], shape=mds[[alt.shape.column.name]])+ scale_colour_hue(alt.color.column.name)+scale_shape(alt.shape.column.name)+opts(aspect.ratio=1) ) 
	( p2 <- qplot(x=Dimension1,y=Dimension3,data=mds ,main="top",color=mds[[alt.color.column.name]], shape=mds[[alt.shape.column.name]])+ scale_colour_hue(alt.color.column.name)+scale_shape(alt.shape.column.name)+opts(aspect.ratio=1) ) 
	( p3 <- qplot(x=Dimension2,y=Dimension3,data=mds ,main="side",color=mds[[alt.color.column.name]], shape=mds[[alt.shape.column.name]])+ scale_colour_hue(alt.color.column.name)+scale_shape(alt.shape.column.name)+opts(aspect.ratio=1) ) 
    print(p1)
    print(p2)
    print(p3)
    dev.off()  
    if(ret){return(list("p1"=p1,"p2"=p2,"p3"=p3)) }
}


###################
## DGEA
###################


# Categorical Linear Fits for DGEA
#-----------------------------------------------------------------------------#
# The following function implements two of the most common designs : 
#  - interaction terms of a two way ANOVA 
#  - one-way ANOVA comparison to a baseline.
# Nothing stops anyone from creating their own fit,fit2 objects. 
# Current reporting function expect fit$variables thoughl
#-----------------------------------------------------------------------------#
# x            : An ExpressionSet. (eg. eset)
# variables    : A character list of Factors from pData(x) used to fit the 
#                model. (eg. c("Tissue","Vaccine"))
# dupcor.block : A character string corresponding to the Factor in pData(x) 
#                used for duplicate correlation mixed modeling. (eg. "Leuka")
# fixed.effect.block : A character string corresponding to the Factor in pData(x) 
#                used fixed-effect blocking. (eg. "Leuka" or "batch")
# output       : a fit object.
dgLimmaCategoricalFit <- function(x,variables,dupcor.block=NULL,fixed.effect.block=NULL)
{
	# Checks
	if(any(!variables %in% varLabels(x))){stop("Some element in variables is not present in the pData")}
	if( (!is.null(dupcor.block)) & (!is.null(fixed.effect.block)) ){warning("Both random and fixed effect blocking. Are you sure of this?")}
	
	# Extract the parameter object from the eset
	param = attr(x,"param") 

	# make sure they are ordered factors
	Class = factor(apply(pData(x)[,variables,drop=FALSE],MARGIN=1,function(ro)paste(ro,collapse="___")))
	if(!is.null(fixed.effect.block)) { # Include blocking variable in the design
		block  = as.factor(pData(x)[[fixed.effect.block]])
		design  <- model.matrix(~ -1+Class+block)
		colnames(design) = gsub('(^Class|^block)','',colnames(design),perl=TRUE) 
	}else{ # No blocking
		design = model.matrix(~-1+Class)	
		colnames(design) = gsub('^Class','',colnames(design),perl=TRUE)	
	}
	rownames(design) = sampleNames(x)
	
	# Run Limmas duplicate correlation for models with correlated residuals within a block
	if(!is.null(dupcor.block))
	{
		print("Running duplicate Correlation....")
		dupcor=duplicateCorrelation(exprs(x),design=design,block=as.factor(         pData(x)[[dupcor.block]]        ))
		fit  <- lmFit(x,block=as.factor(    pData(x)[[dupcor.block]]      ),design,correlation=dupcor$consensus)
		fit$dupcor  = dupcor # Non-standard slot !
		fit$dupcor.block  = dupcor.block # Non-standard slot !
	}else{
		fit  <- lmFit(x,design)
	}  

	# Non-standard slots
	fit$variables = pData(x)[,variables,drop=FALSE]
	fit$TimeStamp = Sys.time()
	fit$author = param$author
	#fit$Class = Class
	
	return(fit)
}




# Fit a group of contrasts to a linear model.
#-----------------------------------------------------------------------------#
# fit         : a fit object, namely the output from dgLimmaCategoricalFit.
# vecCont     : A list of contrasts taken from the output of dgANOVAContrasts.
# gsea.method : A character string refering to the method used for patways
#               analysis. Either "none", "roast" or "romer". 
# gsea.eset   : An ExpressionSet. Default to NULL. Use only if you set a
#               GSEA method.
# output      : A fit2 object, i.e. estimated coefficients for a given set of 
#               contrasts (MArrayLM object).
dgLimmaContrastsFit <- function(fit,vecCont,data.dir=dirData,gsea.method="none",gsea.eset=NULL)
{
	# Check arguments
	if(! gsea.method %in% c("none","roast","romer","AvgLogP")){stop("Invalid GSEA method")} # TODO: romer, wilcoxGST

	# Single-gene analysis	
	cont.matrix  = makeContrasts(contrasts=vecCont, levels=fit$design)
	fit2         = contrasts.fit(fit,cont.matrix)
	fit2         = eBayes(fit2)

	# Non-standard slots in the fit objects
	fit2$TimeStamp = Sys.time()
	fit2$cont.matrix = cont.matrix # who knows we might need this later on
	
	# Gene set analysis
	if(gsea.method %in%  c("roast","romer","AvgLogP"))
	{	
		if(is.null(gsea.eset) ){stop("Unfiltered eset must be provided for roast")}	#check eset
		if(any(fit2$genes$ProbeID!=featureNames(eset))){stop("fit object not consistent with eset")}	#check eset
		gene.sets  = getGeneSets(data.dir,ret.r=TRUE) # Load the gene sets from database
		iset = sapply(gene.sets,function(gs){ fit2$genes$ProbeID %in% gs},simplify=FALSE)
		
		fit2$gsea = sapply(colnames(fit2),function(cont)
		{
			print(paste("Beginning one contrast. ",Sys.time()))
			if(gsea.method=="roast"){
				return( mroast(iset, exprs(gsea.eset), fit2$design, contrast=cont.matrix[,cont], adjust.method="BH") )
			}
			if(gsea.method=="romer"){
				return( romer( sapply(iset,which), exprs(gsea.eset), fit2$design, contrast=cont.matrix[,cont]) )
				#romer(iset,y,design,contrast=ncol(design),array.weights=NULL,block=NULL,correlation,set.statistic="mean",nrot=9999)
			}
			if(gsea.method=="AvgLogP"){
				return( AvgLogP.gsea( isets=iset, lm=fit2,coef=cont ))
				#romer(iset,y,design,contrast=ncol(design),array.weights=NULL,block=NULL,correlation,set.statistic="mean",nrot=9999)
			}
		},simplify=FALSE)
		# gseaF : roast does not support F-test yet? provide wilcoxonGST
	#	F.stat = fit2$F; names(F.stat) = rownames(fit2)
	#	missing.probes = featureNames(x)[! featureNames(x) %in% rownames(fit2)]
	#	missingF = rep(0,times=length(missing.probes));names(missingF) = missing.probes
	#	F.stat = c(F.stat,missingF)
	#	F.stat = F.stat[featureNames(x)]
	#	fit2$gseaF = sapply(iset,function(selected){geneSetTest(selected, F.stat, type="f")}) # TODO: does roast make sense with F stat
	}
	attr(fit2,"gsea.method") = gsea.method
	return(fit2)
} 

# a naivve gsea method: the average -logP
AvgLogP.gsea <- function(isets,lm,coef=1)
{
	ts = abs(lm[,coef]$t[,1])
	scores = data.frame("Mixed"=sapply(isets,function(iset){
	 mean(ts[iset])
	}) )

	return(list("P.Value"=scores ))
	
}


# Create contrast formulas on the fly with dgANOVAContrasts
#-----------------------------------------------------------------------------#
# This is a helper function that automatically creates contrasts of for 
# a two-way anova. 
# eset       : An ExpressionSet. (eg. eset)
# variables  : A character list of Factors from pData(x) used to fit the 
#                model. (eg. c("Protein","Adjuvant"))
# f1         : A character string, representing the first Factor from pData(x) 
#              used to fit the model.
# f1.baseline: A character string representing the baseline to substract from 
#              f1.
# f2         : A character string, representing the second Factor from pData(x) 
#              used to fit the model.
# f2.baseline: A character string representing the baseline to substract from 
#              f2.
# output     : A List of contrasts necessary to input to dgLimmaContrastsFit.
# sep        : factor separator, e.g. as in "Mu__Treated-Mu___Untreated"
dgANOVAContrasts <- function(x,variables
				,f1=NULL,f1.baseline=NULL 
				,f2=NULL,f2.baseline=NULL ,sep="___")
{

	# Checks
	if(any(!variables %in% varLabels(x))){stop("Some element in variables is not present in the pData")}
	if(  ! f1 %in% variables) {stop("First factor not in variables")} 
	if( ! f1.baseline %in% pData(x)[[f1]] ){stop("Invalid baseline level for the first factor")}
	if( ! f1.baseline %in% pData(x)[[f1]] ){stop("Invalid baseline level for the first factor")}
	if( !is.null(f2)){ 
		if(  ! f2 %in% variables) {stop("Second factor not in variables")} 
		if(  ! f2.baseline %in% pData(x)[[f2]] ){stop("Invalid baseline level for the second factor")}
	}

	l = list()

	# First on f1
	# Isolate each by combinations of fixed
	cts = unique(pData(x)[,variables,drop=FALSE])
	rownames(cts) = NULL
	fixed = apply(cts[,setdiff(variables,f1),drop=FALSE],MARGIN=1,function(ro)paste(ro,collapse=sep))
	cts = by(cts,INDICES=fixed,function(df){
	if( ! f1.baseline %in% df[[f1]] | nrow(df)<2){return(NULL)} # Nothing to compare to the baseline
	cls  =  apply(df[,variables,drop=FALSE],MARGIN=1,function(ro)paste(ro,collapse=sep))
		bl.cls = cls[df[[f1]]==f1.baseline]
		return(paste(setdiff(cls,bl.cls),bl.cls,sep="-"))
	})
	cts = as.vector(unlist(cts));names(cts)=NULL
	if(is.null(f2)){return(cts)}
	l[[f1]] = cts

	# ELSE, go on to f2
	cts = unique(pData(x)[,variables,drop=FALSE])
	rownames(cts) = NULL
	fixed = apply(cts[,setdiff(variables,f2),drop=FALSE],MARGIN=1,function(ro)paste(ro,collapse=sep))
	cts = by(cts,INDICES=fixed,function(df){
	if(! f2.baseline %in% df[[f2]] | nrow(df)<2){return(NULL)} # Nothing to compare in this case
	cls  =  apply(df[,variables,drop=FALSE],MARGIN=1,function(ro)paste(ro,collapse=sep))
		bl.cls = cls[df[[f2]]==f2.baseline]
		return(paste(setdiff(cls,bl.cls),bl.cls,sep="-"))
	})
	cts = as.vector(unlist(cts));names(cts)=NULL
	l[[f2]] = cts

	# Now t
	# Now the interaction
	cts = unique(pData(x)[,variables,drop=FALSE])
	rownames(cts) = NULL
	fixed = apply(cts[,setdiff(variables,c(f1,f2)),drop=FALSE],MARGIN=1,function(ro)paste(ro,collapse=sep))
	cts = by(cts,INDICES=fixed,function(df){

		if( length(setdiff(df[[f1]],f1.baseline))==0|length(setdiff(df[[f2]],f2.baseline))==0 ){return(NULL)} # Nothing to compare in this case
		cls  =  apply(df[,variables,drop=FALSE],MARGIN=1,function(ro)paste(ro,collapse=sep))
		sl = c()
		for(f1l in setdiff(df[[f1]],f1.baseline) ){		
			for(f2l in setdiff(df[[f2]],f2.baseline) ){
			if(
			any(df[[f1]]==f1l &  df[[f2]]==f2l)
			&any(df[[f1]]==f1.baseline &  df[[f2]]==f2l)
			&any(df[[f1]]==f1l &  df[[f2]]==f2.baseline)
			&any(df[[f1]]==f1.baseline &  df[[f2]]==f2.baseline )){
					sl	 = c(sl,
					paste('('
						,cls[  df[[f1]]==f1l &  df[[f2]]==f2l ],'-'
						,cls[df[[f1]]==f1.baseline &  df[[f2]]==f2l ],')-('
						,cls[df[[f1]]==f1l &  df[[f2]]==f2.baseline ],'-'
						,cls[df[[f1]]==f1.baseline &  df[[f2]]==f2.baseline ]
						,')',sep="")
					)	
				}
			}
		}
		return(sl)
		})
	cts = as.vector(unlist(cts));names(cts)=NULL
	l[[paste(f1,f2,sep="_int_")]] = cts
	return(l)	
}


##
## dgCompareContrasts: Returns a scatterplot (ggplot2 object) of the fold-changes
## of two contrasts. Text labels are added to the plot according to some p-value
## or fc thresholds.
##
## main : Title of the plot
##
## fitX,fitY: A fit object that holds the contrast for the X,Y axis of the scatterplot.
## These object may well hold other contrasts. Shoud the compared contrasts be in the
## same fit object, simply subset it. Fit objects should have slot fit$genes populated
## a data.frame which has at least column "SYMBOL". Symbols ultimately used for on
## the plot are those from fitX$genes
##
## ctX,ctY : strings, the names of the contrasts as in colnames(fit)
##
## by.x, by.y : name of the column of fit$genes to use for joining. Duplicates are appendend ProbeID instead og beiing collapsed.
##           This implies the need to collapse duplicates prior to joining. Defaults to Probe_Id (ILLXXXX), the nucl. sequence.
##  
## ctXY: Interaction contrats. Not used right now.
##
## adjust.method: p-value adjustment method. Defaults to BH. 
##
## p.value.X,p.value.Y,lfc.X,lfc.Y: p-values and logFC cutoffs. They are used to
##  determine for which gene labels (gene symbols) are added to the plot.
##
## xlab,ylab : Alternative names for contrasts for prettier axis names.
##
## xlim,ylim: x,y limits. Default to -max,+max of logFCs. Ideally use identical limits
##  to get a square plot.
##
## asp: The aspect ratio of the plot. Overrides xlim, ylim.
##
## text.size: numerical value to be passed to geom_text for the gene labels.
##
## gene.color.groups: A list of gene symbols vectors, eg. list("IFN"=c("STAT1","IRF7")
##    ,"RP"=c("RPS20","RSP19")).
##
## smooth : should a loess smooth be added to the plot. Defaults to TRUE
##
## do.venn : Should Vennerable objects be computed and appended to the output
##
## outfile: name for a png where to plot.
##
##
##  VALUE: ggplot2 object.
##
##
## TODO :
##       - Improve coloring method. It sucks.
##       - find a way not to exlcude NAs, like probeID, RefSEq... but then cross-platform problem
##       - add support for not collapsing duplicate symbols
##       - Catch possible errors : too stringeant cutoff...
##       - Maybe merging shoulf be done on the bassis of somethig esle than gene....
dgCompareContrasts <- function(
                               fitX,
                               fitY,
                               ctX = colnames(fitX)[1],
                               ctY = colnames(fitY)[1],
			       by.x        = "Probe_Id",
			       by.y        = by.x, 
                               ctXY = NULL,
                               adjust.method = "BH", 
                               p.value.X = 0.05,
                               p.value.Y = p.value.X,
                               lfc.X = 0,
                               lfc.Y = lfc.X,
                               main = "Comparison",
                               xlab = paste(ctX,' (log2 scale)',sep=""),
                               ylab = paste(ctY,' (log2 scale)',sep=""),
                              # xlim = c( -1*max(max(abs(fitX$coef[,ctX])),max(abs(fitY$coef[,ctY])))  , max(max(abs(fitX$coef[,ctX])),max(abs(fitY$coef[,ctY])))  ),
			       xlim = NULL,
                               ylim = xlim,
			       asp=1,
                               text.size = 1.5,
                               gene.color.groups = NULL,
			       gene.filter=NULL,
                               smooth      = TRUE,
                               do.venn     = TRUE,
                               outFile     = NULL
                                 )
{

   ## Checks
   if(! ctX %in% colnames(fitX))         {stop("ctX arugment is not a contrast in fitX.")}
   if(! ctY %in% colnames(fitY))         {stop("ctY argument is not a contrast in fitY.")}
   if(! by.x %in% colnames(fitX$genes) ) {stop("by.x arugment is not a column name if fitX$genes.")}
   if(! by.y %in% colnames(fitY$genes) ) {stop("by.y arugment is not a column name if fitY$genes.")}
 
   ## Prepare the data data.frame
   data =list( # This puts the info in fit in a more convenient, sortable format
     "X"=  data.frame(
       "logfc"=fitX$coefficients[,ctX]
       ,"fc"= sapply( fitX$coefficients[,ctX]  ,function(z) toFold(z))
       ,"p"=fitX$p.value[,ctX] # pavlue used in sorting required for collapsing
       ,"dt" = decideTests(fitX,adjust.method=adjust.method,lfc=lfc.X,p.value=p.value.X)[,ctX] # decidetest is later used for text label plotting decision
       ,"SYMBOL" = fitX$genes$SYMBOL
	,"ProbeID" = fitX$genes$ProbeID
       ), #,row.names="ProbeID"
     "Y"=  data.frame(
       "logfc"=fitY$coefficients[,ctY]
       ,"fc"= sapply( fitY$coefficients[,ctY]  ,function(z) toFold(z))
       ,"p"=fitY$p.value[,ctY]
       ,"dt" = decideTests(fitY,adjust.method=adjust.method,lfc=lfc.Y,p.value=p.value.Y)[,ctY]
       ,"SYMBOL" = fitY$genes$SYMBOL
	,"ProbeID" = fitY$genes$ProbeID
       )#,row.names="ProbeID"
     )
     # Define the joining columns
     data$X[[by.x]] =  fitX$genes[[by.x]]
     data$Y[[by.y]] =  fitY$genes[[by.y]]

     ## Make sure SYMBOL column is complete, order by p-value and collapse if necessary, print warning
     data = mapply(function(df,by)
      {
       	df = df[order(df$p,decreasing=FALSE),]	# Sort by p.value
	if(any(is.na(df$SYMBOL) | df$SYMBOL=="")) # Check for absent SYMBOLS
	{
		warning("Some elements of the SYMBOL column in fitX/Y are NA or empty. They are replaced by the corresponing ProbeID")
		df$SYMBOL[is.na(df$SYMBOL) | df$SYMBOL==""] =  df$ProbeID[is.na(df$SYMBOL) | df$SYMBOL==""] # just to make sure it is complete
	}
	if( anyDuplicated( df[[by]]) ) # Collpase if necessary
	{
        	df = df[ !duplicated( df[[by]]) ,]
		warning("There were duplicates in fitX/Y[[by.X/Y]] in terms of by.X/Y. Collapsed by the maximum p.value. Watch out for probe effect")
	}
        return(df)
      },data,list(by.x,by.y),SIMPLIFY=FALSE,USE.NAMES=TRUE )

      ## Join!
      print(paste("About to join the two datasets. "
		,length(intersect(data$X[[by.x]],data$Y[[by.y]]))," features are in common, "
		,length(setdiff(data$X[[by.x]],data$Y[[by.y]]))," are unique to X and "
		,length(setdiff(data$Y[[by.y]],data$X[[by.x]]))," are unique to Y.",sep="")  )
      colnames(data$X) = paste("X.",colnames(data$X),sep="")
      colnames(data$Y) = paste("Y.",colnames(data$Y),sep="")
      data = merge(data$X,data$Y,by.x=paste("X.",by.x,sep=""),by.y=paste("Y.",by.y,sep=""),all.x=FALSE,all.y=FALSE)
   
     ## Set SYMBOL as from fitX
     data$SYMBOL = data$X.SYMBOL
     
     
   ## Prepare color column
   if(!is.null(gene.color.groups))
     {
       print(paste(
		length(intersect(unlist(gene.color.groups,recursive=TRUE ),data$SYMBOL))
		," out of ",length(unique(unlist(gene.color.groups,recursive=TRUE )))
		," symbols from color groups map to SYMBOLs of the data.",sep=""))

       data$Group = "Others"
       overlapping = c()
       for(gr in names(gene.color.groups))
         {
           data$Group[ data$SYMBOL %in% gene.color.groups[[gr]] ] = gr
         }
	if(length(gene.color.groups)>=2)
	{
		overlapping = intersect(gene.color.groups[[1]], gene.color.groups[[2]])
		if(length(gene.color.groups)>2)
		{
			for(i in 3:length(gene.color.groups))
        	 	{
		   		overlapping = intersect(overlapping, gene.color.groups[[i]])
        	 	}
		}		
	}
	data$Group[ data$SYMBOL %in% overlapping ] = "Overlapping"
     }
   
   ## Now prepare the plot object
   data$Gene = data$SYMBOL ## These will be the text labels
	if(is.null(gene.filter)) # If text labels should be added to sig genes
	{
   		data$Gene[data$X.dt==0 & data$Y.dt==0 ]  = ""
	}else{
		data$Gene[! data$Gene %in% gene.filter]  = ""
	}# If text labels should to filter only


   p =  ggplot(data=data,aes(x=X.logfc,y=Y.logfc))+  geom_point(aes(),alpha=0.1)+
     geom_abline(aes(),alpha=0.5)+geom_text(aes(label= Gene),alpha=0.75,size=text.size)+
      opts(title = main)+xlab(xlab)+ylab(ylab)
   if(smooth){p=p+geom_smooth(aes(),linetype="dashed",alpha=0.5,color="blue")} ## Add or not the smooth
   if(!is.null(gene.color.groups)){p=p+aes(color=Group)+ scale_colour_brewer(pal = "Set1")}## Add color if applicable
  
 #  if (!is.na(asp)) 
    p <- p + opts(aspect.ratio = asp)
    if (!is.null(xlim)) {p <- p + xlim(xlim)}
    if (!is.null(ylim)) {p <- p + ylim(ylim)}


   ## Hijack ggplot2 object for venn diagram ?
   if(do.venn)
     {
       library(Vennerable)
       p$up   = list("X"=data$SYMBOL[data$X.dt ==  1] ,"Y"=data$SYMBOL[data$Y.dt ==  1]);names(p$up) =c(xlab,ylab)
       p$dwn  = list("X"=data$SYMBOL[data$X.dt == -1] ,"Y"=data$SYMBOL[data$Y.dt == -1]);names(p$dwn)=c(xlab,ylab)
       p$bth  = list("X"=data$SYMBOL[data$X.dt !=  0] ,"Y"=data$SYMBOL[data$Y.dt !=  0]);names(p$bth)=c(xlab,ylab)
       p$VennUp         = Venn(p$up)
       p$VennDown       = Venn(p$dwn)
       p$VennBoth       = Venn(p$bth)
     }


   ## Write to png
   if(!is.null(outFile)){
   png(file=outFile,width=960,height=960) 
   print(p)
   dev.off()
 }
   return(p)
}


####################
## General Functions
####################

# rename contrasts with aliases for report intelligibility purposes.
contrasts2aliases <- function(cts,data.dir=dirData,aliases=getContrastsAliases(data.dir))
{
	names(cts) = cts
	if(length(aliases)==0){return(cts)}
	cts[ cts%in% names(aliases)] = aliases[ cts[ cts %in% names(aliases)] ]
	return(  cts   )
}



# toFold : converts logFC into FC
toFold <- function(x) {
	x[x>=0] =   2^(x[x>=0])
	x[x<0]  =  -2^(-x[x<0])
	return(x)
}



multtest.BH <- function(pvals)
{
    library(multtest)
    adj = mt.rawp2adjp( pvals ,proc="BH")
    adj = adj$adjp[order(adj$index),]
    adj = adj[,"BH"]
	
    return(adj)
}


# TODO : this function needs to be revised before becoming a standard procedure. 
# lm=fit2;output.dir=file.path("comparisons");p.value=0.05;lfc=0.58;adjust.method="none"
my.Venns <-function(lm=fit2,output.dir=file.path("comparisons"),p.value=0.05,lfc=0.58,adjust.method="none")
{
	if(!file.exists(output.dir)){dir.create(output.dir)}
	tops = sapply(colnames(lm),function(cont){
		top = topTable(lm,coef=cont,p.value=p.value,lfc=lfc,number=nrow(lm),adjust.method=adjust.method)
		unique(top$SYMBOL)
		}
		,simplify=FALSE)
	comb = combinations(n=length(tops),r=3,v=names(tops))
	v = list()
	for(i in 1:nrow(comb))
	{
		png(file.path(output.dir,paste(paste(comb[i,],collapse="_VS_"),".png",sep="")),height=600,width=600)
		v[[i]] = Venn(tops[comb[i,]])
		plot(v[[i]])
		dev.off()
	}
	return(v)
}

# A general function to read "masked" tables. Either excel or tab delimited
read.masked.table <- function(fn,sheet=NULL)
{
	if(is.null(sheet)) # It is a tab-delimited file
	{	
		col.mask = strsplit(readLines(fn,n=1),split="\t")[[1]]
		x = read.delim(fn,skip=1,quote="",check.names=FALSE,comment.char="",fill=FALSE,stringsAsFactors=FALSE,colClasses="character")#,colClasses="character"
	}else{ # it is xls
		col.mask = unlist(read.xls(fn,sheet,check.names=FALSE,header=FALSE)[1,])
		x = read.xls(fn,sheet,skip=1,check.names=FALSE,comment.char="",fill=FALSE,stringsAsFactors=FALSE,colClasses="character")#
	}
 	x$hide[is.na(x$hide)] = ""
	x = x[!grepl('\\S',x$hide,perl=TRUE),!grepl('\\S',col.mask,perl=TRUE)] # columns with any non-whitespace character are masked
	x[,ncol(x)] = gsub(' *$','',x[,ncol(x)],perl=TRUE)# This removes the annoying trailing spaces left by read.xls in the last column...
	return(x)
}

# TODO: be able to parse illuminas format directly and hide the controls
getProbesAnn <- function(data.dir = dirData,empty.symbols="ProbeID") # ,controls=FALSE
{
	print("Reading gene annotation from file...")
	probes = read.masked.table(file.path(data.dir,"probes.txt"))
	
	if( any(  ! c("ProbeID","SYMBOL")%in% colnames(probes)) ){stop("Something wrong with probes.txt: ProbeID and/or SYMBOL column missing.")} 

	# make sure ProbeIDs are character, not integers
	probes$ProbeID  = as.character(probes$ProbeID) 
	probes$SYMBOL  = as.character(probes$SYMBOL)
	
       # Perform empty SYMBOL replacement
	probes$SYMBOL[probes$SYMBOL==""] = paste("ProbeID:",probes$ProbeID[probes$SYMBOL==""])
	
	# Assign rownames (primary key) to table
	if(anyDuplicated(probes$ProbeID)){
		stop("Values in ProbeID column of probes.txt do not form a proper primary key (uniqueness required)")
	}
	rownames(probes) = probes$ProbeID

	print("done.")
	return(probes)
}


getSamplesAnn <- function(data.dir = dirData) 
{
	print("Reading sample annotation from worksheet...")
	x = read.masked.table(file.path(data.dir,"worksheet.xls"),sheet=1)
	
	if( any(! c("SampleID")%in% colnames(x)) ){stop("Something wrong with sample annotation: SampleID column absent.")}
	x$SampleID  = as.character(x$SampleID) 

	# Assign rownames (primary key) to table
	if(anyDuplicated(x$SampleID)){
		stop("Values in SampleID column of first tab of worksheet.xls do not form a proper primary key (uniqueness required)")
	}
	rownames(x) = x$SampleID

	print("done.")
	return(x)
}


getArraysAnn <- function(data.dir = dirData) 
{
	print("Reading array annotation from worksheet...")
	x = read.masked.table(file.path(data.dir,"worksheet.xls"),sheet=2)
	
	if( any( ! c("ArrayID","SampleID","outlierFlag")%in% colnames(x)) ){stop("Something wrong with probe annotation in the second tab of worksheet.xls: Missing ArrayID or outlierFlag or SampleID")}
	x$ArrayID  = as.character(x$ArrayID)
	x$SampleID  = as.character(x$SampleID) 
	x$outlierFlag  = as.character(x$outlierFlag) 
	x$outlierFlag[is.na(x$outlierFlag)] = ""

	# Assign rownames (primary key) to table
	if(anyDuplicated(x$ArrayID)){
		stop("Values in ArrayID column of second tab of worksheet.xls do not form a proper primary key (uniqueness required)")
	}
	rownames(x) = x$ArrayID
	
	print("done.")
	return(x)
}

getThirdAnn <- function(data.dir = dirData,third.id="DonorID") 
{
	print("Reading third annotation from worksheet...")
	x = read.masked.table(file.path(data.dir,"worksheet.xls"),sheet=3)
	
	# Check presence of mandatory columns
	if( ! third.id %in% colnames(x) ){
		stop(paste("Something wrong with third annotation in the third tab of worksheet.xls: Column",third.id,"does not exist (Defaults to DonorID)."))}
	x[[third.id]]  = as.character(x[[third.id]])
	
       # Assign rownames (primary key) to table
	if(anyDuplicated(x[[third.id]])){
		stop(paste("Values in column",third.id,"of third tab of worksheet.xls do not form a proper primary key (uniqueness required)"))
	}
	rownames(x) = x[[third.id]]	

	print("done.")
	return(x)
}



# join.by : the column name of Probes.txt to match to the genes in the gene sets. This function will 
# this function returns the pathways either in its originial table format (for reporting),
# but also performs joining which will be used to score calculation later.
#  ProbeIDs could be pre-computed... but then it would complicate matters for the user... not sure
# TODO: return list of lists of probes of same gene SYMBOL... necessary for proper collapsing
getGeneSets <- function(data.dir=dirData,join.by = "SYMBOL",ret.r=FALSE,set.min.size=7)
{
	print("Reading gene sets table...")
	x = read.masked.table(file.path(data.dir,"geneSets.txt")  )#,sheet=1)

	# Check presence of mandatory columns
	if( ! "GeneSetID" %in% colnames(x) ){
		stop(paste("Something wrong with geneSets.txt: Column GeneSetID does not exist."))}
	x[["GeneSetID"]]  = as.character(x[["GeneSetID"]])
	
       # Assign rownames (primary key) to table
	if(anyDuplicated(x[["GeneSetID"]])){
		stop(paste("Values in column GeneSetID of geneSets.txt do not form a proper primary key (uniqueness required)"))
	}
	rownames(x) = x[["GeneSetID"]]	


	print("done.")
	if(!ret.r)
	{
		return(x)
	}else
	{
		obj = x$Members
		names(obj) = x$GeneSetID
		obj = sapply(obj,function(gs){unlist(strsplit(gs,split=","),recursive=TRUE)},simplify=FALSE) # converts vector to list of splitted
        # convert to probeIDS with sym
		probes = getProbesAnn(data.dir)
		probes.genes = probes[[join.by]] # the symbols correspondong to each probes
		names(probes.genes)  = probes$ProbeID  
		print("Mapping gene sets to probes... This might take a minute or two")
		print(Sys.time())
		obj = symbols2indices(obj,probes.genes) # from limma package to map gene set members to probes
		obj = sapply(obj,function(ind){names(probes.genes)[ind]})# CONVERT BACK TO PROBEids
        # Minimum set size
		print("Filtering the corresponding sets of probes by size in terms of mapped genes....")
		size.filter = sapply(obj,function(ps){  length(unique(probes.genes[ps]))>=set.min.size   })
		print(Sys.time())
		print("done")
		return(obj[size.filter])
	}		
}

getContrastsAliases <- function(data.dir=dirData)
{
    #print("Reading contrast aliases...")
	x = read.masked.table(file.path(data.dir,"contrasts_aliases.txt")  )#,sheet=1)
    #print("done.")
	if( any(! c("formula","alias")%in% colnames(x)) ){stop("Something wrong with contrasts aliases file. ")}
    #if(nrow(x)==0){return(x)}
 	formulas  = x$formula
	x = x$alias;names(x) = formulas
	return(x)
}


# getParams : loads default values for optional parameters
#              sets values for all parameters specified in "params.txt"
# file       : parameters file
# returns    : a list of parameters
#-----------------------------------------------------------------------------#
# Alternate  : You can declare a list named param and manually set the 
#              parameters without using a file.
getParams <- function(file="params.txt") {
	
	output <- list()
	
    #default values for optional parameters
	output$iqr = 0.3
	output$bg  = 40
	output$bgCount = 2
	output$targetsFile  = "sample_annotation.txt"
	output$contrastsFile = "contrasts.txt"
	output$normalizationMethod = "quantile"
	output$outFile = "analysis_output.xls"
	
    #dirData = getwd()
	#params = dir(dirData)[grep("arams",dir(dirData))] #read params file
    #params = read.delim(paste(dirData,params,sep="/"), header=F)
	params <- read.delim(file,header=F)
	rownames(params) = params[,1]
    #---------------------------------------------------------------------------
	output$project             <- as.character(params["Project",2])
	output$dirCommon           <- as.character(params["CommonDirectory",2])
	output$geneAnnFile         <- as.character(params["GeneAnnotationFile",2])
	output$version             <- as.character(params["ChipVersion",2])
	if(!is.na(as.character(params["ContrastFile",2]))) {
		output$contrastsFile  <- as.character(params["ContrastFile",2]) }
	if(!is.na(as.character(params["SampleAnnotationFile",2]))) {
		output$targetsFile    <- as.character(params["SampleAnnotationFile",2]) }
	if(!is.na(as.character(params["Normalization",2]))) {
		output$normalizationMethod <- as.character(params["Normalization",2]) }
	if(!is.na(as.numeric(as.character(params["Background",2])))) {
		output$bg        <- as.numeric(as.character(params["Background",2])) }
	if(!is.na(as.character(params["AnalysisOutFile",2]))) {
		output$outFile        <- as.character(params["AnalysisOutFile",2]) }
	if(!is.na(as.character(params["Author",2]))) {
		output$author        <- as.character(params["Author",2]) }	else { 
			output$author = "NA"
		}
	if(output$author == "<Your Initials>" || output$author=="NA") {
		print("Warning! You must specify your initials in : param$author ...")
		print("Enter your initials :")
		output$author = scan(what="character", n=1)
	}
    #-------------------------------------------------------------------------
	output
}


# function to read BeadType files
# All files need to be directly into chips directory.
## TODO : implement within readBeads.
my.readBeadType <- function(chips.dir=file.path("chips"))
{
	
 	 data=lapply( list.files(path=chips.dir,full.names=TRUE,pattern="beadTypeFile.txt"),function(fn)
         {
           mat = as.matrix(read.table(fn,sep=",",header=TRUE,row.names="Illumicode",colClasses="numeric")[,"Mean.GRN",drop=FALSE])
           colnames(mat) = gsub( file.path(chips.dir,""),'',gsub('_beadTypeFile.txt','',fn))
          return(mat)
        })
	if(any(  sapply(data,function(mat){any(rownames(mat)!=rownames(data[[1]]))})   )){stop("Something wrong with rownames")}
        data = do.call(cbind,data)
	eset = new("ExpressionSet",exprs=data)
	return(eset)
}

	
# color.me : Return a color string vector corresponding to x, with a different
#            colour for each type of value provided in x.
color.me <- function(x)
{
	x = factor(x)
	levels(x) = rainbow(nlevels(x)) # rainbow by default
	if(nlevels(x) == 2 )
	{
		levels(x) = c("blue","yellow")
	}
	if(nlevels(x) == 4 )
	{
		levels(x) = c("blue","cyan","yellow","magenta")
	}
	
	return(as.character(x))
}


# 
plotDensitiesMat <- function (mat, subset = c(1:dim(mat)[2]), title = NULL,
legend.cex = 0.95,legend.order =NULL)
{
    val <- mat
    val <- log2(val)
    if(length(subset)==2)  colors <- c("blue","red")
    else  colors <- rainbow(length(subset))
	
    y.max <- c()
    y.min <- c()
    x.max <- c()
    x.min <- c()
    legend.txt <- c()
    for (n in subset) {
        y.max[n] <- max(density(val[, n], na.rm = TRUE)$y)
        y.min[n] <- min(density(val[, n], na.rm = TRUE)$y)
        x.max[n] <- max(density(val[, n], na.rm = TRUE)$x)
        x.min[n] <- min(density(val[, n], na.rm = TRUE)$x)
    }
    ymax <- max(y.max)
    ymin <- min(y.min)
    xmax <- max(x.max)
    xmin <- min(x.min)
	
    for (n in 1:ncol(val)){
		
        if (n == 1){
            par(mgp=c(1,0,0))
            legend.txt <- c(legend.txt, colnames(mat)[n])
			
            plot(density(val[, n], na.rm = TRUE), col = colors[n],
				 main = "",xlab="x",ylab="y",xlim=c(xmin,xmax),ylim=c(ymin,ymax)) }
        else {
			lines(density(val[, n], na.rm = TRUE), col = colors[n])
			legend.txt <- c(legend.txt, colnames(mat)[n]) }
    }
    if (is.null(title))
	title(paste("Density Plot"))
    else title(title)
	
    legend.ncol <- 1
	
    if (length(legend.txt)>25){
		legend.cex  <- 0.7
		legend.ncol <- ceiling(length(legend.txt)/25)
    }
	
    legend(x = "topright", legend = legend.txt, cex = legend.cex,
		   fill = colors, inset = 0.01,bty="n", ncol = legend.ncol)
}


# mergeEsets : merge two ExpressionSets together
# The second set is appended to the first set
# Warning : Please make sure both sets have the same probe number
#           AND same number of phenotypic column in the same order.
# x : first ExpressionSet
# y : second ExpressionSet
# output : merged ExpressionSet
mergeEsets <- function(x,y) {
	output <- x
#y      <- y[featureNames(x),]
	exprs(output) <- cbind(exprs(x),exprs(y)[rownames(exprs(x)),])
	pData(output) <- rbind(pData(x),pData(y))
	output 
}


# variable for versioning system
illupipe.version = "2.0.0.3"

