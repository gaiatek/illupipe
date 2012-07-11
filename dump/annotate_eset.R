# * Introduction *
#
# This function implements the already natural "sample VS chip" paradigm of a microarray experiment. It reads
# the hyb and sample annotations from a 2-table (excel) database and writes it to the pData. It also
# reads the feature Annotation using getAnnotation() and puts it into fData. It also performs a few consistency checks.
# The purpose of this function is to facilitate the correction, uniformization,
# of large merged projects where the annotation is unstable. Navigating, inspecting, sharing between analysts,
# merging, augmenting ,modifying and keeping track of changes in the  the annotation is simply
# user-friendlier using excel than with the pData.
#
#
# * Input *
# 
# - eset : An ExprSet Object with one column names 'ArrayID'. Preexisting pData and fData are overwritten.
#
# - worksheet : xls spreadsheet (the 'master spreadsheet') where the first two sheets hold the sample and hybs annotations.
#   The *first* sheet holds the sample annotation which minimally has one column named 'SampleID' that lists
#   the  (unique) sampleids.
#   The *second* sheet lists the hybridizations. Column 'ArrayID' and 'SampleID' are mandatory.
#
# - f.ann : a data frame to be put in fData.
#
# - (NOT USED) silent.drop.of.hybs : In case th eset contains hybs that are not the the hyb sheet
#             , should they be dropped when merging
#   
# - stricthybs : Usually hybridizations in the raw eset should be the same as in the pData. This is not always
#                the case so  stricthybs=TRUE relaxes this assumption. 
#
#  -prefix : whether or not s. and h. should be added to column names
#
#  - third.id : if not NULL, the name of column in the sample annotation which can be used to cross with a third table found
#               in the third tab of the master sheet. For instance for large projects it is useful to
#               normalize the annotation DB by having a third table which
#               lists attributes of the donors/patients.
#
#
# * Notes *
#
# (1) greek characters like \mu are not supported by read.xls(). Please replace them
# using excel before calling the function. (2) Control samples are usually not given a unique ID by the vgti. You can
# either add an empty line in the sample annotation with a unique id (e.g. 'ambion3') or remove them from the hyb sheet
# alltogether. (3) All samples or hybs in th eset MUST be present in the annotation, but the reverse is not true.
# (4) to avoid problems in excel paste values only, not formulas... (5) this function returns an eset with a fully
# sample/hybs joined annotation.
# 
#
#
# * Recommendations *
#
#  (1) I recommend the following columns in the database :
#
# 'Set' : One setID per sample set or hyb set. Keeping track of the sets makes the decision of what is
#  a batch or not easier. Differetn set ids can be higlihted in different colors.
#
# 'Note' : put this as the first column of the sample and hyb sheets. Note any modification to the annotation, e.g.
#  "tissue changed from BAL to LN" or "inversion with SampleID B4 confirmed by Angela" etc.
#  'outlier.flag' : this column could track which samples are outliers, accompined by an explanation for removal.
#                   This way we could get rid of 'outliers.txt'.
#
#
#  (2) For tracability purposes, I also recommend keeping the original annotation and appending you own,
#  analysis-ready columns, ideally pasting those among the first columns of the sheet.
#  For instance, if the vgti sends a worksheet with the column "TIME POINT (dAYS)", paste a version
#  renamed "Time" as one of the first columns. 
#
#  (3) In my opinion, 'batch' and 'Class' should not be in the worksheet. The reason is that what consitutes a 'batch' or a
#  'Class' is really an analysis decision, not a characteristic of the data itself.
#
#  (4) This function will return the fully joined sample and hyb annotation. Simply subset and re-order to get
#      nice annot table in the report.
#
# * TO DO *
#
# -Second form normalization by adding an optional "Donor table" to database schema, useful
#   for large cohort studies where donor annotation is redundant in the sample table.
# - Try out sourcing the master spreadsheet from Google Docs. Would be great in terms of versioning and sharability.
# - add an option to remove the 's.' and 'h.' which can be annoying.

annotate.eset <- function(
                          eset
                          ,worksheet=file.path(dirData,"worksheets","mastersheet.xls")
                          ,f.ann=getAnnotation()
                          #,silent.drop.of.hybs = FALSE
                          ,stricthybs=TRUE
                          ,prefix=TRUE
                          ,third.id=NULL
                          )
  {

    ## First drop everyting in the pData apart from ArrayID
    sn                  = sampleNames(eset) # to retain the previous ordering of samples
    pData(eset)$ArrayID = sn
    pData(eset)         = pData(eset)[,"ArrayID",drop=FALSE] # drop everything other than ArrayID.
    
    ## Re-cross with masterprint, don't want to re-read the data every time I find a typo in the sheet
    masterprints = read.xls(worksheet, sheet = 2,stringsAsFactors=FALSE) # row name will check for unicity

    ## Rename the field to keep track of whih table they come from
    if(prefix)
      {
        colnames(masterprints) = paste('h',colnames(masterprints),sep=".")
        colnames(masterprints)[colnames(masterprints)=="h.ArrayID"] = "ArrayID" 
        colnames(masterprints)[colnames(masterprints)=="h.SampleID"] = "SampleID" 
      }

    ## Check unicity of hybs
    if(anyDuplicated(masterprints[["ArrayID"]])){stop("Duplicated ArrayID. Check your master hyb sheet")}

    ## Check for completness of Hyb sheet
    if(any( ! sn  %in% masterprints[["ArrayID"]]   )  ){ # & !silent.drop.of.hybs
      stop("Some hybs in eset are not in hyb sheet and I don't like silent dropping of hybridizations.")}

    ## Strict hybs check
    if( !setequal( pData(eset)$ArrayID,masterprints$ArrayID )  & stricthybs) {
      stop("Some hybs in the hyb sheet are not in the eset and strict=TRUE")} 

    ## Now cross
    pdata  = merge(pData(eset),masterprints,by="ArrayID",sort = FALSE)# for every rowname in x, get the lines from y
    rownames(pdata) = pdata$ArrayID
    pData(eset) = pdata[sn,] # they should laready sorted but let's make sure

  
    
    ## Cross with sample information
    mastersamples = read.xls(worksheet, sheet = 1,stringsAsFactors=FALSE) # this sheet is the result of merging all sample desc and hyb sheets
    
    # We rename the field to keep track of whih table they come from
    if(prefix)
      {
        colnames(mastersamples) = paste('s',colnames(mastersamples),sep=".")
        colnames(mastersamples)[colnames(mastersamples)=="s.SampleID"] = "SampleID"
        if(!is.null(third.id)){colnames(mastersamples)[colnames(mastersamples) == paste("s.",third.id,sep="")  ] = third.id }
      }

    ## Here check for unicity and completness of sample annots
    if(anyDuplicated(mastersamples[["SampleID"]])){stop("Duplicated SampleID. Check your master sample sheet")}
    if(   any(! pData(eset)[["SampleID"]]   %in% mastersamples[["SampleID"]]) ){stop("Some SampleID in hyb sheet not in sample sheet. Check for completness.")}

    ## Now cross
    pdata = merge(pData(eset), mastersamples,by="SampleID",sort = FALSE)
    rownames(pdata) = pdata$ArrayID
    pData(eset) = pdata[sn,]

    if( any(colnames(exprs(eset)) != sampleNames(eset)) | any(rownames(pData(eset)) != sampleNames(eset))  )
      {stop("Whoooah. I messed up")}

    ## Next annotate the features
    if( any(! featureNames(eset) %in% rownames(f.ann) )  ) # featureNames must be a strict subset of the annot
      {
        stop("Some of the eset features are not in the supplied annotation")
      }   
    fData(eset) = f.ann[featureNames(eset),]  ## because featureNames could be a subset of annot id


    ## Third table
    if(!is.null(third.id))
      {
        masterdonors             = read.xls(worksheet, sheet = 3,stringsAsFactors=FALSE)
        rownames(masterdonors)   = masterdonors[[third.id]]
        masterdonors[[third.id]] = NULL
        if(prefix){colnames(masterdonors) = paste("d.",colnames(masterdonors),sep="")} ## modifies colnames to fit
        if( any(! pData(eset)[[third.id]] %in% rownames(masterdonors) )  ){stop("Some IDs in sample table not in third table")}
        pData(eset) = cbind( pData(eset), masterdonors[pData(eset)[[third.id]],] )
      }
  
    return(eset)
  }
