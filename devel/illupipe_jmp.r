###############################################################################
# illupipe_jmp.r  v.0.0.0.1                                                   #
#-----------------------------------------------------------------------------#
# Illumina MicroArray Pipeline Toolkit                                        #
# Functions for data import/export of illupipe datasets to/from JMP Genomics  #
#-----------------------------------------------------------------------------#
# @author : Jean-Philippe Goulet                                              #
# @date   : August 25th  2011                                                 #
###############################################################################


# genesColumns : NULL or an index vector representing the column numbers to
#                include into fit$genes, eg. 1:15
importJMP_fit <- function(fn='file.txt', featureNames='Array_Address_Id',
                          genesColumns=NULL, contrastNames=NULL) {
  if(!file.exists(fn)) { stop("Cannot find input file.")}
  temp <- read.delim(fn)
 
  
  fit <- list()
  if(is.null(genesColumns)) {
    i <- grep('Variance',colnames(temp))
    fit$genes <- temp[,1:(i[1] - 1)]
  }
  else {
    fit$genes <- temp[,genesColumns]
  }
  rownames(fit$genes) <- fit$genes[[featureNames]]
  if(featureNames != 'ProbeID') {
    colnames(fit$genes) <- sub(featureNames,'ProbeID',colnames(fit$genes))
  }
  colnames(fit$genes) <- sub('Symbol','SYMBOL',colnames(fit$genes))
  
  fit$p.value <- temp[,grepl('p.Value',colnames(temp))]
  i <- grep('Adjusted',colnames(fit$p.value))
  if (length(i) > 0) {
    fit$fdr <- fit$p.value[,i]
    fit$p.value <- fit$p.value[,1:(i[1] -1)]
  }
  fit$coefficients <- temp[,grepl('Diff',colnames(temp))]
  fit$coefficients <- fit$coefficients[,1:ncol(fit$p.value)]
  fit$t <- temp[,grepl('t.Statistic',colnames(temp))]
  
  rownames(fit$p.value) <- fit$genes$ProbeID
  rownames(fit$fdr) <- fit$genes$ProbeID
  rownames(fit$coefficients) <- fit$genes$ProbeID
  rownames(fit$t) <- fit$genes$ProbeID

  fit$p.value <- apply(fit$p.value,c(1,2),as.numeric)
  if(!is.null(fit$fdr)) { fit$fdr <- apply(fit$fdr,c(1,2),as.numeric) }
  if(!is.null(fit$t)) { fit$t <- apply(fit$t,c(1,2),as.numeric) }
  fit$coefficients <- apply(fit$coefficients,c(1,2),as.numeric)

  if(!is.null(contrastNames) && (length(contrastNames) != ncol(fit$p.value))) {
    stop('Wrong number of contrast names defined by <contrastNames>.')
  }
  else if(!is.null(contrastNames)) {
    colnames(fit$p.value) <- contrastNames
    colnames(fit$fdr) <- contrastNames
    colnames(fit$t) <- contrastNames 
    colnames(fit$coefficients) <- contrastNames 
  }
  fit <- new('MArrayLM', list(coefficients=fit$coefficients, p.value=fit$p.value, genes=fit$genes, fdr=fit$fdr, t=fit$t))
  return(fit)
}

exportJMP_fit <- function() {

}

importJMP_eset <- function() {

}

exportJMP_eset <- function() {

}


