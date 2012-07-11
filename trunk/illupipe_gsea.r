###############################################################################
# illupipe_gsea.r  v.0.55                                                     #
#-----------------------------------------------------------------------------#
# Illumina MicroArray Pipeline Toolkit                                        #
# Functions for pathway analyses of illumina beadArrays transcriptomics       #
#-----------------------------------------------------------------------------#
# @author : Jean-Philippe Goulet                                              #
# @date : September 10th 2011                                                 #
###############################################################################

# gsea.input : Generate input files for GSEA analyses
# method        : 'ranked' (ranked gene list) or 'exp' (samples permutation). 
# lm            : linear model.
# chip.version  : filename used to map probeID to SYMBOL.
# contrasts     : a vector of contrasts.
# adjust.method : method used to adjust p-values for topTable. <deprecated>
# dirCommon     : directory where is located folder 'GeneSetDatabases'.
# dirPath       : Needs to be the absolute path to your gsea folder in order 
#                 for gsea.jar to recognize the file locations.
gsea.input <- function(eset, dirPath=file.path(dirData,"gsea"), 
				dirCommon=param$dirCommon, method="ranked", 
				lm=fit2, contrasts=colnames(lm),  
		    adjust.method="fdr", chip.version="illuminaV3.chip") {

  

  if(method=='exp') dat <- exprs(eset)
  
  topTableList <- list()
  for (i in 1:length(contrasts)) {
        #	topTableList[[i]] <- topTable(lm, coef=i, adjust=adjust.method, 
	#	 number=nrow(lm))
    temp = list()
    temp$t = lm$t[,i]
    temp$ProbeID <- lm$genes$ProbeID
    topTableList[[i]] <- temp
  }

  dirPath_list <- c("C2","C3","C5")
  DB <- c("2", "3", "5")
  
  ##Shell script file//Batch file
  gsea_bat=file(file.path(dirPath, "fileSH"),"w")

	#######################################################################
  for (icont in 1:length(contrasts)){
    ## preparing .gct and .cls files
    arr1 <- strsplit(contrasts[icont],"-")[[1]][1]
    indexarr1 <- which(pData(eset)$Class==arr1)
    arr2 <- strsplit(contrasts[icont],"-")[[1]][2]
    indexarr2 <- which(pData(eset)$Class==arr2)

     #####gsea parameters list and batch commands file
    ParamGsea <- "params_gsea_smal.txt"
    params    <- read.delim(paste(dirPath,ParamGsea,sep="/"), header=F)
    gsea_Params_list <- list()
    gsea_Params_list[[1]] <- file(paste(dirPath, "/",
           strsplit(contrasts[icont], "-")[[1]][1],".VS.", strsplit(contrasts[icont],
           "-")[[1]][2], "_", "C2_Params", ".txt", sep=""), "w")
    gsea_Params_list[[2]] <- file(paste(dirPath, "/", 
           strsplit(contrasts[icont], "-")[[1]][1], ".VS.", 
           strsplit(contrasts[icont], "-")[[1]][2], "_", "C3_Params", ".txt",
 		 sep=""), "w")
    gsea_Params_list[[3]] <- file(paste(dirPath, "/",
	   strsplit(contrasts[icont], "-")[[1]][1], ".VS.",
           strsplit(contrasts[icont], "-")[[1]][2], "_", "C5_Params", ".txt",
		 sep=""), "w")
		
		##.gct file
    if(method=="exp") # more complex contrasts names make this part crash
      {
       	
        data2 <- data.frame(dat[,rownames(pData(eset))[c(indexarr1,indexarr2)]])
        temp <- rep(NA,length(data2[,1]))
        dataF <- data.frame(dimnames(data2)[[1]],temp,data2)
        names(dataF) <- c("NAME","Descrption",dimnames(data2)[[2]])
        Datagct <- file(paste(dirPath,"/", 
                              strsplit(contrasts[icont], "-")[[1]][1],".VS.", 
                              strsplit(contrasts[icont],"-")[[1]][2],".gct",sep=""),"w")
        cat("#1.2",file=Datagct)
        cat("\n",file=Datagct)
        cat(dim(dat)[1],"\t",length(indexarr1)+length(indexarr2),file=Datagct)
        cat("\n",file=Datagct)
        close(Datagct)
        coTY <-read.table(file=paste(dirPath,"/",
                            strsplit(contrasts[icont], "-")[[1]][1], ".VS.",
                            strsplit(contrasts[icont],"-")[[1]][2],".gct",sep=""), header=F, sep="\t")
        write.table(dataF, paste(dirPath, "/",
                                 strsplit(contrasts[icont], "-")[[1]][1], ".VS.",
                                 strsplit(contrasts[icont],"-")[[1]][2], ".gct", sep=""), append=T,
                    row.names=F,col.names=T,sep="\t")
      }
    
    ##.cls file
    gsea_Pheno <- file(paste(dirPath, "/",
                             strsplit(contrasts[icont], "-")[[1]][1], ".VS.", 
                             strsplit(contrasts[icont], "-")[[1]][2], ".cls", sep=""), "w")
    cat(length(indexarr1)+length(indexarr2),2,1,file=gsea_Pheno,sep=" ")
    cat("\n",file=gsea_Pheno)
    cat("#", strsplit(contrasts[icont], "-")[[1]][1],
        strsplit(contrasts[icont], "-")[[1]][2], file=gsea_Pheno,sep=" ")
    cat("\n",file=gsea_Pheno)
    cat(c(rep(1,length(indexarr1)),
          rep(0,length(indexarr2))),file=gsea_Pheno,sep=" ")
    cat("\n",file=gsea_Pheno)
    close(gsea_Pheno)
    #####gsea rankedlist 

    datarank_list <- cbind("#"=topTableList[[icont]]$ProbeID, 
                           " "=topTableList[[icont]]$t)
    names(datarank_list) <- c("#"," ")
    write.table(datarank_list,
                file=paste(dirPath,"/",
                  strsplit(contrasts[icont],"-")[[1]][1],".VS.",
                  strsplit(contrasts[icont],"-")[[1]][2],".rnk",sep=""),
                row.names=F,quote=F,sep="\t")


	
    if(method=="exp"){
      for(igsea in 1:3){
        cat("gmx", "\t", paste(dirCommon, "/GeneSetDatabases/c",
                               DB[igsea], ".all.v2.5.symbols.gmt", sep=""), 
            file=gsea_Params_list[[igsea]])
        cat("\n",file=gsea_Params_list[[igsea]])
        cat("res","\t",paste(dirPath,"/",
                             strsplit(contrasts[icont],"-")[[1]][1],".VS.", 
                             strsplit(contrasts[icont],"-")[[1]][2],".gct",sep=""), 				 file=gsea_Params_list[[igsea]])
        cat("\n",file=gsea_Params_list[[igsea]])
                         #cat("cls", "\t", paste(dirPath, "/",
                                # strsplit(contrasts[icont],"-")[[1]][1], ".VS.", 
				# strsplit(contrasts[icont],"-")[[1]][2],".cls","#", 
				# strsplit(contrasts[icont],"-")[[1]][1],"_versus_", 
				# strsplit(contrasts[icont],"-")[[1]][2],sep=""), 
				# file=gsea_Params_list[[igsea]])
        cat("cls", "\t", paste(dirPath, "/",
                               strsplit(contrasts[icont],"-")[[1]][1], ".VS.", 
                               strsplit(contrasts[icont],"-")[[1]][2],".cls",sep=""), 
            file=gsea_Params_list[[igsea]])
        cat("\n",file=gsea_Params_list[[igsea]])
        cat("chip","\t", paste(dirCommon,
                               "/GeneSetDatabases/", chip.version, sep=""), 
            file=gsea_Params_list[[igsea]])
        cat("\n",file=gsea_Params_list[[igsea]])
        cat("out","\t",paste(dirPath,"/",dirPath_list[igsea],sep=""), 				file=gsea_Params_list[[igsea]])
        cat("\n",file=gsea_Params_list[[igsea]])
        cat("rpt_label","\t", paste(
                                    strsplit(contrasts[icont],"-")[[1]][1],".VS.", 
                                    strsplit(contrasts[icont],"-")[[1]][2], sep=""), 				 file=gsea_Params_list[[igsea]])
        cat("\n",file=gsea_Params_list[[igsea]])
        for( i in 1:length(params[,1])){
          cat(paste(as.character(params[i,1]), "\t",
                    as.character(params[i,2]), sep=""),
              file=gsea_Params_list[[igsea]])
          cat("\n",file=gsea_Params_list[[igsea]])
        }
        close(gsea_Params_list[[igsea]])
        cat("java -cp gsea2.jar -Xmx512m xtools.gsea.Gsea -param_file",
            paste(strsplit(contrasts[icont],"-")[[1]][1],".VS.", 
                  strsplit(contrasts[icont],"-")[[1]][2],"_", 
                  dirPath_list[igsea],"_Params",".txt\n",sep=""), file=gsea_bat)
      }
    }else{
      cat("gmx","\t",paste(dirCommon,"/GeneSetDatabases/c",
                           DB[1],".cp.v3.0.symbols.gmt",sep=""), 
          file=gsea_Params_list[[1]])
      cat("gmx","\t",paste(dirCommon,"/GeneSetDatabases/c",
                           DB[2],".tft.v3.0.symbols.gmt",sep=""), 
          file=gsea_Params_list[[2]])
      cat("gmx","\t",paste(dirCommon,"/GeneSetDatabases/c",
                           DB[3],".bp.v3.0.symbols.gmt",sep=""), 
          file=gsea_Params_list[[3]])
      
      for(igsea in 1:3){
				# changed for MSigDB v.3 (cf.3-10 lines up)
 				#cat("gmx","\t",paste(dirCommon,"/GeneSetDatabases/c",
				# DB[igsea],".v3.0.symbols.gmt.txt",sep=""), 
				# file=gsea_Params_list[[igsea]])
				
        cat("\n",file=gsea_Params_list[[igsea]])
        cat("rnk","\t",paste(dirPath,"/",
                             strsplit(contrasts[icont],"-")[[1]][1],".VS.", 
                             strsplit(contrasts[icont],"-")[[1]][2],".rnk",sep=""), 				 file=gsea_Params_list[[igsea]])
        cat("\n",file=gsea_Params_list[[igsea]])
        cat("chip","\t", paste(dirCommon,
                               "/GeneSetDatabases/", chip.version, sep=""), 
            file=gsea_Params_list[[igsea]])
        cat("\n",file=gsea_Params_list[[igsea]])
        cat("out","\t",paste(dirPath,"/",dirPath_list[igsea],sep=""), 				 file=gsea_Params_list[[igsea]])
        cat("\n",file=gsea_Params_list[[igsea]])
        cat("rpt_label","\t", paste(
                                    strsplit(contrasts[icont],"-")[[1]][1],".VS.", 
                                    strsplit(contrasts[icont],"-")[[1]][2],sep=""),
            file=gsea_Params_list[[igsea]])
        cat("\n",file=gsea_Params_list[[igsea]])
        for  ( i in 1:length(params[,1])){
          cat(paste(as.character(params[i,1]), "\t",
                    as.character(params[i,2]),sep=""),
              file=gsea_Params_list[[igsea]])
          cat("\n",file=gsea_Params_list[[igsea]])
        }
        close(gsea_Params_list[[igsea]])
        cat("java -cp gsea2.jar -Xmx512m xtools.gsea.GseaPreranked -param_file", 
            paste('\"',strsplit(contrasts[icont],"-")[[1]][1],".VS.", 
                  strsplit(contrasts[icont],"-")[[1]][2],"_", dirPath_list[igsea], 
                  "_Params",".txt\"\n",sep=""),file=gsea_bat)
      }
    }
  }
  #cat("\n",file=gsea_bat)
  #cat("exit0",file=gsea_bat)
  close(gsea_bat)
}


# gsea.report : generates a HTML report that summarizes gsea analyses.
# eset    : the corresponding eset for visualization
#     *** You need to set pData(eset)$Class corresponding to contrasts. ***
# dirPath : path to gsea folder.
# dirCommon     : directory where is located folder 'GeneSetDatabases'.
# contrasts : vector of contrasts.
# chip.version  : filename used to map probeID to SYMBOL.
# archive : Logical. Does it save .Rdata files containing the top significant geneSets 
#           with their enriched genes.
# visualize : Logical. Defaults to TRUE : generate visualization plots. cf. gsea.visualize()
gsea.report <- function(eset, dirPath= file.path(dirData,"gsea"),
		dirCommon=param$dirCommon, contrasts=getContrasts(), visualize=T,
		chip.version="illuminaV3.chip", mode="auto", inputFiles=NULL, archive=F) {
	
	vecCont = contrasts
	DB = c("C2", "C3", "C5")
	#method = "Preranked"

	c2_list = dir(file.path(dirPath, "gsea", "C2"))
	c3_list = dir(file.path(dirPath, "gsea", "C3"))
	c5_list = dir(file.path(dirPath, "gsea", "C5"))

	vecCont = unlist( lapply(strsplit(vecCont,split="-"),function(e){paste(e[1:2],collapse="-")} ))	
	# 
	vecCont = sub("-",".VS.", vecCont, fixed=T)
	
	##############################
	## Generate Report Title
	##############################
	#read template 
	temp_file = read.delim(file.path(dirPath,"gsea_template.html")) 

	tmpfic = HTMLInitFile(file.path(dirPath,"gsea"), filename="index",
		CSSFile="report.css", 
		Title="MicroArray Data Gene Set Enrichment Analysis Summary")


	# Report Header with date
	HTML(paste("<table style=\"width: 100%; 
		background-color: rgb(210, 218, 254);\"><tbody><tr>  		<td><span style=\"font-weight: bold;\">Pathways Analysis 
		Report</td><td style=\"text-align: right;\">",
		format(Sys.time(), "%d-%b-%Y"), "</td></tr></tbody></table>"),
		file=tmpfic)
		
	# Report GSEA description (template)	
	for(i in 1:nrow(temp_file)) {
		cat(paste(as.character(temp_file[i,1]),"&nbsp;",sep=""),file=tmpfic, 
			append=T)
	}		
		


	##############################
	## Generate HTML Table Header
	##############################

	HTML("<table cellpadding=\"7\">
	  <tbody>
   	 <tr class=\"head\">
      <td colspan=\"4\" rowspan=\"1\"><big><big>
      <span style=\"font-weight: bold;\">Pathways</span></big></big>
      </td>
    </tr>
    <tr class=\"line1\">
      <td style=\"width: 25%;\"> <span
 	style=\"font-weight: bold;\">Contrast</span></td>
      <td style=\"font-weight: bold; width: 25%; text-align: center;\">C2:
	Canonical Pathways</td>
      <td style=\"width: 25%; text-align: center; font-weight: bold;\">C3:
	Transcription Factors<br>
      </td>
      <td style=\"width: 25%; text-align: center; font-weight: bold;\">C5
	: Biological Processes <br>
      </td>
    </tr>
	",file=tmpfic)
	
	contrastID2 = list()
	contrastID3 = list()
	contrastID5 = list()
	
	line_color = c("line1","line2","line3","line4")
	j=1
	for (i in 1:length(contrasts)){
		contrastID2[[i]] = strsplit(c2_list[grep(vecCont[i],c2_list)],".",
			fixed=T)[[1]][5]
		contrastID3[[i]] = strsplit(c3_list[grep(vecCont[i],c3_list)],".",
			fixed=T)[[1]][5]
		contrastID5[[i]] = strsplit(c5_list[grep(vecCont[i],c5_list)],".",
			fixed=T)[[1]][5]
		j <- j +1
		if(j>4) { j = 1}
	
		#####################################
		## First Column : contrast name + C2
		#####################################
	
		HTML(paste("<tr class=", line_color[j], "><td>", contrasts[i], 
		"</td><td><ul><li><a href=\"visualization_plots/Pathways_C2_", vecCont[i],
		".png\">Summary Figure</a></li><li> <a href=\"C2/", 
		c2_list[grep(vecCont[i],c2_list)], "/gsea_report_for_na_pos_", 
		contrastID2[[i]],".html\">upregulated pathways</a></li><li>
		<a href=\"C2/", c2_list[grep(vecCont[i],c2_list)], "/gsea_report_for_na_neg_",
		contrastID2[[i]], ".html\">downregulated
		pathways</a></li></ul></td>",sep=""), file=tmpfic )
	
	
		#####################################
		## 	Second Column : C3
		#####################################
		HTML(paste("<td><ul><li><a href=\"visualization_plots/Pathways_C3_", vecCont[i],
		".png\">Summary Figure</a></li><li> <a href=\"C3/",
		c3_list[grep(vecCont[i],c3_list)], "/gsea_report_for_na_pos_",
		contrastID3[[i]],".html\">upregulated pathways</a></li><li>
		<a href=\"C3/",  c3_list[grep(vecCont[i],c3_list)], 
		"/gsea_report_for_na_neg_", contrastID3[[i]], 
		".html\">downregulated
		pathways</a></li></ul></td>",sep=""),file=tmpfic)
	 

		#####################################
		## Third Column : C5
		#####################################
		HTML(paste("<td><ul><li><a href=\"visualization_plots/Pathways_C5_", vecCont[i],
		".png\">Summary Figure</a></li><li> <a href=\"C5/",
		c5_list[grep(vecCont[i],c5_list)], "/gsea_report_for_na_pos_",
		contrastID5[[i]],".html\">upregulated pathways</a></li><li>
		<a href=\"C5/", c5_list[grep(vecCont[i],c5_list)], 
		"/gsea_report_for_na_neg_", contrastID5[[i]], 
		".html\">downregulated
		pathways</a></li></ul></td></tr>",sep=""),file=tmpfic)
		 


	}

	## Table Footer
	HTML("  </tbody></table><br><br>", file=tmpfic)

			
	## generate Visualization plots
	gsea.visualize(eset,chip.version=chip.version,mode=mode,inputFiles=inputFiles,
				   contrasts=contrasts, dirPath=file.path(dirPath,"gsea"),dirCommon=dirCommon , archive=archive)

}

# gsea.vizualise : Generate Vizualisation plots for top 5 pathways and
#                  top 5 genes for up and down-regulated pathways.
# eset     : the ExpressionSet to use for heatmaps.
# mode     : Visualization mode -> 'auto' or 'custom'. 
#            'auto' uses systematically top 5 pathways UP and DOWN + top 5 genes
#            from the leading edge subsets. 
#            'custom' uses selected genes and pathways from inputFiles. 
# dirPath  : Path to the C2,C3 and C5 folders (and inputFiles, if mode=='custom')
# dirCommon: Path to the annotations files, where is 'GeneSetDatabases' folder.
# contrasts: Vector of contrasts you want to visualize. for each contrasts,
#            there must be a GSEA report for C2, C3 and C5 in order to work.
# chip.version: What file for probe-to-symbol mapping needs to be used.
#               This file needs to be inside dirPath/GeneSetDatabases/.
# inputFiles: List of filenames containing selected genes and pathways for custom
#             visualization mode. Must be in this order : C2,C3,C5. 
#             See 'dump/C2_MDC_InputPathways-example.txt' for a file example.
# archive : Logical. Does it save .Rdata files containing the top significant geneSets 
#           with their enriched genes.
# output   : PDFs are stored in dirPath (for every contrasts : C2,C3,C5).
gsea.visualize <- function(eset, dirPath=file.path(dirData,"gsea","gsea"), 
	dirCommon=param$dirCommon, contrasts=getContrasts(), 
	chip.version="illuminaV3.chip", mode="auto", inputFiles=NULL, archive=F) {
	
  MyMapFile=read.delim(file.path(dirCommon,"GeneSetDatabases", chip.version))

  filter <- 0
  if(mode == "auto") {
    selectGene <- "Top"
    inputPath <- "none"
  }
  if(mode == "custom") {
    selectGene <- "custom"
    inputPath <- "custom"
    #MyInputFileC2=read.delim(file=paste(dirPath,"/C2.txt",sep=""))
    #MyInputFileC3=read.delim(file=paste(dirPath,"/C3.txt",sep=""))
    #MyInputFileC5=read.delim(file=paste(dirPath,"/C5.txt",sep=""))
    #MyInputFileList=list(MyInputFileC2,MyInputFileC3,MyInputFileC5)
	 MyInputFileList <- as.list(inputFiles)  
  }

  vecCont <- contrasts
  vecCont = unlist( lapply(strsplit(vecCont,split="-"),function(e){paste(e[1:2],collapse="-")} ))
  dirPath_list=c("C2","C3","C5")

  for (iij in 1:length(dirPath_list)) {
	for ( icont in 1:length(vecCont)) {
	  Path.list.all <- list()
	  Path.list.all.bkup <- list()
	  PhenoVecAll <- list()
	  dirGsea <- paste(dirPath,"/",dirPath_list[iij],sep="")
	  setwd(dirGsea)
	  gseafolder <- dir()
	  folder <- gseafolder[grep(paste(strsplit(vecCont[icont],"-")[[1]][1],".VS.",strsplit(vecCont[icont],"-")[[1]][2],sep=""),gseafolder,fixed=T,perl=F)]
	
	  setwd(as.character(paste(dirGsea,"/",folder,sep="")))
	  result= dir()

	  ## import   pathways from the Input file
	  #######################################################################
	  #path_names_pos
	  if(inputPath=="custom") {
	    path_names_na <- MyInputFileList[[iij]][MyInputFileList[[iij]]$contrast
		== vecCont[icont] & MyInputFileList[[iij]]$updown == "down",			"pathways"] 
	  }else{
	    file0na <- result[grep("gsea_report_for_na_neg_",result,perl=T)]
	    file1na <- file0na[grep(".xls",file0na,perl=T)]
	    gsea_report_na <- read.delim(file=file1na,sep="\t")
		if(length(gsea_report_na$NAME) > 40) {
		  path_names_na <- gsea_report_na$NAME[1:40]
	  
		}
		else {
			path_names_na <- gsea_report_na$NAME
		}
	  }
	  if (length(path_names_na)!=0) {
	      Path.list <- list()
		  Path.list.bkup <- list()
	      PhenoVec <- list()
	      for(jh in 1:length(path_names_na)) {
		pathfile <- read.delim(paste(as.character(path_names_na[jh]),
		  ".xls",sep=""),sep="\t")
		yesno <- which(as.character(pathfile[,8])=="Yes")
		Path.list[[jh]] <- as.character(pathfile[,2])[yesno]
		Path.list.bkup[[jh]] <- Path.list[[jh]]
		if(selectGene=="Top"){
		  if((length(Path.list[[jh]])>5)&&(length(path_names_na) >=1)) {
			Path.list[[jh]] <- rev(Path.list[[jh]])[1:5]
		  }
		}
		if(selectGene == "custom") {
			list_genes_na <- MyInputFileList[[iij]][MyInputFileList[[iij]]$contrast
			 == vecCont[icont]&MyInputFileList[[iij]]$pathways == 
			 as.character(path_names_na[jh]),"genes"]   
		   if(length(unlist(strsplit(as.character(list_genes_na),",")))!=0) {  
				Path.list[[jh]] <- unique(c(Path.list[[jh]], 
					unlist(strsplit(as.character(list_genes_na), ","))))
		    }
		}
		names(Path.list)[jh] <- as.character(path_names_na[jh])
		names(Path.list.bkup)[jh] <- as.character(path_names_na[jh])	  
	 	PhenoVec[[jh]] <- strsplit(vecCont[icont],"-")[[1]][2]
		names(PhenoVec)[jh] <- as.character(path_names_na[jh])
	      }
	      Path.list.all=c(Path.list.all,Path.list[1:5][])
		  Path.list.all.bkup=c(Path.list.all.bkup,Path.list.bkup)
		  PhenoVec <- PhenoVec[1:5][]
	      PhenoVecAll=c(PhenoVecAll,PhenoVec)
	  }#end path_names_na !=0



	  #path_names_pos
	  if(inputPath=="custom") {
	    path_names_pos <- MyInputFileList[[iij]][MyInputFileList[[iij]]$contrast == vecCont[icont] & MyInputFileList[[iij]]$updown == "up","pathways"]  
	  }else {
	    file0pos <- result[grep("gsea_report_for_na_pos_", result, perl=T)]
	    file1pos <- file0pos[grep(".xls", file0pos, perl=T)]
	    gsea_report_pos <- read.delim(file=file1pos, sep="\t")
	    path_names_pos  <- gsea_report_pos$NAME[1:40]
	  }
     
	  if(length(path_names_pos)!=0) {
	    Path.list <- list()
		Path.list.bkup <- list()
	    PhenoVec <- list()
	    for(jh in 1:length(path_names_pos)) {
		pathfile <- read.delim(paste(as.character(path_names_pos[jh]),
		  ".xls",sep=""),sep="\t")
		yesno <- which(as.character(pathfile[,8])=="Yes")
		Path.list[[jh]] <- as.character(pathfile[,2])[yesno]
		Path.list.bkup[[jh]] <- Path.list[[jh]]
		if(selectGene=="Top") {
		  
		  if((length(Path.list[[jh]])>5) && (length(path_names_pos) >=1)) {
		    Path.list[[jh]] <- Path.list[[jh]][1:5]
		  }
		}
		if(selectGene == "custom") {
			list_genes_pos <- MyInputFileList[[iij]][MyInputFileList[[iij]]$contrast 
				== vecCont[icont] & MyInputFileList[[iij]]$pathways == 
				as.character(path_names_pos[jh]),"genes"] 
		   if(length(unlist(strsplit(as.character(list_genes_pos),",")))!=0) {  
				Path.list[[jh]] <- unique(c(Path.list[[jh]],
				 unlist(strsplit(as.character(list_genes_pos),","))))
		   }	
		}
		names(Path.list)[jh] <- as.character(path_names_pos[jh])
		names(Path.list.bkup)[jh] <- as.character(path_names_pos[jh])
		PhenoVec[[jh]]	<- strsplit(vecCont[icont],"-")[[1]][1]
		names(PhenoVec)[jh] <- as.character(path_names_pos[jh])
	    }
		Path.list.all <- c(Path.list[1:5][],Path.list.all)
		Path.list.all.bkup=c(Path.list.all.bkup,Path.list.bkup)
		PhenoVec <- PhenoVec[1:5][]
		pathListName = paste("path.list.all.bkup",dirPath_list[iij], vecCont[icont],"Rdata",sep=".")
		if(archive) save(Path.list.all.bkup,file=file.path(dirPath,pathListName))
	    PhenoVecAll <- c(PhenoVec,PhenoVecAll)
	  } # end != 0 length(path_names_pos)

	  # Merge data from na and pos pathways into dataFUD data frame
	  #######################################################################################
	  NPhenoVecAll <- list()
	  for (jd in 1:length(unique(names(PhenoVecAll)))){
	    isa <- which(names(PhenoVecAll)==unique(names(PhenoVecAll))[jd])
		if(length(isa)==1){
		  NPhenoVecAll[[jd]] <- PhenoVecAll[[isa]]
		}else{
		  NPhenoVecAll[[jd]] <- unique(as.character(unlist(PhenoVecAll[isa])))
		}
	  }
	  names(NPhenoVecAll) <- unique(names(PhenoVecAll))
	  NPath.list.all <- list()
	  for (jd in 1:length(unique(names(Path.list.all)))) {
	    isa<- which(names(Path.list.all) == unique(names(Path.list.all))[jd])
	    if(length(isa)==1) {
	      NPath.list.all[[jd]]=Path.list.all[[isa]]
	      NPath.list.all[[jd]]=NPath.list.all[[jd]][!is.na(NPath.list.all[[jd]])]
	    }else{
	      NPath.list.all[[jd]]=unique(as.character(unlist(Path.list.all[isa])))
	      NPath.list.all[[jd]]=NPath.list.all[[jd]][!is.na(NPath.list.all[[jd]])]
	    }
	  }
	  names(NPath.list.all)=unique(names(Path.list.all))

	  len <- length(names(NPath.list.all))
	  all.unique <- unique(unlist(NPath.list.all))
	  Transpose <- matrix(0, ncol=len, nrow=length(all.unique))
	  colnames(Transpose) <- names(NPath.list.all)
	  rownames(Transpose) <- all.unique
	
	  for(ith in 1:len){
	    Transpose[as.character(NPath.list.all[[ith]]), ith] <- 1
	  }
	    sumCount <- apply(Transpose,1,sum)
	    dataFUD <- data.frame(dimnames(Transpose)[[1]],Transpose,sumCount)
	    names(dataFUD) <- c("Genes",dimnames(Transpose)[[2]] ,"Sum")

	    # Mapping Symbol to ProbeId and extraction of expression from eset object into MatExp matrix 
	    #####################################################################################
	    ListPheno <- c(strsplit(vecCont[icont],"-")[[1]][1], strsplit(vecCont[icont],"-")[[1]][2])
	    ## Load eset object and filtering
	    # eset

	    myList <- dataFUD[,1]
	    MyListMap <- list()
	    for (iss in 1:length(myList)) {
		MyListMap[[iss]] <- MyMapFile[toupper(MyMapFile$Gene.Symbol) == 
			as.character(dataFUD[iss,1]),"Probe.Set.ID"]  
	    }
	    names(MyListMap) <- as.character(myList)
	    xx <- MyListMap

	    #xx=as.list(MyIlluminaHumanRef8V3ALIAS2PROBE) 
	    myProbeList <- xx[as.character(myList)]   
	    names(xx) <- toupper(names(xx))
	    myProbeList <- xx[as.character(myList)] 

	    dataFUD <- dataFUD[!is.na(names(myProbeList)),]

	    myProbeList <- myProbeList[!is.na(names(myProbeList))]   
	    myProbeVec <- dimnames(exprs(eset))[[1]][rownames(exprs(eset))%in%as.character(unlist(myProbeList))]
	    expression <- exprs(eset)[myProbeVec,]
	    ##!!!

	    ListListPheno <- numeric()
	    for (id in 1:length(ListPheno)){
		ListListPheno <- c(ListListPheno, which(pData(eset)$Class==ListPheno[id]))
	    }
	    expression <- expression[,rownames(pData(eset))[ListListPheno]]
	    vec1 <- pData(eset)$Class[ListListPheno]


	    MatExp <- matrix(0,nrow=length(names(myProbeList)),ncol=ncol(expression), 
		byrow=T)
	    for (ij in 1:length(names(myProbeList))){
		sumsum <- sum(rownames(expression)%in%as.character(unlist(myProbeList[ij])))
		tem <- expression[rownames(expression)%in%as.character(unlist(myProbeList[ij])),]
		if(sumsum !=0){
			if(class(tem)=="matrix"){
				MatExp[ij,] <- apply(tem ,2,max)
			}else{
				MatExp[ij,] <- tem
			}
		} 
	    }
	    colnames(MatExp) <- vec1
	    rownames(MatExp) <- names(myProbeList)
	    if(length(which(apply(MatExp ,1,sum)==0.0)) >=1) {
		MatExp <- MatExp[-c(as.numeric(which(apply(MatExp ,1,sum)==0.0))),]
	    }

	    dataFUD <- dataFUD[dimnames(MatExp)[[1]],]
	    
	    ## Draw Visualization Plot
	    graphOut(MatExp,dataFUD,dirPath_list[iij],ListPheno,length(Path.list.all),file.path(dirPath,"visualization_plots"))
    }
  }

  setwd(dirData)
					
}


## TODO
# gsea.compare : compare topmost pathways found for the specified contrasts.
# 
gsea.compare <- function(eset, dirPath=paste(dirData,"gsea",sep="/"), 
				dirCommon=param$dirCommon, contrasts=getContrasts(),
				chip.version="illuminaV3.chip") {
	
}

### TODO
# getList : returns expression data from a list of symbols. It converts the 
#           symbols to their corresponding probes within the expression matrix.
#           This symbol list is usually coming from GSEA results.
# expression :  an expression matrix. eg. exprs(eset)
# symbols : a list of gene symbols to retrieve (string list)
#gsea.getList <- function(expression=exprs(eset), symbols="",
#					dirCommon=paste(param$dirCommon, "GeneSet"){
		
#}

## graphOut : Function for vizualisation of GSEA pathways results
##            Created by Ali Filali.
## TODO : cleanup + documentation for pipelinization 
graphOut <- function(MatExp,matdataF,db,pheno,npath,folder){

	if(!file.exists(folder)) {
		dir.create(folder)
	}
	
# Plot Layout
###############################################################################
png(paste(folder,"/Pathways","_",db,"_",paste(pheno[1],".VS.",pheno[2],sep=""),".png",sep=""),width=1800,height=1200,bg="white")
nf=layout(matrix(c(rep(0,npath),1,0,seq(2,npath+3,1),rep(npath+4,npath),npath+5,0),3,npath+2,byrow=TRUE),heights=c(0.3,6,0.3),widths=c(rep(0.03,npath),0.5,0.03))
layout.show(nf)         

#Samples Key
############################################################################
color_fun <- colorRampPalette(c("yellow", "red"))
MatBreaks=data.matrix(c(rep(0,length(which(dimnames(MatExp)[[2]]==pheno[1]))),rep(1,length(which(dimnames(MatExp)[[2]]==pheno[2])))))
breaks <- c(-1,0,1)
op=par(mar=c(0.1,2,0.2,10),lwd=2,cex=1.)
image(1:nrow(MatBreaks),1:ncol(MatBreaks),MatBreaks,col=color_fun(2),breaks=breaks,xlab="p",ylab="pa",tck=0,col.axis=NA,cex.lab=0.1)
#fo1=strwidth(pheno[1])
#fo2=strwidth(pheno[2])
#text(0.2+fo1,1,labels=pheno[1],cex=1.5)
#text(length(which(dimnames(MatExp)[[2]]==pheno[1]))+0.5+fo2,1,labels=pheno[2],cex=1.5)

#loop for pathways membership 
###########################################################################################################################
color_fun <- colorRampPalette(c("gray", "blue"))
vecCol=c(rep(color_fun(2)[1],5),rep(color_fun(2)[2],5))
matdataF=matdataF[,-c(1,dim(matdataF)[[2]])]
for ( i in 1:dim(matdataF)[2]){
MatBreaks=t(data.matrix(matdataF[,i]))
breaks <- c(-1,0,1)
op=par(mar=c(0.1,2,0.1,0.3),lwd=2,cex=1)
image(1:nrow(MatBreaks),1:ncol(MatBreaks),MatBreaks,col=color_fun(2),breaks=breaks,xlab="p",ylab="pa",tck=0,col.axis=NA,cex.lab=0.1)
mtext(names(matdataF)[i],side=2,at=3,line=0.15,adj=0,outer=FALSE,cex=2)
}


# Heatmap MatExp
##################################################################################################################
XY=MatExp
XY=XY[!is.na(dimnames(XY)[[1]]),]
fctz=function(x){(x-mean(x))/sd(x)}
XY=t(apply(XY,1,fctz))
XY=na.omit(XY)
mid=min(XY)
mad=max(XY)
op=par(mar=c(0.1,2,0.1,10),lwd=2,cex=1) #(botom,left,top,right)
image(1:ncol(XY),1:nrow(XY),t(XY),col = greenred(20),tck=0,col.axis=NA,cex.lab=0.1)
mtext(dimnames(XY)[[1]],side=4,at=1:nrow(XY),line=0.5,las=2,cex=2)


#color.key for MatExp
#####################################################################################################################
z <- seq(min(XY), max(XY), length = length(greenred(20)))

breaks=21
breaks <- seq(min(XY), max(XY),
            length = breaks)       
op=par(mar=c(2,0.5,2,3.5),lwd=2,cex=1)
zz=t(matrix(z, ncol = 1))
image(1:nrow(zz),1:ncol(zz),zz, col = greenred(20), breaks = breaks,xlab="",ylab="",tck=0,col.axis=NA,cex.lab=0.1)
mtext(round(z,digits=1),side=4,at=1:ncol(zz),line=0.5,las=2,cex=2)

#Pathway Key
#####################################################################################################################
color_fun <- colorRampPalette(c("gray", "white","blue","white"))
MatBreaks=data.matrix(c(0,1,2,3))
breaks <- c(-1,0,1,2,3)
op=par(mar=c(0.1,2,0.2,0.3),lwd=2,cex=1.)
image(1:nrow(MatBreaks),1:ncol(MatBreaks),MatBreaks,col=color_fun(4),breaks=breaks,xlab="p",ylab="pa",tck=0,col.axis=NA,cex.lab=0.1)
fo1=strwidth("Absent")
fo2=strwidth("Present")
text(1.5+fo1,1,labels=" Absent",cex=2)
text(3.5+fo2,1,labels=" Present",cex=2)

#color.key for Samples
#####################################################################################################################

color_fun <- colorRampPalette(c("yellow", "white","white","red","white","white"))
MatBreaks=data.matrix(c(0,1,2,3,4,5))
breaks <- c(-1,0,1,2,3,4,5)
op=par(mar=c(0.1,2,0.2,10),lwd=2,cex=1.)
image(1:nrow(MatBreaks),1:ncol(MatBreaks),MatBreaks,col=color_fun(6),breaks=breaks,xlab="p",ylab="pa",tck=0,col.axis=NA,cex.lab=0.1)
#axis(side=1,at=1:2,labels=c("Gene abscent","Gene present"))
fo1=strwidth(pheno[1])
fo2=strwidth(pheno[2])
text(1.5+fo1,1,labels=paste(" ",pheno[1],sep=""),cex=2)
text(4.5+fo2,1,labels=paste(" ",pheno[2],sep=""),cex=2)


dev.off()

}

# finalize.gsea.report : Versioning system for GSEA HTML reports.
#-----------------------------------------------------------------------------|
# Once finalized, the html report will be archived with a version tag.
# A file called .details will be hidden within visualization_plots folder with 
# the parameters used for the analysis, along with its description.
# The content will be zipped and named param$project_gsea_version.zip 
# The HTML folder will be copied to html_version as well.              
#-----------------------------------------------------------------------------|
# version : name to be appended to the report.
# author  : name of the analyst who generated it.
# desc    : analysis description to be written in the tracefile and the report.
# source  : Optional. where is the report located (default : 'html')
finalize.gsea.report <- function(version="full", author=param$author,
desc="NA", source=file.path(dirData,"gsea","gsea")) {
	if(file.exists(paste(source,version,sep="_")) || file.exists(
			paste(gsub("[^a-z0-9]",".",param$project, perl=T, 
			ignore.case=T), "_gsea_", version, ".zip", sep=""))) {
		print("Version ", version, " already exist! Try a different version name or delete the previous files.")
		return()
	}
	if(desc == "NA") {
		print("Please enter a detailed description for this GSEA report :")
		desc = scan(what="character", n=1, sep="^")
	}
	author = verify.author(author=author)
	tracefile = file(file.path(dirData,"tracefile.txt"),"a")
	cat(paste(format(Sys.time(), "%Y-%m-%d"),"\t",author,
			  "\tFinalized GSEA report : ", version, ". ",desc,'\n',sep=""), file=tracefile)
	close(tracefile)
	details = file(file.path(source,"visualization_plots",".details"),"w")
	cat("-----------| Analysis Parameters |-----------\n", file=details)
	for(i in 1:length(param)) {
		cat(paste(names(param)[i],"\t",param[i]),"\n",file=details)
	}
	cat("---------------------------------------------\nAuthor: ",
		file=details)
	cat(author, "\n",file=details)
	cat("Analysis version : ", version, file=details)
	cat("\nDescription : ", desc, "\n", file =details)
	cat("\nFinalized on : ", format(Sys.time(), "%Y-%m-%d %X"), 
		", using illupipe for GSEA version ", illupipe.gsea.version, "\n", file =details)
	close(details)
# System Commands -- for UNIX systems only. 
	system(paste("cp -r ", source, " ", source, "_", version, sep=""))
	system(paste("zip -9 -r ", gsub("[^a-z0-9]",".",param$project, perl=T, 
									ignore.case=T), "_gsea_", version, ".zip ", source, sep=""))
}


# gsea.heatmap : Generate heatmaps automatically for GSEA analyses
#                showing all the contrasts within a subfolder. Scaling is automatic for pdfs.
#                It will concatenate each pathways from all analysis in a subfolder, and will
#                aggregate to the heatmap only those with a value passing a certain treshold.
#                Figures are scaled to accomodate for various pathway names length or various
#                number of pathways.
# input.dir : Where to find the subdirectories (dbList)
# dbList : list of subdirectories to scan. eg. c("C2","C3","C5")
# value    : what column to use for filtering. eg. 'NOM.p.value'
# treshold : What minimum value to use for filtering. eg. 0.05
# top.n : <-TODO-> Keep only the top N pathways by report (up or down). eg. 50.
# doPairwise : Logical. do you want pairwise correlation scatterplot mtarices to be generated?
# scale : character ('none', 'row' or 'col'), argument to passe to pheatmap for scaling
#--------------------------------------------------------------------------------------------
# side-effect : Generates pdf files inside current working directory.
gsea.heatmap <- function(input.dir='gsea', dbList=c('C1','C2','C3','C5'),
                         value='FDR.q.val', treshold=0.25, top.n=NULL,
                         keepEmptyContrasts=FALSE, displayScore='NES',
                         doPairwise=TRUE, returnDF = FALSE, scale='none') {
  require(pheatmap)
  require(ggplot2)
  if(returnDF) {output = list()}
  # Loop through all DBs
  for(cx in dbList) {
    contrasts <- dir(file.path(input.dir,cx))
    pathNames <- c()
    results <- list()
    for(contrast in contrasts) {
      pathways  <- dir(file.path(input.dir,cx,contrast))
      #print(pathways)
      pathways <- pathways[grepl('gsea_report_for_na_',pathways)]
      pathways <- pathways[grepl('.xls',pathways)]
      pathDW <- read.delim(file.path(input.dir,cx,contrast,pathways[1]), row.names=1)
      pathUP <- read.delim(file.path(input.dir,cx,contrast,pathways[2]), row.names=1)
      pathUPfiltered <- subset(pathUP, pathUP[[value]] < treshold)
      pathDWfiltered <- subset(pathDW, pathDW[[value]] < treshold)
      ## TODO : subset pathways for top.n
      temp <- unique(union(rownames(pathDWfiltered),rownames(pathUPfiltered)))
      pathNames <- unique(union(temp,pathNames))
      results[[sub('.GseaPreranked..............','',contrast)]] <-
        rbind(pathUP[which(!is.na(pathUP[,displayScore])),3:8],
              pathDW[which(!is.na(pathDW[,displayScore])),3:8])
      
    }
    
    temp <- list()
    i <- 1
    contrasts <- names(results)
    for(pathways in results) {
      if(nrow(pathways) > 0) {
        temp[[names(results)[[i]]]] <- getES(x=pathways,names=pathNames, typeScore=displayScore)
      }
      else if(keepEmptyContrasts){
        temp[[names(results)[[i]]]] <- rep(0,length(pathNames))
      }
      else { contrasts <- contrasts[-i]}
      i <- i + 1
    }
    # test if heatmap is possible (enough rows or col)
    if(length(temp) >= 2  && length(pathNames) > 2) {
      df <- data.frame(matrix(unlist(temp), nrow=length(pathNames)))
      rownames(df) <- pathNames
      colnames(df) <- contrasts
      #colnames(df) <- c('x1','x2','x3','x4','x5','x1b','x2b','x3b','x4b','x5b')
      #print(head(df))
      #print(summary(df))
      #print(dim(df))
      #df <- apply(df,c(1,2),as.numeric)

      names(df) <- colnames(df)
       ## REVERSE some contrasts:
      #df[,!grepl('LTB', contrasts)] <-
      # apply(df[,!grepl('LTB', contrasts)], c(1,2),
      #       function(x) {return(x * -1)})
     
     if(returnDF) {output[[cx]] <- df }
     # Pairwise ScatterPlot showing correlations
      if(doPairwise) {
        pdf(file=file.path(input.dir,
              paste(cx,value,treshold,displayScore,'pairwise','pdf',sep='.')))
        print(plotmatrix(df) + geom_smooth(colour='red'))
      #pairs(df)
        print(ggcorplot(data=df, var_text_size = 4,
                        cor_text_limits = c(5,10)))
        dev.off()
      }
      #df[,-5] <- apply(df[,-5], c(1,2), function(x) {return(x * -1)})
      #colnames(df) <- c('HP - C', 'LP - C', 'P - C', 'F (P - C)', 'HP - LP', 'M (P - C)')
      # Split C2 database
      if(cx == "C2") {
        temp <- df[grepl('REACTOME',rownames(df)),]
        rownames(temp) <- sub('REACTOME_', '', rownames(temp))
        if(ncol(temp) >= 2 && nrow(temp) >= 2) {
          pheatmap(temp,filename=file.path(input.dir,
                  paste(cx,'REACTOME',value,treshold,displayScore,'pdf',sep='.')),
                 width=max(c(7,ceiling(max(nchar(rownames(temp)))/6))),
                 height = max(c(7,ceiling(nrow(temp)/7))), scale=scale)
        }
        
        temp <- df[grepl('KEGG',rownames(df)),]
        if(ncol(temp) >= 2 && nrow(temp) >= 2) {
          rownames(temp) <- sub('KEGG_','',rownames(temp))
          pheatmap(temp,filename=file.path(input.dir,
                        paste(cx,'KEGG',value,treshold,displayScore,'pdf',sep='.')),
                 width=max(c(7,ceiling(max(nchar(rownames(temp)))/6))),
                 height = max(c(7,ceiling(nrow(temp)/7))), scale=scale)
        }
        temp <- df[grepl('BIOCARTA',rownames(df)),]
        if(ncol(temp) >= 2 && nrow(temp) >= 2) {
          rownames(temp) <- sub('BIOCARTA_','',rownames(temp))
          pheatmap(temp,filename=file.path(input.dir,
                  paste(cx,'BIOCARTA',value,treshold,displayScore,'pdf',sep='.')),
               width=max(c(7,ceiling(max(nchar(rownames(temp)))/6))),
               height = max(c(7,ceiling(nrow(temp)/7))), scale=scale)
        }
        df <- df[!grepl('BIOCARTA',rownames(df)),]
        df <- df[!grepl('KEGG',rownames(df)),]
        df <- df[!grepl('REACTOME',rownames(df)),]
      }
      if(ncol(df) >= 2 && nrow(df) >= 2) {
        pheatmap(df,filename=file.path(input.dir,
                  paste(cx,'pathways',value,treshold,displayScore,'pdf',sep='.')),
                 width=max(c(7,ceiling(max(nchar(rownames(df)))/6))),
                 height = max(c(7,ceiling(nrow(df)/7))), scale=scale)
      }
    }
  }
  if(returnDF) { return(output) }
}




# getES : Will return a list of scores, using GSEA outputs and a list of pathway
#         names from which to retrieve that score.
#        x, y : data.frame containing GSEA results for a given contrast.
#               Usually, x will represent the up-regulated pathways and
#               y will represent the down-regulated pathways. x or y can
#               be either one. You can also put down and up regulated pathways
#               in one single data.frame, i.e. <x>.
#       names : list of pathway names to use for output.
#       typeScore : The name of the GSEA output column you want to retrieve.
#----------------------------------------------------------------
# output : a vector of length(<names>) containing <typeScore>
getES <- function(x, y=NULL, names=c1_names, typeScore = "NES") {
  output <- rep(0,length(names))
  names(output) <- names
  dat <- x
  
  for(i in 1:nrow(dat)) {
    if(!is.na(output[rownames(dat)[i]])) {
      if(!is.na(dat[i, typeScore])) {
        if(abs(output[rownames(dat)[i]]) < abs(dat[i, typeScore])) {
          output[rownames(dat)[i]] <- dat[i, typeScore]
        }
      }
    }
  }
  if(!is.null(y)) {
    dat <- y ## Negative Enrichments
    for(i in 1:nrow(dat)) {
      if(!is.na(output[rownames(dat)[i]])) {
        if(!is.na(dat[i, typeScore])) {
          if(abs(output[rownames(dat)[i]]) < abs(dat[i, typeScore])) {
            output[rownames(dat)[i]] <- dat[i, typeScore]
          }
        }
      }
    }
  }
  return(output)
}



# variable for versioning system
illupipe.gsea.version = "0.55"