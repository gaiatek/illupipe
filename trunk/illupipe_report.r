###############################################################################
# illupipe_report.R  v 3.1 Beta                                               #
#-----------------------------------------------------------------------------#
# Illumina MicroArray Pipeline Toolkit - Reporting functions                  #
# Most of the pipeline functions (illupipe) offer the option to generate HTML #
# reports of what has been done during the dat analysis.                      #
# This file contains the HTML-related portions of the pipeline.               #
#-----------------------------------------------------------------------------#
# @author : Jean-Philippe Goulet                                              #
# @date : January 13th 2011                                                   #
###############################################################################



# this function ensures directories and css exist
check.dirs <-function(report.dir=dirReport,data.dir=dirData)
{
	# Prepare where the data will be written
	if(!file.exists(file.path(report.dir))) 
	{
		dir.create(file.path(report.dir))
	}
	files.dir = file.path(report.dir,"files")
	if(!file.exists(files.dir))
	{
			dir.create(files.dir)
	}
	css.dir = file.path(files.dir,'css')
	if(!file.exists(css.dir)) 
	{
		dir.create(css.dir)
		file.copy( file.path(data.dir,'report_parts','report.css'), file.path(css.dir)  )
		file.copy( file.path(data.dir,'report_parts','R2HTML.css'), file.path(css.dir)  )
		file.copy( file.path(data.dir,'report_parts','hwriter.css'), file.path(css.dir)  )
	}
}
# Deal with the css file

# This function fetches html bits need for report contruction
getHtmlBits <-function(data.dir=dirData)
{	
	dir  = file.path(data.dir,'report_parts')  
	html = list(
	"head" = readLines(file.path(dir,'head.html'))
	,"foot" = readLines(file.path(dir,'foot.html'))
	,"annotations" =  readLines(file.path(dir,'annotations.html'))
        ,"removals" =  readLines(file.path(dir,'removals.html'))
	,"preprocessing" =  readLines(file.path(dir,'preprocessing.html'))
	,"exploratory" =  readLines(file.path(dir,'exploratory.html'))
	,"dgea_model_head" =  readLines(file.path(dir,'dgea_model_head.html'))
	,"dgea_model_tail" =  readLines(file.path(dir,'dgea_model_tail.html'))
	)
	return(html)	
}




#project.name="Test" 
#report.dir = dirReport;data.dir=dirData;third.id="DonorID"
#annots.items =NULL;#  list("raw"=TRUE,"samples"=TRUE,"arrays"=TRUE,"probes"=FALSE,"gene.sets"=TRUE)
#removal.items=NULL;#list("removed.arrays"=TRUE,"missing.samples"=TRUE)
#preprocessing.items = FALSE
#exploratory.items=FALSE
#dgea.items=list("Stim, Vaccine, Time"=c("Stim:Vaccine_with_GSEA") ) # A list whose names are model names, and elements are contrast sets names. As in the fits object
#dgea.items=list("Stim, Vaccine, Time"=c("Vaccine"),"Stim"=c("Stim") )
#dgea.items = "all"


create.report <-function(project.name="Unspecified Project Name",report.author="Unspecified author",third.id=NULL
			,report.dir = dirReport, data.dir=dirData,
			annots.items = list("raw"=TRUE,"samples"=TRUE,"arrays"=TRUE,"probes"=TRUE,"gene.sets"=TRUE),
			removal.items=list("removed.arrays"=TRUE,"missing.samples"=TRUE),
			preprocessing.items = TRUE, exploratory.items=TRUE, dgea.items="all",
			dgea.params = list(n.genes.toplists=1000,n.genes.heatmaps=50,heatmaps.alt.text=NULL
						,p.thresh=0.05,q.thresh=0.05,heatmap.bl.subtract=FALSE)
)
{
	check.dirs(report.dir,data.dir)
	
    # Check that basic requirements are fulfiled before running the function

	if(!file.exists(file.path(data.dir,'eset_raw.Rdata'))) { stop("Missing eset_raw.Rdata in your data directory.")}
	if(!file.exists(file.path(data.dir,'eset.Rdata'))) { stop("Missing eset.Rdata in your data directory.")}
	if(!file.exists(file.path(data.dir,'probes.txt'))) { stop("Missing probes.txt in your data directory.")}
	if(!file.exists(file.path(data.dir,'contrasts_aliases.txt'))) { stop("Missing contrasts_aliases.txt in your data directory.")}
	if(!file.exists(file.path(data.dir,'worksheet.xls'))) { stop("Missing worksheet.xls in your data directory.")}
    if(!file.exists(file.path(data.dir,'report_parts'))) { stop("Missing reports_parts templat folder in your data directory.")}
	if(!is.null(annots.items$gene.sets)) {
		if(annots.items$gene.sets&!file.exists(file.path(data.dir,'report_parts'))) { stop("Missing geneSet.txt file in you data directory and you have set gene.sets to TRUE.")}
	}
	# Report Heading
	html = getHtmlBits(data.dir)$head
	html = gsub('TITLE',"Microarray Data Analysis",html) # Replace title
	html = gsub('PROJECT_NAME',project.name,html) # Replace project name
	html = gsub('TIMESTAMP',Sys.time(),html) # timestamp

	# The Annotations section
	if(!is.null(annots.items))
	{
		print("Generating Annotations section...")
    		html = c(  html,report.annotations(report.dir,data.dir,annots.items,third.id=third.id) )
	}
	# The arrays removal section
	if(!is.null(removal.items))
	{
		print("Generating removal section...")
    		html = c(  html, report.removals(report.dir,data.dir,removal.items,third.id=third.id) )
	}
	# The pre-processing section
	if(preprocessing.items)
	{
		print("Generating preprocessing section...")
    		html = c(  html, report.preprocessing(report.dir,data.dir) )
	}
	# Exploratory Analysis
	if(exploratory.items)
	{
		print("Generating exploratory section...")
    		html = c(  html, report.exploratory(report.dir,data.dir) )
	}
	# DGEA
	if(!is.null(dgea.items))
	{
		print("Generating DGEA section...")
    		html = c(  html, report.dgea(report.dir,data.dir,dgea.items,dgea.params,third.id) )
	}

	# ok.. so the numberof chips ncol() is befoere removing technical replicates. So technically you should read the sample
	# annotation. Also you should mention what is done with technical replicates...

	# Finally, add the footer
	html = c(html,getHtmlBits(data.dir)$foot)
	html = c(html,paste("<br><i>Analysis performed by ",report.author,", on ", format(Sys.time(), "%Y-%m-%d %X"),".</i>"))

	writeLines(html,file.path(report.dir,"index.html"))

}


# report.dir = dirReport; data.dir=dirData
report.dgea <- function(report.dir, data.dir, dgea.items, dgea.params, third.id) # Freedom to specifiy the heatmap labels...
{
	files.dir = file.path(report.dir,'files')	

	# Load and subset the fits object to the desired models and contrast sets
	fits = load.fits(data.dir)
	if(is.list(dgea.items))
	{
		fits = fits[names(dgea.items)] # subset to linear model
		if(  any( ! names(dgea.items) %in% names(fits) )){stop("No model of that name. check the names of dgea.items")}
		fits = sapply(names(fits),function(model.name){
			f = fits[[model.name]]
			if(  any( ! dgea.items[[model.name]] %in% names(f$fits2) )){stop("No set of contrasts of that name. Check your contrast sets names of dgea.items")}
			f$fits2 = f$fits2[ dgea.items[[model.name]] ]
			return(f)
			},simplify=FALSE)# subseyt contrast sets
	}

	# Load the eset
	x = load.eset(data.dir,third.id=third.id)

	html = c()
	for(model.name in names(fits))
	{
		print("Beginning one model...")
		f = fits[[model.name]] # subset the fits object. This object has 2 slots: fit and fits2
		f.dir = file.path(files.dir,model.name) # create directories
		dir.create(f.dir, showWarnings = FALSE)

		# Create replication xls and html files
		getReplication(lm=f$fit,output.dir=f.dir, lm.name=model.name,css.link='../css/hwriter.css')
		
		# Get the html and set values
		html = c(html,getHtmlBits(data.dir)$dgea_model_head)# write table head
		html = gsub('MODEL_DIR',file.path("files",model.name),html,fixed=TRUE)
		if(is.null(f$fit$dupcor.block))
		{
			html = gsub('VARIABLE(S)_OF_INTEREST',paste(colnames(f$fit$variables),collapse=", "),html,fixed=TRUE)
		}else{
			html = gsub('VARIABLE(S)_OF_INTEREST',
			paste(paste(colnames(f$fit$variables),collapse=", "),"with",f$fit$dupcor.block,"as a blocked-random effect.")
			,html,fixed=TRUE)
		}

		for(contset.name in names(f$fits2))
		{
			print("Beginning one set of contrasts...")
			contset.dir = file.path(f.dir,contset.name)
			dir.create(contset.dir, showWarnings = FALSE)
			html = c(html, report.contrast.set(data.dir,report.dir,f$fits2[[contset.name]],f$fit,x
								,model.name,contset.name,files.dir, dgea.params))
		}
		html = c(html,getHtmlBits(data.dir)$dgea_model_tail)# write table head# write table tail	
	}	
	return(html)
	# Ok. Begin by writing the data and generating the files
	# gsea or not.	
}

# lm  = fits[[1]]$fit; output.dir="";prefix="lm";lm.name="lm";
# This function takes a fit object and writes the replication to xls and html. The Fit object HAS to have slot "variables".
# lm : MArrayLM object (fit)
# output.dir : where the xls and html should be written
# css.link : the relative location of the css file (hwriter) that should be written in the html file.
# ret : should the data.frame be returned
# value:  writes to output.dir/replication.html/.xls
getReplication <- function(lm,output.dir=".", lm.name="lm",css.link="hwriter.css",ret=FALSE)
{

		# Check variables slot
		if(is.null(lm$variables))
		{stop("Current implementation of getReplication requires lm object to have non-standard slot \"variables\" ")}

		# Create replication xls and html files
		Replication           = as.data.frame(melt(table(lm$variables))) #  data.frame("Replication"=apply(f$fit$design,sum,MARGIN=2)) 
		Replication           = Replication[order(Replication$value,decreasing=TRUE),];
		rownames(Replication) = c(1:nrow(Replication))
		if(ncol(Replication)==2)
		{
			colnames(Replication)[1]=colnames(lm$variables)
		}
		colnames(Replication)[ncol(Replication)] = "Replication"
		WriteXLS("Replication", file.path(output.dir,"replication.xls")
					, row.names =TRUE, AdjWidth = TRUE,BoldHeaderRow = TRUE)
		p=openPage(file.path(output.dir,"replication.html"),link.css=css.link,title=paste("Replication for model",lm.name))
		hwrite(paste("Replication for model",lm.name), p, heading=3,center=TRUE)
		hwrite(Replication,p,row.names=TRUE,row.bgcolor='#ddaaff',row.style=list('font-weight:bold'),br=TRUE)
		closePage(p)
		if(ret){return(Replication)}else{return(NULL)}
}

## lm  = fits[[1]]$fit; output.dir="";lm.name="lm";
# This function writes the content of a fit object to a file with column names formatted for the macro: P. Padj. FC.
# lm : MArrayLM object (fit2)
# output.dir : where the xls should be written
getAnalysisOutput <-function(lm,output.dir=".",ret=FALSE,cont.aliases = colnames(lm))
{	
	pvals    = lm$p.value                            
	colnames(pvals) = paste("P.",cont.aliases,sep="")
	adjpvals = apply(pvals,MARGIN=2,multtest.BH)
	colnames(adjpvals) = paste("Padj.",cont.aliases,sep="")
	fcs      = apply(lm$coefficients,MARGIN=2,toFold)
	colnames(fcs)   = paste("FC.",cont.aliases,sep="")
	analysis.output = data.frame(pvals,adjpvals,fcs,lm$genes,check.names=FALSE,stringsAsFactors=FALSE)
	WriteXLS("analysis.output",file.path(output.dir,"analysis_output.xls")
			,row.names =TRUE, AdjWidth = TRUE,BoldHeaderRow = TRUE)
	if(ret){return(analysis.output)}else{return(NULL)}	
}


# This funcion accepts an MArrayLM object, writes (xls,html) and returns well formatted topTables
# lm  = fits[[1]]$fits2[[1]]; output.dir=".";lm.name="lm";
# lm : A fit2 object
# lm$genes has to have controlled column names SYMBOL and PRobeID
getTopTables <-function(lm,cont.aliases = colnames(lm),output.dir=".",html.n.max=1000,ret=TRUE)
{
	tts   = sapply(colnames(lm),function(cont){  topTable(lm,coef=cont,number=nrow(lm), adjust.method="BH")  },simplify=FALSE)
	tts$F = topTableF(lm,number=nrow(lm), adjust.method="BH")
	tts$F = tts$F[,-1*(1+ncol(lm$genes)):(ncol(lm$genes)+ncol(lm))] # remove annoyinh column in topTableF
	tts   = lapply(tts,function(top){ # This part formats topTable by renaming, reordering and adding columns
			top = top[,c("ProbeID","SYMBOL",setdiff(colnames(top),colnames(lm$genes))
					,setdiff(colnames(lm$genes),c("ProbeID","SYMBOL")))]	
			if("logFC" %in% colnames(top)) 
			{
				top$logFC=toFold(top$logFC);colnames(top) = gsub('logFC','Fold Change',colnames(top))
				top$AveExpr=abs(top[["Fold Change"]]);colnames(top) = gsub('AveExpr','|Fold Change|',colnames(top))
			}
			top$AveExpr=NULL;top$B=NULL			
			return(top)
	})#reorder and remove some coluns
	names(tts) = c(cont.aliases,"F") # Replace by their alias
	for(ct in names(tts)){ # Write to excel AND sortable table
			top = tts[[ct]];
			WriteXLS("top",file.path(output.dir,paste("top_",ct,".xls",sep=""))
					,row.names =FALSE, AdjWidth = TRUE,BoldHeaderRow = TRUE)
			# Round a few columns before outputting to HTML
			#Fold Change	|Fold Change|	t	P.Value	adj.P.Val
			sortable.html.table(top[1:min(html.n.max,nrow(top)),], paste("top_",ct,".html",sep="")
						, output.dir, ct)
		}
	
	if(ret){return(tts)}else{return(NULL)}
}

# This function counts the number of DEGs from an MAArrayLM object. It plots a barplot and writes its results to excel.
# lm  = fits[[1]]$fits2[[1]]; output.dir=".";p=0.05;q=0.05;cont.aliases =colnames(lm)
getCounts <- function(lm,cont.aliases = colnames(lm),output.dir=".",p=0.05,q=0.05,ret=TRUE)
{
	# Count 
	dt.raw      = decideTests(lm, p.value = p, adjust.method="none")
	colnames(dt.raw) = cont.aliases
	dt.adj      = decideTests(lm, p.value = q, adjust.method="BH")
	colnames(dt.adj) = cont.aliases
	topF        = topTableF(lm,adjust.method="BH",number=nrow(lm))  
	n.degsF.raw = sum(topF$P.Value   < p)
	n.degsF.adj = sum(topF$adj.P.Val < q)

	cnts  = list(
	"Raw p-value" = list(
		"BOTH"=data.frame("Contrast" = c( colnames(dt.raw)          ,"F"        )
			  ,"Number of Genes" = c( apply(dt.raw,MARGIN=2,FUN=function(co){sum(co != 0)}) ,n.degsF.raw) 
 			,check.names=FALSE,stringsAsFactors=FALSE)
		,"UP"= data.frame("Contrast" = colnames(dt.raw) 
			         ,"Number of Genes" = apply(dt.raw,MARGIN=2,FUN=function(co){sum(co == 1)})
			,check.names=FALSE,stringsAsFactors=FALSE)
		,"DOWN" = data.frame("Contrast"        = colnames(dt.raw) 
			            ,"Number of Genes" = apply(dt.raw,MARGIN=2,FUN=function(co){sum(co== -1)})
			,check.names=FALSE,stringsAsFactors=FALSE)
	)
	,"Adj. p-value"= list(
		"BOTH"= data.frame("Contrast"=c(colnames(dt.adj) ,"F")
			,"Number of Genes"=c( apply(dt.adj,MARGIN=2,FUN=function(co){sum(co!=  0)}) ,n.degsF.adj) 
			,check.names=FALSE,stringsAsFactors=FALSE)
		,"UP"= data.frame("Contrast"= colnames(dt.adj)
			,"Number of Genes"=apply(dt.adj,MARGIN=2,FUN=function(co){sum(co==  1)})
			,check.names=FALSE,stringsAsFactors=FALSE)
		,"DOWN"= data.frame("Contrast"= colnames(dt.adj)
			,"Number of Genes"=apply(dt.adj,MARGIN=2,FUN=function(co){sum(co== -1)}) 
			,check.names=FALSE,stringsAsFactors=FALSE)
	)
	)
	

		
	# Melt and format for ggplot2
	counts.m = melt(cnts,measure.vars="Number of Genes",id.vars=c("Contrast"))
	colnames(counts.m)[3:5] = c("Number of Genes","Direction","Type")	
	counts.m = counts.m[order(counts.m[["Type"]],counts.m[["Direction"]],counts.m[["Number of Genes"]],decreasing=TRUE),]
	#counts.m = counts.m[counts.m$Direction!="BOTH" & !is.na(counts.m[["Number of Genes"]]),]# This removes Both and F-test from barplot
	counts.m$Direction = factor(counts.m$Direction,ordered=TRUE,levels=c("UP","DOWN","BOTH"))
	#counts.m[["Number of Genes"]][counts.m$Direction=="DOWN"] = -1 * counts.m[["Number of Genes"]][counts.m$Direction=="DOWN"]
	#counts.m$Contrast =  factor(counts.m$Contrast,ordered=TRUE,levels=unique(counts$Contrast)) 	# This reverses down
	
	# Plot and write
	png(file=file.path(output.dir,"ndegs_barplot.png"),height=600,width=(300+300*ncol(lm)/12 ) ) 
	print(qplot(x=Contrast,y=`Number of Genes`,data=counts.m,fill=Direction,geom="bar",position="dodge",facets=Direction~Type)+
		opts(axis.text.x = theme_text(angle = 90,hjust=1)) )
	dev.off()
	WriteXLS("counts.m",file.path(output.dir,"ndegs_barplot.xls"),row.names =FALSE, AdjWidth = TRUE,BoldHeaderRow = TRUE)

	if(ret){return(cnts)}else{return(NULL)}	
}


#report.dir = dirReport; data.dir=dirData; fits=load.fits(data.dir);fit2 = fits[[1]]$fits2[[4]] ;contset.name=names(fits[[1]]$fits2)[4];model.name=names(fits)[1];files.dir = "Report/files";report.dir = dirReport; data.dir=dirData;dgea.params=list(n.genes.toplists=1000,n.genes.heatmaps=50,n.genesets.toplists=100,p.thresh=0.05,q.thresh=0.05);
report.contrast.set <- function(data.dir, report.dir ,fit2, fit,x, model.name, contset.name, files.dir, dgea.params)
{

	####### analysis output file ########################
	print("Writing analysis output...")
	getAnalysisOutput(lm=fit2,output.dir=file.path(files.dir,model.name,contset.name),cont.aliases = contrasts2aliases(colnames(fit2),data.dir))
	
	############ TopTables (Gene lists) ###################
	print("TopTables...")
	tops = 	getTopTables(lm=fit2,cont.aliases = contrasts2aliases(colnames(fit2),data.dir) 
				,output.dir=file.path(files.dir,model.name,contset.name)
				,html.n.max=dgea.params$n.genes.toplists,ret=TRUE)

	############# Number of degs ########################
	print("Number of DEGS...")
	counts = getCounts(lm=fit2,cont.aliases = contrasts2aliases(colnames(fit2),data.dir)
				,output.dir=file.path(files.dir,model.name,contset.name)
				, p=dgea.params$p.thresh, q=dgea.params$q.thresh,ret=TRUE)

	#################### heatmaps ###############################
	print("heatmaps...")
	alt.text = dgea.params$heatmaps.alt.text
	if(is.null(alt.text)){alt.text      = colnames(fit2$variables)}
	if(dgea.params$heatmap.bl.subtract){lm.orig=fit}else{lm.orig=NULL}
	getFitHeatmaps( 
			 x             = x
			,lm          = fit2
		        ,output.dir    = file.path(files.dir,model.name,contset.name)
			,cont.aliases  = contrasts2aliases(colnames(fit2),data.dir)
			,alt.text      = alt.text
			,ngenes        = min(dgea.params$n.genes.heatmaps,nrow(x))   
			,lm.orig       = lm.orig
	)
	#############################################################	
	## GSEA ##	
	if(!is.null(fit2$gsea))
	{
		print("GSEA...")
		gene.sets.db  = getGeneSets(data.dir)
		for(cont in colnames(fit2))
		{
			scores = fit2$gsea[[cont]]$P.Value[,"Mixed",drop=FALSE];colnames(scores) = "P.Value"
			scores = merge(scores,gene.sets.db,all.x=TRUE,all.y=FALSE,by.x=0,by.y="GeneSetID")
			scores = scores[order(scores$P.Value),]
			colnames(scores)[1] = "GeneSetID"
			WriteXLS("scores",file.path(files.dir,model.name,contset.name
							,paste("gsea_",contrasts2aliases(cont,data.dir),".xls")  )
							,row.names =FALSE, AdjWidth = TRUE,BoldHeaderRow = TRUE)
			sortable.html.table(scores[1:min(dgea.params$n.genesets.toplists,nrow(scores)),]
					, paste("gsea_",contrasts2aliases(cont,data.dir),".html",sep="")
					, file.path(files.dir,model.name,contset.name), contrasts2aliases(cont,data.dir))
			
		}
		if(!is.null(fit2$gseaF))
		{
		# Now for the F-test
		scores = as.data.frame(fit2$gseaF);colnames(scores) = "P.Value"
		scores = merge(scores,gene.sets.db,all.x=TRUE,all.y=FALSE,by.x=0,by.y="GeneSetID")
		scores = scores[order(scores$P.Value),]
		colnames(scores)[1] = "GeneSetID"
		WriteXLS("scores",file.path(files.dir,model.name,contset.name
							,paste("gsea_F",".xls")  )
							,row.names =FALSE, AdjWidth = TRUE,BoldHeaderRow = TRUE)
	sortable.html.table(scores[1:min(dgea.params$n.genesets.toplists,nrow(scores)),]
					, paste("gsea_F",".html",sep="")
					,file.path(files.dir,model.name,contset.name), "F")
		}
	}
	
	#############################################################	

	print("Writing HTML...")
	## The actual HTML table : order is presumed to be preserved
	df = list()
	df[[paste("Contrasts of Set<br><b>",contset.name,"</b>",sep="")]] = names(tops)
	df[[paste("Number of DEGs<br>p-value<",dgea.params$p.thresh,"<br>barplot (",hwrite('png',link=paste("files",model.name,contset.name,"ndegs_barplot.png",sep="/")),", ",hwrite('xls',link=paste("files",model.name,contset.name,"ndegs_barplot.xls",sep="/")),")",sep="")]] = counts[["Raw p-value"]]$BOTH[["Number of Genes"]]
	df[[paste("Number of DEGs<br>FDR<",100*dgea.params$q.thresh,"%<br>barplot (",hwrite('png',link=paste("files",model.name,contset.name,"ndegs_barplot.png",sep="/")),", ",hwrite('xls',link=paste("files",model.name,contset.name,"ndegs_barplot.xls",sep="/")),")",sep="")]] = counts[["Adj. p-value"]]$BOTH[["Number of Genes"]]
	
	df[["Gene Lists"]] = sapply(names(tops),function(ctn)
	{
	 paste("(",
	 hwrite('html',link=paste("files",model.name,contset.name,paste("top_",ctn,".html",sep=""),sep="/"))," ,"
	,hwrite('xls',link=paste("files",model.name,contset.name,paste("top_",ctn,".xls",sep=""),sep="/"))
	,")",sep="")
	},simplify=TRUE)

	df[["Heatmaps"]] = sapply(names(tops),function(ctn)
	{
	 paste("(",
	 hwrite('png',link=paste("files",model.name,contset.name,paste("heatmap_",ctn,".png",sep=""),sep="/"))," ,"
	,hwrite('xls',link=paste("files",model.name,contset.name,paste("heatmap_",ctn,".xls",sep=""),sep="/"))
	,")",sep="")
	},simplify=TRUE)
	df[["Heatmaps"]][[length(df[["Heatmaps"]])]] = paste("Fold-Changes ",df[["Heatmaps"]][[length(df[["Heatmaps"]])]],"<br>Expressions ",
	"(",
	 hwrite('png',link=paste("files",model.name,contset.name,"heatmap_F_expressions.png",sep="/"))," ,"
	,hwrite('xls',link=paste("files",model.name,contset.name,"heatmap_F_expressions.xls",sep="/"))
	,")"
	,sep="")

	if(!is.null(fit2$gsea))
	{
		df[["ROAST GSEA"]] = sapply(names(tops),function(ctn)
		{
	 	paste("(",
	 	hwrite('html',link=paste("files",model.name,contset.name,paste("gsea_",ctn,".html",sep=""),sep="/"))," ,"
		,hwrite('xls',link=paste("files",model.name,contset.name,paste("gsea_",ctn,".xls",sep=""),sep="/"))
		,")",sep="")
		},simplify=TRUE)
	}

	print(df) ## REMOVE
	df = do.call(cbind,df)



	
	#p = openPage("Report/test.html")
	#hwrite(df,p)
	#closePage(p)


	# Contrast names
	#cont.names = c( contrasts2aliases(colnames(fit2),data.dir), "F-test")



	# gsea of not
	# maybe good idea to reuse JP topTables
	# Sortable tables
	# (Click for the whole formula)
	# Symbols or top
	# file with everything

	html = hwrite(df,row.names=FALSE,br=TRUE,center=TRUE,row.bgcolor=c('#ffffaa',rep('#f1ecff',times=ncol(fit2)),'#ffbbaa'))


	html = c(html,paste("All contrasts at once (",hwrite("xls"
				,link=paste("files",model.name,contset.name,"analysis_output.xls",sep="/")),")<br><br><br><br><br>",sep=""))

	return(html)
}









# This function creates heatmaps for a fit object : All contrasts and F-test
# TODO : would really be useful to delegate heatmaps.2 job to getheatmap... JP is right.
# the logic shoulb be in calling function.. would allow for instance naming F after contset.
# for each contrasts, the name of the intial coefficient
getFitHeatmaps <- function(x,lm,output.dir,cont.aliases,alt.text=c("SampleID"),ngenes=50,lm.orig=NULL,lfc=0,p.value=1,adjust.method="BH",scale="row")
	{

	cont.matrix = lm$cont.matrix
		
	     for(cont in colnames(lm))
	     {
		# Deal with the baseline subtraction
		e = x
		scale = scale
		symbreaks=TRUE
		if( !is.null(lm.orig)  ) #will bl subtract only if contrast of form A-B
		{
			if(any(! rownames(lm.orig) %in% featureNames(e)) ){stop("rownames of fits$fit not the same as in eset")}
			baseline.coef = gsub('.*-','',cont,perl=TRUE) # Attempts to parse contrast of style A-B
			if(baseline.coef %in% colnames(lm.orig)){
				print("Baseline scaling for the heatmap.")
				exprs(e) = exprs(e) - lm.orig$coefficient[featureNames(e),baseline.coef]
				scale="none" # thus expression values will be centered around bl but not scaled
			}
		}


		# Isolate the arrays
		initial.coefs   = names(which(cont.matrix[,cont]!=0))  # This extracts the Classes involved in a contrast
		sub.design      = lm$design[,initial.coefs,drop=FALSE] # This allows finding arrays involved in the contrast and their class
		involved.arrays = as.data.frame(which(sub.design!=0,arr.ind=TRUE),stringsAsFactors=FALSE)
		involved.arrays$initial.coef = sapply(involved.arrays$col,function(col)colnames(sub.design)[col])
		involved.arrays$color = color.me(involved.arrays$initial.coef) # The colors that will be used for columns
			
		# subset the eset objects to those arrays and add color label to pData
		e = e[,rownames(involved.arrays)] # This reorders exprs et the same time
			
		# Now time to subset on rows according to the level of DE
		top = topTable(lm, coef=cont, number=ngenes, genelist = lm$genes, adjust.method=adjust.method,lfc=lfc,p.value=p.value) # genelist=fit$genes,
		e = e[top$ProbeID,]

	
		# Scale height and width AND margins depending on label width
		labCol = apply(pData(e)[,alt.text,drop=FALSE],MARGIN=1,FUN=function(ro)paste(ro,collapse=" | "))
		labCol.length = max(nchar(labCol))
		labRow = fData(e)$SYMBOL
		labRow.length = max(nchar(labRow))
		pointsize=12
		mar.rows  = 2+labRow.length/1.7#labRow.length
		mar.cols  = 2+labCol.length/1.7#labCol.length
		height = 4+mar.rows/4 +3+nrow(e)/5
		width  = 0+mar.cols/4 +3+ncol(e)/5
		
		# Plot the png		
		png(file.path(output.dir,paste("heatmap_",cont.aliases[cont],".png",sep="")), height=height, width=width
				,pointsize=pointsize, res=150,units="in")
		occo = heatmap.2(exprs(e),main = cont.aliases[cont],col=colorpanel(30, "blue", "black", "red"),cexRow=1,cexCol=1.1, 
					ColSideColors=involved.arrays$color, scale=scale, trace="none",
					Colv=TRUE,symbreaks=symbreaks,labRow=labRow
					,labCol= labCol
					, margin=c(mar.cols,mar.rows)
		) 
		dev.off()

		# write the carpet
		carpet = t(occo$carpet)
		rownames(carpet) =  labRow[occo$rowInd]
		colnames(carpet) =  labCol[occo$colInd]
		carpet = as.data.frame(carpet[nrow(carpet):1,])	
		write.table(carpet,file=file.path(output.dir,paste("heatmap_",cont.aliases[cont],".xls",sep="")), sep="\t")
	    }


	
	# Now the F heatmap
	topF = topTableF(lm, number=ngenes, genelist = lm$genes, adjust.method=adjust.method,p.value=p.value) # genelist=fit$genes,
	if(ncol(lm)>1)
	{
		
		fcs = lm$coefficient[topF$ProbeID,]
		labCol = cont.aliases[colnames(fcs)]
		labCol.length = max(nchar(labCol))
		labRow = lm[rownames(fcs),]$genes$SYMBOL
		labRow.length = max(nchar(labRow))
		pointsize=12
		mar.rows  = 2+labRow.length/1.7
		mar.cols  = 2+labCol.length/1.7
		height =  4+mar.rows/4 +3+nrow(e)/5
		width  =  0+mar.cols/4 +3+ncol(e)/5
	
		png(file.path(output.dir,paste("heatmap_F",".png",sep="")), height=height, width=width
				,pointsize=pointsize, res=150,units="in")
		occo = heatmap.2(fcs,main = "F-test",col=colorpanel(30, "blue", "black", "red"),cexRow=1,cexCol=1.1, 
					scale=scale, trace="none",
					Colv=TRUE,symbreaks=T,labRow=labRow
					,labCol= labCol
					, margin=c(mar.cols,mar.rows)
		) 
		dev.off()
		# write the carpet
		carpet = t(occo$carpet)
		rownames(carpet) =  labRow[occo$rowInd]
		colnames(carpet) =  labCol[occo$colInd]
		carpet = as.data.frame(carpet[nrow(carpet):1,])	
		write.table(carpet,file=file.path(output.dir,paste("heatmap_F",".xls",sep="")), sep="\t")
	}else{
		png(file.path(output.dir,paste("heatmap_F",".png",sep="")) )
		plot.new() 
		dev.off()
		temp = data.frame("Clustering of coefficients does not make sense with only one comparison.")
		WriteXLS("temp",file.path(output.dir,paste("heatmap_F",".xls",sep="")),row.names=TRUE,AdjWidth = TRUE)
	}







	# F *Expression* heatmap
	initial.coefs                 = unique(rownames(which(cont.matrix!=0,arr.ind=TRUE)))
	sub.design                    = lm$design[,initial.coefs,drop=FALSE] # This allows finding arrays involved in the contrast and their class
	involved.arrays               = as.data.frame(which(sub.design!=0,arr.ind=TRUE),stringsAsFactors=FALSE)
	involved.arrays$initial.coefs = sapply(involved.arrays$col,function(col)colnames(sub.design)[col])
	involved.arrays$color         = color.me(involved.arrays$initial.coef) # The colors that will be used for columns
	e                             = x[topF$ProbeID,rownames(involved.arrays)]
	
	# Scale height and width AND margins depending on label width
	labCol = apply(pData(e)[,alt.text,drop=FALSE],MARGIN=1,FUN=function(ro)paste(ro,collapse=" | "))
	labCol.length = max(nchar(labCol))
	labRow = fData(e)$SYMBOL
	labRow.length = max(nchar(labRow))
	pointsize=12
	mar.rows  = 2+labRow.length/1.0#labRow.length
	mar.cols  = 2+labCol.length/1.5#labCol.length
	height =  4+mar.rows/4 +3+nrow(e)/5
	width  =  0+mar.cols/4 +3+ncol(e)/5
		
	# Plot the png		
	png(file.path(output.dir,paste("heatmap_F_expressions.png",sep="")), height=height, width=width
			,pointsize=pointsize, res=150,units="in")
	occo = heatmap.2(exprs(e),main = "",col=colorpanel(30, "blue", "black", "red"),cexRow=1,cexCol=1.1, 
				ColSideColors=involved.arrays$color, scale="row", trace="none",
				Colv=TRUE,symbreaks=TRUE,labRow=labRow
				,labCol= labCol
				, margin=c(mar.cols,mar.rows)
	) 
	dev.off()

	# write the carpet
	carpet = t(occo$carpet)
	rownames(carpet) =  labRow[occo$rowInd]
	colnames(carpet) =  labCol[occo$colInd]
	carpet = as.data.frame(carpet[nrow(carpet):1,])	
	write.table(carpet,file=file.path(output.dir,paste("heatmap_F_expressions.xls",sep="")), sep="\t")

	#X = exprs(e)
	#rownames(X) = fData(e)$SYMBOL
	#colnames(X) = apply(pData(e)[,alt.text,drop=FALSE],MARGIN=1,FUN=function(ro)paste(ro,collapse=" | "))
	#pheatmap(X, scale = "none",annotation = involved.arrays[,"Class",drop=FALSE],
        #       filename = file.path(output.dir,paste("heatmap_F_expressions",".pdf",sep="")) 
	#		)#, width = NA, height = NA)


	return(NULL)
}



# report.dir = dirReport; data.dir=dirData
report.exploratory <- function(report.dir, data.dir)
{
	files.dir = file.path(report.dir,'files')
	file.copy(file.path(data.dir,'exploratory_analysis.pdf'),files.dir)
	html = getHtmlBits(data.dir)$exploratory
	return(html)
}



# report.dir = dirReport; data.dir=dirData
report.preprocessing <- function(report.dir, data.dir)
{
	files.dir = file.path(report.dir,'files')

	x       = load.eset(data.dir)
	arrays  = getArraysAnn(data.dir)
	param   = attr(x,"param")
	html    = getHtmlBits()$preprocessing
	
	# Write the expression values file
	write.exprs(x,file=file.path(files.dir,"expression_values.txt"))

	# Write the paramerter values
	html = gsub('BG_VALUE'         ,param$bg,html,fixed=TRUE)
	html = gsub('BG_COUNT'   ,param$bgCount,html,fixed=TRUE)
	html = gsub('NORM_METHOD',param$normalizationMethod,html,fixed=TRUE)
	html = gsub('IQR_VALUE'  ,param$iqr,html,fixed=TRUE)
	html = gsub('NPROBES'    ,nrow(x),html,fixed=TRUE)
	html = gsub('NARRAYS'    ,nrow(arrays[arrays$outlierFlag=="",]),html,fixed=TRUE) # Before removinf technical replicates...
	if(!is.null(param$batch)){html = gsub('BATCH_METHOD'    ,paste("ComBat with ",param$batch," as a batch and ",paste(param$batch.covariates,collapse=", ")," as covariates",sep=""),html,fixed=TRUE)}else{html = gsub('BATCH_METHOD',"no batch correction",html,fixed=TRUE)}
	return(html)
}


#report.dir = dirReport;data.dir=dirData;third.id="DonorID"; removal.items=list("removed.arrays"=TRUE,"missing.samples"=TRUE)
report.removals <- function(report.dir, data.dir, removal.items, third.id)
{
	files.dir = file.path(report.dir,'files')
	
	html = getHtmlBits(data.dir)$removals

	# Join annotations... and reorder columns
	arrays  = getArraysAnn(data.dir)
	samples = getSamplesAnn(data.dir)
	if(!is.null(third.id)){
		third = getThirdAnn(data.dir);
		samples = merge(third,samples,by=colnames(third)[1],all.x=FALSE,all.y=TRUE)
	}
	arrays = merge(samples,arrays,by="SampleID",all.x=FALSE,all.y=TRUE)
	arrays = arrays[,c("ArrayID","outlierFlag",colnames(arrays)[!colnames(arrays)%in%c("ArrayID","outlierFlag")])]
	samples = samples[,c("SampleID",colnames(samples)[!colnames(samples)%in%c("SampleID")])]
	rownames(arrays) = arrays$ArrayID
	rownames(samples) = samples$SampleID

	#Secondary tables
	successful.arrays =  arrays[arrays$outlierFlag == "",]
	removed.arrays    =  arrays[arrays$outlierFlag != "",]
	missing.samples   =  samples[! samples$SampleID %in% successful.arrays$SampleID,]
	unhybed.samples   =  missing.samples[! missing.samples$SampleID %in% removed.arrays$SampleID,]
	removed.samples   =  missing.samples[  missing.samples$SampleID %in% removed.arrays$SampleID,]

	html = getHtmlBits(data.dir)$removals

	if(removal.items$removed.arrays)
	{
		x = list("removed_arrays" = removed.arrays)
		if(nrow(x$removed_arrays)>0){rownames(x$removed_arrays) = c(1:nrow(x$removed_arrays))}
		WriteXLS("x", file.path(files.dir,"removed_arrays.xls"), row.names =FALSE, AdjWidth = TRUE,BoldHeaderRow = TRUE)
		p=openPage(file.path(files.dir,"removed_arrays.html"),link.css='css/hwriter.css',title="Removed Arrays")
		hwrite('Removed Arrays', p, heading=3,center=TRUE)
		if(nrow(x$removed_arrays)==0){hwrite('None!', p,center=TRUE)}
		if(nrow(x$removed_arrays)>0){hwrite(x$removed_arrays,p,row.names=TRUE,row.bgcolor='#ddaaff',row.style=list('font-weight:bold'),br=TRUE)}
		closePage(p)
		file.copy(file.path(data.dir,"outliers_supporting_figures.pdf"),file.path(report.dir,"files"))
	}else{html = html[!grepl('removed_arrays.xls',html)]}
	if(removal.items$missing.samples)
	{
		x = list("unhybridized_samples" = unhybed.samples,"removed_samples"=removed.samples)
		if(nrow(x$unhybridized_samples)>0){rownames(x$unhybridized_samples) = c(1:nrow(x$unhybridized_samples))}
		if(nrow(x$removed_samples)>0){rownames(x$removed_samples) = c(1:nrow(x$removed_samples))}
		WriteXLS("x", file.path(files.dir,"missing_samples.xls"), row.names =FALSE, AdjWidth = TRUE,BoldHeaderRow = TRUE)
		p=openPage(file.path(files.dir,"missing_samples.html"),link.css='css/hwriter.css',title="Missing Samples")
		hwrite('Unhybridized Samples', p, heading=3,center=TRUE)
		if(nrow(x$unhybridized_samples)==0){hwrite('None!', p,center=TRUE)}
		if(nrow(x$unhybridized_samples)>0){hwrite(x$unhybridized_samples,p,row.names=TRUE,row.bgcolor='#ddaaff',row.style=list('font-weight:bold'),br=TRUE)}
		hwrite('Removed Samples', p, heading=3,center=TRUE,id="rs")
		if(nrow(x$removed_samples)==0){hwrite('None!', p,center=TRUE)}
		if(nrow(x$removed_samples)>0){hwrite(x$removed_samples,p,row.names=TRUE,row.bgcolor='#ddaaff',row.style=list('font-weight:bold'),br=TRUE)}
		closePage(p)
	}else{html = html[!grepl('missing_samples.xls',html)]}

	# unhybridized samples
	# samples missing as a consequence of array removal (hyridizaation could be reattempted)
	# sucessful samples


	# full annotation : the usual samples annotation
	# check if works with empmty data frames !!
	return(html)
}













# This function creates report files and returns the html of the annotation and raw data section
report.annotations <-function(report.dir,data.dir,annots.items, third.id)
{
	files.dir = file.path(report.dir,'files')
	html = 	getHtmlBits(data.dir)$annotations
	if(annots.items$raw){
		x.raw = load.raw(data.dir,third.id=NULL,strictHybs=FALSE,replace.zeroes=NULL)
		write.exprs(x.raw,file=file.path(files.dir,"raw.txt"))
		html = c(html,getHtmlBits(data.dir)$annotlineraw)
	}else{html = html[!grepl('raw.txt',html)]}
	if(annots.items$samples){
		x = list("Samples"=getSamplesAnn(data.dir)) # Write filtes and secondary htmls
		rownames(x$Samples)=c(1:nrow(x$Samples));
		if(!is.null(third.id)){
			x[[third.id]]=getThirdAnn(data.dir)
			rownames(x[[third.id]])=c(1:nrow(x[[third.id]]));
		}# "Donor will not always apply... anyway"
		WriteXLS("x", ExcelFileName = file.path(files.dir,"samples.xls"),row.names =FALSE,AdjWidth = TRUE,BoldHeaderRow = TRUE)
		p=openPage(file.path(files.dir,"samples.html"),link.css='css/hwriter.css',title="Samples")
		hwrite('Samples', p, heading=3,center=TRUE)
		hwrite(x$Samples,p,row.names=TRUE,row.bgcolor='#ddaaff',row.style=list('font-weight:bold'),br=TRUE)
		if(!is.null(third.id)){
			hwrite(third.id, p, heading=3,center=TRUE);
			hwrite(x[[third.id]],p,row.names=TRUE,row.bgcolor='#ddaaff',row.style=list('font-weight:bold'),br=TRUE)
		}
		closePage(p)
	}else{html = html[!grepl('samples.html',html)]}
		if(annots.items$arrays){
			x = list("Arrays"=getArraysAnn(data.dir))
			rownames(x$Arrays)=c(1:nrow(x$Arrays));
			WriteXLS("x", ExcelFileName = file.path(files.dir,"arrays.xls"),row.names =FALSE,AdjWidth = TRUE,BoldHeaderRow = TRUE)	
			p=openPage(file.path(files.dir,"arrays.html"),link.css='css/hwriter.css',title="Arrays")
			hwrite('Arrays', p, heading=3,center=TRUE)
			hwrite(x$Arrays,p,row.names=TRUE,row.bgcolor='#ddaaff',row.style=list('font-weight:bold'),br=TRUE)
			closePage(p)
		}else{html = html[!grepl('arrays.html',html)]}
		if(annots.items$probes){
			x = list("Probes"=getProbesAnn(data.dir))
			WriteXLS("x", ExcelFileName = file.path(files.dir,"probes.xls"),row.names =FALSE,AdjWidth = TRUE,BoldHeaderRow = TRUE)
		}else{html = html[!grepl('probes.xls',html)]}
		if(annots.items$gene.sets){
			x = list("GeneSets"=getGeneSets(data.dir))
			WriteXLS("x", ExcelFileName = file.path(files.dir,"gene_sets.xls"),row.names =FALSE,AdjWidth = TRUE,BoldHeaderRow = TRUE)
			html = c(html,getHtmlBits(data.dir)$annotlinegenesets)
		}else{html = html[!grepl('gene_sets.xls',html)]}
		return(html)
}


################################################################

# finalize.report : Versioning system for HTML reports.
#-----------------------------------------------------------------------------|
# Once finalized, the html report will be archived with a version tag.
# A file called .details will be hidden within html/qc folder witholding 
# the parameters used for the analysis, along with its description.
# The content will be zipped and named param$project_version.zip 
# The HTML folder will be copied to html_version as well.              
#-----------------------------------------------------------------------------|
# version : name to be appended to the report.
# author  : name of the analyst who generated it.
# desc    : analysis description to be written in the tracefile and the report.
# source  : Optional. where is the report located (default : 'html')
## TODO : Adjust the report zipped folder naming to be consistent with 
##        the new naming convention.
finalize.report <- function(version="v0", author=param$author,
desc="NA", source=dirReport) {
	
	if(file.exists(paste(source,version,sep="_")) || file.exists(
		paste(gsub("[^a-z0-9]",".",param$project, perl=T, 
		ignore.case=T), "_report_", version, ".zip", sep=""))) {
		
		print("Version ", version, " already exist! Try a different version 
			  name or delete the previous files.")
		return()
	}
	if(desc == "NA") {
		print("Please enter a detailed description for this report :")
		desc = scan(what="character", n=1, sep="^")
	}
	author = verify.author(author=author)
	tracefile = file(file.path(dirData,"tracefile.txt"),"a")
	cat(paste(format(Sys.time(), "%Y-%m-%d"),"\t",author,
			  "\tFinalized report : ", version, ". ",desc,'\n',sep=""), file=tracefile)
	close(tracefile)
	details = file(file.path(dirData, source,".details"),"w")
	cat("-----------| Analysis Parameters |-----------", file=details)
	for(i in 1:length(param)) {
		cat(paste(names(param)[i],"\t",param[i]),"\n",file=details)
	}
	cat("---------------------------------------------\nAuthor: ",
		file=details)
	cat(author, "\n",file=details)
	cat("Analysis version : ", version, file=details)
	cat("\nDescription : ", desc, "\n", file =details)
	cat("\nFinalized on : ", format(Sys.time(), "%Y-%m-%d %X"), 
		", using illupipe version ", illupipe.version, "\n", file =details)
	close(details)
	# System Commands -- for UNIX systems only. 
	system(paste("zip -9 -r ", gsub("[^a-z0-9]",".",param$project, perl=T, 
		ignore.case=T), "_report_", version, ".zip ", source, sep=""))
}

illupipe.report.version = "3.1 beta"
