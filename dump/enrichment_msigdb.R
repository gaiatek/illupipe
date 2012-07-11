######################################################################
# annotation enrichment.of gene a cluster list                       #
# cmethod is passed to p.adjust                                      #
# all : list of all genes (SYMBOL) on the array                      #
# clusters : list of vectors of SYMBOLs of the genes in each cluster #
# pathways output of                                                 #
######################################################################
gene.annotation.enrichment <- function(all,clusters,pathways,cmethod="BH")
{
	p <- list()
# preprocessing: compile list of all genes in any pathway
	all.genes.pw <- list()
	for (i in 1:length(pathways))
	{
		genes.pw <- c(pathways[i],recursive=TRUE)
		genes.pw <- genes.pw[3:length(genes.pw)]
		all.genes.pw <- c(all.genes.pw,genes.pw)
	}
	all.genes.pw <- unique(all.genes.pw)
    g <- length(intersect(all,all.genes.pw))
	text<-paste("Testing ",length(clusters)," Clusters and",length(all.genes.pw)," Names in ",length(pathways)," Sets.")
	message(text)
# calculate p-values
	for (i in 1:length(clusters))
	{ 
		p_clust <- vector()
		cluster <- unique(clusters[[i]])
		n <- length(intersect(all.genes.pw,cluster))
		for (j in 1:length(pathways))
		{
			genes.pw <- c(pathways[j],recursive=TRUE)
			genes.pw <- genes.pw[3:length(genes.pw)]
			f <- length(intersect(all,genes.pw))
			ov <- length(intersect(genes.pw,cluster))
#			overlap <- c(overlap,ov)
			p_val <- phyper(ov-1, f, g-f, n ,lower.tail=FALSE)
			if (p_val<0) # bug in R 2.7.0
			{
			  text <- paste("Returned negative p-value: phyper(", ov-1, "," , f, "," , g-f, "," ,n, ",lower.tail=FALSE)") 
        warning(text)
			}
			p_clust <- c(p_clust,p_val)
		} 
		p_clust <- p.adjust(as.vector(p_clust),method=cmethod)
		p <- c(p,list(p_clust))
	}
	p
}

# load MSigDb *.gmt files
load.pathways <- function(filename)
{
	path.size <- count.fields(file=filename, sep = "\t", quote = "", skip = 0,blank.lines.skip = TRUE, comment.char = "")
	pathways <- list(scan(file=filename,what="character",sep="\t",quote="",n=path.size[1],strip.white = TRUE,quiet=TRUE))
	for (i in 2:length(path.size))
	{
		pw <- list(scan(file=filename,what="character",sep="\t",quote="",skip=i-1,n=path.size[i],strip.white = TRUE,quiet=TRUE))
		pathways <- c(pathways,pw)
	}
	pathways
}

get.pathnames.clean <- function(enrich,pathways,max_log_p)
{
	pathnames <- vector()
        	for (i in 1:length(enrich))
	{
		d <- log10(c(enrich[i],recursive=TRUE))
		path_indices <- which( (d < max_log_p) & (d != -Inf)  )
		cpath <-vector()
               
		for (j in 1:length(path_indices))
		{
			name <- c(pathways[path_indices[j] ],recursive=TRUE)[1]
                       
			cpath <- c(cpath,name)
                        
		}
		if (length(path_indices)>0){ pathnames <- c(pathnames,cpath)}
                
	}
	pathnames
        
}

get.pvalues.clean <- function(enrich,pathways,max_log_p)
{
	        pvalVec=vector()
	for (i in 1:length(enrich))
	{
		d <- log10(c(enrich[i],recursive=TRUE))
		path_indices <- which( (d < max_log_p) & (d != -Inf)  )
		pp=vector()
		for (j in 1:length(path_indices))
		{
			
                        pval<-c(enrich[[i]][path_indices[j]],recursive=TRUE)  
			
                        pp<-c(pp,pval)
		}
		
                if (length(path_indices)>0){ pvalVec <- c(pvalVec,pp)}
	}
	
        pvalVec
}




get.gene.set<- function(name,pathways)
{
	for (path in pathways)
	{
		if (path[1]==name)
			return(path[3:length(path)])
	}
}

