
## SPIA test
#redo fit with updated JMP output
fit2 <- fit
fit <- importJMP_fit(fn='_Pairwise comparisons within groups_tbuk_tbsa_mal_diffs_V3.txt', 
                     featureNames='ProbeID', contrastNames=c('C-HP','C-TB_UKtest','C-TB_UKtr', 'PTB_tr-LTB_va_SA'))

fit$genes <- fit2$genes
rownames(fit$coefficients) <- rownames(fit2)
rownames(fit$p.value) <- rownames(fit2)
rownames(fit$fdr) <- rownames(fit2)
rownames(fit$t) <- rownames(fit2)
library(annDB,character.only=T)
library(SPIA)


results = decideTests(fit, )
spia.input <- list()


#names(spia.input[[1]]) <- getEG(names(spia.input[[1]]),annDB)

#sort by p-value or FDR

temp <- fit[sort(fit$fdr[,1]),]
allEG <- fit$genes[names(sort(fit$fdr[,1])),'Entrez_Gene_ID']
spia.input[[1]] <- fit$coefficients[names(sort(fit$fdr[,1])),1]
names(spia.input[[1]]) <- allEG

temp <- which(results[names(sort(fit$fdr[,1])),1]!=0)
spia.input[[1]] <- spia.input[[1]][temp]
spia.input[[1]] <- spia.input[[1]][which(!is.na(names(spia.input[[1]])))]
#names(spia.input[[1]]) <- getEG(names(spia.input[[1]]),annDB)
#names(spia.input[[1]]) <- as.vector(fit$genes[temp,'Entrez_Gene_ID'])

spia.res <- list()
#spia.input[[1]] <- spia.input[[1]][!is.na(fit$genes[temp,'Entrez_Gene_ID'])]
#getEG(names(spia.input[[1]]),annDB)
#[which(!is.na(getEG(rownames(fit$genes),annDB)))]
all <- allEG
#all <- getEG(rownames(fit$genes),annDB)
all <- all[!is.na(all)]
all <- as.character(all)
spia.input[[1]] <- spia.input[[1]][!duplicated(names(spia.input[[1]]))]
de <= spia.input[[1]]
spia.res[[1]] <- spia(de = de , all=all, verbose = TRUE, 
                      organism= "hsa", plots=TRUE, nB=2000)


