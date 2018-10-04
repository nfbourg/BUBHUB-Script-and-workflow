#!/usr/bin/Rscript -- Nathanael Bourgeois
preliminary = FALSE
#Set Variables
#formulas = list('Status','AgeAtDeath','RIN','Seq_Batch',
formulas = list(list('AgeAtDeath','RIN','Seq_Batch','Status'))
required_packages <- c("dplyr",'tibble','readr','DESeq2','pheatmap','logistf','ggplot2')

#installed the requried packaged
ipak <- function(pkg){
        new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
        if (length(new.pkg))
                install.packages(new.pkg, dependencies = TRUE)
        sapply(pkg, require, character.only = TRUE)
}


#sub function 
Split_col <- function(cold_cts, str1, str2){
	cold_cts %>% filter((Status == str1) | (Status == str2)) -> filt
	as.data.frame(column_to_rownames(filt,"ID")) -> coldata
	return(coldata)
}

Split_cts <- function(tib_cts,colnames, str1, str2){
	tib_cts %>% select('target_id',colnames) -> cts_filt
	as.data.frame(column_to_rownames(cts_filt,'target_id')) -> rd_cts 
	return(rd_cts)
}

#DESeq counts vs trait 
DESeqtest <- function(rd_cts_t, coldata_t, expVar){
	#datasplit and df conversion
  rd_cts <- as.data.frame(column_to_rownames(rd_cts_t,'target_id'))
  coldata <- as.data.frame(column_to_rownames(coldata_t,"ID"))
	dsgn <- as.formula(paste(' ~ ',expVar, sep=''))
	dds <- DESeqDataSetFromMatrix(countData = rd_cts,colData = coldata, design = as.formula(dsgn))
	#DESeq
	keep <- rowSums(counts(dds)==0)==0 
	dds <- dds[keep,]
	dds <- DESeq(dds)
	res <- results(dds,cooksCutoff=FALSE)
	cf <- resultsNames(dds)[2]
	nm = paste('',expVar,collapse='_')
	name <- paste('R_res_counts_vs_',expVar, '.tsv', sep='')
	writeSeq(res, name)
	return()
}

DESeqAnal <- function(rd_cts, coldata, dsgn, prefiltNum,nm,str1, str2){
	dds <- DESeqDataSetFromMatrix(countData = rd_cts,colData = coldata, design = as.formula(dsgn))
	#prefilt
	keep <- rowSums(counts(dds)==0)==0 
	dds <- dds[keep,]
	#DESeq
	dds <- DESeq(dds,minReplicatesForReplace = Inf)
	res <- results(dds,cooksCutoff=FALSE)
	if(length(resultsNames(dds))==2){
		cf <- resultsNames(dds)[2]
	}
	else{
		cf <- resultsNames(dds)[4]
	}
	res_LFC <- lfcShrink(dds, coef=cf, type="apeglm")
	name <- paste('R_res_',nm,'_',str1,'_vs_',str2, '.tsv', sep='')
	#writeSeq(res, name)
	nameLFC <-  paste('R_deseq_',nm,'_',str1,'_vs_',str2, '.tsv', sep='')
	writeSeq(res_LFC,nameLFC)
	qplot(res$log2FoldChange,
	      geom="histogram",
	      binwidth = .5,  
	      main = "Histogram for Log2FoldChange", 
	      xlab = "Log2FoldChange",
	      ylab = "Counts (or Frequency)",
	      fill=I("dodgerblue1"), 
	      col=I("steelblue4"))
	ggsave(paste(name,'.jpeg',sep=''))
#	return(list(res,res_LFC))
}

DnF_Analysis <- function(rd_cts_t, coldata_t, prefiltNum, str1, str2, expVar){
	#datasplit and df conversion
	coldata <- Split_col(coldata_t,str1,str2)
	rd_cts <- Split_cts(rd_cts_t,rownames(coldata),str1,str2)
	dsgn = ''
	for (trait in expVar){
	  if (length(dsgn)>1){
	    dsgn = paste(dsgn,'+',trait)
	  }
	  else{
	    dsgn = paste('~',trait)
	  }
	}
	resfirth = logAnal(rd_cts, coldata, dsgn)
	#write tsv files
	nm = paste('',expVar,collapse='_',sep='')
	DESeqAnal(rd_cts, coldata, dsgn, prefiltNum, nm,str1,str2)
#	nameFirth <-  paste('R_res_Firth_',nm,'_',str1,'_vs_',str2, '.tsv', sep='')
#	writeTxt(resfirth, expVar, nameFirth, FALSE)
	return()
}

#linear regresion function
linReg <- function(coldata_t, expVar){
  if(expVar == 'Seq_Batch'){
    #	expVar = paste('factor(',expVar,')')
    return
  }
  formu <- as.formula(paste(expVar, ' ~ Status'))
  linMod <- summary(lm(formu, data=coldata_t))
  writeTxt(linMod, expVar, 'linMod.txt', TRUE)
}

logAnal <- function(cnts, coldata, dsgn){
  
}

writeSeq <- function(res, fprefix){
  resOrdered <- res[order(res$pvalue),]
	numP <- sum(res$padj < 0.05, na.rm=TRUE)
	pt05 <<- add_row(pt05,Pair=fprefix,Num_genes=numP)
	as.data.frame(resOrdered) %>% as_tibble(rownames='Gene_id') %>% write_tsv(fprefix)
	txt = paste(fprefix, numP)
	capture.output(cat(txt,'\n',sep=''), file='pt05-pvals.txt',append=TRUE)
	return
}

writeTxt <- function(txt, expVar, fname,APPND){
		title <- paste('------',toString(expVar),'~ Status','------ \n',sep=' ')
		n <- nchar(title)-2
		spcr <- paste(replicate(n, "-"), collapse = "")
		capture.output(cat(spcr, '\n', title, spcr, '\n', sep=''), file='linMod.txt',append=APPND)
		capture.output(txt, file='linMod.txt',append=TRUE)
}


setwd("/projectnb/bubhub/users/nfbourg/cte_risk_variants_mrnaseq/samples/R_stats")
args = commandArgs(trailingOnly=TRUE)
ipak(required_packages)
rd_cts_t = read_delim('./../counts.tsv','\t',col_names = TRUE)
coldata_t = read_delim('./../CTE_sample_info.csv',',')
for (form in formulas){
	expVar = unlist(form)
	if(preliminary==FALSE){
		if (!(is.list(form))){
			if (expVar!= 'Status' && expVar!='Seq_Batch')
				linReg(coldata_t, expVar)
			DESeqtest(rd_cts_t, coldata_t, expVar)
		}
	}
	txt = 'The number of genes with p value below .05 are:'
	pt05 <<- tibble(Pair = '--',Num_genes = '--')
	capture.output(cat(txt,'\n',sep=''), file='pt05-pvals.txt')
	DnF_Analysis(rd_cts_t, coldata_t, 0, 'CNT', 'RHI', expVar)
	DnF_Analysis(rd_cts_t, coldata_t, 0, 'CNT', 'CTE_12', expVar)
	DnF_Analysis(rd_cts_t, coldata_t, 0, 'CNT', 'CTE_34', expVar)
	DnF_Analysis(rd_cts_t, coldata_t, 0, 'RHI', 'CTE_12', expVar)
	DnF_Analysis(rd_cts_t, coldata_t, 0, 'RHI', 'CTE_34', expVar)
	DnF_Analysis(rd_cts_t, coldata_t, 0, 'CTE_12', 'CTE_34', expVar)
  pt05 %>% write_tsv('pt05_cts.txt')
}

