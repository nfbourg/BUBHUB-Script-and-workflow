level =.05

required_packages <- c("dplyr",'tibble','readr','pheatmap','gplots','ggplot2','fgsea')

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos='http://cran.us.r-project.org')
  sapply(pkg, require, character.only = TRUE)
}

readNfilt <- function(path, sffx, level){
  tib <- read_delim(path,'\t',col_names = TRUE)
  nm = paste('l2FC_',sffx,sep='')
  filt <- (tib %>% filter(padj<.05) %>% select('Gene_id','log2FoldChange') %>% rename(!!nm := 'log2FoldChange'))
  return(filt)
}
setwd("/projectnb/bubhub/users/nfbourg/cte_risk_variants_mrnaseq/samples/R_stats")
#setwd("/mnt/bubhub/cte/samples/R_stats/")
#setwd("/mnt/bubhub/samples/R_stats")


ipak(required_packages)

CvM_filt <- readNfilt('R_deseq_AgeAtDeath_RIN_Seq_Batch_Status_CNT_vs_CTE_12.tsv','CvM',level)
CvR_filt <- readNfilt('R_deseq_AgeAtDeath_RIN_Seq_Batch_Status_CNT_vs_CTE_34.tsv','CvR',level)
CvS_filt <- readNfilt('R_deseq_AgeAtDeath_RIN_Seq_Batch_Status_CNT_vs_RHI.tsv','CvS',level)
MvS_filt <- readNfilt('R_deseq_AgeAtDeath_RIN_Seq_Batch_Status_CTE_12_vs_CTE_34.tsv','MvS',level)

tf.venn <- (full <- full_join(CvM_filt, full_join(CvR_filt, full_join(CvS_filt, MvS_filt, by = 'Gene_id'), by = 'Gene_id'), by = 'Gene_id')) %>%
    mutate_at(.vars = vars(starts_with('l2FC')), .funs = funs(!is.na(.))) %>% 
    rename_at(.vars = vars(starts_with('l2FC')), .funs = funs(paste(substr(.,6,12),'_has',sep='')))

venn = venn(tf.venn[2:5])
dev.copy(pdf,'VennPlot.pdf')
dev.off()

gene_list <- (tf.venn %>% select('Gene_id'))
rd_cts <- (read_delim('./../norm_counts.csv',',',col_names = TRUE) %>% rename('Gene_id' = 'target_id'))

trim_rd_cts <- as.data.frame(left_join(gene_list,rd_cts,by = 'Gene_id') %>% filter_all(all_vars(. !=0)) %>% column_to_rownames('Gene_id')) %>% log10() 

#trim_rd_cts_2 <- as.data.frame(left_join(gene_list,rd_cts,by = 'Gene_id') %>% filter(!grepl("210082.2",Gene_id))%>% 
#                                 filter(!grepl("131095",Gene_id))%>% filter(!grepl("251562",Gene_id))%>% 
#                                 filter(!grepl("211459",Gene_id))%>% filter(!grepl("197971",Gene_id))%>% 
#                                 filter(!grepl("198727",Gene_id))%>% filter(!grepl("198763",Gene_id))%>% 
#                                 filter(!grepl("198695",Gene_id))%>% filter(!grepl("198840",Gene_id))%>% 
#                                 filter(!grepl("198804",Gene_id)) %>%filter(!grepl("19888",Gene_id)) %>% 
#                                 column_to_rownames('Gene_id')) %>% log10()
trim_rd_cts <- as.matrix(transform(trim_rd_cts, as.numeric()))
#trim_rd_cts_2 <- as.matrix(transform(trim_rd_cts_2, as.numeric()))
coldata <- read_delim('./../CTE_sample_info.csv',',') %>% select('ID','Status') %>% as.data.frame() %>% column_to_rownames('ID')


pheatmap(trim_rd_cts, cluster_rows= , show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col= coldata,cellheight = 20, cellwidth= 30, filename = 'heatmap_full.pdf')
#pheatmap(trim_rd_cts_2, cluster_rows= , show_rownames=TRUE,
#         cluster_cols=TRUE, annotation_col= coldata,cellheight = 20, cellwidth= 30, filename = '~/Documents/heatmap_trimmed.pdf')

