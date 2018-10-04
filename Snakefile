import numpy as np
import itertools
import os
import pandas as pd

workdir: '../../samples'

csv = 'CTE_sample_info.csv'

SAMPLES = []
READ=[1,2]
PAIR=['U','P']
df = pd.read_csv(csv, delimiter=',')
for id in df['ID']:
	SAMPLES.append(id)
TRAITS = ['Status','AgeAtDeath', 'RIN', 'Seq_Batch', ['AgeAtDeath', 'RIN', 'Seq_Batch','Status']]
dat_cb=[]

#variables
filt = "0"

#datasplit
DATA = ['RHI','CNT','CTE_12','CTE_34']
comb = list(itertools.combinations(DATA,2))
info = pd.read_table('CTE_sample_info.csv', sep = ',')
dat_cbs = []
for cb in comb:
	dat_cb = cb[0] + '_' + cb[1]
	dat_cbs.append(dat_cb)


#indices
kal_ind_loc = "/projectnb/bubhub/bubhub-reference/genome/human/GENCODE/v27/gencode.v27.transcripts.kallisto_index"
star_ind_loc = "/projectnb/bubhub/bubhub-reference/genome/human/GENCODE/v27/GENCODE_v27_star_index"

#rules
rule start:
	input:
		"multiqc_report_latest.html",
		"counts.tsv",
		#"cts_poststats.json",
		"cts_prestats.html",
		#"pca.html",
		#"pca.json",
		"R_stats/vennplot.pdf",
		"R_stats/heatmap_full.pdf",
		expand("{sample}Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
		expand("stats/FirExp{dat}.tsv", dat=dat_cbs),
		expand("stats/deseq_{dat}.tsv", dat=dat_cbs)

rule multiqc:
	input:
		expand("{sample}_{read}{pair}_fastqc.zip", sample=SAMPLES, read=READ, pair=PAIR),
		expand("{sample}Log.final.out", sample=SAMPLES) 
	output:
		"multiqc_report_latest.html"
	shell:
		"export LC_ALL=en_US.utf8;"
		"export LANG=$LC_ALL;"
		"multiqc {input} -n {output} -f"

rule index:
	input:
		"{sample}Aligned.sortedByCoord.out.bam"
	output:
		"{sample}Aligned.sortedByCoord.out.bam.bai"
	shell:
		"samtools index {input}"

rule STAR:
	input:
		s1 = "{sample}_1P.fastq.gz",
		s2 = "{sample}_2P.fastq.gz"
	output:
		"{sample}Log.final.out",
		"{sample}Log.out",
		"{sample}Log.progress.out",
		"{sample}SJ.out.tab",
#		"{sample}_STARtmp/",
		"{sample}Aligned.sortedByCoord.out.bam"
	threads: 16
	shell:
		"STAR --outSAMtype BAM SortedByCoordinate --genomeDir {star_ind_loc} --runThreadN {threads} --readFilesIn {input.s1} {input.s2} --readFilesCommand zcat --outFileNamePrefix ./{wildcards.sample}"

#----Stats Block----#

rule DESeq:
	input:
		"counts_{dat}.tsv",
		"CTE_info_{dat}.csv"
	output:
		o1 = "stats/deseq_{dat}.tsv"
	threads:16
	run:
		cmd = "detk-de deseq2 --cores=16 'counts~ Status + RIN + AgeAtDeath + Seq_Batch' {input} -o {output.o1}"
		print(cmd)
		shell(cmd)

rule firth:
	input:
		"counts_{dat}.tsv",
		"CTE_info_{dat}.csv"
	output:
		o1 = "stats/FirExp{dat}.tsv"
	threads:16
	run:
		cmd = "detk-de firth 'Status ~ counts + RIN + AgeAtDeath + Seq_Batch' {input} -o {output.o1}"
		print(cmd)
		shell(cmd)


rule data_split:
	input:
		"norm_counts.csv"
	output:
                "counts_{dat}.tsv",
		"CTE_info_{dat}.csv"
	run:
		counts = pd.read_table('norm_counts.csv',sep = ',')
		for cb in comb:
		        info_split = info.loc[info['Status'].isin(cb)]
		        pathInfo = 'CTE_info_' + cb[0] + '_' + cb[1] +'.csv'
		        info_split.to_csv(pathInfo, index=False)
		        cols = ['target_id']
		        cols.extend(list(info_split['ID']))
		        counts_split = counts[cols]
			pathC = 'counts_' + cb[0] + '_' + cb[1] +'.tsv'
		        counts_split.to_csv(pathC, sep='\t', index=False)


rule post_stats:
	input:
		"norm_counts.csv"
	output:
		"cts_poststats.json"
	shell:
		"detk-stats coldist --json=tmp1.json {input};cat tmp1.json  >> {output}; rm tmp1.json;"
		"detk-stats coldist --json=tmp1.json --log {input}; cat tmp1.json >> {output}; rm tmp1.json;"
		"detk-stats rowzero --json=tmp1.json {input}; cat tmp1.json >> {output}; rm tmp1.json;"
		"detk-stats entropy --json=tmp1.json {input}; cat tmp1.json >> {output}; rm tmp1.json;"
		#"detk-stats --json=tmp.json pca -m {tsv} -f 'Subject.subject_type' {input}; cat tmp.json >> {output};"

rule norm:
	input:
		"counts_filt.tsv"
	output:
		"norm_counts.csv",
	shell:
		"detk-norm deseq2 {input} -o {output}"

rule filter:
	input:
		"counts_filt0.tsv",
	output:
		"counts_filt.tsv"
	shell:
		"detk-filter 'zeros(all) == 0' {input} -o {output}" 

rule pca_json:
	input:
		"counts_filt0.tsv",
	output:
		"pca.json"
	shell:
		"detk-stats  pca --json={output} --column-data=CTE_sample_info.csv --color-col='Status' {input}"

rule pca_html:
	input:
		"counts_filt0.tsv",
	output:
		"pca.html"
	shell:
		"detk-stats pca --html={output} --column-data=CTE_sample_info.csv --color-col='Status' {input}"

rule pre_stats_json:
	input:
		"counts_filt0.tsv"
	output:
		"cts_prestat.json"
	shell:
		"detk-stats coldist --json=tmp.json {input};cat tmp.json  >> {output}; rm tmp.json;"
		"detk-stats coldist --json=tmp.json --log {input}; cat tmp.json >> {output}; rm tmp.json;"
		#"detk-stats rowzero --json=tmp.json {input}; cat tmp.json >> {output}; rm tmp.json;"
		"detk-stats entropy --json=tmp.json {input}; cat tmp.json >> {output}; rm tmp.json;"
		#"detk-stats pca {in --json=tmp.jsonput}; cat tmp.json >> {output}; rm tmp.json;"

rule pre_stats_html:
	input:
		"counts_filt0.tsv"
	output:
		"cts_prestats.html"
	shell:
		"detk-stats coldist --html=tmp.html {input}; cat tmp.html  >> {output}; rm tmp.html;"
		"detk-stats coldist --html=tmp.html --log {input}; cat tmp.html >> {output}; rm tmp.html;"
		#"detk-stats rowzero --html=tmp.html {input}; cat tmp.html >> {output}; rm tmp.html;"
		"detk-stats entropy --html=tmp.html {input}; cat tmp.html >> {output}; rm tmp.html;"
		#"detk-stats pca {in --json=tmp.jsonput}; cat tmp.json >> {output}; rm tmp.json;"

rule filter0:		
	input:
		"counts.tsv",
	output:
		"counts_filt0.tsv"
	shell:
		"detk-filter 'mean(all) > 0' {input} -o {output}" 
rule R_anal:
	input:
		"norm_counts.csv",
		"R_stats/R_complete.txt"
	output:
		"R_stats/vennplot.pdf",
		"R_stats/heatmap_full.pdf"
	shell:
		"Rscript ../analysis/processing/R_analysis.R"

rule R:
	input:
		"counts.tsv"
	output:
		"R_stats/R_complete.txt"
	threads: 4
	shell:
		"Rscript ../analysis/processing/norm_counts.R;	echo 'True.' > {output}"

#----End Block----#

rule tran_2_gene:
	input:
		"t_counts.tsv"
	output:
		"counts.tsv"
	run:
		filei = open(input[0],'r')
		header = filei.readline()
		c_dict = {}
		for line in filei:
			line = line.split('|')
			gene = line[1]
			l_cts = [float(num) for num in line[-1].split()]
			if gene in c_dict.keys():
				for i in range(len(l_cts)):
					c_dict[gene][i] += l_cts[i]
			else:
				c_dict[gene] = l_cts
		filei.close()
		fileo = open(output[0],'w')
		fileo.writelines(header)
		for key in sorted(c_dict.keys()):
			fileo.writelines(key)
			for num in c_dict[key]:
				line = '\t' + str(int(num)) 
				fileo.write(line) 
			fileo.write('\n')
		fileo.close()
		
rule csvgather:
	input:
		expand("{sample}_kallisto_quant/abundance.tsv", sample = SAMPLES)
	output:
		"t_counts.tsv"
	shell:
		"csvgather -d '\t' -j 0 -f est_counts -t 's:est_counts:{{dir}}:' -t 's:_kallisto_quant::' {input} -o {output}"


rule fastqc:
	input:
		"{sample}_{read}{pair}.fastq.gz",
	output:
		"{sample}_{read}{pair}_fastqc.zip"
#		"{sample}_{read}{pair}_fastqc.html"
	threads: 1
	shell:
		"fastqc -t {threads} {input}" 

rule kallisto:
	input:
		s1 = "{sample}_1P.fastq.gz",
		s2 = "{sample}_2P.fastq.gz"
	output:
#		dir = "{sample}_kallisto_quant",
		file = "{sample}_kallisto_quant/abundance.tsv"
	threads: 16	
	shell:
		"kallisto quant -t {threads} -i {kal_ind_loc} {input.s1} {input.s2} -o '{wildcards.sample}_kallisto_quant'"
		
rule trim:
	output:
		fp = "{sample}_1P.fastq.gz",
		fu = "{sample}_1U.fastq.gz",
		bp = "{sample}_2P.fastq.gz",
		bu = "{sample}_2U.fastq.gz"
	threads: 16
	shell:
		"trimmomatic PE -phred33 {wildcards.sample}*1.fastq.gz {wildcards.sample}*2.fastq.gz {output.fp} {output.fu} {output.bp} {output.bu} -threads {threads} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" 
