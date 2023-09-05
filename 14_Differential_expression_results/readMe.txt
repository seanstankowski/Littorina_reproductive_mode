The scripts in this directory were used in the down stream analysis of the differential expression analysis. 

Expression_and_selection.ipynb:
This Jupyter notebook contains all of the code used to obtain the estimates used to test for the enrichment of DEGs in high Tr regions (Figs 4E and S35, and table).

diff_exp_merged_twisstPopgenwins:
This file contains the majority of the data needed to run  the analyses in Expression_and_selection.ipynb. the columns in the data frame are described below.

scaffold.x	contig ID in the assembly 
start.x	window start bp
end.x	window end bp
mid.x	window midpoint bp
win_size	window size bp
sites	number of SNPs in window
lnL	negative lnL from tree infrence
topoC_count	Twisst weight for the reproduction tree (tr in text) count
topoA_count	Twisst weight for the background tree (tb in text) count
topoB_couunt	Twisst weight for the control tree (tc in text) count
topoC	Twisst weight for the reproduction tree (tr in text) proportion
topoA	Twisst weight for the background tree (tb in text) proportion
topoB	Twisst weight for the control tree (tc in text) proportion
high_topo1	if tr weight > 0.7, == 1
scaffold.1	contig ID in the assembly 
start_check	window start bp
end_check	window end bp
mid_check	window midpoint bp
sites_full	total sites (invariant + variant)
pi_Egglayer	pi in egg-layers
pi_Brooder	pi in live-bearers
dxy_Egglayer_Brooder	dxy between egg-layers and live-bearers
Fst_Egglayer_Brooder	Fstxy between egg-layers and live-bearers
scaffold.2	window start bp
start_check_2	window end bp
end_check_2	window midpoint bp
mid_check_2	window size bp
sites_full_2	total sites (invariant + variant)
pi_arcana	pi in arcana
pi_compressa	pi in compressa
pi_n_sax	pi in northern sax
pi_spain	pi in iberian sax
dxy_arcana_compressa	dxy between species 
dxy_arcana_n_sax	dxy between species 
dxy_arcana_spain	dxy between species 
dxy_compressa_n_sax	dxy between species 
dxy_compressa_spain	dxy between species 
dxy_n_sax_spain	dxy between species 
Fst_arcana_compressa	fst between species
Fst_arcana_n_sax	fst between species
Fst_arcana_spain	fst between species
Fst_compressa_n_sax	fst between species
Fst_compressa_spain	fst between species
Fst_n_sax_spain	fst between species
LG	linkage groups
av	average position on the genetic map
LG_map_position	Lg and map position bound
inv_status	NA, unknown; buffer, in buffer area near boundary; colinear, not inside a known inversion; LGCxxx, ID of inversion that the window falls within.
in_map	1 == achored to the genetic map; NA, not anchored
row	Gene ID in the RNAseq analysis
baseMean.repr	average of the normalized count values, dividing by size factors, taken over all samples
log2FoldChange.repr	log fold change 
lfcSE.repr	standard error of log2foldchange
stat.repr	Wald statistic: the log2FoldChange divided by lfcSE
pvalue.repr	p-value for reproductive tissue
padj.repr	adjusted p-value for reproductive tissue
threshold.repr	NS, not significant; SD significant
baseMean.foot	average of the normalized count values, dividing by size factors, taken over all samples
log2FoldChange.foot	log fold change 
lfcSE.foot	standard error of log2foldchange
stat.foot	Wald statistic: the log2FoldChange divided by lfcSE
pvalue.foot	p-value for foor tissue
padj.foot	adjusted p-value for foot tissue
threshold.foot	NS, not significant; SD significant
differ	differentially expressed in foot, tissue, repro, both, none
seqname	contig ID in the assembly 
gene_start	start position of gene
gene_end	end position of gene 
scaffold.y	contig ID in the assembly 
start.y	window start bp
end.y	window end bp
mid.y	window midpoint bp

The RNA seq results in the above file are for the full analysis with all sequencing pools (see main text). However, we also conducted a reduced analysis. the file 'deseq_results_downsample' contains two columns: 

Row	Gene ID in the RNAseq analysis; same as for the full analysis
differ	differentially expressed in foot, tissue, repro, both, none; same as for the full analysis

Expression_and_selection.ipynb merges this file with the diff_exp_merged_twisstPopgenwins during the analysis. 