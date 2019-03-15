library(GSVA)
library(limma)
library(dplyr)

setwd('/home/grahman/projects/rotation_2019/gsea')

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
	stop('Please provide arguments!')
}

# group file is initially in weird row format
group.file.loc <- args[1]
group <- read.csv(group.file.loc, header=F)
group <- group[c(1,3),]
t.group <- data.frame(t(group))
colnames(t.group) <- c('id', 'group')

meta.file.loc <- '../data/metadata_sid_filt_uniq.csv'
meta.df <- read.csv(meta.file.loc, header=T)

t.group.meta <- inner_join(t.group, meta.df, by = c('id' = 'X'))
t.group.meta <- t.group.meta[,c('id', 'group', 'case_id')]
t.group.meta$case_id <- gsub('-', '.', t.group.meta$case_id)

mat.file <- 'test_exp_mat.csv'
exp.mat <- read.table(mat.file, header=T, sep=',', row.names=1)


ensemble.map.file <- '/home/grahman/data/TCGA_RNA_seq/ensemble_map.txt'
ensemble.map <- read.table(ensemble.map.file, sep='\t', header=F,
						   stringsAsFactors = F)
colnames(ensemble.map) <- c('ensemble', 'hgnc')

gene.file.1 <- '/home/grahman/data/genesets/BIOCARTA_INFLAM_PATHWAY.gmt'
geneset1 <- read.table(gene.file.1, sep='\t', header=F, skip=2,
					   stringsAsFactors = F)
colnames(geneset1) <- 'hgnc'
ensemble.geneset1.df <- semi_join(ensemble.map, geneset1, by = c('hgnc'))
ensemble.geneset1 <- as.character(ensemble.geneset1.df$ensemble)

gene.file.2 <- '/home/grahman/data/genesets/ADAPTIVE_IMMUNE_RESPONSE.gmt'
geneset2 <- read.table(gene.file.2, sep='\t', header=F, skip=2,
					   stringsAsFactors = F)
colnames(geneset2) <- 'hgnc'
ensemble.geneset2.df <- semi_join(ensemble.map, geneset2, by = c('hgnc'))
ensemble.geneset2 <- as.character(ensemble.geneset2.df$ensemble)

geneSets <- list(set1=unlist(ensemble.geneset1.df$ensemble),
				 set2=unlist(ensemble.geneset2.df$ensemble))

# fails when use all the genes? probably because og is filtered??
exp.mat.log <- log10(exp.mat + 0.000001)

out <- gsva(as.matrix(exp.mat.log), geneSets, method=c('ssgsea'), mx.diff=T,
			verbose=T)
write.csv(out, 'out.csv', row.names=F)
write.csv(t.group.meta, 'LGG_group.csv', row.names=F)
