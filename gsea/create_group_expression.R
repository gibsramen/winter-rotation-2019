library(GSVA)
library(dplyr)

# input group file, genesets
# output ssGSEA scores?

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
	stop('Please provide arguments!')
}
setwd('/home/grahman/projects/rotation_2019/gsea')

# group file is initially in weird row format
group.file.loc <- args[1]
group <- read.csv(group.file.loc, header=F)
group <- group[c(1,3),]
t.group <- data.frame(t(group))
colnames(t.group) <- c('id', 'group')

# filter metadata to patients in cluster file
meta.file <- '/home/grahman/projects/rotation_2019/data/metadata_sid_filt_uniq.csv'
meta.df <- read.csv(meta.file)
group.meta.join <- inner_join(t.group, meta.df, by = c("id" = "X"))
group.meta.join <- group.meta.join[,c('id', 'case_id', 'group')]

geneset <- read.delim(args[2])
geneset <- data.frame(geneset[2:nrow(geneset),])
colnames(geneset) <- c('hgnc')

# match ensemble names to hgnc names
ensemble.file.loc <- '/home/grahman/data/TCGA_RNA_seq/ensemble_map.txt'
ensemble.map <- read.table(ensemble.file.loc, sep='\t', skip=2, header=F)
colnames(ensemble.map) <- c('ensemble', 'hgnc')
gene.match <- ensemble.map[ensemble.map$hgnc %in% geneset$hgnc,]
gene.match$ensemble <- as.character(gene.match$ensemble)

# contains file locations for expression data
sample.sheet.file.loc <- '/home/grahman/data/TCGA_RNA_seq/tcga_transcriptome_sample_sheet.csv'
sample.sheet <- read.csv(sample.sheet.file.loc)

group.meta.sample.join <- left_join(group.meta.join, sample.sheet,
									by = c("case_id" = "Case.ID"))
sample.df <- group.meta.sample.join

rna.dirs.loc <- '/home/grahman/data/TCGA_RNA_seq/rna_dirs'
print('Creating expression matrix...')
for (i in c(1:nrow(sample.df))){
	print(i)
	row <- sample.df[i,]
	exp.file <- paste(rna.dirs.loc, 
					  row[['File.ID']], 
					  row[['File.Name']], 
					  sep = '/')
	if (!file.exists(exp.file)){
		next
	}

	# first column is gene names (ensemble)
	# rest of the columns are expression of each tcga barcode (column name)
	this.exp.mat <- read.table(gzfile(exp.file), sep = '\t')
	colnames(this.exp.mat) <- c('ensemble', as.character(row[['case_id']]))
	this.exp.mat$ensemble <- gsub('\\..*', '', as.character(this.exp.mat$ensemble))
	#this.exp.mat <- semi_join(this.exp.mat, gene.match, by = c('ensemble'))
	if (i == 1){
		exp.mat <- this.exp.mat
		colnames(exp.mat)[1] <- c('ensemble')
	}
	else{
		exp.mat <- inner_join(exp.mat, this.exp.mat, by = 'ensemble')
	}
}
print('Finished creating expression matrix!')
write.csv(exp.mat, file = 'test_exp_mat.csv', row.names=F)
