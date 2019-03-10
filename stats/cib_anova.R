library(dplyr)
library(reshape2)

# take in group file and cibersort csv as inputs
# return ANOVA of all microbes/cibersort by group
# probably only return significant values

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2){
	stop('Please provide arguments!')
}

setwd('/home/grahman/projects/rotation_2019/stats')
cluster_dir <- args[1]
cancer_types <- list.files(cluster_dir, full.names = F)
cib_file <- args[2]

print('Loading data frames...')
cib_df <- read.csv(cib_file)

for (cancer in cancer_types){
	print(cancer)
	cancer_dir <- paste(cluster_dir, cancer, sep='/')
	for (i in c(2:10)){
		clust_string <- sprintf('%s_clusters', i)
		print(clust_string)
		num_clust_dir <- paste(cancer_dir, clust_string, sep='/')
		group_file_loc <- paste(num_clust_dir, 'group.csv', sep='/')
		group_df <- read.csv(group_file_loc, header=F)

		t.group_df <- t(group_df)
		colnames(t.group_df) <- c('id', 'index', 'group')
		t.group_df <- t.group_df[,c('id', 'group')]
		t.group_df <- data.frame(t.group_df)
		t.group_df <- t.group_df[order(t.group_df$id),]

		print('Filtering data frames...')
		cib_df_filt <- cib_df[cib_df$feature.id %in% t.group_df$id,]
		cib_group_df <- data.frame(cib_df_filt)
		cib_group_df <- cib_group_df[order(cib_group_df$feature.id),]
		cib_group_df$group <- t.group_df$group
		cib_group_df$feature.id <- NULL
		cib_group_df$CancerType <- NULL
		melt.cib_group_df <- melt(cib_group_df, id.vars = c('group'))

		p.vals <- c()
		for (cell_type in colnames(cib_group_df)){
			if (cell_type == 'group'){ next }
			cell.cib_group_df <- cib_group_df[,c('group', cell_type)]
			colnames(cell.cib_group_df) <- c('group', 'value')
			aov <- aov(value ~ group, data = cell.cib_group_df)
			summ <- summary(aov)
			p.val <- summ[[1]]$"Pr(>F)"[1]
			p.vals <- c(p.vals, p.val)
			tk <- TukeyHSD(aov)
		}

		adj.p.vals <- p.adjust(p.vals, method="fdr")

		cell_types <- colnames(cib_group_df)[!colnames(cib_group_df) %in% c('group')]
		out.df <- data.frame(cbind(cell_types, p.vals, adj.p.vals))
		colnames(out.df) <- c('Cell.Type', 'p-vals', 'FDR-adj-p-vals')
		out.df <- out.df[order(out.df$'p-vals'),]

		out_file <- paste(num_clust_dir, 'cib_anova.csv', sep='/')
		write.csv(out.df, file = out_file, row.names=F)
	}
}
