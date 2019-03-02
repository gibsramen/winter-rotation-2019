library(dplyr)
library(reshape2)

# take in cluster directory
# return ANOVA of all microbes by group
# probably only return significant values

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2){
	stop('Please provide arguments!')
}

cluster_dir <- args[1]
cancer_types <- list.files(cluster_dir, full.names = F)
feature_csv <- args[2]

print('Loading data frames...')
otu_df <- read.csv(feature_csv)

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
		otu_df_filt <- otu_df[otu_df$Unnamed..0 %in% t.group_df$id,]
		otu_group_df <- data.frame(otu_df_filt)
		otu_group_df <- otu_group_df[order(otu_group_df$Unnamed..0),]
		otu_group_df$group <- t.group_df$group
		otu_group_df$Unnamed..0 <- NULL
		melt.otu_group_df <- melt(otu_group_df, id.vars = c('group'))

		p.vals <- c()
		for (otu in colnames(otu_group_df)){
			if (otu == 'group'){ next }
			otu.otu_group_df <- otu_group_df[,c('group', otu)]
			colnames(otu.otu_group_df) <- c('group', 'value')
			aov <- aov(value ~ group, data = otu.otu_group_df)
			summ <- summary(aov)
			p.val <- summ[[1]]$"Pr(>F)"[1]
			p.vals <- c(p.vals, p.val)
			tk <- TukeyHSD(aov)
		}

		adj.p.vals <- p.adjust(p.vals, method = "bonferroni")

		all.otus <- colnames(otu_group_df)[!colnames(otu_group_df) %in% c('group')]
		out.df <- data.frame(cbind(all.otus, p.vals, adj.p.vals))
		colnames(out.df) <- c('OTU', 'p-vals', 'Bonf-adj-p-vals')

		out_file <- paste(num_clust_dir, 'otu_anova.csv', sep='/')
		write.csv(out.df, file = out_file, row.names=F)
	}
}
