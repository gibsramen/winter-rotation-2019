library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
	stop('Please provide cancer type')
}

# run on cluster directory
cluster_dir <- args[1]
cancer_types <- list.files(cluster_dir, full.names = F)
meta_df <- read.csv('../data/metadata_sid_filt.csv')
meta_df <- meta_df[,c('X', 'race', 'age_at_diagnosis', 'gender')]

for (cancer in cancer_types){
	print(cancer)
	cancer_dir <- paste(cluster_dir, cancer, sep='/')
	for (i in c(2:10)){
		clust_string <- sprintf('%s_clusters', i)
		print(clust_string)
		group_file <- paste(cancer_dir, clust_string, 'group.csv', sep='/')
		group <- read.csv(group_file, header=F)

		t.group <- t(group)
		colnames(t.group) <- c('id', 'index', 'group')
		t.group <- t.group[,c('id', 'group')]
		t.group <- data.frame(t.group)

		cancer_meta_df <- meta_df[meta_df$X %in% t.group$id,]
		cancer_meta_df$group <- t.group$group
		cancer_meta_df$X <- NULL

		if (!cancer %in% c('OV', 'CESC', 'TGCT', 'PRAD', 'UCEC')){
			gender_table <- table(cancer_meta_df[,c('group', 'gender')])
			gender_fisher <- fisher.test(gender_table, simulate.p.value=TRUE)
			gender.p.val <- gender_fisher$p.value
		} else{
			gender.p.val <- 'NA'
		}
		race_table <- table(cancer_meta_df[,c('group', 'race')])
		race_fisher <- fisher.test(race_table, simulate.p.value=TRUE)
		race.p.val <- race_fisher$p.value

		age_table <- cancer_meta_df[,c('group', 'age_at_diagnosis')]
		age_aov <- aov(age_at_diagnosis ~ group, age_table)
		summ <- summary(age_aov)
		age.p.val <- summ[[1]]$"Pr(>F)"[1]

		p.vals <- c(race.p.val, age.p.val, gender.p.val)
		labels <- c('race', 'age', 'gender')
		out.df <- cbind(labels, p.vals)
		out.file <- paste(cancer_dir, clust_string, 'cov.csv', sep='/')
		write.csv(out.df, file=out.file, row.names=F)
	}
}
