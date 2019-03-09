library(survival)
library(survminer)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
	stop('Please provide cancer type')
}

# TODO: take arg cluster directory
cluster_dir <- args[1]
cancer_types <- list.files(cluster_dir)

setwd('/home/grahman/projects/rotation_2019/stats')

days_df <- read.table('../clinical/follow_up_counts.txt', sep='\t')
colnames(days_df) <- c('barcode', 'status', 'days')
days_df$status <- as.character(days_df$status)

meta_df <- read.csv('../data/metadata_sid_filt_uniq.csv')
surv_df <- meta_df[,c('X', 'investigation' ,'vital_status_label', 
	'days_to_death', 'new_tumor_event_after_initial_trtmt',
	'primary_therapy_outcome_success_label', 'case_id')]
surv_df$investigation <- sub("TCGA-", '', surv_df$investigation)

cases <- surv_df[,c('case_id'),]
cases <- as.factor(cases)
barcodes <- c()
status_vector <- c()
days_vector <- c()

for (barcode in cases){
	row <- days_df[days_df$barcode == barcode,]
	barcodes <- c(barcodes, barcode)
	status_vector <- c(status_vector, row[[2]])
	days_vector <- c(days_vector, row[[3]])
}

for (cancer in cancer_types){
	print(cancer)
	cancer_dir <- paste(cluster_dir, cancer, sep='/')
	surv_df$status <- status_vector
	surv_df$days <- days_vector
	cancer_surv_df <- surv_df[surv_df$investigation == cancer,]
	cancer_surv_df <- cancer_surv_df[order(cancer_surv_df$X),]
	cancer_surv_df <- cancer_surv_df[,c('X', 'status', 'days')]
	cancer_surv_df$status <- ifelse(cancer_surv_df$status == 'Dead', 1, 0)

	for (i in c(2:10)){
		clust_string <- sprintf('%s_clusters', i)
		print(clust_string)
		num_clust_dir <- paste(cancer_dir, clust_string, sep='/')
		group_file_loc <- paste(num_clust_dir, 'group.csv', sep='/')
		group <- read.csv(group_file_loc, header=F)

		t.group <- t(group)
		colnames(t.group) <- c('id', 'index', 'group')
		t.group <- t.group[,c('id', 'group')]
		t.group <- data.frame(t.group)
		t.group_cancer <- t.group[t.group$id %in% cancer_surv_df$X,]

		cancer_surv_df_i <- cancer_surv_df
		cancer_surv_df_i$group <- t.group_cancer$group
		cancer_surv_df_i <- cancer_surv_df_i[,c('days', 'status', 'group')]

		surv_object <- Surv(time = cancer_surv_df_i$days,
			event=cancer_surv_df_i$status)
		fit1 <- survfit(surv_object ~ group, data=cancer_surv_df_i)

		out_file <- paste(num_clust_dir, 'survival_plot.png', sep='/')
		survp <- ggsurvplot(fit1, pval = TRUE)
		ggsave(file=out_file, print(survp), width=4, height=4, units="in")
	}
}
