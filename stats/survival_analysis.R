library(survival)
library(survminer)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
	stop('Please provide cancer type')
}

# run on input cancer type and input number of clusters
cancer_type <- args[1]
num_clust <- args[2]

project_dir <- '/home/grahman/projects/rotation_2019/cancer_subtypes'
cancer_dir <- paste('cluster_results', cancer_type, sep='/')
clust_string <- sprintf('%s_clusters', num_clust)
group_file <- paste(project_dir, cancer_dir, clust_string, 'group.csv', sep='/')
group <- read.csv(group_file, header=F)

t.group <- t(group)
colnames(t.group) <- c('id', 'index', 'group')
t.group <- t.group[,c('id', 'group')]
t.group <- data.frame(t.group)

days_df <- read.table('../clinical/follow_up_counts.txt', sep='\t')
colnames(days_df) <- c('barcode', 'status', 'days')
days_df$status <- as.character(days_df$status)

meta_df <- read.csv('../data/metadata_sid_filt.csv')
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

surv_df$status <- status_vector
surv_df$days <- days_vector
cancer_surv_df <- surv_df[surv_df$investigation == cancer_type,]
cancer_surv_df <- cancer_surv_df[order(cancer_surv_df$X),]
cancer_surv_df <- cancer_surv_df[,c('X', 'status', 'days')]
cancer_surv_df$status <- ifelse(cancer_surv_df$status == 'Dead', 1, 0)
print(head(cancer_surv_df))

t.group_cancer <- t.group[t.group$id %in% cancer_surv_df$X,]
print(head(t.group_cancer))

surv_object <- Surv(time = cancer_surv_df$days,
					event = cancer_surv_df$status,
					type = 'right')
fit1 <- survfit(surv_object ~ group, data=t.group_cancer)

out_dir <- paste('survival_plots', cancer_type, sep='/')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- paste(out_dir, sprintf('%s_clusters.png', num_clust), sep='/')
png(out_file)
p <- ggsurvplot(fit1, data = t.group_cancer, pval = TRUE)
dev.off()


