#library(CancerSubtypes)
library(survival)
library(survminer)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
	stop('Please provide cancer type')
}

# run on input cancer type and input number of clusters
cancer_type <- args[1]
num_clust <- args[2]

cancer_dir <- paste('cluster_results', cancer_type, sep='/')
clust_string <- sprintf('%s_clusters', num_clust)
group_file <- paste(cancer_dir, clust_string, 'group.csv', sep='/')
group <- read.csv(group_file, header=F)

t.group <- t(group)
colnames(t.group) <- c('id', 'index', 'group')
t.group <- t.group[,c('id', 'group')]
t.group <- data.frame(t.group)

meta_df <- read.csv('../data/metadata_sid_filt.csv')
surv_df <- meta_df[,c('X', 'investigation' ,'vital_status_label', 
					  'days_to_death', 'new_tumor_event_after_initial_trtmt',
					  'primary_therapy_outcome_success_label')]
surv_df$investigation <- sub("TCGA-", '', surv_df$investigation)
cancer_surv_df <- surv_df[surv_df$investigation == cancer_type,]
cancer_surv_df <- cancer_surv_df[,c('X', 'vital_status_label', 'days_to_death')]
cancer_surv_df$vital_status_label <- ifelse(cancer_surv_df$vital_status_label
											== 'Dead', 1, 0)
t.group_cancer <- t.group[t.group$id %in% cancer_surv_df$X,]
cancer_surv_df <- cancer_surv_df[order(cancer_surv_df$X),]

#print(cancer_surv_df$days_to_death)
surv_object <- Surv(time = cancer_surv_df$days_to_death,
					event = cancer_surv_df$vital_status_label)
fit1 <- survfit(surv_object ~ group, data=t.group_cancer)
print(summary(fit1))
png('test.png')
ggsurvplot(fit1, data = t.group_cancer, pval = TRUE)
dev.off()
#print(surv_object)

