library(ggplot2)

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
meta_df <- meta_df[,c('X', 'race', 'age_at_diagnosis', 'gender')]
cancer_meta_df <- meta_df[meta_df$X %in% t.group$id,]
cancer_meta_df$group <- t.group$group
cancer_meta_df$X <- NULL

gender_table <- table(cancer_meta_df[,c('group', 'gender')])
gender_fisher <- fisher.test(gender_table, simulate.p.value=TRUE)
print(gender_fisher)

race_table <- table(cancer_meta_df[,c('group', 'race')])
race_fisher <- fisher.test(race_table, simulate.p.value=TRUE)
print(race_fisher)

age_table <- cancer_meta_df[,c('group', 'age_at_diagnosis')]

age_aov <- aov(age_at_diagnosis ~ group, age_table)
print(summary(age_aov))

print(TukeyHSD(age_aov))
