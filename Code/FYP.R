library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(predictiveFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(h2o)
library(ggplot2)
library(readxl)
library(GenomicRanges)
library(tidyr)
library(dplyr)
library(writexl)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(gridExtra)
library(corrplot)
library(VennDiagram)


data_dir <- "/home/zhen/FYP"

# Create positive and negative labels
m6A_gr <- readRDS(file.path(data_dir, "motif_exon_hg38_gr_atlas2.rds"))
m6A_gr$Mdl_positive_indx <- m6A_gr$Technique_Num >= 2
mature_mRNA <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "tx")
m6A_gr$Mdl_negative_indx <- m6A_gr%over%subsetByOverlaps(mature_mRNA, m6A_gr[m6A_gr$Mdl_positive_indx])
rm(mature_mRNA)
m6A_gr$Mdl_negative_indx[m6A_gr$Mdl_positive_indx] <- FALSE
m6A_gr$Mdl_negative_balanced_indx <- FALSE
set.seed(2318)
m6A_gr$Mdl_negative_balanced_indx[sample(which(m6A_gr$Mdl_negative_indx), sum(m6A_gr$Mdl_positive_indx))] <- TRUE

genomicFeatures_gr <- readRDS(file.path(data_dir, "motif_exon_hg38_gr_gfeatures.rds"))
genome_pos <- mcols(genomicFeatures_gr)[m6A_gr$Mdl_positive_indx,]
genome_neg <- mcols(genomicFeatures_gr)[m6A_gr$Mdl_negative_balanced_indx,]
rm(genomicFeatures_gr)
gc()
library(BSgenome.Hsapiens.UCSC.hg38)
seq_pos <- sequenceDerivedFeatures(m6A_gr[m6A_gr$Mdl_positive_indx]+20, BSgenome.Hsapiens.UCSC.hg38, encoding = "onehot")
seq_neg <- sequenceDerivedFeatures(m6A_gr[m6A_gr$Mdl_negative_balanced_indx]+20, BSgenome.Hsapiens.UCSC.hg38, encoding = "onehot")


library(h2o)
h2o.init()
# h2o.removeAll()
model_matrix <- rbind(cbind(seq_pos, genome_pos), cbind(seq_neg, genome_neg))
model_matrix$Y <- as.factor(rep(c(1,0), each = sum(m6A_gr$Mdl_positive_indx)))
training_frame <- as.h2o(model_matrix)

seq_only_model_xgboost <- h2o.xgboost(x = colnames(seq_pos), y = "Y", training_frame = training_frame, nfolds = 10)
seq_only_model_xgboost 
add_gen_model_xgboost <- h2o.xgboost(x = colnames(model_matrix), y = "Y", training_frame = training_frame, nfolds = 10)
add_gen_model_xgboost

# Calculate AUC
perf_seq_only_model_xgboost <- h2o.performance(seq_only_model_xgboost, training_frame)
perf_add_gen_model_xgboost <- h2o.performance(add_gen_model_xgboost, training_frame)


library(ggplot2)
# Get ROC curve data
roc_seq_only_xgboost <- data.frame(fpr = perf_seq_only_model_xgboost@metrics$thresholds_and_metric_scores$fpr,
                           tpr = perf_seq_only_model_xgboost@metrics$thresholds_and_metric_scores$tpr)
roc_add_gen_xgboost <- data.frame(fpr = perf_add_gen_model_xgboost@metrics$thresholds_and_metric_scores$fpr,
                          tpr = perf_add_gen_model_xgboost@metrics$thresholds_and_metric_scores$tpr)

# Plot ROC curves
auc1 <- ggplot() +
  geom_line(data = roc_seq_only_xgboost, aes(x = fpr, y = tpr, color = "Sequence feature only"), size = 1) +
  geom_line(data = roc_add_gen_xgboost, aes(x = fpr, y = tpr, color = "Sequence + Genomic feature"), size = 1) +
  scale_color_manual(values = c("Sequence feature only" = "red", "Sequence + Genomic feature" = "green")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  annotate("text", x = 0.5, y = 0.75, label = paste0(round(perf_seq_only_model_xgboost@metrics$AUC, 3)), color = "red") +
  annotate("text", x = 0.25, y = 0.98, label = paste0(round(perf_add_gen_model_xgboost@metrics$AUC, 3)), color = "green") +
  labs(color = "Model") + 
  theme_bw()
auc1


## general linear model
seq_only_model_glm <- h2o.glm(y = "Y", x = colnames(seq_pos), training_frame = training_frame, family = "binomial", nfolds = 10) 
seq_only_model_glm 
add_gen_model_glm <- h2o.glm(y = "Y", x = colnames(model_matrix), training_frame = training_frame, family = "binomial", nfolds = 10) 
add_gen_model_glm 

perf_seq_only_model_glm <- h2o.performance(seq_only_model_glm, training_frame)
perf_add_gen_model_glm <- h2o.performance(add_gen_model_glm, training_frame)
roc_seq_only_glm <- data.frame(fpr = perf_seq_only_model_glm@metrics$thresholds_and_metric_scores$fpr,
                                   tpr = perf_seq_only_model_glm@metrics$thresholds_and_metric_scores$tpr)
roc_add_gen_glm <- data.frame(fpr = perf_add_gen_model_glm@metrics$thresholds_and_metric_scores$fpr,
                                  tpr = perf_add_gen_model_glm@metrics$thresholds_and_metric_scores$tpr)
auc2 <- ggplot() +
  geom_line(data = roc_seq_only_glm, aes(x = fpr, y = tpr, color = "Sequence feature only"), size = 1) +
  geom_line(data = roc_add_gen_glm, aes(x = fpr, y = tpr, color = "Sequence + Genomic feature"), size = 1) +
  scale_color_manual(values = c("Sequence feature only" = "red", "Sequence + Genomic feature" = "green")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  annotate("text", x = 0.5, y = 0.7, label = paste0(round(perf_seq_only_model_glm@metrics$AUC, 3)), color = "red") +
  annotate("text", x = 0.25, y = 0.95, label = paste0(round(perf_add_gen_model_glm@metrics$AUC, 3)), color = "green") +
  labs(color = "Model") + 
  theme_bw()
auc2


## random forest
# h2o.removeAll()
n <- nrow(model_matrix) 
set.seed(3038)
idx_rf <- sample(1:n, floor(0.7*n))  
train_data <- model_matrix[idx_rf, ]
test_data <- model_matrix[-idx_rf, ]
training_frame2 <- as.h2o(train_data)
training_frame3 <- as.h2o(test_data)
seq_only_model_rf <- h2o.randomForest(y = "Y", x = colnames(seq_pos), training_frame = training_frame2, nfolds = 10) 
seq_only_model_rf 
add_gen_model_rf <- h2o.randomForest(y = "Y", x = colnames(model_matrix), training_frame = training_frame2, nfolds = 10) 
add_gen_model_rf 

# seq_only_model_rf <- h2o.randomForest(y = "Y", x = colnames(seq_pos), training_frame = training_frame)
# seq_only_model_rf
# add_gen_model_rf <- h2o.randomForest(y = "Y", x = colnames(model_matrix), training_frame = training_frame)
# add_gen_model_rf

perf_seq_only_model_rf <- h2o.performance(seq_only_model_rf, training_frame3)
perf_add_gen_model_rf <- h2o.performance(add_gen_model_rf, training_frame3)
roc_seq_only_rf <- data.frame(fpr = perf_seq_only_model_rf@metrics$thresholds_and_metric_scores$fpr,
                               tpr = perf_seq_only_model_rf@metrics$thresholds_and_metric_scores$tpr)
roc_add_gen_rf <- data.frame(fpr = perf_add_gen_model_rf@metrics$thresholds_and_metric_scores$fpr,
                              tpr = perf_add_gen_model_rf@metrics$thresholds_and_metric_scores$tpr)
auc3 <- ggplot() +
  geom_line(data = roc_seq_only_rf, aes(x = fpr, y = tpr, color = "Sequence feature only"), size = 1) +
  geom_line(data = roc_add_gen_rf, aes(x = fpr, y = tpr, color = "Sequence + Genomic feature"), size = 1) +
  scale_color_manual(values = c("Sequence feature only" = "red", "Sequence + Genomic feature" = "green")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  annotate("text", x = 0.5, y = 0.7, label = paste0(round(perf_seq_only_model_rf@metrics$AUC, 3)),color = "red") +
  annotate("text", x = 0.25, y = 0.9, label = paste0(round(perf_add_gen_model_rf@metrics$AUC, 3)), color = "green") +
  labs(color = "Model") + 
  theme_bw()
auc3


## multilayer perceptron
set.seed(3039)
seq_only_model_mlp <- h2o.deeplearning(x = colnames(seq_pos),
                              y = "Y",
                              training_frame = training_frame,
                              hidden = c(32, 64, 128), 
                              epochs = 10,
                              activation = "Rectifier",
                              input_dropout_ratio = 0.1, nfolds = 10) 
seq_only_model_mlp 
add_gen_model_mlp <- h2o.deeplearning(x = colnames(model_matrix), 
                                       y = "Y",  
                                       training_frame = training_frame, 
                                      hidden = c(32, 64, 128), 
                                       epochs = 10,  
                                       activation = "Rectifier",
                                      input_dropout_ratio = 0.1, nfolds = 10)  
add_gen_model_mlp 

perf_seq_only_model_mlp <- h2o.performance(seq_only_model_mlp, training_frame)
perf_add_gen_model_mlp <- h2o.performance(add_gen_model_mlp, training_frame)
roc_seq_only_mlp <- data.frame(fpr = perf_seq_only_model_mlp@metrics$thresholds_and_metric_scores$fpr,
                              tpr = perf_seq_only_model_mlp@metrics$thresholds_and_metric_scores$tpr)
roc_add_gen_mlp <- data.frame(fpr = perf_add_gen_model_mlp@metrics$thresholds_and_metric_scores$fpr,
                             tpr = perf_add_gen_model_mlp@metrics$thresholds_and_metric_scores$tpr)
auc4 <- ggplot() +
  geom_line(data = roc_seq_only_mlp, aes(x = fpr, y = tpr, color = "Sequence feature only"), size = 1) +
  geom_line(data = roc_add_gen_mlp, aes(x = fpr, y = tpr, color = "Sequence + Genomic feature"), size = 1) +
  scale_color_manual(values = c("Sequence feature only" = "red", "Sequence + Genomic feature" = "green")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  annotate("text", x = 0.5, y = 0.7, label = paste0(round(perf_seq_only_model_mlp@metrics$AUC, 3)), color = "red") +
  annotate("text", x = 0.25, y = 0.95, label = paste0(round(perf_add_gen_model_mlp@metrics$AUC, 3)), color = "green") +
  labs(color = "Model") + 
  theme_bw()
auc4





auc1 <- auc1 + theme(legend.position="none")
auc2 <- auc2 + theme(legend.position="none") 
auc3 <- auc3 + theme(legend.position="none") 
auc4 <- auc4 + theme(legend.position="none") 

auc_plot <- grid.arrange(auc1, auc2, auc3, auc4, ncol = 2)

ggsave("auc_plot.png", auc_plot, dpi = 300, width = 6, height = 6)

# IUPAC_CODE_MAP



# zhang

m6AQTL <- readRDS("m6AQTL.m6APeak_logOR_GC.IP.adjusted_qqnorm.15PCs.fastQTL.nominals.rds")

#Filter by q values
m6AQTL <- m6AQTL[m6AQTL$qvalue < 0.05,]

#Filter by ref&alt lengths
m6AQTL <- m6AQTL[nchar(as.character(m6AQTL$REF))==1 & nchar(as.character(m6AQTL$ALT))==1,]

#Filter by distance
m6AQTL <- m6AQTL[abs(m6AQTL$DIST) <= 20,]

#Create GRanges of peaks
startEND <- gsub(".*:|_.*", "", m6AQTL$PEAK)
peaks <- GRanges(strand = gsub(".*_", "", m6AQTL$PEAK),
                 seqnames = gsub(":.*", "", m6AQTL$PEAK),
                 ranges = IRanges(start = as.numeric(gsub("-.*", "", startEND)),
                                  end = as.numeric(gsub(".*-", "", startEND))))
mcols(peaks) <- DataFrame(m6AQTL[,-1*c(1,2)])

#Save as cisQTL_peaks_gr.rds
saveRDS(peaks, "cisQTL_peaks_gr.rds")

#Extract DRACH motif from peaks
library(BSgenome.Hsapiens.UCSC.hg19)
motifs <- sampleSequence("DRACH", peaks, BSgenome.Hsapiens.UCSC.hg19)
motifs <- motifs - 2
fol <- findOverlaps(motifs, peaks)
motifs$queryHits <- NA
motifs$peakHits <- NA
motifs$queryHits[queryHits(fol)] <- queryHits(fol)
motifs$peakHits[queryHits(fol)] <- subjectHits(fol)
motifs <- motifs[!is.na(motifs$queryHits)]
mcols(motifs) <- cbind(mcols(motifs), mcols(peaks)[motifs$peakHits,])

#Find complement for the reference sequence on negative strands
SNP_gr <- GRanges(seqnames = seqnames(motifs),
                  strand = strand(motifs),
                  ranges = IRanges(start = gsub(".*:", "", motifs$SNP)))
SNP_gr$REF <- as.character(motifs$REF)
indx_neg <- as.vector(strand(SNP_gr) == "-")
SNP_gr$REF[indx_neg] <- as.character(complement(DNAStringSet(SNP_gr$REF[indx_neg])))
SNP_gr$ALT <- as.character(motifs$ALT)
SNP_gr$ALT[indx_neg] <- as.character(complement(DNAStringSet(SNP_gr$ALT[indx_neg])))

#Check if all reference sequences are matching to hg19
all( as.character(DNAStringSet( Views(Hsapiens, SNP_gr) )) == SNP_gr$REF )
motifs$DIST <- distance(motifs, SNP_gr)
indx <- motifs$DIST <=20
motifs <- motifs[indx]
motifs$SNP_gr <- SNP_gr[indx]

#Save as cisQTL_motif_gr.rds
saveRDS(motifs, "cisQTL_motif_gr.rds")
motifs <- readRDS("cisQTL_motif_gr.rds")

#liftOver to hg38
library(rtracklayer)
hg19To38 <- import.chain("hg19ToHg38.over.chain")
motifs_hg38 <- liftOver(motifs, hg19To38)
motifs_hg38 <- motifs_hg38[elementNROWS(motifs_hg38)==1]
motifs_hg38 <- unlist(motifs_hg38)
SNP_hg38 <- liftOver(motifs_hg38$SNP_gr, hg19To38)
all(elementNROWS(SNP_hg38)==1)
SNP_hg38 <- unlist(SNP_hg38)

#Check if all motifs are still DRACH and all SNP are still correct in reference
library(BSgenome.Hsapiens.UCSC.hg38)
indx1 <- vcountPattern("DRACH", DNAStringSet( Views(BSgenome.Hsapiens.UCSC.hg38, motifs_hg38 + 2) ), fixed = FALSE) > 0
indx2 <- as.character(DNAStringSet( Views(BSgenome.Hsapiens.UCSC.hg38, SNP_hg38) )) == SNP_hg38$REF
motifs_hg38$SNP_gr <- SNP_hg38
saveRDS(motifs_hg38[indx1&indx2], "cisQTL_motif_hg38_gr.rds")

## introduce SNP mutations to sequence feature
gr_zhang_m6A <- readRDS("cisQTL_motif_hg38_gr.rds")
# gr_zhang_m6A <- gr_zhang_m6A[abs(gr_zhang_m6A$DIST) <= 5,]
gr_zhang_m6A_snp <- gr_zhang_m6A$SNP_gr
gr_zhang_m6A_snp$mutateTo <- gr_zhang_m6A_snp$ALT
gr_zhang_m6A_snp$ALT <- NULL
gr_zhang_m6A_snp$REF <- NULL

## make mutations to both positive & negative samples
snp_features_zhang <- sequenceDerivedFeatures(gr_zhang_m6A + 20, BSgenome.Hsapiens.UCSC.hg38, encoding = "onehot", mutation = gr_zhang_m6A_snp)

# set.seed(3031)
# gr_zhang_m6A_snp$mutateTo <- sample(c("A","T","C","G"), length(gr_zhang_m6A_snp), replace = TRUE)
# snp_features_zhang <- sequenceDerivedFeatures(gr_zhang_m6A + 20, BSgenome.Hsapiens.UCSC.hg38, encoding = "onehot", mutation = gr_zhang_m6A_snp)

h2o.init()
mutOddsRatio <- function(mdl_obj, mut_obj){
  ref_prob <- as.vector( h2o.predict(mdl_obj, as.h2o(mut_obj$ref))$p1 )
  mut_prob <- as.vector( h2o.predict(mdl_obj, as.h2o(mut_obj$mut))$p1 )
  OR <- (mut_prob/(1-mut_prob))/(ref_prob/(1-ref_prob))
  names(OR) <- rownames(mut_obj$ref)
  return(OR)
}


## Pair up mutation indexes
sub_indx <- rownames(snp_features_zhang$ref ) ==  rownames(snp_features_zhang$mut )
snp_features_zhang$mut <- snp_features_zhang$mut[sub_indx,]
snp_features_zhang$ref <- snp_features_zhang$ref[sub_indx,]

# DRACH_index <- (snp_features_zhang$mut[,"C_19"]!= 1)&
#   (snp_features_zhang$mut[,"C_20"]!= 1 & snp_features_zhang$mut[,"T_20"]!= 1) &
#   (snp_features_zhang$mut[,"A_21"] == 1)&
#   (snp_features_zhang$mut[,"C_22"] == 1)&
#   (snp_features_zhang$mut[,"G_23"]!= 1)
# mean(DRACH_index)


# snp_features_zhang$mut <- snp_features_zhang$mut[DRACH_index,]
# snp_features_zhang$ref <- snp_features_zhang$ref[DRACH_index,]


OR_zhang_seqonly_xgboost <- mutOddsRatio(seq_only_model_xgboost, snp_features_zhang)
OR_zhang_seqonly_glm <- mutOddsRatio(seq_only_model_glm, snp_features_zhang)
OR_zhang_seqonly_rf <- mutOddsRatio(seq_only_model_rf, snp_features_zhang)
OR_zhang_seqonly_mlp <- mutOddsRatio(seq_only_model_mlp, snp_features_zhang)

# OR_zhang_seqonly_rf <- na.omit(OR_zhang_seqonly_rf)
# OR_zhang_seqonly_rf <- OR_zhang_seqonly_rf[is.finite(OR_zhang_seqonly_rf)]
range(OR_zhang_seqonly_xgboost)
range(OR_zhang_seqonly_glm)
range(OR_zhang_seqonly_rf)
range(OR_zhang_seqonly_mlp)

OR_zhang_seqonly_xgboost <- log(OR_zhang_seqonly_xgboost)
OR_zhang_seqonly_glm <- log(OR_zhang_seqonly_glm)
OR_zhang_seqonly_rf <- log(OR_zhang_seqonly_rf)
OR_zhang_seqonly_mlp <- log(OR_zhang_seqonly_mlp)

## Check in-silico m6A mutation results predicted by the sequence model
original_indx <- as.numeric(names(OR_zhang_seqonly_xgboost))
beta <- gr_zhang_m6A$beta[original_indx]


## Retrieve genome derived features
indx <- as.numeric(rownames(snp_features_zhang$ref))
genomicFeatures_gr <- readRDS(file.path(data_dir, "motif_exon_hg38_gr_gfeatures.rds"))
sub_motif_gr <- gr_zhang_m6A[indx]
fol <- findOverlaps(sub_motif_gr, genomicFeatures_gr)
sub_motif_gr <- sub_motif_gr[queryHits(fol)]
sub_snp_gr <- sub_motif_gr$SNP_gr
sub_snp_gr$mutateTo <- sub_snp_gr$ALT
gnm <- as.data.frame(mcols(genomicFeatures_gr[subjectHits(fol)]))
sub_snp_gr$ALT <- NULL
sub_snp_gr$REF <- NULL


## make mutations to both positive & negative samples
snp_features_zhang2 <- sequenceDerivedFeatures(sub_motif_gr + 20, BSgenome.Hsapiens.UCSC.hg38, mutation = sub_snp_gr, encoding = "onehot")

# set.seed(3032)
# sub_snp_gr$mutateTo <- sample(c("A","T","C","G"), length(sub_snp_gr), replace = TRUE)
# snp_features_zhang2 <- sequenceDerivedFeatures(sub_motif_gr + 20, BSgenome.Hsapiens.UCSC.hg38, mutation = sub_snp_gr, encoding = "onehot")




## Pair up mutation indexes
sub_indx <- rownames( snp_features_zhang2$ref ) ==  rownames( snp_features_zhang2$mut )
snp_features_zhang2$ref <- snp_features_zhang2$ref[sub_indx,]
snp_features_zhang2$mut <- snp_features_zhang2$mut[sub_indx,]
snp_features_zhang2$ref <- cbind(snp_features_zhang2$ref, gnm)
snp_features_zhang2$mut <- cbind(snp_features_zhang2$mut, gnm)

## Calculate mutation odds ratio
OR_zhang_addgen_xgboost <- mutOddsRatio(add_gen_model_xgboost, snp_features_zhang2)
OR_zhang_addgen_glm <- mutOddsRatio(add_gen_model_glm, snp_features_zhang2)
OR_zhang_addgen_rf <- mutOddsRatio(add_gen_model_rf, snp_features_zhang2)
OR_zhang_addgen_mlp <- mutOddsRatio(add_gen_model_mlp, snp_features_zhang2)

range(OR_zhang_addgen_xgboost)
range(OR_zhang_addgen_glm)
range(OR_zhang_addgen_rf)
range(OR_zhang_addgen_mlp)

OR_zhang_addgen_xgboost <- log(OR_zhang_addgen_xgboost)
OR_zhang_addgen_glm <- log(OR_zhang_addgen_glm)
OR_zhang_addgen_rf <- log(OR_zhang_addgen_rf)
OR_zhang_addgen_mlp <- log(OR_zhang_addgen_mlp)

## Check in-silico m6A mutation results predicted by the genome + sequence model
original_indx2 <- as.numeric(names(OR_zhang_addgen_xgboost))
beta2 <- gr_zhang_m6A$beta[original_indx2]


length(beta)
length(beta[OR_zhang_seqonly_xgboost>0 & beta>0])
length(beta[beta>0])
length(OR_zhang_seqonly_xgboost[OR_zhang_seqonly_xgboost>0])
fisher1 <- matrix(c(21, 43, 13, 37), nrow = 2) 
colnames(fisher1) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher1) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher1)

length(beta2)
length(beta2[OR_zhang_addgen_xgboost>0 & beta2>0])
length(beta2[beta2>0])
length(OR_zhang_addgen_xgboost[OR_zhang_addgen_xgboost>0])
fisher2 <- matrix(c(17, 46, 8, 41), nrow = 2) 
colnames(fisher2) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher2) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher2)


cor_test_1 <- cor.test(OR_zhang_seqonly_xgboost, beta)
cor_test_2 <- cor.test(OR_zhang_addgen_xgboost, beta2)
text1 <- paste("cor:", round(cor_test_1$estimate, 3), "p-value:", round(cor_test_1$p.value, 3))
text2 <- paste("cor:", round(cor_test_2$estimate, 3), "p-value:", round(cor_test_2$p.value, 3))

p1 <- ggplot(data.frame(OR_zhang_seqonly_xgboost, beta), aes(x = OR_zhang_seqonly_xgboost, y = beta)) + 
  geom_point(color = "green") + 
  xlab("Log odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_zhang_seqonly_xgboost), y = min(beta), label = text1, hjust = 1, vjust = -18) + 
  theme_bw()
p2 <- ggplot(data.frame(OR_zhang_addgen_xgboost, beta2), aes(x = OR_zhang_addgen_xgboost, y = beta2)) + 
  geom_point(color = "red") + 
  xlab("Log odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_zhang_addgen_xgboost), y = min(beta2), label = text2, hjust = 1, vjust = -18) + 
  theme_bw()



original_indx <- as.numeric(names(OR_zhang_seqonly_glm))
beta <- gr_zhang_m6A$beta[original_indx]
original_indx2 <- as.numeric(names(OR_zhang_addgen_glm))
beta2 <- gr_zhang_m6A$beta[original_indx2]


# Fisher's exact test
length(beta)
length(beta[OR_zhang_seqonly_glm>0 & beta>0])
length(beta[beta>0])
length(OR_zhang_seqonly_glm[OR_zhang_seqonly_glm>0])
fisher1 <- matrix(c(21, 43, 17, 33), nrow = 2) 
colnames(fisher1) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher1) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher1)

length(beta2)
length(beta2[OR_zhang_addgen_glm>0 & beta2>0])
length(beta2[beta2>0])
length(OR_zhang_addgen_glm[OR_zhang_addgen_glm>0])
fisher2 <- matrix(c(21, 42, 15, 34), nrow = 2) 
colnames(fisher2) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher2) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher2)


# Correlation test
cor_test_1 <- cor.test(OR_zhang_seqonly_glm, beta)
cor_test_2 <- cor.test(OR_zhang_addgen_glm, beta2)
text1 <- paste("cor:", round(cor_test_1$estimate, 3), "p-value:", round(cor_test_1$p.value, 3))
text2 <- paste("cor:", round(cor_test_2$estimate, 3), "p-value:", round(cor_test_2$p.value, 3))

p3 <- ggplot(data.frame(OR_zhang_seqonly_glm, beta), aes(x = OR_zhang_seqonly_glm, y = beta)) + 
  geom_point(color = "green") + 
  xlab("Log odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_zhang_seqonly_glm), y = min(beta), label = text1, hjust = 1, vjust = -18) + 
  theme_bw()
p4 <- ggplot(data.frame(OR_zhang_addgen_glm, beta2), aes(x = OR_zhang_addgen_glm, y = beta2)) + 
  geom_point(color = "red") + 
  xlab("Log odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_zhang_addgen_glm), y = min(beta2), label = text2, hjust = 1, vjust = -18) + 
  theme_bw()




original_indx <- as.numeric(names(OR_zhang_seqonly_rf))
beta <- gr_zhang_m6A$beta[original_indx]
original_indx2 <- as.numeric(names(OR_zhang_addgen_rf))
beta2 <- gr_zhang_m6A$beta[original_indx2]


length(beta)
length(beta[OR_zhang_seqonly_rf>0 & beta>0])
length(beta[beta>0])
length(OR_zhang_seqonly_rf[OR_zhang_seqonly_rf>0])
fisher1 <- matrix(c(21, 43, 13, 37), nrow = 2) 
colnames(fisher1) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher1) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher1)

length(beta2)
length(beta2[OR_zhang_addgen_rf>0 & beta2>0])
length(beta2[beta2>0])
length(OR_zhang_addgen_rf[OR_zhang_addgen_rf>0])
fisher2 <- matrix(c(17, 46, 10, 39), nrow = 2) 
colnames(fisher2) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher2) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher2)


cor_test_1 <- cor.test(OR_zhang_seqonly_rf, beta)
cor_test_2 <- cor.test(OR_zhang_addgen_rf, beta2)
text1 <- paste("cor:", round(cor_test_1$estimate, 3), "p-value:", round(cor_test_1$p.value, 3))
text2 <- paste("cor:", round(cor_test_2$estimate, 3), "p-value:", round(cor_test_2$p.value, 3))

p5 <- ggplot(data.frame(OR_zhang_seqonly_rf, beta), aes(x = OR_zhang_seqonly_rf, y = beta)) + 
  geom_point(color = "green") + 
  xlab("Log odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_zhang_seqonly_rf), y = min(beta), label = text1, hjust = 1, vjust = -18) + 
  theme_bw()

p6 <- ggplot(data.frame(OR_zhang_addgen_rf, beta2), aes(x = OR_zhang_addgen_rf, y = beta2)) + 
  geom_point(color = "red") + 
  xlab("Log odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_zhang_addgen_rf), y = min(beta2), label = text2, hjust = 1, vjust = -18) + 
  theme_bw()


original_indx <- as.numeric(names(OR_zhang_seqonly_mlp))
beta <- gr_zhang_m6A$beta[original_indx]
original_indx2 <- as.numeric(names(OR_zhang_addgen_mlp))
beta2 <- gr_zhang_m6A$beta[original_indx2]


length(beta)
length(beta[OR_zhang_seqonly_mlp>0 & beta>0])
length(beta[beta>0])
length(OR_zhang_seqonly_mlp[OR_zhang_seqonly_mlp>0])
fisher1 <- matrix(c(32, 32, 20, 30), nrow = 2) 
colnames(fisher1) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher1) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher1)

length(beta2)
length(beta2[OR_zhang_addgen_mlp>0 & beta2>0])
length(beta2[beta2>0])
length(OR_zhang_addgen_mlp[OR_zhang_addgen_mlp>0])
fisher2 <- matrix(c(26, 37, 16, 33), nrow = 2) 
colnames(fisher2) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher2) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher2)


cor_test_1 <- cor.test(OR_zhang_seqonly_mlp, beta)
cor_test_2 <- cor.test(OR_zhang_addgen_mlp, beta2)
text1 <- paste("cor:", round(cor_test_1$estimate, 3), "p-value:", round(cor_test_1$p.value, 3))
text2 <- paste("cor:", round(cor_test_2$estimate, 3), "p-value:", round(cor_test_2$p.value, 3))

p7 <- ggplot(data.frame(OR_zhang_seqonly_mlp, beta), aes(x = OR_zhang_seqonly_mlp, y = beta)) + 
  geom_point(color = "green") + 
  xlab("Log odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_zhang_seqonly_mlp), y = min(beta), label = text1, hjust = 1, vjust = -18) + 
  theme_bw()

p8 <- ggplot(data.frame(OR_zhang_addgen_mlp, beta2), aes(x = OR_zhang_addgen_mlp, y = beta2)) + 
  geom_point(color = "red") + 
  xlab("Log odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_zhang_addgen_mlp), y = min(beta2), label = text2, hjust = 1, vjust = -18) + 
  theme_bw()

cor_plot_zhang <- grid.arrange(p1, p3, p5, p7, p2, p4, p6, p8, ncol = 4)
ggsave("cor_plot_zhang.png", cor_plot_zhang, dpi = 300, width = 10, height = 6)
ggsave("cor_plot_zhang3.png", cor_plot_zhang, dpi = 300, width = 10, height = 6)

cor_plot_zhang_control <- grid.arrange(p1, p3, p5, p7, p2, p4, p6, p8, ncol = 4)
ggsave("cor_plot_zhang_control.png", cor_plot_zhang_control, dpi = 300, width = 10, height = 6)

cor_plot_zhang2 <- grid.arrange(p1, p3, p5, p7, p2, p4, p6, p8, ncol = 4)
ggsave("cor_plot_zhang2.png", cor_plot_zhang2, dpi = 300, width = 10, height = 6)








# xiong
library(readxl)
library(GenomicRanges)
library(dplyr)


# Data preprocessing
data_2 <- read_excel("41588_2021_890_MOESM2_ESM.xlsx", sheet = 2, skip = 1) %>%filter(abs(Dist) <= 20)
data_3 <- read_excel("41588_2021_890_MOESM2_ESM.xlsx", sheet = 3, skip = 1) %>%filter(abs(Dist) <= 20)
data_4 <- read_excel("41588_2021_890_MOESM2_ESM.xlsx", sheet = 4, skip = 1) %>%filter(abs(Dist) <= 20)


m6A_peak_2 <- data.frame(do.call('rbind', strsplit(as.character(data_2$`m6A Peak`), "_")))
m6A_peak_2$chromosome <- m6A_peak_2$X1
m6A_peak_2$start <- as.numeric(m6A_peak_2$X2)
m6A_peak_2$end <- as.numeric(m6A_peak_2$X3)

m6A_peak_3 <- data.frame(do.call('rbind', strsplit(as.character(data_3$`m6A Peak`), "_")))
m6A_peak_3$chromosome <- m6A_peak_3$X1
m6A_peak_3$start <- as.numeric(m6A_peak_3$X2)
m6A_peak_3$end <- as.numeric(m6A_peak_3$X3)

m6A_peak_4 <- data.frame(do.call('rbind', strsplit(as.character(data_4$`m6A Peak`), "_")))
m6A_peak_4$chromosome <- m6A_peak_4$X1
m6A_peak_4$start <- as.numeric(m6A_peak_4$X2)
m6A_peak_4$end <- as.numeric(m6A_peak_4$X3)

library(tidyr)
snp_info_2 <- select(data_2, 2)
snp_info_3 <- select(data_3, 2)
snp_info_4 <- select(data_4, 2)


snp_info_2 <- separate(snp_info_2, col = 1, into = c("chromosome", "position", "ref_base", "mut_base"), sep = ":")
snp_info_3 <- separate(snp_info_3, col = 1, into = c("chromosome", "position", "ref_base", "mut_base"), sep = ":")
snp_info_4 <- separate(snp_info_4, col = 1, into = c("chromosome", "position", "ref_base", "mut_base"), sep = ":")


snp_info_2$position <- as.numeric(snp_info_2$position)
snp_info_3$position <- as.numeric(snp_info_3$position)
snp_info_4$position <- as.numeric(snp_info_4$position)

slope_2 <- select(data_2, 5)
slope_3 <- select(data_3, 5)
slope_4 <- select(data_4, 5)

# Create Granges variables
gr_2 <- GRanges(seqnames = m6A_peak_2$chromosome,
                    ranges = IRanges(start = m6A_peak_2$start, end = m6A_peak_2$end),
                    strand = "*",
                    snp_info = snp_info_2$mut_base,
                    snp_position = snp_info_2$position,
                    snp_chromosome = snp_info_2$chromosome,
                    beta = slope_2$Slope)

gr_3 <- GRanges(seqnames = m6A_peak_3$chromosome,
                    ranges = IRanges(start = m6A_peak_3$start, end = m6A_peak_3$end),
                    strand = "*",
                    snp_info = snp_info_3$mut_base,
                    snp_position = snp_info_3$position,
                    snp_chromosome = snp_info_3$chromosome,
                    beta = slope_3$Slope)

gr_4 <- GRanges(seqnames = m6A_peak_4$chromosome,
                    ranges = IRanges(start = m6A_peak_4$start, end = m6A_peak_4$end),
                    strand = "*",
                    snp_info = snp_info_4$mut_base,
                    snp_position = snp_info_4$position,
                    snp_chromosome = snp_info_4$chromosome,
                    beta = slope_4$Slope)

gr_xiong_m6A <- c(gr_2, gr_3, gr_4)

gr_xiong_m6A_filtered <- predictiveFeatures::sampleSequence("DRACH", gr_xiong_m6A, BSgenome.Hsapiens.UCSC.hg38)
gr_xiong_m6A_filtered

gr_xiong_m6A_filtered$snp_info <- NA
gr_xiong_m6A_filtered$snp_position <- NA
gr_xiong_m6A_filtered$snp_chromosome <- NA
gr_xiong_m6A_filtered$beta <- NA

# loop over the rows in gr_xiong_m6A_filtered
for (i in seq_along(gr_xiong_m6A_filtered)) {
  # get the current row of gr_xiong_m6A_filtered
  gr_sub <- gr_xiong_m6A_filtered[i]
  
  # find overlapping regions in gr_xiong_m6A
  hits <- findOverlaps(gr_sub, gr_xiong_m6A)
  
  # get the SNP info for the overlapping regions
  snp_info <- gr_xiong_m6A[subjectHits(hits)]$snp_info
  snp_position <- gr_xiong_m6A[subjectHits(hits)]$snp_position
  snp_chromosome <- gr_xiong_m6A[subjectHits(hits)]$snp_chromosome
  beta <- gr_xiong_m6A[subjectHits(hits)]$beta
  
  # combine the SNP info into a single string and assign to the new column
  gr_xiong_m6A_filtered$snp_info[i] <- paste(snp_info, collapse = ",")
  gr_xiong_m6A_filtered$snp_position[i] <- paste(snp_position, collapse = ",")
  gr_xiong_m6A_filtered$snp_chromosome[i] <- paste(snp_chromosome, collapse = ",")
  gr_xiong_m6A_filtered$beta[i] <- paste(beta, collapse = ",")
}


gr_xiong_m6A_filtered2 <- GRanges(seqnames = seqnames(gr_xiong_m6A_filtered),
                    ranges = IRanges(start = start(gr_xiong_m6A_filtered) + 2, 
                                     end = end(gr_xiong_m6A_filtered) - 2),
                    strand = "*")




gr_xiong_m6A_filtered2$snp_info <- gr_xiong_m6A_filtered$snp_info
gr_xiong_m6A_filtered2$snp_position <- gr_xiong_m6A_filtered$snp_position
gr_xiong_m6A_filtered2$snp_chromosome <- gr_xiong_m6A_filtered$snp_chromosome
gr_xiong_m6A_filtered2$beta <- gr_xiong_m6A_filtered$beta

# saveRDS(gr_xiong_m6A_filtered2, "gr_xiong_m6A_filtered2_not_filtered.rds")
# saveRDS(gr_xiong_m6A_filtered2, "gr_xiong_m6A_filtered2.rds")
# gr_xiong_m6A_filtered2 <- readRDS("gr_xiong_m6A_filtered2_not_filtered.rds")
gr_xiong_m6A_filtered2 <- readRDS("gr_xiong_m6A_filtered2.rds")


get_first_allele <- function(x) {
  if (length(x) > 1) {
    return(x[1])
  } else {
    return(x)
  }
}


gr_xiong_m6A_filtered2$snp_info <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$snp_info), ",\\s*"), function(x) x[1])
gr_xiong_m6A_filtered2$snp_position <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$snp_position), ",\\s*"), function(x) x[1])
gr_xiong_m6A_filtered2$snp_chromosome <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$snp_chromosome), ",\\s*"), function(x) x[1])
gr_xiong_m6A_filtered2$beta <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$beta), ",\\s*"), function(x) x[1])

gr_xiong_m6A_filtered2$snp_position <- as.numeric(gr_xiong_m6A_filtered2$snp_position)
gr_xiong_m6A_snp <- GRanges(seqnames = gr_xiong_m6A_filtered2$snp_chromosome,
                ranges = IRanges(start = gr_xiong_m6A_filtered2$snp_position, end = gr_xiong_m6A_filtered2$snp_position),
                strand = "*",
                mutateTo = gr_xiong_m6A_filtered2$snp_info)

# saveRDS(gr_xiong_m6A_snp, "gr_xiong_m6A_snp.rds")

all(elementNROWS(gr_xiong_m6A_filtered2)==1)
all(elementNROWS(gr_xiong_m6A_snp)==1)

gr_xiong_m6A_filtered2$beta <- as.numeric(gr_xiong_m6A_filtered2$beta)


# install.packages('sequenceDerivedFeatures')
snp_features_xiong <- sequenceDerivedFeatures(gr_xiong_m6A_filtered2 + 20, BSgenome.Hsapiens.UCSC.hg38, encoding = "onehot", mutation = gr_xiong_m6A_snp)


# set.seed(3033)
# gr_xiong_m6A_snp$mutateTo <- sample(c("A","T","C","G"), length(gr_xiong_m6A_snp), replace = TRUE)
# snp_features_xiong <- sequenceDerivedFeatures(gr_xiong_m6A_filtered2 + 20, BSgenome.Hsapiens.UCSC.hg38, encoding = "onehot", mutation = gr_xiong_m6A_snp)



mutOddsRatio <- function(mdl_obj, mut_obj){
  ref_prob <- as.vector( h2o.predict(mdl_obj, as.h2o(mut_obj$ref))$p1 )
  mut_prob <- as.vector( h2o.predict(mdl_obj, as.h2o(mut_obj$mut))$p1 )
  OR <- (mut_prob/(1-mut_prob))/(ref_prob/(1-ref_prob))
  names(OR) <- rownames(mut_obj$ref)
  return(OR)
}

  
gr_xiong_m6A_filtered2[sub_index]
sub_index <- rownames(snp_features_xiong$ref ) ==  rownames(snp_features_xiong$mut )
snp_features_xiong$ref <- snp_features_xiong$ref[sub_index,]
snp_features_xiong$mut <- snp_features_xiong$mut[sub_index,]


h2o.init()
OR_xiong_seqonly_xgboost <- mutOddsRatio(seq_only_model_xgboost, snp_features_xiong)
OR_xiong_seqonly_glm <- mutOddsRatio(seq_only_model_glm, snp_features_xiong)
OR_xiong_seqonly_rf <- mutOddsRatio(seq_only_model_rf, snp_features_xiong)
OR_xiong_seqonly_mlp <- mutOddsRatio(seq_only_model_mlp, snp_features_xiong)
# names(OR_xiong_seqonly_xgboost)[log(OR_xiong_seqonly_xgboost) > 1]
range(OR_xiong_seqonly_xgboost)
range(OR_xiong_seqonly_glm)
range(OR_xiong_seqonly_rf)
range(OR_xiong_seqonly_mlp)


OR_xiong_seqonly_xgboost <- log(OR_xiong_seqonly_xgboost)
OR_xiong_seqonly_glm <- log(OR_xiong_seqonly_glm)
OR_xiong_seqonly_rf <- log(OR_xiong_seqonly_rf)
OR_xiong_seqonly_mlp <- log(OR_xiong_seqonly_mlp)


original_indx <- as.numeric(names(OR_xiong_seqonly_xgboost))
beta <- as.numeric(gr_xiong_m6A_filtered2$beta[original_indx])



indx <- as.numeric(rownames(snp_features_xiong$ref))
genomicFeatures_gr <- readRDS(file.path(data_dir, "motif_exon_hg38_gr_gfeatures.rds"))
sub_motif_gr <- gr_xiong_m6A_filtered2[indx]
fol <- findOverlaps(sub_motif_gr, genomicFeatures_gr)
sub_motif_gr <- sub_motif_gr[queryHits(fol)]
sub_snp_gr <- GRanges(seqnames = sub_motif_gr$snp_chromosome,
                      ranges = IRanges(start = sub_motif_gr$snp_position, 
                                         end = sub_motif_gr$snp_position),
                      strand = "*",
                      mutateTo = sub_motif_gr$snp_info            )
gnm <- as.data.frame(mcols(genomicFeatures_gr[subjectHits(fol)]))


## make mutations to both positive & negative samples
snp_features_xiong2 <- sequenceDerivedFeatures(sub_motif_gr + 20, BSgenome.Hsapiens.UCSC.hg38, mutation = sub_snp_gr, encoding = "onehot")

# set.seed(3034)
# sub_snp_gr$mutateTo <- sample(c("A","T","C","G"), length(sub_snp_gr), replace = TRUE)
# snp_features_xiong2 <- sequenceDerivedFeatures(sub_motif_gr + 20, BSgenome.Hsapiens.UCSC.hg38, mutation = sub_snp_gr, encoding = "onehot")


sub_indx <- rownames( snp_features_xiong2$ref ) ==  rownames( snp_features_xiong2$mut )
snp_features_xiong2$ref <- snp_features_xiong2$ref[sub_indx,]
snp_features_xiong2$mut <- snp_features_xiong2$mut[sub_indx,]
snp_features_xiong2$ref <- cbind(snp_features_xiong2$ref, gnm)
snp_features_xiong2$mut <- cbind(snp_features_xiong2$mut, gnm)


OR_xiong_addgen_xgboost <- mutOddsRatio(add_gen_model_xgboost, snp_features_xiong2)
OR_xiong_addgen_glm <- mutOddsRatio(add_gen_model_glm, snp_features_xiong2)
OR_xiong_addgen_rf <- mutOddsRatio(add_gen_model_rf, snp_features_xiong2)
OR_xiong_addgen_mlp <- mutOddsRatio(add_gen_model_mlp, snp_features_xiong2)

range(OR_xiong_addgen_xgboost)
range(OR_xiong_addgen_glm)
range(OR_xiong_addgen_rf)
range(OR_xiong_addgen_mlp)
OR_xiong_addgen_rf <- na.omit(OR_xiong_addgen_rf)
OR_xiong_addgen_rf <- OR_xiong_addgen_rf[is.finite(OR_xiong_addgen_rf)]

OR_xiong_addgen_xgboost <- log(OR_xiong_addgen_xgboost)
OR_xiong_addgen_glm <- log(OR_xiong_addgen_glm)
OR_xiong_addgen_rf <- log(OR_xiong_addgen_rf)
OR_xiong_addgen_mlp <- log(OR_xiong_addgen_mlp)


original_indx2 <- as.numeric(names(OR_xiong_addgen_xgboost))
beta2 <- gr_xiong_m6A_filtered2$beta[original_indx2]


length(beta)
length(beta[OR_xiong_seqonly_xgboost>0 & beta>0])
length(beta[beta>0])
length(OR_xiong_seqonly_xgboost[OR_xiong_seqonly_xgboost>0])
fisher1 <- matrix(c(1, 1, 4, 5), nrow = 2) 
colnames(fisher1) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher1) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher1)

length(beta2)
length(beta2[OR_xiong_addgen_xgboost>0 & beta2>0])
length(beta2[beta2>0])
length(OR_xiong_addgen_xgboost[OR_xiong_addgen_xgboost>0])
fisher2 <- matrix(c(1, 1, 0, 5), nrow = 2) 
colnames(fisher2) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher2) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher2)


cor_test_1 <- cor.test(OR_xiong_seqonly_xgboost, beta)
cor_test_2 <- cor.test(OR_xiong_addgen_xgboost, beta2)
text1 <- paste("cor:", round(cor_test_1$estimate, 3), "p-value:", round(cor_test_1$p.value, 3))
text2 <- paste("cor:", round(cor_test_2$estimate, 3), "p-value:", round(cor_test_2$p.value, 3))

p1 <- ggplot(data.frame(OR_xiong_seqonly_xgboost, beta), aes(x = OR_xiong_seqonly_xgboost, y = beta)) + 
  geom_point(color = "green") + 
  xlab("Odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_xiong_seqonly_xgboost), y = min(beta), label = text1, hjust = 1, vjust = -18) + 
  theme_bw()
p2 <- ggplot(data.frame(OR_xiong_addgen_xgboost, beta2), aes(x = OR_xiong_addgen_xgboost, y = beta2)) + 
  geom_point(color = "red") + 
  xlab("Odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_xiong_addgen_xgboost), y = min(beta2), label = text2, hjust = 1, vjust = -18) + 
  theme_bw()



original_indx <- as.numeric(names(OR_xiong_seqonly_glm))
beta <- gr_xiong_m6A_filtered2$beta[original_indx]
original_indx2 <- as.numeric(names(OR_xiong_addgen_glm))
beta2 <- gr_xiong_m6A_filtered2$beta[original_indx2]


length(beta)
length(beta[OR_xiong_seqonly_glm>0 & beta>0])
length(beta[beta>0])
length(OR_xiong_seqonly_glm[OR_xiong_seqonly_glm>0])
fisher1 <- matrix(c(2, 0, 2, 7), nrow = 2) 
colnames(fisher1) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher1) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher1)

length(beta2)
length(beta2[OR_xiong_addgen_glm>0 & beta2>0])
length(beta2[beta2>0])
length(OR_xiong_addgen_glm[OR_xiong_addgen_glm>0])
fisher2 <- matrix(c(1, 1, 1, 4), nrow = 2) 
colnames(fisher2) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher2) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher2)


cor_test_1 <- cor.test(OR_xiong_seqonly_glm, beta)
cor_test_2 <- cor.test(OR_xiong_addgen_glm, beta2)
text1 <- paste("cor:", round(cor_test_1$estimate, 3), "p-value:", round(cor_test_1$p.value, 3))
text2 <- paste("cor:", round(cor_test_2$estimate, 3), "p-value:", round(cor_test_2$p.value, 3))

p3 <- ggplot(data.frame(OR_xiong_seqonly_glm, beta), aes(x = OR_xiong_seqonly_glm, y = beta)) + 
  geom_point(color = "green") + 
  xlab("Odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_xiong_seqonly_glm), y = min(beta), label = text1, hjust = 1, vjust = -18) + 
  theme_bw()

p4 <- ggplot(data.frame(OR_xiong_addgen_glm, beta2), aes(x = OR_xiong_addgen_glm, y = beta2)) + 
  geom_point(color = "red") + 
  xlab("Odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_xiong_addgen_glm), y = min(beta2), label = text2, hjust = 1, vjust = -18) + 
  theme_bw()




original_indx <- as.numeric(names(OR_xiong_seqonly_rf))
beta <- gr_xiong_m6A_filtered2$beta[original_indx]
original_indx2 <- as.numeric(names(OR_xiong_addgen_rf))
beta2 <- gr_xiong_m6A_filtered2$beta[original_indx2]


length(beta)
length(beta[OR_xiong_seqonly_rf>0 & beta>0])
length(beta[beta>0])
length(OR_xiong_seqonly_rf[OR_xiong_seqonly_rf>0])
fisher1 <- matrix(c(1, 1, 5, 1), nrow = 2) 
colnames(fisher1) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher1) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher1)

length(beta2)
length(beta2[OR_xiong_addgen_glm>0 & beta2>0])
length(beta2[beta2>0])
length(OR_xiong_addgen_glm[OR_xiong_addgen_glm>0])
fisher2 <- matrix(c(1, 1, 1, 4), nrow = 2) 
colnames(fisher2) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher2) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher2)


cor_test_1 <- cor.test(OR_xiong_seqonly_rf, beta)
cor_test_2 <- cor.test(OR_xiong_addgen_rf, beta2)
text1 <- paste("cor:", round(cor_test_1$estimate, 3), "p-value:", round(cor_test_1$p.value, 3))
text2 <- paste("cor:", round(cor_test_2$estimate, 3), "p-value:", round(cor_test_2$p.value, 3))

p5 <- ggplot(data.frame(OR_xiong_seqonly_rf, beta), aes(x = OR_xiong_seqonly_rf, y = beta)) + 
  geom_point(color = "green") + 
  xlab("Odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_xiong_seqonly_rf), y = min(beta), label = text1, hjust = 1, vjust = -18) + 
  theme_bw()
p6 <- ggplot(data.frame(OR_xiong_addgen_rf, beta2), aes(x = OR_xiong_addgen_rf, y = beta2)) + 
  geom_point(color = "red") + 
  xlab("Odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_xiong_addgen_rf), y = min(beta2), label = text2, hjust = 1, vjust = -18) + 
  theme_bw()



original_indx <- as.numeric(names(OR_xiong_seqonly_mlp))
beta <- gr_xiong_m6A_filtered2$beta[original_indx]
original_indx2 <- as.numeric(names(OR_xiong_addgen_mlp))
beta2 <- gr_xiong_m6A_filtered2$beta[original_indx2]


length(beta)
length(beta[OR_xiong_seqonly_mlp>0 & beta>0])
length(beta[beta>0])
length(OR_xiong_seqonly_mlp[OR_xiong_seqonly_mlp>0])
fisher1 <- matrix(c(1, 1, 4, 1), nrow = 2) 
colnames(fisher1) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher1) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher1)

length(beta2)
length(beta2[OR_xiong_addgen_mlp>0 & beta2>0])
length(beta2[beta2>0])
length(OR_xiong_addgen_mlp[OR_xiong_addgen_mlp>0])
fisher2 <- matrix(c(1, 1, 2, 3), nrow = 2) 
colnames(fisher2) <- c("m6AQTL gain", "m6AQTL loss") 
rownames(fisher2) <- c("in-silico mutation gain", "in-silico mutation loss")
fisher.test(fisher2)


cor_test_1 <- cor.test(OR_xiong_seqonly_mlp, beta)
cor_test_2 <- cor.test(OR_xiong_addgen_mlp, beta2)
text1 <- paste("cor:", round(cor_test_1$estimate, 3), "p-value:", round(cor_test_1$p.value, 3))
text2 <- paste("cor:", round(cor_test_2$estimate, 3), "p-value:", round(cor_test_2$p.value, 3))

p7 <- ggplot(data.frame(OR_xiong_seqonly_mlp, beta), aes(x = OR_xiong_seqonly_mlp, y = beta)) + 
  geom_point(color = "green") + 
  xlab("Odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_xiong_seqonly_mlp), y = min(beta), label = text1, hjust = 1, vjust = -18) + 
  theme_bw()

p8 <- ggplot(data.frame(OR_xiong_addgen_mlp, beta2), aes(x = OR_xiong_addgen_mlp, y = beta2)) + 
  geom_point(color = "red") + 
  xlab("Odds ratio of ML inference") + 
  ylab("Beta of m6AQTL") + 
  annotate("text", x = max(OR_xiong_addgen_mlp), y = min(beta2), label = text2, hjust = 1, vjust = -18) + 
  theme_bw()

cor_plot_xiong <- grid.arrange(p1, p3, p5, p7, p2, p4, p6, p8, ncol = 4)
ggsave("cor_plot_xiong.png", cor_plot_xiong, dpi = 300, width = 10, height = 6)

cor_plot_xiong_control <- grid.arrange(p1, p3, p5, p7, p2, p4, p6, p8, ncol = 4)
ggsave("cor_plot_xiong_control.png", cor_plot_xiong_control, dpi = 300, width = 10, height = 6)









# RMDisease2/Zhang
RMDisease <- read.csv("m6a_human_associatedSNPs_hg19.csv")
RMDisease_gr <- GRanges(seqnames = RMDisease$seqnames,
                  ranges = IRanges(start = RMDisease$MD_ChromStart, end = RMDisease$MD_ChromEnd),
                  strand = RMDisease$MD_Strand,
                  SNP_ChromStart = RMDisease$SNP_ChromStart,
                  SNP_ChromEnd = RMDisease$SNP_ChromEnd,
                  SNP_Strand = RMDisease$SNP_Strand,
                  ref = RMDisease$ref,
                  alt = RMDisease$alt,
                  AL = RMDisease$AL)


RMDisease_gr_motif <- predictiveFeatures::sampleSequence("DRACH", RMDisease_gr, BSgenome.Hsapiens.UCSC.hg19)

RMDisease_gr_motif$SNP_ChromStart <- NA
RMDisease_gr_motif$SNP_ChromEnd <- NA
RMDisease_gr_motif$SNP_Strand <-NA
RMDisease_gr_motif$ref <- NA
RMDisease_gr_motif$alt <- NA
RMDisease_gr_motif$AL <- NA

for (i in seq_along(RMDisease_gr_motif)) {
  
  gr_sub <- RMDisease_gr_motif[i]

  hits <- findOverlaps(gr_sub, RMDisease_gr)
  
  SNP_ChromStart <- RMDisease_gr[subjectHits(hits)]$SNP_ChromStart
  SNP_ChromEnd <- RMDisease_gr[subjectHits(hits)]$SNP_ChromEnd
  ref <- RMDisease_gr[subjectHits(hits)]$ref
  alt <- RMDisease_gr[subjectHits(hits)]$alt
  
  RMDisease_gr_motif$SNP_ChromStart[i] <- paste(SNP_ChromStart, collapse = ",")
  RMDisease_gr_motif$SNP_ChromEnd[i] <- paste(SNP_ChromEnd, collapse = ",")
  RMDisease_gr_motif$ref[i] <- paste(ref, collapse = ",")
  RMDisease_gr_motif$alt[i] <- paste(alt, collapse = ",")
}

for (i in seq_along(RMDisease_gr_motif)) {
  
  gr_sub <- RMDisease_gr_motif[i]
  
  hits <- findOverlaps(gr_sub, RMDisease_gr)
  
  SNP_Strand <- RMDisease_gr[subjectHits(hits)]$SNP_Strand
  
  RMDisease_gr_motif$SNP_Strand[i] <- paste(SNP_Strand, collapse = ",")

}

for (i in seq_along(RMDisease_gr_motif)) {
  
  gr_sub <- RMDisease_gr_motif[i]
  
  hits <- findOverlaps(gr_sub, RMDisease_gr)
  
  AL <- RMDisease_gr[subjectHits(hits)]$AL
  
  RMDisease_gr_motif$AL[i] <- paste(AL, collapse = ",")
  
}

# saveRDS(RMDisease_gr_motif, "RMDisease_gr_motif.rds")
RMDisease_gr_motif <- readRDS("RMDisease_gr_motif.rds")

get_first_allele <- function(x) {
  if (length(x) > 1) {
    return(x[1])
  } else {
    return(x)
  }
}
RMDisease_gr_motif$SNP_ChromStart <- sapply(strsplit(as.character(RMDisease_gr_motif$SNP_ChromStart), ",\\s*"), function(x) x[1])
RMDisease_gr_motif$SNP_ChromEnd <- sapply(strsplit(as.character(RMDisease_gr_motif$SNP_ChromEnd), ",\\s*"), function(x) x[1])
RMDisease_gr_motif$ref <- sapply(strsplit(as.character(RMDisease_gr_motif$ref), ",\\s*"), function(x) x[1])
RMDisease_gr_motif$alt <- sapply(strsplit(as.character(RMDisease_gr_motif$alt), ",\\s*"), function(x) x[1])
RMDisease_gr_motif$SNP_Strand <- sapply(strsplit(as.character(RMDisease_gr_motif$SNP_Strand), ",\\s*"), function(x) x[1])
RMDisease_gr_motif$AL <- sapply(strsplit(as.character(RMDisease_gr_motif$AL), ",\\s*"), function(x) x[1])
RMDisease_gr_motif <- RMDisease_gr_motif - 2


hg19To38 <- import.chain("hg19ToHg38.over.chain")
RMDisease_motif_hg38 <- liftOver(RMDisease_gr_motif, hg19To38)
RMDisease_motif_hg38 <- unlist(RMDisease_motif_hg38)
all(elementNROWS(RMDisease_motif_hg38)==1)
# saveRDS(RMDisease_motif_hg38, "RMDisease_motif_hg38.rds")
RMDisease_motif_hg38 <- readRDS("RMDisease_motif_hg38.rds")



gr_zhang_m6A <- readRDS("cisQTL_motif_hg38_gr.rds")
RMDisease_motif_hg38 <- readRDS("RMDisease_motif_hg38.rds")

RMDisease_motif_hg38$SNP_ChromStart <- as.numeric(RMDisease_motif_hg38$SNP_ChromStart)
RMDisease_motif_hg38$SNP_ChromEnd <- as.numeric(RMDisease_motif_hg38$SNP_ChromEnd)
RMDisease_snp <- GRanges(seqnames = seqnames(RMDisease_motif_hg38),
                         ranges = IRanges(start = RMDisease_motif_hg38$SNP_ChromStart, end = RMDisease_motif_hg38$SNP_ChromEnd),
                         strand = RMDisease_motif_hg38$SNP_Strand,
                         ref = RMDisease_motif_hg38$ref,
                         alt = RMDisease_motif_hg38$alt,
                         AL = RMDisease_motif_hg38$AL)
RMDisease_motif_hg38$SNP_gr <- RMDisease_snp

overlaps <- findOverlaps(gr_zhang_m6A, RMDisease_motif_hg38)
RMDisease_motif_hg38 <- RMDisease_motif_hg38[subjectHits(overlaps)]
gr_zhang_m6A <- gr_zhang_m6A[queryHits(overlaps)]

RMDisease_snp <- RMDisease_motif_hg38$SNP_gr
hg19To38 <- import.chain("hg19ToHg38.over.chain")
RMDisease_snp_hg38 <- liftOver(RMDisease_snp, hg19To38)
RMDisease_snp_hg38 <- unlist(RMDisease_snp_hg38)
all(elementNROWS(RMDisease_snp_hg38)==1)

gr_zhang_m6A_snp <- gr_zhang_m6A$SNP_gr
gr_zhang_m6A_snp$beta <- gr_zhang_m6A$beta

overlaps2 <- findOverlaps(gr_zhang_m6A_snp, RMDisease_snp_hg38)
RMDisease_snp_hg38 <- RMDisease_snp_hg38[subjectHits(overlaps2)]
gr_zhang_m6A_snp <- gr_zhang_m6A_snp[queryHits(overlaps2)]

RMDisease_snp_hg38$AL <- as.numeric(RMDisease_snp_hg38$AL)
cor.test(gr_zhang_m6A_snp$beta, RMDisease_snp_hg38$AL)





# RMDisease2/Xiong
RMDisease_motif_hg38 <- readRDS("RMDisease_motif_hg38.rds")
RMDisease_motif_hg38$SNP_ChromStart <- as.numeric(RMDisease_motif_hg38$SNP_ChromStart)
RMDisease_motif_hg38$SNP_ChromEnd <- as.numeric(RMDisease_motif_hg38$SNP_ChromEnd)
RMDisease_snp <- GRanges(seqnames = seqnames(RMDisease_motif_hg38),
                         ranges = IRanges(start = RMDisease_motif_hg38$SNP_ChromStart, end = RMDisease_motif_hg38$SNP_ChromEnd),
                         strand = RMDisease_motif_hg38$SNP_Strand,
                         ref = RMDisease_motif_hg38$ref,
                         alt = RMDisease_motif_hg38$alt,
                         AL = RMDisease_motif_hg38$AL)
RMDisease_motif_hg38$SNP_gr <- RMDisease_snp

gr_xiong_m6A_filtered2 <- readRDS("gr_xiong_m6A_filtered2.rds")
get_first_allele <- function(x) {
  if (length(x) > 1) {
    return(x[1])
  } else {
    return(x)
  }
}
gr_xiong_m6A_filtered2$snp_info <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$snp_info), ",\\s*"), function(x) x[1])
gr_xiong_m6A_filtered2$snp_position <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$snp_position), ",\\s*"), function(x) x[1])
gr_xiong_m6A_filtered2$snp_chromosome <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$snp_chromosome), ",\\s*"), function(x) x[1])
gr_xiong_m6A_filtered2$beta <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$beta), ",\\s*"), function(x) x[1])

gr_xiong_m6A_filtered2$snp_position <- as.numeric(gr_xiong_m6A_filtered2$snp_position)
gr_xiong_m6A_snp <- GRanges(seqnames = gr_xiong_m6A_filtered2$snp_chromosome,
                            ranges = IRanges(start = gr_xiong_m6A_filtered2$snp_position, end = gr_xiong_m6A_filtered2$snp_position),
                            strand = "*",
                            mutateTo = gr_xiong_m6A_filtered2$snp_info)

gr_xiong_m6A_filtered2$beta <- as.numeric(gr_xiong_m6A_filtered2$beta)
gr_xiong_m6A_filtered2$SNP_gr <- gr_xiong_m6A_snp


overlaps <- findOverlaps(gr_xiong_m6A_filtered2, RMDisease_motif_hg38)
RMDisease_motif_hg38 <- RMDisease_motif_hg38[subjectHits(overlaps)]
gr_xiong_m6A_filtered2 <- gr_xiong_m6A_filtered2[queryHits(overlaps)]
unique(RMDisease_motif_hg38)

RMDisease_snp <- RMDisease_motif_hg38$SNP_gr
hg19To38 <- import.chain("hg19ToHg38.over.chain")
RMDisease_snp_hg38 <- liftOver(RMDisease_snp, hg19To38)
RMDisease_snp_hg38 <- unlist(RMDisease_snp_hg38)
all(elementNROWS(RMDisease_snp_hg38)==1)

gr_xiong_m6A_snp <- gr_xiong_m6A_filtered2$SNP_gr

overlaps2 <- findOverlaps(gr_xiong_m6A_snp, RMDisease_snp) #无




# RMVar/Zhang

library(readr)
library(openxlsx)


# RMVar <- read_delim("RMVar_Human_basic_info_m6A_hg38.txt", delim = "\t")
# RMVar <- RMVar[!RMVar$sample_peak_id == '-', ]
# write_csv(RMVar, "RMVar.csv")
RMVar <- read_csv("RMVar.csv")

sample_peak_id <- data.frame(do.call('rbind', strsplit(RMVar$sample_peak_id, "[:.,]")))
RMVar$peak_chr <- sample_peak_id$X1
RMVar$peak_start_position <- as.numeric(sample_peak_id$X2)
RMVar$peak_end_position <- as.numeric(sample_peak_id$X4)
RMVar <- RMVar[!RMVar$peak_chr == 'chrMT',]

RMVar_gr <- GRanges(seqnames = RMVar$peak_chr,
                    ranges = IRanges(start = RMVar$peak_start_position, end = RMVar$peak_end_position),
                    strand = RMVar$strand,
                    snp_start = RMVar$snp_start,
                    snp_end = RMVar$snp_end,
                    ref = RMVar$reference_base,
                    alt = RMVar$alterative_base,
                    mf = RMVar$modification_function)

RMVar_gr_motif <- predictiveFeatures::sampleSequence("DRACH", RMVar_gr, BSgenome.Hsapiens.UCSC.hg38)
RMVar_gr_motif <- RMVar_gr_motif - 2

gr_zhang_m6A <- readRDS("cisQTL_motif_hg38_gr.rds")
overlaps <- findOverlaps(gr_zhang_m6A, RMVar_gr_motif)
sub_RMVar_gr_motif <- RMVar_gr_motif[subjectHits(overlaps)]
sub_gr_zhang_m6A <- gr_zhang_m6A[queryHits(overlaps)]
unique(sub_RMVar_gr_motif)

sub_RMVar_gr_motif$snp_start <- NA
sub_RMVar_gr_motif$snp_end <- NA
sub_RMVar_gr_motif$reference_base <- NA
sub_RMVar_gr_motif$alterative_base <- NA
sub_RMVar_gr_motif$mf <- NA


for (i in seq_along(sub_RMVar_gr_motif)) {
  
  gr_sub <- sub_RMVar_gr_motif[i]
  
  hits <- findOverlaps(gr_sub, RMVar_gr)
  
  snp_start <- RMVar_gr[subjectHits(hits)]$snp_start
  snp_end <- RMVar_gr[subjectHits(hits)]$snp_end
  reference_base <- RMVar_gr[subjectHits(hits)]$ref
  alterative_base <- RMVar_gr[subjectHits(hits)]$alt
  mf <- RMVar_gr[subjectHits(hits)]$mf
  
  sub_RMVar_gr_motif$snp_start[i] <- paste(snp_start, collapse = ",")
  sub_RMVar_gr_motif$snp_end[i] <- paste(snp_end, collapse = ",")
  sub_RMVar_gr_motif$reference_base[i] <- paste(reference_base, collapse = ",")
  sub_RMVar_gr_motif$alterative_base[i] <- paste(alterative_base, collapse = ",")
  sub_RMVar_gr_motif$mf[i] <- paste(mf, collapse = ",")
  
}


# saveRDS(sub_RMVar_gr_motif, "sub_RMVar.rds")
sub_RMVar <- readRDS("sub_RMVar.rds")

get_first_allele <- function(x) {
  if (length(x) > 1) {
    return(x[1])
  } else {
    return(x)
  }
}

sub_RMVar$snp_start <- sapply(strsplit(as.character(sub_RMVar$snp_start), ",\\s*"), function(x) x[1])
sub_RMVar$snp_end <- sapply(strsplit(as.character(sub_RMVar$snp_end), ",\\s*"), function(x) x[1])
sub_RMVar$reference_base <- sapply(strsplit(as.character(sub_RMVar$reference_base), ",\\s*"), function(x) x[1])
sub_RMVar$alterative_base <- sapply(strsplit(as.character(sub_RMVar$alterative_base), ",\\s*"), function(x) x[1])
sub_RMVar$mf <- sapply(strsplit(as.character(sub_RMVar$mf), ",\\s*"), function(x) x[1])

sub_RMVar$reference_base <- substr(sub_RMVar$reference_base, start = 1, stop = 1)
sub_RMVar$alterative_base <- substr(sub_RMVar$alterative_base, start = 1, stop = 1)
sub_RMVar$snp_start <- as.numeric(sub_RMVar$snp_start)
sub_RMVar$snp_end <- as.numeric(sub_RMVar$snp_end)

RMVar_snp <- GRanges(seqnames = seqnames(sub_RMVar),
                     ranges = IRanges(start = sub_RMVar$snp_start, end = sub_RMVar$snp_end),
                     strand = "*",
                     mutateTo = sub_RMVar$alterative_base)

RMVar_snp$m6A_site <- sub_RMVar
end(RMVar_snp[width(RMVar_snp) != 1]) <- start(RMVar_snp[width(RMVar_snp) != 1])
width(RMVar_snp) <- 1


gr_zhang_m6A_snp <- sub_gr_zhang_m6A$SNP_gr
gr_zhang_m6A_snp$beta <- sub_gr_zhang_m6A$beta
gr_zhang_m6A_snp$m6A_site <- sub_gr_zhang_m6A


overlaps2 <- findOverlaps(gr_zhang_m6A_snp, RMVar_snp)
sub_RMVar_snp <- RMVar_snp[subjectHits(overlaps2)]
sub_gr_zhang_m6A_snp <- gr_zhang_m6A_snp[queryHits(overlaps2)]






# RMVar/Xiong
RMVar <- read_csv("RMVar.csv")
sample_peak_id <- data.frame(do.call('rbind', strsplit(RMVar$sample_peak_id, "[:.,]")))
RMVar$peak_chr <- sample_peak_id$X1
RMVar$peak_start_position <- as.numeric(sample_peak_id$X2)
RMVar$peak_end_position <- as.numeric(sample_peak_id$X4)
RMVar <- RMVar[!RMVar$peak_chr == 'chrMT',]
RMVar_gr <- GRanges(seqnames = RMVar$peak_chr,
                    ranges = IRanges(start = RMVar$peak_start_position, end = RMVar$peak_end_position),
                    strand = RMVar$strand,
                    snp_start = RMVar$snp_start,
                    snp_end = RMVar$snp_end,
                    ref = RMVar$reference_base,
                    alt = RMVar$alterative_base,
                    mf = RMVar$modification_function)

RMVar_gr_motif <- predictiveFeatures::sampleSequence("DRACH", RMVar_gr, BSgenome.Hsapiens.UCSC.hg38)
RMVar_gr_motif <- RMVar_gr_motif - 2

gr_xiong_m6A_filtered2 <- readRDS("gr_xiong_m6A_filtered2.rds")
get_first_allele <- function(x) {
  if (length(x) > 1) {
    return(x[1])
  } else {
    return(x)
  }
}
gr_xiong_m6A_filtered2$snp_info <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$snp_info), ",\\s*"), function(x) x[1])
gr_xiong_m6A_filtered2$snp_position <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$snp_position), ",\\s*"), function(x) x[1])
gr_xiong_m6A_filtered2$snp_chromosome <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$snp_chromosome), ",\\s*"), function(x) x[1])
gr_xiong_m6A_filtered2$beta <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$beta), ",\\s*"), function(x) x[1])

gr_xiong_m6A_filtered2$snp_position <- as.numeric(gr_xiong_m6A_filtered2$snp_position)
gr_xiong_m6A_snp <- GRanges(seqnames = gr_xiong_m6A_filtered2$snp_chromosome,
                            ranges = IRanges(start = gr_xiong_m6A_filtered2$snp_position, end = gr_xiong_m6A_filtered2$snp_position),
                            strand = "*",
                            mutateTo = gr_xiong_m6A_filtered2$snp_info)

gr_xiong_m6A_filtered2$beta <- as.numeric(gr_xiong_m6A_filtered2$beta)
gr_xiong_m6A_filtered2$SNP_gr <- gr_xiong_m6A_snp


overlaps <- findOverlaps(gr_xiong_m6A_filtered2, RMVar_gr_motif)
sub_RMVar_gr_motif2 <- RMVar_gr_motif[subjectHits(overlaps)]
gr_xiong_m6A_filtered2 <- gr_xiong_m6A_filtered2[queryHits(overlaps)]
unique(sub_RMVar_gr_motif2)


sub_RMVar_gr_motif2$snp_start <- NA
sub_RMVar_gr_motif2$snp_end <- NA
sub_RMVar_gr_motif2$reference_base <- NA
sub_RMVar_gr_motif2$alterative_base <- NA
sub_RMVar_gr_motif2$mf <- NA


for (i in seq_along(sub_RMVar_gr_motif2)) {
  
  gr_sub <- sub_RMVar_gr_motif2[i]
  
  hits <- findOverlaps(gr_sub, RMVar_gr)
  
  snp_start <- RMVar_gr[subjectHits(hits)]$snp_start
  snp_end <- RMVar_gr[subjectHits(hits)]$snp_end
  reference_base <- RMVar_gr[subjectHits(hits)]$ref
  alterative_base <- RMVar_gr[subjectHits(hits)]$alt
  mf <- RMVar_gr[subjectHits(hits)]$mf
  
  sub_RMVar_gr_motif2$snp_start[i] <- paste(snp_start, collapse = ",")
  sub_RMVar_gr_motif2$snp_end[i] <- paste(snp_end, collapse = ",")
  sub_RMVar_gr_motif2$reference_base[i] <- paste(reference_base, collapse = ",")
  sub_RMVar_gr_motif2$alterative_base[i] <- paste(alterative_base, collapse = ",")
  sub_RMVar_gr_motif2$mf[i] <- paste(mf, collapse = ",")
  
}


# saveRDS(sub_RMVar_gr_motif2, "sub_RMVar2.rds")
sub_RMVar2 <- readRDS("sub_RMVar2.rds")

get_first_allele <- function(x) {
  if (length(x) > 1) {
    return(x[1])
  } else {
    return(x)
  }
}

sub_RMVar2$snp_start <- sapply(strsplit(as.character(sub_RMVar2$snp_start), ",\\s*"), function(x) x[1])
sub_RMVar2$snp_end <- sapply(strsplit(as.character(sub_RMVar2$snp_end), ",\\s*"), function(x) x[1])
sub_RMVar2$reference_base <- sapply(strsplit(as.character(sub_RMVar2$reference_base), ",\\s*"), function(x) x[1])
sub_RMVar2$alterative_base <- sapply(strsplit(as.character(sub_RMVar2$alterative_base), ",\\s*"), function(x) x[1])
sub_RMVar2$mf <- sapply(strsplit(as.character(sub_RMVar2$mf), ",\\s*"), function(x) x[1])

sub_RMVar2$reference_base <- substr(sub_RMVar2$reference_base, start = 1, stop = 1)
sub_RMVar2$alterative_base <- substr(sub_RMVar2$alterative_base, start = 1, stop = 1)
sub_RMVar2$snp_start <- as.numeric(sub_RMVar2$snp_start)
sub_RMVar2$snp_end <- as.numeric(sub_RMVar2$snp_end)

RMVar_snp2 <- GRanges(seqnames = seqnames(sub_RMVar2),
                      ranges = IRanges(start = sub_RMVar2$snp_start, end = sub_RMVar2$snp_end),
                      strand = "*",
                      mutateTo = sub_RMVar2$alterative_base)

RMVar_snp2$m6A_site <- sub_RMVar2
end(RMVar_snp2[width(RMVar_snp2) != 1]) <- start(RMVar_snp2[width(RMVar_snp2) != 1])
width(RMVar_snp2) <- 1

gr_xiong_snp <- gr_xiong_m6A_filtered2$SNP_gr
gr_xiong_snp$m6A_site <- gr_xiong_m6A_filtered2

overlaps2 <- findOverlaps(gr_xiong_snp, RMVar_snp2) #无
RMVar_snp2 <- RMVar_snp2[subjectHits(overlaps)]
gr_xiong_snp <- gr_xiong_snp[queryHits(overlaps)] 








# m6AVar/Zhang

# m6AVar <- read_delim("Human_protein_coding.txt", delim = "\t")
# write_csv(m6AVar, "m6AVar2.csv")
# m6AVar <- read_csv("m6AVar2.csv")
# m6AVar <- m6AVar[!is.na(m6AVar$Sample_PeakID), ]
# write_csv(m6AVar, "m6AVar.csv")
m6AVar <- read_csv("m6AVar.csv")

Sample_PeakID <- data.frame(do.call('rbind', strsplit(m6AVar$Sample_PeakID, "[:.;]")))
peak_chr <- Sample_PeakID$X1
peak_start_position <- as.numeric(Sample_PeakID$X2)
peak_end_position <- as.numeric(Sample_PeakID$X4)
m6AVar$Start_position <- as.numeric(m6AVar$Start_position)
m6AVar$End_position <- as.numeric(m6AVar$End_position)

m6AVar_gr <- GRanges(seqnames = peak_chr,
                     ranges = IRanges(start = peak_start_position, end = peak_end_position),
                     strand = m6AVar$Strand,
                     SNP_Chromosome = m6AVar$SNP_Chromosome,
                     Start_position = m6AVar$Start_position,
                     End_position = m6AVar$End_position,
                     Reference_Base = m6AVar$Reference_Base,
                     Alterative_Base = m6AVar$Alterative_Base,
                     m6A_Function = m6AVar$m6A_Function)


m6AVar_gr_motif <- predictiveFeatures::sampleSequence("DRACH", m6AVar_gr, BSgenome.Hsapiens.UCSC.hg19)
m6AVar_gr_motif <- m6AVar_gr_motif - 2

motifs <- readRDS("cisQTL_motif_gr.rds")
overlaps <- findOverlaps(motifs, m6AVar_gr_motif)
sub_gr_zhang_m6A <- motifs[queryHits(overlaps)]
sub_m6AVar_motif <- m6AVar_gr_motif[subjectHits(overlaps)]
unique(sub_gr_zhang_m6A)



sub_m6AVar_motif$SNP_Chromosome <- NA
sub_m6AVar_motif$Start_position <- NA
sub_m6AVar_motif$End_position <- NA
sub_m6AVar_motif$Reference_Base <- NA
sub_m6AVar_motif$Alterative_Base <- NA
sub_m6AVar_motif$m6A_Function <- NA


for (i in seq_along(sub_m6AVar_motif)) {
  
  gr_sub <- sub_m6AVar_motif[i]
  
  hits <- findOverlaps(gr_sub, m6AVar_gr)
  
  SNP_Chromosome <- m6AVar_gr[subjectHits(hits)]$SNP_Chromosome
  Start_position <- m6AVar_gr[subjectHits(hits)]$Start_position
  End_position <- m6AVar_gr[subjectHits(hits)]$End_position
  Reference_Base <- m6AVar_gr[subjectHits(hits)]$Reference_Base
  Alterative_Base <- m6AVar_gr[subjectHits(hits)]$Alterative_Base
  m6A_Function <- m6AVar_gr[subjectHits(hits)]$m6A_Function
  
  sub_m6AVar_motif$SNP_Chromosome[i] <- paste(SNP_Chromosome, collapse = ",")
  sub_m6AVar_motif$Start_position[i] <- paste(Start_position, collapse = ",")
  sub_m6AVar_motif$End_position[i] <- paste(End_position, collapse = ",")
  sub_m6AVar_motif$Reference_Base[i] <- paste(Reference_Base, collapse = ",")
  sub_m6AVar_motif$Alterative_Base[i] <- paste(Alterative_Base, collapse = ",")
  sub_m6AVar_motif$m6A_Function[i] <- paste(m6A_Function, collapse = ",")
  
}


# saveRDS(sub_m6AVar_motif, "sub_m6AVar_motif.rds")
sub_m6AVar_motif <- readRDS("sub_m6AVar_motif.rds")

get_first_allele <- function(x) {
  if (length(x) > 1) {
    return(x[1])
  } else {
    return(x)
  }
}

sub_m6AVar_motif$SNP_Chromosome <- sapply(strsplit(as.character(sub_m6AVar_motif$SNP_Chromosome), ",\\s*"), function(x) x[1])
sub_m6AVar_motif$Start_position <- sapply(strsplit(as.character(sub_m6AVar_motif$Start_position), ",\\s*"), function(x) x[1])
sub_m6AVar_motif$End_position <- sapply(strsplit(as.character(sub_m6AVar_motif$End_position), ",\\s*"), function(x) x[1])
sub_m6AVar_motif$Reference_Base <- sapply(strsplit(as.character(sub_m6AVar_motif$Reference_Base), ",\\s*"), function(x) x[1])
sub_m6AVar_motif$Alterative_Base <- sapply(strsplit(as.character(sub_m6AVar_motif$Alterative_Base), ",\\s*"), function(x) x[1])
sub_m6AVar_motif$m6A_Function <- sapply(strsplit(as.character(sub_m6AVar_motif$m6A_Function), ",\\s*"), function(x) x[1])

hg19To38 <- import.chain("hg19ToHg38.over.chain")
m6AVar_hg38 <- liftOver(sub_m6AVar_motif, hg19To38)
m6AVar_hg38 <- unlist(m6AVar_hg38)
all(elementNROWS(m6AVar_hg38)==1)

sub_zhang_hg38 <- liftOver(sub_gr_zhang_m6A , hg19To38)
sub_zhang_hg38 <- unlist(sub_zhang_hg38 )
all(elementNROWS(sub_zhang_hg38 )==1)

# saveRDS(m6AVar_hg38, "m6AVar_hg38.rds")
m6AVar_hg38 <- readRDS("m6AVar_hg38.rds")

m6AVar_hg38$Start_position <- as.numeric(m6AVar_hg38$Start_position)
m6AVar_hg38$End_position <- as.numeric(m6AVar_hg38$End_position)
m6AVar_snp <- GRanges(seqnames = m6AVar_hg38$SNP_Chromosome,
                      ranges = IRanges(start = m6AVar_hg38$Start_position, end = m6AVar_hg38$End_position),
                      strand = "*",
                      mutateTo = m6AVar_hg38$Alterative_Base)

m6AVar_snp_hg38 <- liftOver(m6AVar_snp, hg19To38)
m6AVar_snp_hg38 <- unlist(m6AVar_snp_hg38)
all(elementNROWS(m6AVar_snp_hg38)==1)

m6AVar_snp_hg38$m6A_site <- m6AVar_hg38

gr_zhang_snp <- sub_gr_zhang_m6A$SNP_gr
gr_zhang_snp_hg38 <- liftOver(gr_zhang_snp, hg19To38)
gr_zhang_snp_hg38 <- unlist(gr_zhang_snp_hg38)
all(elementNROWS(gr_zhang_snp_hg38)==1)

gr_zhang_snp_hg38$m6A_site <- sub_gr_zhang_m6A

overlaps2 <- findOverlaps(gr_zhang_snp_hg38, m6AVar_snp_hg38)
gr_zhang_snp_hg38 <- gr_zhang_snp_hg38[queryHits(overlaps2)]
m6AVar_snp_hg38 <- m6AVar_snp_hg38[subjectHits(overlaps2)]
unique(gr_zhang_snp_hg38)


m6AVar_hg38 <- m6AVar_snp_hg38$m6A_site
m6AVar_snp_hg38$m6A_site <- NULL
sub_zhang_hg38 <- gr_zhang_snp_hg38$m6A_site

snp_features_m6AVar <- sequenceDerivedFeatures(m6AVar_hg38 + 20, BSgenome.Hsapiens.UCSC.hg38, encoding = "onehot", mutation = m6AVar_snp_hg38)
sub_indx <- rownames(snp_features_m6AVar$ref ) ==  rownames(snp_features_m6AVar$mut )
snp_features_m6AVar$mut <- snp_features_m6AVar$mut[sub_indx,]
snp_features_m6AVar$ref <- snp_features_m6AVar$ref[sub_indx,]

DRACH_index <- (snp_features_m6AVar$mut[,"C_19"]!= 1)&
  (snp_features_m6AVar$mut[,"C_20"]!= 1 & snp_features_m6AVar$mut[,"T_20"]!= 1) &
  (snp_features_m6AVar$mut[,"A_21"] == 1)&
  (snp_features_m6AVar$mut[,"C_22"] == 1)&
  (snp_features_m6AVar$mut[,"G_23"]!= 1)
mean(DRACH_index)

snp_features_m6AVar$mut <- snp_features_m6AVar$mut[DRACH_index,]
snp_features_m6AVar$ref <- snp_features_m6AVar$ref[DRACH_index,]

index <- as.numeric(rownames(snp_features_m6AVar$ref))

h2o.init()
OR_m6AVar_seqonly_xgboost <- mutOddsRatio(seq_only_model_xgboost, snp_features_m6AVar)
OR_m6AVar_seqonly_glm <- mutOddsRatio(seq_only_model_glm, snp_features_m6AVar)
OR_m6AVar_seqonly_rf <- mutOddsRatio(seq_only_model_rf, snp_features_m6AVar)
OR_m6AVar_seqonly_mlp <- mutOddsRatio(seq_only_model_mlp, snp_features_m6AVar)

range(OR_m6AVar_seqonly_xgboost)
range(OR_m6AVar_seqonly_glm)
range(OR_m6AVar_seqonly_rf)
range(OR_m6AVar_seqonly_mlp)


OR_m6AVar_seqonly_xgboost <- log(OR_m6AVar_seqonly_xgboost)
OR_m6AVar_seqonly_glm <- log(OR_m6AVar_seqonly_glm)
OR_m6AVar_seqonly_rf <- log(OR_m6AVar_seqonly_rf)
OR_m6AVar_seqonly_mlp <- log(OR_m6AVar_seqonly_mlp)


sub_zhang_hg38[index]$beta <- as.numeric(sub_zhang_hg38[index]$beta)
cor.test(sub_zhang_hg38[index]$beta, OR_m6AVar_seqonly_xgboost)
cor.test(sub_zhang_hg38[index]$beta, OR_m6AVar_seqonly_glm)
cor.test(sub_zhang_hg38[index]$beta, OR_m6AVar_seqonly_rf)
cor.test(sub_zhang_hg38[index]$beta, OR_m6AVar_seqonly_mlp)


genomicFeatures_gr <- readRDS(file.path(data_dir, "motif_exon_hg38_gr_gfeatures.rds"))
fol <- findOverlaps(m6AVar_hg38, genomicFeatures_gr)
sub_motif_gr <- m6AVar_hg38[queryHits(fol)]
sub_snp_gr <- GRanges(seqnames = sub_motif_gr$SNP_Chromosome,
                      ranges = IRanges(start = sub_motif_gr$Start_position, 
                                       end = sub_motif_gr$End_position),
                      strand = "*",
                      mutateTo = sub_motif_gr$Alterative_Base            )
sub_snp_gr <- liftOver(sub_snp_gr, hg19To38)
sub_snp_gr <- unlist(sub_snp_gr)
all(elementNROWS(sub_snp_gr)==1)
gnm <- as.data.frame(mcols(genomicFeatures_gr[subjectHits(fol)]))


snp_features_m6AVar2 <- sequenceDerivedFeatures(sub_motif_gr + 20, BSgenome.Hsapiens.UCSC.hg38, mutation = sub_snp_gr, encoding = "onehot")

sub_indx <- rownames( snp_features_m6AVar2$ref ) ==  rownames( snp_features_m6AVar2$mut )
snp_features_m6AVar2$ref <- snp_features_m6AVar2$ref[sub_indx,]
snp_features_m6AVar2$mut <- snp_features_m6AVar2$mut[sub_indx,]

index2 <- as.numeric(rownames(snp_features_m6AVar2$ref))
gnm <- gnm[index2,]
snp_features_m6AVar2$ref <- cbind(snp_features_m6AVar2$ref, gnm)
snp_features_m6AVar2$mut <- cbind(snp_features_m6AVar2$mut, gnm)


OR_m6AVar_addgen_xgboost <- mutOddsRatio(add_gen_model_xgboost, snp_features_m6AVar2)
OR_m6AVar_addgen_glm <- mutOddsRatio(add_gen_model_glm, snp_features_m6AVar2)
OR_m6AVar_addgen_rf <- mutOddsRatio(add_gen_model_rf, snp_features_m6AVar2)
OR_m6AVar_addgen_mlp <- mutOddsRatio(add_gen_model_mlp, snp_features_m6AVar2)

range(OR_m6AVar_addgen_xgboost)
range(OR_m6AVar_addgen_glm)
range(OR_m6AVar_addgen_rf)
range(OR_m6AVar_addgen_mlp)

OR_m6AVar_addgen_xgboost <- log(OR_m6AVar_addgen_xgboost)
OR_m6AVar_addgen_glm <- log(OR_m6AVar_addgen_glm)
OR_m6AVar_addgen_rf <- log(OR_m6AVar_addgen_rf)
OR_m6AVar_addgen_mlp <- log(OR_m6AVar_addgen_mlp)

sub_zhang_hg38[index2]$beta <- as.numeric(sub_zhang_hg38[index2]$beta)
cor.test(sub_zhang_hg38[index2]$beta, OR_m6AVar_addgen_xgboost)
cor.test(sub_zhang_hg38[index2]$beta, OR_m6AVar_addgen_glm)
cor.test(sub_zhang_hg38[index2]$beta, OR_m6AVar_addgen_rf)
cor.test(sub_zhang_hg38[index2]$beta, OR_m6AVar_addgen_mlp)





# m6AVar/Xiong
gr_xiong_m6A_filtered2 <- readRDS("gr_xiong_m6A_filtered2.rds")
get_first_allele <- function(x) {
  if (length(x) > 1) {
    return(x[1])
  } else {
    return(x)
  }
}
gr_xiong_m6A_filtered2$snp_info <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$snp_info), ",\\s*"), function(x) x[1])
gr_xiong_m6A_filtered2$snp_position <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$snp_position), ",\\s*"), function(x) x[1])
gr_xiong_m6A_filtered2$snp_chromosome <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$snp_chromosome), ",\\s*"), function(x) x[1])
gr_xiong_m6A_filtered2$beta <- sapply(strsplit(as.character(gr_xiong_m6A_filtered2$beta), ",\\s*"), function(x) x[1])

gr_xiong_m6A_filtered2$snp_position <- as.numeric(gr_xiong_m6A_filtered2$snp_position)
gr_xiong_m6A_snp <- GRanges(seqnames = gr_xiong_m6A_filtered2$snp_chromosome,
                            ranges = IRanges(start = gr_xiong_m6A_filtered2$snp_position, end = gr_xiong_m6A_filtered2$snp_position),
                            strand = "*",
                            mutateTo = gr_xiong_m6A_filtered2$snp_info)

gr_xiong_m6A_filtered2$beta <- as.numeric(gr_xiong_m6A_filtered2$beta)
gr_xiong_m6A_filtered2$SNP_gr <- gr_xiong_m6A_snp



hg19To38 <- import.chain("hg19ToHg38.over.chain")
m6AVar_gr_hg38 <- liftOver(m6AVar_gr, hg19To38)
m6AVar_gr_hg38 <- unlist(m6AVar_gr_hg38)

m6AVar_gr_hg38_motif <- predictiveFeatures::sampleSequence("DRACH", m6AVar_gr_hg38, BSgenome.Hsapiens.UCSC.hg38)
m6AVar_gr_hg38_motif <- m6AVar_gr_hg38_motif - 2

overlaps <- findOverlaps(gr_xiong_m6A_filtered2, m6AVar_gr_hg38_motif)
sub_gr_xiong <- gr_xiong_m6A_filtered2[queryHits(overlaps)]
sub_m6AVar_gr_hg38 <- m6AVar_gr_hg38_motif[subjectHits(overlaps)]
unique(sub_gr_xiong)


sub_m6AVar_gr_hg38$SNP_Chromosome <- NA
sub_m6AVar_gr_hg38$Start_position <- NA
sub_m6AVar_gr_hg38$End_position <- NA
sub_m6AVar_gr_hg38$Reference_Base <- NA
sub_m6AVar_gr_hg38$Alterative_Base <- NA
sub_m6AVar_gr_hg38$m6A_Function <- NA


for (i in seq_along(sub_m6AVar_gr_hg38)) {
  
  gr_sub <- sub_m6AVar_gr_hg38[i]
  
  hits <- findOverlaps(gr_sub, m6AVar_gr_hg38)
  
  SNP_Chromosome <- m6AVar_gr[subjectHits(hits)]$SNP_Chromosome
  Start_position <- m6AVar_gr[subjectHits(hits)]$Start_position
  End_position <- m6AVar_gr[subjectHits(hits)]$End_position
  Reference_Base <- m6AVar_gr[subjectHits(hits)]$Reference_Base
  Alterative_Base <- m6AVar_gr[subjectHits(hits)]$Alterative_Base
  m6A_Function <- m6AVar_gr[subjectHits(hits)]$m6A_Function
  
  sub_m6AVar_gr_hg38$SNP_Chromosome[i] <- paste(SNP_Chromosome, collapse = ",")
  sub_m6AVar_gr_hg38$Start_position[i] <- paste(Start_position, collapse = ",")
  sub_m6AVar_gr_hg38$End_position[i] <- paste(End_position, collapse = ",")
  sub_m6AVar_gr_hg38$Reference_Base[i] <- paste(Reference_Base, collapse = ",")
  sub_m6AVar_gr_hg38$Alterative_Base[i] <- paste(Alterative_Base, collapse = ",")
  sub_m6AVar_gr_hg38$m6A_Function[i] <- paste(m6A_Function, collapse = ",")
  
}


# saveRDS(sub_m6AVar_gr_hg38, "sub_m6AVar_gr_hg38.rds")
sub_m6AVar_xiong <- readRDS("sub_m6AVar_gr_hg38.rds")

get_first_allele <- function(x) {
  if (length(x) > 1) {
    return(x[1])
  } else {
    return(x)
  }
}

sub_m6AVar_xiong$SNP_Chromosome <- sapply(strsplit(as.character(sub_m6AVar_xiong$SNP_Chromosome), ",\\s*"), function(x) x[1])
sub_m6AVar_xiong$Start_position <- sapply(strsplit(as.character(sub_m6AVar_xiong$Start_position), ",\\s*"), function(x) x[1])
sub_m6AVar_xiong$End_position <- sapply(strsplit(as.character(sub_m6AVar_xiong$End_position), ",\\s*"), function(x) x[1])
sub_m6AVar_xiong$Reference_Base <- sapply(strsplit(as.character(sub_m6AVar_xiong$Reference_Base), ",\\s*"), function(x) x[1])
sub_m6AVar_xiong$Alterative_Base <- sapply(strsplit(as.character(sub_m6AVar_xiong$Alterative_Base), ",\\s*"), function(x) x[1])
sub_m6AVar_xiong$m6A_Function <- sapply(strsplit(as.character(sub_m6AVar_xiong$m6A_Function), ",\\s*"), function(x) x[1])


sub_m6AVar_xiong$Start_position <- as.numeric(sub_m6AVar_xiong$Start_position)
sub_m6AVar_xiong$End_position <- as.numeric(sub_m6AVar_xiong$End_position)
m6AVar_snp2 <- GRanges(seqnames = sub_m6AVar_xiong$SNP_Chromosome,
                       ranges = IRanges(start = sub_m6AVar_xiong$Start_position, end = sub_m6AVar_xiong$End_position),
                       strand = "*",
                       mutateTo = sub_m6AVar_xiong$Alterative_Base)
m6AVar_snp2_hg38 <- liftOver(m6AVar_snp2, hg19To38)
m6AVar_snp2_hg38 <- unlist(m6AVar_snp2_hg38)

xiong_snp2 <- sub_gr_xiong$SNP_gr

overlaps2 <- findOverlaps(xiong_snp2, m6AVar_snp2_hg38) #无
xiong_snp2 <- xiong_snp2[queryHits(overlaps2)]
m6AVar_snp2_hg38 <- m6AVar_snp2_hg38[subjectHits(overlaps2)]







# m6AVar/RMDisease2 261/26
m6AVar <- read_csv("m6AVar.csv")
Sample_PeakID <- data.frame(do.call('rbind', strsplit(m6AVar$Sample_PeakID, "[:.;]")))
peak_chr <- Sample_PeakID$X1
peak_start_position <- as.numeric(Sample_PeakID$X2)
peak_end_position <- as.numeric(Sample_PeakID$X4)
m6AVar$Start_position <- as.numeric(m6AVar$Start_position)
m6AVar$End_position <- as.numeric(m6AVar$End_position)
m6AVar_gr <- GRanges(seqnames = peak_chr,
                     ranges = IRanges(start = peak_start_position, end = peak_end_position),
                     strand = m6AVar$Strand,
                     SNP_Chromosome = m6AVar$SNP_Chromosome,
                     Start_position = m6AVar$Start_position,
                     End_position = m6AVar$End_position,
                     Reference_Base = m6AVar$Reference_Base,
                     Alterative_Base = m6AVar$Alterative_Base,
                     m6A_Function = m6AVar$m6A_Function)


hg19To38 <- import.chain("hg19ToHg38.over.chain")
m6AVar_gr_hg38 <- liftOver(m6AVar_gr, hg19To38)
m6AVar_gr_hg38 <- unlist(m6AVar_gr_hg38)

m6AVar_gr_motif <- predictiveFeatures::sampleSequence("DRACH", m6AVar_gr_hg38, BSgenome.Hsapiens.UCSC.hg38)
m6AVar_gr_motif <- m6AVar_gr_motif - 2

RMDisease_motif_hg38 <- readRDS("RMDisease_motif_hg38.rds")

overlaps <- findOverlaps(RMDisease_gr_motif, m6AVar_gr_motif)
RMDisease2 <- RMDisease_gr_motif[queryHits(overlaps)]
m6AVar2 <- m6AVar_gr_motif[subjectHits(overlaps)]


m6AVar2$SNP_Chromosome <- NA
m6AVar2$Start_position <- NA
m6AVar2$End_position <- NA


for (i in seq_along(m6AVar2)) {
  
  gr_sub <- m6AVar2[i]
  
  hits <- findOverlaps(gr_sub, m6AVar_gr_hg38)
  
  SNP_Chromosome <- m6AVar_gr_hg38[subjectHits(hits)]$SNP_Chromosome
  Start_position <- m6AVar_gr_hg38[subjectHits(hits)]$Start_position
  End_position <- m6AVar_gr_hg38[subjectHits(hits)]$End_position
  
  m6AVar2$SNP_Chromosome[i] <- paste(SNP_Chromosome, collapse = ",")
  m6AVar2$Start_position[i] <- paste(Start_position, collapse = ",")
  m6AVar2$End_position[i] <- paste(End_position, collapse = ",")
  
}


# saveRDS(m6AVar2, "m6AVar2_RMD.rds")
m6AVar2 <- readRDS("m6AVar2_RMD.rds")

get_first_allele <- function(x) {
  if (length(x) > 1) {
    return(x[1])
  } else {
    return(x)
  }
}

m6AVar2$SNP_Chromosome <- sapply(strsplit(as.character(m6AVar2$SNP_Chromosome), ",\\s*"), function(x) x[1])
m6AVar2$Start_position <- sapply(strsplit(as.character(m6AVar2$Start_position), ",\\s*"), function(x) x[1])
m6AVar2$End_position <- sapply(strsplit(as.character(m6AVar2$End_position), ",\\s*"), function(x) x[1])

m6AVar2$Start_position <- as.numeric(m6AVar2$Start_position)
m6AVar2$End_position <- as.numeric(m6AVar2$End_position)
m6AVar2_snp2 <- GRanges(seqnames = m6AVar2$SNP_Chromosome,
                       ranges = IRanges(start = m6AVar2$Start_position, end = m6AVar2$End_position),
                       strand = "*")
m6AVar2_snp2 <- liftOver(m6AVar2_snp2, hg19To38)
m6AVar2_snp2 <- unlist(m6AVar2_snp2)


RMDisease2$SNP_ChromStart <- as.numeric(RMDisease2$SNP_ChromStart)
RMDisease2$SNP_ChromEnd <- as.numeric(RMDisease2$SNP_ChromEnd)
RMDisease2_snp <- GRanges(seqnames = seqnames(RMDisease2),
                         ranges = IRanges(start = RMDisease2$SNP_ChromStart, end = RMDisease2$SNP_ChromEnd),
                         strand = RMDisease2$SNP_Strand,
                         ref = RMDisease2$ref,
                         alt = RMDisease2$alt,
                         AL = RMDisease2$AL)
RMDisease2_snp <- liftOver(RMDisease2_snp, hg19To38)
RMDisease2_snp <- unlist(RMDisease2_snp)
RMDisease2_snp$m6A_site <- RMDisease2

overlaps2 <- findOverlaps(m6AVar2_snp2, RMDisease2_snp)

RMDisease2_snp <- RMDisease2_snp[subjectHits(overlaps2)]
m6AVar2_snp2 <- m6AVar2_snp2[queryHits(overlaps2)]
RMD_m6AVar_m6A <- RMDisease2_snp$m6A_site
RMD_m6AVar_m6A$SNP <- RMDisease2_snp





# RMDisease2/RMVar 85044/6977
RMVar <- read_csv("RMVar.csv")
sample_peak_id <- data.frame(do.call('rbind', strsplit(RMVar$sample_peak_id, "[:.,]")))
RMVar$peak_chr <- sample_peak_id$X1
RMVar$peak_start_position <- as.numeric(sample_peak_id$X2)
RMVar$peak_end_position <- as.numeric(sample_peak_id$X4)
RMVar <- RMVar[!RMVar$peak_chr == 'chrMT',]
RMVar_gr <- GRanges(seqnames = RMVar$peak_chr,
                    ranges = IRanges(start = RMVar$peak_start_position, end = RMVar$peak_end_position),
                    strand = RMVar$strand,
                    snp_start = RMVar$snp_start,
                    snp_end = RMVar$snp_end,
                    ref = RMVar$reference_base,
                    alt = RMVar$alterative_base,
                    mf = RMVar$modification_function)

RMVar_gr_motif <- predictiveFeatures::sampleSequence("DRACH", RMVar_gr, BSgenome.Hsapiens.UCSC.hg38)
RMVar_gr_motif <- RMVar_gr_motif - 2



RMDisease_motif_hg38 <- readRDS("RMDisease_motif_hg38.rds")
RMDisease_motif_hg38$SNP_ChromStart <- as.numeric(RMDisease_motif_hg38$SNP_ChromStart)
RMDisease_motif_hg38$SNP_ChromEnd <- as.numeric(RMDisease_motif_hg38$SNP_ChromEnd)
RMDisease_snp <- GRanges(seqnames = seqnames(RMDisease_motif_hg38),
                         ranges = IRanges(start = RMDisease_motif_hg38$SNP_ChromStart, end = RMDisease_motif_hg38$SNP_ChromEnd),
                         strand = RMDisease_motif_hg38$SNP_Strand,
                         ref = RMDisease_motif_hg38$ref,
                         alt = RMDisease_motif_hg38$alt,
                         AL = RMDisease_motif_hg38$AL)
RMDisease_motif_hg38$SNP_gr <- RMDisease_snp

overlaps <- findOverlaps(RMVar_gr_motif, RMDisease_motif_hg38)
RMDisease_motif_hg38 <- RMDisease_motif_hg38[subjectHits(overlaps)]
RMvar_RMD <- RMVar_gr_motif[queryHits(overlaps)]


RMvar_RMD$snp_start <- NA

for (i in seq_along(RMvar_RMD)) {
  
  gr_sub <- RMvar_RMD[i]
  
  hits <- findOverlaps(gr_sub, RMVar_gr)
  
  snp_start <- RMVar_gr[subjectHits(hits)]$snp_start
  
  RMvar_RMD$snp_start[i] <- paste(snp_start, collapse = ",")
  
}


# saveRDS(RMvar_RMD, "RMvar_RMD.rds")
RMvar_RMD <- readRDS("RMvar_RMD.rds")

get_first_allele <- function(x) {
  if (length(x) > 1) {
    return(x[1])
  } else {
    return(x)
  }
}
RMvar_RMD$snp_start <- sapply(strsplit(as.character(RMvar_RMD$snp_start), ",\\s*"), function(x) x[1])
RMvar_RMD$snp_start <- as.numeric(RMvar_RMD$snp_start)

RMVar2_snp2 <- GRanges(seqnames = seqnames(RMvar_RMD),
                      ranges = IRanges(start = RMvar_RMD$snp_start, end = RMvar_RMD$snp_start),
                      strand = "*")

RMVar2_snp2$m6A_site <- RMvar_RMD


RMDisease_snp <- RMDisease_motif_hg38$SNP_gr
hg19To38 <- import.chain("hg19ToHg38.over.chain")
RMDisease_snp_hg38 <- liftOver(RMDisease_snp, hg19To38)
RMDisease_snp_hg38 <- unlist(RMDisease_snp_hg38)
all(elementNROWS(RMDisease_snp_hg38)==1)

overlaps2 <- findOverlaps(RMVar2_snp2, RMDisease_snp_hg38)
RMvar_RMD_m6A <- RMVar2_snp2[queryHits(overlaps2)]$m6A_site
RMvar_RMD_m6A$SNP <- RMVar2_snp2[queryHits(overlaps2)]




# m6AVar/RMVar 358259

m6AVar <- read_csv("m6AVar.csv")
Sample_PeakID <- data.frame(do.call('rbind', strsplit(m6AVar$Sample_PeakID, "[:.;]")))
peak_chr <- Sample_PeakID$X1
peak_start_position <- as.numeric(Sample_PeakID$X2)
peak_end_position <- as.numeric(Sample_PeakID$X4)
m6AVar$Start_position <- as.numeric(m6AVar$Start_position)
m6AVar$End_position <- as.numeric(m6AVar$End_position)

m6AVar_gr <- GRanges(seqnames = peak_chr,
                     ranges = IRanges(start = peak_start_position, end = peak_end_position),
                     strand = m6AVar$Strand,
                     SNP_Chromosome = m6AVar$SNP_Chromosome,
                     Start_position = m6AVar$Start_position,
                     End_position = m6AVar$End_position,
                     Reference_Base = m6AVar$Reference_Base,
                     Alterative_Base = m6AVar$Alterative_Base,
                     m6A_Function = m6AVar$m6A_Function)

hg19To38 <- import.chain("hg19ToHg38.over.chain")
m6AVar_gr_hg38 <- liftOver(m6AVar_gr, hg19To38)
m6AVar_gr_hg38 <- unlist(m6AVar_gr_hg38)

m6AVar_gr_motif <- predictiveFeatures::sampleSequence("DRACH", m6AVar_gr_hg38, BSgenome.Hsapiens.UCSC.hg38)
m6AVar_gr_motif <- m6AVar_gr_motif - 2

RMVar <- read_csv("RMVar.csv")
sample_peak_id <- data.frame(do.call('rbind', strsplit(RMVar$sample_peak_id, "[:.,]")))
RMVar$peak_chr <- sample_peak_id$X1
RMVar$peak_start_position <- as.numeric(sample_peak_id$X2)
RMVar$peak_end_position <- as.numeric(sample_peak_id$X4)
RMVar <- RMVar[!RMVar$peak_chr == 'chrMT',]
RMVar_gr <- GRanges(seqnames = RMVar$peak_chr,
                    ranges = IRanges(start = RMVar$peak_start_position, end = RMVar$peak_end_position),
                    strand = RMVar$strand,
                    snp_start = RMVar$snp_start,
                    snp_end = RMVar$snp_end,
                    ref = RMVar$reference_base,
                    alt = RMVar$alterative_base,
                    mf = RMVar$modification_function)

RMVar_gr_motif <- predictiveFeatures::sampleSequence("DRACH", RMVar_gr, BSgenome.Hsapiens.UCSC.hg38)
RMVar_gr_motif <- RMVar_gr_motif - 2

overlaps <- findOverlaps(m6AVar_gr_motif, RMVar_gr_motif)
RMvarm3 <- RMVar_gr_motif[subjectHits(overlaps)]







# three databases (m6AVar & RMD) & RMVar 1
RMVar <- read_csv("RMVar.csv")
sample_peak_id <- data.frame(do.call('rbind', strsplit(RMVar$sample_peak_id, "[:.,]")))
RMVar$peak_chr <- sample_peak_id$X1
RMVar$peak_start_position <- as.numeric(sample_peak_id$X2)
RMVar$peak_end_position <- as.numeric(sample_peak_id$X4)
RMVar <- RMVar[!RMVar$peak_chr == 'chrMT',]
RMVar_gr <- GRanges(seqnames = RMVar$peak_chr,
                    ranges = IRanges(start = RMVar$peak_start_position, end = RMVar$peak_end_position),
                    strand = RMVar$strand,
                    snp_start = RMVar$snp_start,
                    snp_end = RMVar$snp_end,
                    ref = RMVar$reference_base,
                    alt = RMVar$alterative_base,
                    mf = RMVar$modification_function)
RMVar_gr_motif <- predictiveFeatures::sampleSequence("DRACH", RMVar_gr, BSgenome.Hsapiens.UCSC.hg38)
RMVar_gr_motif <- RMVar_gr_motif - 2

overlaps <- findOverlaps(RMD_m6AVar_m6A, RMVar_gr_motif)
dd3 <- RMVar_gr_motif[subjectHits(overlaps)]
dd12 <- RMD_m6AVar_m6A[queryHits(overlaps)]


d12_snp2 <- dd12$SNP


dd3$snp_start <- NA
for (i in seq_along(dd3)) {
  
  gr_sub <- dd3[i]
  
  hits <- findOverlaps(gr_sub, RMVar_gr)
  
  snp_start <- RMVar_gr[subjectHits(hits)]$snp_start
  
  dd3$snp_start[i] <- paste(snp_start, collapse = ",")
  
}

# saveRDS(dd3, "dd3.rds")
dd3 <- readRDS("dd3.rds")

dd3$snp_start <- sapply(strsplit(as.character(dd3$snp_start), ",\\s*"), function(x) x[1])
dd3$snp_start <- as.numeric(dd3$snp_start)

d3_snp2 <- GRanges(seqnames = seqnames(dd3),
                       ranges = IRanges(start = dd3$snp_start, end = dd3$snp_start),
                       strand = "*")
 
overlaps2 <- findOverlaps(d12_snp2_hg38, d3_snp2)





