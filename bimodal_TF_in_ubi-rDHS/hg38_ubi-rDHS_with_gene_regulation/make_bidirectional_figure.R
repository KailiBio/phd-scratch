
# -- Kaili
# This script is for making bidirectional figures.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/ccRE/hg38_rOCR/bidirectional/")

library(ggplot2)

####################
# bidirectional gene type
####################
dat_2end = data.frame(c(2329,439,205,45,23,19,13,13,9,3,3,2,1,1,1,1),
                      c("protein_coding","antisense", "lincRNA", "processed_transcript", 
                        "transcribed_unprocessed_pseudogene", "miRNA", "transcribed_processed_pseudogene",
                        "transcribed_unitary_pseudogene", "misc_RNA", "processed_pseudogene", "sense_intronic",
                        "snoRNA", "non_coding", "ribozyme", "sense_overlapping", "snRNA"))
colnames(dat_2end) = c("num", "genetype")
dat_2end$genetype = factor(dat_2end$genetype)

ggplot(dat_2end, aes(x=genetype, y=num)) + geom_bar(stat = "identity", fill = "#e1e1e1") +
  scale_x_discrete(limits = dat_2end$genetype) +
  theme_minimal() + theme(title=element_text(face="bold", size=12),
                          axis.text.x = element_text(angle=90, hjust=1, face="bold")) +
  labs(title="gene-type for bidirectional pairs that both end overlap ubi-rOCRs") +
  ylab("number of genes") +
  geom_text(aes(label=dat_2end$num), position=position_dodge(width=0.9))
ggsave("gene_type_both_end.pdf", width=9, height = 5)
ggsave("gene_type_both_end.png", width=9, height = 5)


#### all
dat_all = data.frame(c(2760,488,234,47,25,23,15,13,9,3,3,2,1,1,1,1),
                      c("protein_coding","antisense", "lincRNA", "processed_transcript", 
                        "transcribed_unprocessed_pseudogene", "miRNA", "transcribed_processed_pseudogene",
                        "transcribed_unitary_pseudogene", "misc_RNA", "processed_pseudogene", "sense_intronic",
                        "snoRNA", "non_coding", "ribozyme", "sense_overlapping", "snRNA"))
colnames(dat_all) = c("num", "genetype")
dat_all$genetype = factor(dat_all$genetype)

ggplot(dat_all, aes(x=genetype, y=num)) + geom_bar(stat = "identity", fill = "#e1e1e1") +
  scale_x_discrete(limits = dat_all$genetype) +
  theme_minimal() + theme(title=element_text(face="bold", size=12),
                          axis.text.x = element_text(angle=90, hjust=1, face="bold")) +
  labs(title="gene-type for all bidirectional genes that overlaping ubi-rOCRs") +
  ylab("number of genes") +
  geom_text(aes(label=dat_all$num), position=position_dodge(width=0.9))
ggsave("gene_type_all.pdf", width=9, height = 5)
ggsave("gene_type_all.png", width=9, height = 5)
####################
####################