
# -- Kaili
# This script is for making figures for analysing results from Majority Vote.

setwd("/Users/kaili/Dropbox (UMass Medical School)/Project/chr_status/ideas/dhs_ctcf/majority_vote/")

library(ggplot2)
library(gridExtra)
##################
# barplot
##################
ctcfstates_count = data.frame(c(124162, 151852, 126849, 158281, 138743, 118343, 124911, 90895, 150668, 
                                      232832),
      c("forebrain", "midbrain",  "hindbrain", "intestine", "kidney", "liver", "lung", "heart",
        "stomach", "total"))
colnames(ctcfstates_count) = c("count", "tissue")
rownames(ctcfstates_count) = c("forebrain", "midbrain",  "hindbrain", "intestine", "kidney", "liver", "lung", "heart",
                               "stomach", "total")

ggplot(ctcfstates_count, aes(x = tissue, y=count, fill=tissue)) + geom_bar(stat = "identity") +
  theme_minimal() + theme(title=element_text(face="bold", size=12),
                          axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="number of CTCF states in each sample") + ylab("num of CTCF state") +
  geom_text(aes(label=count), vjust=-0.3, size=3.5)
ggsave("9sample_CTCFstates_count.pdf")
ggsave("9sample_CTCFstates_count.png")


###
majority_ctcfstates = data.frame(c(54116, 25533, 18258, 13608, 12006, 11525, 11966, 15006, 70814),
                                 c("1", "2", "3", "4", "5", "6", "7", "8", "9"))
colnames(majority_ctcfstates) = c("count", "num")
rownames(majority_ctcfstates) = c("1", "2", "3", "4", "5", "6", "7", "8", "9")

ggplot(majority_ctcfstates, aes(x = num, y=count)) + geom_bar(stat = "identity") +
  theme_minimal() + theme(title=element_text(face="bold", size=12)) +
  ylab("num of CTCF state") + xlab("number of samples") +
  geom_text(aes(label=count), vjust=-0.3, size=3.5)
ggsave("9sample_CTCFstates_votes.pdf")
ggsave("9sample_CTCFstates_votes.png")


##################
# with peak summit
##################

# liver_14.5
a = cbind(c(54116, 25533, 18258, 13608, 12006, 11525, 11966, 15006, 70814)-c(712,614,644,720,859,1062,1546,2856,38738),
          c("1", "2", "3", "4", "5", "6", "7", "8", "9"), rep("a",9))
b = cbind(c(712,614,644,720,859,1062,1546,2856,38738),
          c("1", "2", "3", "4", "5", "6", "7", "8", "9"), rep("b",9))
dat = data.frame(rbind(a,b))
colnames(dat) = c("num","count","type")
dat$num = as.numeric(as.character(dat$num))

ggplot(dat, aes(x=count, y=num, fill=type)) + geom_bar(stat = "identity") +
  theme_minimal() + theme(title = element_text(face="bold", size=12),
                          legend.title = element_blank()) +
  ylab("num of recalled CTCF peaks") + xlab("number of samples counted as CTCF states") + 
  labs(title="liver_14.5") + 
  scale_fill_manual(values=c("#AEAFAE", "#E73A2F"), 
                    labels = c("without_CTCF-peak-summit","with_CTCF-peak-summit"))
ggsave("9sample_CTCFstates_votes_liver_14.5_summit.pdf")
ggsave("9sample_CTCFstates_votes_liver_14.5_summit.png")

# lung_14.5
a = cbind(c(54116, 25533, 18258, 13608, 12006, 11525, 11966, 15006, 70814)-c(489,542,690,820,998,1346,2148,3583,42791),
          c("1", "2", "3", "4", "5", "6", "7", "8", "9"), rep("a",9))
b = cbind(c(489,542,690,820,998,1346,2148,3583,42791),
          c("1", "2", "3", "4", "5", "6", "7", "8", "9"), rep("b",9))
dat = data.frame(rbind(a,b))
colnames(dat) = c("num","count","type")
dat$num = as.numeric(as.character(dat$num))

ggplot(dat, aes(x=count, y=num, fill=type)) + geom_bar(stat = "identity") +
  theme_minimal() + theme(title = element_text(face="bold", size=12),
                          legend.title = element_blank()) +
  ylab("num of recalled CTCF peaks") + xlab("number of samples counted as CTCF states") + 
  labs(title="lung_14.5") + 
  scale_fill_manual(values=c("#AEAFAE", "#E73A2F"), 
                    labels = c("without_CTCF-peak-summit","with_CTCF-peak-summit"))
ggsave("9sample_CTCFstates_votes_lung_14.5_summit.pdf")
ggsave("9sample_CTCFstates_votes_lung_14.5_summit.png")


##################
# overlapping peaks
##################

# liver_14.5
p = round(c(3259,2640,2504,2376,2529,2756,3294,5663,49214)/c(54116, 25533, 18258, 13608, 12006, 11525, 11966, 15006, 70814),3)*100
a = cbind(c(54116, 25533, 18258, 13608, 12006, 11525, 11966, 15006, 70814)-c(3259,2640,2504,2376,2529,2756,3294,5663,49214),
          c("1", "2", "3", "4", "5", "6", "7", "8", "9"), rep("a",9))
b = cbind(c(3259,2640,2504,2376,2529,2756,3294,5663,49214),
          c("1", "2", "3", "4", "5", "6", "7", "8", "9"), rep("b",9))
dat = data.frame(rbind(a,b))
colnames(dat) = c("num","count","type")
dat$num = as.numeric(as.character(dat$num))

ggplot(dat, aes(x=count, y=num, fill=type)) + geom_bar(stat = "identity") +
  theme_minimal() + theme(title = element_text(face="bold", size=12),
                          legend.title = element_blank()) +
  ylab("num of recalled CTCF peaks") + xlab("number of samples counted as CTCF states") + 
  labs(title="liver_14.5") + 
  scale_fill_manual(values=c("#AEAFAE", "#E73A2F"), 
                    labels = c("non_overlap_CTCF_peaks","overlap_CTCF_peaks")) +
  geom_text(p, vjust=-0.3, size=3.5)

ggsave("9sample_CTCFstates_votes_liver_14.5_peaks.pdf")
ggsave("9sample_CTCFstates_votes_liver_14.5_peaks.png")

# lung_14.5
a = cbind(c(54116, 25533, 18258, 13608, 12006, 11525, 11966, 15006, 70814)-c(3443,2978,2971,2743,2766,3026,3771,5743,51653),
          c("1", "2", "3", "4", "5", "6", "7", "8", "9"), rep("a",9))
b = cbind(c(3443,2978,2971,2743,2766,3026,3771,5743,51653),
          c("1", "2", "3", "4", "5", "6", "7", "8", "9"), rep("b",9))
dat = data.frame(rbind(a,b))
colnames(dat) = c("num","count","type")
dat$num = as.numeric(as.character(dat$num))

ggplot(dat, aes(x=count, y=num, fill=type)) + geom_bar(stat = "identity") +
  theme_minimal() + theme(title = element_text(face="bold", size=12),
                          legend.title = element_blank()) +
  ylab("num of recalled CTCF peaks") + xlab("number of samples counted as CTCF states") + 
  labs(title="lung_14.5") + 
  scale_fill_manual(values=c("#AEAFAE", "#E73A2F"), 
                    labels = c("non_overlap_CTCF_peaks","overlap_CTCF_peaks"))
ggsave("9sample_CTCFstates_votes_lung_14.5_peak.pdf")
ggsave("9sample_CTCFstates_votes_lung_14.5_peak.png")

##################
# peak recall - barplot
##################
peak_recall_liver_14.5 = data.frame(c(50909, 45892, 45136, 31022),
                                 c("total", "IDEAS", "majority_vote", "9-impute-11"))
colnames(peak_recall_liver_14.5) = c("count", "type")
rownames(peak_recall_liver_14.5) = c("total", "IDEAS", "majority_vote", "9-impute-11")

ggplot(peak_recall_liver_14.5, aes(x = type, y=count, fill=type)) + 
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.5)) +
  theme_minimal() + theme(title=element_text(face="bold", size=12)) +
  ylab("num of CTCF peaks") + xlab("") + labs(title="liver_14.5") +
  geom_text(aes(label=count), vjust=-0.3, size=3.5) +
  scale_fill_manual(values=c("#FBB30B", "#E73A2F", "#397AF2", "#737373"))
ggsave("peak_recall_liver_14.5_barplot.pdf")
ggsave("peak_recall_liver_14.5_barplot.png")

###
peak_recall_lung_14.5 = data.frame(c(54953, 51427, 50901, 39297),
                                    c("total", "IDEAS", "majority_vote", "9-impute-11"))
colnames(peak_recall_lung_14.5) = c("count", "type")
rownames(peak_recall_lung_14.5) = c("total", "IDEAS", "majority_vote", "9-impute-11")

ggplot(peak_recall_lung_14.5, aes(x = type, y=count, fill=type)) + 
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.5)) +
  theme_minimal() + theme(title=element_text(face="bold", size=12)) +
  ylab("num of CTCF peaks") + xlab("") + labs(title="lung_14.5") +
  geom_text(aes(label=count), vjust=-0.3, size=3.5) +
  scale_fill_manual(values=c("#FBB30B", "#E73A2F", "#397AF2", "#737373"))
ggsave("peak_recall_lung_14.5_barplot.pdf")
ggsave("peak_recall_lung_14.5_barplot.png")

###
peak_recall_lung_0 = data.frame(c(64742, 58641, 56133, 58558),
                                   c("total", "IDEAS", "majority_vote", "9-impute-11"))
colnames(peak_recall_lung_0) = c("count", "type")
rownames(peak_recall_lung_0) = c("total", "IDEAS", "majority_vote", "9-impute-11")

ggplot(peak_recall_lung_0, aes(x = type, y=count, fill = type)) + 
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.5)) +
  theme_minimal() + theme(title=element_text(face="bold", size=12)) +
  ylab("num of CTCF peaks") + xlab("") + labs(title="lung_0") +
  geom_text(aes(label=count), vjust=-0.3, size=3.5) +
  scale_fill_manual(values=c("#FBB30B", "#E73A2F", "#397AF2", "#737373"))
ggsave("peak_recall_lung_0_barplot.pdf")
ggsave("peak_recall_lung_0_barplot.png")

##################
# Venn
##################
library(VennDiagram)

## liver_14.5
pdf("venn_liver_14.5.pdf")
p = draw.triple.venn(area1 = 45879, area2 = 31017, area3 = 45136, n12 = 30253, n23 = 30558, n13 = 43421, 
                 n123 = 29878, category = c("IDEAS", "9impute11", "majority_vote"), lty = "blank", 
                 fill = c("#E73A2F", "#FBB30B", "#397AF2"), scaled=TRUE)
grid.arrange(gTree(children=p), top="liver_14.5\ntotal 50,909 peaks\n47,678 peaks in the Venn (93.7%)")
dev.off()
grid.newpage()

## lung_14.5
pdf("venn_lung_14.5.pdf")
p = draw.triple.venn(area1 = 51415, area2 = 39289, area3 = 50901, n12 = 38953, n23 = 38664, n13 = 49859, 
                     n123 = 38448, category = c("IDEAS", "9impute11", "majority_vote"), lty = "blank", 
                     fill = c("#E73A2F", "#FBB30B", "#397AF2"), scaled=TRUE)
grid.arrange(gTree(children=p), top="lung_14.5\ntotal 54,953 peaks\n52,577 peaks in the Venn (95.7%)")
dev.off()
##################
# accumalated recall for peak with different count numbers
##################
# liver_14.5
dat = data.frame(c(716, 1337, 1984, 2703, 3564, 4626, 6177, 9032, 47831, 50909),
                 c("=1", "≤2", "≤3", "≤4", "≤5", "≤6", "≤7", "≤8", "≤9", "total_peaks"))
colnames(dat) = c("recall", "count")
rownames(dat) = c("=1", "≤2", "≤3", "≤4", "≤5", "≤6", "≤7", "≤8", "≤9", "total_peaks")

ggplot(dat, aes(x = count, y=recall)) + 
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.5), fill="#AEAFAE") +
  theme_minimal() + theme(title=element_text(face="bold", size=12)) +
  ylab("num of recalled CTCF peaks") + xlab("number of samples counted as CTCF states") + 
  labs(title="liver_14.5") +
  geom_text(aes(label=recall), vjust=-0.3, size=3.5)
ggsave("accumulated_peak_recall_liver_14.5_barplot.pdf")
ggsave("accumulated_peak_recall_liver_14.5_barplot.png")

# lung_14.5
dat = data.frame(c(495, 1046, 1740, 2563, 3562, 4907, 7055, 10637, 53458, 54953),
                 c("=1", "≤2", "≤3", "≤4", "≤5", "≤6", "≤7", "≤8", "≤9", "total_peaks"))
colnames(dat) = c("recall", "count")
rownames(dat) = c("=1", "≤2", "≤3", "≤4", "≤5", "≤6", "≤7", "≤8", "≤9", "total_peaks")

ggplot(dat, aes(x = count, y=recall)) + 
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.5), fill="#AEAFAE") +
  theme_minimal() + theme(title=element_text(face="bold", size=12)) +
  ylab("num of recalled CTCF peaks") + xlab("number of samples counted as CTCF states") + 
  labs(title="lung_14.5") +
  geom_text(aes(label=recall), vjust=-0.3, size=3.5)
ggsave("accumulated_peak_recall_lung_14.5_barplot.pdf")
ggsave("accumulated_peak_recall_lung_14.5_barplot.png")

##################
##################
##################