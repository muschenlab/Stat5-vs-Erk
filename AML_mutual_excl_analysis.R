options(stringsAsFactors = F)
library(tidyverse)

# load mutation data
mut <- readxl::read_excel("Supplementary-Table 2_2019-11-21.xlsx")[,-c(1:4)]

# convert to gene alteration matrix
gam <- as.matrix(mut)

# gene-pathway assignments
genelist <- list(STAT = c("CALR","CBL","JAK1","JAK2","JAK3","KIT","MPL",
                          "PDGFRA","PDGFRB","STAT5A"),
                 ERK = c("RAF1","NRAS","KRAS","PTPN11","NF1","FLT3"))

# check mut counts for each pathway
pathtab <- table(STAT5 = rowSums(gam[,colnames(gam) %in% genelist$STAT]) > 0,
                 ERK = rowSums(gam[,colnames(gam) %in% genelist$ERK]) > 0)
chisq.test(pathtab)

# heatmap
stat_idx <- which(colnames(gam) %in% genelist$STAT)
erk_idx <- which(colnames(gam) %in% genelist$ERK)
plotdat <- gam
plotdat[,erk_idx] <- sapply(plotdat[,erk_idx], function(x) ifelse(x > 0, -1, 0))
plotdat <- plotdat[order(-rowSums(plotdat)), c(stat_idx, erk_idx)]
pheatmap::pheatmap(plotdat, scale = "none", 
                   cluster_rows = F, cluster_cols = F,
                   color = c("firebrick","white","forestgreen"),
                   show_rownames = F)

# run fishers for all pairs 
gene_n <- ncol(gam)
ftest <- lapply(1:(gene_n-1),  function(i) {
  lapply((i+1):gene_n, function(j) {
    f <- fisher.test(gam[,i], gam[,j])
    c <- chisq.test(gam[,i], gam[,j])
    return(data.frame(geneA = colnames(gam)[i],
                      geneB = colnames(gam)[j],
                      nA = colSums(gam, na.rm = T)[i],
                      nB = colSums(gam, na.rm = T)[j],
                      nAB = table(rowSums(gam[,c(i,j)], na.rm = T))["2"],
                      nexpected = c$expected[4],
                      pval = f$p.value,
                      CI_lower = f$conf.int[1],
                      CI_higher = f$conf.int[2],
                      logOR = log2(f$estimate),
                      stringsAsFactors = F))
  }) %>% bind_rows()
}) %>% bind_rows()
ftest$nAB[is.na(ftest$nAB)] <- 0
ftest$ndiff <- ftest$nAB - ftest$nexpected
ftest$FDR <- p.adjust(ftest$pval, "fdr")

# deal with infs & generate co-occurence score
ftest$logOR[ftest$logOR > 10] <- 10
ftest$logOR[ftest$logOR < -10] <- -10
ftest$score <- -log10(ftest$pval) * ftest$logOR
ftest <- ftest[order(ftest$pval, ftest$logOR),]

# assign genes to pathways
ftest$pathwayA <- ifelse(ftest$geneA %in% genelist$STAT, "STAT5", "ERK")
ftest$pathwayB <- ifelse(ftest$geneB %in% genelist$STAT, "STAT5", "ERK")
ftest$pathcomb <- paste0(ftest$pathwayA, "-", ftest$pathwayB)
ftest[ftest$pathcomb == "STAT5-ERK",]$pathcomb <- "ERK-STAT5"
ftest$pathcat <- ifelse(ftest$pathwayA == ftest$pathwayB, "INTRA", "INTER")

# boxplot
plotdat <- ftest
plotdat$freqA <- sapply(plotdat$nA, function(x) x/nrow(gam))
plotdat$freqB <- sapply(plotdat$nB, function(x) x/nrow(gam))
plotdat <- plotdat[-which((plotdat$freqA < .05 | plotdat$freqB < .05) & plotdat$FDR > .5),]
ggplot(plotdat, aes(x = pathcomb, y = logOR, fill = pathcomb)) +
  geom_boxplot(outlier.alpha = 0, col = "#333333") +
  geom_point(position = position_jitterdodge(), shape = 21, size = 1.8) +
  theme_bw(base_size = 16) +
  ylab("Odds Ratio (log2)") +
  xlab("") +
  ylim(c(-10,10)) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#b4000087","#707070","#267f269b"))
ggplot(plotdat, aes(x = pathcat, y = logOR, fill = pathcat)) +
  geom_boxplot(outlier.alpha = 0, col = "#333333") +
  geom_point(position = position_jitterdodge(), shape = 21, size = 1.8) +
  theme_bw(base_size = 16) +
  ylab("Odds Ratio (log2)") +
  xlab("") +
  ylim(c(-10,10)) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#707070", "#dddddd"))

# stats
aovres <- aov(logOR ~ pathcomb, data = plotdat)
TukeyHSD(aovres)
t.test(logOR ~ pathcat, data = plotdat)

# volcano
plotdat <- ftest
plotdat$label <- paste0(plotdat$geneA, "-", plotdat$geneB)
plotdat$freqA <- sapply(plotdat$nA, function(x) x/nrow(gam))
plotdat$freqB <- sapply(plotdat$nB, function(x) x/nrow(gam))
plotdat$size <- rowMeans(plotdat[,c("freqA","freqB")])
ggplot(plotdat, aes(x = ndiff, y = -log10(pval), col = pathcomb, label = label)) +
  geom_hline(yintercept = 0, alpha = .25) +
  geom_vline(xintercept = 0, alpha = .25) +
  geom_point(aes(size = size)) +
  xlab("Odds Ratio (log2)") +
  ylab("Significance (-log10 p-value)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values = c("#b4000099","#60606099","#267f2699")) +
  scale_x_continuous(limits = c(-28,28))

# network
library(igraph)
plotdat <- ftest
plotdat$weight <- 1
plotdat[plotdat$pathcomb == "ERK-STAT5",]$weight <- 0.5
plotdat$dir <- ifelse(plotdat$logOR < 0, "#aa000099", "#00720099")
g <- graph_from_data_frame(plotdat, directed = F)
g <- set.edge.attribute(g, "width", value = abs(plotdat$score)/8)
g <- set.edge.attribute(g, "color", value = plotdat$dir)
g <- set.vertex.attribute(g, "color", value = ifelse(names(V(g)) %in% genelist$STAT,
                                                     "#f2f2f2ff", "#666666ff"))
mut_counts <- colSums(gam)
mut_counts <- mut_counts[match(names(V(g)), names(mut_counts))]
g <- set.vertex.attribute(g, "size", value = log(mut_counts)*5)
g <- set.edge.attribute(g, "weight", value = plotdat$weight)
l <- layout_with_fr(g)
plot(g, vertex.label.color = "black", layout = l)

# pathway level fishers
pam <- cbind(STAT = rowSums(gam[,colnames(gam) %in% genelist$STAT], na.rm = T),
             ERK = rowSums(gam[,colnames(gam) %in% genelist$ERK], na.rm = T))
pam[pam>1] <- 1
pamdf <- as.data.frame(pam)
fres <- fisher.test(pamdf$STAT, pamdf$ERK)
cres <- chisq.test(pamdf$STAT, pamdf$ERK)
pfres <- data.frame(nA = sum(pamdf$STAT),
                    nB = sum(pamdf$ERK),
                    nAB = table(rowSums(pamdf))["2"],
                    nexpected = cres$expected[4],
                    pval = fres$p.value,
                    CI_lower = fres$conf.int[1],
                    CI_higher = fres$conf.int[2],
                    logOR = log2(fres$estimate),
                    stringsAsFactors = F)
pfres$FCdiff <- log2(pfres$nAB / pfres$nexpected)

# permute pathway level
statlen <- length(genelist$STAT[genelist$STAT %in% colnames(gam)])
permres <- lapply(1:1000, function(n) {
  ridx <- sample(colnames(gam), replace = F)
  rpam <- cbind(rSTAT = rowSums(gam[,colnames(gam) %in% ridx[1:statlen]]),
                rERK = rowSums(gam[,colnames(gam) %in% ridx[-c(1:statlen)]]))
  rpam[rpam>1] <- 1
  rpamdf <- as.data.frame(rpam)
  fres <- fisher.test(rpamdf$rSTAT, rpamdf$rERK)
  chires <- chisq.test(rpamdf$rSTAT, rpamdf$rERK)
  return(data.frame(nA = sum(rpamdf$rSTAT),
                    nB = sum(rpamdf$rERK),
                    nAB = table(rowSums(rpamdf))["2"],
                    nexpected = chires$expected[4],
                    pval = fres$p.value,
                    CI_lower = fres$conf.int[1],
                    CI_higher = fres$conf.int[2],
                    logOR = log2(fres$estimate),
                    stringsAsFactors = F))
}) %>% bind_rows()
permres$ndiff <- permres$nAB - permres$nexpected
permres$FCdiff <- log2(permres$nAB / permres$nexpected)
sum(permres$FCdiff < pfres$FCdiff)/length(permres$FCdiff)

# plot permuted vs observed
ggplot(permres, aes(x=FCdiff)) + 
  geom_histogram(alpha = .9, col = "#b3b3b3", binwidth = .01) +
  geom_vline(xintercept = pfres$FCdiff, col = "firebrick") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())
