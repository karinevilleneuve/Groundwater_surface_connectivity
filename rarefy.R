# Load package
library(vegan)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(iNEXT)

ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Number of sequence reads", y = "ASV richness")
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0, show.legend=FALSE)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

#  ------------------ Load data  ------------------ 
asv = read.table("/DADA2_output/raw_asv.csv", sep=",", row.names=1, header=TRUE, check.names=FALSE) 
taxa = read.table("/DADA2_output/raw_taxa.csv", sep=",", row.names=1, header=TRUE) 
meta = read.table("/DADA2_output/raw_meta.csv", sep=",", row.names=1, header=TRUE) 
tree = ape::read.tree("/DADA2_output/raw_tree.tree")

# ------------------ Merge into phyloseq object ------------------ 
ps=phyloseq(otu_table(asv, taxa_are_rows=TRUE), tax_table(as.matrix(taxa)), sample_data(meta), phy_tree(tree))
sums=data.frame(sort(colSums(asv)))
sums[1,1] 

# Remove discharge samples 
ps=subset_samples(ps, Sample_type!="Discharge-Forest")
ps=subset_samples(ps, Sample_type!="Discharge-River")
ps=subset_samples(ps, Sample_type!="Discharge-River-B")

# ------------------ Filter low abundant taxa  ------------------
# One of the reasons to filter in this way is to avoid spending much time analyzing taxa that were only rarely seen. 
# This also turns out to be a useful filter of noise (taxa that are actually just artifacts of the data collection process)
# A step that should probably be considered essential for datasets constructed via heuristic OTU-clustering methods, which are notoriously prone to generating spurious taxa.
# Removing taxa with a relative abundance less than 0.005% as recommended by Bokulich et al., 2013

minTotRelAbun = 5e-5
x = taxa_sums(ps)
keepTaxa =  (x / sum(x)) > minTotRelAbun
prunedSet = prune_taxa(keepTaxa, ps)
sort(sample_sums(prunedSet))
sort(sample_sums(ps))


# ------------------ Generate rarefaction curve ------------------ 

# When differences between library sizes is high (such as 10 fold change), 
# it is recommended to use rarefying. As we observe 10 fold change in 
# the library sizes, we will normalize by rarefying all samples. 


p = ggrare(prunedSet, step = 20, color = "Sample_type", label = "Sample", se = FALSE) # Generate rarefaction curve plot 

p2 = p + geom_vline(xintercept=2530, color= "red", linetype='dashed') + # Add vertical red line and modifying axis limits 
  theme(text = element_text(size=12)) + 
  theme_minimal() +
  labs(color="Habitat") +
  scale_color_manual(values = c("#72B3DA", "royalblue4", "#000000","#A0522D"))

ps_rare=rarefy_even_depth(prunedSet, 2530,rngseed = 112, replace = FALSE, trimOTUs = TRUE, verbose = TRUE) # Rarefy 


# ------------------ Correlate rarefied and un-rarefied communities ------------------ 

asv_r=data.frame(otu_table(ps_rare))
asv=data.frame(otu_table(ps))

# Shannon diversity index correlation between rarefied and non rarefied samples
dShannon=ChaoShannon(asv)
dShannon$rar=diversity(t(asv_r))
# Species richness correlation between rarefied and non rarefied samples
Sr=ChaoRichness(asv)
Sr$rar=ChaoRichness(asv_r)$Observed

# ------------------ Plotting the results  ------------------ 

gshannon=ggplot(dShannon,aes(x=Estimator, y=rar))+
  stat_cor(method = "spearman", digits = 4, aes(label=paste(rr.label,..p.label..,sep='~`,`~'))) +
  geom_point()+
  theme_minimal()+xlab("Non-rarefied Shannon diversity")+ylab("Rarefied Shannon diversity") +
  geom_smooth(method = "lm", se = F, color = "#D6604D") +
  geom_abline(intercept = 0, slope = 1 , lty=2) +
  xlim(0,7) + ylim(0,7) +
  theme(text = element_text(size=12))


grichness <- ggplot(Sr,aes(x=Estimator, y=rar))+
  stat_cor(method = "spearman", digits = 4, aes(label=paste(rr.label,..p.label..,sep='~`,`~')), label.y=4680) +
  geom_point()+
  theme_minimal()+xlab("Non-rarefied ASV richness")+ylab("Rarefied ASV richness") +
  geom_smooth(method = "lm", se = F, color = "#D6604D") +
  geom_abline(intercept = 0, slope = 1 , lty=2) +
  xlim(0,5000) + ylim(0,5000) + 
  theme(text = element_text(size=12))


combined_plots=ggarrange(p2, ggarrange(gshannon, grichness, ncol=2, labels=c("b","c")), nrow=2, labels="a")


# ------------------ Save plot  ------------------ 

ggsave("rarefaction_curve.pdf", plot=combined_plots, width=16.5, height = 18,units="cm", path="/image/")


# ------------------ Save rarefied table  ------------------ 

write.csv(as.data.frame(as(tax_table(ps_rare), "matrix")), file = "/rarefied_taxa.csv")
write.csv(as.data.frame(as(otu_table(ps_rare), "matrix")),file = "/rarefied_asv.csv")
write.csv(as.data.frame(as(sample_data(ps_rare), "matrix")), file="/rarefied_meta.csv")
tree.raw <- phy_tree(ps_rare)
ape::write.tree(tree.raw , file = "/rarefied_tree.tree")
