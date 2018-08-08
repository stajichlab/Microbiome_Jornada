library(ggplot2)
library(ape)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(magrittr)
library(gridExtra)

plot_abundance = function(physeq, ylabn = "",
                          Facet = "Phylum",
                          Color = "Phylum"){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq,
         mapping = aes_string(x = "VegZone", y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + ylab(ylabn) +
    scale_y_log10()
}

meta <- read.table("050517NP515F-mapping2_NP.txt",
                   header=T,sep="\t", row.names=1,
                   stringsAsFactors=FALSE)
head(meta)
treefile = "Sam1_34a.tree.phy"
tree = read.tree(treefile)
otusall <- read.table("Sam1_34a.otu_table.taxonomy.txt",
                      header=T,sep="\t",row.names=1)

otus <- otusall[ , -which(names(otusall) %in% c("Taxonomy"))]
summary(otus)

#otus <- head(otus,n=500L)
otus <- as(as.matrix(otus), "matrix")
OTU = otu_table(otus, taxa_are_rows = TRUE)
sampleData <- sample_data(meta)

rownames(sampleData)
colnames(sampleData)

taxmat <- read.table("Sam1_34a.taxonomy.fix.tab",
                     header=T,sep="\t",row.names=1)
taxmat <- as(as.matrix(taxmat),"matrix")
TAX = tax_table(taxmat)
sample_names(OTU)
sample_names(sampleData)

#TAX <- head(TAX,n=120L)
phylo = phyloseq(OTU, TAX,sampleData,tree)
phylo

ps0 <- subset_taxa(phylo, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps1 = ps0

prevdf = apply(X = otu_table(ps1),
               MARGIN = ifelse(taxa_are_rows(ps1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps1),
                    tax_table(ps1))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))

prevalenceThreshold = 0.10 * nsamples(ps1)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps1)

length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))

ps3 = tax_glom(ps2, "Genus", NArm = TRUE)

total = median(sample_sums(ps3))
standf = function(x, t=total) round(t * (x / sum(x)))
ps3ra = transform_sample_counts(phylo, standf)

theme_set(theme_bw())

gpt <- ps3ra
gpt <- prune_taxa(names(sort(taxa_sums(ps3ra),TRUE)[1:2000]), ps3ra)
#plot_heatmap(gpt)

gpac <- subset_taxa(gpt,Class=="Cyanobacteria")

pdf("Abundance_HeatMap.pdf")
p <- plot_heatmap(gpac, "PCoA", "bray", "VegZone",low="#000033", high="#FF3300")
print(p)

p <- plot_heatmap(gpac, "PCoA", "bray", "VegZone","Genus",low="#000033", high="#FF3300")
print(p)
