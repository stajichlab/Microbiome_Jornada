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
#meta <- meta[which(meta$Crust_type %in% c("LAC","CLC", "Sand")),]
#meta <- meta[which(meta$Site %in% c("GMT","KELSO", "JTNP")),]

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

pdf("Prevalence_Bacteria.pdf")
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps1),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

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

plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2, plotBefore, plotAfter)

phylum <- ps3ra %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "darkgreen","#DA5724", "#508578", "#CD9BCD",
  "orange","#AD6F3B", "#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "blue", "black","purple","yellow","magenta"

)

pdf("Taxoplot_Bacteria.pdf")

xlabstr = "VegZone"
ggplot(phylum, aes(x = factor(VegZone), y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  # Remove x axis title
  theme( axis.text.x = element_text(angle = 60, hjust = 1)) +
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  xlab(xlabstr) + ylab("Relative Abundance (Phyla > 0.5%) \n") +
  ggtitle("Phylum Composition")

class <- ps3ra %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at Class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class)                                      # Sort data frame alphabetically by phylum


ggplot(class, aes(x = factor(VegZone), y = Abundance, fill = Class )) +
  geom_bar(stat = "identity") +
  #scale_fill_brewer(type = "seq", palette = "Paired") +
  #scale_colour_brewer(palette = "Set1") +
  theme( axis.text.x = element_text(angle = 60, hjust = 1)) +
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  xlab(xlabstr) + ylab("Relative Abundance (Class > 2%) \n") +
  ggtitle("Class Composition")

#plot_bar(phylostand, "Phylum", fill="Class",facet_grid=~VegZone)
#plot_bar(phylostand, "Phylum", fill="Class",facet_grid=~Disturbance)

xlabstr = "Sample"
ggplot(class, aes(x = Sample, y = Abundance, fill = Class,facet_grid=VegZone)) +
  geom_bar(stat = "identity") +
  #scale_fill_brewer(type = "seq", palette = "Paired") +
  theme( axis.text.x = element_text(angle = 60, hjust = 1)) +
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    xlab(xlabstr) + ylab("Relative Abundance (Class > 2%) \n") +
  ggtitle("Class Composition ")

xlabstr = "Disturbance"
ggplot(class, aes(x = factor(Disturbance), y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
   #scale_fill_brewer(type = "seq", palette = "Paired") +
  # Remove x axis title
  theme( axis.text.x = element_text(angle = 60, hjust = 1)) +
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    xlab(xlabstr) + ylab("Relative Abundance (Class > 2%) \n") +
  ggtitle("Class Composition")


pdf("PCA_16S.pdf")
GP = ps3ra

tarsamp <- subset_samples(GP,VegZone=="Tarbush")
tarrich <- estimate_richness(tarsamp)
summary(tarrich)
plot_richness(tarsamp,measures=c("Chao1","Shannon"),color="Disturbance")

grasssamp <- subset_samples(GP,VegZone=="Grassland")
grassrich <- estimate_richness(grasssamp)
summary(grassrich)
plot_richness(grasssamp,measures=c("Chao1","Shannon"),color="Disturbance")

playasamp <- subset_samples(GP,VegZone=="Playa")
playarich <- estimate_richness(playasamp)
summary(playarich)
plot_richness(playasamp,measures=c("Chao1","Shannon"),color="Disturbance")


creosamp <- subset_samples(GP,VegZone=="Creosote")
creorich <- estimate_richness(creosamp)
summary(creorich)
plot_richness(creosamp,measures=c("Chao1","Shannon"),color="Disturbance")


plot_richness(GP, x="VegZone", measures=c("Chao1","Shannon","Simpson"),color="Disturbance")
plot_richness(GP,measures=c("Chao1","Shannon"),color="VegZone")

GP.ord <- ordinate(GP, "PCoA", "bray")
GP.ord$Site = factor(GP.ord$VegZone,levels=c("Tarbush","Grassland","Playa","Creosote"))

p1 = plot_ordination(GP, GP.ord, type="taxa", color="Phylum", title="taxa")
p1

p2 = plot_ordination(GP, GP.ord, type="samples", color="VegZone", shape="Disturbance") 
p2 + geom_point(aes(fill=VegZone)) + geom_point(size=3) + ggtitle("Jornada VegZone Samples PCoA")

p2 = plot_ordination(GP, GP.ord, type="samples", color="VegZone", shape="Disturbance") 
p2 + geom_polygon(aes(fill=VegZone)) + geom_point(size=3) + ggtitle("Jornada VegZone Samples PCoA")

GP.ord <- ordinate(GP, "NMDS", "bray")
GP.ord$Site = factor(GP.ord$VegZone,levels=c("Tarbush","Grassland","Playa","Creosote"))

p1 = plot_ordination(GP, GP.ord, type="taxa", color="Phylum", title="taxa")
p1

p2 = plot_ordination(GP, GP.ord, type="samples", color="VegZone", shape="Disturbance") 
p2 + geom_point(aes(fill=VegZone)) + geom_point(size=3) + ggtitle("Jornada VegZone Samples NMDS")

p2 = plot_ordination(GP, GP.ord, type="samples", color="VegZone", shape="Disturbance") 
p2 + geom_polygon(aes(fill=VegZone)) + geom_point(size=3) + ggtitle("Jornada VegZone Samples NMDS")


GP.ord <- ordinate(GP, "PCoA", "unifrac",weighted=TRUE)
GP.ord$Site = factor(GP.ord$VegZone,levels=c("Tarbush","Grassland","Playa","Creosote"))

p1 = plot_ordination(GP, GP.ord, type="taxa", color="Phylum", title="taxa")
p1

p2 = plot_ordination(GP, GP.ord, type="samples", color="VegZone", shape="Disturbance") 
p2 + geom_point(aes(fill=VegZone)) + geom_point(size=3) + ggtitle("Jornada VegZone Samples PCoA Unifrac")

p2 = plot_ordination(GP, GP.ord, type="samples", color="VegZone", shape="Disturbance") 
p2 + geom_polygon(aes(fill=VegZone)) + geom_point(size=3) + ggtitle("Jornada VegZone Samples PCoA Unifrac")