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
         mapping = aes_string(x = "Crust_type", y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + ylab(ylabn) +
    scale_y_log10()
}

meta <- read.table("Biocrust_16S_Mapping_file.txt",
      header=T,sep="\t", row.names=1,
      stringsAsFactors=FALSE)
head(meta)
treefile = "BioCrusts2017Bacteria.tree.phy"
tree = read.tree(treefile)
otusall <- read.table("16S.otu_table.taxonomy.txt",
                   header=T,sep="\t",row.names=1)
meta <- meta[which(meta$Crust_type %in% c("LAC","CLC", "Sand")),]
meta <- meta[which(meta$Site %in% c("GMT","KELSO", "JTNP")),]
#otus <- otusall[ , which(names(otusall) %in%
#     c("NPFeb1SF","NPFeb3SF","NPFeb5SF","NPJan7SF","NPJan8SF","PD.11","PD.5",
#      "NPFeb2SF","NPFeb4SF","NPFeb6SF", "NPJan1SF","NPJan2SF", "NPJan3SF","PD.1",
#      "PD.7","NP45","NPFeb12","NP.2","N.11"))]

otus <- otusall[ , -which(names(otusall) %in% c("Taxonomy"))]
sum(otus)

#otus <- head(otus,n=500L)
otus <- as(as.matrix(otus), "matrix")
OTU = otu_table(otus, taxa_are_rows = TRUE)
sampleData <- sample_data(meta)

rownames(sampleData)
colnames(sampleData)

taxmat <- read.table("BioCrusts2017Bacteria.taxonomy_fix.tab",
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

prevalenceThreshold = 0.05 * nsamples(ps1)
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

xlabstr = "Crust type"
ggplot(phylum, aes(x = factor(Crust_type), y = Abundance, fill = Phylum)) +
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


ggplot(class, aes(x = factor(Crust_type), y = Abundance, fill = Class )) +
  geom_bar(stat = "identity") +
  #scale_fill_brewer(type = "seq", palette = "Paired") +
  #scale_colour_brewer(palette = "Set1") +
  theme( axis.text.x = element_text(angle = 60, hjust = 1)) +
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  xlab(xlabstr) + ylab("Relative Abundance (Class > 2%) \n") +
  ggtitle("Class Composition")

#plot_bar(phylostand, "Phylum", fill="Class",facet_grid=~Crust_type)
#plot_bar(phylostand, "Phylum", fill="Class",facet_grid=~Site)

xlabstr = "Sample"
ggplot(class, aes(x = Sample, y = Abundance, fill = Class,facet_grid=Crust_type)) +
  geom_bar(stat = "identity") +
  #scale_fill_brewer(type = "seq", palette = "Paired") +
  theme( axis.text.x = element_text(angle = 60, hjust = 1)) +
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    xlab(xlabstr) + ylab("Relative Abundance (Class > 2%) \n") +
  ggtitle("Class Composition ")

xlabstr = "Site"
ggplot(class, aes(x = factor(Site), y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
   #scale_fill_brewer(type = "seq", palette = "Paired") +
  # Remove x axis title
  theme( axis.text.x = element_text(angle = 60, hjust = 1)) +
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
    xlab(xlabstr) + ylab("Relative Abundance (Class > 2%) \n") +
  ggtitle("Class Composition")

xlabstr = "Month"
ggplot(class, aes(x = factor(Sampling_month), y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  #scale_fill_brewer(type = "seq", palette = "Paired") +
  # Remove x axis title
  theme( axis.text.x = element_text(angle = 60, hjust = 1)) +
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  xlab(xlabstr) + ylab("Relative Abundance (Class > 2%) \n") +
  ggtitle("Class Composition by Month")

xlabstr = "Humidity"
ggplot(class, aes(x = factor(Humidity), y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  #scale_fill_brewer(type = "seq", palette = "Paired") +
  # Remove x axis title
  theme( axis.text.x = element_text(angle = 60, hjust = 1)) +
  #
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  xlab(xlabstr) + ylab("Relative Abundance (Class > 2%) \n") +
  ggtitle("Class Composition by Humidity")

pdf("PCA_16S.pdf")
GP = ps3ra

clcsamp <- subset_samples(GP,Crust_type=="CLC")
clcrich <- estimate_richness(clcsamp)
summary(clcrich)
plot_richness(clcsamp,measures=c("Chao1","Shannon"),color="Site")

lacsamp <- subset_samples(GP,Crust_type=="LAC")
lacrich <- estimate_richness(lacsamp)
summary(lacrich)
plot_richness(lacsamp,measures=c("Chao1","Shannon"),color="Site")

plot_richness(GP, x="Crust_type", measures=c("Chao1","Shannon","Simpson"),color="Site")
plot_richness(GP,measures=c("Chao1","Shannon"),color="Site")
GP.ord <- ordinate(GP, "PCoA", "bray")
GP.ord$Site = factor(GP.ord$Site,levels=c("KELSO","GMT","CIMA","JTNP"))

p1 = plot_ordination(GP, GP.ord, type="taxa", color="Phylum", title="taxa")
p1

p2 = plot_ordination(GP, GP.ord, type="samples", color="Crust_type", shape="Site") 
p2 + geom_point(aes(fill=Crust_type)) + geom_point(size=3) + ggtitle("Crust Samples PCoA")

p2 = plot_ordination(GP, GP.ord, type="samples", color="Crust_type", shape="Site") 
p2 + geom_polygon(aes(fill=Crust_type)) + geom_point(size=3) + ggtitle("Crust Samples PCoA")

GP.ord <- ordinate(GP, "NMDS", "bray")
GP.ord$Site = factor(GP.ord$Site,levels=c("KELSO","GMT","JTNP"))
p1 = plot_ordination(GP, GP.ord, type="taxa", color="Phylum", title="taxa NMDS")
p1

p2 = plot_ordination(GP, GP.ord, type="samples", color="Crust_type", shape="Site") 
p2 + geom_point(aes(fill=Crust_type)) + geom_point(size=3) + ggtitle("Crust Samples NMDS") 

p2 = plot_ordination(GP, GP.ord, type="samples", color="Crust_type", shape="Site") 
p2 + geom_polygon(aes(fill=Crust_type)) + geom_point(size=3) + ggtitle("Crust Samples NMDS") 


GP.ord <- ordinate(GP, "PCoA", "unifrac",weighted=TRUE)
GP.ord$Site = factor(GP.ord$Site,levels=c("KELSO","GMT","JTNP"))
p1 = plot_ordination(GP, GP.ord, type="taxa", color="Phylum", title="taxa PCoA Unifrac")
p1


p2 = plot_ordination(GP, GP.ord, type="samples", color="Crust_type", shape="Site") 
p2 + geom_point(aes(fill=Crust_type)) + geom_point(size=3) + ggtitle("Crust Samples PCoA Unifrac")


p2 = plot_ordination(GP, GP.ord, type="samples", color="Crust_type", shape="Site") 
p2 + geom_polygon(aes(fill=Crust_type)) + geom_point(size=3) + ggtitle("Crust Samples PCoA Unifrac") 

pdf("16S_test_timeseries.pdf")
GP <- subset_samples(GP,Sampling_month != "HV.August_2017" & Sampling_month != "August_2016" )
newIn <- subset_samples(GP,Site=="JTNP" & Crust_type == "CLC")
newIn.ord <- ordinate(newIn, "PCoA", "bray")
p3 = plot_ordination(newIn, newIn.ord, type="samples", color="Sampling_month", shape="Rain_event") 
p3 + geom_point(aes(fill=Sampling_month)) + geom_point(size=3) + ggtitle("CLC Crust Samples PCoA Unifrac") 

newIn <- subset_samples(GP,Site=="JTNP" & Crust_type == "LAC")
newIn.ord <- ordinate(newIn, "PCoA", "bray")
p3 = plot_ordination(newIn, newIn.ord, type="samples", color="Sampling_month", shape="Rain_event") 
p3 + geom_point(aes(fill=Sampling_month)) + geom_point(size=3) + ggtitle("LAC Crust Samples PCoA Unifrac") 

newIn <- subset_samples(GP,Site=="JTNP" & Crust_type == "Sand")
newIn.ord <- ordinate(newIn, "PCoA", "bray")
p3 = plot_ordination(newIn, newIn.ord, type="samples", color="Sampling_month", shape="Rain_event") 
p3 + geom_point(aes(fill=Sampling_month)) + geom_point(size=3) + ggtitle("Sand Crust Samples PCoA Unifrac") 
