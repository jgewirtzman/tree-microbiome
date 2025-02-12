library(phyloseq)
library(RColorBrewer)
library(ggplot2)
library(microeco)
library(pheatmap)
library(ggpubr)

#### Import Metadata ####
#ddpcr = read.csv(file = "C:/Users/Wyatt/Desktop/UMGC Tree Full Sequences/ddPCR_tree_all_data.csv")
ddpcr = read.csv(file = "C:/Users/Wyatt/Desktop/UMGC Tree Full Sequences/ddPCR_meta_all_data.csv")
water = read.csv(file = "C:/Users/Wyatt/Desktop/UMGC Tree Full Sequences/Tree_Core_Sectioning_Data.csv")
qpcr_16s = read.csv(file = "C:/Users/Wyatt/Desktop/UMGC Tree Full Sequences/qPCR/16s_w_metadata.csv")
qpcr_16s=subset(qpcr_16s, qpcr_16s$Sample.ID!="None")
qpcr_16s=qpcr_16s[,c(3,4,6)]
qpcr_its = read.csv(file = "C:/Users/Wyatt/Desktop/UMGC Tree Full Sequences/qPCR/its_w_metadata.csv")
qpcr_its=subset(qpcr_its, qpcr_its$Sample.ID!="None")
qpcr_its=qpcr_its[,c(3,4,6)]


water$seq_id=toupper(water$seq_id)
ddpcr=merge(ddpcr, water, by = c("core_type","seq_id"), all.x = TRUE)
ddpcr=merge(ddpcr,qpcr_16s, by = c("core_type","seq_id"),  all.x = TRUE)
ddpcr=merge(ddpcr,qpcr_its, by = c("core_type","seq_id"),  all.x = TRUE)
#write.csv(ddpcr, "tree_methane_all_data.csv")

qpcr_bact=ggplot(data = subset(ddpcr, material=="Wood"), aes(x=core_type, y=log10((X16S_per_ul*75/Sample.Mass.Added.to.Tube..mg.*1000*10)+1), col=core_type))+
  stat_summary(size=1)+
  stat_summary(fun = "mean", geom = "text", aes(label = round(..y.., 2)), vjust = -1) + 
  geom_jitter()+
  ylab("Log(Copies/g)")+
  stat_compare_means(label = "p.signif")+
#  facet_grid(.~species, scales="free_x")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
qpcr_bact

ddpcr$core_type=factor(ddpcr$core_type, levels=c("Inner","Outer","Mineral","Organic"), labels = c("Heartwood","Sapwood","Mineral Soil","Organic Soil"))
qpcr_bact=ggplot(data = subset(ddpcr, material=="Wood"), aes(x=log10((X16S_per_ul*75/Sample.Mass.Added.to.Tube..mg.*1000*10)+1), y=species))+
  geom_density()+
  ggridges::geom_density_ridges(aes(fill=species), alpha=0.7, size=.1, rel_min_height = 0.01, scale=0.75)+
  xlab(expression(Log[10](Copies~g^-1)))+
  stat_summary(aes(col=core_type,fill=species), fun="mean", pch=21, size=1.1, alpha=0.6, stroke=1.5)+
  ylab("Tree Species")+
# facet_grid(.~core_type)+
  scale_color_manual(values=c("black","white"))+
  scale_fill_manual(values=species_colors)+
  theme_minimal(base_size = 16)
qpcr_bact


qpcr_bact = ggplot(data = subset(ddpcr, material == "Wood"), aes(x = log10((X16S_per_ul * 75 / Sample.Mass.Added.to.Tube..mg. * 1000 * 10) + 1), y = species)) +
  geom_density() +
  ggridges::geom_density_ridges_gradient(aes(fill = ..x..), alpha = 0.7, size = .1, rel_min_height = 0.001, scale = 1) +
  xlab(expression(Log[10](Copies~g^-1))) +
  stat_summary(aes(col = core_type, fill = ..x..), fun = "mean", pch = 21, size = 1.1, alpha = 0.6, stroke = 1.5) +
  ylab("Tree Species") +
  # facet_grid(.~core_type)+
  scale_color_manual(values = c("black", "white")) +
  scale_fill_viridis_c(option = "B", values = c(0.1,0.35,0.6,0.9)) +
  theme_minimal(base_size = 16)+
  theme(legend.position = "none")
qpcr_bact



qpcr_bact=ggplot(data = subset(ddpcr, material=="Wood"), aes(x=log10((X16S_per_ul*75/Sample.Mass.Added.to.Tube..mg.*1000*10)+1), y=core_type))+
  geom_density()+
  ggridges::geom_density_ridges(aes(fill=core_type), alpha=0.7, size=.1, rel_min_height = 0.01, scale=10)+
  xlab(expression(Log[10](Copies~g^-1)))+
  stat_summary(aes(col=core_type), fun="mean", pch=21, size=1)+
  stat_compare_means()+
  ylab("Tree Species")
  # facet_grid(.~core_type)+
qpcr_bact


bact_comp=ggplot(data=ddpcr, aes(x=log10(arc_loose+bact_loose+1), y=log10(X16S_per_ul*2+1)))+
  geom_point(aes(col=core_type))+
  geom_smooth(method="lm", se=FALSE, col="black")+
  stat_cor()+
  geom_abline(slope=1)+
  ylim(0,6)+
  xlim(0,6)+
  xlab("ddPCR No Plastid")+
  ylab("UMGC qPCR")
bact_comp


qpcr_fun=ggplot(data = subset(ddpcr, material=="Wood"), aes(x=core_type, y=log10((as.numeric(ITS_per_ul)*75*10)+1), col=core_type))+
  stat_summary(size=1)+
  geom_jitter()+
  ylab("Log(Copies/g)")+
  stat_compare_means(label = "p.signif")+
  facet_grid(.~species)
qpcr_fun

cowplot::plot_grid(qpcr_bact, qpcr_fun, nrow =2)

qpcr_fun=ggplot(data = subset(ddpcr, material=="Wood"), aes(x=seq_id, y=log10((as.numeric(ITS_per_ul)*75*10)+1), col=core_type))+
  geom_point()+
  ylab("Log(Copies/g)")+
  facet_grid(.~species, scales="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
qpcr_fun


#### 16S IMPORT #### 

#### NEPHELE into Phyloseq 
tax_tab=read.delim("D:/UMGC Full Tree Sequences/outputs/taxonomy_table.txt", header = TRUE)
otu_tab=read.delim("D:/UMGC Full Tree Sequences/outputs/OTU_table.txt", header=TRUE, row.names = 1)

bastard_tax=otu_tab[,590:596] # Subset taxonomy from end (for QIIME?) and make taxa table
bastard_tax[bastard_tax==""]=NA 
tax_tab_pre=tax_table(bastard_tax) 
taxa_names(tax_tab_pre)=sub("sp","seq",taxa_names(tax_tab_pre)) # tax_table conversts seq## to sp##, need to change back to merge

otu_tab_corr=otu_tab[,1:589] # Remove taxonomy from end of otu table
otu_table_pre=otu_table(otu_tab_corr, taxa_are_rows=TRUE)

phylo_tree=read_tree("D:/UMGC Full Tree Sequences/outputs/phylo/unrooted_tree.nwk")
samp_data=read.delim("D:/UMGC Full Tree Sequences/outputs/tree_16s_mapping_dada2_corrected.txt.no_gz", row.names = 1)
samp_data$RowName=row.names(samp_data)

fasta=read.delim("D:/UMGC Full Tree Sequences/outputs/seq.fasta")
#tree_biom=read_biom("D:/UMGC Full Tree Sequences/outputs/taxa.biom")

#### Merge Metadata 
samp_data$seq_id
samp_data$seq_id = sub("prime","'",samp_data$seq_id) # fix messed up names
samp_data$seq_id = sub("star","*",samp_data$seq_id)
samp_data$seq_id = sub("HM","H",samp_data$seq_id)
samp_data$seq_id[samp_data$ForwardFastqFile=="206_B01_16S_S3_R1_001.fastq"]="RO104"
samp_data$core_type[samp_data$ForwardFastqFile=="206_B01_16S_S3_R1_001.fastq"]="Inner"

samp_data_merged = merge(ddpcr, samp_data, by=c("seq_id","core_type"), all.y = TRUE)

dups=which(duplicated(samp_data_merged$RowName)==TRUE) # Remove duplicates, not sure where they come from
samp_data_merged=samp_data_merged[-c(dups),]

row.names(samp_data_merged)=samp_data_merged$RowName


#### Make PS object 
raw_ps = phyloseq(tax_tab_pre, otu_table_pre, phylo_tree, sam_data(samp_data_merged))


#### Remove Plastids 
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

mitochondria=rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[,5]=="Mitochondria")]
chloroplast=rownames(tax_table(raw_ps))[which(tax_table(raw_ps)[,4]=="Chloroplast")]

badTaxa = c(mitochondria, chloroplast)
no_mito = pop_taxa(raw_ps, badTaxa)


'
sdt = data.table(as(sample_data(no_mito), "data.frame"),
                 TotalReads = sample_sums(no_mito), keep.rownames = TRUE)
setnames(sdt, "rn", "Sano_mitoleID")
pSeqDepth = ggplot(sdt, aes(TotalReads, fill=material)) + geom_histogram() + ggtitle("Sequencing Depth") +geom_vline(xintercept =6000)
pSeqDepth


library(MicrobiotaProcess)
# for reproducibly random number
set.seed(46814)
rareres <- get_rarecurve(obj=no_mito, chunks=50)

p_rare <- ggrarecurve(obj=rareres,
                      indexNames=c("Observe","Chao1","ACE"),) +
  theme(legend.spacing.y=unit(0.01,"cm"),
        legend.text=element_text(size=4))
'



taxa_names(no_mito) <- paste0("ASV", seq(ntaxa(no_mito)))
set.seed(46814)
ps.rare = rarefy_even_depth(no_mito, sample.size = 3500)


#### Transform

ps.ra = transform_sample_counts(ps.rare, function(x) 100 * x/sum(x))
# Removes samples without matching metadata temporarily
ps.ra = subset_samples(ps.ra, core_type=="Inner" | core_type=="Outer" | core_type=="Mineral" | core_type=="Organic")
ps.unk = subset_samples(ps.ra, core_type!="Inner" & core_type!="Outer" & core_type!="Mineral" & core_type!="Organic")







#### ITS IMPORT #### 

tax_tab_its=read.delim("D:/UMGC Full Tree Sequences/ITS_results/outputs/taxonomy_table.txt", header = TRUE)
otu_tab_its=read.delim("D:/UMGC Full Tree Sequences/ITS_results/outputs/OTU_table.txt", header=TRUE, row.names = 1)

bastard_tax_its=otu_tab_its[,590:596] # Subset taxonomy from end (for QIIME?) and make taxa table
bastard_tax_its[bastard_tax_its==""]=NA 
tax_tab_pre_its=tax_table(bastard_tax_its) 
taxa_names(tax_tab_pre_its)=sub("sp","seq",taxa_names(tax_tab_pre_its)) # tax_table conversts seq## to sp##, need to change back to merge

otu_tab_corr_its=otu_tab_its[,1:589] # Remove taxonomy from end of otu table
otu_table_pre_its=otu_table(otu_tab_corr_its, taxa_are_rows=TRUE)

phylo_tree_its=read_tree("D:/UMGC Full Tree Sequences/ITS_results/outputs/phylo/rooted_tree.nwk")
samp_data_its=read.delim("D:/UMGC Full Tree Sequences/ITS_results/outputs/tree_its_mapping_dada2_corrected.txt.no_gz", row.names = 1)
samp_data_its$RowName=row.names(samp_data_its)


fasta_its=read.delim("D:/UMGC Full Tree Sequences/ITS_results/outputs/seq.fasta")
#tree_biom=read_biom("D:/UMGC Full Tree Sequences/outputs/taxa.biom")

#### Merge Metadata 
samp_data_its$seq_id
samp_data_its$seq_id = sub("prime","'",samp_data_its$seq_id) # fix messed up names
samp_data_its$seq_id = sub("star","*",samp_data_its$seq_id)
samp_data_its$seq_id = sub("HM","H",samp_data_its$seq_id)

colnames(samp_data_its)[3]="core_type"

# Not exactly sure why, but some duplicates emerged
samp_data_merged_its = merge(ddpcr, samp_data_its, by=c("seq_id","core_type"), all.y = TRUE)

dups=as.vector(which(duplicated(samp_data_merged_its$ForwardFastqFile)==TRUE)-1)
samp_data_merged_its=samp_data_merged_its[-c(dups),]

row.names(samp_data_merged_its)=samp_data_merged_its$RowName

#### Make PS object 
raw_ps_its = phyloseq(tax_tab_pre_its, otu_table_pre_its, phylo_tree_its, sample_data(samp_data_merged_its))

taxa_names(raw_ps_its) <- paste0("ASV", seq(ntaxa(raw_ps_its)))
set.seed(46814)
ps.its.rare = rarefy_even_depth(raw_ps_its, sample.size = 4000)

#rarecurve(t(as.data.frame(otu_table(raw_ps_its))), step=50, cex=0.5,label = FALSE)

ps.its.ra = transform_sample_counts(ps.its.rare, function(x) 100 * x/sum(x))
# Removes samples without matching metadata temporarily
ps.its.ra = subset_samples(ps.its.ra, core_type=="Inner" | core_type=="Outer" | core_type=="Mineral" | core_type=="Organic")

#### QUVE_16S ####

tax_tab_QUVE_16s=read.delim("D:/UMGC Full Tree Sequences/QUVE/16S/outputs/taxonomy_table.txt", header = TRUE)
otu_tab_QUVE_16s=read.delim("D:/UMGC Full Tree Sequences/QUVE/16S/outputs/OTU_table.txt", header=TRUE, row.names = 1)

bastard_tax_QUVE_16s=otu_tab_QUVE_16s[,61:67] # Subset taxonomy from end (for QIIME?) and make taxa table
bastard_tax_QUVE_16s[bastard_tax_QUVE_16s==""]=NA 
tax_tab_pre_QUVE_16s=tax_table(bastard_tax_QUVE_16s) 
taxa_names(tax_tab_pre_QUVE_16s)=sub("sp","seq",taxa_names(tax_tab_pre_QUVE_16s)) # tax_table converts seq## to sp##, need to change back to merge

otu_tab_corr_QUVE_16s=otu_tab_QUVE_16s[,1:60] # Remove taxonomy from end of otu table
otu_table_pre_QUVE_16s=otu_table(otu_tab_corr_QUVE_16s, taxa_are_rows=TRUE)

phylo_tree_QUVE_16s=read_tree("D:/UMGC Full Tree Sequences/QUVE/16S/outputs/phylo/rooted_tree.nwk")
samp_data_QUVE_16s=read.delim("D:/UMGC Full Tree Sequences/QUVE/16S/outputs/black_oak_mapping_corrected.txt.no_gz", row.names = 1)
samp_data_QUVE_16s$RowName=row.names(samp_data_QUVE_16s)


fasta_QUVE_16s=read.delim("D:/UMGC Full Tree Sequences/QUVE/16S/outputs/seq.fasta")
#tree_biom=read_biom("D:/UMGC Full Tree Sequences/outputs/taxa.biom")

raw_ps_QUVE_16s = phyloseq(tax_tab_pre_QUVE_16s, otu_table_pre_QUVE_16s, phylo_tree_QUVE_16s, sample_data(samp_data_QUVE_16s))

taxa_names(raw_ps_QUVE_16s) <- paste0("ASV", seq(ntaxa(raw_ps_QUVE_16s)))
set.seed(46814)

mitochondria=rownames(tax_table(raw_ps_QUVE_16s))[which(tax_table(raw_ps_QUVE_16s)[,5]=="Mitochondria")]
chloroplast=rownames(tax_table(raw_ps_QUVE_16s ))[which(tax_table(raw_ps_QUVE_16s )[,4]=="Chloroplast")]

badTaxa = c(mitochondria, chloroplast)
no_mito_q16 = pop_taxa(raw_ps_QUVE_16s, badTaxa)

ps.quve.16s.rare = rarefy_even_depth(no_mito_q16, sample.size = 4000)
ps.quve.16s.ra = transform_sample_counts(ps.quve.16s.rare, function(x) 100 * x/sum(x))

rarecurve(t(as.data.frame(otu_table(no_mito_q16))), step=50, cex=0.5)

#### QUVE_ITS ####
tax_tab_QUVE_its=read.delim("/Users/Wyatt/Desktop/UMGC Tree Full Sequences/NEPHELE_Outputs/outputs/taxonomy_table.txt", header = TRUE)
otu_tab_QUVE_its=read.delim("/Users/Wyatt/Desktop/UMGC Tree Full Sequences/NEPHELE_Outputs/outputs/OTU_table.txt", header=TRUE, row.names = 1)

bastard_tax_QUVE_its=otu_tab_QUVE_its[,61:67] # Subset taxonomy from end (for QIIME?) and make taxa table
bastard_tax_QUVE_its[bastard_tax_QUVE_its==""]=NA 
tax_tab_pre_QUVE_its=tax_table(bastard_tax_QUVE_its) 
taxa_names(tax_tab_pre_QUVE_its)=sub("sp","seq",taxa_names(tax_tab_pre_QUVE_its)) # tax_table conversts seq## to sp##, need to change back to merge

otu_tab_corr_QUVE_its=otu_tab_QUVE_its[,1:60] # Remove taxonomy from end of otu table
otu_table_pre_QUVE_its=otu_table(otu_tab_corr_QUVE_its, taxa_are_rows=TRUE)

phylo_tree_QUVE_its=read_tree("/Users/Wyatt/Desktop/UMGC Tree Full Sequences/NEPHELE_Outputs/outputs/phylo/unrooted_tree.nwk")
samp_data_QUVE_its=read.delim("/Users/Wyatt/Desktop/UMGC Tree Full Sequences/NEPHELE_Outputs/outputs/black_oak_mapping_its_corrected.txt.no_gz", row.names = 1)
samp_data_QUVE_its$RowName=row.names(samp_data_QUVE_its)


fasta_QUVE_its=read.delim("/Users/Wyatt/Desktop/UMGC Tree Full Sequences/NEPHELE_Outputs/outputs/seq.fasta")
#tree_biom=read_biom("D:/UMGC Full Tree Sequences/outputs/taxa.biom")

raw_ps_QUVE_its = phyloseq(tax_tab_pre_QUVE_its, otu_table_pre_QUVE_its, phylo_tree_QUVE_its, sample_data(samp_data_QUVE_its))

taxa_names(raw_ps_QUVE_its) <- paste0("ASV", seq(ntaxa(raw_ps_QUVE_its)))
set.seed(46814)
ps.quve.its.rare = rarefy_even_depth(raw_ps_QUVE_its, sample.size = 4000)

ps.quve.its.ra = transform_sample_counts(ps.quve.its.rare, function(x) 100 * x/sum(x))

#### Color Palettes ####
library("RColorBrewer")

coul <- c(brewer.pal(9, "Set3"),brewer.pal(12, "Paired"),rev(brewer.pal(8, "Accent")),brewer.pal(8, "Dark2"),
          brewer.pal(12, "Paired"), brewer.pal(8, "Accent"), brewer.pal(8, "Dark2"))

coul_army <- c("#4d784e","#7C887E", "#6ea171","#91967F", "#e1d798","#BAB79F","#fedc00",
               "#675645","#D3D2CD", "#424756","#A08D83","#726870","#A08D83" ,"#ad1e22",
               "#E89C31","#083248","#DBA858","#ffffff",
               "#4d784e","#7C887E", "#6ea171","#91967F", "#e1d798","#BAB79F",
               "#675645","#ad1e22","#D3D2CD", "#424756","#8C0E0F","#726870", "#8C0E0F",
               "#E89C31","#083248","#DBA858","#ffffff")


color_taxa=function(ps.wood, level){
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  speciesList = unique(tax_table(ps.wood)[,level])
  speciesPalette = getPalette(length(speciesList))
  names(speciesPalette) = speciesList
  return(speciesPalette)
}

color_taxa_meco=function(meco_bacteria, level){
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  speciesList = unique(meco_bacteria$tax_table[,level])
  speciesPalette = getPalette(length(speciesList))
  names(speciesPalette) = speciesList
  return(speciesPalette)
}


assign_species_colors <- function(df) {
  # Combine the color palettes from RColorBrewer
  combined_palette <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))
  
  # Make sure we have enough colors, recycle if needed
  unique_species <- unique(df$species)
  color_count <- length(unique_species)
  
  # If not enough colors, repeat the palette
  if (color_count > length(combined_palette)) {
    combined_palette <- rep(combined_palette, ceiling(color_count / length(combined_palette)))
  }
  
  # Assign a color to each species
  species_colors <- setNames(combined_palette[1:color_count], unique_species)
  
  # Add a color column to the dataframe based on species
  df$color <- species_colors[df$species]
  
  return(df)
}





# Ecological Niches ----

# * ITS ----
ps.its.wood=subset_samples(ps.its.ra)


colnames(tax_table(ps.its.wood))=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")


meco_fungi <- microtable$new(otu_table = as.data.frame(otu_table(ps.its.wood)), 
                             tax_table = noquote(as.data.frame(tax_table(ps.its.wood))),
                             sample_table = as.data.frame(as.matrix(sample_data(ps.its.wood))),
                             phylo_tree = phy_tree(ps.its.wood))


t2 = trans_func$new(meco_fungi)
t2$cal_spe_func(fungi_database = "fungal.traits")
t2$cal_spe_func(fungi_database = "FUNGuild")
FT=as.data.frame(t2$res_spe_func_raw_FungalTraits)
FG=as.data.frame(t2$res_spe_func_raw_funguild)


t2$cal_spe_func_perc(abundance_weighted = FALSE, dec=2)


# * 16S ----

ps.wood=subset_samples(ps.ra)

colnames(tax_table(ps.wood))=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_table(ps.wood)[,1]=paste("k__",tax_table(ps.wood)[,1], sep="")
tax_table(ps.wood)[,2]=paste("p__",tax_table(ps.wood)[,2], sep="")
tax_table(ps.wood)[,3]=paste("c__",tax_table(ps.wood)[,3], sep="")
tax_table(ps.wood)[,4]=paste("o__",tax_table(ps.wood)[,4], sep="")
tax_table(ps.wood)[,5]=paste("f__",tax_table(ps.wood)[,5], sep="")
tax_table(ps.wood)[,6]=paste("g__",tax_table(ps.wood)[,6], sep="")
tax_table(ps.wood)[,7]=paste("s__",tax_table(ps.wood)[,7], sep="")

sample_data(ps.wood)$core_type=factor(sample_data(ps.wood)$core_type, levels=c("Inner","Outer","Mineral","Organic"))

meco_bacteria <- microtable$new(otu_table = as.data.frame(otu_table(ps.wood)), 
                             tax_table = noquote(as.data.frame(tax_table(ps.wood))),
                             sample_table = as.data.frame(as.matrix(sample_data(ps.wood))),
                             phylo_tree = phy_tree(ps.wood))


t16 = trans_func$new(meco_bacteria)
t16$cal_spe_func(prok_database = "FAPROTAX")
FTX=as.data.frame(t16$res_spe_func)
t16$cal_spe_func_perc(abundance_weighted = FALSE, dec=2)

# * QUVE 16S ####
ps.quve.16s.ra=subset_samples(ps.quve.16s.ra, core_type!="NA")

colnames(tax_table(ps.quve.16s.ra))=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_table(ps.quve.16s.ra)[,1]=paste("k__",tax_table(ps.quve.16s.ra)[,1], sep="")
tax_table(ps.quve.16s.ra)[,2]=paste("p__",tax_table(ps.quve.16s.ra)[,2], sep="")
tax_table(ps.quve.16s.ra)[,3]=paste("c__",tax_table(ps.quve.16s.ra)[,3], sep="")
tax_table(ps.quve.16s.ra)[,4]=paste("o__",tax_table(ps.quve.16s.ra)[,4], sep="")
tax_table(ps.quve.16s.ra)[,5]=paste("f__",tax_table(ps.quve.16s.ra)[,5], sep="")
tax_table(ps.quve.16s.ra)[,6]=paste("g__",tax_table(ps.quve.16s.ra)[,6], sep="")
tax_table(ps.quve.16s.ra)[,7]=paste("s__",tax_table(ps.quve.16s.ra)[,7], sep="")

meco_quve_16s <- microtable$new(otu_table = as.data.frame(otu_table(ps.quve.16s.ra)), 
                                tax_table = noquote(as.data.frame(tax_table(ps.quve.16s.ra))),
                                sample_table = as.data.frame(as.matrix(sample_data(ps.quve.16s.ra))),
                                phylo_tree = phy_tree(ps.quve.16s.ra))

#### * QUVE ITS ####
ps.quve.its.ra=subset_samples(ps.quve.its.ra, core_type!="NA")

colnames(tax_table(ps.quve.its.ra))=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
meco_quve_its <- microtable$new(otu_table = as.data.frame(otu_table(ps.quve.its.ra)), 
                             tax_table = noquote(as.data.frame(tax_table(ps.quve.its.ra))),
                             sample_table = as.data.frame(as.matrix(sample_data(ps.quve.its.ra))),
                             phylo_tree = phy_tree(ps.quve.its.ra))


#### Functional  ####

# * ITS ----
#### ```` Metabolisms ####

ps.its.g=tax_glom(ps.its.ra, taxrank=rank_names(ps.its.ra)[6], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
ps.its.g.filt=prune_taxa(taxa_sums(ps.its.g)>110, ps.its.g)

#ps.its.g.filt.wood=subset_samples(ps.its.g.filt, material.x=="Wood")
#no_mito_select=subset_samples(ps.its.g.filt.wood, species=="TSCA"|species=="PIST"|species=="ACSA"|species=="ACRU")
no_mito_select=ps.its.g.filt

# Format FungalTraits data as tax_table to replace existing

FT=replace(FT, FT=="", NA)
FT_tt=tax_table(FT)
row.names(FT_tt)=row.names(FT)

FG=replace(FG, FG=="", NA)
FG_tt=tax_table(FG)
row.names(FG_tt)=row.names(FG)


tax_table(no_mito_select)=FT_tt
colnames(tax_table(no_mito_select))=colnames(FT)

tax_table(no_mito_select)=FG_tt
colnames(tax_table(no_mito_select))=colnames(FG)


no_mito_select_light=subset_samples(no_mito_select, species.x!="QUAL" & species.x!="QUVE" & species.x!="KALA"  & species.x!="PRSE"
                                    & species.x!="CAOV" & species.x!="SAAL")

no_mito_select_light=subset_samples(no_mito_select, material.x=="Wood")

fp=plot_bar(no_mito_select_light, x="seq_id", fill ="primary_lifestyle")+
  geom_bar(stat="identity")+scale_fill_manual(values=coul[12:60])+
  facet_grid(core_type~., scales="free", space = "free") + ylim(0,100) + theme(legend.position = "bottom") + theme(legend.text=element_text(size=12))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
fp


# ```` Alpha ----
ps.its.rare.clean=subset_samples(ps.its.rare, core_type=="Inner" | core_type=="Outer" | core_type=="Mineral" | core_type=="Organic")
alpha_its=plot_richness(ps.its.rare.clean, x="species.x", measures=c("Chao1", "Shannon"), col="core_type", title = "ITS")+
  stat_summary(size=1.2)
alpha_its


# ```` Heatmaps ----
t2$sample_table$core_type=factor(t2$sample_table$core_type, levels=c("Inner","Outer","Mineral","Organic"))

t2$plot_spe_func_perc(add_facet = TRUE)+
  facet_grid(.~t2$sample_table$core_type, scales="free")+
  theme(axis.text.x=element_blank())+
  scale_fill_viridis_c()+
  xlab("Sample Type")

meco_fungi_all_wood <- clone(meco_fungi)
meco_fungi_all_wood$sample_table <- subset(meco_fungi_all_wood $sample_table,  core_type == "Inner" | core_type == "Outer")
meco_fungi_all_wood$tidy_dataset()

fung_core_avg=meco_fungi_all_wood$merge_samples(use_group = "core_type")

tf_c_avg = trans_func$new(fung_core_avg)
tf_c_avg$cal_spe_func(fungi_database = "FUNGuild")
tf_c_avg$cal_spe_func_perc(abundance_weighted = TRUE, dec=2)

tf_c_avg$res_spe_func_perc=tf_c_avg$res_spe_func_perc[c(7,8,11,18,31,44,45,46,51,52,52,54,60,64,66,69,71,72)]

row.names(tf_c_avg$res_spe_func_perc)=c("Heartwood", "Sapwood")

tf_c_avg$plot_spe_func_perc(add_facet = TRUE)+
  #  facet_grid(.~t16_s$sample_table$SampleID, scales="free")+
  geom_text(aes(label = round(value, 2)), color = "white", size = 5) +  # Add labels (numbers) to the heatmap
  scale_fill_gradientn(colors = c("#440154","#3b528b","#21918c", "#5ec962", "#fde725"), na.value = "grey",
                       trans = "log10",  # Logarithmic transformation
                       limits = c(0.01, NA)) +  # Automatically scale to the data range
  theme(axis.text.y = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))



# Environmental heatmaps
meco_fungi_all_wood$sample_table$Percent_Change[(meco_fungi_all_wood$sample_table$Percent_Change)<0]=NA
meco_fungi_all_wood$sample_table$mcra_probe_loose=log10(as.numeric(meco_fungi_all_wood$sample_table$mcra_probe_loose)+1)
meco_fungi_all_wood$sample_table$pmoa_strict=log10(as.numeric(meco_fungi_all_wood$sample_table$pmoa_strict)+1)

t3 <- trans_env$new(dataset = meco_fungi_all_wood, env_cols = 
                      c("mcra_probe_loose", "pmoa_strict", "mmox_strict","ch4_50","ch4_125","ch4_200",
                        "CH4_int", "CO2_int","N2O_int", "O2_int", "co2_50","co2_125","co2_200", "mean_vwc",
                        "mean_st","X16S_per_ul","Percent_Change"), complete_na = TRUE, character2numeric = TRUE)
t3$cal_cor(add_abund_table = tf_c_avg$res_spe_func_perc, cor_method = "spearman", 
           by_group = "core_type")
t3$plot_cor(filter_feature = c("","*"))

t3$cal_cor(use_data = "ta6", p_adjust_method = "none", use_taxa_num = 500, 
           by_group = "core_type")
t3$plot_cor(filter_feature = c("","*","**"))


t3$cal_cor(add_abund_table = tf_c_avg$res_spe_func_perc, cor_method = "spearman", 
           by_group = "core_type")
t3$plot_cor()

t4 <- trans_env$new(dataset = meco_fungi_all_wood, add_data = tf_c_avg$res_spe_func_perc[4:20], complete_na = TRUE, character2numeric = TRUE)
t4$cal_cor(use_data = "Genus", p_adjust_method = "none", use_taxa_num = 100)
t4$plot_cor(filter_feature = c("","*"))

t3$plot_scatterfit(x=tf_c_avg$res_spe_func_perc$Pathotroph, y=log10(t3$data_env$X16S_per_ul+1))


# Taxa heatmap

t3_f <- trans_abund$new(dataset = fung_core_avg, taxrank = "Class", ntaxa = 20)
t3_f$plot_heatmap(facet =fung_core_avg$SampleID ,xtext_keep = FALSE, withmargin = FALSE) 
#  facet_grid(~meco_fungi$sample_table$core_type.y, scales="free")

meco_fungi_all <- clone(meco_fungi)
meco_fungi_all$sample_table <- subset(meco_fungi$sample_table,  core_type == "Inner" | core_type == "Outer" |  core_type == "Mineral" | core_type == "Organic")
meco_fungi_all$tidy_dataset()

meco_fung_merge=meco_fungi_all$merge_samples(use_group = "core_type")
#meco_fung_merge$sample_table$SampleID=factor(meco_fung_merge$sample_table$SampleID, levels = c("Inner","Outer","Mineral","Organic"))
# Create a trans_abund object with the merged dataset
t3_ff <- trans_abund$new(dataset = meco_fung_merge, taxrank = "Class", ntaxa = 23)
t3_ff$data_abund$Sample=factor(t3_ff$data_abund$Sample, 
                            levels = c("Inner", "Outer", "Mineral", "Organic"), 
                            labels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))

# Plot the heatmap with the updated x-axis labels
t3_ff$plot_heatmap(xtext_keep = TRUE, withmargin = TRUE, xtext_size = 16, ytext_size = 16)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "right")+
  theme(legend.title = element_text(size=12))+
  theme(legend.title = element_text())+
  theme(legend.text = element_text(size=12))

#### ```` Abundance ----

meco_fungi_all_wood$sample_table$core_type=factor(meco_fungi_all_wood$sample_table$core_type, levels=c("Inner","Outer","Mineral","Organic"),
                                         labels=c("Heartwood","Sapwood","Soil, Mineral","Soil, Organic"))


t4_f <- trans_abund$new(dataset = meco_fungi, taxrank = "Class", ntaxa = 10, group ="core_type")
t4_f$plot_pie(facet_nrow = 1)+scale_fill_manual(values=coul[12:40])


dataset_merge_f <- meco_fungi_all_wood$merge_samples(use_group = "core_type")
t5_f <- trans_venn$new(dataset_merge_f, ratio = "numratio")
t5_f$plot_venn()



#### ```` Diff ----
t6_f <- trans_diff$new(dataset = meco_fungi_all_wood, method = "lefse", group = "core_type", alpha = 0.01, lefse_subgroup = NULL, 
                       taxa_level = "Class", rm_un=TRUE)

t6_f$plot_diff_bar(threshold = 3.5)
t6_f$res_abund$Group=factor(t6_f$res_abund$Group, 
                          levels = c("Inner", "Outer", "Mineral", "Organic"), 
                          labels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))
t6_f$res_diff$Group=factor(t6_f$res_diff$Group, 
                         levels = c("Inner", "Outer", "Mineral", "Organic"), 
                         labels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))

t6_f$res_diff <- t6_f$res_diff %>%
  filter(!grepl("\\|[^|]*\\|NA", Taxa))

# Print the updated dataframe
#t6_f$res_diff[[2,2]]

t6_f$plot_diff_bar(use_number = 1:30, width = 0.8, group_order = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"), keep_prefix = FALSE, axis_text_y = 14)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  scale_color_manual(values=custom_palette)+
  scale_fill_manual(values=custom_palette)

t6_f$plot_diff_cladogram(use_taxa_num = 100, use_feature_num = 50, clade_label_level = 3)

#### ```` Environment ----

t7_f <- trans_env$new(dataset = meco_fungi, env_cols=c(35:47), complete_na = TRUE)
t7_f <- trans_env$new(dataset = meco_fungi, env_cols=c(22,31,33), complete_na = TRUE)

t7_f$cal_cor(by_group = "core_type", use_data = "other", p_adjust_method = "fdr", other_taxa = t6_f$res_diff$Taxa[1:40])
t7_f$plot_cor()

t7_f$cal_cor(by_group = "core_type",add_abund_table = t2$res_spe_func_perc, cor_method = "spearman")
t7_f$plot_cor()


t7_f$cal_mantel(by_group="core_type", use_measure = "bray")
t7_f$res_mantel
#### ```` Beta ----
set.seed(46814)
t7_ff <- trans_env$new(dataset = meco_fungi_all_wood,  env_cols=c(19, 22, 44, 38, 40, 42), complete_na = TRUE) # ,add_data =tf_c_avg$res_spe_func_perc[4:20])
t7_ff$cal_ordination(method = "CCA", use_measure = "bray", taxa_level = "ta2")
# t1$res_rda is the result list stored in the object
t7_ff$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)
# t1$res_rda_trans is the transformed result for plotting
t7_ff$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
t7_ff$plot_ordination(plot_color = "core_type")

t7_ff$cal_ordination_anova()
t7_ff$res_ordination_terms

t7_ff$cal_ordination_envfit()
t7_ff$res_ordination_envfit

t7_ff$cal_mantel(use_measure = "bray", by_group ="core_type")
t7_ff$res_mantel




# Beta diversity within species
meco_fungi_wood <- clone(meco_fungi)
meco_fungi_wood$sample_table <- subset(meco_fungi_wood$sample_table, species.x != "NA" & core_type == "Outer")
meco_fungi_wood$sample_table <- subset(meco_fungi_wood$sample_table, species.x!="KALA"& species.x!="SAAL" &
                                                                      species.x!="QUAL" &
                                                                      species.x!="QUVE" &
                                                                      species.x!="CAOV" &
                                                                      species.x!="PRSE")
meco_fungi_wood$tidy_dataset()
meco_fungi_wood$cal_betadiv(unifrac = TRUE)


t10_its <- trans_beta$new(dataset = meco_fungi_wood, group = "species.x", measure="unwei_unifrac")
t10_its$cal_group_distance(within_group = TRUE)
#t10$cal_group_distance_diff(method = "anova")
t10_its$plot_group_distance(boxplot_add = "mean")+ylim(0.5,1.1)


# Beta-div phylogeny

t2_w_its <- trans_beta$new(dataset = meco_fungi_wood, measure = "bray", group = "species")
t2_w_its$plot_clustering(group = "core_type", replace_name = c("species.x"), color_values = coul[9:30])

meco_fungi_group <- meco_fungi_wood$merge_samples(use_group = "species.x")
meco_fungi_group$cal_betadiv(unifrac=TRUE)
t3_its_group <- trans_beta$new(dataset = meco_fungi_group, measure = "wei_unifrac")

t3_its_group$plot_clustering(group = "SampleID", replace_name = c("SampleID"), color_values = c("black","black","red","red","red",
                                                                                                "black","black","black","red","black",
                                                                                                "black","black","black","black","red","black"))+
  ggtitle("ITS Sapwood, Weighted Unifrac")
# Phylogeny above

t3_its <- trans_abund$new(dataset = meco_fungi_group, taxrank = "Class", ntaxa = 8, groupmean = "SampleID")

t3_its$plot_radar(values.radar = c("0%", "25%", "50%"), grid.min = 0, grid.mid = 0.25, grid.max = 1)+
  ggtitle("Inner Wood, ITS")

t3_its$plot_donut(label = FALSE, facet_nrow = 4)+
  ggtitle("Inner Wood, ITS")


# Ordination
meco_fungi_wood_both <- clone(meco_fungi)
meco_fungi_wood_both$sample_table <- subset(meco_fungi_wood_both$sample_table, material.x == "Wood" | material.x == "Soil")
meco_fungi_wood_both$tidy_dataset()
meco_fungi_wood_both$cal_betadiv(unifrac = TRUE)


t8_f <- trans_beta$new(dataset = meco_fungi_wood_both, measure = "wei_unifrac")
t8_f$cal_ordination(ordination = "PCoA")
t8_f$res_ordination$scores$core_type=factor(t8_f$res_ordination$scores$core_type, 
                                            levels = c("Inner", "Outer", "Mineral", "Organic"), 
                                            labels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))

custom_palette <- rev(c("#018571", "#80cdc1", "#dfc27d", "#a6611a"))
custom_palette2 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")


t8_f$plot_ordination(plot_color = "core_type", plot_type = c("point","centroid"), centroid_segment_alpha = 0.25)+
  ggtitle("PCoA, WUniFrac: ITS")+
  geom_point(size=4, alpha=0.75)+
  scale_color_manual(values = custom_palette, name="Sample Type")+
  theme_minimal(base_size = 16)


#### Tests

meco_fungi_adonis <- clone(meco_fungi)
meco_fungi_adonis$sample_table <- subset(meco_fungi_adonis$sample_table)
meco_fungi_adonis$tidy_dataset()
meco_fungi_adonis$cal_betadiv(unifrac=TRUE)

t_f <- trans_beta$new(dataset = meco_fungi_adonis, measure = "wei_unifrac")
t_f$sample_table$material.x[is.na(t_f$sample_table$material.x)]="MISC"
t_f$cal_manova(group = "material.x")
t_f$res_manova



#### ```` FEAST ####
ps.g.its.count=tax_glom(ps.its.rare, taxrank=rank_names(ps.its.rare)[5], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
feast_file_its = subset_samples(ps.g.its.count, core_type=="Inner" | core_type=="Outer" | core_type=="Mineral" | core_type=="Organic")
#Should glom to speed up and reduce unknowns

otus_its=as.matrix(as.data.frame(otu_table(feast_file_its)))
metadata_its=as.data.frame(matrix(NA,ncol=4))


### Unpaired Sink/Source Route ###
SampleID=colnames(otus_its)
Env=sample_data(feast_file_its)$core_type
SourceSink=ifelse(sample_data(feast_file_its)$core_type=="Inner","Sink","Source") # Specify here 

count_occurrences <- function(lst) {
  counts <- c()
  seen <- list()
  
  for (value in lst) {
    if (is.null(seen[[value]])) {
      seen[[value]] <- 1
      counts <- c(counts, 1)
    } else {
      seen[[value]] <- seen[[value]] + 1
      counts <- c(counts, seen[[value]])
    }
  }
  
  return(counts)
}
id=count_occurrences(Env)
id[SourceSink=="Source"]=NA


feast_input=as.data.frame(cbind(Env, SourceSink, id))
feast_input$id=as.numeric(feast_input$id)
rownames(feast_input)=SampleID


### Paired Sink/Source Route ###
SampleID=colnames(otus_its)
Env=sample_data(feast_file_its)$core_type
SourceSink=ifelse(sample_data(feast_file_its)$core_type=="Inner","Sink","Source") # Specify here 
Pair_key=sample_data(feast_file_its)$seq_id
feast_input=as.data.frame(cbind(Env, SourceSink, id, Pair_key))


feast_input <- feast_input %>%
  arrange(Pair_key) %>%
  group_by(Pair_key, SourceSink) %>%
  mutate(Group_number = cur_group_id()) %>%
  ungroup() %>%
  group_by(Pair_key) %>%
  mutate(Has_Source = "Source" %in% SourceSink,
         Has_Sink = "Sink" %in% SourceSink,
         Group_number = case_when(Has_Source & Has_Sink ~ Group_number[SourceSink == "Source"][1],
                                  TRUE ~ NA)) %>%
  select(-Has_Source, -Has_Sink)

feast_input$id=as.numeric(feast_input$Group_number)
feast_input=as.data.frame(feast_input)
feast_input=feast_input[,c(1:3)]

#SampleID originally, now species
rownames(feast_input)=SampleID

feast_input=subset(feast_input, is.na(id)==FALSE)
#feast_input=feast_input[c(1:305,315:420),] # Something weird here for sapwood, needs to have last few cut
feast_input=feast_input[c(1:377,379:385),] # Something weird here for heartwood, needs to have last few cut
# Convert lists to vectors
feast_input$Env <- sapply(feast_input$Env, function(x) x[[1]])
feast_input$SourceSink <- sapply(feast_input$SourceSink, function(x) x[[1]])


### Run FEAST


write.table(feast_input, "/Users/Wyatt/Desktop/FEAST/feast_input_paired_its.txt", sep = '\t')
innyy=Load_metadata("/Users/Wyatt/Desktop/FEAST/feast_input_paired_its.txt")

write.table(otus_its, "/Users/Wyatt/Desktop/FEAST/feast_otu_paired_its.txt", sep = '\t')
otuu=Load_CountMatrix("/Users/Wyatt/Desktop/FEAST/feast_otu_paired_its.txt")

#innyy=Load_metadata("/Users/Wyatt/Desktop/metadata_example_multi.txt")
#otuu=Load_CountMatrix("/Users/Wyatt/Desktop/otu_example_multi.txt")

# Different sources flag changed for paired (1) /unpaired (0)
FEAST_output <- FEAST(C = otuu, metadata = innyy, different_sources_flag=1, dir_path = "/Users/Wyatt/Desktop/FEAST/",
                      outfile="its_Outer_paired")

FEAST_tab=read.delim("/Users/Wyatt/Desktop/FEAST/its_Outer_paired_source_contributions_matrix.txt")
colnames(FEAST_tab)=subset(feast_input, Env!="Inner")$Env
colnames(FEAST_tab)[is.na(colnames(FEAST_tab))]="Unknown"

### Plot FEAST

min_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="Mineral")], na.rm = TRUE)
org_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="Organic")], na.rm = TRUE)
unk_sums=(FEAST_tab[,which(colnames(FEAST_tab)=="Unknown")])
out_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="Outer")], na.rm = TRUE)

feast_sum=cbind(min_sums, org_sums, unk_sums, out_sums)
colnames(feast_sum)=c("Mineral","Organic","Unknown","Outer")
feast_sum=as.data.frame(feast_sum)

library(tidyr)

# Add the row names as a separate column
feast_sum$Row_Name <- rownames(feast_sum)

# Transform the dataframe to have two columns: 'Row_Name', 'Value', and 'Column_Name'
transformed_dataframe <- pivot_longer(feast_sum, 
                                      cols = c("Mineral", "Outer", "Organic", "Unknown"), 
                                      names_to = "Column_Name", 
                                      values_to = "Value")

# Set the row names back to the original ones
rownames(transformed_dataframe) <- NULL

rename_data=cbind(as.vector(paste(sample_data(feast_file_its)$RowName,"_",sample_data(feast_file_its)$core_type,sep="")), 
                  as.vector(sample_data(feast_file_its)$species.x))
colnames(rename_data)=c("Row_Name","Species")

transformed_dataframe=merge(transformed_dataframe, rename_data, by="Row_Name")

plot <- ggplot(transformed_dataframe, aes(x = Species, y = Value, fill=Row_Name, group=Column_Name)) +
  #  geom_bar(stat = "identity", position = "dodge_2") +
  facet_grid(Column_Name~., scales="free")+
  geom_jitter(aes(col=Species))+
  geom_hline(yintercept = 0.5)+
  labs(x = "Sample", y = "Percent Contribution", fill = "Contribution") +
  theme(legend.position = "none") +
  #  ylim(0,1) +
  stat_summary(aes(col=Species), size=1.2)
plot

plot + 
  stat_summary(
    fun = function(x) mean(x, na.rm = TRUE), # Change the summary function if needed
    geom = "text",
    aes(label = sprintf("%.2f", after_stat(y)), group = Column_Name), # Use 'after_stat' to refer to summary value
    vjust = -0.5,
    show.legend = FALSE
  )

transformed_dataframe=subset(transformed_dataframe, Species!="NA")

plot <- ggplot(transformed_dataframe, aes(x = Column_Name, y = Value, fill=Species, group=Column_Name)) +
  geom_bar(stat = "identity", position = "dodge2") +
  labs(x = "Sample", y = "Percent Contribution", fill = "Contribution") +
 # theme(legend.position = "none") +
  #  ylim(0,1) +
  facet_grid(.~Column_Name, scale="free_x")+
  stat_summary(size=1, (aes(col=Column_Name)))+
  ggtitle("ITS: Sapwood Origins")+
  scale_fill_manual(values = species_colors_its, name="Species Code")
plot

plot + 
  stat_summary(
    fun = function(x) mean(x, na.rm = TRUE), # Change the summary function if needed
    geom = "text",
    aes(label = sprintf("%.2f", after_stat(y)), group = Column_Name), # Use 'after_stat' to refer to summary value
    vjust = -0.5,
    show.legend = FALSE, 
    size=6
  )



plot <- ggplot(transformed_dataframe, aes(x = Column_Name, y = Value, fill = Species, group = Species)) +
  stat_summary(
    fun = mean,  # Calculate the mean
    geom = "bar",  # Use bar geometry for the mean
    position = position_dodge2(width = 0.9),  # Adjust dodge width if necessary
    na.rm = TRUE  # Remove NA values for mean calculation
  ) +
  labs(x = "Sample", y = "Percent Contribution", fill = "Contribution") +
  ggtitle("ITS: Heartwood Origins") +
  scale_fill_manual(values = species_colors_its, name = "Species Code") +
  facet_grid(. ~ Column_Name, scale = "free_x")

plot

plot + 
  stat_summary(
    fun = function(x) mean(x, na.rm = TRUE), # Change the summary function if needed
    geom = "text",
    aes(label = sprintf("%.2f", after_stat(y)), group = Column_Name), # Use 'after_stat' to refer to summary value
    vjust = -0.5,
    show.legend = FALSE, 
    size=6
  )

# * 16S ----

#### ```` Metabolisms  ####
ps.g=tax_glom(ps.ra, taxrank=rank_names(ps.ra)[5], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
ps.g.filt=prune_taxa(taxa_sums(ps.g)>10, ps.g)
ps.g.filt.wood=subset_samples(ps.g.filt, material=="Wood")
ps.g.filt.wood=subset_samples(ps.g.filt.wood, species.x!="KALA"&
                                species.x!="SAAL" &
                                species.x!="QUAL" &
                                species.x!="QUVE" &
                                species.x!="CAOV" &
                                species.x!="PRSE")



FTX_tt=tax_table(FTX)
row.names(FTX_tt)=row.names(FTX)
tax_table(ps.g.filt.wood)=FTX_tt
colnames(tax_table(ps.g.filt.wood))=colnames(FTX)

# change fill to e.g. "fermentation"  with FAPROTAX
fpt=plot_bar(ps.g.filt.wood, x="seq_id", fill ="plant_pathogen")+
  geom_bar(stat = "identity", )+
  #scale_fill_manual(values = c("#E3A692","#009999"))+
  scale_fill_manual(values=coul)+
  facet_grid(core_type~., scales="free", space = "free") +
  ylim(0,100) + 
  theme(legend.position = "bottom")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
fpt

#methylotrophy
fpt2=plot_bar(ps.g.filt.wood, x="seq_id", fill ="dark_hydrogen_oxidation")+
  geom_bar(stat = "identity", )+
  #scale_fill_manual(values = c("#E3A692","#009999"))+
  scale_fill_manual(values=coul)+
  facet_grid(core_type~species.x, scales="free", space = "free") +
  ylim(0,100) + 
  theme(legend.position = "bottom")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
fpt2

cowplot::plot_grid(fpt, fpt2)

sample_data(ps.g.filt.wood)$mcra10=log10(sample_data(ps.g.filt.wood)$mcra_probe_loose+1)

fpt=plot_bar(ps.g.filt.wood, x="seq_id", fill ="mcra10")+
  geom_bar(stat="identity")+scale_fill_viridis_c()+
  facet_grid(core_type~species.x, scales="free", space = "free") + ylim(0,100) + theme(legend.position = "bottom")



### ```` Alpha ####
ps.rare.clean=subset_samples(ps.rare, core_type=="Inner" | core_type=="Outer" | core_type=="Mineral" | core_type=="Organic")
alpha_16s=plot_richness(ps.rare.clean, x="species.x", measures=c("Chao1", "Shannon"), col="core_type", title = "16S")+
  stat_summary(size=1.2)

shannon_alpha_16s=subset(alpha_16s$data, alpha_16s$data$variable=="Shannon")
sh_al_16s_wood=subset(shannon_alpha_16s, core_type=="Inner")
sh_al_16s_soil=subset(shannon_alpha_16s, core_type=="Mineral")
alpha_16s_merge=merge(sh_al_16s_soil, sh_al_16s_wood, by="seq_id")

alp_cor=ggplot(alpha_16s_merge, aes(x=value.x,
                                     y=value.y, col=species.x.y, group=core_type.x))+
  geom_point()+
  geom_smooth(aes(x=value.x,
                 y=value.y),method = "lm", se=FALSE)+
  xlab("Mineral Soil Alpha Diveersity (Shannon)")+
  ylab("Inner Wood Alpha Diveersity (Shannon)")+
 # facet_wrap(.~species.x.y, scales="free")+
  stat_cor()
alp_cor


meco_bact_alpha <- clone(meco_bacteria)
meco_bact_alpha$sample_table <- subset(meco_bact_alpha$sample_table,  core_type == "Inner" | core_type == "Outer")
meco_bact_alpha$tidy_dataset()

ta_16<- trans_alpha$new(dataset = meco_bact_alpha, group = "core_type")
ta_16$cal_diff(method = "anova", formula = "species.y")
head(ta_16$res_diff)


#### ```` Heatmaps ----
meco_bact_all_wood <- clone(meco_bacteria)
meco_bact_all_wood$sample_table <- subset(meco_bact_all_wood$sample_table,  core_type=="Inner" | core_type=="Outer")
meco_bact_all_wood$tidy_dataset()

bact_species=meco_bact_all_wood$merge_samples(use_group = "core_type")

t16_s = trans_func$new(bact_species)
t16_s$cal_spe_func(prok_database = "FAPROTAX")
t16_s$cal_spe_func_perc(abundance_weighted = TRUE, dec=3)
#t16_s$res_spe_func_perc=t16_s$res_spe_func_perc[c(18,12,7)]

#t16_s$sample_table$core_type=factor(t16_s$sample_table$core_type, levels=c("Inner","Outer","Organic","Mineral"))
row.names(t16_s$res_spe_func_perc)=c("Heartwood", "Sapwood")

data <- t16_s$res_spe_func_perc
filtered_data <- data[, colSums(data >= 0.05) > 0]
t16_s$res_spe_func_perc=filtered_data

t16_s$plot_spe_func_perc(add_facet = TRUE, use_group_list = FALSE) +
  geom_text(aes(label = round(value, 2)), color = "white", size = 5) +  # Add labels (numbers) to the heatmap
  scale_fill_gradientn(colors = c("#440154","#3b528b","#21918c", "#5ec962", "#fde725"), na.value = "grey",
                       trans = "log10",  # Logarithmic transformation
                       limits = c(0.01, NA)) +  # Automatically scale to the data range
  theme(axis.text.y = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))


 # stat_compare_means(col="white")


t16_s$sample_table$species.y[t16_s$sample_table$species.y!="ACSA"]="Other"

func_comp=ggplot(data=fortify(t16_s$sample_table), aes(x=t16_s$sample_table$core_type, y=t16_s$res_spe_func_perc$dark_hydrogen_oxidation))+
  geom_boxplot()+
  stat_compare_means()
func_comp
  

# Environmental heatmaps
meco_bact_all_wood$sample_table$Percent_Change[(meco_bact_all_wood$sample_table$Percent_Change)<0]=NA
meco_bact_all_wood$sample_table$mcra_probe_loose=log10(as.numeric(meco_bact_all_wood$sample_table$mcra_probe_loose)+1)
meco_bact_all_wood$sample_table$pmoa_strict=log10(as.numeric(meco_bact_all_wood$sample_table$pmoa_strict)+1)

t3 <- trans_env$new(dataset = meco_bact_all_wood, env_cols = 
                      c("mcra_probe_loose", "pmoa_strict", "mmox_strict","ch4_50","ch4_125","ch4_200",
                        "CH4_int", "CO2_int","N2O_int", "O2_int", "co2_50","co2_125","co2_200", "mean_vwc",
                        "mean_st","X16S_per_ul"), complete_na = TRUE, character2numeric = TRUE)

t3$cal_cor(add_abund_table = t16_s$res_spe_func_perc, cor_method = "spearman", 
           by_group = "core_type")
t3$plot_cor()

t3$cal_cor(use_data = "Class", p_adjust_method = "none", use_taxa_num = 30,
           by_group = "core_type")
t3$plot_cor(filter_feature = c("","*","**"))

t3$cal_autocor(group="core_type")


#### Taxa Heatmap
t3 <- trans_abund$new(dataset = meco_bacteria, taxrank = "Class", ntaxa = 10)
t3$data_taxanames=c("Clostridia","Bacilli", "Negativicutes","Methanobacteria")

t3$plot_heatmap(xtext_keep = FALSE, withmargin = FALSE) +
  facet_grid(.~meco_bacteria$sample_table$SampleID, scales="free")



meco_bact_all_wood <- clone(meco_bacteria)
meco_bact_all_wood$sample_table <- subset(meco_bact_all_wood$sample_table,  core_type=="Outer")
meco_bact_all_wood$tidy_dataset()

meco_bact_merge=meco_bact_all_wood$merge_samples(use_group = "species.y")
meco_bact_merge$sample_table$SampleID <- factor(meco_bact_merge$sample_table$SampleID, 
                                                levels = c("Inner", "Outer", "Mineral", "Organic"), 
                                                labels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))

# Create a trans_abund object with the merged dataset
t3 <- trans_abund$new(dataset = meco_bact_merge, taxrank = "Class", ntaxa = 23)
t3$data_abund$Sample=factor(t3$data_abund$Sample, 
                            levels = c("Inner", "Outer", "Mineral", "Organic"), 
                            labels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))

# Plot the heatmap with the updated x-axis labels
t3$plot_heatmap(xtext_keep = TRUE, withmargin = TRUE, xtext_size = 16, ytext_size = 16, plot_numbers = TRUE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "right")+
  theme(legend.title = element_text(size=12))+
  theme(legend.title = element_text())+
  theme(legend.text = element_text(size=12))

#### ```` Abundance ----

t4 <- trans_abund$new(dataset = meco_bacteria, taxrank = "Class", ntaxa = 10, group ="core_type")
t4$plot_pie(facet_nrow = 1)+scale_fill_manual(values=coul[12:40])


dataset_merge <- meco_bact_wood$merge_samples(use_group = "core_type")
t5 <- trans_venn$new(dataset_merge, ratio = "numratio")
t5$plot_venn()

meco_bact_wood <- clone(meco_bacteria)
meco_bact_wood$sample_table <- subset(meco_bact_wood$sample_table,  core_type == "Inner" | core_type == "Outer")
meco_bact_wood$tidy_dataset()
meco_bact_group <- meco_bact_wood$merge_samples(use_group = "species.y")
t13 <- trans_venn$new(dataset = meco_bact_group, ratio="seqratio")
t13$plot_venn(petal_plot = TRUE, petal_center_size = 50, petal_r = 1.5, petal_a = 3, petal_move_xy = 3.8, petal_color_center = "#BEBADA")

#### ```` Diff ----
t6 <- trans_diff$new(dataset = meco_bacteria, method = "lefse", group = "core_type", alpha = 0.01, lefse_subgroup = NULL, taxa_level = "Class")
t6$plot_diff_bar(threshold = 3.5)
t6$res_abund$Group=factor(t6$res_abund$Group, 
                          levels = c("Inner", "Outer", "Mineral", "Organic"), 
                          labels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))
t6$res_diff$Group=factor(t6$res_diff$Group, 
                          levels = c("Inner", "Outer", "Mineral", "Organic"), 
                          labels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))
t6$plot_diff_bar(use_number = 1:30, width = 0.8, group_order = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"), keep_prefix = FALSE, axis_text_y = 14)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(size=12))+
  theme(legend.text = element_text(size=12))+
  scale_fill_manual(values=custom_palette)+
  scale_color_manual(values=custom_palette)


t6$plot_diff_cladogram(use_taxa_num = 100, use_feature_num = 50, clade_label_level = 1)

t6 <- trans_diff$new(dataset = meco_bact_wood, method = "lefse", group = "species.y", alpha = 0.01, lefse_subgroup = NULL)
t6$plot_diff_bar(threshold = 4.20)
t6$plot_diff_bar(use_number = 1:30, width = 0.8)



#### ```` Environment ----

t7 <- trans_env$new(dataset = meco_bacteria, env_cols=c(35:47), complete_na = TRUE)
t7 <- trans_env$new(dataset = meco_bacteria, env_cols=c(22,31,33, 63), complete_na = TRUE)

t7$cal_cor(by_group = "core_type", use_data = "other", p_adjust_method = "fdr", other_taxa = t6$res_diff$Taxa[1:40])
t7$cal_cor(by_group = "core_type",p_adjust_method = "fdr", p_adjust_type = "Env",use_data="Or")
t7$plot_cor(filter_feature = c("", "*", "**"))

t7$cal_cor(by_group = "core_type",add_abund_table = t16$res_spe_func_perc, cor_method = "spearman")
t7$plot_cor()







#### ```` Beta ----
t7 <- trans_env$new(dataset = meco_bact_all_wood, env_cols=c( 22, 44, 38, 40, 42), complete_na = TRUE)
t7$cal_ordination(method = "CCA", use_measure = "bray", taxa_level = "Class")
t7$res_ordination_R2
# t1$res_rda is the result list stored in the object
#t7$trans_ordination(adjust_arrow_length = TRUE, max_perc_env = 1.5)
# t1$res_rda_trans is the transformed result for plotting
t7$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
t7$plot_ordination(plot_color = "core_type")

t7$cal_ordination_anova()
t7$res_ordination_terms

t7$cal_ordination_envfit()
t7$res_ordination_envfit

t7$cal_mantel(use_measure = "bray", by_group ="core_type")
t7$res_mantel

t7$ca

### Beta-div within species

meco_bact_species <- clone(meco_bacteria)
meco_bact_species$sample_table <- subset(meco_bact_species$sample_table,  is.na(species.x)==FALSE & core_type=="Outer")
meco_bact_species$tidy_dataset()
meco_bact_species$cal_betadiv(unifrac=TRUE)
t10 <- trans_beta$new(dataset = meco_bact_species, group = "species.x", measure = "unwei_unifrac")
t10$cal_group_distance(within_group = TRUE)
#t10$cal_group_distance_diff(method = "anova")
t10$plot_group_distance(boxplot_add = "mean")+ylim(0.5,1.1)

beta_16s_within=t10$res_group_distance
beta_16s_within$amplicon="16S"
beta_its_within=t10_its$res_group_distance
beta_its_within$amplicon="ITS"

beta_compare=rbind(beta_16s_within, beta_its_within)
beta_compare=subset(beta_compare, species.x!="FAGR")

compbeta=ggplot(data=beta_compare, aes(x=species.x, y= Value, fill=amplicon))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")
compbeta

t10$cal_manova()
t10$res_manova

t10_its$cal_manova()
t10_its$res_manova


#meco_bact_all_wood$cal_betadiv(unifrac=TRUE)
t_core <- trans_beta$new(dataset = meco_bact_all_wood, group = "core_type", measure = "unwei_unifrac")
t_core$cal_group_distance(within_group = TRUE)
#t10$cal_group_distance_diff(method = "anova")
t_core$plot_group_distance(boxplot_add = "mean")+ylim(0.5,1.1)


#### Betadiv phylogenies

library(GUniFrac)
meco_bact_wood <- clone(meco_bacteria)
meco_bact_wood$sample_table <- subset(meco_bact_wood$sample_table, core_type=="Inner" & species.y!="MISC")
meco_bact_wood$tidy_dataset()

# By individual tree
meco_bact_wood$cal_betadiv(unifrac = TRUE)
t2_w_16s <- trans_beta$new(dataset = meco_bact_wood, measure = "unwei_unifrac", group = "species.x")
t2_w_16s$plot_clustering(group = "core_type", replace_name = c("species.x"), color_values = coul[9:30])

# By species
meco_bact_group <- meco_bact_wood$merge_samples(use_group = "species.y")
meco_bact_group$cal_betadiv(unifrac = TRUE)
t3_16s_group <- trans_beta$new(dataset = meco_bact_group, measure = "wei_unifrac")

t3_16s_group$plot_clustering(group = "SampleID", replace_name = c("SampleID"), color_values = 
                               c("red","red","red","black","red",
                                 "black","black","black","black","black",
                                 "black","black","red","red","black","black"))+
  ggtitle("16S Heartwood, Weighted Unifrac")

#### Radar plots

t3_16s <- trans_abund$new(dataset = meco_bact_group, taxrank = "Class", ntaxa = 8, groupmean = "SampleID")

t3_16s$plot_radar(values.radar = c("0%", "25%", "50%"), grid.min = 0, grid.mid = 0.25, grid.max = 1)+
  ggtitle("Inner Wood, 16S")

t3_16s$plot_donut(label = FALSE, facet_nrow = 4, color_values =  )+
  ggtitle("Organic Soil, 16S")



#### Tests

meco_bact_adonis <- clone(meco_bacteria)
meco_bact_adonis$sample_table <- subset(meco_bact_adonis$sample_table)
meco_bact_adonis$tidy_dataset()
meco_bact_adonis$cal_betadiv(unifrac=TRUE)

t_b <- trans_beta$new(dataset = meco_bact_adonis, measure = "wei_unifrac")
t_b$sample_table$material[is.na(t_b$sample_table$material)]="MISC"
t_b$cal_manova(group = "material")
t_b$res_manova



# Ordination

t_b <- trans_beta$new(dataset = meco_bact_adonis, measure = "wei_unifrac")
t_b$cal_ordination(ordination = "PCoA")
t_b$res_ordination$scores$core_type=factor(t_b$res_ordination$scores$core_type, 
                                            levels = c("Inner", "Outer", "Mineral", "Organic"), 
                                            labels = c("Heartwood", "Sapwood", "Mineral Soil", "Organic Soil"))


# Only below this point edit
custom_palette <- rev(c("#018571", "#80cdc1", "#dfc27d", "#a6611a"))

t_b$plot_ordination(plot_color = "core_type", plot_type = c("point","centroid"), centroid_segment_alpha = 0.25)+
  ggtitle("PCoA, WUniFrac: 16S")+
  geom_point(size=4, alpha=0.75)+
  scale_color_manual(values = custom_palette, name="Sample Type")+
  theme_minimal(base_size = 16)




#### ```` FEAST ####
ps.g.count=tax_glom(ps.rare, taxrank=rank_names(ps.rare)[5], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
feast_file = subset_samples(ps.g.count, core_type=="Inner" | core_type=="Outer" | core_type=="Mineral" | core_type=="Organic")
 #Should glom to speed up and reduce unknowns

otus_16s=as.matrix(as.data.frame(otu_table(feast_file)))
metadata_16s=as.data.frame(matrix(NA,ncol=4))


### Unpaired Sink/Source Route ###
SampleID=colnames(otus_16s)
Env=sample_data(feast_file)$core_type
SourceSink=ifelse(sample_data(feast_file)$core_type=="Inner","Sink","Source") # Specify here 

count_occurrences <- function(lst) {
  counts <- c()
  seen <- list()
  
  for (value in lst) {
    if (is.null(seen[[value]])) {
      seen[[value]] <- 1
      counts <- c(counts, 1)
    } else {
      seen[[value]] <- seen[[value]] + 1
      counts <- c(counts, seen[[value]])
    }
  }
  
  return(counts)
}
id=count_occurrences(Env)
id[SourceSink=="Source"]=NA

feast_input=as.data.frame(cbind(Env, SourceSink, id))
feast_input$id=as.numeric(feast_input$id)
rownames(feast_input)=SampleID


### Paired Sink/Source Route ###
SampleID=colnames(otus_16s)
Env=sample_data(feast_file)$core_type
SourceSink=ifelse(sample_data(feast_file)$core_type=="Outer","Sink","Source") # Specify here 
Pair_key=sample_data(feast_file)$seq_id
feast_input=as.data.frame(cbind(Env, SourceSink, id, Pair_key))

"
feast_input <- feast_input %>%
  arrange(Pair_key)

feast_input <- feast_input %>%
  group_by(Pair_key) %>%
  mutate(Group_number = cur_group_id())
"

feast_input <- feast_input %>%
  arrange(Pair_key) %>%
  group_by(Pair_key, SourceSink) %>%
  mutate(Group_number = cur_group_id()) %>%
  ungroup() %>%
  group_by(Pair_key) %>%
  mutate(Has_Source = "Source" %in% SourceSink,
         Has_Sink = "Sink" %in% SourceSink,
         Group_number = case_when(Has_Source & Has_Sink ~ Group_number[SourceSink == "Source"][1],
                                  TRUE ~ NA)) %>%
  select(-Has_Source, -Has_Sink)

feast_input$id=as.numeric(feast_input$Group_number)
feast_input=as.data.frame(feast_input)
feast_input=feast_input[,c(1:3)]

#SampleID originally, now species
rownames(feast_input)=SampleID

feast_input=subset(feast_input, is.na(id)==FALSE)

# Convert lists to vectors
feast_input$Env <- sapply(feast_input$Env, function(x) x[[1]])
feast_input$SourceSink <- sapply(feast_input$SourceSink, function(x) x[[1]])

### Run FEAST

write.table(as.data.frame(feast_input), "/Users/Wyatt/Desktop/FEAST/feast_input_paired.txt", sep = '\t')
innyy=Load_metadata("/Users/Wyatt/Desktop/FEAST/feast_input_paired.txt")

write.table(otus_16s, "/Users/Wyatt/Desktop/FEAST/feast_otu_paired.txt", sep = '\t')
otuu=Load_CountMatrix("/Users/Wyatt/Desktop/FEAST/feast_otu_paired.txt")

#innyy=Load_metadata("/Users/Wyatt/Desktop/metadata_example_multi.txt")
#otuu=Load_CountMatrix("/Users/Wyatt/Desktop/otu_example_multi.txt")

# Different sources flag changed for paired (1) /unpaired (0)
FEAST_output <- FEAST(C = otuu, metadata = innyy, different_sources_flag=1, dir_path = "/Users/Wyatt/Desktop/FEAST/",
                      outfile="16s_Inner_paired")

FEAST_tab=read.delim("/Users/Wyatt/Desktop/FEAST/16s_Inner_paired_source_contributions_matrix.txt")
colnames(FEAST_tab)=subset(feast_input, Env!="Outer")$Env
colnames(FEAST_tab)[is.na(colnames(FEAST_tab))]="Unknown"

### Plot FEAST

create_barplot <- function(data) {
  # Transpose the data frame
  transposed_data <- t(data)
  
  # Create a new column with row names as Sample
  transposed_data <- data.frame(Sample = row.names(transposed_data), transposed_data, row.names = NULL)
  
  # Convert data frame to long format
  melted_data <- reshape2::melt(transposed_data, id.vars = "Sample", variable.name = "Contribution", value.name = "Percent")
  
  # Create the barplot using ggplot2
  plot <- ggplot(melted_data, aes(x = Sample, y = Percent, fill = Contribution)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Sample", y = "Percent Contribution", fill = "Contribution") +
    theme_minimal()+ theme(legend.position = "none") +
    ylim(0,1)
    
  
  return(plot)
}
create_barplot((FEAST_tab))



min_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="Mineral")], na.rm = TRUE)
org_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="Organic")], na.rm = TRUE)
unk_sums=(FEAST_tab[,which(colnames(FEAST_tab)=="Unknown")])
out_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="Inner")], na.rm = TRUE)

feast_sum=cbind(min_sums, org_sums, unk_sums, out_sums)
colnames(feast_sum)=c("Mineral","Organic","Unknown","Inner")
feast_sum=as.data.frame(feast_sum)

library(tidyr)

# Add the row names as a separate column
feast_sum$Row_Name <- rownames(feast_sum)

# Transform the dataframe to have two columns: 'Row_Name', 'Value', and 'Column_Name'
transformed_dataframe <- pivot_longer(feast_sum, 
                                      cols = c("Mineral", "Inner", "Organic", "Unknown"), 
                                      names_to = "Column_Name", 
                                      values_to = "Value")

# Set the row names back to the original ones
rownames(transformed_dataframe) <- NULL

rename_data=cbind(as.vector(paste(sample_data(feast_file)$RowName,"_",sample_data(feast_file)$core_type,sep="")), 
                            as.vector(sample_data(feast_file)$species.x))
colnames(rename_data)=c("Row_Name","Species")

transformed_dataframe=merge(transformed_dataframe, rename_data, by="Row_Name")

plot <- ggplot(transformed_dataframe, aes(x = Species, y = Value, fill=Row_Name, group=Column_Name)) +
#  geom_bar(stat = "identity", position = "dodge_2") +
  facet_grid(Column_Name~., scales="free")+
  geom_jitter(aes(col=Species))+
  labs(x = "Sample", y = "Percent Contribution", fill = "Contribution") +
   theme(legend.position = "none") +
#  ylim(0,1) +
  stat_summary(aes(col=Species), size=1.2)+
  geom_hline(yintercept = 0.5)
plot

plot + 
  stat_summary(
    fun = function(x) mean(x, na.rm = TRUE), # Change the summary function if needed
    geom = "text",
    aes(label = sprintf("%.2f", after_stat(y)), group = Column_Name), # Use 'after_stat' to refer to summary value
    vjust = -0.5,
    show.legend = FALSE
  )




plot <- ggplot(transformed_dataframe, aes(x = Column_Name, y = Value, fill=Row_Name, group=Column_Name)) +
  geom_bar(stat = "identity", position = "dodge_2") +
  labs(x = "Sample", y = "Percent Contribution", fill = "Contribution") +
  theme(legend.position = "none") +
  #  ylim(0,1) +
  stat_summary(aes(col="red"), size=1.2)+
  ggtitle("16S: Inner")
plot

plot + 
  stat_summary(
    fun = function(x) mean(x, na.rm = TRUE), # Change the summary function if needed
    geom = "text",
    aes(label = sprintf("%.2f", after_stat(y)), group = Column_Name), # Use 'after_stat' to refer to summary value
    vjust = -0.5,
    show.legend = FALSE
  )



transformed_dataframe=subset(transformed_dataframe, Species!="NA")

plot <- ggplot(transformed_dataframe, aes(x = Column_Name, y = Value, fill = Species, group = Species)) +
  stat_summary(
    fun = mean,  # Calculate the mean
    geom = "bar",  # Use bar geometry for the mean
    position = position_dodge2(width = 0.9),  # Adjust dodge width if necessary
    na.rm = TRUE  # Remove NA values for mean calculation
  ) +
  labs(x = "Sample", y = "Percent Contribution", fill = "Contribution") +
  ggtitle("16S: Sapwood Origins") +
  scale_fill_manual(values = species_colors, name = "Species Code") +
  facet_grid(. ~ Column_Name, scale = "free_x")

plot

plot + 
  stat_summary(
    fun = function(x) mean(x, na.rm = TRUE), # Change the summary function if needed
    geom = "text",
    aes(label = sprintf("%.2f", after_stat(y)), group = Column_Name), # Use 'after_stat' to refer to summary value
    vjust = -0.5,
    show.legend = FALSE, 
    size=6
  )


#### * QUVE 16S ####

#### ```` Abundance #####
t1_q <- trans_abund$new(dataset = meco_quve_16s, taxrank = "Class", ntaxa = 9, groupmean = "core_type")
t1_q$plot_radar(values.radar = c("0%", "25%", "50%"), grid.min = 0, grid.mid = 0.25, grid.max = 0.75)

t1_q$plot_donut(label = FALSE,facet_nrow = 4)


meco_quve_16s_norm = clone(meco_quve_16s)
meco_quve_16s_norm$sample_table <- subset(meco_quve_16s_norm$sample_table,  
                                          core_type != "inner" | TreatmentGroup=="QUVE200HEART2"
                                          | TreatmentGroup=="QUVE50HEART2"
                                          | TreatmentGroup=="QUVE400HEART2"
                                          | TreatmentGroup=="QUVE600HEART2")
meco_quve_16s_norm$sample_table <- subset(meco_quve_16s_norm$sample_table,  
                                          core_type != "outer" | TreatmentGroup=="QUVE200SAP"
                                          | TreatmentGroup=="QUVE50SAP"
                                          | TreatmentGroup=="QUVE400SAP"
                                          | TreatmentGroup=="QUVE600SAP")

meco_quve_16s_norm$tidy_dataset()

meco_quve_16s$sample_table$core_type_2=meco_quve_16s$sample_table$core_type
meco_quve_16s$sample_table$core_type_2[meco_quve_16s$sample_table$core_type_2=="inner" |
                                         meco_quve_16s$sample_table$core_type_2=="outer"]="wood"     
  
dataset1_q <- meco_quve_16s_norm$merge_samples(use_group = "core_type")
t1 <- trans_venn$new(dataset1_q, ratio="seqratio")
t1$plot_venn(petal_plot = TRUE, petal_center_size = 50, petal_r = 1.5, petal_a = 3, petal_move_xy = 3.8, petal_color_center = "#BEBADA")

#### ```` Alpha ####
alpha_quve16=plot_richness(ps.quve.16s.rare, x="core_type", measures=c("Chao1", "Shannon"), col="core_type", title="16S")+
  stat_summary(size=1.2)
alpha_quve16

#### ```` FEAST #####
library(FEAST)
ps.quve.16s=tax_glom(ps.quve.16s.rare, taxrank=rank_names(ps.quve.16s.rare)[6], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

feast_file=subset_samples(ps.quve.16s, is.na(sample_data(ps.quve.16s)$core_type)==FALSE) #Should glom to speed up and reduce unknowns

otus_q16=as.matrix(as.data.frame(otu_table(feast_file)))
metadata_q16=as.data.frame(matrix(NA,ncol=4))

SampleID=colnames(otus_q16)
Env=sample_data(feast_file)$core_type
SourceSink=ifelse(sample_data(feast_file)$core_type=="inner","Sink","Source")

count_occurrences <- function(lst) {
  counts <- c()
  seen <- list()
  
  for (value in lst) {
    if (is.null(seen[[value]])) {
      seen[[value]] <- 1
      counts <- c(counts, 1)
    } else {
      seen[[value]] <- seen[[value]] + 1
      counts <- c(counts, seen[[value]])
    }
  }
  
  return(counts)
}
id=count_occurrences(Env)
id[SourceSink=="Source"]=NA


feast_input=as.data.frame(cbind(Env, SourceSink, id))
feast_input$id=as.numeric(feast_input$id)
rownames(feast_input)=SampleID


write.table(feast_input, "/Users/Wyatt/Desktop/FEAST/quve_16s_feast_input.txt", sep = '\t')
innyy=Load_metadata("/Users/Wyatt/Desktop/FEAST/quve_16s_feast_input.txt")

write.table(otus_q16, "/Users/Wyatt/Desktop/FEAST/quve_16s_feast_otu.txt", sep = '\t')
otuu=Load_CountMatrix("/Users/Wyatt/Desktop/FEAST/quve_16s_feast_otu.txt")

#innyy=Load_metadata("/Users/Wyatt/Desktop/metadata_example_multi.txt")
#otuu=Load_CountMatrix("/Users/Wyatt/Desktop/otu_example_multi.txt")

FEAST_output <- FEAST(C = otuu, metadata = innyy, different_sources_flag=0, dir_path = "/Users/Wyatt/Desktop/FEAST",
                      outfile="quve_16s_outer")

FEAST_tab=read.delim("/Users/Wyatt/Desktop/FEAST/quve_16s_outer_source_contributions_matrix.txt")
colnames(FEAST_tab)=subset(feast_input, Env!="inner")$Env
colnames(FEAST_tab)[is.na(colnames(FEAST_tab))]="Unknown"


create_barplot <- function(data) {
  # Transpose the data frame
  transposed_data <- t(data)
  
  # Create a new column with row names as Sample
  transposed_data <- data.frame(Sample = row.names(transposed_data), transposed_data, row.names = NULL)
  
  # Convert data frame to long format
  melted_data <- reshape2::melt(transposed_data, id.vars = "Sample", variable.name = "Contribution", value.name = "Percent")
  
  # Create the barplot using ggplot2
  plot <- ggplot(melted_data, aes(x = Sample, y = Percent, fill = Contribution)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Sample", y = "Percent Contribution", fill = "Contribution") +
    theme_minimal()+
    ylim(0,1)+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 90))+
    stat_summary(aes(group=Sample))
  
  return(plot)
}
create_barplot((FEAST_tab))

creat_bp_2<- function(FEAST_tab){
  
  min_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="mineral")], na.rm = TRUE)
  org_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="organic")], na.rm = TRUE)
  unk_sums=(FEAST_tab[,which(colnames(FEAST_tab)=="Unknown")])
  out_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="outer")], na.rm = TRUE)
  branch_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="branch")], na.rm = TRUE)
  bark_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="bark")], na.rm = TRUE)
  coarse_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="coarse")], na.rm = TRUE)
  fine_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="fine")], na.rm = TRUE)
  rot_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="rot")], na.rm = TRUE)
  leaf_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="leaf")], na.rm = TRUE)
  litter_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="litter")], na.rm = TRUE)
  
  
  feast_sum=cbind(min_sums, org_sums, unk_sums, out_sums, branch_sums, bark_sums, coarse_sums, fine_sums,
                  rot_sums, leaf_sums, litter_sums)
  colnames(feast_sum)=c("Mineral","Organic","Unknown","Outer","Branch","Bark","Coarse","Fine","Rot","Leaf",'Litter')
  feast_sum=as.data.frame(feast_sum)
  
  library(tidyr)
  
  # Add the row names as a separate column
  feast_sum$Row_Name <- rownames(feast_sum)
  
  # Transform the dataframe to have two columns: 'Row_Name', 'Value', and 'Column_Name'
  transformed_dataframe <- pivot_longer(feast_sum, 
                                        cols = c("Mineral","Organic","Unknown","Outer","Branch","Bark","Coarse","Fine","Rot","Leaf",'Litter'), 
                                        names_to = "Column_Name", 
                                        values_to = "Value")
  return(transformed_dataframe)
}
q16_inner_result=creat_bp_2(FEAST_tab)

q16_inner_result$Column_Name[q16_inner_result$Column_Name=="Outer"]="Sapwood"
#q16_inner_result$Column_Name[q16_inner_result$Column_Name=="Inner"]="Heartwood"

plot_in<- ggplot(q16_inner_result, aes(x = Column_Name, y = Value, group=Column_Name), col="black") +
  labs(x = "", y = "Contribution") +
  theme(legend.position = "none") +
  geom_jitter(width = .15, size=3) +
  #  ylim(0,1) +
  stat_summary(aes(col=tolower(Column_Name)), size=1.5, alpha=0.85)+
  stat_summary(aes(color = tolower(Column_Name)), fun.data = mean_se, geom = "linerange", size = 1.5, alpha=0.7) + # Thicker SE bars
  theme(axis.text.x = element_text(angle = 45))+
  theme_pubclean(base_size=18)+
  ggtitle("Heartwood: 16S")+
  # scale_color_manual(values="red", na.value = "black")+
  scale_color_viridis_d(option = "viridis")+
  # geom_vline(aes(xintercept = as.numeric(factor(Column_Name))), linetype = "dotted", color = "gray40")+
  theme(legend.position = "none") +
  coord_flip()+
  scale_y_continuous(breaks = seq(0, max(qits_inner_result$Value), by = 0.25)) +
  theme(
    panel.grid.major.y = element_line(color = "grey80", linetype = "solid"),  # Horizontal gridlines
    panel.grid.major.x = element_line(color = "lightgrey", linetype = "dashed", size = 0.1),  # Vertical gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    legend.position = "none"
  )
plot_in

cowplot::plot_grid(plot_in,plot_out)


# ```` Beta ####
meco_quve_16s$cal_betadiv(unifrac = TRUE)
t2_q <- trans_beta$new(dataset = meco_quve_16s, measure = "unwei_unifrac", group = "core_type")
t2_q$plot_clustering(group = "core_type", replace_name = c("core_type"), color_values = coul[9:30])

t2_q$cal_ordination(ordination = "PCoA")

t2_q$res_ordination$scores$core_type=factor(t2_q$res_ordination$scores$core_type, levels = c("mineral","organic","outer","inner","branch","bark","coarse","fine","rot","leaf",'litter'),
                                            labels=c("Mineral","Organic","Sapwood","Heartwood","Branch","Bark","Coarse","Fine","Rot","Leaf",'Litter'))

t2_q$plot_ordination(plot_color = "core_type", plot_type = c("point","centroid"))+
  theme_minimal(base_size=16)+
  theme(legend.title = element_blank())
  


quve_grouped_meco <- meco_quve_16s$merge_samples(use_group = "core_type")
quve_grouped_meco$cal_betadiv(unifrac=TRUE)
t3_q <- trans_beta$new(dataset = quve_grouped_meco, measure = "wei_unifrac")

t3_q$plot_clustering(group = "SampleID", replace_name = c("SampleID"), color_values = coul[9:30])+
  ggtitle("QUVE Tissue Regions")


#### ```` Height variance ####
meco_quve_16s_wood <- clone(meco_quve_16s)
meco_quve_16s_wood$sample_table <- subset(meco_quve_16s_wood$sample_table,  core_type == "inner")
meco_quve_16s_wood$tidy_dataset()
meco_quve_16s_wood$cal_betadiv(unifrac = TRUE)
t2_q_w <- trans_beta$new(dataset = meco_quve_16s_wood, measure = "wei_unifrac")

t2_q_w$cal_manova(group = "Inner.Core.Sample.ID")
t2_q_w$res_manova


quve_grouped_meco_wood <- meco_quve_16s_wood$merge_samples(use_group = "Inner.Core.Sample.ID")
quve_grouped_meco_wood$cal_betadiv(unifrac = TRUE)

t3_q_w <- trans_beta$new(dataset = quve_grouped_meco_wood, measure = "wei_unifrac")
t3_q_w$plot_clustering(group = "SampleID", replace_name = c("SampleID"), color_values = coul[9:30])+
  ggtitle("QUVE Height, Outer")


t3_q_w <- trans_abund$new(dataset=quve_grouped_meco_wood, taxrank = "Class", ntaxa = 6)
t3_q_w$plot_donut(label = FALSE)

meco_quve_16s_wood$cal_betadiv(unifrac=TRUE)
t3_q_wood <- trans_beta$new(dataset = meco_quve_16s_wood, measure = "wei_unifrac")
t3_q_wood$cal_ordination(ordination = "PCoA")
t3_q_wood$plot_ordination(plot_color = "Inner.Core.Sample.ID", plot_type = c("point"))+
  ggtitle("QUVE Outer, 16S")


#### ```` Diff ----
t6_q <- trans_diff$new(dataset = meco_quve_16s, method = "lefse", group = "core_type", alpha = 0.01, lefse_subgroup = NULL)

t6_q$plot_diff_bar(use_number = 1:20, width = 0.8, color_values = coul, simplify_names = TRUE)
t6_q$plot_diff_abund(use_number = 1:5, color_values = coul)
t6_q$plot_diff_cladogram(use_taxa_num = 20, use_feature_num = 50, clade_label_level = 5, color =  coul)


t6_q$res_diff


#### ```` Functional ----

meco_bact_quve_func <- clone(meco_quve_16s)
meco_bact_quve_func$sample_table <- subset(meco_bact_quve_func$sample_table)
meco_bact_quve_func$tidy_dataset()

q_bact_species=meco_bact_quve_func$merge_samples(use_group = "core_type")

t16q_s = trans_func$new(q_bact_species)
t16q_s$cal_spe_func(prok_database = "FAPROTAX")
t16q_s$cal_spe_func_perc(abundance_weighted = TRUE, dec=3)
#t16_s$res_spe_func_perc=t16_s$res_spe_func_perc[c(18,12,7)]

#t16_s$sample_table$core_type=factor(t16_s$sample_table$core_type, levels=c("Inner","Outer","Organic","Mineral"))
#row.names(t16q_s$res_spe_func_perc)=c("Heartwood", "Sapwood")

data <- t16q_s$res_spe_func_perc
filtered_data <- data[, colSums(data >= 0.05) > 0]
t16q_s$res_spe_func_perc=filtered_data

t16q_s$plot_spe_func_perc(add_facet = TRUE, use_group_list = FALSE) +
  geom_text(aes(label = round(value, 2)), color = "white", size = 5) +  # Add labels (numbers) to the heatmap
  scale_fill_gradientn(colors = c("#440154","#3b528b","#21918c", "#5ec962", "#fde725"), na.value = "grey",
                       trans = "log10",  # Logarithmic transformation
                       limits = c(0.05, NA)) +  # Automatically scale to the data range
  theme(axis.text.y = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))




#### * QUVE ITS ####

t1_q_its <- trans_abund$new(dataset = meco_quve_its, taxrank = "Class", ntaxa = 9, groupmean = "core_type")

t1_q_its$plot_radar(values.radar = c("0%", "25%", "50%"), grid.min = 0, grid.mid = 0.25, grid.max = 1)
t1_q_its$plot_donut(label = FALSE)

meco_quve_its$sample_table$core_type_2=meco_quve_its$sample_table$core_type
meco_quve_its$sample_table$core_type_2[meco_quve_its$sample_table$core_type_2=="inner" |
                                         meco_quve_its$sample_table$core_type_2=="outer"]="wood"                                          
  
dataset1_q <- meco_quve_its$merge_samples(use_group = "core_type")
t1 <- trans_venn$new(dataset1_q, ratio =  "seqratio")
t1$plot_venn(petal_plot = TRUE, petal_center_size = 50, petal_r = 1.5, petal_a = 3, petal_move_xy = 3.8, petal_color_center = "#BEBADA")

#### ```` Alpha ####
alpha_quve_its=plot_richness(ps.quve.its.rare, x="core_type", measures=c("Chao1", "Shannon"), col="core_type", title = "ITS")+
  stat_summary(size=1.2)
alpha_quve_its

cowplot::plot_grid(alpha_quve16, alpha_quve_its)
#### ```` FEAST #####
library(FEAST)

ps.quve.its=tax_glom(ps.quve.its.rare, taxrank=rank_names(ps.quve.its.rare)[6], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

feast_file_its=subset_samples(ps.quve.its, is.na(sample_data(ps.quve.its)$core_type)==FALSE) #Should glom to speed up and reduce unknowns

otus_qITS=as.matrix(as.data.frame(otu_table(feast_file_its)))
metadata_qITS=as.data.frame(matrix(NA,ncol=4))

SampleID=colnames(otus_qITS)
Env=sample_data(feast_file_its)$core_type
SourceSink=ifelse(sample_data(feast_file_its)$core_type=="outer","Sink","Source")

count_occurrences <- function(lst) {
  counts <- c()
  seen <- list()
  
  for (value in lst) {
    if (is.null(seen[[value]])) {
      seen[[value]] <- 1
      counts <- c(counts, 1)
    } else {
      seen[[value]] <- seen[[value]] + 1
      counts <- c(counts, seen[[value]])
    }
  }
  
  return(counts)
}
id=count_occurrences(Env)
id[SourceSink=="Source"]=NA


feast_input_its=as.data.frame(cbind(Env, SourceSink, id))
feast_input_its$id=as.numeric(feast_input_its$id)
rownames(feast_input_its)=SampleID


write.table(feast_input_its, "/Users/Wyatt/Desktop/FEAST/quve_its_feast_input_its.txt", sep = '\t')
innyy_its=Load_metadata("/Users/Wyatt/Desktop/FEAST/quve_its_feast_input_its.txt")

write.table(otus_qITS, "/Users/Wyatt/Desktop/FEAST/quve_its_feast_otu_its.txt", sep = '\t')
otuu_its=Load_CountMatrix("/Users/Wyatt/Desktop/FEAST/quve_its_feast_otu_its.txt")

#innyy=Load_metadata("/Users/Wyatt/Desktop/metadata_example_multi.txt")
#otuu=Load_CountMatrix("/Users/Wyatt/Desktop/otu_example_multi.txt")

FEAST_output_its <- FEAST(C = otuu_its, metadata = innyy_its, different_sources_flag=0, dir_path = "/Users/Wyatt/Desktop/FEAST/",
                      outfile="quve_its_inner")

FEAST_tab_its=read.delim("/Users/Wyatt/Desktop/FEAST/quve_its_inner_source_contributions_matrix.txt")
colnames(FEAST_tab_its)=subset(feast_input_its, Env!="outer")$Env
colnames(FEAST_tab_its)[is.na(colnames(FEAST_tab_its))]="Unknown"


create_barplot <- function(data) {
  # Transpose the data frame
  transposed_data <- t(data)
  
  # Create a new column with row names as Sample
  transposed_data <- data.frame(Sample = row.names(transposed_data), transposed_data, row.names = NULL)
  
  # Convert data frame to long format
  melted_data <- reshape2::melt(transposed_data, id.vars = "Sample", variable.name = "Contribution", value.name = "Percent")
  
  # Create the barplot using ggplot2
  plot <- ggplot(melted_data, aes(x = Sample, y = Percent, fill = Contribution)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Sample", y = "Percent Contribution", fill = "Contribution") +
    theme_minimal()
  
  return(plot)
}
create_barplot((FEAST_tab_its))


creat_bp_2<- function(FEAST_tab){
  
  min_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="mineral")], na.rm = TRUE)
  org_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="organic")], na.rm = TRUE)
  unk_sums=(FEAST_tab[,which(colnames(FEAST_tab)=="Unknown")])
  out_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="inner")], na.rm = TRUE)
  branch_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="branch")], na.rm = TRUE)
  bark_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="bark")], na.rm = TRUE)
  coarse_sums=(FEAST_tab[,which(colnames(FEAST_tab_its)=="coarse")])
  fine_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="fine")], na.rm = TRUE)
  rot_sums=(FEAST_tab[,which(colnames(FEAST_tab)=="rot")])
  leaf_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="leaf")], na.rm = TRUE)
  litter_sums=rowSums(FEAST_tab[,which(colnames(FEAST_tab)=="litter")], na.rm = TRUE)
  
  
  feast_sum=cbind(min_sums, org_sums, unk_sums, out_sums, branch_sums, bark_sums, coarse_sums, fine_sums,
                  rot_sums, leaf_sums, litter_sums)
  colnames(feast_sum)=c("Mineral","Organic","Unknown","Inner","Branch","Bark","Coarse","Fine","Rot","Leaf",'Litter')
  feast_sum=as.data.frame(feast_sum)
  
  library(tidyr)
  
  # Add the row names as a separate column
  feast_sum$Row_Name <- rownames(feast_sum)
  
  # Transform the dataframe to have two columns: 'Row_Name', 'Value', and 'Column_Name'
  transformed_dataframe <- pivot_longer(feast_sum, 
                                        cols = c("Mineral","Organic","Unknown","Inner","Branch","Bark","Coarse","Fine","Rot","Leaf",'Litter'), 
                                        names_to = "Column_Name", 
                                        values_to = "Value")
  return(transformed_dataframe)
} # have to manually change inner vs. outer here
qits_inner_result=creat_bp_2(FEAST_tab_its)

#qits_inner_result$Column_Name[qits_inner_result$Column_Name=="Outer"]="Sapwood"
qits_inner_result$Column_Name[qits_inner_result$Column_Name=="Inner"]="Heartwood"

plot_its_out<- ggplot(qits_inner_result, aes(x = Column_Name, y = Value, group=Column_Name), col="black") +
  labs(x = "", y = "Contribution") +
  theme(legend.position = "none") +
  geom_jitter(width = .15, size=3) +
  #  ylim(0,1) +
  stat_summary(aes(col=tolower(Column_Name)), size=1.5, alpha=0.85)+
  stat_summary(aes(color = tolower(Column_Name)), fun.data = mean_se, geom = "linerange", size = 1.5, alpha=0.7) + # Thicker SE bars
  theme(axis.text.x = element_text(angle = 45))+
  theme_pubclean(base_size=18)+
  ggtitle("Sapwood: ITS")+
 # scale_color_manual(values="red", na.value = "black")+
  scale_color_viridis_d(option = "viridis")+
 # geom_vline(aes(xintercept = as.numeric(factor(Column_Name))), linetype = "dotted", color = "gray40")+
  theme(legend.position = "none") +
  coord_flip()+
  scale_y_continuous(breaks = seq(0, max(qits_inner_result$Value), by = 0.25)) +
  theme(
    panel.grid.major.y = element_line(color = "grey80", linetype = "solid"),  # Horizontal gridlines
    panel.grid.major.x = element_line(color = "lightgrey", linetype = "dashed", size = 0.1),  # Vertical gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    legend.position = "none"
  )
plot_its_out

cowplot::plot_grid(plot_its_in,plot_its_out)


create_color_scheme <- function(words, colors) {
  unique_words <- unique(tolower(words))
  
  if (length(unique_words) > length(colors)) {
    stop("Not enough colors for the unique words.")
  }
  
  color_scheme <- setNames(colors[1:length(unique_words)], unique_words)
  return(color_scheme)
}

# Define the words and colors
words <- c("Mineral", "Organic", "Unknown", "Heartwood", "Branch", 
           "Bark", "Coarse", "Fine", "Rot", "Leaf", "Litter", "Sapwood")
colors <- coul[9:30]

# Generate the color scheme
color_scheme <- create_color_scheme(words, colors)




# ```` Beta ####
meco_quve_its$cal_betadiv(unifrac=TRUE)
t2_q_its <- trans_beta$new(dataset = meco_quve_its, measure = "wei_unifrac", group = "core_type")
t2_q_its$plot_clustering(group = "core_type", replace_name = c("core_type"), color_values = coul[9:30])

t2_q_its$cal_ordination(ordination = "PCoA")

t2_q_its$res_ordination$scores$core_type=factor(t2_q_its$res_ordination$scores$core_type, levels = c("mineral","organic","outer","inner","branch","bark","coarse","fine","rot","leaf",'litter'),
                                            labels=c("Mineral","Organic","Sapwood","Heartwood","Branch","Bark","Coarse","Fine","Rot","Leaf",'Litter'))

t2_q_its$plot_ordination(plot_color = "core_type", plot_type = c("point","centroid"))+
  theme_minimal(base_size=16)+
  theme(legend.title = element_blank())


quve_grouped_meco_its <- meco_quve_its$merge_samples(use_group = "core_type")
quve_grouped_meco_its$cal_betadiv(unifrac=TRUE)
t3_q_its <- trans_beta$new(dataset = quve_grouped_meco_its, measure = "wei_unifrac")

t3_q_its$plot_clustering(group = "SampleID", replace_name = c("SampleID"), color_values = coul[9:30])+
  ggtitle("QUVE Tissue Regions:ITS")


#### ```` Height variance ####
meco_quve_its_wood <- clone(meco_quve_its)
meco_quve_its_wood$sample_table <- subset(meco_quve_its_wood$sample_table,  core_type == "inner")
meco_quve_its_wood$tidy_dataset()
meco_quve_its_wood$cal_betadiv(unifrac = TRUE)
t2_q_w <- trans_beta$new(dataset = meco_quve_its_wood, measure = "wei_unifrac")
t2_q_w$cal_manova(group = "Inner.Core.Sample.ID")
t2_q_w$res_manova

quve_grouped_meco_wood <- meco_quve_its_wood$merge_samples(use_group = "Inner.Core.Sample.ID")
quve_grouped_meco_wood$cal_betadiv(unifrac=TRUE)
t3_q_w <- trans_beta$new(dataset = quve_grouped_meco_wood, measure = "wei_unifrac")

t3_q_w$plot_clustering(group = "SampleID", replace_name = c("SampleID"), color_values = coul[9:30])+
  ggtitle("QUVE Height, Inner, ITS")

t3_q_w <- trans_abund$new(dataset=quve_grouped_meco_wood, taxrank = "Class", ntaxa = 6)
t3_q_w$plot_donut(label = FALSE)

meco_quve_its_wood$cal_betadiv(unifrac=TRUE)
t3_q_wood <- trans_beta$new(dataset = meco_quve_its_wood, measure = "wei_unifrac")
t3_q_wood$cal_ordination(ordination = "PCoA")
t3_q_wood$plot_ordination(plot_color = "Inner.Core.Sample.ID", plot_type = c("ellipse","point"))+
  ggtitle("Inner, ITS")

#### ```` Diff ----
t6_q_its <- trans_diff$new(dataset = meco_quve_its, method = "lefse", group = "core_type", alpha = 0.01, lefse_subgroup = NULL)

t6_q_its$plot_diff_bar(use_number = 1:20, width = 0.8, color_values = coul, simplify_names = TRUE)
t6_q_its$plot_diff_abund(use_number = 1:5, color_values = coul)
t6_q_its$plot_diff_cladogram(use_taxa_num = 20, use_feature_num = 50, clade_label_level = 5, color =  coul)
#### ```` Functional ----

meco_fung_quve_func <- clone(meco_quve_its)
meco_fung_quve_func$sample_table <- subset(meco_fung_quve_func$sample_table)
meco_fung_quve_func$tidy_dataset()

q_fung_species=meco_fung_quve_func$merge_samples(use_group = "core_type")

tITSq_s = trans_func$new(q_fung_species)
tITSq_s$cal_spe_func(fungi_database = "fungal.traits")
tITSq_s$cal_spe_func_perc(abundance_weighted = TRUE, dec=3)
#tITS_s$res_spe_func_perc=tITS_s$res_spe_func_perc[c(18,12,7)]

#tITS_s$sample_table$core_type=factor(tITS_s$sample_table$core_type, levels=c("Inner","Outer","Organic","Mineral"))
#row.names(tITSq_s$res_spe_func_perc)=c("Heartwood", "Sapwood")

data <- tITSq_s$res_spe_func_perc
filtered_data <- data[, colSums(data >= 5) > 0]
tITSq_s$res_spe_func_perc=filtered_data

tITSq_s$plot_spe_func_perc(add_facet = TRUE, use_group_list = FALSE) +
  geom_text(aes(label = round(value, 2)), color = "white", size = 5) +  # Add labels (numbers) to the heatmap
  scale_fill_gradientn(colors = c("#440154","#3b528b","#21918c", "#5ec962", "#fde725"), na.value = "grey",
                       trans = "log10",  # Logarithmic transformation
                       limits = c(0.05, NA)) +  # Automatically scale to the data range
  theme(axis.text.y = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))


#### RA Plots ####

# ```` 16S ####
ps.g=tax_glom(ps.ra, taxrank=rank_names(ps.ra)[5], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

ps.g.filt.wood=subset_samples(ps.g, TreatmentGroup=="Wood" & core_type=="Inner")
no_mito_select=subset_samples(ps.g.filt.wood, seq_id=="RM5"| seq_id=="RM502" | seq_id=="RM503")
no_mito_select=prune_taxa(taxa_sums(no_mito_select)>0.5, no_mito_select)

methods_bar_1=plot_bar(no_mito_select, x="seq_id", fill="ta5")+
  geom_bar(stat="identity")+
  scale_fill_viridis_d(name="Phylum")+
  ylab("% Abundance")+
  theme(axis.text.x=element_blank(), #remove x axis labels
         axis.ticks.x=element_blank(),
         axis.title.x=element_blank())+
  scale_fill_manual(values=c(coul_army[c(1,3,5,6,7,16,15,18,19,26)],"black"), name="Phylum")
methods_bar_1

# ```` ITS ####
ps.g=tax_glom(ps.its.ra, taxrank=rank_names(ps.its.ra)[6], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
krs=rownames(tax_table(ps.g)[which(grepl("_Pholiota*",tax_table(ps.g)[,6], ignore.case = TRUE))])
ps.krs=prune_taxa(c(krs),ps.g)

ps.g.filt.wood=subset_samples(ps.g, core_type=="Inner")
#no_mito_select=subset_samples(ps.g.filt.wood, seq_id=="RM5"| seq_id=="RM502" | seq_id=="RM503")
ps.g.filt.wood=prune_taxa(taxa_sums(ps.g.filt.wood)>10, ps.g.filt.wood)

methods_bar_1=plot_bar(ps.krs, x="seq_id", fill="ta6")+
  geom_bar(stat="identity")+
   facet_grid(.~species.x, scales="free_x")+
  ylab("% Abundance")+
  ylim(0,100)+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank())
methods_bar_1




# ```` QUVE 16S ####

ps.quve.16s=tax_glom(no_mito_q16, taxrank=rank_names(no_mito_q16)[6], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
ps.q16s.filt=prune_taxa(taxa_sums(ps.quve.16s)>30, ps.quve.16s)
ps.q16s.filt=subset_samples(ps.q16s.filt, material!="NA")

plot_bar(ps.q16s.filt, x="TreatmentGroup", fill="ta4")+
  geom_bar(stat="identity")+scale_fill_manual(values=coul[10:50])+
  facet_wrap(core_type~material, scales="free_x", ncol=5)

# ```` ITS ####

ps.quve.its=tax_glom(ps.quve.its.ra, taxrank=rank_names(ps.quve.its.ra)[6], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
ps.q.its.filt=prune_taxa(taxa_sums(ps.quve.its)>10, ps.quve.its)
ps.q.its.filt=subset_samples(ps.q.its.filt, material!="NA")

plot_bar(ps.q.its.filt, x="TreatmentGroup", fill="Class")+
  geom_bar(stat="identity")+scale_fill_manual(values=coul[10:50])+
  facet_wrap(core_type~material, scales="free_x", ncol=5)



#### RA Averages ####

wood.inner.no.ra=subset_samples(ps.ra, TreatmentGroup=="Wood" & core_type=="Inner")
ps.species=merge_samples(wood.inner.no.ra, group = sam_data(wood.inner.no.ra)$species.y)
ps.species = transform_sample_counts(ps.species, function(x) 100 * x/sum(x))
sam_data(ps.species)$name=rownames(sam_data(ps.species))
ps.species.g=tax_glom(ps.species, taxrank=rank_names(ps.species)[3], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
ps.species.g.filt=prune_taxa(taxa_sums(ps.species.g)>10, ps.species.g)

species_pal=color_taxa(ps.ra,"ta3")

inner_species=plot_bar(ps.species.g.filt, x="name", fill="ta3")+
  geom_bar(stat="identity", col="black")+scale_fill_manual(values=species_pal)
inner_species



wood.outer.no.ra=subset_samples(ps.ra, TreatmentGroup=="Wood" & core_type=="Outer")
ps.outer.species=merge_samples(wood.outer.no.ra, group = sam_data(wood.outer.no.ra)$species.y)
ps.outer.species = transform_sample_counts(ps.outer.species, function(x) 100 * x/sum(x))
sam_data(ps.outer.species)$name=rownames(sam_data(ps.outer.species))
ps.outer.g=tax_glom(ps.outer.species, taxrank=rank_names(ps.outer.species)[3], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
ps.outer.g.filt=prune_taxa(taxa_sums(ps.outer.g)>10, ps.outer.g)

outer_species=plot_bar(ps.outer.g.filt, x="name", fill="ta3")+
  geom_bar(stat="identity")+scale_fill_manual(values=species_pal)

cowplot::plot_grid(inner_species, outer_species, nrow=2)




#### RA Soil ####


ps.g.filt.soil=subset_samples(ps.g.filt, material=="Soil")
plot_bar(ps.g.filt, x="ForwardFastqFile", fill="ta3")+
  geom_bar(stat="identity") +scale_fill_manual(values=coul)+
  theme(axis.text.x=element_blank())+facet_wrap(.~species.y*core_type, scales="free", nrow=2)




soil.mineral.no.ra=subset_samples(no_mito, TreatmentGroup=="Soil" & core_type=="Mineral")
ps.mineral.species=merge_samples(soil.mineral.no.ra, group = sam_data(soil.mineral.no.ra)$species.y)
ps.mineral.species = transform_sample_counts(ps.mineral.species, function(x) 100 * x/sum(x))
sam_data(ps.mineral.species)$name=rownames(sam_data(ps.mineral.species))
ps.mineral.g=tax_glom(ps.mineral.species, taxrank=rank_names(ps.mineral.species)[5], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
ps.mineral.g.filt=prune_taxa(taxa_sums(ps.mineral.g)>3, ps.mineral.g)

species_pal_soil=color_taxa(ps.mineral.species,"ta3")

mineral_species=plot_bar(ps.mineral.g.filt, x="name", fill="ta3")+
  geom_bar(stat="identity")+scale_fill_manual(values=species_pal_soil)


soil.organic.no.ra=subset_samples(no_mito, TreatmentGroup=="Soil" & core_type=="Organic")
ps.organic.species=merge_samples(soil.organic.no.ra, group = sam_data(soil.organic.no.ra)$species.y)
ps.organic.species = transform_sample_counts(ps.organic.species, function(x) 100 * x/sum(x))
sam_data(ps.organic.species)$name=rownames(sam_data(ps.organic.species))
ps.organic.g=tax_glom(ps.organic.species, taxrank=rank_names(ps.organic.species)[5], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
ps.organic.g.filt=prune_taxa(taxa_sums(ps.organic.g)>3, ps.organic.g)

organic_species=plot_bar(ps.organic.g.filt, x="name", fill="ta3")+
  geom_bar(stat="identity")+scale_fill_manual(values=species_pal_soil)

cowplot::plot_grid(mineral_species, organic_species)


####Taxa Graphs ####
# wyatt, change 'prime' to ' before metadata merge


#### Methanogens, Family
methano=rownames(tax_table(ps.g)[which(grepl("Methano*",tax_table(ps.g)[,5], ignore.case = TRUE))])
methano_rice=rownames(tax_table(ps.g)[which(grepl("Rice Cluster*",tax_table(ps.g)[,5], ignore.case = TRUE))])
methano_bath=rownames(tax_table(ps.g)[which(grepl("Bath*",tax_table(ps.g)[,3], ignore.case = TRUE))])
ps.methano_fam=prune_taxa(c(methano, methano_bath, methano_rice),ps.g)

ps.methano_fam_s=subset_samples(ps.methano_fam, seq_id=="RM5"| seq_id=="RM502" | seq_id=="RM503") #|species.x=="BEAL"|species.x=="PIST")
sam_data(ps.methano_fam_s)$core_type=factor(sam_data(ps.methano_fam_s)$core_type, level=c("Inner","Outer","Mineral","Organic"))

ps.methano_fam_s= prune_taxa(taxa_sums(ps.methano_fam_s)>0.5, ps.methano_fam_s)
ps.methano_fam_s= subset_samples(ps.methano_fam_s, core_type=="Inner") 



#### Abundance stats
ps.methano_fam_ss= subset_samples(ps.methano_fam, core_type=="Inner" | core_type=="Outer") 

mean(otu_table(ps.methano_fam_ss)["ASV25852",])
sd(otu_table(ps.methano_fam_ss)["ASV25852",])
range(otu_table(ps.methano_fam_ss)["ASV25852",])

mean(colSums(otu_table(ps.methano_fam_ss)))
sd(colSums(otu_table(ps.methano_fam_ss)))
range(colSums(otu_table(ps.methano_fam_ss)))

#count(colSums(otu_table(ps.methano_fam_ss))>0)

otu_table(ps.methano_fam)=(otu_table(ps.methano_fam))

ps.methano_fam_ss=subset_taxa(ps.methano_fam_ss, ta5=="Methanomassiliicoccaceae")
# 1. Summarize Taxa Abundance
taxa_abundance <- colSums(otu_table(ps.methano_fam_ss))

# 2. Prepare Data for Analysis
# Extract species information
species_info <- sample_data(ps.methano_fam_ss)$core_type

# Combine taxa abundance with species information
data_for_analysis <- data.frame(species = species_info, abundance = taxa_abundance)

# 3. Statistical Test
# Perform ANOVA (or alternative test)
#anova_result <- aov(abundance ~ species, data = data_for_analysis)
#summary(anova_result)

#data_for_analysis$species[data_for_analysis$species!="ACSA"]="Other"

taxa_t <- ggplot(data_for_analysis, aes(x = species, y = log10(abundance+1))) +
  geom_boxplot()+
  stat_compare_means()
taxa_t






### Abundance graph
sample_data(ps.methano_fam)$core_type=factor(sample_data(ps.methano_fam)$core_type, levels = c("Inner","Outer","Mineral","Organic"))
ps.methano_fam=subset_samples(ps.methano_fam, species.x!="NA")

methano_p_fam=plot_bar(ps.methano_fam, x="seq_id", fill="ta5")  +
  geom_bar(stat="identity")+
  xlab("Plot")+ 
#  ylim(0,55)+
#  facet_grid(.~species.x*material, scales="free_x", space="free")+
  scale_fill_viridis_d(na.value="grey", name="Family", direction = -1)+ 
  facet_grid(core_type~species.x, scales="free", space="free_x")+
#  ggtitle("Archaeal: Methanogens (Family)")+
  ylab("% Abundance")+theme(axis.text.x=element_blank(), #remove x axis labels
                                                                      axis.ticks.x=element_blank(),
                                                                      axis.title.x=element_blank())+
  theme(legend.key.size = unit(0, 'inch'), #change legend key size
        legend.key.height = unit(.1, 'inch'), #change legend key height
        legend.key.width = unit(.150, 'inch'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# theme(aspect.ratio=4/4)+
#  force_panelsizes(cols = unit(1.5, "inch"))
methano_p_fam


#### Output datatable with abundances ###
otu_t <- t(otu_table(ps.methano_fam))
otu_df <- as.data.frame(otu_t)
merged_sample_data <- cbind(sample_data(ps.methano_fam), otu_df)

colnames(merged_sample_data)[81:91]=tax_table(ps.methano_fam)[,5]
colnames(merged_sample_data)[83]="Bathyarchaeia"
write.csv(merged_sample_data,"tree_data_methanogen_group.csv")




library(phyloseq)
library(dplyr)
library(ggplot2)
library(tidyr)

# Convert phyloseq object to data frame
df <- psmelt(ps.methano_fam)

# Sum the abundances of all ASVs per sample
sum_df <- df %>%
  group_by(Sample, species.x, core_type) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  ungroup()

# Calculate average and standard error
summary_df <- sum_df %>%
  group_by(species.x, core_type) %>%
  summarise(
    avg = mean(total_abundance),
    se = sd(total_abundance) / sqrt(n())
  )

# Plot
summary_df$core_type=factor(summary_df$core_type, levels=c("Inner","Outer","Mineral","Organic"),
                            labels=c("Heartwood","Sapwood","Mineral Soil","Organic Soil"))

methano_p_fam_1 <- ggplot(summary_df, aes(x = species.x, y = avg, fill = species.x)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 0.25, position = position_dodge(0.9)) +
 # xlab("Species") + 
  scale_fill_manual(na.value="grey", name="Species", values = species_colors_dd) + 
  ylab("Mean % Abundance Methanogens") +
  facet_wrap(~core_type, ncol=4) + # Separate chart for each core_type
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title.x = element_blank(),
        legend.key.size = unit(0, 'inch'), #change legend key size
        legend.key.height = unit(.1, 'inch'), #change legend key height
        legend.key.width = unit(.150, 'inch'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16),
        legend.position = "none")

methano_p_fam_1






library(phyloseq)
library(dplyr)
library(ggplot2)
library(tidyr)

# Convert phyloseq object to data frame
df <- psmelt(ps.methano_fam)

# Calculate the total abundance for each sample
df <- df %>%
  group_by(Sample, species.x, core_type) %>%
  mutate(Total = sum(Abundance)) %>%
  ungroup()

# Filter out samples with zero total abundance to avoid division by zero
#df <- df %>% filter(Total > 0)

# Calculate the relative abundance of each family-level taxon within each sample
df <- df %>%
  group_by(Sample, species.x, core_type, ta5) %>%
  summarise(RelativeAbundance = sum(Abundance) , .groups = 'drop') %>% #/ Total
  ungroup()

# Calculate the average relative abundance of each family-level taxa for each species
summary_df <- df %>%
  group_by(species.x, core_type, ta5) %>%
  summarise(avg = mean(RelativeAbundance), .groups = 'drop') %>%
  ungroup()

# Plot
summary_df$core_type <- factor(summary_df$core_type, levels = c("Inner", "Outer", "Mineral", "Organic"),
                               labels=c("Heartwood","Sapwood","Mineral Soil","Organic Soil"))
summary_df$ta5[is.na(summary_df$ta5)]="Bathyarchaeia Class"

methano_p_fam_2 <- ggplot(summary_df, aes(x = species.x, y = log10(avg+1), fill = ta5)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("Species") +
  scale_fill_viridis_d(na.value = "grey", name = "Family", direction = -1) +
  ylab(expression(Log[10](Mean~Relative~Abundance))) +
  facet_wrap(~core_type, ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.key.size = unit(0, 'inch'),
        legend.key.height = unit(.1, 'inch'),
        legend.key.width = unit(.150, 'inch'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))
methano_p_fam_2

cowplot::plot_grid(methano_p_fam_1, methano_p_fam_2, ncol=2, align = "v")




# Convert phyloseq object to data frame
df <- psmelt(ps.methano_fam)

# Calculate the total abundance for each sample
df <- df %>%
  group_by(Sample) %>%
  mutate(Total = sum(Abundance)) %>%
  ungroup()

# Filter out samples with zero total abundance to avoid division by zero
df <- df %>% filter(Total > 0)

# Calculate the relative abundance of each family-level taxon within each sample
df <- df %>%
  group_by(Sample, species.x, core_type, ta5) %>%
  summarise(RelativeAbundance = sum(Abundance) / Total, .groups = 'drop') %>%
  ungroup()

# Calculate the average relative abundance of each family-level taxa for each species
summary_df <- df %>%
  group_by(species.x, core_type, ta5) %>%
  summarise(avg = mean(RelativeAbundance), .groups = 'drop') %>%
  ungroup()

# Plot
color_labels <- function(x) {
  sapply(x, function(species) {
    color <- species_colors_dd[species]
    paste0("<span style='color:", color, ";'>", species, "</span>")
  })
}

summary_df$core_type <- factor(summary_df$core_type, levels = c("Inner", "Outer", "Mineral", "Organic"),
                               labels=c("Heartwood","Sapwood","Mineral Soil","Organic Soil"))
summary_df$ta5[is.na(summary_df$ta5)]="Bathyarchaeia Class"

methano_p_fam_2 <- ggplot(summary_df, aes(x = species.x, y = avg, fill = ta5)) +
  geom_bar(stat = "summary") +
  xlab("Species") +
  scale_fill_viridis_d(na.value = "grey", name = "Family", direction = -1) +
  ylab("Mean Relative Abundance") + # Removed the log10 transformation
  facet_wrap(~core_type, ncol = 4) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust=0.5),  # Use ggtext for axis text
        legend.key.size = unit(0, 'inch'),
        legend.key.height = unit(.1, 'inch'),
        legend.key.width = unit(.150, 'inch'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  scale_x_discrete(labels = color_labels)
methano_p_fam_2


cowplot::plot_grid(methano_p_fam_1, methano_p_fam_2, ncol=1, align = "v", axis = "tblr")





















cowplot::plot_grid(methods_bar_1, methano_p_fam, rel_widths = c(1,1.25), labels = "AUTO")
# Susbet just RA and ddPCR target
methano_matrix=as.data.frame(sample_sums(otu_table(ps.methano_fam)))
methano_matrix$mcra_probe_loose=sample_data(ps.methano_fam)$mcra_probe_loose
methano_matrix$mcra_eva_loose=sample_data(ps.methano_fam)$mcra_eva_loose
methano_matrix$type=sample_data(ps.methano_fam)$core_type
methano_matrix$PC=as.numeric(sample_data(ps.methano_fam)$ch4_50)
colnames(methano_matrix)[1]="RA"
methano_matrix_abundant=subset(methano_matrix, RA!=0 & type!="MISC")

probe_ra=ggplot(data=methano_matrix_abundant,aes(y=log(mcra_probe_loose+1, base=10),
                                               x=log(PC, base=10), col=type))+
  geom_point(size=2.5)+
 # ggtitle("mcrA Probe Validation")+
  ylab("Log10(ddPCR Methanogen Abundance)")+
  xlab("Log10(CH4 Flux 50cm)")+
  geom_smooth(method="lm", col="black")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
            digits = 2)+
#  stat_regline_equation()+
  scale_color_manual(values=coul[c(28,30,34,32)], name="Material")+
  facet_grid(.~type, scales="free_x")
probe_ra





#### Methanotrophs, Family
methyl=rownames(tax_table(ps.g)[which(grepl("Methyl*",tax_table(ps.g)[,5], ignore.case = TRUE))])
beij=rownames(tax_table(ps.g)[which(grepl("Beij*",tax_table(ps.g)[,5], ignore.case = TRUE))])
#methano_bath=rownames(tax_table(no_mito)[which(grepl("Bath*",tax_table(no_mito)[,3], ignore.case = TRUE))])
ps.methyl_fam=prune_taxa(c(methyl),ps.g)

methyl_p_fam=plot_bar(ps.methyl_fam, x="species.x", fill="ta5")  +
  geom_bar(stat="identity")+
  xlab("Plot")+ 
  scale_fill_manual(values=coul, na.value="grey")+ 
  facet_wrap(core_type~., scales="free_x")+
  ggtitle("Methanotrophs (Family)")+ylab("% Abundance")+theme(axis.text.x=element_blank(), #remove x axis labels
                                                                      axis.ticks.x=element_blank(),
                                                                      axis.title.x=element_blank())+
  facet_grid(core_type~species.x, scales="free", space="free")+
  theme(legend.key.size = unit(0, 'inch'), #change legend key size
        legend.key.height = unit(.1, 'inch'), #change legend key height
        legend.key.width = unit(.150, 'inch'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=16),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#  force_panelsizes(cols = unit(1.5, "inch"))
methyl_p_fam




# Convert phyloseq object to data frame
df <- psmelt(ps.methyl_fam)

# Sum the abundances of all ASVs per sample
sum_df <- df %>%
  group_by(Sample, species.x, core_type) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  ungroup()

# Calculate average and standard error
summary_df <- sum_df %>%
  group_by(species.x, core_type) %>%
  summarise(
    avg = mean(total_abundance),
    se = sd(total_abundance) / sqrt(n())
  )

# Plot
summary_df$core_type=factor(summary_df$core_type, levels=c("Inner","Outer","Mineral","Organic"))

methyl_p_fam <- ggplot(summary_df, aes(x = species.x, y = avg, fill = species.x)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 0.25, position = position_dodge(0.9)) +
  xlab("Species") + 
  scale_fill_viridis_d(na.value="grey", name="Species", direction = -1, option = 3) + 
  stat_summary(fun = "mean", geom = "text", aes(label = round(..y.., 2)), vjust = -1) +
  ylab("Mean % Abundance Methanotrophs") +
  facet_wrap(~core_type, ncol=2) + # Separate chart for each core_type
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0, 'inch'), #change legend key size
        legend.key.height = unit(.1, 'inch'), #change legend key height
        legend.key.width = unit(.150, 'inch'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=16))

methyl_p_fam



ps.methyl_fam=subset_samples(ps.methyl_fam, species.x!="NA")
# Convert phyloseq object to data frame
df <- psmelt(ps.methyl_fam)

# Calculate the total abundance for each sample
df <- df %>%
  group_by(Sample, species.x, core_type) %>%
  mutate(Total = sum(Abundance)) %>%
  ungroup()

# Filter out samples with zero total abundance to avoid division by zero
#df <- df %>% filter(Total > 0)

# Calculate the relative abundance of each family-level taxon within each sample
df <- df %>%
  group_by(Sample, species.x, core_type, ta5) %>%
  summarise(RelativeAbundance = sum(Abundance) , .groups = 'drop') %>%
  ungroup()

# Calculate the average relative abundance of each family-level taxa for each species
summary_df <- df %>%
  group_by(species.x, core_type, ta5) %>%
  summarise(avg = mean(RelativeAbundance), .groups = 'drop') %>%
  ungroup()

# Plot
summary_df$core_type <- factor(summary_df$core_type, levels = c("Inner", "Outer", "Mineral", "Organic"))

methyl_p_fam_2 <- ggplot(summary_df, aes(x = species.x, y = log10(avg+1), fill = ta5)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("Species") +
  scale_fill_viridis_d(na.value = "grey", name = "Family", direction = -1, option = "C") +
  ylab(expression(Log[10](Mean~Relative~Abundance))) +
  facet_wrap(~core_type, ncol = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.key.size = unit(0, 'inch'),
        legend.key.height = unit(.1, 'inch'),
        legend.key.width = unit(.150, 'inch'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))
methyl_p_fam_2











methyl_matrix=as.data.frame(sample_sums(otu_table(ps.methyl_fam)))
methyl_matrix$pmoa_loose=sample_data(ps.methyl_fam)$pmoa_loose
methyl_matrix$pmoa_loose=sample_data(ps.methyl_fam)$pmoa_loose
methyl_matrix$type=sample_data(ps.methyl_fam)$TreatmentGroup
colnames(methyl_matrix)[1]="RA"
methyl_matrix_abundant=subset(methyl_matrix, RA!=0 & type!="MISC")

primer_ra=ggplot(data=methyl_matrix_abundant,aes(y=log(RA, base=10),
                                                 x=log(pmoa_loose+1, base=10), col=type))+
  geom_point(size=2.5)+
  ggtitle("mcrA Probe Validation")+
  ylab("Log10(Total methylgen RA)")+
  xlab("Log10(ddPCR Copies/rxn)")+
  geom_smooth(method="lm", col="black")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 0.75, digits = 2)+
  #  stat_regline_equation()+
  scale_color_manual(values=coul[c(28,30)], name="Material")
primer_ra



#### Ordinations ####
# ```` 16S ####

#ps.wood=subset_samples(ps.w, TreatmentGroup=="Wood")
#ord.nmds.bray <- ordinate(ps.wood, method="NMDS", distance="bray")

sam_data(ps.rare)$mcra10=log10(sam_data(ps.rare)$mcra_probe_loose+1)
sam_data(ps.rare)$mmox10=log10(sam_data(ps.rare)$mmox_loose+1)
sam_data(ps.rare)$Percent_Change=as.numeric(sam_data(ps.rare)$Percent_Change)
sam_data(ps.rare)$Percent_Change=ifelse(sam_data(ps.rare)$Percent_Change==1 | sam_data(ps.rare)$Percent_Change<0,NA,sam_data(ps.rare)$Percent_Change)
sam_data(ps.rare)$Soft_Hard=ifelse(sam_data(ps.rare)$species.x=="TSCA"|sam_data(ps.rare)$species.x=="PIST","Soft","Hard")
sam_data(ps.rare)$CO2_int=log10(sam_data(ps.rare)$CO2_int)
sam_data(ps.rare)$CH4_int=log10(sam_data(ps.rare)$CH4_int)
sam_data(ps.rare)$O2_int=log10(sam_data(ps.rare)$O2_int)
sam_data(ps.rare)$N2O_int=log10(sam_data(ps.rare)$N2O_int)


#p=plot_ordination(ps.wood, ord.nmds.bray, color="species.x",title="Bray NMDS (16S Universal)")
#p+geom_point(size=3, aes(shape=sam_data(ps.wood)$Soft_Hard))+scale_color_viridis_d()

set.seed(46814)

set3_palette <- brewer.pal(12, "Set3")
dark2_palette <- brewer.pal(8, "Dark2")

# We make sure not to repeat colors from Set3 that might also appear in Dark2
combined_palette <- c(set3_palette, setdiff(dark2_palette, set3_palette))
species_colors <- setNames(rep(combined_palette, length.out = length(unique(sam_data(ps.rare)$species.x))), unique(sam_data(ps.rare)$species.x))

no_mito_wood=subset_samples(ps.rare, core_type=="Inner" & species.x!="NA")
bact_umap=umap(as.matrix(distance(no_mito_wood, method = "unifrac")))


sam_data(no_mito_wood)$Percent_Change=as.numeric(sam_data(no_mito_wood)$Percent_Change)
sam_data(no_mito_wood)$Percent_Change=ifelse(sam_data(no_mito_wood)$Percent_Change==1 | sam_data(no_mito_wood)$Percent_Change<0,NA,sam_data(no_mito_wood)$Percent_Change)
sam_data(no_mito_wood)$Soft_Hard=ifelse(sam_data(no_mito_wood)$species.x=="TSCA"|sam_data(no_mito_wood)$species.x=="PIST","Soft","Hard")






library(ggplot2)
library(dplyr)
create_umap_scatterplot <- function(umap_output, color_data) {
  umap_df <- as.data.frame(umap_output$layout)
  umap_df$color_variable <- color_data

  
  # Calculate centroids and standard deviations for each color group (species)
  centroids <- umap_df %>%
    group_by(color_variable) %>%
    summarise(
      sd_V1 = sd(V1),
      sd_V2 = sd(V2),
      V1 = mean(V1),
      V2 = mean(V2)
    )
  
  # Calculate a combined standard deviation for the size (e.g., Euclidean distance)
  centroids <- centroids %>%
    mutate(size = sqrt(sd_V1^2 + sd_V2^2))
  
  # Print the centroids to check their values
  print(centroids)
  
  "
  # Normalize sizes to a specified range
  size_range <- c(3, 10)  # Adjust size range as needed
  centroids$size <- rescale(centroids$size, to = size_range)
  centroids$size = centroids$size*10
  "
  
  # Create the scatter plot with centroids
  ggplot(umap_df, aes(x = V1, y = V2)) +
    geom_point(aes(col = color_variable), size = 5, alpha=0.9) +
    #    geom_point(data = centroids, aes(x = V1, y = V2, col = "black"), size=6, 
    #               alpha = 1, shape=4)+
    #    geom_point(data = centroids, aes(x = V1, y = V2, col = color_variable), size=4, 
    #             alpha = 1, shape=18)+
    scale_color_manual(name = "Species Code") +
    theme_minimal(base_size = 16) +
    guides(color = guide_legend(ncol = 2, byrow = TRUE))
}

heart_map=create_umap_scatterplot(bact_umap, sam_data(no_mito_wood)$species.x)+
  scale_color_manual(values = species_colors, name="Species Code")+
  ggtitle("UMAP 16S Sapwood (UniFrac)")+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.position = "bottom")
#  scale_color_manual(values=species_colors_dd)



# Define the get_legend function
get_legend <- function(myggplot) {
  tmp <- ggplotGrob(myggplot)
  leg <- gtable::gtable_filter(tmp, "guide-box")
  return(leg)
}

# Assuming sap_map and heart_map are your ggplot objects

# Create a shared legend plot
shared_legend <- get_legend(heart_map)


# Remove the legend from both plots
sap_map <- sap_map + theme(legend.position = "none")
heart_map <- heart_map + theme(legend.position = "none")


# Combine the plots
grid.arrange(
  arrangeGrob(sap_map, heart_map, nrow = 1),
  shared_legend,
  ncol = 2,
  widths = c(5, 1) # Adjust widths as necessary
)


# Example usage:
# Assuming umap_output is the result from UMAP, color_data is the species information, 
# and shape_data is the shape variable.

# create_umap_scatterplot(umap_output, sam_data(no_mito_wood)$species.x, sam_data(no_mito_wood)$shape_var)
# ```` ITS ####

species_colors_its <- setNames(rep(combined_palette, length.out = length(unique(sam_data(ps.its.rare)$species.x))), unique(sam_data(ps.its.rare)$species.x))

ps.its.wood=subset_samples(ps.its.rare, core_type=="Inner" & species.x != "NA")
#ord.nmds.bray <- ordinate(ps.its.wood, method="PCoA", distance="bray")

sam_data(ps.its.wood)$mcra10=log10(sam_data(ps.its.wood)$mcra_probe_loose+1)
sam_data(ps.its.wood)$mmox10=log10(sam_data(ps.its.wood)$mmox_loose+1)
sam_data(ps.its.wood)$Percent_Change=as.numeric(sam_data(ps.its.wood)$Percent_Change)
sam_data(ps.its.wood)$Percent_Change=ifelse(sam_data(ps.its.wood)$Percent_Change==1 | sam_data(ps.its.wood)$Percent_Change<0,NA,sam_data(ps.its.wood)$Percent_Change)
sam_data(ps.its.wood)$Soft_Hard=ifelse(sam_data(ps.its.wood)$species.x=="TSCA"|sam_data(ps.its.wood)$species.x=="PIST","Soft","Hard")
sam_data(ps.its.wood)$CO2_int=log10(sam_data(ps.its.wood)$CO2_int)
sam_data(ps.its.wood)$CH4_int=log10(sam_data(ps.its.wood)$CH4_int)
sam_data(ps.its.wood)$O2_int=log10(sam_data(ps.its.wood)$O2_int)
sam_data(ps.its.wood)$N2O_int=log10(sam_data(ps.its.wood)$N2O_int)

p=plot_ordination(ps.its.wood, ord.nmds.bray, color="Percent_Change",title="Bray NMDS (ITS)")
p+geom_point(size=3)

library(umap)

its_umap=umap(as.matrix(distance(ps.its.wood, method = "unifrac")))



library(ggplot2)
library(dplyr)
create_umap_scatterplot <- function(umap_output, color_data) {
  umap_df <- as.data.frame(umap_output$layout)
  umap_df$color_variable <- color_data
  
  
  # Calculate centroids and standard deviations for each color group (species)
  centroids <- umap_df %>%
    group_by(color_variable) %>%
    summarise(
      sd_V1 = sd(V1),
      sd_V2 = sd(V2),
      V1 = mean(V1),
      V2 = mean(V2)
    )
  
  # Calculate a combined standard deviation for the size (e.g., Euclidean distance)
  centroids <- centroids %>%
    mutate(size = sqrt(sd_V1^2 + sd_V2^2))
  
  # Print the centroids to check their values
  print(centroids)
  
  "
  # Normalize sizes to a specified range
  size_range <- c(3, 10)  # Adjust size range as needed
  centroids$size <- rescale(centroids$size, to = size_range)
  centroids$size = centroids$size*10
  "
  
  # Create the scatter plot with centroids
  ggplot(umap_df, aes(x = V1, y = V2)) +
    geom_point(aes(col = color_variable), size = 5, alpha=0.9) +
    #    geom_point(data = centroids, aes(x = V1, y = V2, col = "black"), size=6, 
    #               alpha = 1, shape=4)+
    #    geom_point(data = centroids, aes(x = V1, y = V2, col = color_variable), size=4, 
    #             alpha = 1, shape=18)+
    scale_color_manual(name = "Species Code") +
    theme_minimal(base_size = 16) +
    guides(color = guide_legend(ncol = 2, byrow = TRUE))
}

heart_map_its=create_umap_scatterplot(its_umap, sam_data(ps.its.wood)$species.x)+
  scale_color_manual(values = species_colors, name="Species Code")+
  ggtitle("UMAP ITS Heartwood (UniFrac)")+
  guides(color=guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.position = "bottom") 



# Define the get_legend function
get_legend <- function(myggplot) {
  tmp <- ggplotGrob(myggplot)
  leg <- gtable::gtable_filter(tmp, "guide-box")
  return(leg)
}

# Assuming sap_map and heart_map are your ggplot objects

# Create a shared legend plot
shared_legend <- get_legend(heart_map_its)


# Remove the legend from both plots
sap_map_its <- sap_map_its + theme(legend.position = "none")
heart_map_its <- heart_map_its + theme(legend.position = "none")


# Combine the plots
grid.arrange(
  arrangeGrob(sap_map_its, heart_map_its, nrow = 1),
  shared_legend,
  ncol = 2,
  widths = c(5, 1) # Adjust widths as necessary
)

##### ```` QUVE #### 


ord.nmds.bray <- ordinate(no_mito_q16, method="PCoA", distance="bray")
p=plot_ordination(ps.quve.16s.rare, ord.nmds.bray, color="core_type",title="Bray NMDS (QUVE)")
p+geom_point(aes(size=sample_data(ps.quve.16s.rare)$Inner.Core.Sample.ID))+
  scale_color_manual(values=coul[6:25])

#### Tree phylogeny ####

library(phytools)

S.PhyloMaker<-function (tree, spList, nodes, output.spList = T, scenarios = c("S1", "S2", "S3")) 
{
  options(scipen=999)
  tree0 <- tree
  spList[sapply(spList, is.factor)] <- lapply(spList[sapply(spList, is.factor)], as.character)
  if (any(duplicated(spList$species))) 
  {
    warning("Duplicated species detected and removed.")
    print(spList$species[duplicated(spList$species)])
  }
  spList <- spList[!duplicated(spList$species), ]
  spList.original <- spList
  spList$species <- gsub(" ", "_", spList$species)
  spList$species <- gsub("(^[[:alpha:]])", "\\U\\1", spList$species, perl = TRUE)
  spList$genus <- gsub("(^[[:alpha:]])", "\\U\\1", spList$genus, perl = TRUE)
  spList$family <- gsub("(^[[:alpha:]])", "\\U\\1", spList$family, perl = TRUE)
  rnN <- data.frame(node.label = paste("N", 1:length(tree$node.label), sep = ""), oriN = tree$node.label, stringsAsFactors = FALSE)
  nodes[,c("level","family","genus","rn","bn","taxa")]<-lapply(nodes[,c("level","family","genus","rn","bn","taxa")], as.character)
  tree$node.label <- paste("N", 1:length(tree$node.label), sep = "")
  kk<-c()
  for (i in 1:length(tree$tip.label)) {
    kk<-c(kk,substring(tree$tip.label[i],1,gregexpr("_",tree$tip.label[i])[[1]][1]-1))
  }
  m<-data.frame(num=1:length(kk),genus=kk,species=tree$tip.label)
  m<-merge(m,nodes[,c("genus","family")])
  mX <- m
  m <- m[,c("genus","family")]
  m <- m[!duplicated(m$genus),]
  dimnames(m)[[2]][2] <- "family_in_PhytoPhylo"
  m<-m[,c("genus","family_in_PhytoPhylo")]
  m0 <- spList[!duplicated(spList$genus), c("genus","family")]
  dimnames(m0)[[2]][2] <- "family_in_spList"
  mm<-merge(m0, m)
  g<-mm[which(is.na(match(paste(mm$genus,mm$family_in_spList,sep="_"),paste(mm$genus,mm$family_in_PhytoPhylo,sep="_")))),]
  if (dim(g)[1]>0)
  {
    print("Taxonomic classification not consistent between spList and PhytoPhylo.")
    print(g) 
  }
  add.tip <- spList[which(is.na(match(spList$species, tree$tip.label))), ]
  status <- rep("match(prune)", dim(spList)[1])
  status[which(is.na(match(spList$species, tree$tip.label)))] <- "match(add)"
  if (dim(add.tip)[1] == 0 & length(na.omit(match(spList$species, tree$tip.label))) == 0)
    stop("Incorrect format of species list.")
  if (length(setdiff(spList$species, tree0$tip.label)) == 0 & length(na.omit(match(spList$species, tree$tip.label))) > 0)
  {
    print("There is no species needs to be added, all the species are pruned from PhytoPhylo.")
    splis <- spList.original
    treeX <- drop.tip(tree0, setdiff(tree0$tip.label, splis$species))
    splis$status <- "match(prune)"
    phylo0 <- list(Scenario.1 = NULL, Scenario.2 = NULL, Scenario.3 = NULL, Species.list = splis)
    if ("S1" %in% scenarios) {phylo0$Scenario.1 <- treeX}
    if ("S2" %in% scenarios) {phylo0$Scenario.2 <- treeX}
    if ("S3" %in% scenarios) {phylo0$Scenario.3 <- treeX}
    phylo0[sapply(phylo0, is.null)] <- NULL
    return(phylo0)
    stop()
  }
  add.tip$sort <- ""
  add.tip$sort[which(!is.na(match(add.tip$genus, nodes[nodes$level == "G", ]$genus)))] <- "G1"
  add.tip$sort[which(is.na(match(add.tip$genus, nodes[nodes$level == "G", ]$genus)) & !is.na(match(add.tip$family, nodes[nodes$level == "F", ]$family)))] <- "F1"
  add.tip$sort[add.tip$sort == "F1"][duplicated(add.tip[add.tip$sort == "F1", ]$genus)] <- "F2"
  a <- which(add.tip$sort == "")
  if (length(a) > 0)
  {  
    print(paste("Note:", length(a), "taxa unmatch:",sep=" "))
    print(add.tip$species[a])
    status[match(add.tip$species[a], spList$species)] <- "unmatch"         
  }  
  spList.original$status <- status
  if ("S1" %in% scenarios) {
    t1 <- tree
    rnN1<-rnN
    nG <- nodes[nodes$level == "G", ]
    nF <- nodes[nodes$level == "F", ]
    data <- add.tip[add.tip$sort == "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        n <- match(data$family[i], nF$family)
        g <- nF$gen.n[n]
        s <- nF$sp.n[n]
        if (g == 1 & s == 1) {                                                                          
          num <- grep(nF$taxa[n], t1$tip.label)
          len <- t1$edge.length[match(num, t1$edge[, 2])]
          t1 <- bind.tip(t1, tip.label = data$species[i], edge.length = len, where = num, position = len)
          nF$gen.n[n] <- g + 1
          nF$sp.n[n] <- s + 1
          num <- grep(nF$taxa[n], t1$tip.label)
          t1$node.label[match(t1$edge[match(num, t1$edge[, 2]), 1], unique(t1$edge[, 1]))] <- paste("NN", t1$Nnode + 1, sep = "")
          rnN1$node.label[match(nF$bn[n], rnN1$node.label)]<- paste("NN", t1$Nnode + 1, sep = "")
          nF$bn[n] <- paste("NN", t1$Nnode + 1, sep = "")
          nF$bn.bl[n] <- len
        }             
        else {
          num <- unique(t1$edge[, 1])[match(nF$bn[n], t1$node.label)]
          len <- nF$bn.bl[n]
          t1 <- bind.tip(t1, tip.label = data$species[i], edge.length = len, where = num)
        }
      }
    }
    data <- add.tip[add.tip$sort != "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        n <- grep(paste(data$genus[i], "_", sep = ""), t1$tip.label)
        if (length(n) == 1) {
          num <- n
          len <- t1$edge.length[match(num, t1$edge[, 2])]
          t1 <- bind.tip(t1, tip.label = data$species[i], edge.length = len, where = num, position = len)
        }
        if (length(n) > 1) {
          num <- fastMRCA(t1, t1$tip.label[min(n)], t1$tip.label[max(n)])
          len <- fastDist(t1, t1$tip.label[min(n)], t1$tip.label[max(n)])/2
          t1 <- bind.tip(t1, tip.label = data$species[i], edge.length = len, where = num)
        }   
      }
    }
    toDrop <- setdiff(1:length(t1$tip.label), which(!is.na(match(t1$tip.label, spList$species))))
    t1 <- drop.tip(t1, tip = toDrop)
    Re <- which(!is.na(match(t1$node.label, rnN1$node.label)))
    noRe <- which(is.na(match(t1$node.label, rnN1$node.label)))
    t1$node.label[Re] <- rnN1$oriN[match(t1$node.label, rnN1$node.label)[Re]]
    t1$node.label[noRe] <- ""
  }
  else {
    t1 <- NULL
  }
  if ("S2" %in% scenarios) {
    t2 <- tree
    rnN2<-rnN
    nG <- nodes[nodes$level == "G", ]
    nF <- nodes[nodes$level == "F", ]
    data <- add.tip[add.tip$sort == "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        n <- match(data$family[i], nF$family)
        g <- nF$gen.n[n]
        s <- nF$sp.n[n]
        if (g == 1 & s == 1) {
          num <- grep(nF$taxa[n], t2$tip.label)
          len <- t2$edge.length[match(num, t2$edge[,2])] * sample((1:99)/100,1)
          t2 <- bind.tip(t2, tip.label = data$species[i], edge.length = len, where = num, position = len)    
          nF$gen.n[n] <- g + 1
          nF$sp.n[n] <- s + 1
          num <- grep(data$species[i], t2$tip.label)
          t2$node.label[match(t2$edge[match(num, t2$edge[, 2]), 1], unique(t2$edge[,1]))] <- paste("NN", t2$Nnode + 1, sep = "")
          rnN2$node.label[match(nF$bn[n], rnN2$node.label)]<- paste("NN", t2$Nnode + 1, sep = "")
          nF$bn[n] <- paste("NN", t2$Nnode + 1, sep = "")
          nF$bn.bl[n] <- len
        }
        else {
          num <- unique(t2$edge[, 1])[match(nF$bn[n], t2$node.label)]
          len <- t2$edge.length[match(num,t2$edge[,2])] * sample((1:99)/100,1)
          t2 <- bind.tip(t2, tip.label = data$species[i], edge.length = len, where = num, position = len)
          nF$gen.n[n] <- g + 1
          nF$sp.n[n] <- s + 1
          num <- grep(data$species[i], t2$tip.label)
          t2$node.label[match(t2$edge[match(num, t2$edge[, 2]), 1], unique(t2$edge[, 1]))] <- paste("NN", t2$Nnode + 1, sep = "")
          rnN2$node.label[match(nF$bn[n], rnN2$node.label)]<- paste("NN", t2$Nnode + 1, sep = "")
          nF$bn[n] <- paste("NN", t2$Nnode + 1, sep = "")
          nF$bn.bl[n] <- nF$bn.bl[n]+len
        }
      }
    }
    data <- add.tip[add.tip$sort != "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        n <- grep(paste(data$genus[i], "_", sep = ""), t2$tip.label)
        if (length(n) == 1) {
          num <- n
          len <- t2$edge.length[match(num, t2$edge[, 2])] * sample((1:99)/100,1)
          t2 <- bind.tip(t2, tip.label = data$species[i], edge.length = len, where = num, position = len)
        }
        if (length(n) > 1) {
          num <- sample(n,1)
          len <- t2$edge.length[match(num, t2$edge[, 2])] * sample((1:99)/100,1)
          t2 <- bind.tip(t2, tip.label = data$species[i], edge.length = len, where = num, position = len)
        }   
      }
    }
    toDrop <- setdiff(1:length(t2$tip.label), which(!is.na(match(t2$tip.label, spList$species))))
    t2 <- drop.tip(t2, tip = toDrop)
    Re <- which(!is.na(match(t2$node.label, rnN2$node.label)))
    noRe <- which(is.na(match(t2$node.label, rnN2$node.label)))
    t2$node.label[Re] <- rnN2$oriN[match(t2$node.label, rnN2$node.label)[Re]]
    t2$node.label[noRe] <- ""
  }
  else {
    t2 <- NULL
  }
  if ("S3" %in% scenarios) {
    t3 <- tree
    rnN3<-rnN
    nG <- nodes[nodes$level == "G", ]
    nF <- nodes[nodes$level == "F", ]
    data <- add.tip[add.tip$sort == "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        n <- match(data$family[i], nF$family)
        g <- nF$gen.n[n]
        s <- nF$sp.n[n]
        if (g == 1 & s == 1) {                                                                          
          num <- grep(nF$taxa[n], t3$tip.label)
          len <- t3$edge.length[match(num, t3$edge[, 2])] * (2/3)
          t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = num, position = len)
          nF$gen.n[n] <- g + 1
          nF$sp.n[n] <- s + 1
          num <- grep(nF$taxa[n], t3$tip.label)
          t3$node.label[match(t3$edge[match(num, t3$edge[, 2]), 1], unique(t3$edge[, 1]))] <- paste("NN", t3$Nnode + 1, sep = "")
          rnN3$node.label[match(nF$bn[n], rnN3$node.label)]<- paste("NN", t3$Nnode + 1, sep = "")
          nF$bn[n] <- paste("NN", t3$Nnode + 1, sep = "")
          nF$bn.bl[n] <- len
        }
        if (g == 1 & s > 1) {
          num <- unique(t3$edge[, 1])[match(nF$bn[n], t3$node.label)]
          if ((2/3)*nF$rn.bl[n] <= nF$bn.bl[n])  { len <-(nF$rn.bl[n]-nF$bn.bl[n])/2 }
          if ((2/3)*nF$rn.bl[n] > nF$bn.bl[n])   { len <-nF$rn.bl[n]*2/3-nF$bn.bl[n] }
          t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = num, position = len)
          nF$gen.n[n] <- g + 1
          nF$sp.n[n] <- s + 1
          num <- unique(t3$edge[, 1])[match(nF$bn[n], t3$node.label)]
          t3$node.label[match(t3$edge[match(num, t3$edge[, 2]), 1], unique(t3$edge[, 1]))] <- paste("NN", t3$Nnode + 1, sep = "")
          rnN3$node.label[match(nF$bn[n], rnN3$node.label)]<- paste("NN", t3$Nnode + 1, sep = "")
          nF$bn[n] <- paste("NN", t3$Nnode + 1, sep = "")
          nF$bn.bl[n] <- len+nF$bn.bl[n]
        }
        if (g > 1) {
          num <- unique(t3$edge[, 1])[match(nF$bn[n], t3$node.label)]
          len <- nF$bn.bl[n]
          t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = num)
        }
      }
    }
    data <- add.tip[add.tip$sort != "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        n <- grep(paste(data$genus[i], "_", sep = ""), t3$tip.label)
        if (length(n) == 1) {
          len <- t3$edge.length[match(n, t3$edge[, 2])]/2
          t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = n, position = len)
          nG$sp.n[match(data$genus[i], nG$genus)] <- length(n) + 1
          nu <- grep(paste(data$genus[i], "_", sep = ""), t3$tip.label)
          num <- fastMRCA(t3, t3$tip.label[nu[1]], t3$tip.label[nu[2]])
          t3$node.label[match(num, unique(t3$edge[, 1]))] <- paste("NN", t3$Nnode + 1,  sep = "")
          rnN3$node.label[match(nG$bn[n], rnN3$node.label)]<- paste("NN", t3$Nnode + 1, sep = "")
          nG$bn[match(data$genus[i], nG$genus)] <- paste("NN", t3$Nnode + 1, sep = "")
          nG$bn.bl[match(data$genus[i], nG$genus)] <- len
        }
        if (length(n) > 1) {
          num <- fastMRCA(t3, t3$tip.label[min(n)], t3$tip.label[max(n)])
          len <- as.numeric(branching.times(t3))[match(num, unique(t3$edge[, 1]))]
          t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = num)
        }
      }
    }
    toDrop <- setdiff(1:length(t3$tip.label), which(!is.na(match(t3$tip.label, spList$species))))
    t3 <- drop.tip(t3, tip = toDrop)
    Re <- which(!is.na(match(t3$node.label, rnN3$node.label)))
    noRe <- which(is.na(match(t3$node.label, rnN3$node.label)))
    t3$node.label[Re] <- rnN3$oriN[match(t3$node.label, rnN3$node.label)[Re]]
    t3$node.label[noRe] <- ""
  }
  else {
    t3 <- NULL
  }
  if (output.spList == FALSE) 
    spList <- NULL
  phylo <- list(Scenario.1 = t1, Scenario.2 = t2, Scenario.3 = t3, Species.list = spList.original)
  phylo[sapply(phylo, is.null)] <- NULL
  return(phylo)
}

example<-read.csv("/Users/Wyatt/Documents/S.PhyloMaker-master/example.splist.csv",header=T) # read in the example species list.
samples= read.csv("/Users/Wyatt/Desktop/UMGC Tree Full Sequences/phylo_trees_selections.csv", header=T)

phylo<-read.tree("/Users/Wyatt/Documents/PhytoPhylo.tre") # read in the megaphylogeny.
nodes<-read.csv(file = "/Users/Wyatt/Documents/S.PhyloMaker-master/nodes.csv",header=T) # read in the nodes information of the megaphylogeny.

result<-S.PhyloMaker(spList=samples, tree=phylo, nodes=nodes) # run the function S.PhyloMaker.

str(result) # the structure of the ouput of S.PhyloMaker.

phytools::plotTree(result$Scenario.3,fsize=1.2,ftype="i")




library(ape)
library(phytools)

# Assuming 'result_tree' is the tree object obtained from running S.PhyloMaker
result_tree <- result$Scenario.3  # For example, using scenario 3

# Calculate the pairwise phylogenetic distances between all species
dist_matrix <- cophenetic.phylo(result_tree)

# Function to find the top two closest relatives for each species
find_closest_relatives <- function(dist_matrix) {
  closest_relatives <- apply(dist_matrix, 1, function(row) {
    # Order distances, take the second and third since the first is the distance to itself (which is 0)
    ordered <- order(row, decreasing = FALSE)
    names(row)[ordered[2:3]]
  })
  
  return(closest_relatives)
}

# Get the closest relatives for all species
closest_relatives <- find_closest_relatives(dist_matrix)

list(closest_relatives)




















#### Misc ####
vanilla=subset_samples(raw_ps, RowName=="WMAVan.16S.S101")
ps.van=vanilla

#ps.van=tax_glom(vanilla, taxrank=rank_names(ps.ra)[4], NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
ps.van.select=prune_taxa(taxa_sums(ps.van)>100, ps.van)

methods_bar_1=plot_bar(ps.van.select, x="seq_id", fill="ta6")+
  geom_bar(stat="identity")+
  scale_fill_viridis_d(name="Phylum")+
  ylab("% Abundance")+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank())
methods_bar_1


# Create the scatter plot
p <- ggplot(ddpcr, aes(x=log10(CH4_int+1), y=log10(ch4_200+1))) +
  geom_point() +  # Add points
  labs(x="log10(ch4_int + 1)", y="log10(ch4_200)") +  # Label axes
  theme_minimal()+
  stat_cor()+
  ylim(-0.1,0.3)+
  facet_grid(.~core_type)

# Display the plot
print(p)


