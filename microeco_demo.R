library(file2meco)
library(microeco)
library(magrittr)
library(ggplot2)
library(paletteer)
library(ggradar)
#--------------------------------------------
# 1. Set input data
#--------------------------------------------
abund_file_path <- system.file("extdata", "example_kraken2_merge.txt", package="file2meco")
sample_file_path <- system.file("extdata", "example_metagenome_sample_info.tsv", package="file2meco")
match_file_path <- system.file("extdata", "example_metagenome_match_table.tsv", package="file2meco")
data <- mpa2meco(abund_file_path, sample_table = sample_file_path, match_table = match_file_path, rel = TRUE)

#--------------------------------------------
# 2. Display abundances
#--------------------------------------------
taxa_of_interest <- c("Bradyrhizobium", "Psuedomonas")
# Standard barplot
abundance <- trans_abund$new(dataset = data, taxrank = "Genus", ntaxa = 20)
abundance_plot <- abundance$plot_bar(others_color = "grey70", facet = "Group", xtext_keep = FALSE, legend_text_italic = FALSE, color_values = paletteer_d("ggthemes::Classic_20"))
dev.off()
print(abundance_plot)
# Grouped-mean barplot
abundance_mean <- trans_abund$new(dataset = data, taxrank = "Genus", ntaxa = 20, groupmean = "Group")
grouped_mean<- abundance_mean$plot_bar(others_color = "grey70", xtext_keep = FALSE, legend_text_italic = FALSE, color_values = paletteer_d("ggthemes::Classic_20"))
grouped_mean <- grouped_mean + theme_classic() + theme(axis.title.y = element_text(size = 18))
print(grouped_mean)
# Abundance distribution between groups
grouped_distrib_box <- abundance$plot_box(group = "Group", xtext_angle = 30)
print(grouped_distrib_box)
# Heatmap visualisation
heatmap <- abundance$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
heatmap <- heatmap + theme(axis.text.y = element_text(face = 'italic'))
print(heatmap)
# Radar visualisation
abundance_mean <- trans_abund$new(dataset = data, taxrank = "Genus", ntaxa = 20, groupmean = "Group")
radar <- abundance_mean$plot_radar(values.radar = c("0%", "25%", "50%"), grid.min = 0, grid.mid = 0.25, grid.max = 0.5)
print(radar)
# Number of unique/shared taxa per group
merged <- data$merge_samples("Group")
# tmp is a new microtable object
# create trans_venn object
merged_venn_d <- trans_venn$new(merged, ratio = "seqratio")
venn <- merged_venn_d$plot_venn()
print(venn)

#--------------------------------------------
# 3. Calculate alpha diversity
#--------------------------------------------
# Alpha diversity (within samples)
alpha_diversity <- trans_alpha$new(dataset = data, group = "Group")
# Test we use is dependent on no. groups 
alpha_diversity$cal_diff(method = "t.test")
alpha_diversity$res_diff %<>% base::subset(Significance != "ns") # Remove non-significant label
# Plot Shannon diversity
shannon <- alpha_diversity$plot_alpha(measure = "Shannon", y_increase = 0.3, add = "dotplot")
print(shannon)
# Plot inverse Simpson
inv_simpson <- alpha_diversity$plot_alpha(measure="InvSimpson", y_increase = 0.3, add = "dotplot") +
  ylab("Inverse Simpson")
print(inv_simpson)

#--------------------------------------------
# 4. Calculate beta diversity
#--------------------------------------------
# Beta diversity (between samples)
data$cal_betadiv(unifrac = F)
beta_diversity <- trans_beta$new(dataset = data, group = "Group", measure = "bray")
# Generate PCoA
beta_diversity$cal_ordination(method = "PCoA")
beta_pcoa <- beta_diversity$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
# Generate Bray-Curtis plots
beta_diversity$cal_group_distance(within_group = TRUE)
beta_diversity$cal_group_distance_diff(method = "t.test")
beta_diversity$res_group_distance_diff %<>% base::subset(Significance != "ns") # Remove non-significant label
bray_curtis <- beta_diversity$plot_group_distance(add = "dotplot")
print(bray_curtis)

#--------------------------------------------
# 5. Identify enrichments
#--------------------------------------------
lefse <- trans_diff$new(dataset = data, method = "lefse", group = "Group", alpha = 0.01, lefse_subgroup = NULL)
