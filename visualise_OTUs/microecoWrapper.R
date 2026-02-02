library(file2meco)
library(microeco)
library(magrittr)
library(ggplot2)
library(paletteer)
library(ggradar)
library(retroPal)
library(scales) # Log-transform
library(optparse)
#--------------------------------------------
# 1. Set input data
#--------------------------------------------
option_list <- list(
  make_option(c("-a", "--abundance"), type = "character",
              help = "Abundance file in MPA format"),
  make_option(c("-s", "--sample"), type = "character",
              help = "Sample file (CSV file: row names equal sample ID's in abundance file, followed by categorical information"),
  make_option(c("--categorical"), type = "character",
              help = "Column name of categorical variable for splitting"),
  make_option(c("--taxon"), type = "character",
              default = "Genus",
              help = "Taxon level for analysis. Default: Genus\n
              Must choose one of Domain, Kingdom, Phylum, Class, order, Family, Genus or Species"),
  make_option(c("--cat-order"), type = "character",dest = "cat_order",
              help = "Comma separated level order of categorical variables [e.g. Indoor,Outdoor]"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output directory")
)


taxon_level <- list("Domain" = "d__",
                    "Kingdom"="k__",
                    "Phylum"="p__",
                    "Class"="c__",
                    "Order"="o__",
                    "Family"="f__",
                    "Genus"="g__",
                    "Species"="s__"
                    )

opt <- parse_args(OptionParser(option_list = option_list))

abund_file_path <- opt$abundance
sample_file_path <- opt$sample
categorical_variable <- opt$categorical
output_folder <- opt$output
taxon <- opt$taxon
#match_file_path <- system.file("extdata", "example_metagenome_match_table.tsv", package="file2meco")
data <- mpa2meco(abund_file_path, sample_table = sample_file_path, rel = TRUE, use_level = taxon_level[[taxon]])

# Set order
if (is.null(opt$`cat_order`)) {
  cat_order <- unique(data$sample_table[[categorical_variable]])
} else {
  cat_order <- trimws(unlist(strsplit(opt$cat_order, ",")))
}

data$sample_table[[categorical_variable]] %<>% factor(., levels = cat_order)
# Set colourscheme
cat_cols <- get_retro_pal("mid4")[1:length(cat_order)]
names(cat_cols) <- cat_order
cat_outline <- get_retro_pal("dark5")[1:length(cat_order)]
names(cat_outline) <- cat_order

#--------------------------------------------
# 2. Display abundances
#--------------------------------------------
# Standard barplot
abundance <- trans_abund$new(dataset = data, taxrank = taxon, ntaxa = 20)
abundance_plot <- abundance$plot_bar(others_color = "grey70", facet = categorical_variable, xtext_keep = FALSE, legend_text_italic = TRUE, color_values = get_retro_pal("categorical")[1:20])
abundance_plot <- abundance_plot + guides(fill = guide_legend(ncol = 1)) + pub_theme(x_axis_labels = TRUE, x_axis_rotation = -90) 
ggsave(file=file.path(output_folder, "relative_abundance.svg"), plot=abundance_plot, width=16, height=16, units="cm", bg="white")

# Grouped-mean barplot
abundance_mean <- trans_abund$new(dataset = data, taxrank = taxon, ntaxa = 20, groupmean = categorical_variable)
grouped_mean<- abundance_mean$plot_bar(others_color = "grey70", xtext_keep = FALSE, legend_text_italic = FALSE, color_values = get_retro_pal("categorical")[1:20])
grouped_mean <- grouped_mean + guides(fill = guide_legend(ncol = 1)) + pub_theme(x_axis_labels=TRUE)
ggsave(file=file.path(output_folder, "relative_abundance.grouped_mean.svg"), plot=grouped_mean, width=16, height=16, units="cm", bg="white")

# Abundance distribution between groups
grouped_distrib_box <- abundance$plot_box(group = categorical_variable, color_values = cat_outline) 
grouped_distrib_box <- grouped_distrib_box + scale_fill_manual(values=cat_cols) + guides(fill=guide_legend(title="Group")) + pub_theme(x_axis_labels=T, x_axis_rotation = 30, x_axis_italic = TRUE)
grouped_distrib_box@layers$geom_segment$aes_params$colour <- grouped_distrib_box@layers$geom_segment$data$fill
ggsave(file=file.path(output_folder, "relative_abundance.grouped_distributions.svg"), plot=grouped_distrib_box, width=16, height=16, units="cm", bg="white")

# Heatmap visualisation
heatmap <- abundance$plot_heatmap(facet = categorical_variable, xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10), legend_title = "Relative\nabundance (%)")
heatmap <- heatmap + 
  pub_theme(y_axis_italic = TRUE, drop_grid = TRUE) +
  scale_fill_gradientn(colours = c(get_retro_pal("brown")), trans = "pseudo_log", na.value = "grey90", breaks = c(0, 1, 2, 5, 10, 25, 50, 75), limits = c(0, max(abundance$data_abund$Abundance)))
ggsave(file=file.path(output_folder, "relative_abundance.heatmap.svg"), plot=heatmap, width=16, height=16, units="cm", bg="white")

# Number of unique/shared taxa per group
merged <- data$merge_samples(categorical_variable)
# tmp is a new microtable object
# create trans_venn object
# Set colour scheme
cat_cols <- get_retro_pal("mid4")[1:length(cat_order)]
names(cat_cols) <- cat_order
cat_outline <- get_retro_pal("dark5")[1:length(cat_order)]
names(cat_outline) <- cat_order
merged_venn_d <- trans_venn$new(merged, ratio = "seqratio")
venn_col_order <- cat_cols[merged_venn_d$res_names]
venn <- merged_venn_d$plot_venn(color_circle = venn_col_order, petal_plot = F, text_name_size = 6 / .pt, text_size = 5 / .pt) +
  ylab("") + xlab("") +
  pub_theme(y_axis_labels = FALSE, drop_grid = TRUE, x_axis_labels = FALSE)
ggsave(file=file.path(output_folder, "shared_taxa_venn.svg"), plot=venn, width=8, height=8, units="cm", bg="white")
#--------------------------------------------
# 3. Calculate alpha diversity
#--------------------------------------------
# Alpha diversity (within samples)
alpha_diversity <- trans_alpha$new(dataset = data, group = categorical_variable)
# Test we use is dependent on no. groups 
alpha_diversity$cal_diff(method = "KW")
#alpha_diversity$res_diff %<>% base::subset(Significance != "ns") # Remove non-significant label
# Plot Shannon diversity
shannon <- alpha_diversity$plot_alpha(measure = "Shannon", y_increase = 0.3, add = "dotplot", color_values = cat_outline, point_size = 1) 
shannon@layers$geom_dotplot$geom_params$dotsize <- 0.25
shannon@layers$geom_boxplot$aes_params$fill <- get_retro_pal("mid3")[1:length(cat_order)]
shannon <- shannon +
  pub_theme(x_axis_labels = TRUE)
ggsave(file=file.path(output_folder, "alpha.shannon.svg"), plot=shannon, width=8, height=8, units="cm", bg="white")
# Plot inverse Simpson
inv_simpson <- alpha_diversity$plot_alpha(measure="InvSimpson", y_increase = 0.3, add = "dotplot", color_values = cat_outline, point_size=1) +
  ylab("Inverse Simpson") +
  pub_theme(x_axis_labels = TRUE)
inv_simpson@layers$geom_dotplot$geom_params$dotsize <- 0.25
inv_simpson@layers$geom_boxplot$aes_params$fill <- get_retro_pal("mid3")[1:length(cat_order)]
ggsave(file=file.path(output_folder, "alpha.inverse_simpson.svg"), plot=inv_simpson, width=8, height=8, units="cm", bg="white")
#--------------------------------------------
# 4. Calculate beta diversity
#--------------------------------------------
# Beta diversity (between samples)
data$cal_betadiv(unifrac = F)
beta_diversity <- trans_beta$new(dataset = data, group = categorical_variable, measure = "bray")
# Generate PCoA
beta_diversity$cal_ordination(method = "PCoA")
beta_pcoa <- beta_diversity$plot_ordination(plot_color = categorical_variable, color_values = cat_cols, plot_shape = categorical_variable, plot_type = c("point", "ellipse"), point_size = 1) +
  guides(colour=guide_legend(title="Group"), shape=guide_legend(title="Group"), fill=guide_legend(title="Group")) +
  pub_theme(x_axis_labels = TRUE)
ggsave(file=file.path(output_folder, "beta.pcoa.svg"), plot=beta_pcoa, width=8, height=8, units="cm", bg="white")
# Generate Bray-Curtis plots
beta_diversity$cal_group_distance(within_group = TRUE)
print("DEBUG: pre-KW beta-div")
beta_diversity$cal_group_distance_diff(method = "KW")
print("DEBUG: post-KW beta-div")
#beta_diversity$res_group_distance_diff %<>% base::subset(Significance != "ns") # Remove non-significant label
bray_curtis <- beta_diversity$plot_group_distance(add = "dotplot",plot_group_order = cat_order, color_values = cat_outline, point_size=1) +
  pub_theme(x_axis_labels = TRUE)

bray_curtis@layers$geom_dotplot$geom_params$dotsize <- 0.25
bray_curtis@layers$geom_boxplot$aes_params$fill <- get_retro_pal("mid3")[1:length(cat_order)]
ggsave(file=file.path(output_folder, "beta.bray_curtis.svg"), plot=bray_curtis, width=8, height=8, units="cm", bg="white")


#--------------------------------------------
# 5. Identify enrichments
#--------------------------------------------
lefse <- tryCatch(
  {
    trans_diff$new(dataset = data, method = "lefse", group = categorical_variable, alpha = 0.01, lefse_subgroup = NULL)  
},
  error = function(e) {
    message("No significant taxa in LEfSe")
    NULL
  }
)

if(!is.null(lefse)){
  lefse_lda <- lefse$plot_diff_bar(use_number = 1:30, width = 0.8, threshold = 4, color_group_map = T, color_values = cat_cols) +
    pub_theme(x_axis_labels = T)
  ggsave(file=file.path(output_folder, "lefse.lda.svg"), plot=lefse_lda, width=8, height=8, units="cm", bg="white")
  
  lefse_ra <- lefse$plot_diff_abund(plot_type = "barerrorbar", coord_flip = T, threshold=4, color_values = cat_cols) +
    guides(colour=guide_legend(title="Group"), shape=guide_legend(title="Group"), fill=guide_legend(title="Group")) +
    pub_theme(x_axis_labels = T)
  ggsave(file=file.path(output_folder, "lefse.threshold.svg"), plot=lefse_ra, width=10, height=8, units="cm", bg="white")
}



