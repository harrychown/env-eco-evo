library(ape)
library(dplyr)
library(tidyr)
library(pez)
library(geiger)
library(lme4)
library(mvtnorm)
library(rstan)
library(bayesplot)
library(ggplot2)
library(patchwork)

# For each group, find MRCA and rename clade
collapseTree <- function(input_tree, key){
  output_tree <- input_tree
  # Retain only keys that are in current tree
  key <- key[which(key[,1] %in% input_tree$tip.label),]
  for (grp in unique(key[,2])) {
    tips <- key[,1][key[,2] == grp]
    mrca_node <- getMRCA(output_tree, tips)
    subtree_tips <- extract.clade(output_tree, mrca_node)$tip.label
    
    # Drop all but one representative (first tip)
    output_tree <- drop.tip(output_tree, setdiff(subtree_tips, tips[1]))
    
    # Rename the remaining tip
    output_tree$tip.label[output_tree$tip.label == tips[1]] <- grp
  }
  # Drop anything that hasn't been collapsed
  output_tree <- drop.tip(output_tree, output_tree$tip.label[-match(unique(key[,2]), output_tree$tip.label)])
  return(output_tree)
}


generateInputData <- function(input_data, species_var, env_var, site_var = "site_name"){
  colnames(input_data)[which(colnames(input_data) == species_var)] <- "species"
  trait_per_species <- input_data[,c(site_var, "species", trait_var)]
  # Calculate trait per species
  trait_per_species <- trait_per_species %>%  
    group_by(species) %>%
    filter(!is.na(species)) %>% 
    summarise(
      trait = mean(.data[[trait_var]])
    )
  trait <- trait_per_species$trait
  names(trait) <- trait_per_species$species
  # Generate binary presence/absence information
  tmp <- input_data %>% mutate(presence = 1) %>%
    complete(.data[[site_var]], species, fill = list(presence = 0))
  presence_info <- unique(tmp[,c(site_var,"species","presence")]) %>% 
    left_join(trait_per_species, by = c("species")) 
  colnames(presence_info)[1] <- "site"
  # Add numerical site/species information
  num_species <- 1:length(unique(presence_info$species))
  names(num_species) <- unique(presence_info$species)
  presence_info$num_species <- num_species[presence_info$species]
  num_site <- 1:length(unique(presence_info$site))
  names(num_site) <- unique(presence_info$site)
  presence_info$num_site <- num_site[presence_info$site]
  env_info <- input_data %>% 
    distinct(.data[[site_var]], .data[[env_var]])
  colnames(env_info) <- c("site", "env")
  traits_env_occ <- presence_info %>% 
    left_join(env_info, by=c("site")) %>%
    filter(!is.na(species))
  return(list("trait" = trait,
              "occ" = traits_env_occ))
}

plotLatent <- function(model_output, trait_input, type = "basic" ){
  if(type == "basic"){
    effects <- ranef(model_output)
    latent.traits <- setNames(effects$species$env,  rownames(effects$species))
    latent.traits <- latent.traits[names(trait_input)]
  }
  else{
    latent.traits <- apply(rstan::extract(model_output)[[type]], 2, mean)
    names(latent.traits) <- names(trait_input)
    latent.traits <- latent.traits[names(trait_input)]
  }
  c <- cor.test(latent.traits, trait_input)
  return(list("latent" = latent.traits,
              "actual" = trait_input,
              "c" = c))
}

assessPerformance <- function(model_output, y,  model_params, env_name){
  # Generate a basic summary of parameters
  basic_summary <- as.data.frame(summary(model_output, pars=model_params)$summary)
  basic_summary <- cbind(rep(env_name, nrow(basic_summary)), rownames(basic_summary), basic_summary)
  colnames(basic_summary)[1:2] <- c("env_var", "posterior")
  # Predictive power
  yrep <- rstan::extract(model_output, pars = "yrep")$yrep
  posterior_checks_barplot <- ppc_bars(y, yrep[1:100, ])
  posterior_checks_stat <- ppc_stat(y, yrep, stat = "mean")
  ppcMean <- c(mean(y), mean(yrep))
  names(ppcMean) <- c("y", "yrep")
  # Generate trace plots
  posterior <- as.array(model_output)
  trace_plots <- mcmc_trace(posterior, pars = model_params)
  return(list("basicSummary" = basic_summary,
              "ppcBar" = posterior_checks_barplot,
              "ppcStat" = posterior_checks_stat,
              "trace" = trace_plots,
              "ppcMean" = ppcMean))
}



runBasicGLMM <- function(input_data){
  basic_GLMM <- glmer(presence ~ (env|species), family=binomial, data=input_data$occ)
  basic_GLMM.latent <- plotLatent(basic_GLMM, input_data$trait, type = "basic")
  return(list("model" = basic_GLMM,
              "latentTraits" = basic_GLMM.latent))
}


runBayes <- function(input_data, env_var){
  stan_code <- "
// Tell rstan what our data are
data{
int Ntotal;              // Number of rows in our dataset
int Nspp;                // Number of species
int presence[Ntotal];    // The response variable (presence/absence)
real env[Ntotal];        // An explanatory variable (environment)
int spp[Ntotal];         // An explanatory variable (species)
}
// The parameters/coefficients we want to estimate and report back
parameters{
vector[Nspp] spp_int;  // Each species' overall occupancy
vector[Nspp] spp_slp;      // Each species' environmental response
}
// Some calculations rstan is going to perform internally to help model-fitting
transformed parameters{
vector[Ntotal] predictions;   // Make a variable to hold our model predictions
// Loop over all our input data and specify our model, which is:
// a species' intercept (overall occupancy) +  env response x the environment
for (i in 1:Ntotal)
predictions[i] = spp_int[spp[i]] + spp_slp[spp[i]]*env[i];
}
// Fit our model to our predictions
model{
// Species' occupancies and responses are drawn from uninformative priors
spp_int ~ normal(0, 0.5);
spp_slp ~ normal(0, 0.5);
// The model itself: our presences are drawn from our predictions
presence ~ bernoulli_logit(predictions);
}
// Extract predictions
generated quantities {
int yrep[Ntotal];
for (n in 1:Ntotal)
  yrep[n] = bernoulli_logit_rng(predictions[n]);
}
"
bayesModel <- stan(model_code=stan_code,
                           data=list(Ntotal=nrow(input_data$occ), Nspp=length(unique(input_data$occ$num_species)),
                                     presence=input_data$occ$presence,
                                     env=input_data$occ$env, spp=input_data$occ$num_species),
                           iter=1000, chains=4, seed=123456, cores = 8,
                   control = list(adapt_delta = 0.99))
# Update chains to 4
# Update iterations to 1000

bayesModel.latent <- plotLatent(bayesModel, input_data$trait, type = "spp_slp")
bayesModel.performance <- assessPerformance(pglmm, input_data$occ$presence,  c("mean_int", "mean_slp"), env_var)

return(list("model" = bayesModel,
            "latentTraits" = bayesModel.latent,
            "performance"=bayesModel.performance))
}

runPGLMM <- function(input_data, input_tree, env_var){
  stan_code <- "
// Function to transform a phylogenetic VCV according to Pagel's Lambda
// - and multiply through by overall sigma
functions {
matrix lambda_vcv(matrix vcv, real lambda, real sigma){
matrix[rows(vcv),cols(vcv)] local_vcv; // Make a local copy of the VCV to modify
local_vcv = vcv * lambda;              // Lambda transforms are just a multiplier
for(i in 1:rows(local_vcv))            // ... but we do have to loop over the matrix
local_vcv[i,i] = vcv[i,i];           // ... and make the diagonal the same as it was before
return(local_vcv * sigma);             // Return the transformed matrix x sigma (overall variation)
}
}
data {
int Ntotal;
int Nspp;
int presence[Ntotal];
real env[Ntotal];
int spp[Ntotal];
//
matrix[Nspp,Nspp]Vphy;     // Give the phylogeny as data
}
parameters{
vector[Nspp] spp_int;
vector[Nspp] spp_slp;
// Lambda transforms for the intercepts and slopes
real<lower=0> lam_int;     // (priors --> cannot be negative)
real<lower=0> lam_slp;
// Coefficients for the NON-phylogenetically-derived variance in model terms
real<lower=0> null_int;    // (priors --> cannot be negative)
real<lower=0> null_slp;
// Coefficients for the mean intercepts/slopes
real mean_int;
real mean_slp;
}
transformed parameters{
vector[Ntotal] predictions;
for (i in 1:Ntotal)
predictions[i] = spp_int[spp[i]] + spp_slp[spp[i]]*env[i];
}
model{
// Specify priors
mean_int ~ normal(0,0.5);
mean_slp ~ normal(0,0.5);
lam_int ~ normal(0,0.5);
lam_slp ~ normal(0,0.5);
null_int ~ normal(0,0.5);
null_slp ~ normal(0,0.5);

// Now we draw our species coefficients, incorporating the lambda-transformed phylogeny
spp_int ~ multi_normal(rep_vector(mean_int,Nspp), lambda_vcv(Vphy,lam_int,null_int));
spp_slp ~ multi_normal(rep_vector(mean_slp,Nspp), lambda_vcv(Vphy,lam_slp,null_slp));
//
presence ~ bernoulli_logit(predictions);
}
// Extract predictions
generated quantities {
int yrep[Ntotal];
for (n in 1:Ntotal)
  yrep[n] = bernoulli_logit_rng(predictions[n]);
}
"
pglmm <- stan(model_code=stan_code,
                             data=list(Ntotal=nrow(input_data$occ),  Nspp=length(unique(input_data$occ$num_species)),
                                       presence=input_data$occ$presence,
                                       env=input_data$occ$env, spp=input_data$occ$num_species, Vphy=vcv(input_tree)),
              iter=1000, chains=4, seed=123456, cores = 8,
              control = list(adapt_delta = 0.99))
# Update chains to 4
# Update iterations to 1000
# Assess model performance
pglmm.latent <- plotLatent(pglmm, input_data$trait, type = "spp_slp")
pglmm.performance <- assessPerformance(pglmm, input_data$occ$presence,  c("mean_int", "mean_slp", "lam_int", "lam_slp"), env_var)
return(list("model" = pglmm,
            "latentTraits" = pglmm.latent,
            "performance" = pglmm.performance))
}

#--------------------------------------------
# 1. Set inputs and parameters
#--------------------------------------------

tree_filename <- "/home/harry/Documents/imperial/pglmm/nerc/data/NERC_w_outgroup-021025-bs100.raxml.support"
traits_filename <- "/home/harry/Documents/imperial/pglmm/nerc/data/trait_information.csv"
env_filename <- "/home/harry/Documents/imperial/nerc_year1_summary/data/metadata_w_climate.antifungal-filled.csv"
group_filename <- "/home/harry/Documents/imperial/pglmm/nerc/data/tip_group_key.csv"
subgroup_filename <- "/home/harry/Documents/imperial/pglmm/nerc/data/tip_subgroup_key.csv"
output_directory <- "/home/harry/Documents/imperial/pglmm/nerc/results"
trait_var <- "resistance"


#--------------------------------------------
# 2. Load and clean phylogeny
#--------------------------------------------
raw_tree <- read.tree(tree_filename)
raw_tree <- root(raw_tree, outgroup = "C250.final_snps_WGS", resolve.root = TRUE)

id_to_drop <- c("C250", "Sample_78-", "Sample_115-", "Sample_75-", "Sample_11-", "Sample_141-")
tip_to_drop <- raw_tree$tip.label[grepl(paste(id_to_drop, collapse = "|"), raw_tree$tip.label)]
tree <- drop.tip(raw_tree, tip_to_drop)

#--------------------------------------------
# 3. Load and merge trait/environmental data
#--------------------------------------------
raw_traits <- read.csv(traits_filename)
raw_env <- read.csv(env_filename)
env <- raw_env[raw_env$sample_type == "air", c(2,4,10:ncol(raw_env))]
traits_env <- unique(merge(raw_traits, env, by = "site_name"))

# Map sequences to tree tip labels
seq_pattern <- paste0(traits_env$seq_id, ".final_snps")
traits_env$seq_name <- sapply(seq_pattern, function(p) {
  match <- tree$tip.label[grepl(p, tree$tip.label, fixed = TRUE)]
  if (length(match) == 0) "0" else match[1]
})

traits_env <- traits_env[traits_env$seq_name != "0",]
traits_env <- traits_env[!is.na(traits_env$landcover),]
# Update tree to match known traits
tree <- drop.tip(tree, tree$tip.label[-match(unique(traits_env$seq_name), tree$tip.label)])


#--------------------------------------------
# 4. Load and merge group/subgroup data
#--------------------------------------------
group_info <- read.csv(group_filename)
subgroup_info <- read.csv(subgroup_filename)
colnames(group_info)[1] <-  "seq_name"
colnames(subgroup_info)[1] <-  "seq_name"

traits_env <- traits_env %>%
  left_join(group_info, by = "seq_name") %>% 
  left_join(subgroup_info, by = "seq_name") 

#--------------------------------------------
# 4. Generate group/subgroup phylogenies
#--------------------------------------------
group_tree <- collapseTree(tree, group_info)
subgroup_tree <- collapseTree(tree, subgroup_info)

env_var_list <- c("Cyproconazole", "Prochloraz", "Pyraclostrobin")

#--------------------------------------------
# 5. Generate model inputs
#--------------------------------------------


for(environmental_variable in env_var_list){
  
  individual_data <-  generateInputData(traits_env, "seq_name", environmental_variable)
  group_data <-  generateInputData(traits_env, "group.x", environmental_variable)
  subgroup_data <-  generateInputData(traits_env, "group.y", environmental_variable)
  
  #--------------------------------------------
  # 5. Run basic GLMM
  #--------------------------------------------
  #basic_individual <- runBasicGLMM(individual_data)
  #basic_group <- runBasicGLMM(group_data)
  #basic_subgroup <- runBasicGLMM(subgroup_data)
  
  #--------------------------------------------
  # 5. Run initial Bayesian model
  #--------------------------------------------
  #bayes_individual <- runBayes(individual_data)
  #bayes_group <- runBayes(group_data)
  #bayes_subgroup <- runBayes(subgroup_data)
  
  #--------------------------------------------
  # 5. Run PGLMM
  #--------------------------------------------
  pglmm_individual <- runPGLMM(individual_data, tree, environmental_variable)
  pglmm_group <- runPGLMM(group_data, group_tree, environmental_variable)
  pglmm_subgroup <- runPGLMM(subgroup_data, subgroup_tree, environmental_variable)
  
  #--------------------------------------------
  # 6. Generate overview across grouped models
  #--------------------------------------------
  plot_types <- c("ppcBar", "ppcStat", "trace")
  for(i in plot_types){
    p1 <- pglmm_individual$performance[[i]] + ggtitle("Individual")
    p2 <- pglmm_group$performance[[i]] + ggtitle("Group")
    p3 <- pglmm_subgroup$performance[[i]] + ggtitle("Subgroup")
    p_all <- p1 + p2 + p3
    outplot_name <- paste(c(output_directory, "/", i, ".", environmental_variable, ".jpeg"), sep="", collapse="")
    ggsave(file=outplot_name, plot=p_all, width=45, height=30, units="cm", bg="white")
  }
  
  pglmm_combined <- rbind(pglmm_individual$performance$basicSummary, 
                          pglmm_group$performance$basicSummary, 
                          pglmm_subgroup$performance$basicSummary)
  sub_division <- c()
  for(i in c("individual", "group", "subgroup")){
    sub_division <- c(sub_division, rep(i, nrow(pglmm_combined) / 3))
  }
  pglmm_combined$sub_division <- sub_division
  outtab_name <- paste(c(output_directory, "/",environmental_variable, ".csv"), sep="", collapse="")
  write.csv(pglmm_combined, outtab_name)
  
  # Save modelling results
  individual_out <- paste(c(output_directory, "/individual.", environmental_variable, ".rds"), sep="", collapse="")
  saveRDS(pglmm_individual, file = individual_out)
  
  group_out <- paste(c(output_directory, "/group.", environmental_variable, ".rds"), sep="", collapse="")
  saveRDS(pglmm_group, file = group_out)
  
  subgroup_out <- paste(c(output_directory, "/subgroup.", environmental_variable, ".rds"), sep="", collapse="")
  saveRDS(pglmm_subgroup, file = subgroup_out)
  
  }
