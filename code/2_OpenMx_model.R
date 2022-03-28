# Fit a bivariate GCSM model for vocal tract heritability estimates
# on the dataset with all DZ twins included
#
# Copyright (C) 2019-2020, Emily Jennings (adapted from Conor Dolan)
# checked and modified by Dan Dediu (in interaction with Conor Dolan), 2021-2022
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


# Assuming the current directory is this script's directory!

run_start_time <- Sys.time(); # start time: takes about 40 minutes on a Core i9-10900 (10 cores) with 32Gb RAM and abut 30 mins on a Ryzen 3700X (8 cores) with 64GB RAM


## Libraries ####
library(stringr);
library(dplyr);
library(parallel);
library(doParallel);
library(foreach);
library(psych);
library(OpenMx);


# Notation conventions: 
# - cervicalSpine_CS2H_atlas_c2base_distance.1_1 is phenotype.rater1_twin1. All phenotypes follow this notation. So:
# - phenotype.1_1 = phenotype.rater1_twin1
# - phenotype.2_1 = phenotype.rater2_twin1
# - phenotype.1_2 = phenotype.rater1_twin2
# - phenotype.2_2 = phenotype.rater2_twin2


# Logs:
if( !dir.exists("../logs") ) dir.create("../logs", showWarnings=FALSE);


## Various choices:
choice_openmx_npsol            <- FALSE; # should we use the NPSOL optimizer for OpenMX?
choice_do_exploratory_analysis <- FALSE; # should we run the exploratory analysis?
choice_z_score_phenotypes      <- TRUE;  # should we z-score the phenotypes?
choice_include_age2            <- TRUE;  # should we include the quadratic effect of age?
choice_include_age_x_sex       <- FALSE; # should we include the interaction of age and sex?
choice_include_icv_x_sex       <- FALSE; # should we include the interaction of ICV and sex?

## Run sequentially or in parallel? ####
# Create the worker cluster:
run_parallel <- TRUE; # run sequentially or in parallel?
if( run_parallel )
{
  no_cores <- parallel::detectCores(); if( no_cores > 1 ) no_cores <- no_cores-1; # try to use almost all cores on the system but feel free to override!
}

## NPSOL optimiser?
if( choice_openmx_npsol )
{
  # Make sure NSPOL is actually available (may need to install the non-CRAN build as per ?omxGetNPSOL()
  mxOption(NULL, "Default optimizer", "NPSOL"); # use NPSOL as the default optimiser
}


## OpenMX extend fitting function ####
# Try even harder to fit an OpenMX model and possibly compare it to another one (already fitted), and return NULL and warn if this still fails...
fit_openmx_model_and_compare <- function(model1, model2_fit=NULL, # the model to fit and possibly the one to compare to (already fitted)
                                         model1_embedded_in_model2=FALSE, # which model contains which (important for comparison)?
                                         extraTries=20, # the number of tries to send to mxTryHard
                                         intervals=TRUE, # compute intervals?
                                         get_summary=FALSE, # return the summary as well?
                                         silent=TRUE, debug_file=NULL)
{
  ret_val <- list("fit"=NULL, "summary"=NULL, "comparison"=NULL); # the return value
  
  ret_val$fit <- NULL;
  try(ret_val$fit <- OpenMx::mxTryHard(model1, extraTries=extraTries, intervals=intervals), silent=silent); # try to fit the model
  if( is.null(ret_val$fit) || inherits(ret_val$fit, "try-error") ) next; # fitting seems to have failed: try again if possible!
  
  # Apparently fitting succeeded:
  if( get_summary ) ret_val$summary <- summary(ret_val$fit, verbose=TRUE); # return the summary
  
  # Model comparison (if the case):
  if( !is.null(model2_fit) )
  {
    ret_val$comparison <- NULL;
    if( model1_embedded_in_model2 ) # the order of model comparison
    {
      try(ret_val$comparison <- OpenMx::mxCompare(model2_fit, ret_val$fit), silent=silent); # try to do model comparison
    } else
    {
      try(ret_val$comparison <- OpenMx::mxCompare(ret_val$fit, model2_fit), silent=silent); # try to do model comparison
    }
  }
  
  if( is.null(ret_val$fit) || inherits(ret_val$fit, "try-error") ){ warning(paste0("The model did not coverge:",ret_val$fit,"\n")); return (NULL); }
  if( !is.null(model2_fit) && is.null(ret_val$comparison) ) warning("The model comparison did not coverge!\n");
  return (ret_val);
}


## Load data and adjust variables ####
if( !file.exists("../data/intermediate/data_for_SEM.csv") ) source("./1_data_preparation.R", echo=FALSE); # run the data preparation first...
SEMdata <- read.csv("../data/intermediate/data_for_SEM.csv");
SEMdata$X <- NULL; # rownames are not needed
SEMdata <- SEMdata %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance.1_1":"general_GAPD_allPointSet_matlabProcDistance.2_1"),
                              c("cervicalSpine_CS2H_atlas_c2base_distance.1_2":"general_GAPD_allPointSet_matlabProcDistance.2_2"),
                              c("Sex_1":"ZYGMZDZ"),"FISNumber_1","FISNumber_2");

# z-score all phenotypic measures for the whole dataset:
if( choice_z_score_phenotypes ) # should we z-score the phenotypic measures?
{
  phenotype_names <- unique(stringr::str_sub(names(SEMdata)[grep("\\.[[:digit:]]_[[:digit:]]$", names(SEMdata))], 1, -5)); # get the phenotype names
  for( p in phenotype_names )
  {
    x.names  <- grep(p, names(SEMdata), fixed=TRUE); # phenotype names
    x.values <- SEMdata[,x.names]; # phenotype values concatenated
    
    # z-score both twins for rater 1:
    x.rater <- grep(".1_",names(x.values),fixed=TRUE);
    x.valrater <- as.numeric(unlist(x.values[,x.rater])); # phenotype values concatenated
    x.valrater <- scale(x.valrater, center=TRUE, scale=TRUE); # z-score the values
    x.values[,x.rater] <- as.data.frame(matrix(x.valrater, ncol=length(x.rater), byrow=FALSE)); # put them back
    
    # z-score both twins for rater 2:
    x.rater <- grep(".2_",names(x.values),fixed=TRUE);
    x.valrater <- as.numeric(unlist(x.values[,x.rater])); # phenotype values concatenated
    x.valrater <- scale(x.valrater, center=TRUE, scale=TRUE); # z-score the values
    x.values[,x.rater] <- as.data.frame(matrix(x.valrater, ncol=length(x.rater), byrow=FALSE)); # put them back
    
    # Put the z-scores values back:
    SEMdata[,x.names] <- x.values; # put them back
  }
}

# Read the phenotype info:
phenotype_names <- read.csv("../data/input/phenotype_names_decriptions.csv", stringsAsFactors=FALSE);


# Split into MZ and DZ data:
dmz_notmodel <- subset(SEMdata, ZYGMZDZ=="MZ");
ddz_notmodel <- subset(SEMdata, ZYGMZDZ=="DZ");

# Remove periods from column names so that mxData will work:
names(dmz_notmodel) <- gsub("\\.", "", names(dmz_notmodel));
names(ddz_notmodel) <- gsub("\\.", "", names(ddz_notmodel));



# Select data needed for the model:
dmz <- dmz_notmodel %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance1_1":"general_GAPD_allPointSet_matlabProcDistance1_1"),
                               c("cervicalSpine_CS2H_atlas_c2base_distance2_1":"general_GAPD_allPointSet_matlabProcDistance2_1"), 
                               c("cervicalSpine_CS2H_atlas_c2base_distance1_2":"general_GAPD_allPointSet_matlabProcDistance1_2"),
                               c("cervicalSpine_CS2H_atlas_c2base_distance2_2":"general_GAPD_allPointSet_matlabProcDistance2_2"),
                               "Sex_1","z_Age_at_MRI_1","z_Age_at_MRI_SQR_1","z_ICV_new_1","Sex_2","z_Age_at_MRI_2","z_Age_at_MRI_SQR_2", 
                               "z_ICV_new_2", "FISNumber_1", "FISNumber_2", "FamilyNumber", "mult_id_fam");
ddz <- ddz_notmodel %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance1_1":"general_GAPD_allPointSet_matlabProcDistance1_1"),
                               c("cervicalSpine_CS2H_atlas_c2base_distance2_1":"general_GAPD_allPointSet_matlabProcDistance2_1"), 
                               c("cervicalSpine_CS2H_atlas_c2base_distance1_2":"general_GAPD_allPointSet_matlabProcDistance1_2"),
                               c("cervicalSpine_CS2H_atlas_c2base_distance2_2":"general_GAPD_allPointSet_matlabProcDistance2_2"),
                               "Sex_1","z_Age_at_MRI_1","z_Age_at_MRI_SQR_1","z_ICV_new_1","Sex_2","z_Age_at_MRI_2","z_Age_at_MRI_SQR_2",
                               "z_ICV_new_2", "FISNumber_1", "FISNumber_2", "FamilyNumber", "mult_id_fam");

# Recode sex numerically (male=0, female=1):
dmz$Sex_1 <- ifelse(dmz$Sex_1 == "male", 0, 1);
ddz$Sex_1 <- ifelse(ddz$Sex_1 == "male", 0, 1);
dmz$Sex_2 <- ifelse(dmz$Sex_2 == "male", 0, 1);
ddz$Sex_2 <- ifelse(ddz$Sex_2 == "male", 0, 1);

# Rename columns (use shorter names):
colnames(dmz)[colnames(dmz)=="Sex_1"]              <- "sex1";
colnames(ddz)[colnames(ddz)=="Sex_1"]              <- "sex1";
colnames(dmz)[colnames(dmz)=="Sex_2"]              <- "sex2";
colnames(ddz)[colnames(ddz)=="Sex_2"]              <- "sex2";
colnames(dmz)[colnames(dmz)=="z_Age_at_MRI_1"]     <- "age1";
colnames(ddz)[colnames(ddz)=="z_Age_at_MRI_1"]     <- "age1";
colnames(dmz)[colnames(dmz)=="z_Age_at_MRI_2"]     <- "age2";
colnames(ddz)[colnames(ddz)=="z_Age_at_MRI_2"]     <- "age2";
colnames(dmz)[colnames(dmz)=="z_Age_at_MRI_SQR_1"] <- "satm1";
colnames(ddz)[colnames(ddz)=="z_Age_at_MRI_SQR_1"] <- "satm1";
colnames(dmz)[colnames(dmz)=="z_Age_at_MRI_SQR_2"] <- "satm2";
colnames(ddz)[colnames(ddz)=="z_Age_at_MRI_SQR_2"] <- "satm2";
colnames(dmz)[colnames(dmz)=="z_ICV_new_1"]        <- "icv1";
colnames(ddz)[colnames(ddz)=="z_ICV_new_1"]        <- "icv1";
colnames(dmz)[colnames(dmz)=="z_ICV_new_2"]        <- "icv2";
colnames(ddz)[colnames(ddz)=="z_ICV_new_2"]        <- "icv2";

# Replace NAs with 0s (this makes sure that only the twin with the missing data is ignored but not the other twin that may contain useful info)
# Full explanation:
# Let's assume a case with variables X, Y and Z, and covariates Age and Sex. 
# In OpenMx the covariates are treated as fixed regressors, which are not allowed to be missing, because the model is fitted conditional on the covariates.
# Here is a complete case -> no issues:
#   X  Y Z age sex
#   2 5  2  19   1
# Here is a case with missing phenotype Z -> no problem, the stimation works fine:
#   X Y  Z age sex
#   2 5 NA  19   1
# Missing Sex -> this is a problem, as OpenMx will not accept this type of missingness (if you do not know the sex, you lose the case!):
#   X Y  Z age sex
#   2 5  2  19  NA
# More precisely, here is what is going on (1=twin 1 ; 2 = twin 2); suppose
#   X1 Y1 Z1 age1 sex1   X2 Y2 Z2 age2 sex2
#    3  4  5   19    0   NA NA NA   NA   NA
# because of NA age2 and NA sex2, we would lose this twin pair, but that would be a waste as we do have the twin 1 data.
# Solution: recode age2 and sex2 as have any value (e.g., 0):
#   X1 Y1 Z1 age1 sex1   X2 Y2 Z2 age2 sex2
#    3  4  5   19    0   NA NA NA    0    0
# With this, we can retain the case, and the "made-up" values 0 and 0 have no effect, because the twin 2 phenotypes are missing anyway, and will be ignored.
# However, it is important to check that all the phenotypes are missing (X2 Y2 Z2).
dmz$sex2[is.na(dmz$sex2)]   <- 0;
ddz$sex2[is.na(ddz$sex2)]   <- 0;
dmz$age2[is.na(dmz$age2)]   <- 0;
ddz$age2[is.na(ddz$age2)]   <- 0;
dmz$satm2[is.na(dmz$satm2)] <- 0;
ddz$satm2[is.na(ddz$satm2)] <- 0;
dmz$icv2[is.na(dmz$icv2)]   <- 0;
ddz$icv2[is.na(ddz$icv2)]   <- 0;


# The names and number of phenotypic measures:
varnames <- colnames(dmz %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance1_1":"general_GAPD_allPointSet_matlabProcDistance1_1")));
n_vars <- length(varnames); stopifnot(n_vars == 150); # there should be 150 of them!


## MZ and DZ correlations ####
rmz <- foreach(i=1:n_vars) %do% cor(dmz[,c(i, i+n_vars, i+2*n_vars, i+3*n_vars)],use = "complete.obs");
names(rmz) <- colnames(dmz %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance1_1":"general_GAPD_allPointSet_matlabProcDistance1_1")));
rdz <- foreach(i=1:n_vars) %do% cor(ddz[,c(i, i+n_vars, i+2*n_vars, i+3*n_vars)],use = "complete.obs");
names(rdz) <- colnames(ddz %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance1_1":"general_GAPD_allPointSet_matlabProcDistance1_1")));

## MZ and DZ covariances ####
smz <- foreach(i=1:n_vars) %do% cov(dmz[,c(i, i+n_vars, i+2*n_vars, i+3*n_vars)],use = "complete.obs");
names(smz) <- colnames(dmz %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance1_1":"general_GAPD_allPointSet_matlabProcDistance1_1")));
sdz <- foreach(i=1:n_vars) %do% cov(ddz[,c(i, i+n_vars, i+2*n_vars, i+3*n_vars)],use = "complete.obs");
names(sdz) <- colnames(ddz %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance1_1":"general_GAPD_allPointSet_matlabProcDistance1_1")));

# All the phenotypic measures:
selVars <- colnames(dmz %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance1_1":"general_GAPD_allPointSet_matlabProcDistance1_1"),
                                   c("cervicalSpine_CS2H_atlas_c2base_distance2_1":"general_GAPD_allPointSet_matlabProcDistance2_1"), 
                                   c("cervicalSpine_CS2H_atlas_c2base_distance1_2":"general_GAPD_allPointSet_matlabProcDistance1_2"),
                                   c("cervicalSpine_CS2H_atlas_c2base_distance2_2":"general_GAPD_allPointSet_matlabProcDistance2_2")))


## Checks using the  MZ and DZ "raw" correlations ####
raw_checks <- lapply(1:n_vars, function(i)
  {
    # Select the MZ and DZ data for phenotype i:
    dmz_tmp <- dmz[,i + n_vars*(0:3)]; ddz_tmp <- ddz[,i + n_vars*(0:3)]; 
    
    # Standard deviations, means and correlations:
    ss_mz     <- corFiml(dmz_tmp, covar=TRUE);
    ss_dz     <- corFiml(ddz_tmp, covar=TRUE);
    stsd_mz   <- sqrt(diag(ss_mz$cov));
    stsd_dz   <- sqrt(diag(ss_dz$cov));
    str_mz    <- ss_mz$cor;
    str_dz    <- ss_dz$cor;
    stmean_mz <- ss_mz$mean;
    stmean_dz <- ss_dz$mean;
    
    # ACE or ADE?
    model_type_1 <- ifelse(str_mz[3,1] > 2*str_dz[3,1], "ADE", "ACE");
    model_type_2 <- ifelse(str_mz[4,2] > 2*str_dz[4,2], "ADE", "ACE");
    
    # Negative correlations?
    neg_cor_mz <- (str_mz[3,1] < 0 || str_mz[4,2] < 0);
    neg_cor_dz <- (str_dz[3,1] < 0 || str_dz[4,2] < 0);
    
    # Correlations and mean-differences between raters:
    draters <- data.frame("rater1"=c(dmz_tmp[,1], dmz_tmp[,3], ddz_tmp[,1], ddz_tmp[,3]),  # rater 1, MZ+DZ
                          "rater2"=c(dmz_tmp[,2], dmz_tmp[,4], ddz_tmp[,2], ddz_tmp[,4])); # rater 2, MZ+DZ
    r_raters <- cor.test(draters$rater1, draters$rater2); rho_raters <- cor.test(draters$rater1, draters$rater2, method="spearman");
    t_raters <- t.test(draters$rater1, draters$rater2, paired=TRUE); u_raters <- wilcox.test(draters$rater1, draters$rater2, paired=TRUE);
    
    # Correlations between twins:
    dtwins_mz <- data.frame("twin1_MZ"=c(dmz_tmp[,1], dmz_tmp[,2]),  # MZ twin 1 for both raters
                            "twin2_MZ"=c(dmz_tmp[,3], dmz_tmp[,4])); # MZ twin 2 for both raters
    dtwins_dz <- data.frame("twin1_DZ"=c(ddz_tmp[,1], ddz_tmp[,2]),  # DZ twin 1 for both raters
                            "twin2_DZ"=c(ddz_tmp[,3], ddz_tmp[,4])); # DZ twin 2 for both raters
    r_mz <- cor.test(dtwins_mz$twin1_MZ, dtwins_mz$twin2_MZ); rho_mz <- cor.test(dtwins_mz$twin1_MZ, dtwins_mz$twin2_MZ, method="spearman");
    r_dz <- cor.test(dtwins_dz$twin1_DZ, dtwins_dz$twin2_DZ); rho_sz <- cor.test(dtwins_dz$twin1_DZ, dtwins_dz$twin2_DZ, method="spearman");
    
    # Return value:
    return (list("df"=data.frame("phenotype_number"=i, "phentype_name"=stringr::str_sub(colnames(dmz[i]), 1, -4),
                       "mz_mean1_1"=stmean_mz[1], "mz_mean2_1"=stmean_mz[2], "mz_mean1_2"=stmean_mz[3], "mz_mean2_2"=stmean_mz[4], 
                       "dz_mean1_1"=stmean_dz[1], "dz_mean2_1"=stmean_dz[2], "dz_mean1_2"=stmean_dz[3], "dz_mean2_2"=stmean_dz[4], 
                       "mz_stds1_1"=stsd_mz[1],   "mz_stds2_1"=stsd_mz[2],   "mz_stds1_2"=stsd_mz[3],   "mz_stds2_2"=stsd_mz[4], 
                       "dz_stds1_1"=stsd_dz[1],   "dz_stds2_1"=stsd_dz[2],   "dz_stds1_2"=stsd_dz[3],   "dz_stds2_2"=stsd_dz[4], 
                       "model_type_rater_1"=model_type_1, "model_type_rater_2"=model_type_2, 
                       "negative_corr_mz"=neg_cor_mz, "negative_corr_dz"=neg_cor_dz,
                       "mean_rater_1"=mean(draters$rater1, na.rm=TRUE), "mean_rater_2"=mean(draters$rater2, na.rm=TRUE), 
                       "sd_rater_1"=sd(draters$rater1, na.rm=TRUE), "sd_rater_2"=sd(draters$rater2, na.rm=TRUE), 
                       "raters_r"=r_raters$estimate, "raters_r_p"=r_raters$p.value,
                       "raters_rho"=rho_raters$estimate, "raters_rho_p"=rho_raters$p.value,
                       "raters_t"=t_raters$statistic, "raters_t_df"=t_raters$parameter, "raters_t_p"=t_raters$p.value, "raters_t_mean_diff"=t_raters$estimate,
                       "raters_u"=u_raters$statistic, "raters_u_p"=u_raters$p.value, 
                       "rMZ"=r_mz$estimate, "rMZ_p"=r_mz$p.value, "rDZ"=r_dz$estimate, "rDZ_p"=r_dz$p.value,
                       row.names=NULL),
                 "stmean_mz"=stmean_mz, "stmean_dz"=stmean_dz,
                 "str_mz"=str_mz, "str_dz"=str_dz));
  });
raw_checks_df <- do.call(rbind, lapply(raw_checks, function(x) x$df));





## Set up the clusters for running in parallel if so requested ####
if( run_parallel )
{
  cl <- parallel::makeCluster(no_cores);
  doParallel::registerDoParallel(cl);
} else
{
  foreach::registerDoSEQ();
}



## How do we model the two raters? ####
# We can imagine 4 main options: 
# RM_imiv. the two raters have the same means and error variances: this basically assumes that they are identical, but has the smallest number of parameters
# RM_dmiv. the two raters have different means but identical error variances: this assumes that the raters have different "true scores" but identical error rates
# RM_imdv. the two raters have identical means but different error variances: this assumes that the raters have the same "true scores" but different error rates
# RM_dmdv. the two raters have different means and error variances: this assumes that the raters may be different in both respects, but has the largest number of parameters
if( choice_do_exploratory_analysis ) # This is part of the exploratory analysis and should not be normally run...
{
  if( !choice_z_score_phenotypes ) # This should be run on the raw measures!
  {
    # NB: please run with choice_z_score_phenotypes = FALSE to see the actual differences in means, but with choice_z_score_phenotypes=TRUE for the analysis!!! 
    # Let's look at the inter-rater correlations and mean differences across the MZ and DZ twins:
    raw_checks_df <- raw_checks_df[ order(raw_checks_df$raters_r, decreasing=FALSE), ];
    raw_checks_df <- raw_checks_df[ order(abs(raw_checks_df$raters_t_mean_diff), decreasing=TRUE), ];
    jpeg(file="../figures/rater_models.jpg", width=4*3, height=2*3, units="in", res=150, quality=80);
    par(mfrow=c(2,4));
    hist(raw_checks_df$raters_r, main="Histogram of Pearson's r", xlab="Inter-rater Pearson's r");
    hist(raw_checks_df$raters_rho, main="Histogram of Spearman's rho", xlab="Inter-rater Spearman's rho");
    hist(raw_checks_df$raters_t_mean_diff, main="Histogram of mean differences", xlab="Inter-rater mean differences");
    plot(raw_checks_df$raters_t_mean_diff, raw_checks_df$raters_r, main="Pearson's r vs mean difference", xlab="Inter-rater mean differences", ylab="Inter-rater Pearson's r");
    plot(raw_checks_df$mean_rater_1, raw_checks_df$mean_rater_2, main="Mean rater 1 vs mean rater 2", xlab="Mean rater 1", ylab="Mean rater 2");
    plot(raw_checks_df$mean_rater_1, raw_checks_df$raters_t_mean_diff, main="Mean rater 1 vs mean differences", xlab="Mean rater 1", ylab="Inter-rater mean differences");
    plot(raw_checks_df$mean_rater_2, raw_checks_df$raters_t_mean_diff, main="Mean rater 2 vs mean differences", xlab="Mean rater 2", ylab="Inter-rater mean differences");
    par(mfrow=c(1,1));
    dev.off();
    # -> so, it seems like the two raters tend to have the same overall mean
    # -> these suggest that the models with the same "true score" might be ok, so RM_dmiv and RM_dmdv can be dropped
    # -> the choice between identical or different error variances is hard to make a priori
  }
  
  if( choice_z_score_phenotypes ) # This should be run on the z-scored measures !
  {
    # So let's compare the models with 1 and 2 rater error variances (RM_imiv vs RM_imdv):
    # Please note that we are using these comparisons for decide on the best model across all phenotypes and it is therefore not double-dipping (e.g., we don't decide for each phenotype individually in a way that affects their specific results)
    # For detailed explanations about these models and the OpenMX code please see the actual fitting of the models below:
    cat("Compare the 1 vs 2 rater error variance models...\n");
    results <- foreach(i=1:n_vars, .packages=c("OpenMx"), .errorhandling="stop") %dopar% 
      {
        progress_filename <- ifelse(run_parallel, paste0("../logs/progress_phenotype_",i,"_of_",n_vars,".txt"), "");
        cat("Fitting ", varnames[i],"(",i,"of",n_vars,")...\n", file=progress_filename, append=FALSE); # print progress
        
        # Note: we use the estimated values from the raw data as starting points...
        if( sum(raw_checks_df$phenotype_number == i, na.rm=TRUE) != 1 ) stop(paste0("Cannot find the raw estimates for phenotype ",i,"!"));
        rawd <- raw_checks_df[ raw_checks_df$phenotype_number == i, ];
        
        # We do this comparison on the factor model with variance decomposition:
        # N.B., we use implicit starting values, as they work pretty well...
        
        # MZ 2x2 covariance matrix:
        Phmz <- mxMatrix(type="Symm", nrow=2, ncol=2,
                         free  =TRUE,
                         values=matrix(c(1.0,   0.6,
                                         0.6,   1.0),  2, 2, byrow=TRUE),
                         labels=matrix(c("v1",  "cmz",
                                         "cmz", "v1"), 2, 2, byrow=TRUE),
                         name="SMZ1");
        # DZ 2x2 covariance matrix:
        Phdz <- mxMatrix( type="Symm", nrow=2, ncol=2,
                          free  =TRUE,
                          values=matrix(c(1.0,   0.3,
                                          0.3,   1.0),  2, 2, byrow=TRUE),
                          labels=matrix(c("v1",  "cdz",
                                          "cdz", "v1"), 2, 2, byrow=TRUE),
                          name="SDZ1");
        
        # Factor loadings (relates observed phenotypes to latent phenotypes; both fixed to 1.0):
        LY   <- mxMatrix(type="Full", nrow=4, ncol=2,
                         free  =matrix(c(FALSE, FALSE,
                                         FALSE, FALSE,
                                         FALSE, FALSE,
                                         FALSE, FALSE), 4, 2, byrow=TRUE),
                         values=matrix(c(1,     0,
                                         1,     0,
                                         0,     1,
                                         0,     1), 4, 2, byrow=TRUE),
                         labels=matrix(c("rr1", NA,
                                         "rr2", NA,
                                         NA,    "rr1",
                                         NA,    "rr2"), 4, 2, byrow=TRUE),
                         name="LY");
        # Error (initialize with the estimated value):
        TE   <- mxMatrix(type="Diag", nrow=4, ncol=4, free=TRUE, labels=c("ve1", "ve2", "ve1", "ve2"), name="TE"); # different error variances for the two raters
        
        # Matrices to calculate A, C and E variance components:
        covMZ <- mxAlgebra(LY %*% SMZ1 %*% t(LY) + TE, name="SMZ");
        covDZ <- mxAlgebra(LY %*% SDZ1 %*% t(LY) + TE, name="SDZ");
        corMZ <- mxAlgebra(cov2cor(SMZ1), name="RMZ");
        corDZ <- mxAlgebra(cov2cor(SDZ1), name="RDZ");
        
        # The regression model:
        # Means (no reason to assume that intercepts vary with zygosity):
        meanMZ <- mxMatrix(type="Full", nrow=1, ncol=2, free=TRUE, 
                           labels=c("mb01", "mb01"),
                           name="MeanMZ");
        meanDZ <- mxMatrix(type="Full", nrow=1, ncol=2, free=TRUE, 
                           labels=c("mb01", "mb01"),
                           name="MeanDZ");
        
        #
        defSex1    <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.sex1"),  name="Sex1");
        defSex2    <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.sex2"),  name="Sex2");
        defAge1    <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age1"),  name="Age1");
        defAge2    <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age2"),  name="Age2");
        defAgeSqr1 <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.satm1"), name="AgeSqr1");
        defAgeSqr2 <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.satm2"), name="AgeSqr2");
        defIcv1    <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.icv1"),  name="ICV1");
        defIcv2    <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.icv2"),  name="ICV2");
        
        # Regression coefficients sex, age, age^2, and icv, and for interaction terms:
        pathBsex    <- mxMatrix(type="Full", nrow=1, ncol=2, free=TRUE, values=0.0, labels=c("b11","b11"), name="bsex");
        pathBage    <- mxMatrix(type="Full", nrow=1, ncol=2, free=TRUE, values=0.0, labels=c("b22","b22"), name="bage");
        pathBsatm   <- mxMatrix(type="Full", nrow=1, ncol=2, free=choice_include_age2, values=0.0, labels=c("b33","b33"), name="bsatm");
        pathBicv    <- mxMatrix(type="Full", nrow=1, ncol=2, free=TRUE, values=0.0, labels=c("b44","b44"), name="bicv");
        pathBagesex <- mxMatrix(type="Full", nrow=1, ncol=2, free=choice_include_age_x_sex, values=0.0, labels=c("b55","b55"), name="bagesex");
        pathBicvsex <- mxMatrix(type="Full", nrow=1, ncol=2, free=choice_include_icv_x_sex, values=0.0, labels=c("b66","b66"), name="bicvsex");
        
        # Mean corrected for covariates:
        correctedMeanMZ <- mxAlgebra(expression=
                                       cbind(
                                         MeanMZ + (Sex1 %*% bsex) + (Age1 %*% bage) + (AgeSqr1 %*% bsatm) + (ICV1 %*% bicv) + (Age1 %*% Sex1 %*% bagesex) + (ICV1 %*% Sex1 %*% bicvsex),
                                         MeanMZ + (Sex2 %*% bsex) + (Age2 %*% bage) + (AgeSqr2 %*% bsatm) + (ICV2 %*% bicv) + (Age2 %*% Sex2 %*% bagesex) + (ICV2 %*% Sex2 %*% bicvsex)),
                                     name="correctedMeanMZ");
        correctedMeanDZ <- mxAlgebra(expression=
                                       cbind(
                                         MeanDZ + (Sex1 %*% bsex) + (Age1 %*% bage) + (AgeSqr1 %*% bsatm) + (ICV1 %*% bicv) + (Age1 %*% Sex1 %*% bagesex) + (ICV1 %*% Sex1 %*% bicvsex),
                                         MeanDZ + (Sex2 %*% bsex) + (Age2 %*% bage) + (AgeSqr2 %*% bsatm) + (ICV2 %*% bicv) + (Age2 %*% Sex2 %*% bagesex) + (ICV2 %*% Sex2 %*% bicvsex)),
                                     name="correctedMeanDZ");
        
        # Data:
        data_mz <- mxData(observed=dmz, type="raw");
        data_dz <- mxData(observed=ddz, type="raw");
        
        # Expectation objects for Multiple Groups:
        exp_mz <- mxExpectationNormal(covariance="SMZ", means="correctedMeanMZ",
                                      dimnames=c(colnames(dmz)[i], colnames(dmz)[i+n_vars], colnames(dmz)[i+2*n_vars], colnames(dmz)[i+3*n_vars]));
        exp_dz <- mxExpectationNormal(covariance="SDZ", means="correctedMeanDZ",
                                      dimnames=c(colnames(ddz)[i], colnames(ddz)[i+n_vars], colnames(ddz)[i+2*n_vars], colnames(ddz)[i+3*n_vars]));
        
        #
        bits <- c(defSex1, defSex2, defAge1, defAge2, defAgeSqr1, defAgeSqr2, defIcv1, defIcv2,
                  pathBsex, pathBage, pathBsatm, pathBicv, pathBagesex, pathBicvsex);
        pars <- list(Phmz, Phdz, LY, TE);
        
        # ML function and model:
        funML <- mxFitFunctionML();
        model_mz <- mxModel(meanMZ, correctedMeanMZ, bits, pars, covMZ, corMZ, data_mz, funML, exp_mz, name="MZ");
        model_dz <- mxModel(meanDZ, correctedMeanDZ, bits, pars, covDZ, corDZ, data_dz, funML, exp_dz, name="DZ");
        
        # Combine Groups:
        multi <- mxFitFunctionMultigroup(c("MZ","DZ"));
        model_fm <- mxModel("Phenotypic", pars, model_mz, model_dz, funML, multi);
        
        # Fit model:
        cat("Fitting Phenotypic model [model_fm]...\n", file=progress_filename, append=TRUE);
        fit <- fit_openmx_model_and_compare(model_fm, get_summary=TRUE);
        fit_fm <- fit$fit; summary_acde <- fit$summary;
        
        # Force the two rater error variances to be equal:
        model_fm_eqvars <- omxSetParameters(fit_fm, free=TRUE, 
                                            labels   =c("ve1",  "ve2"), 
                                            newlabels=c("ve12", "ve12"),
                                            #values=stTE,
                                            name="Phenotypic (same rater error variance)");
        model_fm_eqvars <- omxAssignFirstParameters(model_fm_eqvars); # deal with the "The free parameter 've12' has been assigned multiple starting values!" error...
        cat("Fitting Phenotypic model with single rater error variance [model_fm_eqvars]...\n", file=progress_filename, append=TRUE);
        fit <- fit_openmx_model_and_compare(model_fm_eqvars, fit_fm, model1_embedded_in_model2=TRUE, get_summary=TRUE); 
        fit_fm_eqvars <- fit$fit; summary_acde_eqvars <- fit$summary; cmp_freevars_eqvars <- fit$comparison; # comparison between two rater errors and same rater error
        if( is.null(cmp_freevars_eqvars) ) cmp_freevars_eqvars <- list("AIC"=c(NA,NA), "p"=c(NA,NA)); # allow cde to continue...
        
        # Return results:
        return (data.frame(
          # Phenotype info:
          "phenotype_number"    =i, 
          "phenotype_name"      =stringr::str_sub(varnames[i], end=-4), 
          "phenotype_name_short"=phenotype_names$short.names[ toupper(stringr::str_sub(phenotype_names$long.names.rater1, end=-3)) == toupper(stringr::str_sub(varnames[i], end=-4)) ],
          
          # Compare the 1 vs 2 rater error variances model:
          "AIC_1r"=cmp_freevars_eqvars$AIC[2], "AIC_2r"=cmp_freevars_eqvars$AIC[1], 
          "deltaAIC"=(cmp_freevars_eqvars$AIC[2] - cmp_freevars_eqvars$AIC[1]), 
          "p"=cmp_freevars_eqvars$p[2],
          
          # Model fit diagnostics:
          "statusCode_2r"=summary_acde_eqvars$statusCode, "statusCode_1r"=summary_acde$statusCode
        ));
      }
    results_df <- do.call(rbind, results);
    
    unlink("../logs/progress_phenotype_*.txt"); # clean the logs....
    
    
    # Save results to file:
    write.csv(results_df, "../data/intermediate/SEM_results_2vs1_raters.csv", row.names=FALSE);
    
    
    # Save the explanation of the results:
    cat(paste0('
# This document explains the content of the `SEM_results_2vs1_raters.csv` file.\n
\n
*N.B.* this document uses GitHub-flavoured Markdown, in particular \\<sup\\>x\\</sup\\> for superscript <sup>x</sup> and \\<sub\\>x\\</sub\\> for subscript <sub>x</sub>.
The columns are:\n
\n
- **phenotype\\_name**: the full phenotype name\n
- **phenotype\\_number**: the phenotype number\n
- **phenotype\\_name\\_short**: the phenotype short name\n
- **AIC_1r** = AIC(1 rate)\n
- **AIC_2r** = AIC(2 rates)\n
- **deltaAIC** = AIC(2 rates) - AIC(1 rate): it is negative if 2 rates are better than 1 rate, and positive otherwise\n
- **p**: the *p*-value of the likelihood test of the difference between the two models\n
      ', 
               "Script took ", difftime(Sys.time(), run_start_time, units="mins"), " minutes to run...\n"), 
        file="../data/intermediate/SEM_results_2vs1_raters_explanations.md", append=FALSE);
    
    
    # Decide the best rater model to use:
    cat("N.B. ΔAIC = AIC(2 rates model) - AIC(1 rate model):\n", file="../data/intermediate/SEM_results_2vs1_raters.txt", append=FALSE);
    
    cat("\nFine-grained ΔAIC guidelines (10 and 2 AIC points):\n", file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes clearly supporting 2 rates (ΔAIC < -10) = ", (x <- sum(results_df$deltaAIC < -10, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes suggestively supporting 2 rates (-10 <= ΔAIC < -2) = ", (x <- sum(results_df$deltaAIC >= -10 & results_df$deltaAIC < -2, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes where 2 raters and 1 rate are equivalent (-2 <= ΔAIC <= 2) = ", (x <- sum(results_df$deltaAIC >= -2 & results_df$deltaAIC <= 2, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes suggestively supporting 1 rate (2 < ΔAIC <= 10) = ", (x <- sum(results_df$deltaAIC > 2 & results_df$deltaAIC <= 10, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes clearly supporting 1 rate1 (ΔAIC > 10) = ", (x <- sum(results_df$deltaAIC > 10, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    
    cat("\nBlunt ΔAIC guidelines (2 AIC points):\n", file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes supporting 2 rates (ΔAIC < -2) = ", (x <- sum(results_df$deltaAIC < -2, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes where 2 raters and 1 rate are equivalent (-2 <= ΔAIC <= 2) = ", (x <- sum(results_df$deltaAIC >= -2 & results_df$deltaAIC <= 2, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes supporting 1 rate (2 < ΔAIC) = ", (x <- sum(results_df$deltaAIC > 2, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    
    cat("\nJust negative and positive ΔAIC:\n", file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes with ΔAIC < 0 = ", (x <- sum(results_df$deltaAIC < 0, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes with ΔAIC > 0 = ", (x <- sum(results_df$deltaAIC > 0, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    
    cat("\nLikelihood ratio test:\n", file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes with p < 0.05 = ", (x <- sum(results_df$p < 0.05, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes with p < 0.05 and ΔAIC < 0 = ", (x <- sum(results_df$p < 0.05 & results_df$deltaAIC < 0, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    cat(paste0("Phenotypes with p < 0.05 and ΔAIC > 0 = ", (x <- sum(results_df$p < 0.05 & results_df$deltaAIC > 0, na.rm=TRUE)), " (", round(x/nrow(results_df)*100,1), "%)\n"), file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
    
    cat("\n", file="../data/intermediate/SEM_results_2vs1_raters.txt", append=TRUE);
  }
  
  # Stop the clusters:
  if( run_parallel )
  {
    parallel::stopCluster(cl);
  }
  
  
  # Based on overall AIC score comparisons, it seems the RM_imiv (1 mean and 1 error variance for the two raters) is by far the best model!
  stop("This is not the real analysis, so stop!\n");
}



## Fit the genetic SEM models ####
cat(">>> Fit the genetic SEM models...\n");
results <- foreach(i=1:n_vars, .packages=c("OpenMx"), .errorhandling="stop") %dopar% 
  {
    progress_filename <- ifelse(run_parallel, paste0("../logs/progress_phenotype_",i,"_of_",n_vars,".txt"), "");
    cat("Fitting ", varnames[i],"(",i,"of",n_vars,")...\n", file=progress_filename, append=FALSE); # print progress
    
    # Note: we use the estimated values from the raw data as starting points...
    if( sum(raw_checks_df$phenotype_number == i, na.rm=TRUE) != 1 ) stop(paste0("Cannot find the raw estimates for phenotype ",i,"!"));
    rawd <- raw_checks_df[ raw_checks_df$phenotype_number == i, ];
    
    
    ## Phenotypic model ####
    cat(">>> Phenotypic models:\n\n", file=progress_filename, append=TRUE);
    comments <- "OK";
    
    
    ##
    ## Explanations (refer to the code below) ##
    ##
    ## We define the MZ and DZ covariance matrices. 
    ## We want to estimate these as unstructured covariance matrices and to obtain ML (Maximum Likelihood) estimates of these summary statistics. 
    ## These covariance matrices are 4x4 and are parameterized by estimating the standard deviations (`sm11`, `sm12`, `sm21`, `sm22`, `sd11`, `sd12`, `sd21`, and `sd22`) 
    ## and the correlation matrices (`Rmz` and `Rdz`). 
    ## The actual covariance matrices are calculated as `Sdmz %*% Rmz %*% t(Sddz)` and `Sddz %*% Rdz %*% t(Sddz)`. 
    ## There are no constraints here, we do not assume that the MZ sd's equal the DZ sd's or that they are equal between twin 1 and twin 2 members of a twin pair. 
    ## The value `mb01` is intercept of the regression of the phenotypes 1 and 2 on the fixed regressors. 
    ## Note that this intercept is the same over twin members and zygosity. 
    ## The intercepts are called `MeanMZ` and `MeanDZ`, respectively.
    ##
    ## We define a phenotypic rater model. 
    ## The associated path diagram is shown in figure SEM_path_01.jpg in the same folder as this script.
    ## In the syntax, `rr1`, `rr2` correspond to the loadings of the phenotypes Ph11 and Ph21 on the latent phenotype f1 (see path diagram). 
    ## In the code these parameters are fixed to 1 (see the definition of `LY`).
    ## In addition we define the covariance matrices of f1 and f2 in the MZ and the DZ twins. 
    ## The MZ and DZ variances are labeled `v1`, so they are constrained to be equal. 
    ## The MZ and DZ covariance (cov(f1,f2) in the path model) are labeled differently, as these are certainly not expected to be equal, hence: `cmz` and `cdz`.
    ## In addition we define the residual variances, i.e., the variances of the phenotypes Ph11 and Ph21 (and of course Ph12 and Ph22) that are not explained by the latent phenotype f1 ( and f2). 
    ## These variances may be interpreted as measurement (rater) error. 
    ## They are labeled `ve12` and are assumed equal over twin members and raters. 
    ##
    
    # MZ 2x2 covariance matrix:
    Phmz <- mxMatrix(type="Symm", nrow=2, ncol=2,
                     free  =TRUE,
                     values=matrix(c(1.0,   0.6,
                                     0.6,   1.0),  2, 2, byrow=TRUE),
                     labels=matrix(c("v1",  "cmz",
                                     "cmz", "v1"), 2, 2, byrow=TRUE),
                     name="SMZ1");
    # DZ 2x2 covariance matrix:
    Phdz <- mxMatrix( type="Symm", nrow=2, ncol=2,
                      free  =TRUE,
                      values=matrix(c(1.0,   0.3,
                                      0.3,   1.0),  2, 2, byrow=TRUE),
                      labels=matrix(c("v1",  "cdz",
                                      "cdz", "v1"), 2, 2, byrow=TRUE),
                      name="SDZ1");
    
    # Factor loadings (relates observed phenotypes to latent phenotypes; both fixed to 1.0):
    LY   <- mxMatrix(type="Full", nrow=4, ncol=2,
                     free  =matrix(c(FALSE, FALSE,
                                     FALSE, FALSE,
                                     FALSE, FALSE,
                                     FALSE, FALSE), 4, 2, byrow=TRUE),
                     values=matrix(c(1,     0,
                                     1,     0,
                                     0,     1,
                                     0,     1), 4, 2, byrow=TRUE),
                     labels=matrix(c("rr1", NA,
                                     "rr2", NA,
                                     NA,    "rr1",
                                     NA,    "rr2"), 4, 2, byrow=TRUE),
                     name="LY");
    # Error (use default init values, as they seem to work well):
    TE   <- mxMatrix(type="Diag", nrow=4, ncol=4, free=TRUE, labels=c("ve12", "ve12", "ve12", "ve12"), name="TE"); # same error variance for the two raters
    
    
    ##
    ## Explanations (refer to the code below) ##
    ##
    ## We compute the phenotypic covariance matrices. 
    ## For example: 
    ##    LY %*% SMZ1 %*% t(LY) + TE
    ## which gives us the parameterization of the MZ covariance matrices, based on the present model (see the path diagram SEM_path_01.jpg). 
    ## If using `mxPath()` or `Mplus`, one would specify this model using indvidual equations; e.g., in `Mplus` it would be something like (for the MZ group): 
    ##    f1 by Ph11* (rr1);
    ##    f1 by Ph21* (rr2);
    ##    f2 by Ph12* (rr1);
    ##    f2 by Ph22* (rr2);
    ##    f1 with f2* (cmz);
    ##    f1* (v1);
    ##    f2* (v1);
    ##    Ph11* (ve1);
    ##    Ph21* (ve2);
    ##    Ph12* (ve1);
    ##    Ph22* (ve2);
    ## `Mplus` forms matrices based on this specification and ultimately uses these model matrices to calculate the expected covariance matrices, as we do in the code below. 
    ##
    
    # Matrices to calculate A, C and E variance components:
    covMZ <- mxAlgebra(LY %*% SMZ1 %*% t(LY) + TE, name="SMZ");
    covDZ <- mxAlgebra(LY %*% SDZ1 %*% t(LY) + TE, name="SDZ");
    corMZ <- mxAlgebra(cov2cor(SMZ1), name="RMZ");
    corDZ <- mxAlgebra(cov2cor(SDZ1), name="RDZ");
    
    
    ##
    ## Explanations (refer to the code below) ##
    ##
    ## We introduce the fixed regressors; this code is the same as for the previous model. 
    ##
    
    # The regression model:
    # Means (no reason to assume that intercepts vary with zygosity):
    means_init_vals <- mean(c(raw_checks[[i]]$stmean_mz, raw_checks[[i]]$stmean_dz)); # average of both MZ and DZ twins across both raters
    meanMZ <- mxMatrix(type="Full", nrow=1, ncol=2, free=TRUE, 
                       values=means_init_vals,
                       labels=c("mb01", "mb01"),
                       name="MeanMZ");
    meanDZ <- mxMatrix(type="Full", nrow=1, ncol=2, free=TRUE, 
                       values=means_init_vals,
                       labels=c("mb01", "mb01"),
                       name="MeanDZ");
    
    #
    defSex1    <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.sex1"),  name="Sex1");
    defSex2    <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.sex2"),  name="Sex2");
    defAge1    <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age1"),  name="Age1");
    defAge2    <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age2"),  name="Age2");
    defAgeSqr1 <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.satm1"), name="AgeSqr1");
    defAgeSqr2 <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.satm2"), name="AgeSqr2");
    defIcv1    <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.icv1"),  name="ICV1");
    defIcv2    <- mxMatrix(type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.icv2"),  name="ICV2");
    
    # Regression coefficients sex, age, age^2, and icv, and for interaction terms:
    pathBsex    <- mxMatrix(type="Full", nrow=1, ncol=2, free=TRUE, values=0.0, labels=c("b11","b11"), name="bsex");
    pathBage    <- mxMatrix(type="Full", nrow=1, ncol=2, free=TRUE, values=0.0, labels=c("b22","b22"), name="bage");
    pathBsatm   <- mxMatrix(type="Full", nrow=1, ncol=2, free=choice_include_age2, values=0.0, labels=c("b33","b33"), name="bsatm");
    pathBicv    <- mxMatrix(type="Full", nrow=1, ncol=2, free=TRUE, values=0.0, labels=c("b44","b44"), name="bicv");
    pathBagesex <- mxMatrix(type="Full", nrow=1, ncol=2, free=choice_include_age_x_sex, values=0.0, labels=c("b55","b55"), name="bagesex");
    pathBicvsex <- mxMatrix(type="Full", nrow=1, ncol=2, free=choice_include_icv_x_sex, values=0.0, labels=c("b66","b66"), name="bicvsex");
    
    # Mean corrected for covariates:
    correctedMeanMZ <- mxAlgebra(expression=
                                   cbind(
                                     MeanMZ + (Sex1 %*% bsex) + (Age1 %*% bage) + (AgeSqr1 %*% bsatm) + (ICV1 %*% bicv) + (Age1 %*% Sex1 %*% bagesex) + (ICV1 %*% Sex1 %*% bicvsex),
                                     MeanMZ + (Sex2 %*% bsex) + (Age2 %*% bage) + (AgeSqr2 %*% bsatm) + (ICV2 %*% bicv) + (Age2 %*% Sex2 %*% bagesex) + (ICV2 %*% Sex2 %*% bicvsex)),
                                 name="correctedMeanMZ");
    correctedMeanDZ <- mxAlgebra(expression=
                                   cbind(
                                     MeanDZ + (Sex1 %*% bsex) + (Age1 %*% bage) + (AgeSqr1 %*% bsatm) + (ICV1 %*% bicv) + (Age1 %*% Sex1 %*% bagesex) + (ICV1 %*% Sex1 %*% bicvsex),
                                     MeanDZ + (Sex2 %*% bsex) + (Age2 %*% bage) + (AgeSqr2 %*% bsatm) + (ICV2 %*% bicv) + (Age2 %*% Sex2 %*% bagesex) + (ICV2 %*% Sex2 %*% bicvsex)),
                                 name="correctedMeanDZ");
    
    
    ##
    ## Explanations (refer to the code below) ##
    ##
    ## We define the data sets and the likelihood function, and we assemble the model into a single `OpenMx` model called `model_fm`.
    ##
    
    # Data:
    data_mz <- mxData(observed=dmz, type="raw");
    data_dz <- mxData(observed=ddz, type="raw");
    
    # Expectation objects for Multiple Groups:
    exp_mz <- mxExpectationNormal(covariance="SMZ", means="correctedMeanMZ",
                                  dimnames=c(colnames(dmz)[i], colnames(dmz)[i+n_vars], colnames(dmz)[i+2*n_vars], colnames(dmz)[i+3*n_vars]));
    exp_dz <- mxExpectationNormal(covariance="SDZ", means="correctedMeanDZ",
                                  dimnames=c(colnames(ddz)[i], colnames(ddz)[i+n_vars], colnames(ddz)[i+2*n_vars], colnames(ddz)[i+3*n_vars]));
    
    #
    bits <- c(defSex1, defSex2, defAge1, defAge2, defAgeSqr1, defAgeSqr2, defIcv1, defIcv2,
              pathBsex, pathBage, pathBsatm, pathBicv, pathBagesex, pathBicvsex);
    pars <- list(Phmz, Phdz, LY, TE);
    
    # ML function and model:
    funML <- mxFitFunctionML();
    model_mz <- mxModel(meanMZ, correctedMeanMZ, bits, pars, covMZ, corMZ, data_mz, funML, exp_mz, name="MZ");
    model_dz <- mxModel(meanDZ, correctedMeanDZ, bits, pars, covDZ, corDZ, data_dz, funML, exp_dz, name="DZ");
    
    # Combine Groups:
    multi <- mxFitFunctionMultigroup(c("MZ","DZ"));
    model_fm <- mxModel("Phenotypic", pars, model_mz, model_dz, funML, multi);
    
    # Fit model:
    cat("Fitting Phenotypic model [model_fm]...\n", file=progress_filename, append=TRUE);
    fit <- fit_openmx_model_and_compare(model_fm, get_summary=TRUE);
    fit_fm <- fit$fit; summary_acde <- fit$summary;

    # Summary statistics (use the 1 rater error variance model):
    SMZ1 <- fit_fm$MZ$SMZ1$values; colnames(SMZ1) <- rownames(SMZ1) <- c("twin_1", "twin_2"); # conditional MZ latent phenotype
    SDZ1 <- fit_fm$DZ$SDZ1$values; colnames(SDZ1) <- rownames(SDZ1) <- c("twin_1", "twin_2"); # conditional DZ latent phenotype
    # Sanity checks:
    if( any(eigen(SMZ1)$values < 0 ) || any(eigen(SDZ1)$values < 0 ) ){ warning("Phenotypic covariance matrix has negative eigenvalue!\n"); }
    rMZ1 <- fit_fm$MZ$RMZ$result[1,2]; # MZ latent phenotype correlation
    rDZ1 <- fit_fm$DZ$RDZ$result[1,2]; # DZ latent phenotype correlation
    tv <- diag((fit_fm$MZ$SMZ$result)); # total variance
    ev <- diag(fit_fm$MZ$TE$values); # error variance
    err_vR2 <- ev/tv; # proportion of rater variance
    
    
    ##
    ## Explanations (refer to the code below) ##
    ##
    ## Here, the ACE or ADE decomposition of the latent cov(f1,f2) is carried out. 
    ## In the phenotypic model we estimated the covariance matrix of f1 and f2 in the MZ and the DZ twins; now we are going to model them. 
    ## The path diagram is in figure SEM_path_02.jpg.
    ##
    
    # ACE/ADE model (heuristic decision between ACE and ADE models based on r_{MZ} and r_{DZ}):
    if( is.na(rMZ1) || is.na(rDZ1) )
    {
      comments <- "Twin correlations estimation from phenotypic model failed: falling back to raw correlations.";
      genetic_model <- ifelse((rawd$rMZ < 2*rawd$rDZ), # r_{MZ} vs 2r_{DZ}
                              "ACE",  # ACE model
                              "ADE"); # ADE model
    } else
    {
      genetic_model <- ifelse((rMZ1 < 2*rDZ1), # r_{MZ} vs 2r_{DZ}
                               "ACE",  # ACE model
                               "ADE"); # ADE model
      
    }
    
    
    ##
    ## Explanations (refer to the code below) ##
    ##
    ## In the phenotypic model we estimated the covariance matrix of f1 and f2 in the MZ and the DZ twins. 
    ## These are modeled using the code below.
    ## That is:
    ##
    ##            phenotypic model	  ACE model
    ##    SMZ1 = 	v1	  cmz	        =	a^2 + c^2 + e^2	    a^2 + c^2 		
    ##            cmz	  v1		        a^2 + c^2 		      a^2 + c^2 + e^2
    ##
    ## and:
    ##
    ##            phenotypic model	  ACE model
    ##    SDZ1 = 	v1	  cmz	        =	a^2 + c^2 + e^2	    .5a^2 + c^2 		
    #             cmz	  v1		        .5a^2 + c^2 		    a^2 + c^2 + e^2
    ##
    ## The parameters a, c and e are defined in the code below (as `A`, `CD` and `E`, respectively).
    ##
    
    # Variance decomposition:
    A  <- mxMatrix(type="Lower", nrow=1, ncol=1, free=TRUE,  values=c(.5), labels=c("a11"), name="a" );
    CD <- mxMatrix(type="Lower", nrow=1, ncol=1, free=TRUE,  values=c(.5), labels=c("c11"), name="cd"); # C or D
    E  <- mxMatrix(type="Lower", nrow=1, ncol=1, free=TRUE,  values=c(.5), labels=c("e11"), name="e" );
    
    
    ##
    ## Explanations (refer to the code below) ##
    ##
    ## We calculate the variances and the total phenotypic variance.
    ##
    
    vA  <- mxAlgebra(a * a,      name="A");
    vCD <- mxAlgebra(cd * cd,    name="CD");
    vE  <- mxAlgebra(e * e,      name="E")
    vPh <- mxAlgebra(A + CD + E, name="Ph");
    
    
    ##
    ## Explanations (refer to the code below) ##
    ##
    ## With these in place, the latents `SMZ1` and `SDZ1` are formed, and used to calculate the phenotypic covariance matrices. 
    ##
    
    # MZ 2x2 covariance matrix:
    Phmz <- mxAlgebra(expression=
                        cbind(rbind(Ph,     A + CD),
                              rbind(A + CD, Ph    )),
                      name="SMZ1");
    # DZ 2x2 covariance matrix:
    if( genetic_model == "ACE" ) # explicitly define the Phdz matrix focusing of the model to allow %dopar% to work correctly
    {
      Phdz <- mxAlgebra(expression=
                          cbind(rbind(Ph,              0.5*A + 1.00*CD),
                                rbind(0.5*A + 1.00*CD, Ph             )),
                        name="SDZ1");
    } else # ADE
    {
      Phdz <- mxAlgebra(expression=
                          cbind(rbind(Ph,              0.5*A + 0.25*CD),
                                rbind(0.5*A + 0.25*CD, Ph             )),
                        name="SDZ1");
    }
    
    
    ##
    ## Explanations (refer to the code below) ##
    ##
    ## We calculate the phenotypic covariance matrices, and the latent correlation matrices (`RMZ` and `RDZ`).
    ##
    
    # Matrices to calculate A, C and E variance components:
    covMZ <- mxAlgebra(LY %*% SMZ1 %*% t(LY) + TE, name="SMZ");
    covDZ <- mxAlgebra(LY %*% SDZ1 %*% t(LY) + TE, name="SDZ");
    corMZ <- mxAlgebra(cov2cor(SMZ1), name="RMZ");
    corDZ <- mxAlgebra(cov2cor(SDZ1), name="RDZ");
    
    
    ##
    ## Explanations (refer to the code below) ##
    ##
    ## We calculate various statistics and informative representations of the results, as explained in the annotation of the code. 
    ##
    
    # Important things to estimate:
    stcomp <- mxAlgebra(expression=
                          cbind(
                            # The proportions of *standardized* latent phenotypic variances A, CD, A+CD and E:
                            A/Ph,      # standardized A variance, aka the narrow-sense heritability: var(A) / var(latent phenotype)                    
                            CD/Ph,     # standardized C or D variance:  var(CD) / var(latent phenotype)
                            (A+CD)/Ph, # if CD==C, total familiality, else then this is the broad-sense heritability
                            E/Ph,      # standardized unshared environmental variance: var(E)/var(latent phenotype)
                            
                            # The proportions of observed rater 1 phenotypic variances A, CD, A+CD and E, and rater variances:
                            A/SMZ[1,1],     # as A/Ph above but as proportion of the observed phenotypic variance var(phenotype), not var(latent phenotype) 
                            CD/SMZ[1,1],    # as above CD/Ph but as proportion of the observed phenotypic variance var(phenotype), not var(latent phenotype) 
                            (A+CD)/SMZ[1,1],# as (A+CD)/Ph above but as proportion of the observed phenotypic variance var(phenotype), not var(latent phenotype) 
                            E/SMZ[1,1],     # as E/Ph above but as proportion of the observed phenotypic variance var(phenotype), not var(latent phenotype)      
                            TE[1,1]/SMZ[1,1],   # the standardized rater error variance, i.e., (1 - rater reliability) 
                            (A+CD+E)/SMZ[1,1]#, # the rater reliability: {(A+CD+E)/SMZ[1,1]} = 1 - (TE[1,1]/SMZ[1,1])
                            
                            ## Same as above but for rater 2 (only of interest if LY[2,1] and LY[4,2] are estimated, which is not the case for us here):
                            #(A*(LY[2,1]^2))/SMZ[2,2], 
                            #(CD*(LY[2,1]^2))/SMZ[2,2], 
                            #((A+CD)*(LY[2,1]^2))/SMZ[2,2], 
                            #(E*(LY[2,1]^2))/SMZ[2,2], 
                            #TE[1,1]/SMZ[2,2],
                            #(A+CD+E)/SMZ[2,2]
                          ),
                        name="stvc");
    
    
    ##
    ## Explanations (refer to the code below) ##
    ##
    ## Finally, the model is assembled and run.
    ##
    
    # Expectation objects for Multiple Groups:
    exp_mz <- mxExpectationNormal(covariance="SMZ", means="correctedMeanMZ",
                                  dimnames=c(colnames(dmz)[i], colnames(dmz)[i+n_vars], colnames(dmz)[i+2*n_vars], colnames(dmz)[i+3*n_vars]));
    exp_dz <- mxExpectationNormal(covariance="SDZ", means="correctedMeanDZ",
                                  dimnames=c(colnames(ddz)[i], colnames(ddz)[i+n_vars], colnames(ddz)[i+2*n_vars], colnames(ddz)[i+3*n_vars]));
    
    #
    bits <- c(defSex1, defSex2, defAge1, defAge2, defAgeSqr1, defAgeSqr2, defIcv1, defIcv2,
              pathBsex, pathBage, pathBsatm, pathBicv, pathBagesex, pathBicvsex);
    pars <- list(Phmz, Phdz, LY, TE, A, CD, E, vA, vCD, vE, vPh);
    
    # ML function and model:
    funML <- mxFitFunctionML();
    CI    <- mxCI(c("ACorDE.A", "ACorDE.CD", "ACorDE.E", "MZ.stvc"));
    model_mz <- mxModel(CI, stcomp, meanMZ, correctedMeanMZ, bits, pars, covMZ, corMZ, data_mz, funML, exp_mz, name="MZ");
    model_dz <- mxModel(meanDZ, correctedMeanDZ, bits, pars, covDZ, corDZ, data_dz, funML, exp_dz, name="DZ");
    
    # Combine Groups:
    multi <- mxFitFunctionMultigroup(c("MZ","DZ"));
    
    # Model with variance decomposition:
    model_acde_vardecomp <- mxModel("ACorDE", pars, model_mz, model_dz, funML, CI, multi);
    
    # Run models:
    cat("Fitting ACorDE model [model_acde_vardecomp]...\n", file=progress_filename, append=TRUE);
    fit <- fit_openmx_model_and_compare(model_acde_vardecomp, get_summary=TRUE);
    fit_acde_vardecomp <- fit$fit; summary_acde_vardecomp <- fit$summary;

    
    # Test Significance of A (by removing it from the model):
    model_noA <- omxSetParameters(fit_acde_vardecomp, labels=c("a11"), free=FALSE, values=0, name="NoA");
    cat("Fitting ACorDE model without A [model_noA]...\n", file=progress_filename, append=TRUE);
    fit <- fit_openmx_model_and_compare(model_noA, fit_acde_vardecomp, model1_embedded_in_model2=TRUE, get_summary=TRUE);
    fit_noA <- fit$fit;  summary_noA <- fit$summary; comp_A_vs_noA <- fit$comparison;
    
    # Test Significance of CD (by removing it from the model):
    model_noCD <- omxSetParameters(fit_acde_vardecomp, labels=c("c11"), free=FALSE, values=0, name="NoCD");
    cat("Fitting ACorDE model without CD [model_noCD]...\n", file=progress_filename, append=TRUE);
    fit <- fit_openmx_model_and_compare(model_noCD, fit_acde_vardecomp, model1_embedded_in_model2=TRUE, get_summary=TRUE);
    fit_noCD <- fit$fit;  summary_noCD <- fit$summary; comp_CD_vs_noCD <- fit$comparison;
    
    # Test Significance of E (by removing it from the model):
    model_noE <- omxSetParameters(fit_acde_vardecomp, labels=c("e11"), free=FALSE, values=0, name="NoE");
    cat("Fitting ACorDE model without E [model_noE]...\n", file=progress_filename, append=TRUE);
    fit <- fit_openmx_model_and_compare(model_noE, fit_acde_vardecomp, model1_embedded_in_model2=TRUE, get_summary=TRUE);
    fit_noE <- fit$fit;  summary_noE <- fit$summary; comp_E_vs_noE <- fit$comparison;
    
    
    # Return results:
    cat("Assembling return value...\n", file=progress_filename, append=TRUE);
    cis <- summary_acde_vardecomp$CI; # use a shorter name for the CIs
    return (data.frame(
      # Phenotype info:
      "phenotype_number"    =i, 
      "phenotype_name"      =stringr::str_sub(varnames[i], end=-4), 
      "phenotype_name_short"=phenotype_names$short.names[ toupper(stringr::str_sub(phenotype_names$long.names.rater1, end=-3)) == toupper(stringr::str_sub(varnames[i], end=-4)) ],

      # Fit measures:
      "numObs"=summary_acde_vardecomp$numObs, "numParams"=summary_acde_vardecomp$estimatedParameters, 
      "observedStatistics"=summary_acde_vardecomp$observedStatistics, "degreesOfFreedom"=summary_acde_vardecomp$degreesOfFreedom,
      "AIC"=summary_acde_vardecomp$AIC.Mx, "BIC"=summary_acde_vardecomp$BIC.Mx,
      "fit"=summary_acde_vardecomp$fit, #"CFI"=summary_acde_vardecomp$CFI, "TLI"=summary_acde_vardecomp$TLI, "RMSEA"=summary_acde_vardecomp$RMSEA, 
      "statusCode"=summary_acde_vardecomp$statusCode,
      
      # Relevant parameter estimates and 95%CIs:
      # Which genetic model is this?
      "genetic_model"=genetic_model, 
      "CD_meaning"=ifelse(genetic_model == "ACE", "C", "D"), 
      "H2_meaning"=ifelse(genetic_model == "ACE", "familiality", "broad sense heritability"),
      "twin_correlations"=ifelse(is.na(rMZ1) || is.na(rDZ1), "raw", "corrected"),
      "twin_corr_rMZ"=ifelse(is.na(rMZ1) || is.na(rDZ1), rawd$rMZ, rMZ1), "twin_corr_rDZ"=ifelse(is.na(rMZ1) || is.na(rDZ1), rawd$rDZ, rDZ1), 
      
      # Non-standardized A, C or D (as per CD_meaning), and E:
      "A"=cis["ACorDE.A[1,1]", "estimate"], "A_CiL"=cis["ACorDE.A[1,1]", "lbound"], "A_CiH"=cis["ACorDE.A[1,1]", "ubound"], 
      "A_deltaAIC"=(comp_A_vs_noA$AIC[2] - comp_A_vs_noA$AIC[1]), "A_p"=comp_A_vs_noA$p[2], # deltaAIC <0 iff NoA better than A
      "CD"=cis["ACorDE.CD[1,1]", "estimate"], "CD_CiL"=cis["ACorDE.CD[1,1]", "lbound"], "CD_CiH"=cis["ACorDE.CD[1,1]", "ubound"], 
      "CD_deltaAIC"=(comp_CD_vs_noCD$AIC[2] - comp_CD_vs_noCD$AIC[1]), "CD_p"=comp_CD_vs_noCD$p[2],
      "E"=cis["ACorDE.E[1,1]", "estimate"], "E_CiL"=cis["ACorDE.E[1,1]", "lbound"], "E_CiH"=cis["ACorDE.E[1,1]", "ubound"], 
      "E_deltaAIC"=(comp_E_vs_noE$AIC[2] - comp_E_vs_noE$AIC[1]), "E_p"=comp_E_vs_noE$p[2],
      
      # Latent h2, cd2, H2 (as per H2_meaning) and e2, i.e., the standardized variance terms corrected for rater effects:
      "h2" =cis["MZ.stvc[1,1]", "estimate"], "h2_CiL" =cis["MZ.stvc[1,1]", "lbound"], "h2_CiH" =cis["MZ.stvc[1,1]", "ubound"],
      "cd2"=cis["MZ.stvc[1,2]", "estimate"], "cd2_CiL"=cis["MZ.stvc[1,2]", "lbound"], "cd2_CiH"=cis["MZ.stvc[1,2]", "ubound"],
      "H2" =cis["MZ.stvc[1,3]", "estimate"], "H2_CiL" =cis["MZ.stvc[1,3]", "lbound"], "H2_CiH" =cis["MZ.stvc[1,3]", "ubound"],
      "e2" =cis["MZ.stvc[1,4]", "estimate"], "e2_CiL" =cis["MZ.stvc[1,4]", "lbound"], "e2_CiH" =cis["MZ.stvc[1,4]", "ubound"],
      
      # Rater effects: A, C or D (as per CD_meaning), E and the rater effect (R):
      "rh2"  =cis["MZ.stvc[1,5]", "estimate"], "rh2_CiL"  =cis["MZ.stvc[1,5]", "lbound"], "rh2_CiH"  =cis["MZ.stvc[1,5]", "ubound"],
      "rcd2" =cis["MZ.stvc[1,6]", "estimate"], "rcd2_CiL" =cis["MZ.stvc[1,6]", "lbound"], "rcd2_CiH" =cis["MZ.stvc[1,6]", "ubound"],
      "rH2"  =cis["MZ.stvc[1,7]", "estimate"], "rH2_CiL"  =cis["MZ.stvc[1,7]", "lbound"], "rH2_CiH"  =cis["MZ.stvc[1,7]", "ubound"],
      "re2"  =cis["MZ.stvc[1,8]", "estimate"], "re2_CiL"  =cis["MZ.stvc[1,8]", "lbound"], "re2_CiH"  =cis["MZ.stvc[1,8]", "ubound"],
      "rerr2"=cis["MZ.stvc[1,9]", "estimate"], "rerr2_CiL"=cis["MZ.stvc[1,9]", "lbound"], "rerr2_CiH"=cis["MZ.stvc[1,9]", "ubound"],
      
      # Rater reliability (1 - standardized rater error variance):
      "rreliability"=cis["MZ.stvc[1,10]", "estimate"], "rreliability_CiL"=cis["MZ.stvc[1,10]", "lbound"], "rreliability_CiH"=cis["MZ.stvc[1,10]", "ubound"],
      
      # Comments and/or errors:
      "Comments"=comments
    ));
  }
results_df <- do.call(rbind, results);

# Add the raw MZ and DZ correlations:
results_df <- merge(results_df, raw_checks_df[,c("phentype_name", "rMZ", "rDZ")], by.x="phenotype_name", by.y="phentype_name", all.x=TRUE);

# "Classic" h2 estimate based on raw correlations:
results_df$classic_h2 <- 2*(results_df$rMZ - results_df$rDZ);

## Save results to file ####
write.csv(results_df, "../data/intermediate/SEM_results.csv", row.names=FALSE);


# Save the explanation of the results:
cat(paste0('
# This document explains the content of the `SEM_results.csv` file.\n
\n
*N.B.* this document uses GitHub-flavoured Markdown, in particular \\<sup\\>x\\</sup\\> for superscript <sup>x</sup> and \\<sub\\>x\\</sub\\> for subscript <sub>x</sub>.
The columns are:\n
\n
## General:\n
- **phenotype\\_name**: the full phenotype name\n
- **phenotype\\_number**: the phenotype number\n
- **phenotype\\_name\\_short**: the phenotype short name\n
\n
## For the genetic decomposition model:\n
- **numObs**, **numParams**, **observedStatistics**, **degreesOfFreedom**, **AIC**, **BIC**, **fit**, **statusCode**: number of observations, of parameters, observed statistics, degrees of freedom, AIC, BIC, fit and error code (if any)\n
- **genetic\\_model**, **CD\\_meaning**, and **H2\\_meaning**: given that we cannot estimate simultaneously the *C* and *D* components, we must chose between the two using the heuristic [if *r*<sub>MZ</sub> < 2*r*<sub>DZ</sub> then we estimate *C* otherwise we estimate *D*]; this, in turn, determines the fitted genetic model (*ACE* or *ADE*), the meaning of the "CD" (short for "C or D") parameter (*C* or *D*), and the meaning of the "H2" (actually *H*<sup>2</sup>) estimate ("familiality" for *C* and "broad sense heritability" for *D*, respectively);\n
- Corrected or raw twin correlations used to decide between the ADE and ACE models:\n
  + **twin_correlations**: the raw or the corrected twin correlations were used?\n
  + **twin_corr_rMZ**, **twin_corr_rDZ**: the MZ and DZ correlations\n
- for each of **A**, **CD** and **E**:
  + the point estimate (i.e., of *A*, (*C* or *D*), and *E*, respectively)\n
  + **\\_CiL**, **\\_CiH**: the lower and higher (upper) bounds of the 95% confidence intervals (CIs) of the corresponding point estimate (i.e., of *A*, (*C* or *D*), and *E*, respectively); `NA` should probably be interpreted as 0.0\n
  + **\\_deltaAIC** and **\\_p**: the results of the comparison between the model with the component (e.g., *A*) and without it; the AIC difference is negative if including the component fits the data better than when excluding it; the *p* is significant if including the component and excluding it are significantly different\n
- The latent (standardized)  estimates of relevant parameters, i.e., the standardized variance terms corrected for rater effects; for all, the point estimate and the lower and upper 95% CIs:\n
  + **h2**: the narrow heritability *h*<sup>2</sup>\n
  + **cd2**: *C*<sup>2</sup> of *D*<sup>2</sup>, respectively\n
  + **H2**: the familiality or broad sense heritability *H*<sup>2</sup>, respectively\n
  + **e2**: *E*<sup>2</sup>\n
- The rater-specific estimates of various parameters (with 95% CIs):
  + **rh2**: the narrow heritability *h*<sup>2</sup>\n
  + **rcd2**: *C*<sup>2</sup> of *D*<sup>2</sup>, respectively\n
  + **rH2**: the familiality or broad sense heritability *H*<sup>2</sup>, respectively\n
  + **re2**: *E*<sup>2</sup>\n
  + **rerr2**: the standardized rater error variance, i.e., (1 - rater reliability)\n
- The rater reliability:\n
  + **rreliability**: rater reliability (1 - standardized rater error variance)\n
- **Comments**: any issues chosing the genetic model?
\n
## For the "classic" (i.e., non-SEM) heritability estimates:\n
- **rMZ** and **rDZ**: the raw (i.e., uncorrected) observed phenotypic correlations between monozygotic *r*<sub>MZ</sub> and dizygotic *r*<sub>DZ</sub> twins\n
- **classic\\_h2**: the "classic" narrow heritability estimate *h*<sup>2</sup> = 2(*r*<sub>MZ</sub> - *r*<sub>DZ</sub>)\n
\n
      ', 
           "Script took ", difftime(Sys.time(), run_start_time, units="mins"), " minutes to run...\n"), 
    file="../data/intermediate/SEM_results_explanations.md", append=FALSE);



## END ####
cat("Took ", difftime(Sys.time(), run_start_time, units="mins"), " minutes...\n");

# Stop the clusters:
if( run_parallel )
{
  parallel::stopCluster(cl);
}

unlink("../logs/progress_phenotype_*.txt"); # clean the logs....

