# Fit a bivariate GCSM model for vocal tract heritability estimates
# on the dataset with all DZ twins included, using lavaan
#
# Copyright (C) 2022, Dan Dediu
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

# Check if the analysis can be run at all...
if( !file.exists("../data/input/measures.csv") )
{
  stop("Please note that for running the analysis, the file '../data/input/measures.csv' must first be obtained from the Netherlands Twin Register (https://tweelingenregister.vu.nl/information_for_researchers/working-with-ntr-data). This file cannot be made public as it contains potentially identifying information about the participants ('../data/input/measures-STUB.csv' is an empty stub to show the structure of this file)...");
}


run_start_time <- Sys.time(); # start time: takes about 40 minutes on a Core i9-10900 (10 cores) with 32Gb RAM and abut 30 mins on a Ryzen 3700X (8 cores) with 64GB RAM


## Libraries ####
library(stringr);
library(foreach);
library(dplyr);
library(psych);
library(lavaan);
library(lavaanPlot);


# Notation conventions: 
# - cervicalSpine_CS2H_atlas_c2base_distance.1_1 is phenotype.rater1_twin1. All phenotypes follow this notation. So:
# - phenotype.1_1 = phenotype.rater1_twin1
# - phenotype.2_1 = phenotype.rater2_twin1
# - phenotype.1_2 = phenotype.rater1_twin2
# - phenotype.2_2 = phenotype.rater2_twin2


## Load data and adjust variables ####
if( !file.exists("../data/intermediate/data_for_SEM.csv") ) source("./1_data_preparation.R", echo=FALSE); # run the data preparation first...
SEMdata <- read.csv("../data/intermediate/data_for_SEM.csv");
SEMdata$X <- NULL; # rownames are not needed
SEMdata <- SEMdata %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance.1_1":"general_GAPD_allPointSet_matlabProcDistance.2_1"),
                              c("cervicalSpine_CS2H_atlas_c2base_distance.1_2":"general_GAPD_allPointSet_matlabProcDistance.2_2"),
                              c("Sex_1":"ZYGMZDZ"),"FISNumber_1","FISNumber_2");

# z-score all phenotypic measures for the whole dataset:
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



## Fit the genetic SEM models ####
cat(">>> Fit the genetic SEM models...\n");
results <- foreach(i=1:n_vars, .packages=c("OpenMx"), .errorhandling="stop") %do% # fast enough for single threading...
  {
    cat("Fitting ", varnames[i],"(",i,"of",n_vars,")...\n"); # print progress

    # Arrange the data in the format (where T=twin member, and R=rater):
    # Zygosity  T1R1  T1R2  T2R1  T2R2  covariates...
    d <- rbind(data.frame("zyg"="MZ", 
                          "T1R1"=dmz[,i], "T1R2"=dmz[,i+n_vars], "T2R1"=dmz[,i+2*n_vars], "T2R2"=dmz[,i+3*n_vars],
                          "T1sex"=dmz$sex1, "T1age"=dmz$age1, "T1age2"=dmz$satm1, "T1icv"=dmz$icv1,
                          "T2sex"=dmz$sex2, "T2age"=dmz$age2, "T2age2"=dmz$satm2, "T2icv"=dmz$icv2),
               data.frame("zyg"="DZ", "T1R1"=ddz[,i], "T1R2"=ddz[,i+n_vars], "T2R1"=ddz[,i+2*n_vars], "T2R2"=ddz[,i+3*n_vars],
                          "T1sex"=ddz$sex1, "T1age"=ddz$age1, "T1age2"=ddz$satm1, "T1icv"=ddz$icv1,
                          "T2sex"=ddz$sex2, "T2age"=ddz$age2, "T2age2"=ddz$satm2, "T2icv"=ddz$icv2));
    
    # Residualize the phenotype on the covariates (trying to include the covariates in the SEM model does not work):
    d$T1R1rez <- residuals(lm(T1R1 ~ T1sex + T1age + T1age2 + T1icv, data=d, na.action=na.exclude));
    d$T1R2rez <- residuals(lm(T1R2 ~ T1sex + T1age + T1age2 + T1icv, data=d, na.action=na.exclude));
    d$T2R1rez <- residuals(lm(T2R1 ~ T2sex + T2age + T2age2 + T2icv, data=d, na.action=na.exclude));
    d$T2R2rez <- residuals(lm(T2R2 ~ T2sex + T2age + T2age2 + T2icv, data=d, na.action=na.exclude));
    
    # ACE:
    ace.model  <- " # ACE model (using reziduals after controllng for covariates):
    
                    # Raters:
                    T1 =~ c(r1,r1)*T1R1rez + c(r2,r2)*T1R2rez # twin 1 is measured by both raters
                    T2 =~ c(r1,r1)*T2R1rez + c(r2,r2)*T2R2rez # twin 2 is measured by both raters

                    # Latents:
                    A1 =~ NA*T1 + c(a,a)*T1 # A for twin 1
                    A2 =~ NA*T2 + c(a,a)*T2 # A for twin 2 
                    C1 =~ NA*T1 + c(c,c)*T1 # C for twin 1 
                    C2 =~ NA*T2 + c(c,c)*T2 # C for twin 2 

                    # Variances:
                    A1 ~~ 1*A1 # A has variance 1
                    A2 ~~ 1*A2 # A has variance 1
                    C1 ~~ 1*C1 # C has variance 1
                    C2 ~~ 1*C2 # C has variance 1
                    
                    T1 ~~ c(e2,e2)*T1 # E is the residual variance of the phenotype
                    T2 ~~ c(e2,e2)*T2 # E is the residual variance of the phenotype
                    
                    # Covariances
                    A1 ~~ c(1,.5)*A2 # for A, the correlation is 1 for MZ twins and 0.5 for DZ twins
                    A1 ~~ 0*C1 + 0*C2 # A and C are uncorrelated
                    A2 ~~ 0*C1 + 0*C2 # A and C are uncorrelated 
                    C1 ~~ c(1,1)*C2 # C is correlated 1 regardless of twin status
                                    # by lavaan default, E is uncorrelated with A and C, so no need to model it
                    ";
    ace.fit <- sem(ace.model, data=d, group="zyg", se="robust", optim.method="BFGS");
    #fitMeasures(ace.fit, c("chisq", "df", "pvalue", "aic", "bic", "rmsea", "srmr")); summary(ace.fit , standardized=TRUE); 
    #lavaanPlot(model=ace.fit, node_options=list(shape="box", fontname="Helvetica"), coefs=TRUE, covs=TRUE, stars="latent");

    # ADE:
    ade.model  <- " # ADE model (using reziduals after controllng for covariates):
    
                    # Raters:
                    T1 =~ c(r1,r1)*T1R1rez + c(r2,r2)*T1R2rez # twin 1 is measured by both raters
                    T2 =~ c(r1,r1)*T2R1rez + c(r2,r2)*T2R2rez # twin 2 is measured by both raters

                    # Latents:
                    A1 =~ NA*T1 + c(a,a)*T1 # A for twin 1
                    A2 =~ NA*T2 + c(a,a)*T2 # A for twin 2 
                    D1 =~ NA*T1 + c(d,d)*T1 # D for twin 1 
                    D2 =~ NA*T2 + c(d,d)*T2 # D for twin 2 

                    # Variances:
                    A1 ~~ 1*A1 # A has variance 1
                    A2 ~~ 1*A2 # A has variance 1
                    D1 ~~ 1*D1 # D has variance 1
                    D2 ~~ 1*D2 # D has variance 1
                    
                    T1 ~~ c(e2,e2)*T1 # E is the residual variance of the phenotype
                    T2 ~~ c(e2,e2)*T2 # E is the residual variance of the phenotype
                    
                    # Covariances
                    A1 ~~ c(1,.5)*A2 # for A, the correlation is 1 for MZ twins and 0.5 for DZ twins
                    A1 ~~ 0*D1 + 0*D2 # A and D are uncorrelated
                    A2 ~~ 0*D1 + 0*D2 # A and D are uncorrelated 
                    D1 ~~ c(1,.25)*D2 # non-additive genetic effects D are correlated 1 for MZ twins and .25 for DZ
                                      # by lavaan default, E is uncorrelated with A and D, so no need to model it
                    ";
    ade.fit <- sem(ade.model, data=d, group="zyg", se="robust", optim.method="BFGS");
    #fitMeasures(ace.fit, c("chisq", "df", "pvalue", "aic", "bic", "rmsea", "srmr")); summary(ace.fit , standardized=TRUE); 
    #lavaanPlot(model=ade.fit, node_options=list(shape="box", fontname="Helvetica"), coefs=TRUE, covs=TRUE, stars="latent");

    # ACE or ADE? Use AIC for this comparison:
    if( (deltaAIC_ACE_ADE <- AIC(ace.fit) - AIC(ade.fit)) <= 0 )
    {
      # ACE:
      genetic_model <- "ACE"; model_fit <- ace.fit;
      
      # Test A:
      ce.model   <- " # CE model (using reziduals after controllng for covariates):
    
                    # Raters:
                    T1 =~ c(r1,r1)*T1R1rez + c(r2,r2)*T1R2rez # twin 1 is measured by both raters
                    T2 =~ c(r1,r1)*T2R1rez + c(r2,r2)*T2R2rez # twin 2 is measured by both raters

                    # Latents:
                    #A1 =~ NA*T1 + c(a,a)*T1 # A for twin 1
                    #A2 =~ NA*T2 + c(a,a)*T2 # A for twin 2 
                    C1 =~ NA*T1 + c(c,c)*T1 # C for twin 1 
                    C2 =~ NA*T2 + c(c,c)*T2 # C for twin 2 

                    # Variances:
                    #A1 ~~ 1*A1 # A has variance 1
                    #A2 ~~ 1*A2 # A has variance 1
                    C1 ~~ 1*C1 # C has variance 1
                    C2 ~~ 1*C2 # C has variance 1
                    
                    T1 ~~ c(e2,e2)*T1 # E is the residual variance of the phenotype
                    T2 ~~ c(e2,e2)*T2 # E is the residual variance of the phenotype
                    
                    # Covariances
                    #A1 ~~ c(1,.5)*A2 # for A, the correlation is 1 for MZ twins and 0.5 for DZ twins
                    #A1 ~~ 0*C1 + 0*C2 # A and C are uncorrelated
                    #A2 ~~ 0*C1 + 0*C2 # A and C are uncorrelated 
                    C1 ~~ c(1,1)*C2 # C is correlated 1 regardless of twin status
                                    # by lavaan default, E is uncorrelated with A and C, so no need to model it
                    ";
      ce.fit <- sem(ce.model, data=d, group="zyg", se="robust", optim.method="BFGS");
      #fitMeasures(ce.fit, c("chisq", "df", "pvalue", "aic", "bic", "rmsea", "srmr")); summary(ace.fit , standardized=TRUE); 
      #lavaanPlot(model=ce.fit, node_options=list(shape="box", fontname="Helvetica"), coefs=TRUE, covs=TRUE, stars="latent");
      deltaAIC_A <- AIC(ace.fit) - AIC(ce.fit); LRtest_A <- lavTestLRT(ace.fit, ce.fit);
      
      # Test C:
      ae.model   <- " # AE model (using reziduals after controllng for covariates):
    
                    # Raters:
                    T1 =~ c(r1,r1)*T1R1rez + c(r2,r2)*T1R2rez # twin 1 is measured by both raters
                    T2 =~ c(r1,r1)*T2R1rez + c(r2,r2)*T2R2rez # twin 2 is measured by both raters

                    # Latents:
                    A1 =~ NA*T1 + c(a,a)*T1 # A for twin 1
                    A2 =~ NA*T2 + c(a,a)*T2 # A for twin 2 
                    #C1 =~ NA*T1 + c(c,c)*T1 # C for twin 1 
                    #C2 =~ NA*T2 + c(c,c)*T2 # C for twin 2 

                    # Variances:
                    A1 ~~ 1*A1 # A has variance 1
                    A2 ~~ 1*A2 # A has variance 1
                    #C1 ~~ 1*C1 # C has variance 1
                    #C2 ~~ 1*C2 # C has variance 1
                    
                    T1 ~~ c(e2,e2)*T1 # E is the residual variance of the phenotype
                    T2 ~~ c(e2,e2)*T2 # E is the residual variance of the phenotype
                    
                    # Covariances
                    A1 ~~ c(1,.5)*A2 # for A, the correlation is 1 for MZ twins and 0.5 for DZ twins
                    #A1 ~~ 0*C1 + 0*C2 # A and C are uncorrelated
                    #A2 ~~ 0*C1 + 0*C2 # A and C are uncorrelated 
                    #C1 ~~ c(1,1)*C2 # C is correlated 1 regardless of twin status
                                    # by lavaan default, E is uncorrelated with A and C, so no need to model it
                    ";
      ae.fit <- sem(ae.model, data=d, group="zyg", se="robust", optim.method="BFGS");
      #fitMeasures(ae.fit, c("chisq", "df", "pvalue", "aic", "bic", "rmsea", "srmr")); summary(ace.fit , standardized=TRUE); 
      #lavaanPlot(model=ae.fit, node_options=list(shape="box", fontname="Helvetica"), coefs=TRUE, covs=TRUE, stars="latent");
      deltaAIC_other <- AIC(ace.fit) - AIC(ae.fit); LRtest_other <- lavTestLRT(ace.fit, ae.fit);
    } else
    {
      # ADE:
      genetic_model <- "ADE"; model_fit <- ade.fit;

      # Test A:
      de.model   <- " # DE model (using reziduals after controllng for covariates):
    
                    # Raters:
                    T1 =~ c(r1,r1)*T1R1rez + c(r2,r2)*T1R2rez # twin 1 is measured by both raters
                    T2 =~ c(r1,r1)*T2R1rez + c(r2,r2)*T2R2rez # twin 2 is measured by both raters

                    # Latents:
                    #A1 =~ NA*T1 + c(a,a)*T1 # A for twin 1
                    #A2 =~ NA*T2 + c(a,a)*T2 # A for twin 2 
                    D1 =~ NA*T1 + c(d,d)*T1 # D for twin 1 
                    D2 =~ NA*T2 + c(d,d)*T2 # D for twin 2 

                    # Variances:
                    #A1 ~~ 1*A1 # A has variance 1
                    #A2 ~~ 1*A2 # A has variance 1
                    D1 ~~ 1*D1 # D has variance 1
                    D2 ~~ 1*D2 # D has variance 1
                    
                    T1 ~~ c(e2,e2)*T1 # E is the residual variance of the phenotype
                    T2 ~~ c(e2,e2)*T2 # E is the residual variance of the phenotype
                    
                    # Covariances
                    #A1 ~~ c(1,.5)*A2 # for A, the correlation is 1 for MZ twins and 0.5 for DZ twins
                    #A1 ~~ 0*D1 + 0*D2 # A and D are uncorrelated
                    #A2 ~~ 0*D1 + 0*D2 # A and D are uncorrelated 
                    D1 ~~ c(1,.25)*D2 # non-additive genetic effects D are correlated 1 for MZ twins and .25 for DZ
                                      # by lavaan default, E is uncorrelated with A and D, so no need to model it
                    ";
      de.fit <- sem(de.model, data=d, group="zyg", se="robust", optim.method="BFGS");
      #fitMeasures(de.fit, c("chisq", "df", "pvalue", "aic", "bic", "rmsea", "srmr")); summary(ace.fit , standardized=TRUE); 
      #lavaanPlot(model=de.fit, node_options=list(shape="box", fontname="Helvetica"), coefs=TRUE, covs=TRUE, stars="latent");
      deltaAIC_A <- AIC(ade.fit) - AIC(de.fit); LRtest_A <- lavTestLRT(ade.fit, de.fit);
      
      # Test D:
      ae.model   <- " # AE model (using reziduals after controllng for covariates):
    
                    # Raters:
                    T1 =~ c(r1,r1)*T1R1rez + c(r2,r2)*T1R2rez # twin 1 is measured by both raters
                    T2 =~ c(r1,r1)*T2R1rez + c(r2,r2)*T2R2rez # twin 2 is measured by both raters

                    # Latents:
                    A1 =~ NA*T1 + c(a,a)*T1 # A for twin 1
                    A2 =~ NA*T2 + c(a,a)*T2 # A for twin 2 
                    #D1 =~ NA*T1 + c(d,d)*T1 # D for twin 1 
                    #D2 =~ NA*T2 + c(d,d)*T2 # D for twin 2 

                    # Variances:
                    A1 ~~ 1*A1 # A has variance 1
                    A2 ~~ 1*A2 # A has variance 1
                    #D1 ~~ 1*D1 # D has variance 1
                    #D2 ~~ 1*D2 # D has variance 1
                    
                    T1 ~~ c(e2,e2)*T1 # E is the residual variance of the phenotype
                    T2 ~~ c(e2,e2)*T2 # E is the residual variance of the phenotype
                    
                    # Covariances
                    A1 ~~ c(1,.5)*A2 # for A, the correlation is 1 for MZ twins and 0.5 for DZ twins
                    #A1 ~~ 0*D1 + 0*D2 # A and D are uncorrelated
                    #A2 ~~ 0*D1 + 0*D2 # A and D are uncorrelated 
                    #D1 ~~ c(1,.25)*D2 # non-additive genetic effects D are correlated 1 for MZ twins and .25 for DZ
                                      # by lavaan default, E is uncorrelated with A and D, so no need to model it
                    ";
      ae.fit <- sem(ae.model, data=d, group="zyg", se="robust", optim.method="BFGS");
      #fitMeasures(ae.fit, c("chisq", "df", "pvalue", "aic", "bic", "rmsea", "srmr")); summary(ace.fit , standardized=TRUE); 
      #lavaanPlot(model=ae.fit, node_options=list(shape="box", fontname="Helvetica"), coefs=TRUE, covs=TRUE, stars="latent");
      deltaAIC_other <- AIC(ade.fit) - AIC(ae.fit); LRtest_other <- lavTestLRT(ade.fit, ae.fit);
    }

    model_fit_measures <- fitMeasures(model_fit, c("chisq", "df", "pvalue", "aic", "bic", "rmsea", "srmr"));
    model_fit_summary  <- standardizedSolution(model_fit, type="std.all", level=0.95);
    return (data.frame(
      # Phenotype info:
      "phenotype_number"    =i, 
      "phenotype_name"      =stringr::str_sub(varnames[i], end=-4), 
      "phenotype_name_short"=phenotype_names$short.names[ toupper(stringr::str_sub(phenotype_names$long.names.rater1, end=-3)) == toupper(stringr::str_sub(varnames[i], end=-4)) ],
      
      # Which genetic model is this?
      "genetic_model"=genetic_model, 
      "CD_meaning"=ifelse(genetic_model == "ACE", "C", "D"), 
      "H2_meaning"=ifelse(genetic_model == "ACE", "familiality", "broad sense heritability"),
      "deltaAIC_ACE_ADE"=deltaAIC_ACE_ADE,
      
      # Fit measures:
      "fit_Chisq"=model_fit_measures["chisq"], "fit_Chisq.df"=model_fit_measures["df"], "fit_Chisq.p"=model_fit_measures["pvalue"], 
      "fit_AIC"=model_fit_measures["aic"], "fit_BIC"=model_fit_measures["bic"], "fit_RMSEA"=model_fit_measures["rmsea"], "fit_SRMR"=model_fit_measures["srmr"], 
      
      # Relevant parameter estimates and 95%CIs:
      # a^2 (or h^2):
      "a2"=model_fit_summary$est.std[ model_fit_summary$label == "a" ][1]^2, 
      "a2.p"=model_fit_summary$pvalue[ model_fit_summary$label == "a" ][1],
      "a2.CiL"=model_fit_summary$ci.lower[ model_fit_summary$label == "a" ][1]^2,
      "a2.CiH"=model_fit_summary$ci.upper[ model_fit_summary$label == "a" ][1]^2,
      "a2.deltaAIC"=deltaAIC_A, "a2.LRT.Chisq"=LRtest_A[2,"Chisq diff"], "a2.LRT.Chisq.df"=LRtest_A[2,"Df diff"], "a2.LRT.Chisq.p"=LRtest_A[2,"Pr(>Chisq)"], # model comparison
      # cd^2:
      "cd2"=model_fit_summary$est.std[ model_fit_summary$label == ifelse(genetic_model=="ACE","c","d") ][1]^2, 
      "cd2.p"=model_fit_summary$pvalue[ model_fit_summary$label == ifelse(genetic_model=="ACE","c","d") ][1],
      "cd2.CiL"=model_fit_summary$ci.lower[ model_fit_summary$label == ifelse(genetic_model=="ACE","c","d") ][1]^2,
      "cd2.CiH"=model_fit_summary$ci.upper[ model_fit_summary$label == ifelse(genetic_model=="ACE","c","d") ][1]^2,
      "cd2.deltaAIC"=deltaAIC_other, "cd2.LRT.Chisq"=LRtest_other[2,"Chisq diff"], "cd2.LRT.Chisq.df"=LRtest_other[2,"Df diff"], "cd2.LRT.Chisq.p"=LRtest_other[2,"Pr(>Chisq)"], # model comparison
      # e^2:
      "e2"=model_fit_summary$est.std[ model_fit_summary$label == "e2" ][1], 
      "e2.p"=model_fit_summary$pvalue[ model_fit_summary$label == "e2" ][1],
      "e2.CiL"=model_fit_summary$ci.lower[ model_fit_summary$label == "e2" ][1],
      "e2.CiH"=model_fit_summary$ci.upper[ model_fit_summary$label == "e2" ][1]
    ));
  }
results_df <- do.call(rbind, results);

## Save results to file ####
write.csv(results_df, "../data/intermediate/lavaan_results.csv", row.names=FALSE);


# Save the explanation of the results:
cat(paste0('
# This document explains the content of the `lavaan_results.csv` file.\n
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
- **genetic\\_model**, **CD\\_meaning**, and **H2\\_meaning**: given that we cannot estimate simultaneously the *C* and *D* components, we must chose between the two using model compaison with AIC (given by **deltaAIC_ACE_ADE** = AIC(ACE) - AIC(ADE)); this, in turn, determines the fitted genetic model (*ACE* or *ADE*), the meaning of the "CD" (short for "C or D") parameter (*C* or *D*), and the meaning of the "H2" (actually *H*<sup>2</sup>) estimate ("familiality" for *C* and "broad sense heritability" for *D*, respectively);\n
- The fit measures of the full model (ACE or ADE):\n
  + **fit_Chisq**, **fit_Chisq.df** and **fit_Chisq.p**: the Chi squared test fit; **fit_AIC** and **fit_BIC** are AIC and BIC of the model, and **fit_RMSEA** the model\'s RMSEA\n
- The latent (standardized) estimates of relevant parameters, i.e., the standardized variance terms corrected for rater effects; for all, the point estimate, the lower and upper 95% CIs, and the *p*-value; for *a* and *cd*, we also give the model compairson with the model without the component as **deltaAIC** = AIC(model with component) - AIC(model without component) and as LR test\'s Chi squared test **LRT.Chisq**, **LRT.Chisq.df** and **LRT.Chisq.p**:\n
  + **a2**: the narrow heritability *h*<sup>2</sup>\n
  + **cd2**: *C*<sup>2</sup> of *D*<sup>2</sup>, respectively\n
  + **e2**: *E*<sup>2</sup>\n
\n
      ', 
           "Script took ", difftime(Sys.time(), run_start_time, units="mins"), " minutes to run...\n"), 
    file="../data/intermediate/lavaan_results_explanations.md", append=FALSE);



## END ####
cat("Took ", difftime(Sys.time(), run_start_time, units="mins"), " minutes...\n");

