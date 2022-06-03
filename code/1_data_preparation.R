# Data cleaning and preparation for the VT-HERITABILITY project
#
# Copyright (C) 2018-2021, Emily Jennings 
# checked and modified by Dan Dediu, 2022
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


# Figures:
if( !dir.exists("../figures/") ) dir.create("../figures/", showWarnings=FALSE);

# Intermediate data:
if( !dir.exists("../data/intermediate/") ) dir.create("../data/intermediate/", showWarnings=FALSE);


## Libraries ####
library(car);
library(dplyr);
library(lsr);
library(effsize);
library(moments);
library(corrplot);
library(mets);


## Load data and adjust variables ####

# Load the input data:
raw.data <- read.csv("../data/input/measures.csv", stringsAsFactors=FALSE);

# Relabelling sex with male=1 and female=2; from now on use "Sex" to obtain numeric values as code for the sexes:
raw.data$Sex <- factor(raw.data$sex, levels = c(1,2), labels = c("male", "female"));
# Relabelling twzyg with MZM=1, DZM=2, MZF=3, DZF=4, DOSmf=5, DOSfm=6:
raw.data$TWZYG <- factor(raw.data$twzyg, levels=c(1,2,3,4,5,6), labels=c("MZM", "DZM", "MZF", "DZF", "DOSmf", "DOSfm"));
# Relabelling Study with ADHD=1, OCD=2, Depression=3, EMIF=6, OBES=7:
raw.data$STUDY <- factor(raw.data$Study, levels=c(1,2,3,6,7), labels=c("ADHD", "OCS", "Depression", "Aging", "Obesity"));
# Recode twyzyg (3=1) (1=1) (2=2) (4 through 6=2) into zygmzdz:
raw.data$zygmzdz <- car::recode(raw.data$twzyg, "1=1; 3=1; 2=2; 4:6=2");
# Adding a column so that zygmzdz has a character label:
raw.data$ZYGMZDZ <- factor(raw.data$zygmzdz, levels = c(1,2), labels = c("MZ", "DZ"));
# Standardization of age at time of MRI:
raw.data$z_Age_at_MRI <- scale(raw.data$Age_MRI_1);
# Square standardized age:
raw.data$z_Age_at_MRI_SQR <- ((raw.data$z_Age_at_MRI)^2);
# str(raw.data[374:423]); # Check that everything is fine....


## MZ/DZ status ####

# Number of cases with no known MZ/DZ status:
nrow(raw.data[which(is.na(raw.data$ZYGMZDZ=='MZ')),]); # 39
# Select cases with no known MZ/DZ status:
noMZDZstatus <- raw.data[which(is.na(raw.data$ZYGMZDZ)),];
# Number of cases with no known MZ/DZ status that are parents:
nrow(raw.data[which(noMZDZstatus$bioprnt==1),]); #12
# Number of cases with no known MZ/DZ status that are siblings:
nrow(raw.data[which(noMZDZstatus$biosib==1),]); # 13
# Number of cases with no data on MZ/DZ status available (so twins with no MZ/DZ status):
nobioprntstatus <- noMZDZstatus[which(noMZDZstatus$bioprnt==-1),]; # not parents
nobiosibstatus <- nobioprntstatus[which(nobioprntstatus$biosib==-1),]; # also not sibs
# Of those that are not sibs or parents...
nrow(nobiosibstatus[which(nobiosibstatus$multiple_type==2),]); # 5 cases belong to twin pairs with no MZ/DZ status (and there are two sets of twin pairs and 1 twin without its co-twin)
nrow(nobiosibstatus[which(nobiosibstatus$multiple_type==0),]); # 6 cases are set as not belonging to any multiple set
nrow(nobiosibstatus[which(nobiosibstatus$multiple_type==3),]); # 3 cases are triplets (two from the same family)

# Remove cases for which MZ/DZ status is unknown:
raw.data <- raw.data[which(complete.cases(raw.data$ZYGMZDZ)),];

# Number of non-duplicate cases (is 632, so there are 6 duplicate cases (638-632=6):
length(unique(raw.data$FISNumber)); # 632


## Duplicate cases ####

# Extract duplicate cases FIS numbers:
duplicates <- raw.data$FISNumber[duplicated(raw.data$FISNumber)];
# View the duplicates cases:
#raw.data[which(raw.data$FISNumber==514351512444),]; # twin in OCD (54 years old) and EMIF (60.9 years old) studies; removed EMIF study data
#raw.data[which(raw.data$FISNumber==514351512566),]; # twin in OCD (54 years old) and EMIF (60.9 years old) studies (Dupe 1's twin); removed EMIF study data
#raw.data[which(raw.data$FISNumber==514351956244),]; # twin in OCD (56 years old) and EMIF (62.8 years old) studies; removed EMIF study data
#raw.data[which(raw.data$FISNumber==514351956966),]; # twin in OCD (56 years old) and EMIF (62.8 years old) studies (Dupe 3's twin); removed EMIF study data
#raw.data[which(raw.data$FISNumber==514366843656),]; # twin in OCD (56 years old) and EMIF (65.5 years old) studies; removed EMIF study data
#raw.data[which(raw.data$FISNumber==514366843834),]; # twin in OCD (56 years old) and EMIF (65.5 years old) studies (Dupe 5's twin); removed EMIF study data

# Remove duplicate cases:
raw.data <- distinct(raw.data, FISNumber, .keep_all= TRUE);

# Total number of cases after removing duplicates and cases for which MZ/DZ status was unknown:
nrow(raw.data); # 632


## ICV ####

# Number of cases with missing ICV:
nrow(raw.data[which(is.na(raw.data$ICV_new)),]); # 7
#View(raw.data[which(is.na(raw.data$ICV_new)),]);

# Method: average ICV for cases of the same sex and within +1 and -1 year of the age of the case with missing ICV:
raw.data$ICV_new <- as.numeric(raw.data$ICV_new);

# Split the file by sex:
male <- raw.data[which(raw.data$Sex=="male"),];
nrow(male[which(is.na(male$ICV_new)),]); # 3 male cases with missing ICV
female <- raw.data[which(raw.data$Sex=="female"),];
nrow(female[which(is.na(female$ICV_new)),]); # 4 female cases with missing ICV

#View(male[which(is.na(male$ICV_new)),]);

males36ish <- male[which(male$Age_MRI_1>35),];
males36ish <- males36ish[which(males36ish$Age_MRI_1<37),];
raw.data$ICV_new[137] <- mean(males36ish$ICV_new, na.rm=TRUE);
#View(raw.data[which(raw.data$FamilyNumber==22260),])
raw.data$ICV_new[136]; # check against twin

males33ish <- male[which(male$Age_MRI_1>32),];
males33ish <- males33ish[which(males33ish$Age_MRI_1<34),];
raw.data$ICV_new[164] <- mean(males33ish$ICV_new, na.rm = TRUE);
#View(raw.data[which(raw.data$FamilyNumber==23653),]);
raw.data$ICV_new[165]; # check against twin

males22ish <- male[which(male$Age_MRI_1>21),];
males22ish <- males22ish[which(males22ish$Age_MRI_1<23),];
raw.data$ICV_new[604] <- mean(males22ish$ICV_new, na.rm = TRUE);
#View(raw.data[which(raw.data$FamilyNumber==27987),]);
raw.data$ICV_new[603]; # check against twin

#View(female[which(is.na(female$ICV_new)),]);

females24ish <- female[which(female$Age_MRI_1>23),];
females24ish <- females24ish[which(females24ish$Age_MRI_1<25),];
raw.data$ICV_new[114] <- mean(females24ish$ICV_new, na.rm = TRUE);
raw.data$ICV_new[115] <- mean(females24ish$ICV_new, na.rm = TRUE);
# can't compare with twin because both twins have missing ICV

females35ish <- female[which(female$Age_MRI_1>34),];
females35ish <- females35ish[which(females35ish$Age_MRI_1<36),];
raw.data$ICV_new[284] <- mean(females35ish$ICV_new, na.rm = TRUE);
#View(raw.data[which(raw.data$FamilyNumber==12921),]);
raw.data$ICV_new[283]; # check against twin

females93ish <- female[which(female$Age_MRI_1>90),];
females93ish <- females93ish[which(females93ish$Age_MRI_1<96),];
raw.data$ICV_new[480] <- mean(females93ish$ICV_new, na.rm = TRUE);
#View(raw.data[which(raw.data$FamilyNumber==41435),]); # no twin available to check against

nrow(raw.data[which(is.na(raw.data$ICV_new)),]); # 0

# Standardization of ICV:
raw.data$z_ICV_new <- scale(raw.data$ICV_new);



## Select phenotype measures ####

# Split by landmarkers 1 and 2; excludes initial PCAs, but not Procrustes:
phenotypic_measures <- raw.data  %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance.1":"pharynx_PWTP_midsagittalTrace_matlabProcDistance.2"),
                                            c("hardPalate_HNSL_nsa_nsp_distance.1":"hardPalate_HMSP_midsagittalTrace_matlabProcDistance.2"),
                                            c("hardPalate_HCCP_coronalCanineTrace_matlabProcDistance.1":"hardPalate_HCCP_coronalCanineTrace_matlabProcDistance.2"),
                                            c("hardPalate_HCPP_coronalPM2Trace_matlabProcDistance.1":"hardPalate_HCPP_coronalPM2Trace_matlabProcDistance.2"),
                                            c("hardPalate_HCMP_coronalM2Trace_matlabProcDistance.1":"hardPalate_HCMP_coronalM2Trace_matlabProcDistance.2"),
                                            c("softPalate_SPUN_uvula_nsp_distance.1":"dentition_DDAP_dentalArchTrace_matlabProcDistance.2"),
                                            c("denition_DAIS_incisionSupApexSupLine_nsaNspLine_angle.1":"general_GRPD_registrationPointSet_matlabProcDistance.2"),
                                            c("general_GAPD_allPointSet_matlabProcDistance.1":"general_GAPD_allPointSet_matlabProcDistance.2"),
                                            "FISNumber");
#names(phenotypic_measures);

genetics_info <- raw.data %>% select(c("FISNumber":"z_ICV_new"));
#names(genetics_info);

# Read the phenotype info (name, description...):
phenotype_metainfo <- read.csv("../data/input/phenotype_names_decriptions.csv", stringsAsFactors=FALSE);

# Replace all phenotypic measures' outliers with NAs:
phenotypic_measures_noFIS <- phenotypic_measures %>% select(-"FISNumber");

# Exclude all outliers column-wise (i.e., any point farther than 3 SD from the column mean):
new_phenotypic_measures <- phenotypic_measures_noFIS;
for( i in 1:ncol(phenotypic_measures_noFIS) )
{
  new_phenotypic_measures[,i] <- (function(x) ifelse(!is.na(x) & (abs(x - mean(x,na.rm=TRUE)) <= 3*sd(x,na.rm=TRUE)), x, NA))(new_phenotypic_measures[,i]);
}
nrow(new_phenotypic_measures); # 632
new_phenotypic_measures.list <- as.list(new_phenotypic_measures);
#View(new_phenotypic_measures.list); # shows that no cases were thrown out (there are still 632 cases per PM), but just that the outliers were replaced with NAs
new_phenotypic_measures$FISNumber <- phenotypic_measures$FISNumber;

# Reorder columns so that FISNumber is the first one:
measures <- new_phenotypic_measures[,c(ncol(phenotypic_measures), 1:(ncol(phenotypic_measures)-1))];

# Merge the measures and the geentic info:
combo <- merge(measures, genetics_info, by="FISNumber");


## Descriptive stats ####
# (skip to Reshaping data section if descriptive statistics have already been done)

## Zygosity ####
# Number of MZ cases:
nrow(raw.data[which(raw.data$ZYGMZDZ=='MZ'),]);
# Number of DZ cases:
nrow(raw.data[which(raw.data$ZYGMZDZ=='DZ'),]);
# Number of DZ same-sex cases:
(n_DZ_samesex <- nrow(raw.data[which(raw.data$TWZYG=='DZM'),]) + nrow(raw.data[which(raw.data$TWZYG=='DZF'),]));
# Number of DZ DOS cases:
(n_DOS <- nrow(raw.data[which(raw.data$TWZYG=='DOSfm'),]) + nrow(raw.data[which(raw.data$TWZYG=='DOSmf'),]));

# table for zygosity number (n) before grouped by MZ or DZ:
(twzyg_bar <- table(raw.data$TWZYG));
#t able for twin zygosity number (n) when split by MZ or DZ:
(zyg_bar <- table(raw.data$ZYGMZDZ));


# Studies ####
# Studies from which MRI data was obtained:
(studies_bar <- table(raw.data$STUDY));

# Number of cases per study:
ADHDcases <- raw.data[which(raw.data$STUDY=="ADHD"),];
OCDcases <- raw.data[which(raw.data$STUDY=="OCS"),];
Depressioncases <- raw.data[which(raw.data$STUDY=="Depression"),];
EMIFcases <- raw.data[which(raw.data$STUDY=="Aging"),];
OBEScases <- raw.data[which(raw.data$STUDY=="Obesity"),];
# Check that they add up:
(nrow_total <- nrow(ADHDcases) + nrow(OCDcases) + nrow(Depressioncases) + nrow(EMIFcases) + nrow(OBEScases));

# Table of age range per study:
studysummaries <- rbind(summary(ADHDcases$Age_MRI_1), summary(OCDcases$Age_MRI_1), summary(Depressioncases$Age_MRI_1), summary(EMIFcases$Age_MRI_1), summary(OBEScases$Age_MRI_1));
studysummaries.df <- as.data.frame(studysummaries);
studysummaries.df$`1st Qu.` <- NULL;
studysummaries.df$Median <- NULL;
studysummaries.df$`3rd Qu.` <- NULL;
studysummaries.df$SD <- rbind(sd(ADHDcases$Age_MRI_1, na.rm=TRUE), sd(OCDcases$Age_MRI_1, na.rm=TRUE), sd(Depressioncases$Age_MRI_1, na.rm=TRUE), sd(EMIFcases$Age_MRI_1, na.rm=TRUE), sd(OBEScases$Age_MRI_1, na.rm=TRUE));
studysummaries.df$n <- rbind(nrow(ADHDcases), nrow(OCDcases), nrow(Depressioncases), nrow(EMIFcases), nrow(OBEScases));
rownames(studysummaries.df) <- c("ADHD", "OCS", "Depression", "Aging", "Obesity");
colnames(studysummaries.df) <- c("Min. age", "mean age", "Max. age", "SD", "N");
studysummaries.df$`mean age` <- format(round(studysummaries.df$`mean age`, digits=2), nsmall=2);
studysummaries.df$SD <- format(round(studysummaries.df$SD, digits=2), nsmall=2);
studysummaries.df; #View(studysummaries.df);


## Sex-specific stats ####

# Split the file by sex:
males <- raw.data[which(raw.data$Sex=="male"),];
females <- raw.data[which(raw.data$Sex=="female"),];
# Number of males:
nrow(males); # 223
# Number of females
nrow(females); # 409

# Combining male and female boxplots into one image:
attach(raw.data);
# cex.axis determines the size of the "Males" and "Females" labels.
# Make this mfrow=c(2,4) when the PRSs are done because there will be 8 panels then rather than 4
slices <- c(223, 409);
lbls <- c("males", "females");
pct <- round(slices/sum(slices)*100);
lbls <- paste(lbls, pct); # add percents to labels 
lbls <- paste(lbls,"%",sep=""); # add % to labels 
if( FALSE ) # don't actually plot it!
{
  tiff(file="../figures/descriptives_1.tiff", width=2*4, height=2*4, units="in", res=600, compression="lzw");
  par(mar=c(5,5,5,3)+0.1,mfrow=c(2,2));
  sex.name <- c("Males", "Females");
  pie(slices,labels = lbls, col=(c("light blue", "pink")), main="Percentage of \nM & F cases", cex.main=1.5);
  mtext("A", side=3, line=2, adj=-0.3, cex=1.2);
  boxplot(ICV_new~Sex, data=raw.data, 
          notch=FALSE, col=(c("light blue", "pink")), 
          main="M & F \nIntracranial volume (ICV)", cex.main=1.5, 
          xlab="Sex", ylab=expression('ICV value (mm'^3*')'), cex.lab=1.3, cex.axis=1.3, xaxt='n');
  axis(1, at=1:2, labels=sex.name, cex.axis=1.2, lwd.ticks=0);
  mtext("B", side=3, line=2, adj=-0.3, cex=1.2);
  boxplot(Age_MRI_1~Sex, data=raw.data, 
          notch=FALSE, col=(c("light blue", "pink")), 
          main="M & F Age at MRI", cex.main=1.5, 
          xlab="Sex", ylab="Age", cex.lab=1.3, cex.axis=1.3, xaxt='n');
  axis(1, at=1:2, labels=sex.name, cex.axis=1.2, lwd.ticks=0);
  mtext("C", side=3, line=2, adj=-0.3, cex=1.2);
  dev.off();
}

# Descriptives table sorted by males and females:
descriptives_M <- rbind(summary(males$Age_MRI_1), summary(males$ICV_new));
descriptives_F <- rbind(summary(females$Age_MRI_1), summary(females$ICV_new));
descriptives_M.df <- as.data.frame(descriptives_M);
descriptives_F.df <- as.data.frame(descriptives_F);
rownames(descriptives_M.df) <- c("Age", "ICV (mm^3)");
rownames(descriptives_F.df) <- c("Age ", "ICV (mm^3) ");
descriptives_M.df$SD <- rbind(sd(males$Age_MRI_1, na.rm=TRUE), sd(males$ICV_new, na.rm=TRUE));
descriptives_F.df$SD <- rbind(sd(females$Age_MRI_1, na.rm=TRUE), sd(females$ICV_new, na.rm=TRUE));
descriptives_1.df <- rbind(descriptives_M.df,descriptives_F.df);
colnames(descriptives_1.df) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "SD");
descriptives_1.df$Mean <- format(round(descriptives_1.df$Mean, digits=2), nsmall=2);
descriptives_1.df$SD <- format(round(descriptives_1.df$SD, digits=2), nsmall=2);

# Welch's t-test comparing males and females for age and ICV:
(ttest_ICV <- t.test(males$ICV_new, females$ICV_new));
(ttest_Age_MRI_1 <- t.test(males$Age_MRI_1, females$Age_MRI_1));
# Start combining into one table:
tvalues_extr <- rbind(ttest_ICV$statistic, ttest_Age_MRI_1$statistic);
inf_table.df <- as.data.frame(tvalues_extr);
# Add df to table:
inf_table.df$df <- rbind(ttest_ICV$parameter, ttest_Age_MRI_1$parameter);
#A dd CI to table:
inf_table.df$CI <- rbind(ttest_ICV$conf.int, ttest_Age_MRI_1$conf.int);
# Add p-values to table:
inf_table.df$pvalues <- rbind(ttest_ICV$p.value, ttest_Age_MRI_1$p.value);
# Add effect size to the table:
inf_table.df$CohensD <- rbind((cohensD(males$ICV_new, females$ICV_new, method = "unequal")), (cohensD(males$Age_MRI_1, females$Age_MRI_1, method = "unequal")));
# Round all numbers:
inf_table.df$t <- format(round(inf_table.df$t, digits=2), nsmall=2);
inf_table.df$df <- format(round(inf_table.df$df, digits=2), nsmall=2);
inf_table.df$CI <- format(round(inf_table.df$CI, digits=2), nsmall=2);
inf_table.df$CohensD <- format(round(inf_table.df$CohensD, digits=2), nsmall=2);
# Convert the dataframe to matrix:
inferentials_table <- as.matrix(inf_table.df);
rownames(inferentials_table) <- c("ICV", "Age at MRI");
colnames(inferentials_table) <- c("t", "df", "CI lower", "CI upper", "p", "Cohen's d");
inferentials_table; #View(inferentials_table)


## Phenotypic measures sorted by domain and type ####

names_list <- mat.or.vec(16,1);
names_list.df <- as.data.frame(names_list);
colnames(names_list.df) <- "% of phenotypic measures";
rownames(names_list.df) <- c("cervical", "dentition", "general", "hard palate", "hyoid", "larynx", "mandible", "oral", "pharynx", "skull", "soft palate", "angle", "curvature", "distance", "ratio", "Procrustes dist.");
names_list.df[,1] <- c(
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$domain=="cervical"),])/nrow(phenotype_metainfo))*100),
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$domain=="dentition"),])/nrow(phenotype_metainfo))*100), 
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$domain=="general"),])/nrow(phenotype_metainfo))*100), 
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$domain=="hard palate"),])/nrow(phenotype_metainfo))*100), 
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$domain=="hyoid"),])/nrow(phenotype_metainfo))*100), 
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$domain=="larynx"),])/nrow(phenotype_metainfo))*100), 
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$domain=="mandible"),])/nrow(phenotype_metainfo))*100), 
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$domain=="oral"),])/nrow(phenotype_metainfo))*100), 
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$domain=="pharynx"),])/nrow(phenotype_metainfo))*100), 
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$domain=="skull"),])/nrow(phenotype_metainfo))*100), 
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$domain=="soft palate"),])/nrow(phenotype_metainfo))*100),
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$type=="angle"),])/nrow(phenotype_metainfo))*100),
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$type=="curvature"),])/nrow(phenotype_metainfo))*100),
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$type=="distance"),])/nrow(phenotype_metainfo))*100),
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$type=="ratio"),])/nrow(phenotype_metainfo))*100),
  ((nrow(phenotype_metainfo[which(phenotype_metainfo$type=="Procrustes dist."),])/nrow(phenotype_metainfo))*100));
names_list.df$`% of phenotypic measures` <- format(round(names_list.df$`% of phenotypic measures`, digits=2), nsmall=2);



## Detecting outliers ####
if( FALSE ) ## DON'T RUN!!!
{
  # Assumptions of linear regression to test:
  # 1) The Y-values (or errors, "e" (our residual values, the difference between the observed and predicted y-values))
  #    are independent --> validity of this assumption is determined by knowledge of the study design or data collection
  # Check the following by examining the model residuals or errors
  # 2) The Y-values can be expressed as a linear function of the X variable
  # - if plotting the residuals (y-axis) vs the fitted values (lm() for any x variable as a predictor; x-axis),
  #   then the red line should be fairly flat if the linearity assumption is met
  # 3) Variation of observations around the regression line (the residual SE (aka the standard deviation of the residuals)) is constant (aka homoscedastic)
  # - if the variation is constant, then there should also be no pattern of the dots plotted in the residuals vs. fitted plot: the dots
  #   should look like a cloud of dots
  # 4) For a given value of X, Y values (or the error) are normally distributed
  # - use a Q-Q plot to check normality of the Y values: y-axis = standardized residuals, x-axis = theoretical quantiles (i.e. lm(CS2H ~ age))
  #   if the residuals are normally distributed, then they should fall roughly along a diagonal line
  
  # Note that the following have been done on both the averages (using dataframe "co" (which is generated in the section below in "Identifying Interaction Terms"))
  # and on the data for each rater (saved to a separate folder)
  
  # Check for linearity, homoscedasticity, and normality of each phenotypic phenotypic measure with age, age^2, ICV, and Sex:
  a <- colnames(measures); # when done per rater
  
  # Residual vs fitted and Q-Q plots with all measures and all covariates
  lm_rater1 <- lapply(combo[,c(2:301)], function(x) lm(x ~ z_Age_at_MRI + z_Age_at_MRI_SQR + z_ICV_new + Sex, data=combo));
  pdf("../figures/outlier_detection__phens_vs_IVs__fitted_vs_residuals.pdf");
  for(d in 1:length(measures)) plot(lm_rater1[[d]], main=paste0(a[[d]])); # one plot per page
  dev.off();
  
  # z age:
  lm_age <- lapply(combo[,c(2:301)], function(x) lm(x ~ z_Age_at_MRI, data=combo));
  pdf("../figures/outlier_detection__phens_vs_age__fitted_vs_residuals.pdf")
  for(d in 1:length(measures)) plot(lm_age[[d]], main=paste0(a[[d]]));
  dev.off();
  # some pharynx measures, hard palate angles, and hard palate, dentition, and 
  # general Procrustes distances were rather heavily tailed and non-linear, but  
  # most everything else met the assumptions
  
  # (z age)^2:
  lm_age2 <- lapply(combo[,c(2:301)], function(x) lm(x ~ z_Age_at_MRI_SQR, data=combo));
  pdf("../figures/outlier_detection__phens_vs_age2__fitted_vs_residuals.pdf");
  for(d in 1:length(measures)) plot(lm_age2[[d]], main=paste0(a[[d]]));
  dev.off();
  # most are not linearly related, as makes sense because the predictor is the square of a variable that was linearly related;
  # most are normally distributed; the pharynx wall measures were more heavily tailed than all of the other measures
  
  # ICV:
  lm_ICV <- lapply(combo[,c(2:301)], function(x) lm(x ~ z_ICV_new, data=combo));
  pdf("../figures/outlier_detection__phens_vs_ICV__fitted_vs_residuals.pdf");
  for(d in 1:length(measures)) plot(lm_ICV[[d]], main=paste0(a[[d]]));
  dev.off();
  # all phenotypes normally distributed except some pharynx wall measures, hard palate angle measures, 
  # hard palate and dentition Procrustes distance measures, and general distance measures
  
  # Sex:
  lm_Sex <- lapply(combo[,c(2:301)], function(x) lm(x ~ Sex, data=combo));
  pdf("../figures/outlier_detection__phens_vs_Sex__fitted_vs_residuals.pdf");
  for(d in 1:length(measures)) plot(lm_Sex[[d]], main=paste0(a[[d]]));
  dev.off();
  # once again, pharynx wall measures and HDMP were heavily tailed, but the rest were normally distributed; all measure are linear, 
  # but the distinction between males and females is very clear 
  
  # Safe to delete:
  cat("The outlier_detection__*.pdf files are safe to delete if needed.\n", file="../figures/outlier_detection__are_safe_to_delete.txt", append=FALSE);
}


## Interaction terms & stepwise regressions ####

# The covariates are: ICV_new, Age_MRI_1 and Sex:
covariates <- combo %>% select(c("z_Age_at_MRI", "z_Age_at_MRI_SQR","z_ICV_new"));

# look for interactions between ICV_new (x-axis) * Sex (likely to exist because ICV is significantly larger in males)
#                               Age_MRI_1 * Sex (likely to exist because women live longer than men in the Netherlands, where our data are from)
# don't look for interactions between ICV_new (x-axis) * Age_MRI_1 because it's unlikely to occur due to ICV plateauing at 10 years of age, which is an age before our sample age minimum

# Correlations among covariates that are continuous variables:
# age is not normally distributed (see q-q plots below), but ICV is, so use a non-parametric test?
# no, stick with parametric tests because the sample we have is large (n = 632), the Dutch population is pretty normally 
# distributed when it comes to age (https://www.statista.com/statistics/276710/age-distribution-in-the-netherlands/)
# and the equal environment assumption in the classical twin model does not necessarily hold anyway...
# and Pearson's correlation does not assume normality anyway https://stats.stackexchange.com/questions/3730/pearsons-or-spearmans-correlation-with-non-normal-data
# https://en.wikipedia.org/wiki/Correlation_and_dependence

qqnorm(combo$z_ICV_new); qqline(combo$z_ICV_new, col="red"); # ICV is actually fairly normally distributed
qqnorm(combo$z_Age_at_MRI); qqline(combo$z_Age_at_MRI, col="red"); # age is does not look normally distributed...
qqnorm(combo$z_Age_at_MRI_SQR); qqline(combo$z_Age_at_MRI_SQR, col="red"); # ...so I don't expect age^2 to be either

# Skewness:
skewness(combo$z_ICV_new) # distribution is approximately symmetric
skewness(combo$z_Age_at_MRI) # distribution is approximately symmetric still, even though it doesn't look normal
skewness(combo$z_Age_at_MRI_SQR) # distribution is highly skewed, but that's expected when it's the square of another distribution

# ICV by Age:
plot(combo$z_Age_at_MRI, combo$z_ICV_new);

# Correlations:
(S_correlations <- cor(covariates, covariates, method="kendall"));
(P_correlations <- cor(covariates, covariates, method="pearson"));
P_pvalues <- cor.mtest(covariates, method="pearson");
P_pvalues_mat <- as.matrix(P_pvalues$p);
p.adjust((P_pvalues_mat), method = "holm");
# the continuous covariates are not highly correlated (no matter if using Pearson or Spearman's correlation) and so they can be included in the multiple regression equation together

# The continuous variable age is not normally distributed (see Q-Q plots above), so should use a non-parametric test to evaluate the relationships 
## between the continuous and categorical variables (Wilcoxon-Mann-Whitney U test); however, we'll use the parametric version (Welch's t-test) 
## because we are using all parametric tests for the reasons described above

# Wilcoxon-Mann-Whitney U test for whether or not ICV differs between males and females in our population:
# MW_zICV_sex <- wilcox.test(z_ICV_new ~ sex, data = combo)
# WMW for whether or not age differs between males and females in our population:
# MW_zAge_sex <- wilcox.test(z_Age_at_MRI ~ sex, data = combo)
# WMW for whether or not age^2 differs between males and females in our population: 
# MW_zAgesqr_sex <- wilcox.test(z_Age_at_MRI_SQR ~ sex, data = combo) # but it really shouldn't because standardised age itself didn't....  

# the t-tests:
# Welch's t-test comparing males and females for standardized age and standardized ICV:
(ttest_zICV <- t.test(males$z_ICV_new, females$z_ICV_new));
(ttest_zAge_MRI_1 <- t.test(males$z_Age_at_MRI, females$z_Age_at_MRI));
# Combine into one table:
ztvalues_extr <- rbind(ttest_zICV$statistic, ttest_zAge_MRI_1$statistic);
z_inf_table.df <- as.data.frame(ztvalues_extr);
z_inf_table.df$df <- rbind(ttest_zICV$parameter, ttest_zAge_MRI_1$parameter);
z_inf_table.df$CI_lower <- rbind(ttest_zICV$conf.int[[1]], ttest_zAge_MRI_1$conf.int[[1]]);
z_inf_table.df$CI_upper <- rbind(ttest_zICV$conf.int[[2]], ttest_zAge_MRI_1$conf.int[[2]]);
z_inf_table.df$pvalues <- rbind(ttest_ICV$p.value, ttest_Age_MRI_1$p.value);
z_inf_table.df$CohensD <- rbind((cohen.d(males$z_ICV_new, females$z_ICV_new)), (cohen.d(males$z_Age_at_MRI, females$z_Age_at_MRI)));
z_inf_table.df[1,7] <- as.numeric(z_inf_table.df$CohensD[1,3]);
z_inf_table.df[2,7] <- as.numeric(z_inf_table.df$CohensD[2,3]);
colnames(z_inf_table.df)[7] <- "cohensD";
z_inf_table.df$CohensD <- NULL;
z_inferentials_table <- as.matrix(z_inf_table.df);
rownames(z_inferentials_table) <- c("z ICV", "z Age at MRI");
colnames(z_inferentials_table) <- c("t", "df", "CI lower", "CI upper", "p","Cohens D");
z_inferentials_table.df <- as.data.frame(z_inferentials_table);
z_inferentials_table.df; #View(z_inferentials_table.df)

# Adjust p-values:
tt_pvals <- mat.or.vec(1,2);
names(tt_pvals) <- c("z ICV", "z age");
tt_pvals[,1] <- z_inferentials_table.df[1,"p"];
tt_pvals[,2] <- z_inferentials_table.df[2,"p"];
(tt_p <- p.adjust(tt_pvals, method = "holm"));
# Like above with the non-standardized data, there is a significant difference in ICV between the sexes, but not a significant difference in age between the sexes


## VIF including the interaction terms ####

# Do the following to streamline the stepwise regressions:
# make a dataframe for non-PM info: already exists as "genetics"
# create averages for PMs out of rater's scores:
measures$FISNumber <- NULL;
# Measures created by rater 1:
phenotypic_measures_rater1 <- measures[ , c(TRUE,FALSE)];
names(phenotypic_measures_rater1);
# Measures created by rater 2:
phenotypic_measures_rater2 <- measures[ ,c(FALSE,TRUE)];
names(phenotypic_measures_rater2);
# Checks:
all( substring(names(phenotypic_measures_rater1), nchar(names(phenotypic_measures_rater1))-1) == ".1" );
all( substring(names(phenotypic_measures_rater2), nchar(names(phenotypic_measures_rater2))-1) == ".2" );
nrow(phenotypic_measures_rater1) == nrow(phenotypic_measures_rater2);
n <- ncol(phenotypic_measures_rater1);

PMs_byRater <- cbind(phenotypic_measures_rater1, phenotypic_measures_rater2);

PMs_averages <- mat.or.vec(632,n);
PMs_averages.df <- as.data.frame(PMs_averages);
# NAs will be kept for instances where one rater does not have a value supplied
for (i in 1:n)
{
  PMs_averages.df[[i]] <- (PMs_byRater[[i]] + PMs_byRater[[i+n]])/2;
}
colnames(PMs_averages.df) <- substring(names(phenotypic_measures_rater1), 1, nchar(names(phenotypic_measures_rater1))-2);
#nrow(PMs_averages.df[which(is.na(PMs_averages.df$cervicalSpine_CS2H_atlas_c2base_distance)),])

# Combine averages with genetics info:
co <- cbind(PMs_averages.df, genetics_info);


# VIF analysis:
myvarc <- colnames(co %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance":"general_GAPD_allPointSet_matlabProcDistance")));  # get the names of the phenotypic variables
myvars <- as.list(myvarc) # make the names into a list

# Compare changes in VIF values as predictors are individually removed from the equations:
dvar <- co %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance":"general_GAPD_allPointSet_matlabProcDistance")) # data for all phenotypic measures

# Start will all meaningful IVs and interactions:
vifs <- do.call(rbind, lapply(names(dvar), function(v) vif(lm(as.formula(paste0(v, "~ z_Age_at_MRI + z_Age_at_MRI_SQR + Sex + z_ICV_new + (z_Age_at_MRI:Sex) + (z_Age_at_MRI_SQR:Sex) + (z_ICV_new:Sex)")), data=co)))); 
rownames(vifs) <- names(dvar); 
vifs; #View(vifs);
sum(vifs > 5.0); # z_Age_at_MRI_SQR:Sex has VIFs > 0 -> remove it!

vifs <- do.call(rbind, lapply(names(dvar), function(v) vif(lm(as.formula(paste0(v, "~ z_Age_at_MRI + z_Age_at_MRI_SQR + Sex + z_ICV_new + (z_Age_at_MRI:Sex) + (z_ICV_new:Sex)")), data=co)))); 
rownames(vifs) <- names(dvar); 
vifs; #View(vifs);
sum(vifs > 5.0); # much better: none of the VIF values exceed 5!
# -> so, the predictors to keep are z_Age_at_MRI, z_Age_at_MRI_SQR, Sex, z_ICV_new, (z_Age_at_MRI:Sex) and (z_ICV_new:Sex)


## Stepwise regressions ####

library(MASS);
full.models <- lapply(names(dvar), function(v) lm(as.formula(paste0(v, "~ Sex + z_Age_at_MRI + z_Age_at_MRI_SQR + z_ICV_new + (z_Age_at_MRI:Sex) + (z_ICV_new:Sex)")), data = co));
names(full.models) <- names(dvar);
stepwise_models <- lapply(full.models, function(m) MASS::stepAIC(m, direction = "both", na.rm=TRUE));
detach("package:MASS"); # so that I can use select() from the dplyr package again

# Table of covariates kept in the final models per phenotypic PM in stepAIC():
All_covariates <- c("Sexfemale", "z_Age_at_MRI", "z_Age_at_MRI_SQR", "z_ICV_new", "Sexfemale:z_Age_at_MRI", "Sexfemale:z_ICV_new");
out <- lapply(names(dvar), function(x1) lapply(All_covariates, function(x2) as.numeric(x2 %in% names(eval(parse(text=paste("stepwise_models",x1,"coefficients",sep="$")))))));
finalStepCovariates <- matrix(unlist(out), ncol = 6, byrow = TRUE);
colnames(finalStepCovariates) <- c("Sex", "z age", "(z age)^2", "z ICV", "Sex:z_Age_at_MRI","Sex:zICV");
rownames(finalStepCovariates) <- names(dvar);
finalStepCovariates.df <- as.data.frame(finalStepCovariates);
# Patterns: number of VT measures for which each covariate was important:
# Sex
sum(finalStepCovariates.df$Sex) / nrow(finalStepCovariates.df) * 100; # 92%, 138
# Age
sum(finalStepCovariates.df$`z age`) / nrow(finalStepCovariates.df) * 100; # 81%, 12
# Age^2
sum(finalStepCovariates.df$`(z age)^2`) / nrow(finalStepCovariates.df) * 100; # 55%, 83
# ICV
sum(finalStepCovariates.df$`z ICV`) / nrow(finalStepCovariates.df) * 100; # 67%, 101
# Age:Sex
sum(finalStepCovariates.df$`Sex:z_Age_at_MRI`) / nrow(finalStepCovariates.df) * 100; # 31%, 46
# ICV:Sex
sum(finalStepCovariates.df$`Sex:zICV`) / nrow(finalStepCovariates.df) * 100; # 32%, 48
# -> Age and Sex come out as the covariates relevant for most DVs, and Age^2 and ICV as relevant for 1/2-1/3 of DVs...


# Extract the p-values from the stepwise regressions:
finalStepCovariates_pvalues <- mat.or.vec(ncol(dvar),6);
colnames(finalStepCovariates_pvalues) <- c("Sexfemale", "z_Age_at_MRI", "z_Age_at_MRI_SQR", "z_ICV_new", "Sexfemale:z_Age_at_MRI", "Sexfemale:z_ICV_new");
rownames(finalStepCovariates_pvalues) <- names(dvar);
# Table of covariates kept in the final model for per phenotypic PM in stepAIC(), by p-value:
out1 <- lapply(names(dvar), function(x1) lapply(All_covariates, function(x2) ifelse(x2 %in% names(eval(parse(text=paste("stepwise_models",x1,"coefficients",sep="$")))), 1, NA) )); 
names(out1) <- names(dvar);
out.1A <- out1;
out.1B <- lapply(names(dvar), function(x) parse(text=paste("out.1A$",x,"[which(!is.na(out1$",x,"),arr.ind=T)] <- summary(stepwise_models$",x,")$coefficients[2:length(stepwise_models$",x,"$coefficients),4]",sep="")))
for (i in 1:length(out.1B)) if (sum(!is.na(out1[[i]])) > 0) {eval(out.1B[[i]])};
finalStepCovariates_pvals_table <- matrix(unlist(out.1A), ncol = 6, byrow = TRUE);
colnames(finalStepCovariates_pvals_table) <- c("Sexfemale", "z_Age_at_MRI", "z_Age_at_MRI_SQR", "z_ICV_new", "Sexfemale:z_Age_at_MRI", "Sexfemale:z_ICV_new");
rownames(finalStepCovariates_pvals_table) <- names(dvar);
finalStepCovariates_pvals_table.df <- as.data.frame(finalStepCovariates_pvals_table);

# Check how many IV and interaction terms are needed:
finalStepCovariates_pvals_table.df$SexPgr.5    <- ifelse(finalStepCovariates_pvals_table.df$Sexfemale<=0.05, paste("yes"), paste("no"));
finalStepCovariates_pvals_table.df$AgePgr.5    <- ifelse(finalStepCovariates_pvals_table.df$z_Age_at_MRI<=0.05, paste("yes"), paste("no"));
finalStepCovariates_pvals_table.df$AgeSqrPgr.5 <- ifelse(finalStepCovariates_pvals_table.df$z_Age_at_MRI_SQR<=0.05, paste("yes"), paste("no"));
finalStepCovariates_pvals_table.df$ICVPgr.5    <- ifelse(finalStepCovariates_pvals_table.df$z_ICV_new<=0.05, paste("yes"), paste("no"));
finalStepCovariates_pvals_table.df$AgeSexPgr.5 <- ifelse(finalStepCovariates_pvals_table.df$`Sexfemale:z_Age_at_MRI`<=0.05, paste("yes"), paste("no"));
finalStepCovariates_pvals_table.df$SexICVPgr.5 <- ifelse(finalStepCovariates_pvals_table.df$`Sexfemale:z_ICV_new`<=0.05, paste("yes"), paste("no"));

# Sex:
sum(finalStepCovariates_pvals_table.df$SexPgr.5=='yes', na.rm=TRUE) / nrow(finalStepCovariates_pvals_table.df) * 100; # 81%, 121
# Age:
sum(finalStepCovariates_pvals_table.df$AgePgr.5=='yes', na.rm=TRUE) / nrow(finalStepCovariates_pvals_table.df) * 100; # 65%, 97
# Age^2:
sum(finalStepCovariates_pvals_table.df$AgeSqrPgr.5=='yes', na.rm=TRUE) / nrow(finalStepCovariates_pvals_table.df) * 100; # 47%, 71
# ICV:
sum(finalStepCovariates_pvals_table.df$ICVPgr.5=='yes', na.rm=TRUE) / nrow(finalStepCovariates_pvals_table.df) * 100; # 43%, 64
# Age:Sex:
sum(finalStepCovariates_pvals_table.df$AgeSexPgr.5=='yes', na.rm=TRUE) / nrow(finalStepCovariates_pvals_table.df) * 100; # 11%, 17
# ICV:Sex:
sum(finalStepCovariates_pvals_table.df$SexICVPgr.5=='yes', na.rm=TRUE) / nrow(finalStepCovariates_pvals_table.df) * 100; # 9%, 14


# Does including all terms in each VT's full model equation greatly make the equation less predictive?
# I need to know if the model fits better for each PM with only the final stepwise model IVs present versus all IVs present:
mcod <- mat.or.vec(length(stepwise_models), 2);
mcod.df <- as.data.frame(mcod);
colnames(mcod.df) <- c("mcod M1", "mcod M2");
for(i in (1:n)) mcod.df[i,1]<- summary(stepwise_models[[i]])$adj.r.squared;
for(i in (1:n)) mcod.df[i,2]<- summary(full.models[[i]])$adj.r.squared;
mcod.df; #View(mcod.df);
sum(mcod.df$`mcod M1` < mcod.df$`mcod M2`); mcod.df[ mcod.df$`mcod M1` < mcod.df$`mcod M2`, ];
# -> the multiple coefficients of determination (mcod's) are not really different (and not different at all for some VTs) when all terms are included...

# Thus, we will include Age and Sex, but also the interaction terms, as these are important for at least one PM...



## Reshape the data for GCSM ####

# Put the males in the first columns (so that SEM is done correctly):
# 1) reshape by sorting by FamilyNumber and multi_id_fam:
#    there is a family (FamilyNumber=14331) that has two sets of DZ twins, born 5 years apart, so use idcombine=TRUE to make sure that the two identifers are used together 
#    ("FamilyNumber" and "mult_id_fam" can prevent those twins from being grouped as a group of four under the same FamilyNumber):
twinwide <- fast.reshape(combo,
                         id=c("FamilyNumber", "mult_id_fam"), 
                         idcombine=TRUE, factor=TRUE, sep = "_", labelnum=FALSE, 
                         varying=c("FISNumber", as.character(colnames(combo[2:301])), "ICV_new", "birthorder", "sex","Sex","z_ICV_new","z_Age_at_MRI","Age_MRI_1","z_Age_at_MRI_SQR"));
#View(twinwide);

# 2) reorder columns:
twinwide_new <- twinwide %>% select(c("FISNumber_1":"general_GAPD_allPointSet_matlabProcDistance.2_1"),
                                    c("FISNumber_2":"general_GAPD_allPointSet_matlabProcDistance.2_2"),
                                    "Sex_1", "z_Age_at_MRI_1", "z_Age_at_MRI_SQR_1", "z_ICV_new_1",
                                    "Sex_2", "z_Age_at_MRI_2", "z_Age_at_MRI_SQR_2", "z_ICV_new_2",
                                    "Study", c("PedigreeNumber":"ZYGMZDZ"));
names(twinwide_new); # check that the reordering is correct

# 3) split up the dataframe into two dataframes (for same-sex MZ/DZ pairs and one for DOS pairs) 
#    so that DOS pairs with males in the second columns can have males moved to the first columns.
# Selecting same-sex cases with MZ or DZ status:
MZ_DZ_cases <- twinwide_new[which(twinwide_new$twzyg <= 4),];
# Selecting DOS cases:
DOS <- twinwide_new[which(twinwide_new$twzyg >= 5),];

DOS_males_Sex1 <- DOS[which(DOS$Sex_1=='male'),];
DOS_males_Sex2 <- DOS[which(DOS$Sex_2=='male'),]; # these are the ones I want to move to the first columns
# DOS single cases:
DOS_NAforSex2 <- DOS[is.na(DOS$Sex_2),];

# DOS single cases not already taken care of by DOS_males_sex1 (aka the ones where sex 1 = female):
DOS_NAforSex2_female <- DOS_NAforSex2[which(DOS_NAforSex2$Sex_1=='female'),];

# Check that there are still 62 DOS pairs: 
(total_n <- nrow(DOS_males_Sex1) + nrow(DOS_males_Sex2) + nrow(DOS_NAforSex2_female)); # 62

# Reorder columns so that males and single cases are in the first columns...
DOS_males_Sex2_reordered <- DOS_males_Sex2 %>% select(c("FISNumber_2":"general_GAPD_allPointSet_matlabProcDistance.2_2"),
                                                      c("FISNumber_1":"general_GAPD_allPointSet_matlabProcDistance.2_1"),
                                                      c("Sex_2":"z_ICV_new_2"),c("Sex_1":"z_ICV_new_1"),
                                                      c("Study":"ZYGMZDZ"));
colnames(DOS_males_Sex2_reordered)[1:301]   <- c(names(DOS_males_Sex2)[1:301])   # relabelled because now twin 2 is twin 1...
colnames(DOS_males_Sex2_reordered)[302:602] <- c(names(DOS_males_Sex2)[302:602]) # ....and twin 1 is twin 2, making all males twin 1...
colnames(DOS_males_Sex2_reordered)[603:606] <- c(names(DOS_males_Sex2)[603:606]) # ....and need to therefore also apply the twin 1 label...
colnames(DOS_males_Sex2_reordered)[607:610] <- c(names(DOS_males_Sex2)[607:610]) # ....and the twin 2 label to the appropriate columns

# Putting the reorderings all back together:
DOS_reordered <- rbind(DOS_males_Sex1,           # DOS males first
                       DOS_males_Sex2_reordered, # DOS originally with males second, but now males first as well
                       DOS_NAforSex2_female);    # DOS females

# 4) putting MZ/DZ cases and DOS cases back into the same dataframe:
recombined <- rbind(MZ_DZ_cases,    # MZ and same-sex DZ
                    DOS_reordered); # DOS twins (males first)
nrow(twinwide_new) == nrow(recombined); # same number of pairs+single cases as in twinwide_new

# Remove the genetics info:
twin1 <- recombined %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance.1_1":"general_GAPD_allPointSet_matlabProcDistance.2_1"));
twin2 <- recombined %>% select(c("cervicalSpine_CS2H_atlas_c2base_distance.1_2":"general_GAPD_allPointSet_matlabProcDistance.2_2"));

# TWIN 1:
# Rater 1:
twin1_rater1 <- twin1[ , c(TRUE,FALSE)];
# Rater 2:
twin1_rater2 <- twin1[ , c(FALSE,TRUE)];

# TWIN 2:
# Rater 1:
twin2_rater1 <- twin2[ , c(TRUE,FALSE)];
# Rater 2:
twin2_rater2 <- twin2[ , c(FALSE,TRUE)];

# for running the rater-based model:
rerecombined <- cbind(recombined[1],                     # FISNumber_1
                      twin1_rater1, twin1_rater2,        # twin 1, rater 1 then rater 2
                      recombined[302],                   # FISNumber_2
                      twin2_rater1, twin2_rater2,        # twin 2, rater 1 then rater 2
                      recombined[603:ncol(recombined)]); # reminder of the info...


## 5) save to csv for later use:
write.csv(rerecombined, "../data/intermediate/data_for_SEM.csv");


# Checks:
# Number of MZ pairs;
nrow(recombined[which(recombined$zygmzdz==1),]) # 249

# Total number of DZ pairs:
nrow(recombined[which(recombined$zygmzdz==2),]) # 91

# Number of DZ (same-sex) pairs:
nrow(recombined[which(recombined$twzyg==2),]) # 5 DZM
nrow(recombined[which(recombined$twzyg==4),]) # 24 DZF = 29 DZ same-sex total

# Number of DOS pairs: 
nrow(recombined[which(recombined$twzyg>=5),]) # 62

# Table for zygosity number when pairs are split into MZF, MZM, DZF, DZM, DOSmf:
(twzyg_bar_pairs <- table(recombined$TWZYG));

# Table for twin zygosity number per twin pair (MZ or DZ):
(zyg_bar_pairs <- table(recombined$ZYGMZDZ));



# Main figure for all descriptive statistics:
if( FALSE ) # don't actually plot it!
{
  tiff(file="../figures/descriptives_2.tiff", width=2*4, height=2*4, units="in", res=600, compression="lzw");
  #mar= is for margins: mar=c(bottom, left, top, right)
  par(mar=c(5,7,5,3)+0.1,mfrow=c(2,2));
  barplot(studies_bar, main="Distribution of twins by MRI study", 
          ylab="Number of twins", cex.lab=1.3, 
          cex.main=1.4, las=2,cex.axis=1.3);
  mtext("A", side=3, line=2, adj=-0.30, cex=1.2);
  hist(raw.data$ICV_new, main="Intracranial volume (ICV)", xlab=expression('ICV (mm'^3*')'), 
       ylab="Number of twins", breaks=12,cex.main=1.4,cex.axis=1.3,cex.lab=1.3);
  mtext("B", side=3, line=2, adj=-0.30, cex=1.2);
  hist(raw.data$Age_MRI_1, main="Frequency of age at time of MRI", 
       xlab="Age", ylab="Number of twins", breaks=12,cex.main=1.4,cex.axis=1.3,cex.lab=1.3);
  mtext("C", side=3, line=2, adj=-0.30, cex=1.2);
  barplot(twzyg_bar_pairs, main="Twin pair zygosity", las=2,cex.main=1.4,cex.axis=1.3,cex.lab=1.3,ylab="Number of twin pairs", ylim=c(0,200));
  mtext("D", side=3, line=2, adj=-0.30, cex=1.2);
  dev.off();
}






