
# This document explains the content of the `SEM_results.csv` file.



*N.B.* this document uses GitHub-flavoured Markdown, in particular \<sup\>x\</sup\> for superscript <sup>x</sup> and \<sub\>x\</sub\> for subscript <sub>x</sub>.
The columns are:



## General:

- **phenotype\_name**: the full phenotype name

- **phenotype\_number**: the phenotype number

- **phenotype\_name\_short**: the phenotype short name



## For the genetic decomposition model:

- **numObs**, **numParams**, **observedStatistics**, **degreesOfFreedom**, **AIC**, **BIC**, **fit**, **statusCode**: number of observations, of parameters, observed statistics, degrees of freedom, AIC, BIC, fit and error code (if any)

- **genetic\_model**, **CD\_meaning**, and **H2\_meaning**: given that we cannot estimate simultaneously the *C* and *D* components, we must chose between the two using the heuristic [if *r*<sub>MZ</sub> < 2*r*<sub>DZ</sub> then we estimate *C* otherwise we estimate *D*]; this, in turn, determines the fitted genetic model (*ACE* or *ADE*), the meaning of the "CD" (short for "C or D") parameter (*C* or *D*), and the meaning of the "H2" (actually *H*<sup>2</sup>) estimate ("familiality" for *C* and "broad sense heritability" for *D*, respectively);

- Corrected or raw twin correlations used to decide between the ADE and ACE models:

  + **twin_correlations**: the raw or the corrected twin correlations were used?

  + **twin_corr_rMZ**, **twin_corr_rDZ**: the MZ and DZ correlations

- for each of **A**, **CD** and **E**:
  + the point estimate (i.e., of *A*, (*C* or *D*), and *E*, respectively)

  + **\_CiL**, **\_CiH**: the lower and higher (upper) bounds of the 95% confidence intervals (CIs) of the corresponding point estimate (i.e., of *A*, (*C* or *D*), and *E*, respectively); `NA` should probably be interpreted as 0.0

  + **\_deltaAIC** and **\_p**: the results of the comparison between the model with the component (e.g., *A*) and without it; the AIC difference is negative if including the component fits the data better than when excluding it; the *p* is significant if including the component and excluding it are significantly different

- The latent (standardized)  estimates of relevant parameters, i.e., the standardized variance terms corrected for rater effects; for all, the point estimate and the lower and upper 95% CIs:

  + **h2**: the narrow heritability *h*<sup>2</sup>

  + **cd2**: *C*<sup>2</sup> of *D*<sup>2</sup>, respectively

  + **H2**: the familiality or broad sense heritability *H*<sup>2</sup>, respectively

  + **e2**: *E*<sup>2</sup>

- The rater-specific estimates of various parameters (with 95% CIs):
  + **rh2**: the narrow heritability *h*<sup>2</sup>

  + **rcd2**: *C*<sup>2</sup> of *D*<sup>2</sup>, respectively

  + **rH2**: the familiality or broad sense heritability *H*<sup>2</sup>, respectively

  + **re2**: *E*<sup>2</sup>

  + **rerr2**: the standardized rater error variance, i.e., (1 - rater reliability)

- The rater reliability:

  + **rreliability**: rater reliability (1 - standardized rater error variance)

- **Comments**: any issues chosing the genetic model?


## For the "classic" (i.e., non-SEM) heritability estimates:

- **rMZ** and **rDZ**: the raw (i.e., uncorrected) observed phenotypic correlations between monozygotic *r*<sub>MZ</sub> and dizygotic *r*<sub>DZ</sub> twins

- **classic\_h2**: the "classic" narrow heritability estimate *h*<sup>2</sup> = 2(*r*<sub>MZ</sub> - *r*<sub>DZ</sub>)



      Script took 42.7410683194796 minutes to run...
