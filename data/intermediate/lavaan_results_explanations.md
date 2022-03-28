
# This document explains the content of the `lavaan_results.csv` file.



*N.B.* this document uses GitHub-flavoured Markdown, in particular \<sup\>x\</sup\> for superscript <sup>x</sup> and \<sub\>x\</sub\> for subscript <sub>x</sub>.
The columns are:



## General:

- **phenotype\_name**: the full phenotype name

- **phenotype\_number**: the phenotype number

- **phenotype\_name\_short**: the phenotype short name



## For the genetic decomposition model:

- **genetic\_model**, **CD\_meaning**, and **H2\_meaning**: given that we cannot estimate simultaneously the *C* and *D* components, we must chose between the two using model compaison with AIC (given by **deltaAIC_ACE_ADE** = AIC(ACE) - AIC(ADE)); this, in turn, determines the fitted genetic model (*ACE* or *ADE*), the meaning of the "CD" (short for "C or D") parameter (*C* or *D*), and the meaning of the "H2" (actually *H*<sup>2</sup>) estimate ("familiality" for *C* and "broad sense heritability" for *D*, respectively);

- The fit measures of the full model (ACE or ADE):

  + **fit_Chisq**, **fit_Chisq.df** and **fit_Chisq.p**: the Chi squared test fit; **fit_AIC** and **fit_BIC** are AIC and BIC of the model, and **fit_RMSEA** the model's RMSEA

- The latent (standardized) estimates of relevant parameters, i.e., the standardized variance terms corrected for rater effects; for all, the point estimate, the lower and upper 95% CIs, and the *p*-value; for *a* and *cd*, we also give the model compairson with the model without the component as **deltaAIC** = AIC(model with component) - AIC(model without component) and as LR test's Chi squared test **LRT.Chisq**, **LRT.Chisq.df** and **LRT.Chisq.p**:

  + **a2**: the narrow heritability *h*<sup>2</sup>

  + **cd2**: *C*<sup>2</sup> of *D*<sup>2</sup>, respectively

  + **e2**: *E*<sup>2</sup>



      Script took 10.5611515959104 minutes to run...
