
# This document explains the content of the `final_heritability_results.csv` file.



*N.B.* this document uses GitHub-flavoured Markdown, in particular \<sup\>x\</sup\> for superscript <sup>x</sup> and \<sub\>x\</sub\> for subscript <sub>x</sub>.
The columns are:



## General:

- **short\_name**, **domain** and **type**: phenotype info



## For the genetic decomposition model:

- **genetic\_model**, **CD\_meaning** and **H2\_meaning**: given that we cannot estimate simultaneously the *C* and *D* components, we must chose between the two using the heuristic [if *r*<sub>MZ</sub> < 2*r*<sub>DZ</sub> then we estimate *C* otherwise we estimate *D*]; this, in turn, determines the fitted genetic model (*ACE* or *ADE*), and the meaning of the "CD" (short for "C or D") parameter (*C* or *D*) and of "H2" (as familiality or broad-sense heritability);

- for each of **h2**, **cd2** and **e2**:
  + the point estimate (i.e., of *h*^2^, (*c*^2^ or *d*^2^), and *e*^2^, respectively)

  + **\_CiL**, **\_CiH**: the lower and higher (upper) bounds of the 95% confidence intervals (CIs) of the corresponding point estimate

  + **\_p** and **\_p_adj**: the results of the comparison between the model with the component and without it (the *p* is significant if including the component and excluding it are significantly different), as well as after Holm's multiple testing correction

  + **H2\_signif**: for *H*^2^ and *F*^2^ respectively there is no model comparison and we check instead of 0 falls within the 95%CI or not to judge the significance of this component



## Observed twin correlations:

- **rMZ** and **rDZ**: the phenotypic correlations between monozygotic *r*<sub>MZ</sub> and dizygotic *r*<sub>DZ</sub> twins



## For inter-rater agreement ICC(C,1):

- **ICC**, **ICC\_CiL** and **ICC\_CiH**: the point estimate and 95%CI

- **ICC\_0.75** and **ICC\_0.75\_strict**: if the ICC(C,1) is greater than 75% and if the low 95%CI is greater than 75%, respectively



## For inter-rater agreement from the GSEM model agr(GSEM) we ive the same as for ICC(C,1) above



## The strength of evidence: for each relevant component we report:

- **\_star**: is the component significantly different from 0?

- **\_star_c**: is the component significantly different from 0 after Holm's multiple testing correction?

- **\_gr**: is the component greater than the agreed threshold?

- **\_grgr**: is the lower 95%CI bound greater than the agreed threshold?



## For inter-rater agreement ICC(C,1):

- **\_plus**: is the agreement greater than the agreed threshold?

- **\_plusplus**: is the lower 95%CI bound greater than the agreed threshold?



## These are combined to produce classes of strength of evidence (please note that they mean different things for diferent components as described in the analysis report, but that the lower the number the greater the strength of evidence):

- **h2\_class**, **c2\_class**, and **d2\_class**: roman numerals, the lower the better, with I providing very strong evidence and V no evidence at all



      