---
title: "app_home"
author: "J.S. Huisman"
date: "12/6/2019"
output: html_document
---



Quantifying plasmid conjugation just got easier! This app is designed to make the estimation and reporting of plasmid conjugation rates from liquid mating cultures easier, more accurate, and more comparable. All data analysis functionality can also be found in the [accompanying R package](https://github.com/JSHuisman/conjugator), and the methods are described in more detail in the [preprint](https://www.biorxiv.org/content/10.1101/2020.03.09.980862v1). 
<br/>


#### Analyse experimental data
The tab **Analyse experimental data** allows you to upload your own data from liquid culture mating experiments and compute the corresponding plasmid conjugation rates. The two recommended methods to estimate these rates are **the Simonsen end-point formula** (Simonsen *et al*, J.Gen. Microbiol, 1990), applicable if all strains involved in conjugation grow at the same rate, and **the ASM end-point formula**, which relaxes this assumption (Huisman *et al*, 2020). 

Both of these formulae will cease to be accurate once either the contribution of transconjugants to the overall conjugation becomes substantial, or the recipient dynamics are dominated by conjugation rather than growth. Based on the relative timing of these events, we derived the **critical time** within which the conjugation rate estimates remain valid. 
<br/>


#### Simulate population dynamics 
The tab **Simulate population dynamics** allows users to simulate bacterial population dynamics under different models.
We currently support the following models:
- **SM: the Simonsen model** (Simonsen *et al*, J.Gen. Microbiol, 1990), which assumes equal growth rates of donor, recipient and transconjugant populations, and equal transfer rates from donor to recipient and from transconjugant to recipient.
- **ESM: the extended Simonsen model**, which allows for population-specific growth and conjugation rates.
- **ASM: the approximate extended Simonsen model**, which simplifies the extended Simonsen model to a case with infinite resources.
<br/>


#### Contact and citation 
This app was created by Jana S. Huisman, in collaboration with Fabienne Benz, Sarah J.N. Duxbury, J. Arjan G.M. de Visser, Alex R. Hall, Egil A.J. Fischer, and Sebastian Bonhoeffer.

For questions or suggestions concerning the app and the implemented methods please contact:<br/> jana.huisman [at] env.ethz.ch .<br/> 

If the conjugator app or R package helped you in your work, please cite the manuscript:<br/>
[Jana S. Huisman, Fabienne Benz, Sarah J.N. Duxbury, J. Arjan G.M. de Visser, Alex R. Hall, Egil A.J. Fischer, Sebastian Bonhoeffer. Estimating the rate of plasmid transfer in liquid mating cultures. *BioRxiv* 2020.](https://www.biorxiv.org/content/10.1101/2020.03.09.980862v1)
