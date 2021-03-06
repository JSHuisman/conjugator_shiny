---
title: "analyse_data_intro"
author: "J.S.Huisman"
date: "1/5/2020"
output: html_document
---



### Analyse experimental data
Here you can analyse experimental data from plasmid conjugation experiments and compute the corresponding conjugation rates.

Once you upload data under the "DRT" tab, you will be able to select which conjugation rate estimates to calculate. We recommended using the Simonsen end-point formula ("SM"), if all strains involved in conjugation grow at the same rate, and the ASM end-point formula ("ASM"), if this is not the case. The other methods are heuristics, which require less information on the growth rate of individual strains but are also more prone to biases. Mouse over the "i" next to each estimate's name to get information on the required columns in the experimental data frame.

These methods to estimate conjugation rates will generally cease to be accurate once the experiment is measured past the *minimal critical time*. At that time either the contribution of transconjugants to the overall conjugation has become substantial, or the recipient dynamics are dominated by conjugation rather than growth. We provide estimates of the critical time to help determine whether the measured conjugation rate estimates are to be trusted. In general the calculation of the critical time requires an estimate of the transconjugant to recipient conjugation rate $\gamma_T$, measured in a separate TRT experiment (for more information see [Huisman *et al.*, 2020](https://www.biorxiv.org/content/10.1101/2020.03.09.980862v1)). However, even if this data is not provided, we show how the critical times are expected to vary as a function of the ratio between $\gamma_D$ and $\gamma_T$.
