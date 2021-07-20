---
title: "simulation_plot_explanation"
author: "J.S. Huisman"
date: "12/6/2019"
output: html_document
---



The left plot shows the population dynamics of $D, R$ and $T$. The middle plot shows the Simonsen estimate of the growth rate of the total mating population. Dashed lines indicate the actual values of $\psi_D, \psi_R$ and $\psi_T$ (selected with the input sliders). The right plot compares two estimates of the conjugation rate: $\gamma_{max}$ from the Simonsen end-point formula and $\gamma_{Dmax}$ based on the ASM end-point method (independent of the model used for simulation). Dashed lines indicate the actual values of $\gamma_D$ and $\gamma_T$ (selected with the input sliders). 

The critical times $t_{crit1}, t_{crit2}, t_{crit3}$ indicate the time points at which the assumptions underlying the ASM cease to be valid. To get an accurate estimate of the conjugation rate, one should ideally measure the mating populations before the first critical time. These critical times and the time at which stationary phase is reached ($t_{stat}$) can be toggled on and off with the 'Display critical times' button. The end of the simulation is given by the 'measurement time'.
