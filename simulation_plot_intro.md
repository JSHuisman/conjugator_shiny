---
title: "simulation_plot_intro"
author: "J.S. Huisman"
date: "12/6/2019"
output: html_document
---



### Simulate bacterial population dynamics
On this page one can simulate the dynamics of three bacterial populations: donors ($D$), recipients ($R$) and transconjugants ($T$). The populations compete for the same resources ($C$), as they grow with growth rates ($\psi_{D/R/T}$). Donors conjugate with recipients at rate $\gamma_D$, and transconjugants conjugate with recipients with rate $\gamma_T$. 
Each time the parameters are changed (using the sliders), or a new model is selected, the dynamics are simulated anew.

Three population dynamic models can be selected: 
- **SM: the Simonsen model**, which assumes the growth and conjugation rate are the same for all populations (i.e. $\psi_{max} = \psi_{D} = \psi_{R} = \psi_{T}$, and $\gamma_{max} = \gamma_{D} = \gamma_{T}$), 
- **ESM: the extended Simonsen model**, which allows for population specific growth and conjugation rates, 
- **ASM: the approximate extended Simonsen model**, which simplifies the extended Simonsen model to a case with infinite resources.
