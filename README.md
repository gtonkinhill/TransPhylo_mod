# TransPhylo

# Introduction

This is a modified version of TransPhylo, a software package that can reconstruct infectious disease transmission using genomic data. The input is a dated phylogeny, where leaves correspond to pathogens sampled from the known infected hosts. The main output is a transmission tree which indicates who infected whom, including the potential existence of unsampled individuals who may have acted as missing transmission links. TransPhylo works by colouring the branches of the phylogeny using a separate colour for each host, sampled or not. Each section of the tree  coloured in a unique colour represents the pathogen evolution happening within a distinct host. Changes of colours on branches therefore correspond to transmission events from one host to another.

For example, in the tree below, the outbreak started with the unsampled host 8, who transmitted to sampled host 4, who transmitted to unsampled host 3, who transmitted to both sampled hosts 1 and 2. Host 8 also transmitted to unsampled host 7 who transmitted to both sampled hosts 5 and 6.

<img src="https://raw.githubusercontent.com/wiki/xavierdidelot/TransPhylo/example.png" width="300">

For a more formal description of TransPhylo, see the following preprint:

Didelot, Fraser, Gardy and Colijn (2016)
Genomic infectious disease epidemiology in partially sampled and ongoing outbreaks
http://biorxiv.org/content/early/2016/07/22/065334

# Installation

You can install TransPhylo mod in R using the following command:

`devtools::install_github("gtonkinhill/TransPhylo_mod", auth_token="37064ad364627e23377a492cf790285642256a23")`

# Documentation

The TransPhylo tutorial can be downloaded [here](https://raw.githubusercontent.com/wiki/xavierdidelot/TransPhylo/TransPhylo-Tutorial.pdf). This tutorial is also available within the R package as a vignette. The TransPhylo reference manual can be downloaded [here](https://raw.githubusercontent.com/wiki/xavierdidelot/TransPhylo/TransPhylo-RefMan.pdf).

