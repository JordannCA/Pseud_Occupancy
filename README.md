# Code and Data for "Hidden declines in a seemingly secure Australian frog species"

**Authors:** Jordann Crawford-Ash, Stephanie Pulsford, Will Osborne, Danswell Starrs, Geoffrey W. Heard, Ben C. Scheele 

**Corresponding author:** jordann.crawford-ash@anu.edu.au 

## Description

This repository contains the code and data used in the paper: "Hidden declines in a seemingly secure Australian frog species"
Submitted to Ecology and Evolution for review, October 2025. 

## Methodology

We compiled >1,000 pre-1990 occurrence records for the Pseudophryne bibronii complex across south-eastern Australia and selected 70 historically occupied sites spanning low–high elevations in the ACT, NSW and VIC. At each site, we conducted standardized 15-minute audio-visual surveys within 50 m of breeding habitat up to six times over the 2023–2024 breeding seasons, using call playback to elicit calls and recording detection and abundance in ordinal categories, along with weather and habitat variables and counts of co-occurring Crinia signifera. We measured site attributes (elevation, canopy cover, moisture, disturbance, fire severity) from field observations and spatial datasets, processed all data in R, and analysed occupancy and abundance with single-season occupancy models (using the unmarked package in R) and cumulative-link mixed models (using the ordinal package in R). A subset of adult males was swabbed for Batrachochytrium dendrobatidis and infection load quantified using qPCR. Further detail on the methodology used is provided in the methods section of the manuscript. 

## Structure
- `data/` – analysis-ready datasets
- `code/` – R scripts for analysis and figures  

## Requirements

- All analysis done in R version 4.4.1.
- Language: R 4.3 or above
- all packages required to run listed in code files. 
