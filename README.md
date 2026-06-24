# GINAMO : Galaxy tools and workflow
GINAMO (Genetic Indicators for NAture MOnitoring) focuses on developing best practices for estimating effective population size from genetic data and evaluating genetic indicators based on non-genetic data. You can find more information about the GINAMO project here: https://ginamo.org/

This repository contains the code for the tool developped by GINAMO, which are required to calculate genetic diversity indicators on Galaxy.

These tool are integrated into Galaxy workflows on Galaxy Ecology Platform (https://ecology.usegalaxy.eu/) and on Galaxy Europe (https://usegalaxy.eu/). 

This workflows have been designed for the following data type : 
-	SNPs on VCF file format
-	SSR/microsatellites on « tabular » format
-	Proxy data (*to be determined*)

Here’s how to choose which workflow(s) to use based on your data: 
[add decision tree]

You can find all the workflows in :
- Galaxy Ecology > Workflow > Workflow Public > search GINAMO
- Galaxy Europe > Workflow > Workflow Public > search GINAMO

Workflows names, description and links : 
-	**GINAMO : VCF filtering** : This workflow quality filters a VCF file. Filters will be applied in the following order : Genotype quality, Read depth, biallelic SNPs only, loci with missing data, minor allele count, individuals with missing data and heterozygosity. For each filter, you can select the parameter value. (https://ecology.usegalaxy.eu/published/workflow?id=e5cacefc738764f9)
-	**GINAMO : From SNPs to genetic EBVs** : This workflow splits VCF file into individual population file and computes genetic essentiel biodiversity variables (EBVs), including diversity, inbreeding, differentiation and effective population size. (https://ecology.usegalaxy.eu/published/workflow?id=8f5c84fb0d286050)
-	**GINAMO : SSR filtering** : This workflow quality filters microsatellite data on tabular format. Filters will be applied in the following order : individuals with missing data, loci with missing data and null alleles.  (*add links*)
-	**GINAMO: From SSRs to genetic EBVs** : This workflow computes genetic essential biodiversity variables (EBVs), including diversity, inbreeding, differentiation and effective population size using microsatellite data. (https://ecology.usegalaxy.eu/published/workflow?id=5dc5e056c5f4e0da)
-	**GINAMO : Population delineation and genetic clustering** : This workflow enables the delineation of populations using genetic data. It performs a dAPC and a estimates admixture coefficients using sparse Non-Negative Matrix Factorization algorithms. If you already have a preliminary population delineation, this allows you to test you populations using pairwise Fst. (*add links*)

(*add guidelines to use correctly Galaxy Ecology and how to use the workflows*)

(*add links to Galaxy Guidelines and future GINAMO workflows guidelines -- decision tree*)
