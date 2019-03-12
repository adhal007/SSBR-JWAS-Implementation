# SSBR-JWAS-Implementation
For genomic estimated breeding value(EBV) and Genome wide association studies(GWAS) using JWAS and XSim. The workflow includes 
1) Simulating SNP data using XSim. Files generated: pedigree.txt, phenotype.txt, genotype.txt, map.txt and ID.txt 
2) Statistical Analysis for EBV and GWAS using JWAS
3) Includes XSimPreProcess package for automating XSim -> JWAS analysis.

## Steps for adding XSimPreProcess Module for XSim package.
1. push!(LOAD_PATH, "/Path/To/My/Module/") to add XSimPreProcess Module to julia. Requires this push for module use everytime julia is reopned
2. Putting this statement in the file ~/.julia/config/startup.jl will extend LOAD_PATH on every Julia startup.
3. Generate XSim Data for multiple generations and pre-process the data for use in SSBR for JWAS

