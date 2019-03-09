# SSBR-JWAS-Implementation
For genomic estimated breeding value simulations using JWAS and XSim, the workflow shows pre-processing for XSim generated
pedigree, phenotype and genotype files. Pre-processing of the simulated datasets is required to do BLUP and Bayesian analaysis
using JWAS.

## Steps
1. push!(LOAD_PATH, "/Path/To/My/Module/") to add XSimPreProcess Module to julia. Requires this push for module use everytime julia is reopned
2. Putting this statement in the file ~/.julia/config/startup.jl will extend LOAD_PATH on every Julia startup.
3. Generate XSim Data for multiple generations and pre-process the data for use in SSBR for JWAS

## Functionalities
1. XSimPreProcess Module Pre-processes XSim package's simulated phenotypic and genotypic data into a format compatible for JWAS package
2. Generating training and testing sets in 70:30 ratio 
3. 5 processed output files required for JWAS a) Phenotype b) Genotype c) Pedigree(if single step/optional) d)observation ID's (for EBV prediction) e) mapfile(for hypothesis testing) in the Processed-Data folder
4. 3 raw input files from XSim in the Data folder
5. A python script to generate multiple folders for multiple simulations with XSim
