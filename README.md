# SSBR-JWAS-Implementation
For genomic estimated breeding value(EBV) and Genome wide association studies(GWAS) using JWAS and XSim. The workflow includes 
1) Simulating SNP data using XSim. Files generated: pedigree.txt, phenotype.txt, genotype.txt, map.txt and ID.txt <br/>
2) Statistical Analysis for EBV and GWAS using JWAS <br/>
3) Includes XSimPreProcess package for automating XSim -> JWAS analysis. <br/>
4) To use GWAS function for WPPA in JWAS package run ```Pkg.add(PackageSpec(name=“JWAS”, rev=“master”))``` and then import JWAS.misc by ```using JWAS.misc ```
5) ```Test-Split.jl``` splits your input files into training and testing sets for Single Step Analysis on JWAS

## Steps for adding XSimPreProcess Module for XSim package.
1. ```push!(LOAD_PATH, "/Path/To/My/Module/")``` to add XSimPreProcess Module to julia. Requires this push for module use everytime julia is reopned
2. Copy/paste the command in the file ```~/.julia/config/startup.jl``` will extend LOAD_PATH on every Julia startup.

## Jupyter-Examples.
Usage of XSimPreProcess module for SSBR-analysis in this folder

## Updated Commandline Implementation 
ARGS = 5 

ARG[1] = Input Genotype Path
ARG[2] = Processed Data Path
ARG[3] = Analysis Data Path 
ARG[4] = Number of QTL to set for analysis 
ARG[5] = Number of Markers to include for Phenotype Simulation 

1. Make 3 Directories - 1 for input data, 1 for processed data, 1 for analysis 
2. Files inside Input Data -  Genotype file, XSim-Map-File
3. Files inside Processed Data - Processed genotype, phenotype, pedigree, mapfile, rowIDs and other files as per pre-processing script 
4. Files inside Analysis Data - Output files of JWAS run 
5. Run command as follows ```julia Phenotype-SSBR-Simulation.jl .\Directory1// .\Directory2// .\Directory3// 5 1 ```

