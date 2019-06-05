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
