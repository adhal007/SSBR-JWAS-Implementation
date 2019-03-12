### Load training data for model training 
### Generate/Build model on phenotrain with 5600 genotypes 
### Predict EBV for phenotest with 2400 genotypes
using Pkg
Pkg.add("JWAS")
#Pkg.add("JWAS.Datasets")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("LinearAlgebra")
Pkg.add("DelimitedFiles")
#Pkg.add("JWAS.misc")
Pkg.add(PackageSpec(name="JWAS", rev="master"))

using JWAS,DataFrames,CSV, LinearAlgebra, DelimitedFiles, JWAS.misc

for i in 1:99 
	newmapfile = "/home/adhal/home/Projects/SSBR4/XSim-Data$i/Processed-Data/map_XSim_Raw.txt";
	phenofile  = "/home/adhal/home/Projects/SSBR4/XSim-Data$i/Processed-Data/pheno_processed.txt" #Datasets.dataset
	genofile   = "/home/adhal/home/Projects/SSBR4/XSim-Data$i/Processed-Data/geno_processed.txt"#Datasets.dataset
	pedfile = "/home/adhal/home/Projects/SSBR4/XSim-Data$i/Processed-Data/pedigree_processed.txt"
	obs_ID = readdlm("/home/adhal/home/Projects/SSBR4/XSim-Data$i/Processed-Data/Obs-ID.txt")

	phenotypes1 = CSV.read(phenofile,delim = ',',header=true)
	pedigree=get_pedigree(pedfile,separator=",",header=true)

	#var(phenotypes1[:, 2])

	model_equation1  ="y = intercept";
	R      = 1.0
	model1 = build_model(model_equation1,R);

	G3 =1.0
	g = add_genotypes(model1,genofile,G3, separator=',', header=true);

	outputEBV(model1, obs_ID);

	outputMCMCsamples(model1)

	out2=runMCMC(model1, phenotypes1, Pi=0.95,estimatePi=true,methods="BayesC",single_step_analysis=true,
		         pedigree=pedigree, chain_length=100,output_samples_frequency=100,outputEBV=true, output_samples_file="/home/adhal/home/Projects/SSBR4/XSim-Data$i/MCMC_Samples");

	keys(out2)
	
	#EBV_out = open("/home/adhal/home/Projects/SSBR4/XSim-Data$i/Ebv_results.txt", "w")
	#writedlm(EBV_out, out2["EBV_y"])
	#close(EBV_out)

	## WPPA 
	mrkr_file = "/home/adhal/home/Projects/SSBR4/XSim-Data$i/MCMC_Samples_marker_effects_y.txt";
	g = GWAS(mrkr_file, newmapfile, model1, header=true, window_size="1 MB", threshold=0.05)
	writedlm("/home/adhal/home/Projects/SSBR4/XSim-Data$i/WPPA-results.txt", g)
end
