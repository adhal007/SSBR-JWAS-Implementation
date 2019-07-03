using Pkg
Pkg.add("DelimitedFiles")
Pkg.add("JWAS")
Pkg.add("CSV")
Pkg.add("XSim")
Pkg.add("Distributions")
Pkg.add("LinearAlgebra")
Pkg.add("DataFrames")
Pkg.add("Statistics")
Pkg.add("StatsBase")
#### Adding Packages to Julia Module

using DelimitedFiles, JWAS, JWAS.misc, CSV, Distributions, LinearAlgebra, DataFrames, Statistics, StatsBase, XSim

#### Importing the libraries
## path = ARGS[1]
## include(path*"XSimPreProcess.jl")

#### Add Pre-Process Functions for Phenotype, Genotype, Map, Pedigree and IDs
##################################################################################################
using Pkg
function test_split(y::Vector{Float64}, g::Matrix{Int64}, ped::Matrix{Int64})
    ## Create Phenotype and Genotype for BayesC and Single Step 
    rows_t = length(y)
    geno_ind_ssbr = vcat(rand(collect(5600:8400), 600), rand(collect(8401:11200),600), rand(collect(11201:14000), 300))
    pheno_ind = collect(1:rows_t);
    println(geno_ind_ssbr)
    println("space")
    println(pheno_ind)

    y_SSBR = y[pheno_ind]
    y_BR = y[geno_ind_ssbr]
    g_SSBR = g[geno_ind_ssbr, :];
    g_BR = g[geno_ind_ssbr, :];
    return y_SSBR, g_SSBR, y_BR,  g_BR, geno_ind_ssbr
end

## Pre-Process Pedigree file Ip: "[1, 2, 3]" to "1,2,3"
function ped_ssbr_process(pedfile)
  ped_size = size(pedfile, 1)
  ped_Array = Array{String}(undef, ped_size + 1);
  ped_Array[1] = "ID,Sire,Dam";
  f = Int64.(pedfile)
  g = string.(f)
  g = "a".*(g)
  for j in 2:ped_size+1
    ped_Array[j] = g[j-1, 1]*","*g[j-1, 2]*","*g[j-1, 3]
  end
  return ped_Array
end

##################################################################################
## rowID extraction for phenotype and genotypes is different
## Using only 1 function for obtaining phenotype IDs
function pheno_row_IDs(ped_Array)
	final_ids = Array{String}(undef, size(ped_Array,1)-1) ## since 1st row is "ID,Sire,Dam")
	ID_list = split.(ped_Array, ",")
	for i in 1:size(ped_Array,1)-1
	     final_ids[i] = ID_list[i][1]
	end
	return final_ids
end

function geno_row_IDs(ped_Array, geno_ind_ssbr)
    	final_ids = Array{String}(undef, size(geno_ind_ssbr,1)); ## since 1st row is "ID,Sire,Dam")
	ped_Array = ped_Array[2:end];
	ID_list = split.(ped_Array[geno_ind_ssbr, 1])
	for i in 1:size(ped_Array,1)-1
	     final_ids[i] = ID_list[i][1]
	end
	return final_ids
end

## give a string path to the raw mapfile input generated by XSim
function process_map(rawmapfile, nmarkers)
    mapfile = readdlm(rawmapfile)
    m = mapfile;
    new_map = Array{Any}(undef, nmarkers+1);

    for i in 1:(nmarkers+1)
        if i == 1
            new_map[i] = ""
            for j in 1:4
                if j >= 2
                    new_map[i] *= ","
                    new_map[i] *= m[i, j]
            
                else
                    new_map[i] *= m[i, j]
                end
            end
        else
            new_map[i] = ""
            for j in 1:4
                if j >= 2
                    new_map[i] *= ","
                    m[i, j] = string(m[i, j])
                    new_map[i] *= m[i, j]
                else
                    m[i, j] = string(m[i, j])
                    new_map[i] *= m[i, j]
                end
            end
        end
    end
    
    return new_map
end



###########################################################################################
## GENOTYPE FILE FORMATTING
function geno_ssbr_format(genofile, rowIDs)
    geno_copy = genofile
    n_obs, n_markers = size(genofile)
    n_obs += 1
    n_markers += 1
    ## hcat rowIDs
    geno_cat1 = Array{Any}(undef, n_obs, n_markers);
    geno_cat1[2:end, :] = hcat(rowIDs, geno_copy)
    println("genotype file concatenated with ID's: ", "a1, 1, 0, 2")
    ## add markerIDs and vcat
    empty_l = Array{String}(undef, n_markers-1)
    empty_l[1] = "ID"
    empty_l[2:end] = repeat(["m"], n_markers-2)
    
    for i in 2:n_markers-1
        empty_l[i] = empty_l[i] * "$(i-1)"
    end
    len = n_markers-1
    println("$len markerIDs created")
    ## Header split for genotype file
    H_split = ""
    for j in 1:n_markers
        H_split = H_split*empty_l[j]
        if j == n_markers
            H_split = H_split
        else
            H_split = H_split*","
        end
    end
    println("markerIDs split into:", "ID, m1, m2, m3", H_split )
    ## return the H_split for 1st column of geno_formatte
    ## Main split for genotype file 
    f = x -> Int(x);
    geno_string = Array{String}(undef, n_obs);
    @time for j in 2:n_obs
        v = geno_cat1[j, 2:end]
        println()
        println(v)
        v = f.(v)
        v = string(v)
        v = split(v, "")
        println(v)
        len = length(v)
        println()
        println(len)
        z = ""
        for i in 2:3:len-1
            if i !== len-1
                z = z*","*v[i]
            else
                z = z
            end
        end
        geno_string[j] = rowIDs[j-1]*z
    end
    println("Matrix split into", "a1, 0, 1, 2")
    geno_string[1] = H_split
    return  geno_string
end

######################################################################################
## Do similar conversion for genotypes except keep the values from 2:end - {Float}
## PHENOTYPE FILE FORMATTING
function pheno_ssbr_format(phenofile, rowIDs)
    
    if isa(phenofile, String) == true         
        phenofile = readdlm(phenofile);
    else
        phenofile = phenofile;
        end 
    nrows = length(phenofile)
    ## concatenate the rowIDs
    pheno_cat = hcat(rowIDs, phenofile[2:end])
    pheno_array = Array{String}(undef, nrows)
    pheno_array[1] = "ID,y"
    for j in 2:nrows
        m1 = pheno_cat[j-1, 2];
        m2 = string(m1);
        pheno_array[j] = pheno_cat[j-1, 1]*","*m2
    end
    return pheno_array
end

########################################################################################
## End of Functions for Processing Genotype, IDs, Phenotype, Pedigree and chromosome map
########################################################################################

### Phenotype Simulation ###############################################################
### ARGS = 1: genotype path, 2: Processed path, 3: Analysis Path 4: number of QTL 5: number of genotype markers

geno_path =  ARGS[1]
Mfull = readdlm(geno_path*"Geno_Small.txt", '\t');
nQTL = ARGS[4]
h2 = 0.5
varg = 1
geno_reduce = ARGS[5]
nInd =size(Mfull,1)
nMarkers = Int64(round(size(Mfull,2)/parse(Float64,geno_reduce)))
QTL_effects = randn(parse(Int64,nQTL))
nMarkers_lst = collect(1:nMarkers)
QTLpos = sample(nMarkers_lst,parse(Int64,nQTL))

#### save QTL pos for ROC curve
fQTL = open(geno_path*"QTLpositions.txt", "w");
writedlm(fQTL, QTLpos)
close(fQTL)

## Simulate Breeding values
BV = Mfull[:,QTLpos]*QTL_effects
BV = BV/std(BV)*sqrt(varg)
vare = (1-h2)/h2*varg
pheno = BV+randn(nInd)*vare

#### save Phenotypes
Phenotypes = DataFrame();
Phenotypes[:y] = pheno
fpheno = open(geno_path*"Phenotypes.txt", "w");
CSV.write(fpheno, Phenotypes)
close(fpheno)


## Set Processed Data Path
procpath = ARGS[2]

## pheno_setup
f_pheno = open(geno_path*"Phenotypes_XSim_Raw.txt", "w");
Ped = readdlm(geno_path*"XSim_Ped.txt", '\t');         
CSV.write(f_pheno, Phenotypes)
close(f_pheno)

## Get variance for phenotypes
v = var(pheno)
f_var = open(geno_path*"PhenotypicVariance.txt", "w")
writedlm(f_var, v)
close(f_var)


## Get Phenotypes for SSBR and BR for same set of genotypes
pheno_ssbr,  geno_ssbr, pheno_br, geno_br, geno_ind_ssbr = test_split(Float64.(pheno), Int.(Mfull[:,1:nMarkers]), Int.(Ped));

## Save raw phenotypes for SSBR
Phenotypes = DataFrame();
Phenotypes[:y] = pheno_ssbr;
CSV.write(procpath*"pheno_ssbr.txt", Phenotypes)

## preprocess and save pedigree file
ped_out_1 = ped_ssbr_process(Ped);
ped_out_process = open(procpath*"pedigree_processed.txt", "w");
writedlm(ped_out_process, ped_out_1);
close(ped_out_process)

## pre-process ObsID for phenotypes
obs_ID = pheno_row_IDs(ped_out_1);
println(size(obs_ID, 1))
obs_out = open(procpath*"Obs_ID.txt", "w")
writedlm(obs_out, obs_ID)
close(obs_out)

## pre-process phenotype
pheno_ssbr = readdlm(procpath*"pheno_ssbr.txt")
pheno_processed = pheno_ssbr_format(pheno_ssbr, obs_ID);
println(size(pheno_processed, 1))
pheno_process_out = open(procpath*"pheno_processed.txt", "w");
writedlm(pheno_process_out, pheno_processed);
close(pheno_process_out)
 
## pre-process genotype IDs and pre-process genotype
geno_ssbr_IDs = geno_row_IDs(ped_out_1, geno_ind_ssbr);
geno_out = geno_ssbr_format(geno_ssbr, geno_ssbr_IDs);
println(size(geno_out, 1))
geno_name = open(procpath*"geno_processed.txt", "w");
writedlm(geno_name, geno_out)
close(geno_name)
####################################################################################
## End of Pre-processing for Raw Files from XSim
####################################################################################

## input file reading consideration for JWAS-SSBR :
## read io streams of 1) phenotype file 2) Obs_ID files for ID's to be tested
## io streams of 1) genotype file 2) map file 3) pedigree file to keep as iostream
## process map file using XSimPreProcess since it is generated from XSim (array format different from JWAS)

phenofile  = procpath*"pheno_processed.txt" #Datasets.dataset
phenotypes1 = CSV.read(phenofile,delim = ',',header=true)
obs_ID = readdlm(procpath*"Obs_ID.txt")
mapfile = geno_path*("map_XSim_Raw.txt");
genofile   = procpath*"geno_processed.txt"#Datasets.dataset
pedfile = procpath*"pedigree_processed.txt"


## processing and reading mapfile as formatted iostream
new_map = process_map(mapfile, 45611)
processed_map = open(procpath*"new_map_file_processed.txt", "w");
writedlm(processed_map, new_map)
close(processed_map)
newmapfile = procpath*"new_map_file_processed.txt"


## Run JWAS ########################################################################
analysis_path = ARGS[3]
pedigree = get_pedigree(pedfile, separator=",",header=true)
model_equation1  ="y = intercept";
R      = 1.0
model1 = build_model(model_equation1,R);
G3 =1.0
g = add_genotypes(model1,genofile,G3, separator=',', header=true);
outputEBV(model1, obs_ID);
outputMCMCsamples(model1)
out2=runMCMC(model1, phenotypes1, Pi=0.95,estimatePi=true,methods="BayesC",single_step_analysis=true,
		                  pedigree=pedigree, chain_length=5000,output_samples_frequency=100,outputEBV=true, output_samples_file=analysis_path);
keys(out2)

mrkr_file = analysis_path*"_marker_effects_y.txt";
g = GWAS(mrkr_file, newmapfile, model1, header=true, window_size="1.0 MB", threshold=1/3100)
writedlm(analysis_path*"WPPA-results-1.txt", g);

