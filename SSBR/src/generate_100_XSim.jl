using Pkg
### Add Pkgs if using on remote server
using XSim, JWAS, CSV, DelimitedFiles, DataFrames, Statistics, StatsBase, LinearAlgebra
### Introduce for loop for generating multiple XSim datasets.
### Functions ##########################################################################
function ped_ssbr_convert_1(pedfile)
    ## pedfile input in only one format
    ped_Array = Array{String}(undef, 14001);
    ped_Array[1] = "ID,Sire,Dam";
    temp_arr = Array{String}(undef, 3)
    
    ## convert the pedigree array from Float64 to Int and append the header to the file
    ## For Pedigree file IDs, Sires and Dams must be string in format
    f = (x) -> Int(x);
    @time for i in 1:14000
        y = pedfile[i, :]
        y = f.(y)
        y = string(y)
        z = ""
        s = split(y, "")
        l = length(s)
        for j in 1:length(s)
            if j !== 1 && j !== length(s) && y[j] !== ' '
             z = z*y[j]
            end
        end
    # this converts "[1, 2, 3]" to "1, 2, 3"
        ped_Array[i + 1] = z
    end
    return ped_Array
end

function ped_ssbr_append_char(ped_Array)
    ## split the file into components
    new_Arr = Array{String}(undef, 14001);

    for j in 2:14001
        x = split(string(ped_Array[j]), ",")
        x = "a" .* x
        #println(x)
        new_Arr[j] = x[1]*","*x[2]*","*x[3]
    end
    new_Arr[1] = "ID,sire,dam"
    return new_Arr
end

function pheno_row_IDs_train(phenofile, pedfile)
    no_rows = size(phenofile, 1)
    rowIDs = Array{Any}(undef, no_rows);
    rowIDs[2:end] = pedfile[1:(no_rows)-1, 1];
    rowIDs[1] = "ID";
    f = (x) -> Int(x);
    for i in 2:no_rows
        rowIDs[i] = f.(rowIDs[i])
    end

    for j in 2:no_rows
        rowIDs[j] = "a"*string(rowIDs[j])
    end
    #iter_size = no_rows
    return rowIDs[2:end]
end

function pheno_row_IDs_test(phenofile, pedfile)
    no_rows = size(phenofile, 1)
    rowIDs = Array{Any}(undef, no_rows+1);
    rowIDs[2:end] = pedfile[1:(no_rows), 1];
    rowIDs[1] = "ID";
    f = (x) -> Int(x);
    for i in 2:no_rows+1
        rowIDs[i] = f.(rowIDs[i])
    end

    for j in 2:no_rows+1
        rowIDs[j] = "a"*string(rowIDs[j])
    end
    #iter_size = no_rows
    return rowIDs[2:end]
end

## rowID extraction for phenotype and genotypes is different
function geno_row_IDs(genofile, pedfile)
    no_rows = size(genofile, 1) + 1
    rowIDs = Array{Any}(undef, no_rows);
    rowIDs[2:end] = pedfile[1:(no_rows-1), 1];
    rowIDs[1] = "ID";
    f = (x) -> Int(x);
    for i in 2:no_rows
        rowIDs[i] = f.(rowIDs[i])
    end

    for j in 2:no_rows
        rowIDs[j] = "a"*string(rowIDs[j])
    end
    #iter_size = no_rows
    return rowIDs[2:end]
end

###################################################################################
## Do similar conversion for genotypes except keep the values from 2:end - {Float}
## create a copy of the phenofile
## PHENOTYPE FILE FORMATTING
function pheno_ssbr_format(phenofile, rowIDs)
    #rowIDs = row_IDs(phenofile);
    pheno_copy = phenofile;
    ## concatenate the rowIDs
    pheno_cat = hcat(rowIDs, pheno_copy[2:end])
    pheno_array = Array{String}(undef, 4201)
    pheno_array[1] = "ID,y"
    for j in 2:4201
        m1 = pheno_cat[j-1, 2];
        m2 = string(m1);

        pheno_array[j] = pheno_cat[j-1, 1]*","*m2
    end
    return pheno_array
end

## GENOTYPE FILE FORMATTING
## use the genotype file and do same transformation for genotypes
## genotype script


function geno_ssbr_format(genofile, rowIDs)
    geno_copy = genofile;
    
    ## hcat rowIDs
    geno_cat1 = Array{Any}(undef, 5881, 3001);
    geno_cat1[2:5881, :] = hcat(rowIDs, geno_copy)
    println("genotype file concatenated with ID's: ", "a1, 1, 0, 2")
    
    ## add markerIDs and vcat
    empty_l = Array{String}(undef, 3000)
    empty_l[1] = "ID"
    empty_l[2:3000] = repeat(["m"], 2999)

    for i in 2:3000
        empty_l[i] = empty_l[i] * "$(i-1)"
    end
    
    len = 2999
    println("$len markerIDs created")
    ## Header split for genotype file
    H_split = ""
    for j in 1:3000
        H_split = H_split*empty_l[j]        
        if j == 3000
            H_split = H_split
        else
            H_split = H_split*","
        end
    end
    println("markerIDs split into:", "ID, m1, m2, m3", H_split )
    ## return the H_split for 1st column of geno_formatted    
    ## Main split for genotype file 
    ## convert the float to Int and then string split 
    f = x -> Int(x);
    geno_string = Array{String}(undef, 5881);
    @time for j in 2:5881
        v = geno_cat1[j, 2:end]
        #println(v)
        v = f.(v)
        v = string(v)
        v = split(v, "")
        len = length(v)
        z = ""
        #println(v)
        for i in 1:len
            if i !== 1 && i !== 2 && i !== 3 && i !== len
                z = z*v[i]
            end
        end
        geno_string[j] = rowIDs[j-1]*","*z
    end
    println("Matrix split into", "a1, 0, 1, 2")
    geno_string[1] = H_split
    return  geno_string
end

function generate_XSim_data()
    for i in 1:99
        #numChr,numLoci,chrLength,mutRate = 2,[1,2],[0.1,0.2],0.0
        chrLength= 1.0  #length of each chromosome
        numChr   = 2    #number of chromosomes
        nmarkers = 1500   #number of loci for each chromosome
        nQTL     = 25    #number of QTL for each chromosomefects,mutRate);
        mutRate = 0.000025
        #geneFreq = fill(0.5,nmarkers)
        build_genome(numChr,chrLength,nmarkers,nQTL)

        # generate 4 generations

        #generation 0
        #generate founders
        popSizeFounder = 100
        sires = sampleFounders(popSizeFounder);
        dams  = sampleFounders(popSizeFounder);

        #random mating
        ngen,popSize = 1,2800
        sires1,dams1,gen1 = sampleRan(popSize, ngen, sires, dams);


        ngen, popSize=1, 2800
        sires2,dams2,gen2 = sampleRan(popSize, ngen, sires1, dams1);

        ngen, popSize=1, 2800
        sires3,dams3,gen3 = sampleRan(popSize, ngen, sires2, dams2);

        ngen, popSize= 1, 2800
        sires4,dams4,gen4 = sampleRan(popSize, ngen, sires3, dams3);


        ngen, popSize = 1, 2800
        sires5, dams5, gen5 = sampleRan(popSize, ngen, sires4, dams4);

        # #collect the animals data
        animals5 = concatCohorts(sires5,dams5);
        animals4 = concatCohorts(sires4, dams4);
        animals3 = concatCohorts(sires3, dams3);
        animals2 = concatCohorts(sires2, dams2);
        animals1 = concatCohorts(sires1, dams1);
        animals0 = concatCohorts(sires, dams);

        #get genotypes
        M5 = getOurGenotypes(animals5);
        M4 = getOurGenotypes(animals4);
        M3 = getOurGenotypes(animals3);

        M = vcat(M3, M4, M5);

        #get phenotypes
        Phen1 = getOurPhenVals(animals1, 0.1);
        var1 = var(Phen1)
        Phen2 = getOurPhenVals(animals2, 0.5*var1);
        Phen3 = getOurPhenVals(animals3, 0.5*var1);
        Phen4 = getOurPhenVals(animals4, 0.5*var1);
        Phen5 = getOurPhenVals(animals5, 0.5*var1);
        Phen = vcat(Phen1, Phen2, Phen3, Phen4, Phen5);


        #get pedigree
        Ped5 = getPedigree(animals5);
        Ped4 = getPedigree(animals4);
        Ped3 = getPedigree(animals3);
        Ped2 = getPedigree(animals2);
        Ped1 = getPedigree(animals1)
        Ped = vcat(Ped1, Ped2, Ped3, Ped4, Ped5);

        #write out
        Phenotypes = DataFrame();
        Phenotypes[:y] = Phen;

        f_pheno = open("/home/adhal/home/Projects/SSBR3/XSim-Data$i/Phenotypes_XSim_Raw$i.txt", "w")
        f_geno = open("/home/adhal/home/Projects/SSBR3/XSim-Data$i/Genotypes_XSim_Raw$i.txt", "w")
        f_ped = open("/home/adhal/home/Projects/SSBR3/XSim-Data$i/Pedigree_XSim_Raw$i.txt", "w")
        outputCatData("/home/adhal/home/Projects/SSBR3/XSim-Data$i/Processed-Data/map_XSim_Raw.txt")        
        CSV.write(f_pheno, Phenotypes)
        writedlm(f_geno, M)
        writedlm(f_ped, Ped)
        close(f_pheno)
        close(f_geno)
        close(f_ped)

        #Phenotypes = readdlm("/home/adhal/home/Projects/SSBR3/XSim-Data$i/Phenotypes_XSim_Raw$i.txt");
        #M = readdlm("/home/adhal/home/Projects/SSBR3/XSim-Data$i/Genotypes_XSim_Raw$i.txt");
        #Ped = readdlm("/home/adhal/home/Projects/SSBR3/XSim-Data$i/Pedigree_XSim_Raw$i.txt");

        ## Making trainingandtest
        ## Make sure to conserve orderintheconcatenation

        ## pheno_setup
        ## pheno training setup
       
	pheno_train_ng = Phenotypes[1:3920, :];
        pheno_train_g = Phenotypes[5601:11480, :];
        pheno_train = vcat(pheno_train_ng, pheno_train_g);
	CSV.write("/home/adhal/home/Projects/SSBR3/XSim-Data$i/pheno_train.txt", pheno_train)

        ## pheno testing setup
        pheno_test_ng = Phenotypes[3921:5600, :];
        pheno_test_g = Phenotypes[11481:14000, :];
        pheno_test = vcat(pheno_test_ng, pheno_test_g)

        ## geno_setup
        geno_train_g = M[1:5880, :]
        geno_test_g = M[5881:8400, :];

        ## Split pedigree file for extracting rowIDs for the simualtion
        ## ped_training setup
        ped_train_ng = Ped[1:3920, :];
        ped_train_g = Ped[5601:11480, :];
        ped_train = vcat(ped_train_ng, ped_train_g);

        ## ped_testing setup
        ped_test_ng = Ped[3921:5600, :];
        ped_test_g = Ped[11481:14000, :];
        ped_test = vcat(ped_test_ng, ped_test_g);

        ## preprocess and save pedigree file
        ped_out_1 = ped_ssbr_convert_1(Ped);
        ped_out = ped_ssbr_append_char(ped_out_1);
        ped_out_process = open("/home/adhal/home/Projects/SSBR3/XSim-Data$i/Processed-Data/pedigree_processed.txt", "w");
        writedlm(ped_out_process, ped_out);
        close(ped_out_process)

        ## get ObsID
        obs_ID = pheno_row_IDs_test(pheno_test, ped_test);
        obs_out = open("/home/adhal/home/Projects/SSBR3/XSim-Data$i/Processed-Data/Obs-ID.txt", "w")

        ## get genotype pre-processed
        geno_train_IDs = geno_row_IDs(geno_train_g, ped_train_g);
	geno_out = geno_ssbr_format(geno_train_g, geno_train_IDs);
        geno_name = open("/home/adhal/home/Projects/SSBR3/XSim-Data$i/Processed-Data/geno_processed.txt", "w");
        writedlm(geno_name, geno_out)
        close(geno_name)

        ## get pheno-preprocessed
	pheno_train = readdlm("/home/adhal/home/Projects/SSBR3/XSim-Data$i/pheno_train.txt")
	## pheno_train IDs
	pheno_train_IDs = pheno_row_IDs_train(pheno_train, ped_train);

	## phenotype test IDs
	pheno_test_IDs = pheno_row_IDs_test(pheno_test, ped_test);

        pheno_train_processed = pheno_ssbr_format(pheno_train, pheno_train_IDs);
        ## write out pheno_train to output
        pheno_train_process_out = open("/home/adhal/home/Projects/SSBR3/XSim-Data$i/Processed-Data/pheno_processed.txt", "w");
        writedlm(pheno_train_process_out, pheno_train_processed);
        close(pheno_train_process_out)
 
	end
   
end

generate_XSim_data()


