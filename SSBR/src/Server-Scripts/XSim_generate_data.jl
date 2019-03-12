## Script for running on Server for multiple dataset JWAS analysis

using Pkg
### Add Pkgs if using on remote server
using XSim, JWAS, CSV, DelimitedFiles, DataFrames, Statistics, StatsBase, LinearAlgebra, XSimPreProces

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
