using DelimitedFiles
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

###################################################################################
## Do similar conversion for genotypes except keep the values from 2:end - {Float}
## create a copy of the phenofile
## PHENOTYPE FILE FORMATTING
function pheno_ssbr_format(phenofile, rowIDs)
    
    if isa(phenofile, String) == true         
        phenofile = readdlm(phenofile);
    else
        phenofile = phenofile;
        end 
    nrows = length(phenofile) - 1
    ## concatenate the rowIDs
    pheno_cat = hcat(rowIDs, phenofile[2:end])
    pheno_array = Array{String}(undef, nrows+1)
    pheno_array[1] = "ID,y"
    for j in 2:(nrows+1)
        m1 = pheno_cat[j-1, 2];
        m2 = string(m1);

        pheno_array[j] = pheno_cat[j-1, 1]*","*m2
    end
    return pheno_array
end