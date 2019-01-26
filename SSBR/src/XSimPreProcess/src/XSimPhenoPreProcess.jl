
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

###################################################################################
## Do similar conversion for genotypes except keep the values from 2:end - {Float}
## create a copy of the phenofile
## PHENOTYPE FILE FORMATTING
function pheno_ssbr_format(phenofile, rowIDs)
    #rowIDs = row_IDs(phenofile);
    pheno_copy = phenofile;

#     ## extract rowIDs and convert to Int
#     rowIDs = Array{Any}(undef, 14001);
#     rowIDs[2:end] = pedfile[:, 1];
#     rowIDs[1] = "ID";
#     f = (x) -> Int(x);
#     for i in 2:14001
#         rowIDs[i] = f.(rowIDs[i])
#     end

# ## for Ids in phenotype, genotype and pedigree add a string "a"
#     for i in 2:14001
#         rowIDs[i] = "a"*string(rowIDs[i])
#     end

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
