
#################################################################################
## GENOTYPE FILE FORMATTING
## use the genotype file and do same transformation for genotypes
## genotype script
function geno_ssbr_format(genofile, rowIDs, n_markers)
    geno_copy = genofile;
    #pedfile = readdlm("../../XSim-data\\15800_ped_train_rowID.txt");


    ## hcat rowIDs
    geno_cat1 = Array{Any}(undef, 5881, n_markers+1);
    geno_cat1[2:5881, :] = hcat(rowIDs, geno_copy)
    println("genotype file concatenated with ID's: ", "a1, 1, 0, 2")

    ## add markerIDs and vcat
    empty_l = Array{String}(undef, n_markers)
    empty_l[1] = "ID"
    empty_l[2:n_markers] = repeat(["m"], n_markers-1)

    for i in 2:n_markers
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
    ## return the H_split for 1st column of geno_formatted


    ## Main split for genotype file
    ## Main Split
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
        #println(z)
#         geno_string[j] = string(rowIDs[j])
#         #println(geno_string[j])
#         geno_string[j] = geno_string[j]*","*z
#         #println(geno_string[j])
        geno_string[j] = rowIDs[j-1]*","*z
    end
    println("Matrix split into", "a1, 0, 1, 2")
    geno_string[1] = H_split
    return  geno_string
end

###########################################################
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
