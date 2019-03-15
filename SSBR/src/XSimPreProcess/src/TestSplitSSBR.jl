using Pkg
function test_split(y::Vector{Float64}, g::Matrix{Int64}, ped::Matrix{Int64}, ratio1::Float64, ratio2::Float64)
    ## ratio1 is for training testing split
    ## ratio 2 is for NG:G split 
    ## NG:G must be >= 0.3 for good SSBR results assuming same number in each generation
    
    rows_t = length(y);
    rows_ng = convert(Int64, (ratio2)*rows_t);
    rows_g = convert(Int64, (1.0 - ratio2)*rows_t);
    rows_ng_tr = convert(Int64, (1.0 - ratio1)*rows_ng);
    rows_g_tr = convert(Int64, (1-0 - ratio1)*rows_ng);
    
    ## number of rows in each category 
    ## splitting works as such (1:rows_t) -> (1:rows_ng, rows_ng+1:rows_t) 
    ## (1:rows_ng_tr, row_ng_tr+1:rows_ng)::(tr_ng, test_ng)
    ## (rows_ng+1:rows_ng+rows_g_tr, rows_ng+rows_g_tr+1:rows_t)
    ## training and testing indices for pheno and pedigree 
    
    tr_ind = vcat(collect(1:(rows_ng_tr)), collect((rows_ng + 1):(rows_ng+rows_g_tr))) ;
    test_ind = vcat(collect((rows_ng_tr+1):rows_ng), collect((rows_ng+rows_g_tr + 1):rows_t));
    
    ## these are indices for phenotype and pedigree training and testing
    y_train = y[tr_ind]
    y_test = y[test_ind]
    ped_train = ped[tr_ind]
    ped_test = ped[test_ind]
            
    ## for genotype split
    geno_rows, geno_cols = size(g);
    geno_rows_tr = convert(Int64, (1-ratio1)*geno_rows);
    g_train = g[1:geno_rows_tr, :];
    g_test = g[1:(geno_rows - geno_rows_tr), :];
    
    return y_train, ped_train, g_train, y_test, ped_test, g_test
end