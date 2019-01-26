## make this into a module using
module XSimPreProcess
using Pkg
    #using DelimitedFiles, CSV, LinearAlgebra
import Base.include
import Base.show
import Base

include("XSimPedPreProcess.jl")
include("XSimPhenoPreProcess.jl")
include("XSimGenoPreProcess.jl")
export ped_ssbr_convert_1, ped_ssbr_append_char
export pheno_row_IDs_test, pheno_row_IDs_train, pheno_ssbr_format
export geno_ssbr_format, geno_row_IDs

end
