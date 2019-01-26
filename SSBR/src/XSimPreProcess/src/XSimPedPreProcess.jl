function ped_ssbr_convert_1(pedfile)
  ped_size = size(pedfile, 1)
  ped_Array = Array{String}(undef, ped_size + 1);
  ped_Array[1] = "ID,Sire,Dam";
  temp_arr = Array{String}(undef, 3)

  ## convert the pedigree array from Float64 to Int and append the header to the file
  ## For Pedigree file IDs, Sires and Dams must be string in format

  f = (x) -> Int(x);
  @time for i in 1:ped_size
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
    ped_size = size(ped_Array, 1)
    new_Arr = Array{String}(undef, ped_size);

    for j in 2:ped_size
        x = split(string(ped_Array[j]), ",")
        x = "a" .* x
        #println(x)
        new_Arr[j] = x[1]*","*x[2]*","*x[3]
    end
    new_Arr[1] = "ID,sire,dam"
    return new_Arr
end
