include("WormPyrochlore.jl")
include("wormsanitychecks.jl")
include("W3W4functions.jl")
res=RK_MC_with_L2M2([8,8,8],500000,10000)
filename="L2M2_"*ARGS[1]*".txt"
open(filename, "w") do io
    writedlm(io,res)
end