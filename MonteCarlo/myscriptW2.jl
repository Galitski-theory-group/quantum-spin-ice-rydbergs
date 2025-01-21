include("WormPyrochlore.jl")
include("wormsanitychecks.jl")
configin=initializebycopying([8,8,8])
intmatrixfourier=loadfourierint("4X8X8X8W2Fourier.txt")
res=(RKenergyMC(configin,intmatrixfourier,[8,8,8],400))[1]
filename=ARGS[1]
open(filename, "w") do io
    writedlm(io,res)
end