using Plots
using LinearAlgebra
using SparseArrays
using LightGraphs
using GraphPlot
using DisplayAs
using DelimitedFiles
using FFTW
using Random
using Dates
using Base.Threads

#=
Conventions used:
Periodic boundary conditions unless mentioned
The diamond lattice is an fcc lattice with a 2-site basis at each fcc site. We call them A and B respectively.
Cartesian coordinates r_A =(0,0,0)+ r_l and r_B=(1/2,1/2,1/2)*(2/√3)+r_l.
r_l for a lattice site l specified by (l₁,l₂,l₃) has Cartesian coordinates r_l= (2/√3)*(l₂+l₃,l₃+l₁,l₁+l₂).
Note that our fcc primitive basis vectors are
a₁=(2/√3)(0,1,1)
a₂=(2/√3)(1,0,1)
a₃=(2/√3)(1,1,0)

Update 10/14/22: Spins are now Ising (±1 and not 0 or 1)
=#

function nextdiamsite(currsite,AvsB,bonddir,Nlist)
    #bonddir is 0,1,2 or 3
    #AvsB is of the currsite.
    nextsite=[currsite[1],currsite[2],currsite[3]]
    if bonddir!=0
        nextsite[bonddir]=mod(nextsite[bonddir]-1-AvsB,Nlist[bonddir])+1
    end
    return [nextsite,-AvsB]
end

function fetchspinsite(currsite,nextsite,bonddir,AvsB,spinconfig)
    #AvsB is for currsite not nextsite.
    # We use bonddir+1, because the 0,1,2,3 sites of the tetrahedron are labelled as 1,2,3,4 in Julia.
    if AvsB==-1
        return [[bonddir+1,nextsite[1],nextsite[2],nextsite[3]],spinconfig[bonddir+1,nextsite[1],nextsite[2],nextsite[3]]]
    elseif AvsB==1
        return [[bonddir+1,currsite[1],currsite[2],currsite[3]],spinconfig[bonddir+1,currsite[1],currsite[2],currsite[3]]]
    end
end


function findshortflippableloop(spinconfig::Array{Int64,4},Nlist::Array{Int64,1})
    for i in 1:10000000
        startover=0
        σAB=[0,0] #Initialize

        currsite=mod.(rand(Int,3),Nlist)+[1,1,1]
        AvsB=1 # 1 for A sublattice and -1 for B

        bonddir=mod(rand(Int),4)
        looplistdiam=[[currsite,AvsB]]
        bondslist=[bonddir]

        newsite=(nextdiamsite(currsite,AvsB,bonddir,Nlist))[1]
        AvsB=-AvsB
        σAB[1]=(fetchspinsite(currsite,newsite,bonddir,1,spinconfig))[2]#σA (σ for a bond pointing from A to B)
        σAB[2]=-σAB[1]#σB
        #This sets the requirement that any bond towards an A diamond site should have the pyrochlore site spin for that bond to equal sA≡σAB[1]. Similarly for B.
        pyrsiteslist=[(fetchspinsite(currsite,newsite,bonddir,1,spinconfig))[1]]
        prevbonddir=bonddir
        currsite=newsite
        push!(looplistdiam,[currsite,AvsB])

        for it in 1:100000000
            if startover==1
                break
            end

            availables=[]
            for b in 0:3
                if prevbonddir!=b
                    newsite=(nextdiamsite(currsite,AvsB,b,Nlist))[1]
                    pyrsitespin=fetchspinsite(currsite,newsite,b,AvsB,spinconfig)
                    if pyrsitespin[2]==σAB[((-AvsB+1)÷2)+1] #ensures flippability
                        push!(availables,[newsite,b,pyrsitespin[1]])
                    end
                end
            end
            if length(availables)==0
                #println("Deadend! Start over")
                startover=1
                break
            end
            randint=mod(rand(Int),length(availables))+1
            newsite=(availables[randint])[1]

            if ([newsite,-AvsB] in looplistdiam)==true
                closed=1
            else
                closed=0
            end
            currsite=newsite
            AvsB=-AvsB
            prevbonddir=(availables[randint])[2]
            
            push!(looplistdiam,[currsite,AvsB])
            push!(bondslist,prevbonddir)
            push!(pyrsiteslist,((availables[randint])[3]))

            if closed==1
                newstart=findfirst(isequal([newsite,AvsB]),looplistdiam)
                if newstart!=1
                    #println("Snipping out the first $(newstart-1) nodes.")
                    for k in 1:newstart-1
                        popfirst!(looplistdiam)
                        popfirst!(bondslist)
                        popfirst!(pyrsiteslist)
                    end
                end
                return [looplistdiam,bondslist,pyrsiteslist]
            end
        end
    end
end

function enumerateconfigs(config0,Nlist,Nmax)
    configbag=[config0]
    #println(typeof(configbag))
    L=1
    #=
    loopbag=(newfindflippableloop(config0,Nlist))[1]
    confignew=loopflip_config(config0,loop,Nlist)
    push!(configbag,copy(confignew))
    L+=1
    configprev=confignew=#
    configprev=config0
    confignew=configprev #Initialize. Don't take this seriously
    loopbag=findshortflippableloop(configprev,Nlist) #initialize
    repeats=0
    ind=0
    succ=0
    cond=0
    
    while repeats!=Nmax
        while cond==0
            if ind%2==0
                loopbag=findshortflippableloop(configprev,Nlist)
            elseif ind%2==1
                loopbag=findLongflippableloop(configprev,Nlist)
            end
            firstsite=(loopbag[3])[1]
            if winding(loopbag[2],((loopbag[1])[1])[2],configprev[firstsite[1],firstsite[2],firstsite[3],firstsite[4]])==[0,0,0]
                cond=1
            end
        end
        cond=0
        confignew=newloopflip_config(configprev,loopbag[3],Nlist)
        for i in 1:L
            if confignew==configbag[i]
                repeats+=1
                if repeats%2000==0
                    println("Oh oh $(repeats)th repeat")
                end
                succ=0
                configprev=confignew
                break
            end
            if i==L
                push!(configbag,copy(confignew))
                configprev=confignew
                succ=1
            end
        end
        if succ==1 
            L+=1
        end
        ind+=1
    end
    println("Total $L")
    return [L,configbag[L÷2]]
end

function nextdiamsite(currsite,AvsB,bonddir,Nlist)
    #bonddir is 0,1,2 or 3
    #AvsB is of the currsite.
    nextsite=[currsite[1],currsite[2],currsite[3]]
    if bonddir!=0
        nextsite[bonddir]=mod(nextsite[bonddir]-1-AvsB,Nlist[bonddir])+1
    end
    return [nextsite,-AvsB]
end

function fetchspinsite(currsite,nextsite,bonddir,AvsB,spinconfig)
    #AvsB is for currsite not nextsite.
    # We use bonddir+1, because the 0,1,2,3 sites of the tetrahedron are labelled as 1,2,3,4 in Julia.
    if AvsB==-1
        return [[bonddir+1,nextsite[1],nextsite[2],nextsite[3]],spinconfig[bonddir+1,nextsite[1],nextsite[2],nextsite[3]]]
    elseif AvsB==1
        return [[bonddir+1,currsite[1],currsite[2],currsite[3]],spinconfig[bonddir+1,currsite[1],currsite[2],currsite[3]]]
    end
end

function findshortflippableloop(spinconfig,Nlist)
    for i in 1:10000000
        startover=0
        σAB=[0,0] #Initialize

        currsite=mod.(rand(Int,3),Nlist)+[1,1,1]
        AvsB=1 # 1 for A sublattice and -1 for B

        bonddir=mod(rand(Int),4)
        looplistdiam=[[currsite,AvsB]]
        bondslist=[bonddir]

        newsite=(nextdiamsite(currsite,AvsB,bonddir,Nlist))[1]
        AvsB=-AvsB
        σAB[1]=(fetchspinsite(currsite,newsite,bonddir,1,spinconfig))[2]#σA (σ for a bond pointing from A to B)
        σAB[2]=-σAB[1]#σB
        #This sets the requirement that any bond towards an A diamond site should have the pyrochlore site spin for that bond to equal sA≡σAB[1]. Similarly for B.
        pyrsiteslist=[(fetchspinsite(currsite,newsite,bonddir,1,spinconfig))[1]]
        prevbonddir=bonddir
        currsite=newsite
        push!(looplistdiam,[currsite,AvsB])

        for it in 1:100000000
            if startover==1
                break
            end

            availables=[]
            for b in 0:3
                if prevbonddir!=b
                    newsite=(nextdiamsite(currsite,AvsB,b,Nlist))[1]
                    pyrsitespin=fetchspinsite(currsite,newsite,b,AvsB,spinconfig)
                    if pyrsitespin[2]==σAB[((-AvsB+1)÷2)+1] #ensures flippability
                        push!(availables,[newsite,b,pyrsitespin[1]])
                    end
                end
            end
            if length(availables)==0
                #println("Deadend! Start over")
                startover=1
                break
            end
            randint=mod(rand(Int),length(availables))+1
            newsite=(availables[randint])[1]

            if ([newsite,-AvsB] in looplistdiam)==true
                closed=1
            else
                closed=0
            end
            currsite=newsite
            AvsB=-AvsB
            prevbonddir=(availables[randint])[2]
            
            push!(looplistdiam,[currsite,AvsB])
            push!(bondslist,prevbonddir)
            push!(pyrsiteslist,((availables[randint])[3]))

            if closed==1
                newstart=findfirst(isequal([newsite,AvsB]),looplistdiam)
                if newstart!=1
                    #println("Snipping out the first $(newstart-1) nodes.")
                    for k in 1:newstart-1
                        popfirst!(looplistdiam)
                        popfirst!(bondslist)
                        popfirst!(pyrsiteslist)
                    end
                end
                return [looplistdiam,bondslist,pyrsiteslist]
            end
        end
    end
end

function initializebycopying(Nlist)
    #Nlist should be an array of even numbers
    #filename="4X2X2X2config.txt"
    #configsmallflat=readdlm(filename, '\t', Float64, '\n',comments=true)
    #configsmall=reshape(configsmallflat,(4,2,2,2))
    config=fill(0,4,Nlist[1],Nlist[2],Nlist[3])
    #println(typeof(config))
    configpipi0=fill(0,4,2,2,2);
    for k in 1:2
        for j in 1:2
            for i in 1:2
                configpipi0[1,i,j,k]=1*(1-2*mod(i+j-2,2))
                configpipi0[2,i,j,k]=1*(1-2*mod(i+j-2,2))
                configpipi0[3,i,j,k]=-1*(1-2*mod(i+j-2,2))
                configpipi0[4,i,j,k]=-1*(1-2*mod(i+j-2,2))
            end
        end
    end
    configsmall=(enumerateconfigs(configpipi0,[2,2,2],3000))[2]
    for k in 1:Nlist[3]
        for j in 1:Nlist[2]
            for i in 1:Nlist[1]
                for l in 1:4
                    config[l,i,j,k]=configsmall[l,(mod(i-1,2)+1),(mod(j-1,2)+1),(mod(k-1,2)+1)]
                end
            end
        end
    end
    return config
end

function findLongflippableloop(spinconfig,Nlist)
    for i in 1:10000000
        startover=0
        σAB=[0,0] #Initialize

        currsite=mod.(rand(Int,3),Nlist)+[1,1,1]
        AvsB=1 # 1 for A sublattice and -1 for B

        bonddir=mod(rand(Int),4)
        looplistdiam=[[currsite,AvsB]]
        bondslist=[bonddir]

        newsite=(nextdiamsite(currsite,AvsB,bonddir,Nlist))[1]
        AvsB=-AvsB
        σAB[1]=(fetchspinsite(currsite,newsite,bonddir,1,spinconfig))[2]#σA (σ for a bond pointing from A to B)
        σAB[2]=-σAB[1]#σB
        #This sets the requirement that any bond towards an A diamond site should have the pyrochlore site spin for that bond to equal sA≡σAB[1]. Similarly for B.
        pyrsiteslist=[(fetchspinsite(currsite,newsite,bonddir,1,spinconfig))[1]]
        prevbonddir=bonddir
        currsite=newsite
        push!(looplistdiam,[currsite,AvsB])

        for it in 1:100000000
            if startover==1
                break
            end

            availables=[]
            for b in 0:3
                if prevbonddir!=b
                    newsitefull=nextdiamsite(currsite,AvsB,b,Nlist)
                    newsite=newsitefull[1]
                    pyrsitespin=fetchspinsite(currsite,newsite,b,AvsB,spinconfig)
                    if pyrsitespin[2]==σAB[((-AvsB+1)÷2)+1] #ensures flippability
                        if ((newsitefull in looplistdiam)==false)|| (newsitefull==looplistdiam[1])
                            push!(availables,[newsite,b,pyrsitespin[1]])
                        end
                    end
                end
            end
            if length(availables)==0
                #println("Deadend! Start over")
                startover=1
                break
            end
            randint=mod(rand(Int),length(availables))+1
            newsite=(availables[randint])[1]

            if ([newsite,-AvsB]==looplistdiam[1])  
                closed=1
            else
                closed=0
            end
            currsite=newsite
            AvsB=-AvsB
            prevbonddir=(availables[randint])[2]
            
            push!(looplistdiam,[currsite,AvsB])
            push!(bondslist,prevbonddir)
            push!(pyrsiteslist,((availables[randint])[3]))

            if closed==1
                return [looplistdiam,bondslist,pyrsiteslist]
            end
        end
    end
end

function newloopflip_config(spinconfig,looppyr,Nlist)
    newconfig=copy(spinconfig)
    L=length(looppyr)
    indarray=[0,0,0,0]
    for i in 1:L
        indarray=[(looppyr[i])[1],(looppyr[i])[2],(looppyr[i])[3],(looppyr[i])[4]]
        newconfig[indarray[1],indarray[2],indarray[3],indarray[4]]=-spinconfig[indarray[1],indarray[2],indarray[3],indarray[4]]
    end
    return newconfig
end

function winding(bondslist,AvsBinit,σinit)
    AvsB=AvsBinit #AvsB of starting point of bond and not end point
    windnum=[0,0,0]
    for i in 1:length(bondslist)
        if bondslist[i]!=0
            windnum[bondslist[i]]+=-AvsB
        end
        AvsB=-AvsB
    end
    return windnum*AvsBinit*σinit
end

function CheckFlippable(loopbag,config,Nlist)
    σ=0 #Initialize
    L=length(loopbag[1])-1
    if L%2==1
        return 0
    end
    looplistpyr=[]
    for i in 1:L
        push!(looplistpyr,(fetchspinsite((loopbag[1])[i][1],(loopbag[1])[i+1][1],(loopbag[2])[i],(loopbag[1])[i][2],config))[2])
        if i!=1 
            if looplistpyr[i]==looplistpyr[i-1]
                return 0
            end
        end
    end
    return 1
end

function SignedCheckFlippable(loopbag,config,Nlist)
    σ=0 #Initialize
    L=length(loopbag[1])-1
    if L%2==1
        return 0
    end
    looplistpyr=[]
    for i in 1:L
        push!(looplistpyr,(fetchspinsite((loopbag[1])[i][1],(loopbag[1])[i+1][1],(loopbag[2])[i],(loopbag[1])[i][2],config))[2])
        if i!=1 
            if looplistpyr[i]==looplistpyr[i-1]
                return 0
            end
        end
    end
    return looplistpyr[1]
end

function fetchhexplaquette(fccsite,orientation,Nlist)
    loop1=fill([[0,0,0],1],7)
    loop2=fill([[0,0,0],1],7)
    loop3=fill([[0,0,0],1],7)
    loop4=fill([[0,0,0],1],7)
    bondlist1=zeros(Int,6)
    bondlist2=zeros(Int,6)
    bondlist3=zeros(Int,6)
    bondlist4=zeros(Int,6)

    #7 because first and last elements are the same.
    if orientation==1
        loop1[1]=[fccsite,1]
        bondlist1[1]=0
        loop1[2]=[fccsite,-1]
        bondlist1[2]=1
        loop1[3]=[(mod.(fccsite+[-1,-1,-1]+[1,0,0],Nlist)+[1,1,1]),1]
        bondlist1[3]=2
        loop1[4]=[(mod.(fccsite+[-1,-1,-1]+[1,-1,0],Nlist)+[1,1,1]),-1]
        bondlist1[4]=0
        loop1[5]=[(mod.(fccsite+[-1,-1,-1]+[1,-1,0],Nlist)+[1,1,1]),1]
        bondlist1[5]=1
        loop1[6]=[(mod.(fccsite+[-1,-1,-1]+[0,-1,0],Nlist)+[1,1,1]),-1]
        bondlist1[6]=2
        loop1[7]=loop1[1]
        return [loop1,bondlist1]    
    elseif orientation==2
        loop2[1]=[fccsite,1]
        bondlist2[1]=0
        loop2[2]=[fccsite,-1]
        bondlist2[2]=3
        loop2[3]=[(mod.(fccsite+[-1,-1,-1]+[0,0,1],Nlist)+[1,1,1]),1]
        bondlist2[3]=2
        loop2[4]=[(mod.(fccsite+[-1,-1,-1]+[0,-1,1],Nlist)+[1,1,1]),-1]
        bondlist2[4]=0
        loop2[5]=[(mod.(fccsite+[-1,-1,-1]+[0,-1,1],Nlist)+[1,1,1]),1]
        bondlist2[5]=3
        loop2[6]=[(mod.(fccsite+[-1,-1,-1]+[0,-1,0],Nlist)+[1,1,1]),-1]
        bondlist2[6]=2
        loop2[7]=loop2[1]
        return [loop2,bondlist2]
    elseif orientation==3
        loop3[1]=[fccsite,1]
        bondlist3[1]=0
        loop3[2]=[fccsite,-1]
        bondlist3[2]=3
        loop3[3]=[(mod.(fccsite+[-1,-1,-1]+[0,0,1],Nlist)+[1,1,1]),1]
        bondlist3[3]=1
        loop3[4]=[(mod.(fccsite+[-1,-1,-1]+[-1,0,1],Nlist)+[1,1,1]),-1]
        bondlist3[4]=0
        loop3[5]=[(mod.(fccsite+[-1,-1,-1]+[-1,0,1],Nlist)+[1,1,1]),1]
        bondlist3[5]=3
        loop3[6]=[(mod.(fccsite+[-1,-1,-1]+[-1,0,0],Nlist)+[1,1,1]),-1]
        bondlist3[6]=1
        loop3[7]=loop3[1]
        return [loop3,bondlist3]
    elseif orientation==4
        loop4[1]=[fccsite,1]
        bondlist4[1]=1
        loop4[2]=[(mod.(fccsite+[-1,-1,-1]+[-1,0,0],Nlist)+[1,1,1]),-1]
        bondlist4[2]=3
        loop4[3]=[(mod.(fccsite+[-1,-1,-1]+[-1,0,1],Nlist)+[1,1,1]),1]
        bondlist4[3]=2
        loop4[4]=[(mod.(fccsite+[-1,-1,-1]+[-1,-1,1],Nlist)+[1,1,1]),-1]
        bondlist4[4]=1
        loop4[5]=[(mod.(fccsite+[-1,-1,-1]+[0,-1,1],Nlist)+[1,1,1]),1]
        bondlist4[5]=3
        loop4[6]=[(mod.(fccsite+[-1,-1,-1]+[0,-1,0],Nlist)+[1,1,1]),-1]
        bondlist4[6]=2
        loop4[7]=loop4[1]
        return [loop4,bondlist4]
    end
end

function CountFlippableHex(config,Nlist)
    count=0
    for k in 1:Nlist[3]
        for j in 1:Nlist[2]
            for i in 1:Nlist[1]
                #looplistset=fetchhexplaquette([i,j,k],Nlist)
                #println("$i,$j,$k:")
                for a in 1:4      
                    #println("$a: $((looplistset[a])[1])")
                    if CheckFlippable(fetchhexplaquette([i,j,k],a,Nlist),config,Nlist)==1
                        count+=1
                    end
                end
            end
        end
    end
    return count
end

function runMCfreely(config,Nlist,Nsweep)
    configcurr=copy(config)
    count=0
    #count=CountFlippableHex(configcurr,Nlist)
    firstsite=[1,1,1,1] #initialize
    loopbag=findshortflippableloop(configcurr,Nlist) #initialize
    Ndiamsites=Nlist[1]*Nlist[2]*Nlist[3]
    cond=0
    for i in 1:((Nsweep)*Ndiamsites)
        while cond==0
            loopbag=findshortflippableloop(configcurr,Nlist)
            firstsite=(loopbag[3])[1]
            if winding(loopbag[2],((loopbag[1])[1])[2],configcurr[firstsite[1],firstsite[2],firstsite[3],firstsite[4]])==[0,0,0]
                cond=1
            end
        end
        cond=0

        configcurr=newloopflip_config(configcurr,loopbag[3],Nlist)
    end
    return configcurr
end


function RKnflipMC(config,Nlist,Nsweep,Nskip)
    configcurr=copy(config)
    count=0
    #count=CountFlippableHex(configcurr,Nlist)
    firstsite=[1,1,1,1] #initialize
    loopbag=findshortflippableloop(configcurr,Nlist) #initialize
    Ndiamsites=Nlist[1]*Nlist[2]*Nlist[3]
    cond=0
    totalskips=Nskip*Ndiamsites
    for i in 1:((Nsweep+Nskip)*Ndiamsites)
        while cond==0
            loopbag=findshortflippableloop(configcurr,Nlist)
            firstsite=(loopbag[3])[1]
            if winding(loopbag[2],((loopbag[1])[1])[2],configcurr[firstsite[1],firstsite[2],firstsite[3],firstsite[4]])==[0,0,0]
                cond=1
            end
        end
        cond=0

        configcurr=newloopflip_config(configcurr,loopbag[3],Nlist)
        if i>totalskips
            if i%Ndiamsites==0
                count+=CountFlippableHex(configcurr,Nlist)
            end
        end
    end
    count=count/(Nsweep)
    return count/(Ndiamsites*4)
end

function YYcorrSnap(config,Nlist,relativeFCCseparation,orientation1,orientation2)
    Nunitcells=Nlist[1]*Nlist[2]*Nlist[3]
    count=0
    fcc1=[1,1,1]
    fcc2=[1,1,1]
    for k in 1:Nlist[3]
        for j in 1:Nlist[2]
            for i in 1:Nlist[1]
                fcc1=[i,j,k]
                fcc2=mod.(fcc1+[-1,-1,-1]+relativeFCCseparation,Nlist)+[1,1,1]
                count+=-SignedCheckFlippable(fetchhexplaquette(fcc1,orientation1,Nlist),config,Nlist)*SignedCheckFlippable(fetchhexplaquette(fcc2,orientation2,Nlist),config,Nlist)
            end
        end
    end
    return count/Nunitcells
end

function tabulateYYcorrSnap(config,Nlist)
    counttable=zeros((Nlist[1]÷2)+1,Nlist[2],Nlist[3],4,4)
    #The relative separation is (i-1,j-1,k-1).
    #i.e. The first loop is fcc site r, orientation l. The second loop is fcc site r+(i-1,j-1,k-1), orientation m.
    for m in 1:4
        for l in 1:4
            for k in 1:Nlist[3]
                for j in 1:Nlist[2]
                    for i in 1:(Nlist[1]÷2)+1
                        counttable[i,j,k,l,m]=YYcorrSnap(config,Nlist,[i-1,j-1,k-1],l,m)
                    end
                end
            end
        end
    end
    return counttable
end

function RKyyCorrMC(config::Array{Int64,4},Nlist::Array{Int64,1},relativeFCCseparation::Array{Int64,1},orientation1::Int64,orientation2::Int64,Nsweep::Int64,Nskip::Int64)
    configcurr=copy(config)
    corr=0
    #count=CountFlippableHex(configcurr,Nlist)
    firstsite=[1,1,1,1] #initialize
    @time loopbag=findshortflippableloop(configcurr,Nlist) #initialize
    Ndiamsites=Nlist[1]*Nlist[2]*Nlist[3]
    cond=0
    totalskips=Nskip*Ndiamsites
    for i in 1:((Nsweep+Nskip)*Ndiamsites)
        while cond==0
            @time loopbag=findshortflippableloop(configcurr,Nlist)
            firstsite=(loopbag[3])[1]
            if winding(loopbag[2],((loopbag[1])[1])[2],configcurr[firstsite[1],firstsite[2],firstsite[3],firstsite[4]])==[0,0,0]
                cond=1
            end
        end
        cond=0

        configcurr=newloopflip_config(configcurr,loopbag[3],Nlist)
        if i>totalskips
            if i%Ndiamsites==0
                corr+=YYcorrSnap(configcurr,Nlist,relativeFCCseparation,orientation1,orientation2)
            end
        end
    end
    corr=corr/(Nsweep)
    return corr
end

function RKyyCorrTableMC(config,Nlist,Nsweep,Nskip)
    configcurr=copy(config)
    counttable=zeros((Nlist[1]÷2)+1,Nlist[2],Nlist[3],4,4)
    #count=CountFlippableHex(configcurr,Nlist)
    firstsite=[1,1,1,1] #initialize
    loopbag=findshortflippableloop(configcurr,Nlist) #initialize
    Ndiamsites=Nlist[1]*Nlist[2]*Nlist[3]
    cond=0
    totalskips=Nskip*Ndiamsites
    for i in 1:((Nsweep+Nskip)*Ndiamsites)
        while cond==0
            loopbag=findshortflippableloop(configcurr,Nlist)
            firstsite=(loopbag[3])[1]
            if winding(loopbag[2],((loopbag[1])[1])[2],configcurr[firstsite[1],firstsite[2],firstsite[3],firstsite[4]])==[0,0,0]
                cond=1
            end
        end
        cond=0

        configcurr=newloopflip_config(configcurr,loopbag[3],Nlist)
        if i>totalskips
            if i%Ndiamsites==0
                counttable+=tabulateYYcorrSnap(configcurr,Nlist)
            end
        end
    end
    counttable=counttable/(Nsweep)
    return counttable
end