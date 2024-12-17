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
function nextlatticesite(currsite,AvsB,bonddir,Nlist)
    #bonddir is 0,1,2 or 3
    #AvsB is of the currsite.
    nextsite=[currsite[1],currsite[2],currsite[3]]
    if bonddir!=0
        nextsite[bonddir]=mod(nextsite[bonddir]-1-AvsB,Nlist[bonddir])+1
    end
    return nextsite
end

function fetchspin(currsite,nextsite,bonddir,AvsB,spinconfig)
    #AvsB is for the new diamond site of the (potential) bond, not the old.
    # We use bonddir+1, because the 0,1,2,3 sites of the tetrahedron are labelled as 1,2,3,4 in Julia.
    if AvsB==1
        return spinconfig[bonddir+1,nextsite[1],nextsite[2],nextsite[3]]
    elseif AvsB==-1
        return spinconfig[bonddir+1,currsite[1],currsite[2],currsite[3]]
    end
end

function fetchhexplaquette(fccsite,Nlist)
    loop1=fill([[0,0,0],1],7)
    loop2=fill([[0,0,0],1],7)
    loop3=fill([[0,0,0],1],7)
    loop4=fill([[0,0,0],1],7)
    bondlist1=zeros(Int,6)
    bondlist2=zeros(Int,6)
    bondlist3=zeros(Int,6)
    bondlist4=zeros(Int,6)

    #7 because first and last elements are the same.
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

    return [[loop1,bondlist1],[loop2,bondlist2],[loop3,bondlist3],[loop4,bondlist4]]
end

function newfetchhexplaquette(fccsite,orientation,Nlist)
    if orientation==1
        loop1=fill([[0,0,0],1],7)
        #7 because first and last elements are the same.
        loop1[1]=[fccsite,1]
        loop1[2]=[fccsite,-1]
        loop1[3]=[(mod.(fccsite+[-1,-1,-1]+[1,0,0],Nlist)+[1,1,1]),1]
        loop1[4]=[(mod.(fccsite+[-1,-1,-1]+[1,-1,0],Nlist)+[1,1,1]),-1]
        loop1[5]=[(mod.(fccsite+[-1,-1,-1]+[1,-1,0],Nlist)+[1,1,1]),1]
        loop1[6]=[(mod.(fccsite+[-1,-1,-1]+[0,-1,0],Nlist)+[1,1,1]),-1]
        loop1[7]=loop1[1]
        
        return loop1
    elseif orientation==2
        loop2=fill([[0,0,0],1],7)
        loop2[1]=[fccsite,1]
        loop2[2]=[fccsite,-1]
        loop2[3]=[(mod.(fccsite+[-1,-1,-1]+[0,0,1],Nlist)+[1,1,1]),1]
        loop2[4]=[(mod.(fccsite+[-1,-1,-1]+[0,-1,1],Nlist)+[1,1,1]),-1]
        loop2[5]=[(mod.(fccsite+[-1,-1,-1]+[0,-1,1],Nlist)+[1,1,1]),1]
        loop2[6]=[(mod.(fccsite+[-1,-1,-1]+[0,-1,0],Nlist)+[1,1,1]),-1]
        loop2[7]=loop2[1]

        return loop2
    elseif orientation==3
        loop3=fill([[0,0,0],1],7)
        loop3[1]=[fccsite,1]
        loop3[2]=[fccsite,-1]
        loop3[3]=[(mod.(fccsite+[-1,-1,-1]+[0,0,1],Nlist)+[1,1,1]),1]
        loop3[4]=[(mod.(fccsite+[-1,-1,-1]+[-1,0,1],Nlist)+[1,1,1]),-1]
        loop3[5]=[(mod.(fccsite+[-1,-1,-1]+[-1,0,1],Nlist)+[1,1,1]),1]
        loop3[6]=[(mod.(fccsite+[-1,-1,-1]+[-1,0,0],Nlist)+[1,1,1]),-1]
        loop3[7]=loop3[1]

        return loop3
    elseif orientation==4
        loop4=fill([[0,0,0],1],7)
        loop4[1]=[fccsite,1]
        loop4[2]=[(mod.(fccsite+[-1,-1,-1]+[-1,0,0],Nlist)+[1,1,1]),-1]
        loop4[3]=[(mod.(fccsite+[-1,-1,-1]+[-1,0,1],Nlist)+[1,1,1]),1]
        loop4[4]=[(mod.(fccsite+[-1,-1,-1]+[-1,-1,1],Nlist)+[1,1,1]),-1]
        loop4[5]=[(mod.(fccsite+[-1,-1,-1]+[0,-1,1],Nlist)+[1,1,1]),1]
        loop4[6]=[(mod.(fccsite+[-1,-1,-1]+[0,-1,0],Nlist)+[1,1,1]),-1]
        loop4[7]=loop4[1]

        return loop4
    end
end

function checkflippable(looplist,config,Nlist)
    σ=0 #Initialize
    L=length(looplist)-1
    result=0
    for i in 1:L+1
        lp=looplist[i][1]
        p=looplist[i][2]
        lq=looplist[mod(i,L)+1][1]
        q=looplist[mod(i,L)+1][2]

        if p==1
            if lp==lq
                if i!=1
                    if config[1,lp[1],lp[2],lp[3]]!=-σ
                        return 0
                    end
                end
                σ=config[1,lp[1],lp[2],lp[3]]
            elseif lp==[lq[1]%Nlist[1]+1,lq[2],lq[3]]
                if i!=1
                    if config[2,lp[1],lp[2],lp[3]]!=-σ
                        return 0
                    end
                end
                σ=config[2,lp[1],lp[2],lp[3]]
            elseif lp==[lq[1],lq[2]%Nlist[2]+1,lq[3]]
                if i!=1
                    if config[3,lp[1],lp[2],lp[3]]!=-σ
                        return 0
                    end
                end
                σ=config[3,lp[1],lp[2],lp[3]]
            elseif lp==[lq[1],lq[2],lq[3]%Nlist[3]+1]
                if i!=1
                    if config[4,lp[1],lp[2],lp[3]]!=-σ
                        return 0
                    end
                end
                σ=config[4,lp[1],lp[2],lp[3]]
            end
        elseif p==-1
            if lq==lp
                if i!=1
                    if config[1,lq[1],lq[2],lq[3]]!=-σ
                        return 0
                    end
                end
                σ=config[1,lq[1],lq[2],lq[3]]    
            elseif lq==[lp[1]%Nlist[1]+1,lp[2],lp[3]]
                if i!=1
                    if config[2,lq[1],lq[2],lq[3]]!=-σ
                        return 0
                    end
                end
                σ=config[2,lq[1],lq[2],lq[3]] 
            elseif lq==[lp[1],lp[2]%Nlist[2]+1,lp[3]]
                if i!=1
                    if config[3,lq[1],lq[2],lq[3]]!=-σ
                        return 0
                    end
                end
                σ=config[3,lq[1],lq[2],lq[3]]
            elseif lq==[lp[1],lp[2],lp[3]%Nlist[3]+1]
                if i!=1
                    if config[4,lq[1],lq[2],lq[3]]!=-σ
                        return 0
                    end
                end
                σ=config[4,lq[1],lq[2],lq[3]]
            end
        end
    end
    return 1
end

function signedcheckflippable(looplist,config,Nlist)
    σ=0 #Initialize
    L=length(looplist)-1
    result=0
    for i in 1:L+1
        lp=looplist[i][1]
        p=looplist[i][2]
        lq=looplist[mod(i,L)+1][1]
        q=looplist[mod(i,L)+1][2]

        if p==1
            if lp==lq
                if i!=1
                    if config[1,lp[1],lp[2],lp[3]]!=-σ
                        return 0
                    end
                end
                σ=config[1,lp[1],lp[2],lp[3]]
            elseif lp==[lq[1]%Nlist[1]+1,lq[2],lq[3]]
                if i!=1
                    if config[2,lp[1],lp[2],lp[3]]!=-σ
                        return 0
                    end
                end
                σ=config[2,lp[1],lp[2],lp[3]]
            elseif lp==[lq[1],lq[2]%Nlist[2]+1,lq[3]]
                if i!=1
                    if config[3,lp[1],lp[2],lp[3]]!=-σ
                        return 0
                    end
                end
                σ=config[3,lp[1],lp[2],lp[3]]
            elseif lp==[lq[1],lq[2],lq[3]%Nlist[3]+1]
                if i!=1
                    if config[4,lp[1],lp[2],lp[3]]!=-σ
                        return 0
                    end
                end
                σ=config[4,lp[1],lp[2],lp[3]]
            end
        elseif p==-1
            if lq==lp
                if i!=1
                    if config[1,lq[1],lq[2],lq[3]]!=-σ
                        return 0
                    end
                end
                σ=config[1,lq[1],lq[2],lq[3]]    
            elseif lq==[lp[1]%Nlist[1]+1,lp[2],lp[3]]
                if i!=1
                    if config[2,lq[1],lq[2],lq[3]]!=-σ
                        return 0
                    end
                end
                σ=config[2,lq[1],lq[2],lq[3]] 
            elseif lq==[lp[1],lp[2]%Nlist[2]+1,lp[3]]
                if i!=1
                    if config[3,lq[1],lq[2],lq[3]]!=-σ
                        return 0
                    end
                end
                σ=config[3,lq[1],lq[2],lq[3]]
            elseif lq==[lp[1],lp[2],lp[3]%Nlist[3]+1]
                if i!=1
                    if config[4,lq[1],lq[2],lq[3]]!=-σ
                        return 0
                    end
                end
                σ=config[4,lq[1],lq[2],lq[3]]
            end
        end
    end
    return config[looplist[1][1],looplist[1][2][1],looplist[1][2][2],looplist[1][2][3]]
end


function findflippableloop(spinconfig,Nlist)
    # Nlist is of the form [N1,N2,N3]
    for i in 1:10000000
        # The outermost loop is to allow deadends.
        startover=0
        σAB=[0,0] #Initialize

        init=mod.(rand(Int,3),Nlist)+[1,1,1]
        AvsB=2*mod(rand(Int),2)-1 # 1 for A sublattice and -1 for B

        looplist=[[init,AvsB]]
        bonddir=mod(rand(Int),4)
        second=nextlatticesite(init,AvsB,bonddir,Nlist)
        AvsB=-AvsB
        prevdir=bonddir
        push!(looplist,[second,AvsB])
        currsite=second
        if AvsB==1
            σAB[1]=fetchspin(init,second,bonddir,1,spinconfig)#σA
            σAB[2]=-σAB[1]#σB
        elseif AvsB==-1
            σAB[2]=fetchspin(init,second,bonddir,-1,spinconfig)
            σAB[1]=-σAB[2]
        end #This block sets the requirement that any bond towards an A diamond site should have the pyrochlore site spin for that bond to equal sA≡σAB[1]. Similarly for B.

        for it in 1:100000000
            #This is the loop through different bonds.
            if startover==1
                break
            end
            
            jobdone=0
            extra=0
            randint=mod(rand(Int),3)+1
            while jobdone==0
                if extra==3
                    #println("Deadend! Start over")
                    startover=1
                    break
                end
                bonddir=mod(prevdir+mod(randint-1+extra,3)+1,4)
                nextsite=nextlatticesite(currsite,AvsB,bonddir,Nlist)
                if fetchspin(currsite,nextsite,bonddir,-AvsB,spinconfig)==σAB[((AvsB+1)÷2)+1] # Ensures flippability.
                    if ([nextsite,-AvsB] in looplist)==false
                        AvsB=-AvsB
                        prevdir=bonddir
                        push!(looplist,[nextsite,AvsB])
                        currsite=nextsite
                        jobdone=1
                    else
                        newstart=findfirst(isequal([nextsite,-AvsB]),looplist)
                        AvsB=-AvsB
                        prevdir=bonddir
                        push!(looplist,[nextsite,AvsB])
                        currsite=nextsite
                        jobdone=1
                        if newstart!=1
                            #println("Snipping out the first $(newstart-1) nodes.")
                            for k in 1:newstart-1
                                popfirst!(looplist)
                            end
                        end
                        return looplist
                    end
                else
                    #println("Wrong turn at [$currsite,$AvsB]")
                    extra+=1
                end
            end
        end
        
    end    
end

function findlongflippableloop(spinconfig,Nlist)
    # Nlist is of the form [N1,N2,N3]
    for i in 1:10000000
        # The outermost loop is to allow deadends.
        startover=0
        σAB=[0,0] #Initialize

        init=mod.(rand(Int,3),Nlist)+[1,1,1]
        AvsB=2*mod(rand(Int),2)-1 # 1 for A sublattice and -1 for B

        looplist=[[init,AvsB]]
        bonddir=mod(rand(Int),4)
        second=nextlatticesite(init,AvsB,bonddir,Nlist)
        AvsB=-AvsB
        prevdir=bonddir
        push!(looplist,[second,AvsB])
        currsite=second
        if AvsB==1
            σAB[1]=fetchspin(init,second,bonddir,1,spinconfig)#nA
            σAB[2]=-σAB[1]#nB
        elseif AvsB==-1
            σAB[2]=fetchspin(init,second,bonddir,-1,spinconfig)
            σAB[1]=-σAB[2]
        end #This block sets the requirement that any bond towards an A diamond site should have the pyrochlore site occupation for that bond to equal nA≡σAB[1]. Similarly for B.

        for it in 1:100000000
            #This is the loop through different bonds.
            if startover==1
                break
            end
            
            jobdone=0
            extra=0
            randint=mod(rand(Int),3)+1
            while jobdone==0
                if extra==3
                    #println("Deadend! Start over")
                    startover=1
                    break
                end
                bonddir=mod(prevdir+mod(randint-1+extra,3)+1,4)
                nextsite=nextlatticesite(currsite,AvsB,bonddir,Nlist)
                if fetchspin(currsite,nextsite,bonddir,-AvsB,spinconfig)==σAB[((AvsB+1)÷2)+1] # Ensures flippability.
                    if ([nextsite,-AvsB] in looplist)==false #Ensures non-visitedness
                        AvsB=-AvsB
                        prevdir=bonddir
                        push!(looplist,[nextsite,AvsB])
                        currsite=nextsite
                        jobdone=1
                    elseif [nextsite,-AvsB]==looplist[1]
                        AvsB=-AvsB
                        prevdir=bonddir
                        push!(looplist,[nextsite,AvsB])
                        currsite=nextsite
                        jobdone=1
                        return looplist
                    else
                        #println("Wrong turn at [$currsite,$AvsB] or reached a previously visited site")
                        extra+=1
                    end
                else
                    #println("Wrong turn at [$currsite,$AvsB] or reached a previously visited site")
                    extra+=1
                end
            end
        end 
    end    
end

function cartesiandisplacement(δl,a,b,Nlist)
    #In units of (2/√3)
    δl_new=[δl[1],δl[2],δl[3]] #Initialize
    for i in 1:3
        if δl_new[i]<-(Nlist[i]/2)
            δl_new[i]+=Nlist[i]
        elseif δl_new[i]>Nlist[i]/2
            δl_new[i]+=(-(Nlist[i]))
        end
    end
    displacement=[δl_new[2]+δl_new[3],δl_new[3]+δl_new[1],δl_new[1]+δl_new[2]]
    # if a==b do nothing
    if a==0 && b==1
        displacement+=[0,0.5,0.5]
    elseif a==0 && b==2
        displacement+=[0.5,0,0.5]
    elseif a==0 && b==3
        displacement+=[0.5,0.5,0]
    elseif a==1 && b==0
        displacement+=-[0,0.5,0.5]
    elseif a==2 && b==0
        displacement+=-[0.5,0,0.5]
    elseif a==3 && b==0
        displacement+=-[0.5,0.5,0]

    elseif a==1 && b==2
        displacement+=[0.5,-0.5,0]
    elseif a==2 && b==1
        displacement+=-[0.5,-0.5,0]
    elseif a==2 && b==3
        displacement+=[0,0.5,-0.5]
    elseif a==3 && b==2
        displacement+=-[0,0.5,-0.5]
    elseif a==3 && b==1
        displacement+=[-0.5,0,0.5]
    elseif a==1 && b==3
        displacement+=-[-0.5,0,0.5]
    end
    return displacement
end

function emuvector(mu)
    if mu==0
        return [0.5,0.5,0.5]
    elseif mu==1
        return [0.5,-0.5,-0.5]
    elseif mu==2
        return [-0.5,0.5,-0.5]
    elseif mu==3
        return [-0.5,-0.5,0.5]
    end
end

function vdwpotentialpbc(δl,a,b,Nlist,V,α=6) 
    #Returns the potential between [l,a]=[[l1,l2,l3],a] and [l',b]=[[l'1,l'2,l'3],b] such that l-l'=δl
    #Does not include nearest neighbour interactions
    displacement=cartesiandisplacement(δl,a,b,Nlist)
    distsquare=sum(abs2.(displacement))
    if distsquare<=0.5
        return 0
    else 
        return ((V/2)*(2*distsquare)^(-α/2))
    end
end

function vdwMatrixReal(Nlist,V,α=6)
    interactionmatrix=zeros(4,4,Nlist[1],Nlist[2],Nlist[3])
    for k in 1:Nlist[3]
        for j in 1:Nlist[2]
            for i in 1:Nlist[1]
                for a2 in 0:3
                    for a1 in 0:3
                        δl=[i-1,j-1,k-1] #Julia fft conventions
                        interactionmatrix[a1+1,a2+1,i,j,k]=vdwpotentialpbc(δl,a1,a2,Nlist,V,α)
                    end
                end
            end
        end
    end
    return interactionmatrix
end

function vdwMatrixFourier(Nlist,V,α=6)
    interactionmatrix=vdwMatrixReal(Nlist,V,α)
    return fft(interactionmatrix,[3,4,5])
end

function dipolepotentialpbc(δl,a,b,Nlist,D)
    displacement=cartesiandisplacement(δl,a,b,Nlist)
    distsquare=sum(abs2.(displacement))
    
    if distsquare<=0.5
        return 0
    else
        num=0
        en=0
        if a==b
            num+=1
        else
            num+=-1/3
        end
        num+=-3*(dot(emuvector(a),displacement)*dot(emuvector(b),displacement))/(distsquare*0.75)
        en=((D/2)*num)*(2*distsquare)^(-3/2)
        return en
    end
end

function dipoleMatrixReal(Nlist,D)
    interactionmatrix=zeros(4,4,Nlist[1],Nlist[2],Nlist[3])
    for k in 1:Nlist[3]
        for j in 1:Nlist[2]
            for i in 1:Nlist[1]
                for a2 in 0:3
                    for a1 in 0:3
                        δl=[i-1,j-1,k-1] #Julia fft conventions
                        interactionmatrix[a1+1,a2+1,i,j,k]=dipolepotentialpbc(δl,a1,a2,Nlist,D)
                    end
                end
            end
        end
    end
    return interactionmatrix
end

function dipoleMatrixFourier(Nlist,D)
    interactionmatrix=dipoleMatrixReal(Nlist,D)
    return fft(interactionmatrix,[3,4,5])
end

function energycost_lr(spinconfig,interactionfourier,Nlist)
    spinconfig_fft=(1/(sqrt(Nlist[1]*Nlist[2]*Nlist[3])))*fft(spinconfig,[2,3,4])
    FourierEnergy=0+0im
    #FourierEnergy+=(V/2)*icerulecost(spinconfig,Nlist)
    for k in 1:Nlist[3]
        for j in 1:Nlist[3]
            for i in 1:Nlist[3]
                FourierEnergy+=0.25*((spinconfig_fft[:,i,j,k])')*interactionfourier[:,:,i,j,k]*spinconfig_fft[:,i,j,k]
            end
        end
    end
    return real(FourierEnergy)
end

function energycost_sanitycheck_lr(spinconfig,interactionreal,Nlist)
    energy=0
    #energy+=(V/2)*icerulecost(spinconfig,Nlist)
    for ka in 1:Nlist[3]
        for ja in 1:Nlist[2]
            for ia in 1:Nlist[1]
                for a in 0:3
                    for kb in 1:Nlist[3]
                        for jb in 1:Nlist[2]
                            for ib in 1:Nlist[1]
                                for b in 0:3
                                    m1=mod(ia-ib,Nlist[1])+1
                                    m2=mod(ja-jb,Nlist[2])+1
                                    m3=mod(ka-kb,Nlist[3])+1
                                    energy+=0.25*interactionreal[a+1,b+1,m1,m2,m3]*spinconfig[a+1,ia,ja,ka]*spinconfig[b+1,ib,jb,kb]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return energy
end

function loopflip_config(spinconfig,looplist,Nlist)
    newconfig=copy(spinconfig)
    L=length(looplist)-1
    for i in 1:L
        lp=looplist[i][1]
        p=looplist[i][2]
        lq=looplist[mod(i,L)+1][1]
        q=looplist[mod(i,L)+1][2]

        if p==1
            if lp==lq
                newconfig[1,lp[1],lp[2],lp[3]]=-newconfig[1,lp[1],lp[2],lp[3]]
            elseif lp==[lq[1]%Nlist[1]+1,lq[2],lq[3]]
                newconfig[2,lp[1],lp[2],lp[3]]=-newconfig[2,lp[1],lp[2],lp[3]]
            elseif lp==[lq[1],lq[2]%Nlist[2]+1,lq[3]]
                newconfig[3,lp[1],lp[2],lp[3]]=-newconfig[3,lp[1],lp[2],lp[3]]
            elseif lp==[lq[1],lq[2],lq[3]%Nlist[3]+1]
                newconfig[4,lp[1],lp[2],lp[3]]=-newconfig[4,lp[1],lp[2],lp[3]]
            end
        elseif p==-1
            if lq==lp
                newconfig[1,lq[1],lq[2],lq[3]]=-newconfig[1,lq[1],lq[2],lq[3]]
            elseif lq==[lp[1]%Nlist[1]+1,lp[2],lp[3]]
                newconfig[2,lq[1],lq[2],lq[3]]=-newconfig[2,lq[1],lq[2],lq[3]]
            elseif lq==[lp[1],lp[2]%Nlist[2]+1,lp[3]]
                newconfig[3,lq[1],lq[2],lq[3]]=-newconfig[3,lq[1],lq[2],lq[3]]
            elseif lq==[lp[1],lp[2],lp[3]%Nlist[3]+1]
                newconfig[4,lq[1],lq[2],lq[3]]=-newconfig[4,lq[1],lq[2],lq[3]]
            end
        end
    end
    return newconfig
end

function siteflip_config(config,Nlist,Vnn)
    l=[mod(rand(Int),Nlist[1])+1,mod(rand(Int),Nlist[2])+1,mod(rand(Int),Nlist[3])+1]
    i=mod(rand(Int),4)+1
    lnext=nextlatticesite(l,1,i-1,Nlist)
    newconfig=copy(config)
    eninit=(Vnn/2)*((tet_charge(newconfig,[l,1],Nlist))^2 +(tet_charge(newconfig,[lnext,-1],Nlist))^2)

    newconfig[i,l[1],l[2],l[3]]=-newconfig[i,l[1],l[2],l[3]]

    enfinal=(Vnn/2)*((tet_charge(newconfig,[l,1],Nlist))^2 +(tet_charge(newconfig,[lnext,-1],Nlist))^2)
   
    return (newconfig,(enfinal-eninit))
end

function tet_charge(spinconfig,diamsite,Nlist)
    l=diamsite[1] #fcc site
    i=diamsite[2] #AvsB
    stot=0
    if i==1
        for v in 1:4
            stot+=spinconfig[v,l[1],l[2],l[3]]
        end
    elseif i==-1
        stot+=spinconfig[1,l[1],l[2],l[3]]
        stot+=spinconfig[2,l[1]%Nlist[1]+1,l[2],l[3]]
        stot+=spinconfig[3,l[1],l[2]%Nlist[2]+1,l[3]]
        stot+=spinconfig[4,l[1],l[2],l[3]%Nlist[3]+1]
    end
    return 0.5*stot
end

function icerulecost(spinconfig,Nlist)
    qtot=0
    for k in 1:Nlist[3]
        for j in 1:Nlist[2]
            for i in 1:Nlist[1]
                qa=tet_charge(spinconfig,[[i,j,k],1],Nlist)
                #=if qa!=0
                    println("Ice-rule violated at $([i,j,k]),1")
                end=#
                qb=tet_charge(spinconfig,[[i,j,k],-1],Nlist)
                #=if qb!=0
                    println("Ice-rule violated at $([i,j,k]),-1")
                end=#
                qtot+=qa^2 + qb^2
            end
        end
    end
    return qtot
end

function MCInitialize(config0,Nlist,runs,looponly)
    config=config0
    if runs==0
        return config
    end
    for i in 1:runs
        if i%2==1 || (i%2==0 && looponly==1)
            loop=findflippableloop(config,Nlist)
            config=loopflip_config(config,loop,Nlist)
        elseif (i%2==0 && looponly==0) 
            config=siteflip_config(config,Nlist,1.0)[1]
        end
        
    end
    return config
end

function LoopAnneal_linear(config0,Nlist,Vnn,interactionfourier,T,jumbleruns,annealruns,settleruns,interval,looponly=0,quick=1,α=6)
    #interactionfourier=vdwMatrixFourier(Nlist,V,α)
    en0=energycost_lr(config0,interactionfourier,Nlist)
    if jumbleruns!=0
        config=MCInitialize(config0,Nlist,jumbleruns,looponly)
    else
        config=config0
    end
    en_curr_lr=energycost_lr(config,interactionfourier,Nlist)

    en_ice=0.5*Vnn*icerulecost(config,Nlist)
    en_curr=en_curr_lr+en_ice
    if quick==0
        en_list=[en_curr]
    end
    
    delen_ice=0

    InitialT=T+(abs(en_curr-en0)*16)/(Nlist[1]*Nlist[2]*Nlist[3])
    T_dyn=InitialT
    for i in 1:(annealruns+settleruns)
        if i<=annealruns
            T_dyn = InitialT+(i-1)*(T-InitialT)/(annealruns)
        else
            if T==0
                if rand()<=0.2
                    T_dyn=InitialT
                else
                    T_dyn=T
                end
            else
                T_dyn=T
            end
        end
        if i%2==1 || (i%2==0 && looponly==1)
            #println("Now: lr_en=$en_curr_lr, ice_en1=$en_ice, ice_en2=$(0.5*V*icerulecost(config,Nlist)), en=$en_curr")
            if i%4==1
                loop=findlongflippableloop(config,Nlist)
            else
                loop=findflippableloop(config,Nlist)
            end
            newconfig=loopflip_config(config,loop,Nlist)
            en_new_lr=energycost_lr(newconfig,interactionfourier,Nlist)
            delen_ice=0
            #println("Found loop move: lr_new=$en_new_lr, delice=$delen_ice, en_ice=$(0.5*V*icerulecost(newconfig,Nlist))")

        elseif (i%2==0 && looponly==0)
            #println("Now: lr_en=$en_curr_lr, ice_en1=$en_ice, ice_en2=$(0.5*V*icerulecost(config,Nlist)), en=$en_curr")
            (newconfig,delen_ice)=siteflip_config(config,Nlist,Vnn)
            en_new_lr=energycost_lr(newconfig,interactionfourier,Nlist)
            #println("Found flip move: lr_new=$en_new_lr, delice=$delen_ice, en_ice=$(0.5*V*icerulecost(newconfig,Nlist))")
        end
        en_new=en_new_lr+en_ice+delen_ice

        accept=MCaccept(en_curr,en_new,T)
        if accept==1
            #println("Move accepted \n")
            #println("Move accepted. Old energy=$en_curr, old lr=$en_curr_lr, en_ice=$en_ice, new lr=$en_new_lr , New = $en_new, ice change = $delen_ice")
            en_curr=en_new
            en_ice+=delen_ice
            en_curr_lr=en_new_lr

            config=newconfig
        end
        if quick==0
            if i%interval==0
                push!(en_list,en_curr)
                #if i%(50*interval)==0
                    #println("$i steps done.")
                #end
            end
        end
    end
    if quick==0
        return (config,en_list)
    elseif quick==1
        return (config)
    end
end

function MCaccept(en_curr,en_new,T)
    if en_new<=en_curr
        accept=1
    else
        if T==0
            accept=0
        else
            p=exp(-(en_new-en_curr)/T)
            r=rand()
            if r<=p
                accept=1
            else
                accept=0
            end
        end
    end
end

function LoopMC(config0,Nlist,Vnn,interactionfourier,Temperature,jumbleruns,annealruns,settleruns,runs,interval,looponly=0,quick=1,α=6)
    #interactionfourier=vdwMatrixFourier(Nlist,V,α)
    
    if quick==1
        config=LoopAnneal_linear(config0,Nlist,Vnn,interactionfourier,Temperature,jumbleruns,annealruns,settleruns,interval,looponly,1,α)
    elseif quick==0
        config=LoopAnneal_linear(config0,Nlist,Vnn,interactionfourier,Temperature,jumbleruns,annealruns,settleruns,interval,looponly,0,α)[1]
    end
    en_curr_lr=energycost_lr(config,interactionfourier,Nlist)
    en_ice=0.5*Vnn*icerulecost(config,Nlist)
    en_curr=en_curr_lr+en_ice

    #en_list=[en_curr] 
    if quick==0
        en_list=[en_curr]
    end

    if quick==1
        Nsamp=1
        En_tot=en_curr
        EnSq_tot=en_curr^2
    end

    delen_ice=0
    accept=0 #Intiaialize
    for i in 1:runs
        if i%2==1 || (i%2==0 && looponly==1)
            if i%4==1
                loop=findlongflippableloop(config,Nlist)
            else
                loop=findflippableloop(config,Nlist)
            end
            newconfig=loopflip_config(config,loop,Nlist)
            en_new_lr=energycost_lr(newconfig,interactionfourier,Nlist)
            delen_ice=0

        elseif (i%2==0 && looponly==0)
            (newconfig,delen_ice)=siteflip_config(config,Nlist,Vnn)
            en_new_lr=energycost_lr(newconfig,interactionfourier,Nlist)
        end
        en_new=en_new_lr+en_ice+delen_ice

        accept=MCaccept(en_curr,en_new,Temperature)
        if accept==1
            en_curr=en_new
            en_ice+=delen_ice
            en_curr_lr=en_new_lr

            config=newconfig
        end
        if i%interval==0
            if quick==0
                push!(en_list,en_curr)
                if i%(50*interval)==0
                    #println("$i steps done.")
                end
            elseif quick==1
                Nsamp+=1
                En_tot+=en_curr
                EnSq_tot+=en_curr^2
                #push!(en_list,en_curr) 
            end 
        end
    end
    if quick==0
        return (config,en_list)
    elseif quick==1
        Emean=En_tot/Nsamp
        Cv=(EnSq_tot/Nsamp - (Emean)^2)/(Temperature^2)
        return (Cv,Emean,config)
    end
    
end

function mean(arr)
    return sum(arr)/length(arr)
end

function variance(arr)
    length(arr)=Lab
    return mean(arr.^2)-(mean(arr))^2
end

function temploop(config0,Nlist,V,Vnn,Tmin,Tmax,Tmid,coarse,fine,jumbleruns,annealruns,settleruns,runs,interval,looponly=1,α=6)
    interactionfourier=vdwMatrixFourier(Nlist,V);
    filename="Cv_size"*string(Nlist[1])*"_"*string(today())*".txt"
    open(filename, "w") do io
        write(io,"# $Nlist,V=$V, Vnn=$Vnn, jumbleruns=$jumbleruns, annealruns=$annealruns, settleruns=$settleruns, runs=$runs, interval=$interval, looponly=$looponly \n")
    
        @threads for i in 1:(fine+coarse)
                if i<=fine 
                    T=Tmin+(i-1)*(Tmid-Tmin)/(fine-1)
                else
                    T=Tmid+(i-fine)*(Tmax-Tmid)/coarse
                end
            
                (Cv,Emean,config)=LoopMC(config0,Nlist,Vnn,interactionfourier,T,jumbleruns,annealruns,settleruns,runs,interval,1,1)
                writedlm(io, [T Cv])
        end
    end
end

function makeTarray(Tmin,Tmid,Tmax,n1,n2,n3=14)
    ntot=n1+n2+n3
    betamin=1/Tmax
    betamid=1/Tmid
    betamax=1/Tmin
    betaarray=zeros(ntot)
    Tarray=zeros(ntot)
    for i in 1:n3
        betaarray[i]=betamin*(1.2)^(i-1)
        #println("$(1/betaarray[i])")
    end
    for i in 1:n2
        betaarray[n3+i]=betaarray[n3]+((betamid-betaarray[n3])*i/n2)
        #println("$(1/betaarray[i+n3])")
    end
    for i in 1:n1
        betaarray[n3+n2+i]=1/(Tmid-((Tmid-Tmin)*i/n1))
        #println("$(1/betaarray[i+n3+n2])")
    end
    for i in 1:ntot
        Tarray[ntot+1-i]=1/betaarray[i]
    end
    return Tarray
end


function SeqTempSweepMC(config0,interactionfourier,Nlist,Vnn,Tarray,settleruns,runs,interval,looponly=0)
    NTemp=length(Tarray)
    config=config0
    filename="Seq_Cv_size"*string(Nlist[1])*"_"*string(today())*".txt"
    open(filename, "w") do io
        write(io,"# $Nlist, Vnn=$Vnn, settleruns=$settleruns, runs=$runs, interval=$interval, looponly=$looponly \n")
    end
        for i in 1:NTemp
            T=Tarray[NTemp+1-i]  
            (Cv,Emean,config)=LoopMC(config,Nlist,Vnn,interactionfourier,T,0,0,settleruns,runs,interval,looponly,1)
            open(filename, "a") do io
                writedlm(io, [T Cv])
            end
            println("$i: T=$T, icerulecost=$(icerulecost(config,Nlist)), Cv=$Cv, Emean=$Emean")
        end
    return config
end

function countflippablehex(config,Nlist)
    count=0
    for k in 1:Nlist[3]
        for j in Nlist[2]
            for i in Nlist[1]
                looplistset=fetchhexplaquette([i,j,k],Nlist)
                for a in 1:4
                    if checkflippable(looplistset[a],config,Nlist)==1
                        count+=1
                    end
                end
            end
        end
    end
    return count
end

function loadfourierint(filename)
    #filename="4X8X8X8InteractionFourier.txt"
    #interactionfourier=zeros(Complex{Float64},4,4,8,8,8)
    interactionfourierflat=readdlm(filename, '\t', Complex{Float64}, '\n',comments=true)
    interactionfourier=reshape(interactionfourierflat,(4,4,8,8,8))
    return interactionfourier
end

#April 11 2023 New function
#function monopolegauge

#function RKplaquettecorr(fccsite1,orientation1,fccsite2,orientation2)