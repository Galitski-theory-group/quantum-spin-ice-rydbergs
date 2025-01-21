function W3W4vdwcalc(config,Nlist,intfourier)
    en_W3=0
    en_W4=0
    ki_1=0
    ki_2=0
    ki_3=0 #initialization
    config_fft=(1/(sqrt(Nlist[1]*Nlist[2]*Nlist[3])))*fft(config,[2,3,4])
    for k1_3 in 1:Nlist[3]
        for k1_2 in 1:Nlist[2]
            for k1_1 in 1:Nlist[1]
                for k2_3 in 1:Nlist[3]
                    for k2_2 in 1:Nlist[2]
                        for k2_1 in 1:Nlist[1]
                            for k3_3 in 1:Nlist[3]
                                for k3_2 in 1:Nlist[2]
                                    for k3_1 in 1:Nlist[1]
                                        for a_k1 in 0:3
                                            for a_k2 in 0:3
                                                for a_k3 in 0:3
                                                    for a_i in 0:3
                                                        ki_1=mod(-(k1_1+k2_1+k3_1-3),Nlist[1])+1
                                                        ki_2=mod(-(k1_2+k2_2+k3_2-3),Nlist[2])+1
                                                        ki_3=mod(-(k1_3+k2_3+k3_3-3),Nlist[3])+1
                                                        en_W3+=intfourier[a_i+1,a_k1+1,k1_1,k1_2,k1_3]*intfourier[a_i+1,a_k2+1,k2_1,k2_2,k2_3]*intfourier[a_i+1,a_k3+1,k3_1,k3_2,k3_3]*config_fft[a_k1+1,k1_1,k1_2,k1_3]*config_fft[a_k2+1,k2_1,k2_2,k2_3]*config_fft[a_k3+1,k3_1,k3_2,k3_3]*config_fft[a_i+1,ki_1,ki_2,ki_3]
                                                        for a_k4 in 0:3
                                                            en_W4+=intfourier[a_i+1,a_k1+1,k1_1,k1_2,k1_3]*intfourier[a_i+1,a_k2+1,k2_1,k2_2,k2_3]*intfourier[a_i+1,a_k3+1,k3_1,k3_2,k3_3]*intfourier[a_i+1,a_k4+1,ki_1,ki_2,ki_3]*config_fft[a_k1+1,k1_1,k1_2,k1_3]*config_fft[a_k2+1,k2_1,k2_2,k2_3]*config_fft[a_k3+1,k3_1,k3_2,k3_3]*config_fft[a_i+1,ki_1,ki_2,ki_3]
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return ([real(en_W3)/32,real(en_W4)/64]./(Nlist[1]*Nlist[2]*Nlist[3]))
end

function RK_MC_with_W3W4(Nlist,Nsweep,interval)
    #The load function already assumes Nlist=[8,8,8]
    configcurr=initializebycopying(Nlist)
    HLRfourier=loadfourierint("4X8X8X8InteractionFourier.txt")
    W2fourier=loadfourierint("4X8X8X8W2Fourier.txt")
    
    nflipcount=0
    HLRtot=0
    W2tot=0
    W3tot=0
    W4tot=0

    nfliplist=zeros(Nsweep÷interval)
    enlist=zeros(Nsweep÷interval)
    W2list=zeros(Nsweep÷interval)
    W3list=zeros(Nsweep÷interval)
    W4list=zeros(Nsweep÷interval)

    #println("1: en=$HLRtot")

    firstsite=[1,1,1,1] #initialize
    loopbag=findshortflippableloop(configcurr,Nlist) #initialize
    Ndiamsites=Nlist[1]*Nlist[2]*Nlist[3]
    cond=0
    for i in 1:Nsweep
        for j in 1:Ndiamsites
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

        nflipcount+=CountFlippableHex(configcurr,Nlist)
        HLRtot+=energycost_lr(configcurr,HLRfourier,Nlist)
        W2tot+=energycost_lr(configcurr,W2fourier,Nlist)
        if i%interval==0
            res=W3W4vdwcalc(configcurr,Nlist,HLRfourier)
            W3tot+=res[1]
            W4tot+=res[2]
            println("$(i): HLRavg=$(HLRtot/i), W2avg=$(W2tot/i), W3avg=$(W3tot/(i÷interval)), W4avg=$(W4tot/(i÷interval))")
        end
    end
    nflip=(nflipcount/(Nlist[1]*Nlist[2]*Nlist[3]*4))/Nsweep
    HLRavg=HLRtot/(Nsweep)
    W2avg=W2tot/Nsweep
    W3avg=W3tot/(Nsweep÷interval)
    W4avg=W4tot/(Nsweep÷interval)
    return [nflip,HLRavg,W2avg,W3avg,W4avg]
end

function RK_MC_with_L2M2(Nlist,Nsweep,interval)
    #The load function already assumes Nlist=[8,8,8]
    configcurr=initializebycopying(Nlist)
    HLRfourier=loadfourierint("4X8X8X8InteractionFourier.txt")
    L2fourier=loadfourierint("4X8X8X8L2vdWFourier.txt")
    
    L2tot=0
    M2tot=0

    #println("1: en=$HLRtot")

    firstsite=[1,1,1,1] #initialize
    loopbag=findshortflippableloop(configcurr,Nlist) #initialize
    Ndiamsites=Nlist[1]*Nlist[2]*Nlist[3]
    cond=0
    for i in 1:Nsweep
        for j in 1:Ndiamsites
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

        L2tot+=energycost_lr(configcurr,L2fourier,Nlist)
        if i%interval==0
            M2tot+=M2vdwcalc(configcurr,Nlist,HLRfourier)
            println("$(i): L2avg=$(L2tot/i), M2avg=$(M2tot/(i÷interval))")
        end
    end
    L2avg=L2tot/Nsweep
    M2avg=M2tot/(Nsweep÷interval)
    return [L2avg,M2avg]
end

function neighbdata(d,μ1)
    μ2=mod(μ1+(mod(d-1,3)+1),4)
    x=[0,0,0] #Initialize
    if d<=3
        x=[0,0,0]
    end
    if d>3
        if μ2!=0
            x[μ2]+=1
        end
        if μ1!=0
            x[μ1]+=-1
        end
    end
    return [x,μ2]
end

function M2vdwcalc(config,Nlist,intfourier)
    en_M2=0
    
    phase=0.86602540378+0.5im; ang=0; #initialize
    k5=[0,0,0.1] #initialize
    config_fft=(1/(sqrt(Nlist[1]*Nlist[2]*Nlist[3])))*fft(config,[2,3,4])
    for k1_3 in 1:Nlist[3]
        for k1_2 in 1:Nlist[2]
            for k1_1 in 1:Nlist[1]
                for k2_3 in 1:Nlist[3]
                    for k2_2 in 1:Nlist[2]
                        for k2_1 in 1:Nlist[1]
                            for d in 1:6
                                for a_k1 in 0:3
                                    x=[0,0,0] #initialize
                                    a_k4=mod(a_k1+(mod(d-1,3)+1),4)
                                    if d<=3
                                        x=[0,0,0]
                                    end
                                    if d>3
                                        if a_k4!=0
                                            x[a_k4]+=1
                                        end
                                        if a_k1!=0
                                            x[a_k1]+=-1
                                        end
                                    end
                                    k5=-2*π*[(k1_1+k2_1-2)/Nlist[1],(k1_2+k2_2-2)/Nlist[2],(k1_3+k2_3-2)/Nlist[3]]
                                    phase=exp(1im*dot(x,k5))
                                    for k3_3 in 1:Nlist[3]
                                        for k3_2 in 1:Nlist[2]
                                            for k3_1 in 1:Nlist[1]
                                                for a_k2 in 0:3
                                                    for a_k3 in 0:3
                                                        k4_1=mod(-(k1_1+k2_1+k3_1-3),Nlist[1])+1
                                                        k4_2=mod(-(k1_2+k2_2+k3_2-3),Nlist[2])+1
                                                        k4_3=mod(-(k1_3+k2_3+k3_3-3),Nlist[3])+1
                                                        en_M2+=real(phase*intfourier[a_k1+1,a_k2+1,k2_1,k2_2,k2_3]*intfourier[a_k4+1,a_k3+1,k3_1,k3_2,k3_3]*config_fft[a_k1+1,k1_1,k1_2,k1_3]*config_fft[a_k2+1,k2_1,k2_2,k2_3]*config_fft[a_k4+1,k4_1,k4_2,k4_3]*config_fft[a_k3+1,k3_1,k3_2,k3_3])
                                                        #en_M2+=real(intfourier[a_k1+1,a_k2+1,k2_1,k2_2,k2_3]*intfourier[a_k4+1,a_k3+1,k3_1,k3_2,k3_3]*config_fft[a_k1+1,k1_1,k1_2,k1_3]*config_fft[a_k2+1,k2_1,k2_2,k2_3]*config_fft[a_k4+1,k4_1,k4_2,k4_3]*config_fft[a_k3+1,k3_1,k3_2,k3_3])
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return en_M2/(Nlist[1]*Nlist[2]*Nlist[3])
end
#=
function testM2vdwcalc(config,Nlist,intfourier)
    en_M2=0
    res=[[0,0,0],0] #initialize
    config_fft=(1/(sqrt(Nlist[1]*Nlist[2]*Nlist[3])))*fft(config,[2,3,4])
    for k1_3 in 1:Nlist[3]
        for k1_2 in 1:Nlist[2]
            for k1_1 in 1:Nlist[1]
                for k2_3 in 1:Nlist[3]
                    for k2_2 in 1:Nlist[2]
                        for k2_1 in 1:Nlist[1]
                            for k3_3 in 1:Nlist[3]
                                for k3_2 in 1:Nlist[2]
                                    for k3_1 in 1:Nlist[1]
                                        for a_k1 in 0:3
                                            for a_k2 in 0:3
                                                for a_k3 in 0:3
                                                    for d in 1:6
                                                        res=neighbdata(d,a_k1)
                                                        x=res[1]
                                                        a_k4=res[2]
                                                        k5=-2*π*[(k1_1+k2_1-2)/Nlist[1],(k1_2+k2_2-2)/Nlist[2],(k1_3+k2_3-2)/Nlist[3]]
                                                        phase=exp(1im*dot(x,k5))
                                                        k4_1=mod(-(k1_1+k2_1+k3_1-3),Nlist[1])+1
                                                        k4_2=mod(-(k1_2+k2_2+k3_2-3),Nlist[2])+1
                                                        k4_3=mod(-(k1_3+k2_3+k3_3-3),Nlist[3])+1
                                                        en_M2+=phase*intfourier[a_k1+1,a_k2+1,k2_1,k2_2,k2_3]*intfourier[a_k4+1,a_k3+1,k3_1,k3_2,k3_3]*config_fft[a_k1+1,k1_1,k1_2,k1_3]*config_fft[a_k2+1,k2_1,k2_2,k2_3]*config_fft[a_k4+1,k4_1,k4_2,k4_3]*config_fft[a_k3+1,k3_1,k3_2,k3_3]
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return en_M2/(Nlist[1]*Nlist[2]*Nlist[3])
end
=#
function L2vdwRealSpace(Nlist,vdwintmatrix)
    L2vdwmatrix=zeros(4,4,Nlist[1],Nlist[2],Nlist[3])
    for k in 1:Nlist[3]
        for j in 1:Nlist[2]
            for i in 1:Nlist[1]
                for a1 in 0:3
                    for a2 in 0:3
                        for n in 1:Nlist[3]
                            for m in 1:Nlist[2]
                                for l in 1:Nlist[1]
                                    for b1 in 0:3
                                        for d in 1:6
                                            res=neighbdata(d,b1)
                                            x=res[1]
                                            b2=res[2]

                                            p1=mod(l-i,Nlist[1])+1
                                            q1=mod(m-j,Nlist[2])+1
                                            r1=mod(n-k,Nlist[3])+1
                                            p2=mod(l-1+x[1],Nlist[1])+1
                                            q2=mod(m-1+x[2],Nlist[2])+1
                                            r2=mod(n-1+x[3],Nlist[3])+1
                                            L2vdwmatrix[a1+1,a2+1,i,j,k]+=vdwintmatrix[b1+1,a1+1,p1,q1,r1]*vdwintmatrix[b2+1,a2+1,p2,q2,r2]
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    println("Real space done")
    return L2vdwmatrix
end