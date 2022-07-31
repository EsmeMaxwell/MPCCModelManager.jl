
# TODO rename testdense



function mm_spec_testdense_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(5, 2, 3, 2, 3, 2, 0)
    return dimspec
end



function mm_spec_testdense_defn(x, pr, ps)

    f = x[1]^2 * x[2]^3 * x[3]^4 * x[4]^5 * x[5]^6 * ( pr[1] + pr[2] )^2
            + pr[1]^2 * pr[2]^3 * exp(-2*x[1] - 3*x[2] - 5*x[3] - 7*x[4] - 11*x[5])

    ce = Vector{Num}(undef, 2)
    ce[1] = ( x[1] + x[2] + x[3] + x[4] + x[5] )^2 * sin(pr[1]) * cos(pr[2])
    ce[2] = pr[1] * ( x[1] + x[2] )^3 + pr[2] * ( x[3] + x[4] + x[5] )^2

    ci = Vector{Num}(undef, 3)
    ci[1] = x[1] + x[2] * ( pr[1] + pr[2] )^2
    ci[2] = x[3] + x[4] * x[5] * ( pr[1] - pr[2] )^2
    ci[3] = ( x[2] - x[1] ) * ( x[4] - x[3] ) * x[5]^2 * ( pr[2] - pr[1] )^2

    F = Matrix{Num}(undef, 3, 2)
    F[1,1] = ( x[1] + x[2] )^2
    F[1,2] = ( x[3] + x[4] + x[5] ) * ( pr[1] - pr[2] )^2
    F[2,1] = ( x[3] + x[4] ) * ( pr[1] - pr[2] )
    F[2,2] = ( x[4] - x[1] )
    F[3,1] = ( 1 + x[5] )^7
    F[3,2] = ( 1 + pr[1] + pr[2] ) * sin(x[1])

    return MPCCDefinition(f, ce, ci, F)
end


function mm_spec_testdense_nzmask()::MPCCModelNZMask
    nz_mask = MPCCModelNZMask(
        # ce
        [],
        # ci
        [],
        # F
        [ 
            # F1
            [ ],
            # F2
            [ ] 
        ] )    
    return nz_mask
end



function mm_spec_testdense_testvectors()

    x_val = Vector{Rational}(undef, 5)
    x_val = [4//3, 27//10, 7//2, 1//13, 5]

    pr_val = Vector{Rational}(undef, 2)
    pr_val = [1//7, 23//10]

    ps_val = Int64[]

    # f:
    f_chk_eval = (12167*exp(-16342/195))/49000 + 15667766415/11881376

    # ce:
    ce_chk_eval = Vector{Float64}(undef, 2)
    ce_chk_eval[1] = (6046681*cos(23/10)*sin(1/7))/38025
    ce_chk_eval[2] = 1425923221/7985250

    # ci:
    ci_chk_eval = Vector{Float64}(undef, 3)
    ci_chk_eval[1] = 2564521/147000
    ci_chk_eval[2] = 67391/12740
    ci_chk_eval[3] = -83200849/152880

    # F:
    F_chk_eval = Matrix{Float64}(undef, 3, 2)
    F_chk_eval[1, 1] = 14641/900
    F_chk_eval[1, 2] = 5084623/127400
    F_chk_eval[2, 1] = -14043/1820
    F_chk_eval[2, 2] = -49/39
    F_chk_eval[3, 1] = 279936
    F_chk_eval[3, 2] = (241*sin(4/3))/70



    # gradf:
    gradf_chk_eval = Vector{Float64}(undef, 5)
    gradf_chk_eval[1] = 47003299245/23762752 - (12167*exp(-16342/195))/24500
    gradf_chk_eval[2] = 8704314675/5940688 - (36501*exp(-16342/195))/49000
    gradf_chk_eval[3] = 2238252345/1485172 - (12167*exp(-16342/195))/9800
    gradf_chk_eval[4] = 78338832075/913952 - (12167*exp(-16342/195))/7000
    gradf_chk_eval[5] = 9400659849/5940688 - (133837*exp(-16342/195))/49000

    # ce:
    gradce_chk_eval = Matrix{Float64}(undef, 5, 2)
    gradce_chk_eval[1, 1] = (4918*cos(23/10)*sin(1/7))/195
    gradce_chk_eval[1, 2] = 14641/2100
    gradce_chk_eval[2, 1] = (4918*cos(23/10)*sin(1/7))/195
    gradce_chk_eval[2, 2] = 14641/2100
    gradce_chk_eval[3, 1] = (4918*cos(23/10)*sin(1/7))/195
    gradce_chk_eval[3, 2] = 5129/130
    gradce_chk_eval[4, 1] = (4918*cos(23/10)*sin(1/7))/195
    gradce_chk_eval[4, 2] = 5129/130
    gradce_chk_eval[5, 1] = (4918*cos(23/10)*sin(1/7))/195
    gradce_chk_eval[5, 2] = 5129/130

    # gradci:
    gradci_chk_eval = Matrix{Float64}(undef, 5, 3)
    gradci_chk_eval[1, 1] = 1
    gradci_chk_eval[1, 2] = 0
    gradci_chk_eval[1, 3] = 2029289/5096
    gradci_chk_eval[2, 1] = 29241/4900
    gradci_chk_eval[2, 2] = 0
    gradci_chk_eval[2, 3] = -2029289/5096
    gradci_chk_eval[3, 1] = 0
    gradci_chk_eval[3, 2] = 1
    gradci_chk_eval[3, 3] = -934841/5880
    gradci_chk_eval[4, 1] = 0
    gradci_chk_eval[4, 2] = 22801/980
    gradci_chk_eval[4, 3] = 934841/5880
    gradci_chk_eval[5, 1] = 0
    gradci_chk_eval[5, 2] = 22801/63700
    gradci_chk_eval[5, 3] = -83200849/382200

    # gradF:
    gradF_chk_eval = Matrix{Vector{Float64}}(undef, 3, 2)
    gradF_chk_eval[1, 1] = Vector{Float64}(undef, 5)
    gradF_chk_eval[1, 1][1] = 121/15
    gradF_chk_eval[1, 1][2] = 121/15
    gradF_chk_eval[1, 1][3] = 0
    gradF_chk_eval[1, 1][4] = 0
    gradF_chk_eval[1, 1][5] = 0
    gradF_chk_eval[1, 2] = Vector{Float64}(undef, 5)
    gradF_chk_eval[1, 2][1] = 0
    gradF_chk_eval[1, 2][2] = 0
    gradF_chk_eval[1, 2][3] = 22801/4900
    gradF_chk_eval[1, 2][4] = 22801/4900
    gradF_chk_eval[1, 2][5] = 22801/4900
    gradF_chk_eval[2, 1] = Vector{Float64}(undef, 5)
    gradF_chk_eval[2, 1][1] = 0
    gradF_chk_eval[2, 1][2] = 0
    gradF_chk_eval[2, 1][3] = -151/70
    gradF_chk_eval[2, 1][4] = -151/70
    gradF_chk_eval[2, 1][5] = 0
    gradF_chk_eval[2, 2] = Vector{Float64}(undef, 5)
    gradF_chk_eval[2, 2][1] = -1
    gradF_chk_eval[2, 2][2] = 0
    gradF_chk_eval[2, 2][3] = 0
    gradF_chk_eval[2, 2][4] = 1
    gradF_chk_eval[2, 2][5] = 0
    gradF_chk_eval[3, 1] = Vector{Float64}(undef, 5)
    gradF_chk_eval[3, 1][1] = 0
    gradF_chk_eval[3, 1][2] = 0
    gradF_chk_eval[3, 1][3] = 0
    gradF_chk_eval[3, 1][4] = 0
    gradF_chk_eval[3, 1][5] = 326592
    gradF_chk_eval[3, 2] = Vector{Float64}(undef, 5)
    gradF_chk_eval[3, 2][1] = (241*cos(4/3))/70
    gradF_chk_eval[3, 2][2] = 0
    gradF_chk_eval[3, 2][3] = 0
    gradF_chk_eval[3, 2][4] = 0
    gradF_chk_eval[3, 2][5] = 0



    # hessf:
    hessf_chk_eval = Matrix{Float64}(undef, 5, 5)
    hessf_chk_eval[1, 1] = (12167*exp(-16342/195))/12250 + 141009897735/95051008
    hessf_chk_eval[1, 2] = (36501*exp(-16342/195))/24500 + 26112944025/11881376
    hessf_chk_eval[1, 3] = (12167*exp(-16342/195))/4900 + 6714757035/2970344
    hessf_chk_eval[1, 4] = (12167*exp(-16342/195))/3500 + 235016496225/1827904
    hessf_chk_eval[1, 5] = (133837*exp(-16342/195))/24500 + 28201979547/11881376
    hessf_chk_eval[2, 1] = (36501*exp(-16342/195))/24500 + 26112944025/11881376
    hessf_chk_eval[2, 2] = (109503*exp(-16342/195))/49000 + 1611910125/1485172
    hessf_chk_eval[2, 3] = (36501*exp(-16342/195))/9800 + 1243473525/742586
    hessf_chk_eval[2, 4] = (36501*exp(-16342/195))/7000 + 43521573375/456976
    hessf_chk_eval[2, 5] = (401511*exp(-16342/195))/49000 + 5222588805/2970344
    hessf_chk_eval[3, 1] = (12167*exp(-16342/195))/4900 + 6714757035/2970344
    hessf_chk_eval[3, 2] = (36501*exp(-16342/195))/9800 + 1243473525/742586
    hessf_chk_eval[3, 3] = (12167*exp(-16342/195))/1960 + 959251005/742586
    hessf_chk_eval[3, 4] = (12167*exp(-16342/195))/1400 + 11191261725/114244
    hessf_chk_eval[3, 5] = (133837*exp(-16342/195))/9800 + 1342951407/742586
    hessf_chk_eval[4, 1] = (12167*exp(-16342/195))/3500 + 235016496225/1827904
    hessf_chk_eval[4, 2] = (36501*exp(-16342/195))/7000 + 43521573375/456976
    hessf_chk_eval[4, 3] = (12167*exp(-16342/195))/1400 + 11191261725/114244
    hessf_chk_eval[4, 4] = (12167*exp(-16342/195))/1000 + 78338832075/17576
    hessf_chk_eval[4, 5] = (133837*exp(-16342/195))/7000 + 47003299245/456976
    hessf_chk_eval[5, 1] = (133837*exp(-16342/195))/24500 + 28201979547/11881376
    hessf_chk_eval[5, 2] = (401511*exp(-16342/195))/49000 + 5222588805/2970344
    hessf_chk_eval[5, 3] = (133837*exp(-16342/195))/9800 + 1342951407/742586
    hessf_chk_eval[5, 4] = (133837*exp(-16342/195))/7000 + 47003299245/456976
    hessf_chk_eval[5, 5] = (1472207*exp(-16342/195))/49000 + 9400659849/5940688

    # hessce:
    hessce_chk_eval = Vector{Matrix{Float64}}(undef, 2)
    hessce_chk_eval[1] = Matrix{Float64}(undef, 5, 5)
    hessce_chk_eval[1][1, 1] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][1, 2] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][1, 3] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][1, 4] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][1, 5] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][2, 1] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][2, 2] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][2, 3] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][2, 4] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][2, 5] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][3, 1] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][3, 2] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][3, 3] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][3, 4] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][3, 5] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][4, 1] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][4, 2] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][4, 3] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][4, 4] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][4, 5] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][5, 1] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][5, 2] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][5, 3] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][5, 4] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[1][5, 5] = 2*cos(23/10)*sin(1/7)
    hessce_chk_eval[2] = Matrix{Float64}(undef, 5, 5)
    hessce_chk_eval[2][1, 1] = 121/35
    hessce_chk_eval[2][1, 2] = 121/35
    hessce_chk_eval[2][1, 3] = 0
    hessce_chk_eval[2][1, 4] = 0
    hessce_chk_eval[2][1, 5] = 0
    hessce_chk_eval[2][2, 1] = 121/35
    hessce_chk_eval[2][2, 2] = 121/35
    hessce_chk_eval[2][2, 3] = 0
    hessce_chk_eval[2][2, 4] = 0
    hessce_chk_eval[2][2, 5] = 0
    hessce_chk_eval[2][3, 1] = 0
    hessce_chk_eval[2][3, 2] = 0
    hessce_chk_eval[2][3, 3] = 23/5
    hessce_chk_eval[2][3, 4] = 23/5
    hessce_chk_eval[2][3, 5] = 23/5
    hessce_chk_eval[2][4, 1] = 0
    hessce_chk_eval[2][4, 2] = 0
    hessce_chk_eval[2][4, 3] = 23/5
    hessce_chk_eval[2][4, 4] = 23/5
    hessce_chk_eval[2][4, 5] = 23/5
    hessce_chk_eval[2][5, 1] = 0
    hessce_chk_eval[2][5, 2] = 0
    hessce_chk_eval[2][5, 3] = 23/5
    hessce_chk_eval[2][5, 4] = 23/5
    hessce_chk_eval[2][5, 5] = 23/5

    # hessci:
    hessci_chk_eval = Vector{Matrix{Float64}}(undef, 3)
    hessci_chk_eval[1] = Matrix{Float64}(undef, 5, 5)
    hessci_chk_eval[1][1, 1] = 0
    hessci_chk_eval[1][1, 2] = 0
    hessci_chk_eval[1][1, 3] = 0
    hessci_chk_eval[1][1, 4] = 0
    hessci_chk_eval[1][1, 5] = 0
    hessci_chk_eval[1][2, 1] = 0
    hessci_chk_eval[1][2, 2] = 0
    hessci_chk_eval[1][2, 3] = 0
    hessci_chk_eval[1][2, 4] = 0
    hessci_chk_eval[1][2, 5] = 0
    hessci_chk_eval[1][3, 1] = 0
    hessci_chk_eval[1][3, 2] = 0
    hessci_chk_eval[1][3, 3] = 0
    hessci_chk_eval[1][3, 4] = 0
    hessci_chk_eval[1][3, 5] = 0
    hessci_chk_eval[1][4, 1] = 0
    hessci_chk_eval[1][4, 2] = 0
    hessci_chk_eval[1][4, 3] = 0
    hessci_chk_eval[1][4, 4] = 0
    hessci_chk_eval[1][4, 5] = 0
    hessci_chk_eval[1][5, 1] = 0
    hessci_chk_eval[1][5, 2] = 0
    hessci_chk_eval[1][5, 3] = 0
    hessci_chk_eval[1][5, 4] = 0
    hessci_chk_eval[1][5, 5] = 0
    hessci_chk_eval[2] = Matrix{Float64}(undef, 5, 5)
    hessci_chk_eval[2][1, 1] = 0
    hessci_chk_eval[2][1, 2] = 0
    hessci_chk_eval[2][1, 3] = 0
    hessci_chk_eval[2][1, 4] = 0
    hessci_chk_eval[2][1, 5] = 0
    hessci_chk_eval[2][2, 1] = 0
    hessci_chk_eval[2][2, 2] = 0
    hessci_chk_eval[2][2, 3] = 0
    hessci_chk_eval[2][2, 4] = 0
    hessci_chk_eval[2][2, 5] = 0
    hessci_chk_eval[2][3, 1] = 0
    hessci_chk_eval[2][3, 2] = 0
    hessci_chk_eval[2][3, 3] = 0
    hessci_chk_eval[2][3, 4] = 0
    hessci_chk_eval[2][3, 5] = 0
    hessci_chk_eval[2][4, 1] = 0
    hessci_chk_eval[2][4, 2] = 0
    hessci_chk_eval[2][4, 3] = 0
    hessci_chk_eval[2][4, 4] = 0
    hessci_chk_eval[2][4, 5] = 22801/4900
    hessci_chk_eval[2][5, 1] = 0
    hessci_chk_eval[2][5, 2] = 0
    hessci_chk_eval[2][5, 3] = 0
    hessci_chk_eval[2][5, 4] = 22801/4900
    hessci_chk_eval[2][5, 5] = 0
    hessci_chk_eval[3] = Matrix{Float64}(undef, 5, 5)
    hessci_chk_eval[3][1, 1] = 0
    hessci_chk_eval[3][1, 2] = 0
    hessci_chk_eval[3][1, 3] = 22801/196
    hessci_chk_eval[3][1, 4] = -22801/196
    hessci_chk_eval[3][1, 5] = 2029289/12740
    hessci_chk_eval[3][2, 1] = 0
    hessci_chk_eval[3][2, 2] = 0
    hessci_chk_eval[3][2, 3] = -22801/196
    hessci_chk_eval[3][2, 4] = 22801/196
    hessci_chk_eval[3][2, 5] = -2029289/12740
    hessci_chk_eval[3][3, 1] = 22801/196
    hessci_chk_eval[3][3, 2] = -22801/196
    hessci_chk_eval[3][3, 3] = 0
    hessci_chk_eval[3][3, 4] = 0
    hessci_chk_eval[3][3, 5] = -934841/14700
    hessci_chk_eval[3][4, 1] = -22801/196
    hessci_chk_eval[3][4, 2] = 22801/196
    hessci_chk_eval[3][4, 3] = 0
    hessci_chk_eval[3][4, 4] = 0
    hessci_chk_eval[3][4, 5] = 934841/14700
    hessci_chk_eval[3][5, 1] = 2029289/12740
    hessci_chk_eval[3][5, 2] = -2029289/12740
    hessci_chk_eval[3][5, 3] = -934841/14700
    hessci_chk_eval[3][5, 4] = 934841/14700
    hessci_chk_eval[3][5, 5] = -83200849/1911000

    # hessF:
    hessF_chk_eval = Matrix{Matrix{Float64}}(undef, 3, 2)
    hessF_chk_eval[1, 1] = Matrix{Float64}(undef, 5, 5)
    hessF_chk_eval[1, 1][1, 1] = 2
    hessF_chk_eval[1, 1][1, 2] = 2
    hessF_chk_eval[1, 1][1, 3] = 0
    hessF_chk_eval[1, 1][1, 4] = 0
    hessF_chk_eval[1, 1][1, 5] = 0
    hessF_chk_eval[1, 1][2, 1] = 2
    hessF_chk_eval[1, 1][2, 2] = 2
    hessF_chk_eval[1, 1][2, 3] = 0
    hessF_chk_eval[1, 1][2, 4] = 0
    hessF_chk_eval[1, 1][2, 5] = 0
    hessF_chk_eval[1, 1][3, 1] = 0
    hessF_chk_eval[1, 1][3, 2] = 0
    hessF_chk_eval[1, 1][3, 3] = 0
    hessF_chk_eval[1, 1][3, 4] = 0
    hessF_chk_eval[1, 1][3, 5] = 0
    hessF_chk_eval[1, 1][4, 1] = 0
    hessF_chk_eval[1, 1][4, 2] = 0
    hessF_chk_eval[1, 1][4, 3] = 0
    hessF_chk_eval[1, 1][4, 4] = 0
    hessF_chk_eval[1, 1][4, 5] = 0
    hessF_chk_eval[1, 1][5, 1] = 0
    hessF_chk_eval[1, 1][5, 2] = 0
    hessF_chk_eval[1, 1][5, 3] = 0
    hessF_chk_eval[1, 1][5, 4] = 0
    hessF_chk_eval[1, 1][5, 5] = 0
    hessF_chk_eval[1, 2] = Matrix{Float64}(undef, 5, 5)
    hessF_chk_eval[1, 2][1, 1] = 0
    hessF_chk_eval[1, 2][1, 2] = 0
    hessF_chk_eval[1, 2][1, 3] = 0
    hessF_chk_eval[1, 2][1, 4] = 0
    hessF_chk_eval[1, 2][1, 5] = 0
    hessF_chk_eval[1, 2][2, 1] = 0
    hessF_chk_eval[1, 2][2, 2] = 0
    hessF_chk_eval[1, 2][2, 3] = 0
    hessF_chk_eval[1, 2][2, 4] = 0
    hessF_chk_eval[1, 2][2, 5] = 0
    hessF_chk_eval[1, 2][3, 1] = 0
    hessF_chk_eval[1, 2][3, 2] = 0
    hessF_chk_eval[1, 2][3, 3] = 0
    hessF_chk_eval[1, 2][3, 4] = 0
    hessF_chk_eval[1, 2][3, 5] = 0
    hessF_chk_eval[1, 2][4, 1] = 0
    hessF_chk_eval[1, 2][4, 2] = 0
    hessF_chk_eval[1, 2][4, 3] = 0
    hessF_chk_eval[1, 2][4, 4] = 0
    hessF_chk_eval[1, 2][4, 5] = 0
    hessF_chk_eval[1, 2][5, 1] = 0
    hessF_chk_eval[1, 2][5, 2] = 0
    hessF_chk_eval[1, 2][5, 3] = 0
    hessF_chk_eval[1, 2][5, 4] = 0
    hessF_chk_eval[1, 2][5, 5] = 0
    hessF_chk_eval[2, 1] = Matrix{Float64}(undef, 5, 5)
    hessF_chk_eval[2, 1][1, 1] = 0
    hessF_chk_eval[2, 1][1, 2] = 0
    hessF_chk_eval[2, 1][1, 3] = 0
    hessF_chk_eval[2, 1][1, 4] = 0
    hessF_chk_eval[2, 1][1, 5] = 0
    hessF_chk_eval[2, 1][2, 1] = 0
    hessF_chk_eval[2, 1][2, 2] = 0
    hessF_chk_eval[2, 1][2, 3] = 0
    hessF_chk_eval[2, 1][2, 4] = 0
    hessF_chk_eval[2, 1][2, 5] = 0
    hessF_chk_eval[2, 1][3, 1] = 0
    hessF_chk_eval[2, 1][3, 2] = 0
    hessF_chk_eval[2, 1][3, 3] = 0
    hessF_chk_eval[2, 1][3, 4] = 0
    hessF_chk_eval[2, 1][3, 5] = 0
    hessF_chk_eval[2, 1][4, 1] = 0
    hessF_chk_eval[2, 1][4, 2] = 0
    hessF_chk_eval[2, 1][4, 3] = 0
    hessF_chk_eval[2, 1][4, 4] = 0
    hessF_chk_eval[2, 1][4, 5] = 0
    hessF_chk_eval[2, 1][5, 1] = 0
    hessF_chk_eval[2, 1][5, 2] = 0
    hessF_chk_eval[2, 1][5, 3] = 0
    hessF_chk_eval[2, 1][5, 4] = 0
    hessF_chk_eval[2, 1][5, 5] = 0
    hessF_chk_eval[2, 2] = Matrix{Float64}(undef, 5, 5)
    hessF_chk_eval[2, 2][1, 1] = 0
    hessF_chk_eval[2, 2][1, 2] = 0
    hessF_chk_eval[2, 2][1, 3] = 0
    hessF_chk_eval[2, 2][1, 4] = 0
    hessF_chk_eval[2, 2][1, 5] = 0
    hessF_chk_eval[2, 2][2, 1] = 0
    hessF_chk_eval[2, 2][2, 2] = 0
    hessF_chk_eval[2, 2][2, 3] = 0
    hessF_chk_eval[2, 2][2, 4] = 0
    hessF_chk_eval[2, 2][2, 5] = 0
    hessF_chk_eval[2, 2][3, 1] = 0
    hessF_chk_eval[2, 2][3, 2] = 0
    hessF_chk_eval[2, 2][3, 3] = 0
    hessF_chk_eval[2, 2][3, 4] = 0
    hessF_chk_eval[2, 2][3, 5] = 0
    hessF_chk_eval[2, 2][4, 1] = 0
    hessF_chk_eval[2, 2][4, 2] = 0
    hessF_chk_eval[2, 2][4, 3] = 0
    hessF_chk_eval[2, 2][4, 4] = 0
    hessF_chk_eval[2, 2][4, 5] = 0
    hessF_chk_eval[2, 2][5, 1] = 0
    hessF_chk_eval[2, 2][5, 2] = 0
    hessF_chk_eval[2, 2][5, 3] = 0
    hessF_chk_eval[2, 2][5, 4] = 0
    hessF_chk_eval[2, 2][5, 5] = 0
    hessF_chk_eval[3, 1] = Matrix{Float64}(undef, 5, 5)
    hessF_chk_eval[3, 1][1, 1] = 0
    hessF_chk_eval[3, 1][1, 2] = 0
    hessF_chk_eval[3, 1][1, 3] = 0
    hessF_chk_eval[3, 1][1, 4] = 0
    hessF_chk_eval[3, 1][1, 5] = 0
    hessF_chk_eval[3, 1][2, 1] = 0
    hessF_chk_eval[3, 1][2, 2] = 0
    hessF_chk_eval[3, 1][2, 3] = 0
    hessF_chk_eval[3, 1][2, 4] = 0
    hessF_chk_eval[3, 1][2, 5] = 0
    hessF_chk_eval[3, 1][3, 1] = 0
    hessF_chk_eval[3, 1][3, 2] = 0
    hessF_chk_eval[3, 1][3, 3] = 0
    hessF_chk_eval[3, 1][3, 4] = 0
    hessF_chk_eval[3, 1][3, 5] = 0
    hessF_chk_eval[3, 1][4, 1] = 0
    hessF_chk_eval[3, 1][4, 2] = 0
    hessF_chk_eval[3, 1][4, 3] = 0
    hessF_chk_eval[3, 1][4, 4] = 0
    hessF_chk_eval[3, 1][4, 5] = 0
    hessF_chk_eval[3, 1][5, 1] = 0
    hessF_chk_eval[3, 1][5, 2] = 0
    hessF_chk_eval[3, 1][5, 3] = 0
    hessF_chk_eval[3, 1][5, 4] = 0
    hessF_chk_eval[3, 1][5, 5] = 326592
    hessF_chk_eval[3, 2] = Matrix{Float64}(undef, 5, 5)
    hessF_chk_eval[3, 2][1, 1] = -(241*sin(4/3))/70
    hessF_chk_eval[3, 2][1, 2] = 0
    hessF_chk_eval[3, 2][1, 3] = 0
    hessF_chk_eval[3, 2][1, 4] = 0
    hessF_chk_eval[3, 2][1, 5] = 0
    hessF_chk_eval[3, 2][2, 1] = 0
    hessF_chk_eval[3, 2][2, 2] = 0
    hessF_chk_eval[3, 2][2, 3] = 0
    hessF_chk_eval[3, 2][2, 4] = 0
    hessF_chk_eval[3, 2][2, 5] = 0
    hessF_chk_eval[3, 2][3, 1] = 0
    hessF_chk_eval[3, 2][3, 2] = 0
    hessF_chk_eval[3, 2][3, 3] = 0
    hessF_chk_eval[3, 2][3, 4] = 0
    hessF_chk_eval[3, 2][3, 5] = 0
    hessF_chk_eval[3, 2][4, 1] = 0
    hessF_chk_eval[3, 2][4, 2] = 0
    hessF_chk_eval[3, 2][4, 3] = 0
    hessF_chk_eval[3, 2][4, 4] = 0
    hessF_chk_eval[3, 2][4, 5] = 0
    hessF_chk_eval[3, 2][5, 1] = 0
    hessF_chk_eval[3, 2][5, 2] = 0
    hessF_chk_eval[3, 2][5, 3] = 0
    hessF_chk_eval[3, 2][5, 4] = 0
    hessF_chk_eval[3, 2][5, 5] = 0



    # fdp:
    fdp_chk_eval = Vector{Float64}(undef, 2)
    fdp_chk_eval[1] = (12167*exp(-16342/195))/3500 + 3206852775/2970344
    fdp_chk_eval[2] = (1587*exp(-16342/195))/4900 + 3206852775/2970344

    # cedp:
    cedp_chk_eval = Vector{Vector{Float64}}(undef, 2)
    cedp_chk_eval[1] = Vector{Float64}(undef, 2)
    cedp_chk_eval[1][1] = (6046681*cos(1/7)*cos(23/10))/38025
    cedp_chk_eval[1][2] = 1771561/27000
    cedp_chk_eval[2] = Vector{Float64}(undef, 2)
    cedp_chk_eval[2][1] = -(6046681*sin(1/7)*sin(23/10))/38025
    cedp_chk_eval[2][2] = 49729/676

    # cidp:
    cidp_chk_eval = Vector{Vector{Float64}}(undef, 2)
    cidp_chk_eval[1] = Vector{Float64}(undef, 3)
    cidp_chk_eval[1][1] = 4617/350
    cidp_chk_eval[1][2] = -151/91
    cidp_chk_eval[1][3] = 550999/1092
    cidp_chk_eval[2] = Vector{Float64}(undef, 3)
    cidp_chk_eval[2][1] = 4617/350
    cidp_chk_eval[2][2] = 151/91
    cidp_chk_eval[2][3] = -550999/1092

    # Fdp:
    Fdp_chk_eval = Vector{Matrix{Float64}}(undef, 2)
    Fdp_chk_eval[1] = Matrix{Float64}(undef, 3, 2)
    Fdp_chk_eval[1][1, 1] = 0
    Fdp_chk_eval[1][1, 2] = -33673/910
    Fdp_chk_eval[1][2, 1] = 93/26
    Fdp_chk_eval[1][2, 2] = 0
    Fdp_chk_eval[1][3, 1] = 0
    Fdp_chk_eval[1][3, 2] = sin(4/3)
    Fdp_chk_eval[2] = Matrix{Float64}(undef, 3, 2)
    Fdp_chk_eval[2][1, 1] = 0
    Fdp_chk_eval[2][1, 2] = 33673/910
    Fdp_chk_eval[2][2, 1] = -93/26
    Fdp_chk_eval[2][2, 2] = 0
    Fdp_chk_eval[2][3, 1] = 0
    Fdp_chk_eval[2][3, 2] = sin(4/3)



    # gradfdp:
    gradfdp_chk_eval = Vector{Vector{Float64}}(undef, 2)
    gradfdp_chk_eval[1] = Vector{Float64}(undef, 5)
    gradfdp_chk_eval[1][1] = 9620558325/5940688 - (12167*exp(-16342/195))/1750
    gradfdp_chk_eval[1][2] = 1781584875/1485172 - (36501*exp(-16342/195))/3500
    gradfdp_chk_eval[1][3] = 458121825/371293 - (12167*exp(-16342/195))/700
    gradfdp_chk_eval[1][4] = 16034263875/228488 - (12167*exp(-16342/195))/500
    gradfdp_chk_eval[1][5] = 1924111665/1485172 - (133837*exp(-16342/195))/3500
    gradfdp_chk_eval[2] = Vector{Float64}(undef, 5)
    gradfdp_chk_eval[2][1] = 9620558325/5940688 - (1587*exp(-16342/195))/2450
    gradfdp_chk_eval[2][2] = 1781584875/1485172 - (4761*exp(-16342/195))/4900
    gradfdp_chk_eval[2][3] = 458121825/371293 - (1587*exp(-16342/195))/980
    gradfdp_chk_eval[2][4] = 16034263875/228488 - (1587*exp(-16342/195))/700
    gradfdp_chk_eval[2][5] = 1924111665/1485172 - (17457*exp(-16342/195))/4900

    # gradcedp:
    gradcedp_chk_eval = Vector{Matrix{Float64}}(undef, 2)
    gradcedp_chk_eval[1] = Matrix{Float64}(undef, 5, 2)
    gradcedp_chk_eval[1][1, 1] = (4918*cos(1/7)*cos(23/10))/195
    gradcedp_chk_eval[1][1, 2] = 14641/300
    gradcedp_chk_eval[1][2, 1] = (4918*cos(1/7)*cos(23/10))/195
    gradcedp_chk_eval[1][2, 2] = 14641/300
    gradcedp_chk_eval[1][3, 1] = (4918*cos(1/7)*cos(23/10))/195
    gradcedp_chk_eval[1][3, 2] = 0
    gradcedp_chk_eval[1][4, 1] = (4918*cos(1/7)*cos(23/10))/195
    gradcedp_chk_eval[1][4, 2] = 0
    gradcedp_chk_eval[1][5, 1] = (4918*cos(1/7)*cos(23/10))/195
    gradcedp_chk_eval[1][5, 2] = 0
    gradcedp_chk_eval[2] = Matrix{Float64}(undef, 5, 2)
    gradcedp_chk_eval[2][1, 1] = -(4918*sin(1/7)*sin(23/10))/195
    gradcedp_chk_eval[2][1, 2] = 0
    gradcedp_chk_eval[2][2, 1] = -(4918*sin(1/7)*sin(23/10))/195
    gradcedp_chk_eval[2][2, 2] = 0
    gradcedp_chk_eval[2][3, 1] = -(4918*sin(1/7)*sin(23/10))/195
    gradcedp_chk_eval[2][3, 2] = 223/13
    gradcedp_chk_eval[2][4, 1] = -(4918*sin(1/7)*sin(23/10))/195
    gradcedp_chk_eval[2][4, 2] = 223/13
    gradcedp_chk_eval[2][5, 1] = -(4918*sin(1/7)*sin(23/10))/195
    gradcedp_chk_eval[2][5, 2] = 223/13

    # gradcidp:
    gradcidp_chk_eval = Vector{Matrix{Float64}}(undef, 2)
    gradcidp_chk_eval[1] = Matrix{Float64}(undef, 5, 3)
    gradcidp_chk_eval[1][1, 1] = 0
    gradcidp_chk_eval[1][1, 2] = 0
    gradcidp_chk_eval[1][1, 3] = -67195/182
    gradcidp_chk_eval[1][2, 1] = 171/35
    gradcidp_chk_eval[1][2, 2] = 0
    gradcidp_chk_eval[1][2, 3] = 67195/182
    gradcidp_chk_eval[1][3, 1] = 0
    gradcidp_chk_eval[1][3, 2] = 0
    gradcidp_chk_eval[1][3, 3] = 6191/42
    gradcidp_chk_eval[1][4, 1] = 0
    gradcidp_chk_eval[1][4, 2] = -151/7
    gradcidp_chk_eval[1][4, 3] = -6191/42
    gradcidp_chk_eval[1][5, 1] = 0
    gradcidp_chk_eval[1][5, 2] = -151/455
    gradcidp_chk_eval[1][5, 3] = 550999/2730
    gradcidp_chk_eval[2] = Matrix{Float64}(undef, 5, 3)
    gradcidp_chk_eval[2][1, 1] = 0
    gradcidp_chk_eval[2][1, 2] = 0
    gradcidp_chk_eval[2][1, 3] = 67195/182
    gradcidp_chk_eval[2][2, 1] = 171/35
    gradcidp_chk_eval[2][2, 2] = 0
    gradcidp_chk_eval[2][2, 3] = -67195/182
    gradcidp_chk_eval[2][3, 1] = 0
    gradcidp_chk_eval[2][3, 2] = 0
    gradcidp_chk_eval[2][3, 3] = -6191/42
    gradcidp_chk_eval[2][4, 1] = 0
    gradcidp_chk_eval[2][4, 2] = 151/7
    gradcidp_chk_eval[2][4, 3] = 6191/42
    gradcidp_chk_eval[2][5, 1] = 0
    gradcidp_chk_eval[2][5, 2] = 151/455
    gradcidp_chk_eval[2][5, 3] = -550999/2730

    # gradFdp:
    gradFdp_chk_eval = Vector{Matrix{Vector{Float64}}}(undef, 2)
    gradFdp_chk_eval[1] = Matrix{Vector{Float64}}(undef, 3, 2)
    gradFdp_chk_eval[1][1, 1] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[1][1, 1][1] = 0
    gradFdp_chk_eval[1][1, 1][2] = 0
    gradFdp_chk_eval[1][1, 1][3] = 0
    gradFdp_chk_eval[1][1, 1][4] = 0
    gradFdp_chk_eval[1][1, 1][5] = 0
    gradFdp_chk_eval[1][1, 2] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[1][1, 2][1] = 0
    gradFdp_chk_eval[1][1, 2][2] = 0
    gradFdp_chk_eval[1][1, 2][3] = -151/35
    gradFdp_chk_eval[1][1, 2][4] = -151/35
    gradFdp_chk_eval[1][1, 2][5] = -151/35
    gradFdp_chk_eval[1][2, 1] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[1][2, 1][1] = 0
    gradFdp_chk_eval[1][2, 1][2] = 0
    gradFdp_chk_eval[1][2, 1][3] = 1
    gradFdp_chk_eval[1][2, 1][4] = 1
    gradFdp_chk_eval[1][2, 1][5] = 0
    gradFdp_chk_eval[1][2, 2] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[1][2, 2][1] = 0
    gradFdp_chk_eval[1][2, 2][2] = 0
    gradFdp_chk_eval[1][2, 2][3] = 0
    gradFdp_chk_eval[1][2, 2][4] = 0
    gradFdp_chk_eval[1][2, 2][5] = 0
    gradFdp_chk_eval[1][3, 1] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[1][3, 1][1] = 0
    gradFdp_chk_eval[1][3, 1][2] = 0
    gradFdp_chk_eval[1][3, 1][3] = 0
    gradFdp_chk_eval[1][3, 1][4] = 0
    gradFdp_chk_eval[1][3, 1][5] = 0
    gradFdp_chk_eval[1][3, 2] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[1][3, 2][1] = cos(4/3)
    gradFdp_chk_eval[1][3, 2][2] = 0
    gradFdp_chk_eval[1][3, 2][3] = 0
    gradFdp_chk_eval[1][3, 2][4] = 0
    gradFdp_chk_eval[1][3, 2][5] = 0
    gradFdp_chk_eval[2] = Matrix{Vector{Float64}}(undef, 3, 2)
    gradFdp_chk_eval[2][1, 1] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[2][1, 1][1] = 0
    gradFdp_chk_eval[2][1, 1][2] = 0
    gradFdp_chk_eval[2][1, 1][3] = 0
    gradFdp_chk_eval[2][1, 1][4] = 0
    gradFdp_chk_eval[2][1, 1][5] = 0
    gradFdp_chk_eval[2][1, 2] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[2][1, 2][1] = 0
    gradFdp_chk_eval[2][1, 2][2] = 0
    gradFdp_chk_eval[2][1, 2][3] = 151/35
    gradFdp_chk_eval[2][1, 2][4] = 151/35
    gradFdp_chk_eval[2][1, 2][5] = 151/35
    gradFdp_chk_eval[2][2, 1] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[2][2, 1][1] = 0
    gradFdp_chk_eval[2][2, 1][2] = 0
    gradFdp_chk_eval[2][2, 1][3] = -1
    gradFdp_chk_eval[2][2, 1][4] = -1
    gradFdp_chk_eval[2][2, 1][5] = 0
    gradFdp_chk_eval[2][2, 2] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[2][2, 2][1] = 0
    gradFdp_chk_eval[2][2, 2][2] = 0
    gradFdp_chk_eval[2][2, 2][3] = 0
    gradFdp_chk_eval[2][2, 2][4] = 0
    gradFdp_chk_eval[2][2, 2][5] = 0
    gradFdp_chk_eval[2][3, 1] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[2][3, 1][1] = 0
    gradFdp_chk_eval[2][3, 1][2] = 0
    gradFdp_chk_eval[2][3, 1][3] = 0
    gradFdp_chk_eval[2][3, 1][4] = 0
    gradFdp_chk_eval[2][3, 1][5] = 0
    gradFdp_chk_eval[2][3, 2] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[2][3, 2][1] = cos(4/3)
    gradFdp_chk_eval[2][3, 2][2] = 0
    gradFdp_chk_eval[2][3, 2][3] = 0
    gradFdp_chk_eval[2][3, 2][4] = 0
    gradFdp_chk_eval[2][3, 2][5] = 0

    eval_test_pt = MPCCPointEval{Float64}(
        f_chk_eval,
        ce_chk_eval,
        ci_chk_eval,
        F_chk_eval,
        gradf_chk_eval,
        gradce_chk_eval,
        gradci_chk_eval,
        gradF_chk_eval,
        hessf_chk_eval,
        hessce_chk_eval,
        hessci_chk_eval,
        hessF_chk_eval,
        fdp_chk_eval,
        cedp_chk_eval,
        cidp_chk_eval,
        Fdp_chk_eval,
        gradfdp_chk_eval,
        gradcedp_chk_eval,
        gradcidp_chk_eval,
        gradFdp_chk_eval
    )

    return [ MPCCModelTestVector{Float64, Float64, Float64}(x_val, pr_val, ps_val, eval_test_pt) ]

end



function mm_spec_testdense_knownsols()
    return Vector{Function}(undef, 0)
end


function mm_spec_testdense_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t, 2t ]),
            (-1e9, 1e9),
            "Standard"    
        )

    return  Vector{MPCCParameterisationDefn}( [ defn1 ] )
end

