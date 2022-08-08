
# TODO rename testsparse



function mm_spec_testsparse_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(5, 2, 3, 2, 3, 2, 0)
    return dimspec
end



function mm_spec_testsparse_defn(x, pr, ps)

    f= x[1]^2 * x[2]^3 * x[3]^4 * ( pr[1] )^2 + pr[2]^3 * exp( - 7*x[4] - 11*x[5])

    ce = Vector{Num}(undef, 2)
    ce[1] = [ x[1]^2 * sin(pr[1]^2) ]
    ce[2] = [ pr[1] * ( x[3] + x[5] ) + pr[2] * x[4]^2 ]
    
    ci = Vector{Num}(undef, 3)
    ci[1] = [ x[1] + x[2] * ( pr[1] + pr[2] ) ]
    ci[2] = [ x[3] + x[4] * ( pr[1] - pr[2] )^2 ] 
    ci[3] = [ ( x[2] - x[1] ) * ( x[4] - x[3] ) * x[5]^2 * ( pr[2] - pr[1] ) ]

    F = Matrix{Num}(undef, 3, 2)
    F[1,1] = [ x[1] ]
    F[1,2] = [ x[2] ]
    F[2,1] = [ x[3] * pr[2] ]
    F[2,2] = [ -x[1] ]
    F[3,1] = [ ( 1 + x[5] )^7 ]
    F[3,2] = [ ( 1 + pr[1] + pr[2] ) * sin(x[1]) ]

    return MPCCDefinition(f, ce, ci, F)
end


function mm_spec_testsparse_nzmask()::MPCCModelNZMask
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



function mm_spec_testsparse_testvectors()

    x_val = Vector{Rational}(undef, 5)
    x_val = [4//3, 27//10, 7//2, 1//13, 5]

    pr_val = Vector{Rational}(undef, 2)
    pr_val = [1//7, 23//10]

    ps_val = Int64[]

    # f:
    f_chk_eval = (12167*exp(-722/13))/1000 + 107163/1000

    # ce:
    ce_chk_eval = Vector{Float64}(undef, 2)
    ce_chk_eval[1] = (16*sin(1/49))/9
    ce_chk_eval[2] = 7263/5915

    # ci:
    ci_chk_eval = Vector{Float64}(undef, 3)
    ci_chk_eval[1] = 16651/2100
    ci_chk_eval[2] = 245751/63700
    ci_chk_eval[3] = -550999/2184

    # F:
    F_chk_eval = Matrix{Float64}(undef, 3, 2)
    F_chk_eval[1, 1] = 4/3
    F_chk_eval[1, 2] = 27/10
    F_chk_eval[2, 1] = 161/20
    F_chk_eval[2, 2] = -4/3
    F_chk_eval[3, 1] = 279936
    F_chk_eval[3, 2] = (241*sin(4/3))/70



    # gradf:
    gradf_chk_eval = Vector{Float64}(undef, 5)
    gradf_chk_eval[1] = 321489/2000
    gradf_chk_eval[2] = 11907/100
    gradf_chk_eval[3] = 15309/125
    gradf_chk_eval[4] = -(85169*exp(-722/13))/1000
    gradf_chk_eval[5] = -(133837*exp(-722/13))/1000

    # ce:
    gradce_chk_eval = Matrix{Float64}(undef, 5, 2)
    gradce_chk_eval[1, 1] = (8*sin(1/49))/3
    gradce_chk_eval[1, 2] = 0
    gradce_chk_eval[2, 1] = 0
    gradce_chk_eval[2, 2] = 0
    gradce_chk_eval[3, 1] = 0
    gradce_chk_eval[3, 2] = 1/7
    gradce_chk_eval[4, 1] = 0
    gradce_chk_eval[4, 2] = 23/65
    gradce_chk_eval[5, 1] = 0
    gradce_chk_eval[5, 2] = 1/7

    # gradci:
    gradci_chk_eval = Matrix{Float64}(undef, 5, 3)
    gradci_chk_eval[1, 1] = 1
    gradci_chk_eval[1, 2] = 0
    gradci_chk_eval[1, 3] = 67195/364
    gradci_chk_eval[2, 1] = 171/70
    gradci_chk_eval[2, 2] = 0
    gradci_chk_eval[2, 3] = -67195/364
    gradci_chk_eval[3, 1] = 0
    gradci_chk_eval[3, 2] = 1
    gradci_chk_eval[3, 3] = -6191/84
    gradci_chk_eval[4, 1] = 0
    gradci_chk_eval[4, 2] = 22801/4900
    gradci_chk_eval[4, 3] = 6191/84
    gradci_chk_eval[5, 1] = 0
    gradci_chk_eval[5, 2] = 0
    gradci_chk_eval[5, 3] = -550999/5460

    # gradF:
    gradF_chk_eval = Matrix{Vector{Float64}}(undef, 3, 2)
    gradF_chk_eval[1, 1] = Vector{Float64}(undef, 5)
    gradF_chk_eval[1, 1][1] = 1
    gradF_chk_eval[1, 1][2] = 0
    gradF_chk_eval[1, 1][3] = 0
    gradF_chk_eval[1, 1][4] = 0
    gradF_chk_eval[1, 1][5] = 0
    gradF_chk_eval[1, 2] = Vector{Float64}(undef, 5)
    gradF_chk_eval[1, 2][1] = 0
    gradF_chk_eval[1, 2][2] = 1
    gradF_chk_eval[1, 2][3] = 0
    gradF_chk_eval[1, 2][4] = 0
    gradF_chk_eval[1, 2][5] = 0
    gradF_chk_eval[2, 1] = Vector{Float64}(undef, 5)
    gradF_chk_eval[2, 1][1] = 0
    gradF_chk_eval[2, 1][2] = 0
    gradF_chk_eval[2, 1][3] = 23/10
    gradF_chk_eval[2, 1][4] = 0
    gradF_chk_eval[2, 1][5] = 0
    gradF_chk_eval[2, 2] = Vector{Float64}(undef, 5)
    gradF_chk_eval[2, 2][1] = -1
    gradF_chk_eval[2, 2][2] = 0
    gradF_chk_eval[2, 2][3] = 0
    gradF_chk_eval[2, 2][4] = 0
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
    hessf_chk_eval[1, 1] = 964467/8000
    hessf_chk_eval[1, 2] = 35721/200
    hessf_chk_eval[1, 3] = 45927/250
    hessf_chk_eval[1, 4] = 0
    hessf_chk_eval[1, 5] = 0
    hessf_chk_eval[2, 1] = 35721/200
    hessf_chk_eval[2, 2] = 441/5
    hessf_chk_eval[2, 3] = 3402/25
    hessf_chk_eval[2, 4] = 0
    hessf_chk_eval[2, 5] = 0
    hessf_chk_eval[3, 1] = 45927/250
    hessf_chk_eval[3, 2] = 3402/25
    hessf_chk_eval[3, 3] = 13122/125
    hessf_chk_eval[3, 4] = 0
    hessf_chk_eval[3, 5] = 0
    hessf_chk_eval[4, 1] = 0
    hessf_chk_eval[4, 2] = 0
    hessf_chk_eval[4, 3] = 0
    hessf_chk_eval[4, 4] = (596183*exp(-722/13))/1000
    hessf_chk_eval[4, 5] = (936859*exp(-722/13))/1000
    hessf_chk_eval[5, 1] = 0
    hessf_chk_eval[5, 2] = 0
    hessf_chk_eval[5, 3] = 0
    hessf_chk_eval[5, 4] = (936859*exp(-722/13))/1000
    hessf_chk_eval[5, 5] = (1472207*exp(-722/13))/1000

    # hessce:
    hessce_chk_eval = Vector{Matrix{Float64}}(undef, 2)
    hessce_chk_eval[1] = Matrix{Float64}(undef, 5, 5)
    hessce_chk_eval[1][1, 1] = 2*sin(1/49)
    hessce_chk_eval[1][1, 2] = 0
    hessce_chk_eval[1][1, 3] = 0
    hessce_chk_eval[1][1, 4] = 0
    hessce_chk_eval[1][1, 5] = 0
    hessce_chk_eval[1][2, 1] = 0
    hessce_chk_eval[1][2, 2] = 0
    hessce_chk_eval[1][2, 3] = 0
    hessce_chk_eval[1][2, 4] = 0
    hessce_chk_eval[1][2, 5] = 0
    hessce_chk_eval[1][3, 1] = 0
    hessce_chk_eval[1][3, 2] = 0
    hessce_chk_eval[1][3, 3] = 0
    hessce_chk_eval[1][3, 4] = 0
    hessce_chk_eval[1][3, 5] = 0
    hessce_chk_eval[1][4, 1] = 0
    hessce_chk_eval[1][4, 2] = 0
    hessce_chk_eval[1][4, 3] = 0
    hessce_chk_eval[1][4, 4] = 0
    hessce_chk_eval[1][4, 5] = 0
    hessce_chk_eval[1][5, 1] = 0
    hessce_chk_eval[1][5, 2] = 0
    hessce_chk_eval[1][5, 3] = 0
    hessce_chk_eval[1][5, 4] = 0
    hessce_chk_eval[1][5, 5] = 0
    hessce_chk_eval[2] = Matrix{Float64}(undef, 5, 5)
    hessce_chk_eval[2][1, 1] = 0
    hessce_chk_eval[2][1, 2] = 0
    hessce_chk_eval[2][1, 3] = 0
    hessce_chk_eval[2][1, 4] = 0
    hessce_chk_eval[2][1, 5] = 0
    hessce_chk_eval[2][2, 1] = 0
    hessce_chk_eval[2][2, 2] = 0
    hessce_chk_eval[2][2, 3] = 0
    hessce_chk_eval[2][2, 4] = 0
    hessce_chk_eval[2][2, 5] = 0
    hessce_chk_eval[2][3, 1] = 0
    hessce_chk_eval[2][3, 2] = 0
    hessce_chk_eval[2][3, 3] = 0
    hessce_chk_eval[2][3, 4] = 0
    hessce_chk_eval[2][3, 5] = 0
    hessce_chk_eval[2][4, 1] = 0
    hessce_chk_eval[2][4, 2] = 0
    hessce_chk_eval[2][4, 3] = 0
    hessce_chk_eval[2][4, 4] = 23/5
    hessce_chk_eval[2][4, 5] = 0
    hessce_chk_eval[2][5, 1] = 0
    hessce_chk_eval[2][5, 2] = 0
    hessce_chk_eval[2][5, 3] = 0
    hessce_chk_eval[2][5, 4] = 0
    hessce_chk_eval[2][5, 5] = 0

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
    hessci_chk_eval[2][4, 5] = 0
    hessci_chk_eval[2][5, 1] = 0
    hessci_chk_eval[2][5, 2] = 0
    hessci_chk_eval[2][5, 3] = 0
    hessci_chk_eval[2][5, 4] = 0
    hessci_chk_eval[2][5, 5] = 0
    hessci_chk_eval[3] = Matrix{Float64}(undef, 5, 5)
    hessci_chk_eval[3][1, 1] = 0
    hessci_chk_eval[3][1, 2] = 0
    hessci_chk_eval[3][1, 3] = 755/14
    hessci_chk_eval[3][1, 4] = -755/14
    hessci_chk_eval[3][1, 5] = 13439/182
    hessci_chk_eval[3][2, 1] = 0
    hessci_chk_eval[3][2, 2] = 0
    hessci_chk_eval[3][2, 3] = -755/14
    hessci_chk_eval[3][2, 4] = 755/14
    hessci_chk_eval[3][2, 5] = -13439/182
    hessci_chk_eval[3][3, 1] = 755/14
    hessci_chk_eval[3][3, 2] = -755/14
    hessci_chk_eval[3][3, 3] = 0
    hessci_chk_eval[3][3, 4] = 0
    hessci_chk_eval[3][3, 5] = -6191/210
    hessci_chk_eval[3][4, 1] = -755/14
    hessci_chk_eval[3][4, 2] = 755/14
    hessci_chk_eval[3][4, 3] = 0
    hessci_chk_eval[3][4, 4] = 0
    hessci_chk_eval[3][4, 5] = 6191/210
    hessci_chk_eval[3][5, 1] = 13439/182
    hessci_chk_eval[3][5, 2] = -13439/182
    hessci_chk_eval[3][5, 3] = -6191/210
    hessci_chk_eval[3][5, 4] = 6191/210
    hessci_chk_eval[3][5, 5] = -550999/27300

    # hessF:
    hessF_chk_eval = Matrix{Matrix{Float64}}(undef, 3, 2)
    hessF_chk_eval[1, 1] = Matrix{Float64}(undef, 5, 5)
    hessF_chk_eval[1, 1][1, 1] = 0
    hessF_chk_eval[1, 1][1, 2] = 0
    hessF_chk_eval[1, 1][1, 3] = 0
    hessF_chk_eval[1, 1][1, 4] = 0
    hessF_chk_eval[1, 1][1, 5] = 0
    hessF_chk_eval[1, 1][2, 1] = 0
    hessF_chk_eval[1, 1][2, 2] = 0
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
    fdp_chk_eval[1] = 750141/500
    fdp_chk_eval[2] = (1587*exp(-722/13))/100

    # cedp:
    cedp_chk_eval = Vector{Vector{Float64}}(undef, 2)
    cedp_chk_eval[1] = Vector{Float64}(undef, 2)
    cedp_chk_eval[1][1] = (32*cos(1/49))/63
    cedp_chk_eval[1][2] = 17/2
    cedp_chk_eval[2] = Vector{Float64}(undef, 2)
    cedp_chk_eval[2][1] = 0
    cedp_chk_eval[2][2] = 1/169

    # cidp:
    cidp_chk_eval = Vector{Vector{Float64}}(undef, 2)
    cidp_chk_eval[1] = Vector{Float64}(undef, 3)
    cidp_chk_eval[1][1] = 27/10
    cidp_chk_eval[1][2] = -151/455
    cidp_chk_eval[1][3] = 18245/156
    cidp_chk_eval[2] = Vector{Float64}(undef, 3)
    cidp_chk_eval[2][1] = 27/10
    cidp_chk_eval[2][2] = 151/455
    cidp_chk_eval[2][3] = -18245/156

    # Fdp:
    Fdp_chk_eval = Vector{Matrix{Float64}}(undef, 2)
    Fdp_chk_eval[1] = Matrix{Float64}(undef, 3, 2)
    Fdp_chk_eval[1][1, 1] = 0
    Fdp_chk_eval[1][1, 2] = 0
    Fdp_chk_eval[1][2, 1] = 0
    Fdp_chk_eval[1][2, 2] = 0
    Fdp_chk_eval[1][3, 1] = 0
    Fdp_chk_eval[1][3, 2] = sin(4/3)
    Fdp_chk_eval[2] = Matrix{Float64}(undef, 3, 2)
    Fdp_chk_eval[2][1, 1] = 0
    Fdp_chk_eval[2][1, 2] = 0
    Fdp_chk_eval[2][2, 1] = 7/2
    Fdp_chk_eval[2][2, 2] = 0
    Fdp_chk_eval[2][3, 1] = 0
    Fdp_chk_eval[2][3, 2] = sin(4/3)



    # gradfdp:
    gradfdp_chk_eval = Vector{Vector{Float64}}(undef, 2)
    gradfdp_chk_eval[1] = Vector{Float64}(undef, 5)
    gradfdp_chk_eval[1][1] = 2250423/1000
    gradfdp_chk_eval[1][2] = 83349/50
    gradfdp_chk_eval[1][3] = 214326/125
    gradfdp_chk_eval[1][4] = 0
    gradfdp_chk_eval[1][5] = 0
    gradfdp_chk_eval[2] = Vector{Float64}(undef, 5)
    gradfdp_chk_eval[2][1] = 0
    gradfdp_chk_eval[2][2] = 0
    gradfdp_chk_eval[2][3] = 0
    gradfdp_chk_eval[2][4] = -(11109*exp(-722/13))/100
    gradfdp_chk_eval[2][5] = -(17457*exp(-722/13))/100

    # gradcedp:
    gradcedp_chk_eval = Vector{Matrix{Float64}}(undef, 2)
    gradcedp_chk_eval[1] = Matrix{Float64}(undef, 5, 2)
    gradcedp_chk_eval[1][1, 1] = (16*cos(1/49))/21
    gradcedp_chk_eval[1][1, 2] = 0
    gradcedp_chk_eval[1][2, 1] = 0
    gradcedp_chk_eval[1][2, 2] = 0
    gradcedp_chk_eval[1][3, 1] = 0
    gradcedp_chk_eval[1][3, 2] = 1
    gradcedp_chk_eval[1][4, 1] = 0
    gradcedp_chk_eval[1][4, 2] = 0
    gradcedp_chk_eval[1][5, 1] = 0
    gradcedp_chk_eval[1][5, 2] = 1
    gradcedp_chk_eval[2] = Matrix{Float64}(undef, 5, 2)
    gradcedp_chk_eval[2][1, 1] = 0
    gradcedp_chk_eval[2][1, 2] = 0
    gradcedp_chk_eval[2][2, 1] = 0
    gradcedp_chk_eval[2][2, 2] = 0
    gradcedp_chk_eval[2][3, 1] = 0
    gradcedp_chk_eval[2][3, 2] = 0
    gradcedp_chk_eval[2][4, 1] = 0
    gradcedp_chk_eval[2][4, 2] = 2/13
    gradcedp_chk_eval[2][5, 1] = 0
    gradcedp_chk_eval[2][5, 2] = 0

    # gradcidp:
    gradcidp_chk_eval = Vector{Matrix{Float64}}(undef, 2)
    gradcidp_chk_eval[1] = Matrix{Float64}(undef, 5, 3)
    gradcidp_chk_eval[1][1, 1] = 0
    gradcidp_chk_eval[1][1, 2] = 0
    gradcidp_chk_eval[1][1, 3] = -2225/26
    gradcidp_chk_eval[1][2, 1] = 1
    gradcidp_chk_eval[1][2, 2] = 0
    gradcidp_chk_eval[1][2, 3] = 2225/26
    gradcidp_chk_eval[1][3, 1] = 0
    gradcidp_chk_eval[1][3, 2] = 0
    gradcidp_chk_eval[1][3, 3] = 205/6
    gradcidp_chk_eval[1][4, 1] = 0
    gradcidp_chk_eval[1][4, 2] = -151/35
    gradcidp_chk_eval[1][4, 3] = -205/6
    gradcidp_chk_eval[1][5, 1] = 0
    gradcidp_chk_eval[1][5, 2] = 0
    gradcidp_chk_eval[1][5, 3] = 3649/78
    gradcidp_chk_eval[2] = Matrix{Float64}(undef, 5, 3)
    gradcidp_chk_eval[2][1, 1] = 0
    gradcidp_chk_eval[2][1, 2] = 0
    gradcidp_chk_eval[2][1, 3] = 2225/26
    gradcidp_chk_eval[2][2, 1] = 1
    gradcidp_chk_eval[2][2, 2] = 0
    gradcidp_chk_eval[2][2, 3] = -2225/26
    gradcidp_chk_eval[2][3, 1] = 0
    gradcidp_chk_eval[2][3, 2] = 0
    gradcidp_chk_eval[2][3, 3] = -205/6
    gradcidp_chk_eval[2][4, 1] = 0
    gradcidp_chk_eval[2][4, 2] = 151/35
    gradcidp_chk_eval[2][4, 3] = 205/6
    gradcidp_chk_eval[2][5, 1] = 0
    gradcidp_chk_eval[2][5, 2] = 0
    gradcidp_chk_eval[2][5, 3] = -3649/78

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
    gradFdp_chk_eval[1][1, 2][3] = 0
    gradFdp_chk_eval[1][1, 2][4] = 0
    gradFdp_chk_eval[1][1, 2][5] = 0
    gradFdp_chk_eval[1][2, 1] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[1][2, 1][1] = 0
    gradFdp_chk_eval[1][2, 1][2] = 0
    gradFdp_chk_eval[1][2, 1][3] = 0
    gradFdp_chk_eval[1][2, 1][4] = 0
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
    gradFdp_chk_eval[2][1, 2][3] = 0
    gradFdp_chk_eval[2][1, 2][4] = 0
    gradFdp_chk_eval[2][1, 2][5] = 0
    gradFdp_chk_eval[2][2, 1] = Vector{Float64}(undef, 5)
    gradFdp_chk_eval[2][2, 1][1] = 0
    gradFdp_chk_eval[2][2, 1][2] = 0
    gradFdp_chk_eval[2][2, 1][3] = 1
    gradFdp_chk_eval[2][2, 1][4] = 0
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




function mm_spec_testsparse_knownsols()
    return Vector{Function}(undef, 0)
end






