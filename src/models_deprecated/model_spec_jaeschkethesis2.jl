

function mm_spec_jaeschkethesis2_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(5, 0, 0, 3, 7, 9, 0)
    # x = [ c_A c_B c_C F_A F_B]
    # p = [ k_1 k_2 (-\Delta H_1) (-\Delta H_2) c_{A,in} c_{B,in} V F_{max}, q_{max} ]
    return dimspec
end


function mm_spec_jaeschkethesis2_defn(x, pr, ps)

    f = ( (-1) * (x[4] + x[5])^2 * (x[3]^2) ) / (x[4] * pr[5])

    ce = Vector{Num}(undef, 3)
    ce[1] = (x[4] * pr[5]) - ((x[4] + x[5]) * x[1]) - (pr[1] * x[1] * x[2] * pr[7])
    ce[2] = (x[5] * pr[6]) - ((x[4] + x[5]) * x[2]) - (pr[1] * x[1] * x[2] * pr[7]) - (2 * pr[2] * x[2]^2 * pr[7])
    ce[3] = -((x[4] + x[5]) * x[3]) + (pr[1] * x[1] * x[2] * pr[7])

    ci = Vector{Num}(undef, 7)
    ci[1] = x[1] + 1e-15
    ci[2] = x[2] + 1e-15
    ci[3] = x[3] + 1e-15
    ci[4] = x[4] + 1e-15
    ci[5] = x[5] + 1e-15
    ci[6] = pr[8] - (x[4] + x[5])
    ci[7] = pr[9] - (pr[1] * x[1] * x[2] * pr[7] * pr[3]) - (2 * pr[2] * x[2]^2 * pr[7] * pr[4])


    F = Matrix{Num}(undef, 0, 0)
    # F[1,1] = x[1]
    # F[2,1] = x[2]

    return MPCCDefinition(f, ce, ci, F)
end



function mm_spec_jaeschkethesis2_nzmask()::MPCCModelNZMask
    nz_mask = MPCCModelNZMask(
        # ce
        [],
        # ci
        [   Set{Int64}([1, 2]),
            Set{Int64}([1]) ],
        # F
        [ 
            # F1
            [   Set{Int64}([1]) ],
            # F2
            [   Set{Int64}([2]) ] 
        ] )    
    return nz_mask
end


function mm_spec_jaeschkethesis2_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function mm_spec_jaeschkethesis2_knownsols()

    return [ ]
end




function mm_spec_jaeschkethesis2_parameterisations(t)

    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t, 0.014, 7e4, 5e4, 2, 1.5, 500, 22, 1e6 ]),
            (0.3, 1.5),
            "Standard"
        )

    return Vector{MPCCParameterisationDefn}( [ defn1 ] )
end

