

function mm_spec_dev1_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(4, 2, 3, 2, 3, 2, 1)
    return dimspec
end


function mm_spec_dev1_defn(x, pr, ps)

    f = exp(-x[1]+x[2]+2x[3]-5x[4]) + x[1]^2*x[4]^3*pr[1]^2 + x[2]*pr[2]^2 + x[3]*pr[1]*pr[2]*ps[1]

    ce = Vector{Num}(undef, 2)
    ce[1] = (x[1]-x[2])^2 - pr[1]*x[2]*x[3]^3
    ce[2] = (x[4]-x[1])^3 - pr[1]^2*x[4]

    ci = Vector{Num}(undef, 3)
    ci[1] = sin(0.5*pi*x[1])*pr[2]^3 + cos(x[1]^2*x[2]^3*x[3]^3*x[4]^4)
    ci[2] = pr[1]*(x[3]-2x[2])^4
    ci[3] = -x[4]^3*pr[2]^2 + x[1]^2

    F = Matrix{Num}(undef, 3, 2)
    F[1,1] = x[1]^3*x[2]^3*pr[1]^2
    F[2,1] = x[2]
    F[3,1] = x[3]
    F[1,2] = (x[1]+x[4])^5
    F[2,2] = pr[1]^2*x[3]^2
    F[3,2] = pr[1]*pr[2]*(x[1]-x[2])^2

    return MPCCDefinition(f, ce, ci, F)
end



function mm_spec_dev1_nzmask()::MPCCModelNZMask
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


function mm_spec_dev1_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function mm_spec_dev1_knownsols()

    return [ ]
end



function mm_spec_dev1_parameterisations(t)

    return Vector{MPCCParameterisationDefn}( [ ] )
end
