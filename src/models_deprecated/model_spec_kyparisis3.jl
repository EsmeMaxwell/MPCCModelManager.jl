

function mm_spec_kyparisis3_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(2, 0, 0, 0, 3, 2, 0)
    return dimspec
end


function mm_spec_kyparisis3_defn(x, pr, ps)

    f = (1/2)*(x[1]-pr[1])^2 + (1/2)*(x[2]-pr[2])^2 + x[1] + x[2]

    ce = Vector{Num}()

    ci = Vector{Num}(undef, 3)
    ci[1] = x[1]
    ci[2] = x[2]
    ci[3] = x[1] + x[2] - pr[1] - pr[2]

    F = Matrix{Num}(undef, 0, 0)
    # F[1,1] = x[1]
    # F[2,1] = x[2]

    return MPCCDefinition(f, ce, ci, F)
end



function mm_spec_kyparisis3_nzmask()::MPCCModelNZMask
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


function mm_spec_kyparisis3_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function mm_spec_kyparisis3_knownsols()

    function local_kyparisis3_knownsol(pr::Vector{R}, ps::Vector{Int64}) where {R <: Real}
        @assert ( 2 == length(pr) && 0 == length(ps) )
        e1 = pr[1]
        e2 = pr[2]
        
        if ( e1 >= 0 && e2 >= 0 )
            x_sol = (e1, e2)
        elseif ( e1 >= 0 && e2 < 0 )
            x_sol = (max(e1-1, e1 + e2, 0), 0)
        elseif ( e1 < 0 && e2 >= 0 )
            x_sol = (0, max(e2-1, e1 + e2, 0))
        else
            x_sol = (0, 0)
        end
        f_sol = (1/2)*(x_sol[1]-e1)^2 + (1/2)*(x_sol[2]-e2)^2 + x_sol[1] + x_sol[2]

        return (x_sol, f_sol)
    end

    return [ local_kyparisis3_knownsol ]
end




function mm_spec_kyparisis3_parameterisations(t)

    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ 3.0, t ]),
            (-3.0, 3.0),
            "Vertical, rising"
        )

    return Vector{MPCCParameterisationDefn}( [ defn1 ] )
end

