

function mm_spec_kj3_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(2, 1, 2, 0, 0, 1, 0)
    return dimspec
end


function mm_spec_kj3_defn(x, pr, ps)

    f = (x[1] - pr[1])^2 + (x[2] - pr[1])^2

    ce = Vector{Num}()

    ci = Vector{Num}()

    F = Matrix{Num}(undef, 2, 1)
    F[1,1] = x[1]
    F[2,1] = x[2]

    return MPCCDefinition(f, ce, ci, F)
end



function mm_spec_kj3_nzmask()::MPCCModelNZMask
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


function mm_spec_kj3_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function mm_spec_kj3_knownsols()

    function local_kj3_knownsol_1(pr::Vector{R}, ps::Vector{Int64}) where {R <: Real}
        @assert ( 1 == length(pr) && 0 == length(ps) )
        t = pr[1]

        if ( t <= 0 )
            x_sol = [0, 0]
            f_sol = 2(-t)^2
        else
            x_sol = [t, 0]
            f_sol = (x_sol[1] - t)^2 + (x_sol[2] - t)^2
        end

        return (x_sol, f_sol)
    end

    function local_kj3_knownsol_2(pr::Vector{R}, ps::Vector{Int64}) where {R <: Real}
        @assert ( 1 == length(pr) && 0 == length(ps) )
        t = pr[1]

        if ( t <= 0 )
            x_sol = [0, 0]
            f_sol = 2(-t)^2
        else
            x_sol = [0, t]
            f_sol = (x_sol[1] - t)^2 + (x_sol[2] - t)^2
        end

        return (x_sol, f_sol)
    end

    return [ local_kj3_knownsol_1, local_kj3_knownsol_2 ]
end



function mm_spec_kj3_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t ]),
            (-5.0, 5.0),
            "Standard"    
        )

    return  Vector{MPCCParameterisationDefn}( [ defn1 ] )
end

