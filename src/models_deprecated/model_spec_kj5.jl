

function mm_spec_kj5_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(3, 1, 2, 0, 2, 1, 0)
    return dimspec
end



function mm_spec_kj5_defn(x, pr, ps)

    f = x[1] + x[2] - (1-pr[1]) * x[3]

    ce = Vector{Num}()

    ci = Vector{Num}(undef, 2)
    ci[1] = 4x[1] - x[3]
    ci[2] = 4x[2] - x[3]

    F = Matrix{Num}(undef, 2, 1)
    F[1,1] = x[1]
    F[2,1] = x[2]

    return MPCCDefinition(f, ce, ci, F)
end



function mm_spec_kj5_nzmask()::MPCCModelNZMask
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


function mm_spec_kj5_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function mm_spec_kj5_knownsols()

    function local_kj5_knownsol(pr::Vector{R}, ps::Vector{Int64}) where {R <: Real}
        @assert ( 1 == length(pr) && 0 == length(ps) )

        t = pr[1]
        if ( t < 0 || t > 1)
            error("t out of range, must be in [0,1]")
        end

        x_sol = [0, 0, 0]
        f_sol = 0

        return (x_sol, f_sol)
    end

    return [ local_kj5_knownsol ]
end

