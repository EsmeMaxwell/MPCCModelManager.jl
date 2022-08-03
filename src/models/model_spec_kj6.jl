

function mm_spec_kj6_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(2, 1, 2, 0, 2, 1, 0)
    return dimspec
end


function mm_spec_kj6_defn(x, pr, ps)
    f = exp(-x[1]+x[2])

    ce = Vector{Num}()

    ci = Vector{Num}(undef, 2)
    ci[1] = ((x[1]-2)^2 + (x[2]+1)^2 -2pr[1] -6)
    ci[2] = (1 - x[1])

    F = Matrix{Num}(undef, 2, 1)
    F[1,1] = x[1]
    F[2,1] = x[2]

    return MPCCDefinition(f, ce, ci, F)
end



# function mm_spec_kj6_nzmask()::MPCCModelNZMask
#     nz_mask = MPCCModelNZMask(
#         # ce
#         [],
#         # ci
#         [   Set{Int64}([1, 2]),
#             Set{Int64}([1]) ],
#         # F
#         [ 
#             # F1
#             [   Set{Int64}([1]) ],
#             # F2
#             [   Set{Int64}([2]) ] 
#         ] )    
#     return nz_mask
# end


function mm_spec_kj6_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function mm_spec_kj6_knownsols()

    function local_kj6_knownsol(pr::Vector{R}, ps::Vector{Int64}) where {R <: Real}
        @assert ( 1 == length(pr) && 0 == length(ps) )
        t = pr[1]
        if ( t <= -2 )
            x_sol = [1.0, 0.0]
            f_sol = exp(-1)
        elseif ( t <= -(1//2) && t > -2 )
            x_sol = [ (-sqrt(5 + 2t) + 2), 0]
            f_sol = exp(-x_sol[1])
        else
            x_sol = [ 0, (sqrt(6 + 2*t - (0-2)^2) - 1) ]
            f_sol = exp(x_sol[2])
        end

        return (x_sol, f_sol)
    end

    return [ local_kj6_knownsol ]
end



function mm_spec_kj6_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t ]),
            (-5.0, 2.0),
            "Standard"    
        )

    return  Vector{MPCCParameterisationDefn}( [ defn1 ] )
end
