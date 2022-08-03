

function mm_spec_empty_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(2, 0, 0, 0, 0, 1, 0)
    return dimspec
end


function mm_spec_empty_defn(x, pr, ps)
    f = (-x[1]+x[2])^2

    ce = Vector{Num}(undef, 0)
    ci = Vector{Num}(undef, 0)
    F = Matrix{Num}(undef, 0, 0)

    return MPCCDefinition(f, ce, ci, F)
end




function mm_spec_empty_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function mm_spec_empty_knownsols()
    return Vector{Function}(undef, 0)
end



function mm_spec_empty_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t ]),
            (-1.0, 2.0),
            "Standard"    
        )

    return  Vector{MPCCParameterisationDefn}( [ defn1 ] )
end
