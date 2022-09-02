

function model_pmpcc_empty_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(2, 0, 0, 0, 0, 1, 0)
    return dimspec
end


function model_pmpcc_empty_defn(x, pr, ps)
    label_model = "Empty PMPCC model for testing"

    f = (-x[1]+x[2])^2
    label_f = "Dummy objective fn"

    ce = Vector{Num}(undef, 0)
    ci = Vector{Num}(undef, 0)
    F = Matrix{Num}(undef, 0, 0)

    label_ce = Vector{String}([])
    label_ci = Vector{String}([])
    label_F = Vector{String}([])

    return MPCCDefinition(f, ce, ci, F, label_model, label_f, label_ce, label_ci, label_F)
end




function model_pmpcc_empty_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function model_pmpcc_empty_knownsols()
    return Vector{Function}(undef, 0)
end



function model_pmpcc_empty_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t ]),
            (-1.0, 2.0),
            "Standard"    
        )

    return  Vector{MPCCParameterisationDefn}( [ defn1 ] )
end





function model_pmpcc_empty_build()
    dimspec = model_pmpcc_empty_dimspec()

    (x, pr, ps, t) = mpccmodel_build_sym_nums(dimspec)
    defn = model_pmpcc_empty_defn(x, pr, ps)

    pdefns = model_pmpcc_empty_parameterisations(t)
    testvectors = model_pmpcc_empty_testvectors()
    knownsols = model_pmpcc_empty_knownsols()

    model_cfg = mpccmodel_construct_config(dimspec, defn; pdefns = pdefns, testvectors = testvectors, knownsols = knownsols)

    return model_cfg
end

