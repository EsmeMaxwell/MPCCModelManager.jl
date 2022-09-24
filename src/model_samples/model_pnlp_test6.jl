


function model_pnlp_test6_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(3, 0, 0, 1, 2, 1, 0)
    return dimspec
end



function model_pnlp_test6_defn(x, pr, ps)
    label_model = "Active Set Test 6 (NLP only); 3d sphere CE, and nearby flat regions as CI"

    f = exp(0.1*x[3])
    label_f = "Exp in z"

    ce = Vector{Num}(undef, 1)
    ce[1] = x[1]^2 + x[2]^2 + x[3]^3 - 1
    
    label_ce = Vector{String}([])

    ci = Vector{Num}(undef, 2)
    ci[1] = x[3] -exp(-x[1]) - x[2]
    ci[2] = x[3] + pr[1]

    label_ci = Vector{String}([
            "Mildly curved 'base' (active)",
            "Plane (inactive)"
        ])

    F = Matrix{Num}(undef, 0, 0)
    label_F = Vector{String}([])

    return MPCCDefinition(f, ce, ci, F, label_model, label_f, label_ce, label_ci, label_F)
end




function model_pnlp_test6_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function model_pnlp_test6_knownsols()
    return Vector{Function}([ ])
end



function model_pnlp_test6_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t ]),
            (0.1, 1e-6),
            "Standard"    
        )

    return Vector{MPCCParameterisationDefn}( [ defn1 ] )
end






function model_pnlp_test6_build()
    dimspec = model_pnlp_test6_dimspec()

    (x, pr, ps, t) = mpccmodel_build_sym_nums(dimspec)
    defn = model_pnlp_test6_defn(x, pr, ps)

    pdefns = model_pnlp_test6_parameterisations(t)
    testvectors = model_pnlp_test6_testvectors()
    knownsols = model_pnlp_test6_knownsols()

    model_cfg = mpccmodel_construct_config(dimspec, defn; pdefns = pdefns, testvectors = testvectors, knownsols = knownsols)

    return model_cfg
end

