

# c1 = -x + 2.0 - e1;
# c2 = -x + 2.0 + e1;

function model_pnlp_test3_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(2, 0, 0, 0, 2, 1, 0)
    return dimspec
end



function model_pnlp_test3_defn(x, pr, ps)
    label_model = "Active Set Test 3 (NLP only); Rosenbrock"

    f = (1.0-x[1])^2 + 100.0*(x[2]-x[1]^2)^2
    label_f = "Rosenbrock objective"

    ce = Vector{Num}()
    label_ce = Vector{String}([])

    ci = Vector{Num}(undef, 2)
    ci[1] = x[2] + x[1] - 2.0 + pr[1]
    ci[2] = x[2] + x[1] - 2.0 - pr[1]

    label_ci = Vector{String}([
            "Straight line (active)",
            "Straight line (inactive)"
        ])

    F = Matrix{Num}(undef, 0, 0)
    label_F = Vector{String}([])

    return MPCCDefinition(f, ce, ci, F, label_model, label_f, label_ce, label_ci, label_F)
end




function model_pnlp_test3_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function model_pnlp_test3_knownsols()

    return Vector{Function}([ ])
end



function model_pnlp_test3_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t ]),
            (0.05, 0.0),
            "Standard"    
        )

    return  Vector{MPCCParameterisationDefn}( [ defn1 ] )
end






function model_pnlp_test3_build()
    dimspec = model_pnlp_test3_dimspec()

    (x, pr, ps, t) = mpccmodel_build_sym_nums(dimspec)
    defn = model_pnlp_test3_defn(x, pr, ps)

    pdefns = model_pnlp_test3_parameterisations(t)
    testvectors = model_pnlp_test3_testvectors()
    knownsols = model_pnlp_test3_knownsols()

    model_cfg = mpccmodel_construct_config(dimspec, defn; pdefns = pdefns, testvectors = testvectors, knownsols = knownsols)

    return model_cfg
end

