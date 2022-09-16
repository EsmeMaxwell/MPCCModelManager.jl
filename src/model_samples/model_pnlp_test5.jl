

# c1 = 1 -exp(-x.^2 / 2*e1^2);
# c2 = e2*x + e3;

function model_pnlp_test5_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(2, 0, 0, 0, 2, 3, 0)
    return dimspec
end



function model_pnlp_test5_defn(x, pr, ps)
    label_model = "Active Set Test 5 (NLP only); contracting Gaussian with straight line"

    f = exp(0.01*(x[1]^2 + x[2]^2))
    label_f = "Exp of lin objective"

    ce = Vector{Num}()
    label_ce = Vector{String}([])

    ci = Vector{Num}(undef, 2)
    ci[1] = x[2] - 1 + exp(-x[1]^2 / 2*pr[1]^2)
    ci[2] = x[2] - pr[2]*x[1] - pr[3]

    label_ci = Vector{String}([
            "Gaussian (active)",
            "Straight line (active)"
        ])

    F = Matrix{Num}(undef, 0, 0)
    label_F = Vector{String}([])

    return MPCCDefinition(f, ce, ci, F, label_model, label_f, label_ce, label_ci, label_F)
end




function model_pnlp_test5_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function model_pnlp_test5_knownsols()
    return Vector{Function}([ ])
end



function model_pnlp_test5_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t, 0.1, 0.05 ]),
            (10.0, 1000.0),
            "Standard"    
        )

    return Vector{MPCCParameterisationDefn}( [ defn1 ] )
end






function model_pnlp_test5_build()
    dimspec = model_pnlp_test5_dimspec()

    (x, pr, ps, t) = mpccmodel_build_sym_nums(dimspec)
    defn = model_pnlp_test5_defn(x, pr, ps)

    pdefns = model_pnlp_test5_parameterisations(t)
    testvectors = model_pnlp_test5_testvectors()
    knownsols = model_pnlp_test5_knownsols()

    model_cfg = mpccmodel_construct_config(dimspec, defn; pdefns = pdefns, testvectors = testvectors, knownsols = knownsols)

    return model_cfg
end

