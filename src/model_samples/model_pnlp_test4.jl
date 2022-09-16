
# c1 = 0*x + 1;
# c2 = x + 1;

# c10 = 2*x + 1 -e1;
# c11 = 3*x + 1 -3*e1;
# c12 = 4*x + 1 -5*e1;

# c20 = -0.5*x + 1 -e2;
# c21 = -x + 1 -e2;
# c22 = -7*x + 1 -3*e2;



function model_pnlp_test4_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(2, 0, 0, 0, 8, 2, 0)
    return dimspec
end



function model_pnlp_test4_defn(x, pr, ps)
    label_model = "Active Set Test 4 (NLP only)"

    f = exp(0.01*(x[1]^2 + x[2]^2))
    label_f = "Exp quadratic objective "

    ce = Vector{Num}()
    label_ce = Vector{String}([])

    ci = Vector{Num}(undef, 8)
    ci[1] = x[2] - 1
    ci[2] = x[2] - x[1] - 1
    ci[3] = x[2] - 2x[1] - 1 + pr[1]
    ci[4] = x[2] - 3x[1] - 1 + 3pr[1]
    ci[5] = x[2] - 4x[1] - 1 + 5pr[1]
    ci[6] = x[2] + 0.5x[1] -1 + pr[2]
    ci[7] = x[2] + x[1] -1 + pr[2]
    ci[8] = x[2] + 7x[1] -1 + 3pr[2]
    label_ci = Vector{String}([
            "Straight line (active)",
            "Straight line (active)",
            "Straight line (inactive)",
            "Straight line (inactive)",
            "Straight line (inactive)",
            "Straight line (inactive)",
            "Straight line (inactive)",
            "Straight line (inactive)"
        ])

    F = Matrix{Num}(undef, 0, 0)
    label_F = Vector{String}([])

    return MPCCDefinition(f, ce, ci, F, label_model, label_f, label_ce, label_ci, label_F)
end




function model_pnlp_test4_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function model_pnlp_test4_knownsols()

    function local_test4_knownsol_1(pr::Vector{R}, ps::Vector{Int64}) where {R <: Real}
        x_sol = Vector{R}([0.0, 1.0])
        f_sol = exp(0.01*(x_sol[1].^2 + x_sol[2].^2))
        return (x_sol, f_sol)
    end

    return [ local_test4_knownsol_1 ]
end



function model_pnlp_test4_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t, t ]),
            (0.1, 0.0),
            "Standard"    
        )

    return  Vector{MPCCParameterisationDefn}( [ defn1 ] )
end






function model_pnlp_test4_build()
    dimspec = model_pnlp_test4_dimspec()

    (x, pr, ps, t) = mpccmodel_build_sym_nums(dimspec)
    defn = model_pnlp_test4_defn(x, pr, ps)

    pdefns = model_pnlp_test4_parameterisations(t)
    testvectors = model_pnlp_test4_testvectors()
    knownsols = model_pnlp_test4_knownsols()

    model_cfg = mpccmodel_construct_config(dimspec, defn; pdefns = pdefns, testvectors = testvectors, knownsols = knownsols)

    return model_cfg
end

