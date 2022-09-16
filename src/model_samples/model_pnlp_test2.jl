

# exp(0.01*(v^2 - (1/2)*v*w + w^2)) + (1/17)*w*sin(3*(v + w));

# % constraints
# c1 = (x+0.5).^2 + 2;
# c2 = cos(3*x) + 1.5;
# c3 = (x+0.5).^4 + 2.0 - e1;


function model_pnlp_test2_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(2, 0, 0, 0, 3, 1, 0)
    return dimspec
end



function model_pnlp_test2_defn(x, pr, ps)
    label_model = "Active Set Test 2 (NLP only)"

    f = exp(0.01*(x[1]^2 - (1/2)*x[1]*x[2] + x[2]^2)) + (1/17)*x[2]*sin(3*(x[1] + x[2]))
    label_f = "Mental objective"

    ce = Vector{Num}()
    label_ce = Vector{String}([])

    ci = Vector{Num}(undef, 3)
    ci[1] = x[2] - (x[1]+0.5)^2 - 2
    ci[2] = x[2] - cos(3x[1]) - 1.5
    ci[3] = x[2] - (x[1]+0.5)^4 - 2.0 + pr[1]
    label_ci = Vector{String}([
            "Quadratic",
            "Cos",
            "Quartic"
        ])

    F = Matrix{Num}(undef, 0, 0)
    label_F = Vector{String}([])

    return MPCCDefinition(f, ce, ci, F, label_model, label_f, label_ce, label_ci, label_F)
end




function model_pnlp_test2_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function model_pnlp_test2_knownsols()

    # function local_test2_knownsol_1(pr::Vector{R}, ps::Vector{Int64}) where {R <: Real}
    #     x_sol = Vector{R}([0.0, 1.0])
    #     f_sol = exp(0.01*(x_sol[1].^2 + x_sol[2].^2))
    #     return (x_sol, f_sol)
    # end

    return Vector{Function}([ ])
end



function model_pnlp_test2_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t ]),
            (0.005, 0.0),
            "Standard"    
        )

    return  Vector{MPCCParameterisationDefn}( [ defn1 ] )
end






function model_pnlp_test2_build()
    dimspec = model_pnlp_test2_dimspec()

    (x, pr, ps, t) = mpccmodel_build_sym_nums(dimspec)
    defn = model_pnlp_test2_defn(x, pr, ps)

    pdefns = model_pnlp_test2_parameterisations(t)
    testvectors = model_pnlp_test2_testvectors()
    knownsols = model_pnlp_test2_knownsols()

    model_cfg = mpccmodel_construct_config(dimspec, defn; pdefns = pdefns, testvectors = testvectors, knownsols = knownsols)

    return model_cfg
end

