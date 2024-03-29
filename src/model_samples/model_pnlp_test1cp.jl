

# constraints
# c1 = x.^2 + 1;
# c2 = 0*x + 1;
# c3 = exp(x);
# c4 = -x.^4 + 1;
# c5 = (1/4)*x.*(x-2).*(x+2) + 1;

# c10 = 0*x + 1 - e1;
# c11 = 5*x + 1 - e2;

# fc20 = @(r, s) (1/3^2)*(r*cos(pi/3) + (s-1)*sin(pi/3)).^2 + (r*sin(pi/3)-(s-1)*cos(pi/3)).^2 - (1/16);

function model_pnlp_test1cp_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(2, 0, 0, 0, 8, 2, 0)
    return dimspec
end



function model_pnlp_test1cp_defn(x, pr, ps)
    label_model = "Active Set Test 1 CP mod (NLP only)"

    f = exp(0.01*(x[1]^2 + x[2]^2 -sin(x[1]*x[2])))
    label_f = "Exp quadratic slanted objective"

    ce = Vector{Num}()
    label_ce = Vector{String}([])

    ci = Vector{Num}(undef, 8)
    ci[1] = 1 - x[1]^2 + x[2]
    ci[2] = x[2] - 1
    ci[3] = exp(x[1])
    ci[4] = x[1]^4 + x[2] - 1
    ci[5] = x[2] -(1/4)*x[1]*(x[1]-2)*(x[1]+2) - 1
    ci[6] = x[2] + 1 -pr[1]
    ci[7] = 5*x[2] + 1 -pr[2]
    ci[8] = (1/16) - (1/3^2)*(x[1]*cos(pi/3) + (x[2]-1)*sin(pi/3)).^2 - (x[1]*sin(pi/3)-(x[2]-1)*cos(pi/3)).^2
    label_ci = Vector{String}([
            "Quadratic (active)",
            "Horizonal line (active)",
            "Exp (active)",
            "-x^4 (active)",
            "Cubic fn (active)",
            "Horizontal line (not active for epsilon > 0)",
            "Straight diagonal line (not active for epsilon > 0)",
            "Enclosing rotated ellipse"
        ])

    F = Matrix{Num}(undef, 0, 0)
    label_F = Vector{String}([])

    return MPCCDefinition(f, ce, ci, F, label_model, label_f, label_ce, label_ci, label_F)
end




function model_pnlp_test1cp_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function model_pnlp_test1cp_knownsols()

    function local_test1cp_knownsol_1(pr::Vector{R}, ps::Vector{Int64}) where {R <: Real}
        x_sol = Vector{R}([0.0, 1.0])
        f_sol = exp(0.01*(x_sol[1].^2 + x_sol[2].^2))
        return (x_sol, f_sol)
    end

    return [ local_test1cp_knownsol_1 ]
end



function model_pnlp_test1cp_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t, t ]),
            (0.1, 0.0),
            "Standard"    
        )

    return  Vector{MPCCParameterisationDefn}( [ defn1 ] )
end






function model_pnlp_test1cp_build()
    dimspec = model_pnlp_test1cp_dimspec()

    (x, pr, ps, t) = mpccmodel_build_sym_nums(dimspec)
    defn = model_pnlp_test1cp_defn(x, pr, ps)

    pdefns = model_pnlp_test1cp_parameterisations(t)
    testvectors = model_pnlp_test1cp_testvectors()
    knownsols = model_pnlp_test1cp_knownsols()

    model_cfg = mpccmodel_construct_config(dimspec, defn; pdefns = pdefns, testvectors = testvectors, knownsols = knownsols)

    return model_cfg
end

