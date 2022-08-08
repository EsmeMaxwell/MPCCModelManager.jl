

function model_pmpcc_kj3_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(2, 1, 2, 0, 0, 1, 0)
    return dimspec
end



function model_pmpcc_kj3_defn(x, pr, ps)
    label_model = "Example 6 from Kungurtsev and Jaeschke"

    f = (x[1] - pr[1])^2 + (x[2] - pr[1])^2
    label_f = "Quadratic objective with minimum dependent on parameter"    

    ce = Vector{Num}()
    label_ce = Vector{String}([])

    ci = Vector{Num}()
    label_ci = Vector{String}([])

    F = Matrix{Num}(undef, 2, 1)
    F[1,1] = x[1]
    F[2,1] = x[2]
    label_F = Vector{String}([ "x1 and x2 axes complement" ])

    return MPCCDefinition(f, ce, ci, F, label_model, label_f, label_ce, label_ci, label_F)
end




function model_pmpcc_kj3_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function model_pmpcc_kj3_knownsols()

    function local_kj3_knownsol_1(pr::Vector{R}, ps::Vector{Int64}) where {R <: Real}
        @assert ( 1 == length(pr) && 0 == length(ps) )
        t = pr[1]

        if ( t <= 0 )
            x_sol = [0, 0]
            f_sol = 2(-t)^2
        else
            x_sol = [t, 0]
            f_sol = (x_sol[1] - t)^2 + (x_sol[2] - t)^2
        end

        return (x_sol, f_sol)
    end

    function local_kj3_knownsol_2(pr::Vector{R}, ps::Vector{Int64}) where {R <: Real}
        @assert ( 1 == length(pr) && 0 == length(ps) )
        t = pr[1]

        if ( t <= 0 )
            x_sol = [0, 0]
            f_sol = 2(-t)^2
        else
            x_sol = [0, t]
            f_sol = (x_sol[1] - t)^2 + (x_sol[2] - t)^2
        end

        return (x_sol, f_sol)
    end

    return [ local_kj3_knownsol_1, local_kj3_knownsol_2 ]
end



function model_pmpcc_kj3_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t ]),
            (-5.0, 5.0),
            "Standard"    
        )

    return  Vector{MPCCParameterisationDefn}( [ defn1 ] )
end






function model_pmpcc_kj3_build()
    dimspec = model_pmpcc_kj3_dimspec()

    (x, pr, ps, t) = mpccmodel_build_sym_nums(dimspec)
    defn = model_pmpcc_kj3_defn(x, pr, ps)

    pdefns = model_pmpcc_kj3_parameterisations(t)
    testvectors = model_pmpcc_kj3_testvectors()
    knownsols = model_pmpcc_kj3_knownsols()

    model_cfg = mpccmodel_construct_config(dimspec, defn; pdefns = pdefns, testvectors = testvectors, knownsols = knownsols)

    return model_cfg
end

