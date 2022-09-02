

using Symbolics

function model_pmpcc_kj6_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(2, 1, 2, 0, 2, 1, 0)
    return dimspec
end


function model_pmpcc_kj6_defn(x, pr, ps)
    label_model = "Example 6 from Kungurtsev and Jaeschke"

    f = exp(-x[1]+x[2])
    label_f = "Exp of x1 and x2 so descent in SE direction"

    ce = Vector{Num}()
    label_ce = Vector{String}([])

    ci = Vector{Num}(undef, 2)
    ci[1] = ((x[1]-2)^2 + (x[2]+1)^2 -2pr[1] -6)
    ci[2] = (1 - x[1])
    label_ci = Vector{String}([
            "Circle that comes in from right",
            "Static vertical constraint"
        ])

    F = Matrix{Num}(undef, 2, 1)
    F[1,1] = x[1]
    F[2,1] = x[2]
    label_F = Vector{String}([ "x1 and x2 axes complement" ])

    return MPCCDefinition(f, ce, ci, F, label_model, label_f, label_ce, label_ci, label_F)
end



function model_pmpcc_kj6_parameterisations(t)
    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t ]),
            (-5.0, 2.0),
            "Standard"    
        )

    return Vector{MPCCParameterisationDefn}( [ defn1 ] )
end



function model_pmpcc_kj6_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function model_pmpcc_kj6_knownsols()

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




function model_pmpcc_kj6_build()
    dimspec = model_pmpcc_kj6_dimspec()

    (x, pr, ps, t) = mpccmodel_build_sym_nums(dimspec)
    defn = model_pmpcc_kj6_defn(x, pr, ps)

    pdefns = model_pmpcc_kj6_parameterisations(t)
    testvectors = model_pmpcc_kj6_testvectors()
    knownsols = model_pmpcc_kj6_knownsols()

    model_cfg = mpccmodel_construct_config(dimspec, defn; pdefns = pdefns, testvectors = testvectors, knownsols = knownsols)

    return model_cfg
end

