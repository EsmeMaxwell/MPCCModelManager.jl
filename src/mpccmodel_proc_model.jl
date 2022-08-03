


"""
    mpccmodel_build_fn_from_defn(dimspec, defnnum, x, pr, ps)

Assemble `f`, `ce`, `ci`, and `F` Julia native functions from Symbolic
definition and returns both non-indexed and indexed variants. Generally no need
to call directly as this is called from `mpccmodel_load_defn_from_file()`.

**NOTE:** To understand how this works, one should also play with a few examples
of`build_function()` in Symbolics.jl.  The output from `build_function()` is a
2-element vector whereby the first element is the standard function adn the
second is a mutating version of the same function. Hence the formulation in the
code here.

# Arguments
- `dimspec::MPCCDimSpec`: The dimension specifications for the proram.
- `defnnum::MPCCDefinition`: The Symbolics.jl Num definitions of the
  expressions.
- `x`, `pr`, `ps`: Symbolics.jl Num variables for x, pr, ps.

"""
function mpccmodel_build_fn_from_defn(dimspec::MPCCDimSpec, defnnum::MPCCDefinition, x, pr, ps)
    
    @unpack n, q, l, me, mi = dimspec

    # Before building any functions, we need to setup the Num for Fq (product down columns)
    Fq = Vector{Num}(undef, q)
    for lp_q=1:q
        Fq[lp_q] = prod(defnnum.F[:, lp_q])
    end

    # Objective fn
    f_fn = build_function([ defnnum.f ], x, pr, ps; expression=Val{false}, linenumbers=false)

    # Whole vector / array functions
    ce_fn = build_function(defnnum.ce, x, pr, ps; expression=Val{false}, linenumbers=false)
    ci_fn = build_function(defnnum.ci, x, pr, ps; expression=Val{false}, linenumbers=false)
    if ( l == 0 || q == 0 )
        # Workaround for wee bug in Symbolics when trying to compile empty matrix
        # it produces keech. So we hand it an empty vector instead.
        F_fn = build_function(Vector{Num}([]), x, pr, ps; expression=Val{false}, linenumbers=false)
    else
        F_fn = build_function(defnnum.F, x, pr, ps; expression=Val{false}, linenumbers=false)
    end
    Fq_fn = build_function(Fq, x, pr, ps; expression=Val{false}, linenumbers=false)

    # Indexed functions (need to separate std from mutating versions with [1] and [2])
    ce_i_fn = Vector{Tuple{Function, Function}}(undef, me)
    for lp_ce=1:me
        ce_i_fn[lp_ce] = build_function([ defnnum.ce[lp_ce] ], x, pr, ps; expression=Val{false}, linenumbers=false)
    end
    ce_i_fn_std = [ ce_i_fn[lp_ce][1] for lp_ce=1:me ]
    ce_i_fn_mut = [ ce_i_fn[lp_ce][2] for lp_ce=1:me ]

    ci_i_fn = Vector{Tuple{Function, Function}}(undef, mi)
    for lp_ci=1:mi
        ci_i_fn[lp_ci] = build_function([ defnnum.ci[lp_ci] ], x, pr, ps; expression=Val{false}, linenumbers=false)
    end
    ci_i_fn_std = [ ci_i_fn[lp_ci][1] for lp_ci=1:mi ]
    ci_i_fn_mut = [ ci_i_fn[lp_ci][2] for lp_ci=1:mi ]
    
    F_i_fn = Matrix{Tuple{Function, Function}}(undef, l, q)
    for lp_q=1:q
        for lp_l=1:l
            F_i_fn[lp_l, lp_q] = build_function([ defnnum.F[lp_l, lp_q] ], x, pr, ps; expression=Val{false}, linenumbers=false)
        end
    end
    F_i_fn_std = [ F_i_fn[lp_l, lp_q][1] for lp_l=1:l, lp_q=1:q ]
    F_i_fn_mut = [ F_i_fn[lp_l, lp_q][2] for lp_l=1:l, lp_q=1:q ]

    Fq_i_fn = Vector{Tuple{Function, Function}}(undef, q)
    for lp_q=1:q
        Fq_i_fn[lp_q] = build_function([ Fq[lp_q] ], x, pr, ps; expression=Val{false}, linenumbers=false)
    end
    Fq_i_fn_std = [ Fq_i_fn[lp_q][1] for lp_q=1:q ]
    Fq_i_fn_mut = [ Fq_i_fn[lp_q][2] for lp_q=1:q ]

    
    return MPCCFunctions(
            f_fn[1], f_fn[2],
            ce_fn[1], ce_fn[2],
            ce_i_fn_std, ce_i_fn_mut,
            ci_fn[1], ci_fn[2],
            ci_i_fn_std, ci_i_fn_mut,
            F_fn[1], F_fn[2],
            F_i_fn_std, F_i_fn_mut,
            Fq_fn[1], Fq_fn[2],
            Fq_i_fn_std, Fq_i_fn_mut
        )
end




"""
    mpccmodel_build_parameterisation(defns, t)

Build parameterisation from `defns` and Symbolics Num `t`, returns struct with
Julia native functions including derivatives.

"""
function mpccmodel_build_parameterisation(defns::Vector{MPCCParameterisationDefn}, t::Num)

    len_defns = length(defns)

    fns = Vector{MPCCParameterisationFunctions}(undef, len_defns)
    bf_results = Vector{Any}(undef, len_defns)
    bf_dt_results = Vector{Any}(undef, len_defns)
    sym_prdt = Vector{Vector{Num}}(undef, len_defns)
    for lp_i=1:len_defns
        # Calc derivatives
        sym_prdt[lp_i] = Symbolics.derivative(defns[lp_i].pr, t)

        # Compile them
        bf_results[lp_i] = Symbolics.build_function(defns[lp_i].pr, t; expression=Val{false})
        bf_dt_results[lp_i] = Symbolics.build_function(sym_prdt[lp_i], t; expression=Val{false})
        fns[lp_i] = MPCCParameterisationFunctions(
                bf_results[lp_i][1],
                bf_results[lp_i][2],
                bf_dt_results[lp_i][1],
                bf_dt_results[lp_i][2]
            )
    end

    return MPCCParameterisations(t, defns, fns)
end





function mpccmodel_load_defn_from_file(model_id::String)

    # TODO some sanity checking on model_id

    # Filename of the model we want to import, include() it
    src_filename = joinpath(@__DIR__, "models", "model_spec_$(model_id).jl")
    println("Loading model_spec_$(model_id).jl")
    include(src_filename)

    # We need to call two specific functions in the included file.
    # Need to setup these function variables from strings then we'll later use invokelatest
    src_fn_dimspec = "mm_spec_" * model_id * "_dimspec"
    src_fn_defn = "mm_spec_" * model_id * "_defn"
    # src_fn_nzmask = "mm_spec_" * model_id * "_nzmask"
    src_fn_testvectors = "mm_spec_" * model_id * "_testvectors"
    src_fn_knownsols = "mm_spec_" * model_id * "_knownsols"
    src_fn_param_defns = "mm_spec_" * model_id * "_parameterisations"
    fn_dimspec = getfield(@__MODULE__, Symbol(src_fn_dimspec))
    fn_defn = getfield(@__MODULE__, Symbol(src_fn_defn))
    # fn_nzmask = getfield(@__MODULE__, Symbol(src_fn_nzmask))
    fn_testvectors = getfield(@__MODULE__, Symbol(src_fn_testvectors))
    fn_knownsols = getfield(@__MODULE__, Symbol(src_fn_knownsols))
    fn_param_defns = getfield(@__MODULE__, Symbol(src_fn_param_defns))

    # We're getting the dimspec and setting up Symbolics.jl variables
    dimspec = Base.invokelatest(fn_dimspec)
    @variables x[1:dimspec.n] pr[1:dimspec.r] ps[1:dimspec.s] t
    
    # Grab the model definition, expressed in Symbolics.jl Num
    defn = Base.invokelatest(fn_defn, x, pr, ps)

    # Grab the param definitions    
    param_defns = Base.invokelatest(fn_param_defns, t)

    # Build callable functions for the basic model: f, ce, ci, F
    fns = mpccmodel_build_fn_from_defn(dimspec, defn, x, pr, ps)

    # Build parameterisations
    parameterisations = mpccmodel_build_parameterisation(param_defns, t)

    # nzmask = Base.invokelatest(fn_nzmask)
    testvectors = Base.invokelatest(fn_testvectors)
    knownsols = Base.invokelatest(fn_knownsols)

    return MPCCModelConfig(
        Symbolics.scalarize(x), Symbolics.scalarize(pr), Symbolics.scalarize(ps),
        dimspec,
        defn,
        fns,
        # nzmask,
        testvectors,
        knownsols,
        parameterisations
    )
end



# Because JuMP is a bit shit, we need to construct scalar functions for it
# We fix the parameter set here too.
function mpccmodel_build_fixed_jump_fns(model::MPCCModelConfig, pr0::Vector{R}, ps0::Vector{Int64}) where {R <: Real}

    # Create suitable objective function for use in JuMP, fixing parameters
    jump_f(y::T...) where {T <: Real} = model.fns.f(collect(y), pr0, ps0)[1]  # This would otherwise return a vector!

    # Create suitable penalty objective function for use in JuMP, fixing parameters
    jump_f_pen = function(ρ::R, y::T...) where {R <: Real, T <: Real}
        x = collect(y)

        # Initial f value
        f_val = model.fns.f(collect(y), pr0, ps0)[1]

        # Add penalty
        dot_prod = zero(T)
        for lp_q=1:model.dimspec.q
            cnstr_prod = one(T)
            for lp_l=1:model.dimspec.l
                cnstr_prod *= abs(model.fns.F_i[lp_l, lp_q](collect(y), pr0, ps0)[1])
            end
            dot_prod += cnstr_prod
        end
        
        return (f_val + ρ*dot_prod)
    end

    # Create suitable constraint functions for use in JuMP, fixing parameters
    jump_ce = Vector{Function}(undef, model.dimspec.me)
    for lp_me=1:model.dimspec.me
        jump_ce[lp_me] = function (y::T...) where {T <: Real}
            return model.fns.ce_i[lp_me](collect(y), pr0, ps0)[1]
        end
    end

    # Create suitable constraint functions for use in JuMP, fixing parameters
    jump_ci = Vector{Function}(undef, model.dimspec.mi)
    for lp_mi=1:model.dimspec.mi
        jump_ci[lp_mi] = function (y::T...) where {T <: Real}
            return model.fns.ci_i[lp_mi](collect(y), pr0, ps0)[1]
        end
    end
    
    # Create suitable constraint functions for use in JuMP, fixing parameters
    jump_F = Matrix{Function}(undef, model.dimspec.l, model.dimspec.q)
    for lp_q=1:model.dimspec.q
        for lp_l=1:model.dimspec.l
            jump_F[lp_l, lp_q] = function (y::T...) where {T <: Real}
                return model.fns.F_i[lp_l, lp_q](collect(y), pr0, ps0)[1]
            end
        end
    end

    return MPCCJumpFixedFunctions(jump_f, jump_f_pen, jump_ce, jump_ci, jump_F)
end
