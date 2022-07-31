

# Should take a specification of indexing and return the sparsity and functions.  User may then call multiple times.

# 20220726 TODO : for some inexplicable reason, there still type instability in hessci, etc, despite explicit declaration in comprehension.
# Had tested similar ideas on the REPL and it works fine, so no idea why it's messed-up here. But it's enough to get empty spare arrays
# printing in the REPL.


function mpccmodel_setup_symdiff_sparse(config::MPCCModelConfig, Tf::Type)

    # Tf == "Type force"
    @assert Tf <: Number "Tf should be a subtype of Number"
    @assert isconcretetype(Tf) "Tf must be a concrete type"
    
    # NOTE indexed versions will be slower as they call, sequentially, a vector of
    # already built functions, in sequence.

    @unpack dimspec, defn, fns = config
    @unpack n, q, l, me, mi, r, s = config.dimspec

    s_x = config.x
    s_pr = config.pr
    s_ps = config.ps

    # TEMPORARY FIX 202207 until Symbolics.jl fixes sparsity detection code for Hessians.
    # We must substitute scalar variables in place of the pr[i], then do the calcs.
    # Otherwise it won't work

    @assert (r == 1) "Temporary fix restricts to one parameter"
    @variables p1
    dict_subst = Dict([ s_pr[1] => p1 ])
    fixup_defn_f_num = Symbolics.substitute(defn.f, dict_subst)
    fixup_defn_ce_num = Symbolics.substitute(defn.ce, dict_subst)
    fixup_defn_ci_num = Symbolics.substitute(defn.ci, dict_subst)
    fixup_defn_F_num = Symbolics.substitute(defn.F, dict_subst)


    # We need to setup the num arrays for the sparse jacobian and hessians
    # Note that although SparseMatrixCSC is used consistently below, most of these
    # are vectors; Symbolics spits out the mtrx type and can't be bothered converting

    # grad entries
    gradf_sparse_num = Symbolics.sparsejacobian([ fixup_defn_f_num ], s_x)
    gradce_sparse_num = Symbolics.sparsejacobian(fixup_defn_ce_num, s_x)
    gradci_sparse_num = Symbolics.sparsejacobian(fixup_defn_ci_num, s_x)
    gradF_sparse_num = Matrix{SparseMatrixCSC}(undef, l, q)
    for lp_q=1:q
        for lp_l=1:l
            gradF_sparse_num[lp_l, lp_q] = Symbolics.sparsejacobian([ fixup_defn_F_num[lp_l, lp_q] ], s_x)
        end
    end

    sparsity_gradf = Symbolics.jacobian_sparsity([ fixup_defn_f_num ], s_x)
    sparsity_gradce = Symbolics.jacobian_sparsity(fixup_defn_ce_num, s_x)
    sparsity_gradci = Symbolics.jacobian_sparsity(fixup_defn_ci_num, s_x)
    sparsity_gradF = Matrix{SparseMatrixCSC}(undef, l, q)
    for lp_q=1:q
        for lp_l=1:l
            sparsity_gradF[lp_l, lp_q] = Symbolics.jacobian_sparsity([ fixup_defn_F_num[lp_l, lp_q] ], s_x)
        end
    end

    println("gradf: ", gradf_sparse_num)
    println("gradce: ", gradce_sparse_num)
    println("gradci: ", gradci_sparse_num)
    println("gradF: ", gradF_sparse_num)
    println("sparsity_gradf: ", sparsity_gradf)
    println("sparsity_gradce: ", sparsity_gradce)
    println("sparsity_gradci: ", sparsity_gradci)
    println("sparsity_gradF: ", sparsity_gradF)
    println("")

    # hess entries
    # Second TEMPORARY FIX: hessian doesn't work well in Symbolics.jl, so we
    # simulate by doing two jacobian calcs.

    jacf_dense_num = Symbolics.jacobian([ fixup_defn_f_num ], s_x)
    hessf_sparse_num = Symbolics.sparsejacobian(vec(jacf_dense_num), s_x)
    sparsity_hessf = Symbolics.jacobian_sparsity(vec(jacf_dense_num), s_x)
    # hessf_sparse_num = Symbolics.sparsehessian(defn.f, s_x)

    hessce_sparse_num = Vector{SparseMatrixCSC}(undef, me)
    sparsity_hessce = Vector{SparseMatrixCSC}(undef, me)
    for lp_ce=1:me
        jacce_dense_num = Symbolics.jacobian([ fixup_defn_ce_num[lp_ce] ], s_x)
        hessce_sparse_num[lp_ce] = Symbolics.sparsejacobian(vec(jacce_dense_num), s_x)        
        sparsity_hessce[lp_ce] = Symbolics.jacobian_sparsity(vec(jacce_dense_num), s_x)        
        # hessce_sparse_num[lp_ce] = Symbolics.sparsehessian(defn.ce[lp_ce], s_x)
    end

    hessci_sparse_num = Vector{SparseMatrixCSC}(undef, mi)
    sparsity_hessci = Vector{SparseMatrixCSC}(undef, mi)
    for lp_ci=1:mi
        jacci_dense_num = Symbolics.jacobian([ fixup_defn_ci_num[lp_ci] ], s_x)
        hessci_sparse_num[lp_ci] = Symbolics.sparsejacobian(vec(jacci_dense_num), s_x)
        sparsity_hessci[lp_ci] = Symbolics.jacobian_sparsity(vec(jacci_dense_num), s_x)
        # hessci_sparse_num[lp_ci] = Symbolics.sparsehessian(defn.ci[lp_ci], s_x)
    end

    hessF_sparse_num = Matrix{SparseMatrixCSC}(undef, l, q)
    sparsity_hessF = Matrix{SparseMatrixCSC}(undef, l, q)
    for lp_q=1:q
        for lp_l=1:l
            jacF_dense_num = Symbolics.jacobian([ fixup_defn_F_num[lp_l, lp_q] ], s_x)
            hessF_sparse_num[lp_l, lp_q] = Symbolics.sparsejacobian(vec(jacF_dense_num), s_x)
            sparsity_hessF[lp_l, lp_q] = Symbolics.jacobian_sparsity(vec(jacF_dense_num), s_x)
            # hessF_sparse_num[lp_l, lp_q] = Symbolics.sparsehessian(defn.F[lp_l, lp_q], s_x)
        end
    end

    println("hessf: ", hessf_sparse_num)
    println("hessce: ", hessce_sparse_num)
    println("hessci: ", hessci_sparse_num)
    println("hessF: ", hessF_sparse_num)
    println("sparsity_hessf: ", sparsity_hessf)
    println("sparsity_hessce: ", sparsity_hessce)
    println("sparsity_hessci: ", sparsity_hessci)
    println("sparsity_hessF: ", sparsity_hessF)
    println("")

    # dp entries  TEMP edit to math issues above re Symbolics.jl
    fdp_sparse_num = Vector{Num}(undef, r)
    for lp_r=1:r
        fdp_sparse_num[lp_r] = Symbolics.derivative(fixup_defn_f_num, p1; simplify=true)

        # fdp_sparse_num[lp_r] = Symbolics.derivative(defn.f, s_pr[lp_r]; simplify=true)
    end
    cedp_sparse_num = Vector{Vector{Num}}(undef, r)
    for lp_r=1:r        
        cedp_sparse_num[lp_r] = Symbolics.derivative(fixup_defn_ce_num, p1; simplify=true)
        # cedp_sparse_num[lp_r] = Symbolics.derivative(defn.ce, s_pr[lp_r]; simplify=true)
    end
    cidp_sparse_num = Vector{Vector{Num}}(undef, r)
    for lp_r=1:r   
        cidp_sparse_num[lp_r] = Symbolics.derivative(fixup_defn_ci_num, p1; simplify=true)
        # cidp_sparse_num[lp_r] = Symbolics.derivative(defn.ci, s_pr[lp_r]; simplify=true)
    end
    Fdp_sparse_num = Vector{Matrix{Num}}(undef, r)
    for lp_r=1:r        
        Fdp_sparse_num[lp_r] = Symbolics.derivative(fixup_defn_F_num, p1; simplify=true)
        # Fdp_sparse_num[lp_r] = Symbolics.derivative(defn.F, s_pr[lp_r]; simplify=true)
    end

    println("fdp: ", fdp_sparse_num)
    println("cedp: ", cedp_sparse_num)
    println("cidp: ", cidp_sparse_num)
    println("Fdp: ", Fdp_sparse_num)


    # graddp entries 
    gradfdp_sparse_num = Vector{SparseMatrixCSC}(undef, r)
    sparsity_gradfdp =  Vector{SparseMatrixCSC}(undef, r)
    for lp_r=1:r
        gradfdp_sparse_num[lp_r] = Symbolics.sparsejacobian([fdp_sparse_num[lp_r]], s_x)
        sparsity_gradfdp[lp_r] = Symbolics.jacobian_sparsity([fdp_sparse_num[lp_r]], s_x)
    end
    gradcedp_sparse_num = Vector{SparseMatrixCSC}(undef, r)
    sparsity_gradcedp = Vector{SparseMatrixCSC}(undef, r)
    for lp_r=1:r
        gradcedp_sparse_num[lp_r] = Symbolics.sparsejacobian( cedp_sparse_num[lp_r], s_x )
        sparsity_gradcedp[lp_r] = Symbolics.jacobian_sparsity( cedp_sparse_num[lp_r], s_x )
    end
    gradcidp_sparse_num = Vector{SparseMatrixCSC}(undef, r)
    sparsity_gradcidp = Vector{SparseMatrixCSC}(undef, r)
    for lp_r=1:r
        gradcidp_sparse_num[lp_r] = Symbolics.sparsejacobian( cidp_sparse_num[lp_r], s_x )
        sparsity_gradcidp[lp_r] = Symbolics.jacobian_sparsity( cidp_sparse_num[lp_r], s_x )
    end
    gradFdp_sparse_num = Vector{Matrix{SparseMatrixCSC}}(undef, r)
    sparsity_gradFdp = Vector{Matrix{SparseMatrixCSC}}(undef, r)
    for lp_r=1:r
        gradFdp_sparse_num[lp_r] = Matrix{SparseMatrixCSC}(undef, l, q)
        sparsity_gradFdp[lp_r] = Matrix{SparseMatrixCSC}(undef, l, q)
        for lp_q=1:q
            for lp_l=1:l
                gradFdp_sparse_num[lp_r][lp_l, lp_q] = Symbolics.sparsejacobian( [ Fdp_sparse_num[lp_r][lp_l, lp_q] ], s_x )
                sparsity_gradFdp[lp_r][lp_l, lp_q] = Symbolics.jacobian_sparsity( [ Fdp_sparse_num[lp_r][lp_l, lp_q] ], s_x )
            end
        end
    end

    println("gradfdp: ", gradfdp_sparse_num)
    println("gradcedp: ", gradcedp_sparse_num)
    println("gradcidp: ", gradcidp_sparse_num)
    println("gradFdp: ", gradFdp_sparse_num)
    println("sparsity_gradfdp: ", sparsity_gradfdp)
    println("sparsity_gradcedp: ", sparsity_gradcedp)
    println("sparsity_gradcidp: ", sparsity_gradcidp)
    println("sparsity_gradFdp: ", sparsity_gradFdp)

    # Compiled versions
    gradf_sparse_fn = build_function(gradf_sparse_num, s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false)
    gradce_sparse_fn = build_function(gradce_sparse_num, s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false)
    gradci_sparse_fn = build_function(gradci_sparse_num, s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false)
    gradF_sparse_fn = [ build_function(gradF_sparse_num[lp_l, lp_q], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false) for lp_l=1:l, lp_q=1:q ]
    # gradF_sparse_fn = Matrix(undef, l, q)       # TODO specify type here    
    # for lp_q=1:q
    #     for lp_l=1:l
    #         gradF_sparse_fn[lp_l, lp_q] = build_function(gradF_sparse_num[lp_l, lp_q], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false)
    #     end
    # end

    hessf_sparse_fn = build_function(hessf_sparse_num, s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false)
    hessce_sparse_fn = [ build_function(hessce_sparse_num[lp_me], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false) for lp_me=1:me ]
    hessci_sparse_fn = [ build_function(hessci_sparse_num[lp_mi], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false) for lp_mi=1:mi ]
    hessF_sparse_fn = [ build_function(hessF_sparse_num[lp_l, lp_q], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false) for lp_l=1:l, lp_q=1:q ]

    fdp_sparse_fn= [ build_function(fdp_sparse_num[lp_r], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false) for lp_r=1:r ]
    cedp_sparse_fn= [ build_function(cedp_sparse_num[lp_r], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false) for lp_r=1:r ]
    cidp_sparse_fn= [ build_function(cidp_sparse_num[lp_r], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false) for lp_r=1:r ]
    Fdp_sparse_fn= [ build_function(Fdp_sparse_num[lp_r], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false) for lp_r=1:r ]

    gradfdp_sparse_fn = [ build_function(gradfdp_sparse_num[lp_r], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false) for lp_r=1:r ]
    gradcedp_sparse_fn = [ build_function(gradcedp_sparse_num[lp_r], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false) for lp_r=1:r ]
    gradcidp_sparse_fn = [ build_function(gradcidp_sparse_num[lp_r], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false) for lp_r=1:r ]
    gradFdp_sparse_fn = [ [ build_function(gradcidp_sparse_num[lp_r][lp_l, lp_q], s_x, s_pr, s_ps; expression=Val{false}, linenumbers=false) for lp_l=1:l, lp_q=1:q ] for lp_r=1:r ]

    # f, ce, ci, F all return usual dense vectors/matrices because they are, well, dense.

    # Full f functions (really just converting a 1 element vector to a scalar)
    capture_f = config.fns.f
    function local_f(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        # return fns.f(x, pr, ps)[1]
        return capture_f(x, pr, ps)[1]
    end



    # Full ce functions
    function local_ce(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}   # ::Vector{S}
        return fns.ce(x, pr, ps)
    end




    # Full ci functions
    function local_ci(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}  # ::Vector{S}
        return fns.ci(x, pr, ps)
    end



    # Full F functions
    function local_F(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}  # ::Matrix{S}
        return fns.F(x, pr, ps)
    end



    # gradf functions
    function local_gradf(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return gradf_sparse_fn[1](x, pr, ps)
    end


    # gradce functions
    function local_gradce(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return gradce_sparse_fn[1](x, pr, ps)
    end


    # gradci functions
    function local_gradci(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return gradci_sparse_fn[1](x, pr, ps)
    end
    

    # gradF functions
    function local_gradF(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}                
        return mm_sym_sp_gradF(dimspec, gradF_sparse_fn, x, pr, ps)
    end


    # hessf functions
    function local_hessf(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return hessf_sparse_fn[1](x, pr, ps)
    end


    # hessce functions
    function local_hessce(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}                
        return mm_sym_sp_hessce(dimspec, hessce_sparse_fn, x, pr, ps)
    end


    # hessci functions
    function local_hessci(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return mm_sym_sp_hessci(dimspec, Tf, hessci_sparse_fn, x, pr, ps)
    end


    # hessF functions
    function local_hessF(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return mm_sym_sp_hessF(dimspec, hessF_sparse_fn, x, pr, ps)
    end


    # fdp functions
    function local_fdp(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return mm_sym_sp_fdp(dimspec, fdp_sparse_fn, x, pr, ps)
    end


    # cedp functions
    function local_cedp(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return mm_sym_sp_cedp(dimspec, cedp_sparse_fn, x, pr, ps)
    end


    # cedp functions
    function local_cidp(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return mm_sym_sp_cidp(dimspec, cidp_sparse_fn, x, pr, ps)
    end


    # Fdp functions
    function local_Fdp(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return mm_sym_sp_Fdp(dimspec, Fdp_sparse_fn, x, pr, ps)
    end


    # gradfdp functions
    function local_gradfdp(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return mm_sym_sp_gradfdp(dimspec, gradfdp_sparse_fn, x, pr, ps)
    end


    # gradcedp functions
    function local_gradcedp(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return mm_sym_sp_gradcedp(dimspec, gradcedp_sparse_fn, x, pr, ps)
    end


    # gradcedp functions
    function local_gradcidp(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return mm_sym_sp_gradcidp(dimspec, gradcidp_sparse_fn, x, pr, ps)
    end


    # gradFdp functions
    function local_gradFdp(x::Vector{S}, pr::Vector{T}, ps::Vector{Int64}) where {S <: Real, T <: Real}
        return mm_sym_sp_gradFdp(dimspec, gradFdp_sparse_fn, x, pr, ps)
    end

    return MPCCModelSparseSym(
        config,
        sparsity_gradf,
        sparsity_gradce,
        sparsity_gradci,
        sparsity_gradF,
        sparsity_hessf,
        sparsity_hessce,
        sparsity_hessci,
        sparsity_hessF,
        sparsity_gradfdp,
        sparsity_gradcedp,
        sparsity_gradcidp,
        sparsity_gradFdp,
        gradf_sparse_num,
        gradce_sparse_num,
        gradci_sparse_num,
        gradF_sparse_num,
        hessf_sparse_num,
        hessce_sparse_num,
        hessci_sparse_num,
        hessF_sparse_num,
        fdp_sparse_num,
        cedp_sparse_num,
        cidp_sparse_num,
        Fdp_sparse_num,
        gradfdp_sparse_num,
        gradcedp_sparse_num,
        gradcidp_sparse_num,
        gradFdp_sparse_num,
        local_f, missing,
        local_ce, missing,
        local_ci, missing,
        local_F, missing,
        local_gradf, missing,
        local_gradce, missing,
        local_gradci, missing,
        local_gradF, missing,
        local_hessf, missing,
        local_hessce, missing,
        local_hessci, missing,
        local_hessF, missing,
        local_fdp, missing,
        local_cedp, missing,
        local_cidp, missing,
        local_Fdp, missing,
        local_gradfdp, missing,
        local_gradcedp, missing,
        local_gradcidp, missing,
        local_gradFdp, missing

    )

end


# NOTE don't restrict the ::AbstractMatrix or AbstractVector since Symbolics makes Tuples of RGFs


function mm_sym_sp_gradF(dimspec::MPCCDimSpec, fn_gradF::AbstractMatrix, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @unpack n, l, q = dimspec
    gradF = [ fn_gradF[lp_l, lp_q][1](x, pr, ps) for lp_l in 1:l, lp_q in 1:q ]        # Call each (l,q) function, the [1] is the non-mutating version
    return gradF
end



function mm_sym_sp_hessce(dimspec::MPCCDimSpec, fn_hessce::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @unpack me = dimspec
    hessce = [ fn_hessce[lp_me][1](x, pr, ps) for lp_me in 1:me ]
    for lp_me=1:me
        nt_types = sp_getparamtypes(hessce[lp_me])
        if !(nt_types.T <: Number)
            hessce[lp_me] = SparseMatrixCSC{Number, nt_types.S}(hessce[lp_me])
        end
    end 
    return hessce
end


function mm_sym_sp_hessci(dimspec::MPCCDimSpec, Tf::Type, fn_hessci::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @unpack mi = dimspec
    hessci = SparseMatrixCSC{Tf, Int64}[ fn_hessci[lp_mi][1](x, pr, ps) for lp_mi in 1:mi ]
    return hessci
end


function mm_sym_sp_hessF(dimspec::MPCCDimSpec, fn_hessF::AbstractMatrix, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @unpack n, l, q = dimspec
    hessF = [ fn_hessF[lp_l, lp_q][1](x, pr, ps) for lp_l in 1:l, lp_q in 1:q ]        # Call each (l,q) function, the [1] is the non-mutating version
    for lp_q=1:q
        for lp_l=1:l
            nt_types = sp_getparamtypes(hessF[lp_l, lp_q])
            println("at $lp_l, $lp_q")
            println(nt_types)
            if !(nt_types.T <: Number)
                println("fixing...")
                hessF[lp_l, lp_q] = SparseMatrixCSC{Number, nt_types.S}(hessF[lp_l, lp_q])
            end
        end
    end 
    return hessF
end



function mm_sym_sp_fdp(dimspec::MPCCDimSpec, fn_fdp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @unpack r = dimspec
    fdp = [ fn_fdp[lp_r][1](x, pr, ps) for lp_r in 1:r ]
    return fdp
end


function mm_sym_sp_cedp(dimspec::MPCCDimSpec, fn_cedp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @unpack r = dimspec
    cedp = [ fn_cedp[lp_r][1](x, pr, ps) for lp_r in 1:r ]
    return cedp
end


function mm_sym_sp_cidp(dimspec::MPCCDimSpec, fn_cidp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @unpack r = dimspec
    cidp = [ fn_cidp[lp_r][1](x, pr, ps) for lp_r in 1:r ]
    return cidp
end


function mm_sym_sp_Fdp(dimspec::MPCCDimSpec, fn_Fdp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @unpack r = dimspec
    Fdp = [ fn_Fdp[lp_r][1](x, pr, ps) for lp_r in 1:r ]
    return Fdp
end



function mm_sym_sp_gradfdp(dimspec::MPCCDimSpec, fn_gradfdp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @unpack r = dimspec
    gradfdp = [ fn_gradfdp[lp_r][1](x, pr, ps) for lp_r in 1:r ]
    return gradfdp
end


function mm_sym_sp_gradcedp(dimspec::MPCCDimSpec, fn_gradcedp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @unpack r = dimspec
    gradcedp = [ fn_gradcedp[lp_r][1](x, pr, ps) for lp_r in 1:r ]
    return gradcedp
end


function mm_sym_sp_gradcidp(dimspec::MPCCDimSpec, fn_gradcidp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @unpack r = dimspec
    gradcidp = [ fn_gradcidp[lp_r][1](x, pr, ps) for lp_r in 1:r ]
    return gradcidp
end


function mm_sym_sp_gradFdp(dimspec::MPCCDimSpec, fn_gradFdp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @unpack l, q, r = dimspec
    gradFdp = [ [ fn_gradFdp[lp_r][lp_l, lp_q][1](x, pr, ps) for lp_l=1:l, lp_q=1:q ] for lp_r in 1:r ]
    return gradFdp
end


# gradfdp_sparse_num = Vector{SparseMatrixCSC}(undef, r)
# for lp_r=1:r
#     gradfdp_sparse_num[lp_r] = Symbolics.sparsejacobian([fdp_sparse_num[lp_r]], s_x)
# end
# gradcedp_sparse_num = Vector{SparseMatrixCSC}(undef, r)
# for lp_r=1:r
#     gradcedp_sparse_num[lp_r] = Symbolics.sparsejacobian( cedp_sparse_num[lp_r], s_x )
# end
# gradcidp_sparse_num = Vector{SparseMatrixCSC}(undef, r)
# for lp_r=1:r
#     gradcidp_sparse_num[lp_r] = Symbolics.sparsejacobian( cidp_sparse_num[lp_r], s_x )
# end
# gradFdp_sparse_num = Vector{Matrix{SparseMatrixCSC}}(undef, r)
# for lp_r=1:r
#     gradFdp_sparse_num[lp_r] = Matrix{SparseMatrixCSC}(undef, l, q)
#     for lp_q=1:q
#         for lp_l=1:l
#             gradFdp_sparse_num[lp_r][lp_l, lp_q] = Symbolics.sparsejacobian( [ Fdp_sparse_num[lp_r][lp_l, lp_q] ], s_x )
#         end
#     end
# end

