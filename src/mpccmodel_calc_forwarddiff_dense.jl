


# Finish checcking type stability....






# RESOLVED (faster as it is, and type stable) Should the ::Matrix{P} annotation be on the rhs, or move to the lhs?


# --------------------------- older stuff?

# TODO: have just removed some converts() because the promote_type obviously doesn't work when input empty

# Note: Might need to relax the type on mutation arguments for the mm_ functions, not the exported functions though.
# Might also need to relax the typing on arguments for the local functions in mm_ functions.

# Need to ensure convert() is used consistently

#------------------------------

    




"""
    mpccmodel_setup_forwarddiff_dense(config::MPCCModelConfig)


Accepts model config and returns a struct with all the nice Julia functions to
calculate Jacobians, Hessians, and parametric derivatives using
ForwardDiff.jl.

Both standard and mutating functions are constructed. Most of the functions
accept an index as the last argument to specify a subset of the constraints.
For F constraints, this should be a vector of either tuple pairs or
CartesianIndexes.
"""
function mpccmodel_setup_forwarddiff_dense(config::MPCCModelConfig)

    @unpack dimspec, fns = config

    @unpack f, f!, ce, ce!, ce_i, ce_i!, ci, ci!, ci_i, ci_i!, F, F!, F_i, F_i!, Fq, Fq!, Fq_i, Fq_i! = fns

    # Full f functions
    function local_f(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return f(x, pr, ps)[1]
        # return fns.f(x, pr, ps)[1]
    end

    function local_f!(out_f::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        # fns.f!(@view(out_f[1]), x, pr, ps)
        f!(out_f, x, pr, ps)
        return nothing
    end
  

    # Full ce functions
    function local_ce(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}   
        return ce(x, pr, ps)
    end

    function local_ce!(out_ce::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        # fns.ce!(@view(out_ce[:]), x, pr, ps)
        ce!(out_ce, x, pr, ps)
        return nothing
    end


    # Indexed ce functions
    function local_ce(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64}) where {S <: Real, T <: Real} 
        return mm_fd_dn_ce_i(dimspec, ce_i, x, pr, ps, idxs_ce)
    end
    
    function local_ce!(out_ce_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_ce_i!(out_ce_i, dimspec, ce_i!, x, pr, ps, idxs_ce)
        return nothing
    end


    # Full ci functions
    function local_ci(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}  
        return ci(x, pr, ps)
        # return fns.ci(x, pr, ps)
    end

    function local_ci!(out_ci::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        # fns.ci!(@view(out_ci[:]), x, pr, ps)
        ci!(out_ci, x, pr, ps)
        return nothing
    end


    # Indexed ci functions
    function local_ci(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64}) where {S <: Real, T <: Real}  
        return mm_fd_dn_ci_i(dimspec, ci_i, x, pr, ps, idxs_ci)
    end
    
    function local_ci!(out_ci_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_ci!(out_ci_i, dimspec, ci_i!, x, pr, ps, idxs_ci)
        return nothing
    end


    # Full F functions
    function local_F(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}  # ::Matrix{S}
        # return convert(Matrix{promote_type(S, T)}, fns.F(x, pr, ps))
        return F(x, pr, ps)
    end

    function local_F!(out_F::AbstractMatrix, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        # fns.F!(@view(out_F[:,:]), x, pr, ps)
        F!(out_F, x, pr, ps)
        return nothing
    end


    # Indexed F functions; NOTE, these yeild a vector in the order that the (l,q) index tuples were in!
    function local_F(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}}) where {S <: Real, T <: Real}  
        return mm_fd_dn_F_i(dimspec, F_i, x, pr, ps, idxs_F)
    end

    function local_F(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{Tuple{Int64,Int64}}) where {S <: Real, T <: Real}  
        cart_idxs_F = map(x -> CartesianIndex{2}(x), idxs_F)
        return mm_fd_dn_F_i(dimspec, F_i, x, pr, ps, cart_idxs_F)
    end

    function local_F!(out_F_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_F_i!(out_F_i, dimspec, F_i!, x, pr, ps, idxs_F)
        return nothing
    end    

    function local_F!(out_F_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{Tuple{Int64,Int64}})::Nothing where {S <: Real, T <: Real}
        cart_idxs_F = map(x -> CartesianIndex{2}(x), idxs_F)
        mm_fd_dn_F_i!(out_F_i, dimspec, F_i!, x, pr, ps, cart_idxs_F)
        return nothing
    end


    # Full Fq functions
    function local_Fq(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}  
        # return convert(Vector{promote_type(S, T)}, fns.Fq(x, pr, ps))
        return Fq(x, pr, ps)
    end

    function local_Fq!(out_Fq::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        # fns.Fq!(@view(out_Fq[:]), x, pr, ps)
        Fq!(out_Fq, x, pr, ps)
        return nothing
    end


    # Indexed Fq functions
    function local_Fq(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_Fq::AbstractVector{Int64}) where {S <: Real, T <: Real}  
        return mm_fd_dn_Fq_i(dimspec, Fq_i, x, pr, ps, idxs_Fq)
    end
    
    function local_Fq!(out_Fq_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_Fq::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_Fq!(out_Fq_i, dimspec, Fq_i!, x, pr, ps, idxs_Fq)
        return nothing
    end



    # gradf functions
    function local_gradf(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}  
        return mm_fd_dn_gradf(dimspec, f, x, pr, ps)
    end

    function local_gradf!(out_gradf::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradf!(out_gradf, dimspec, f, x, pr, ps)
        return nothing
    end


    # Full gradce functions
    function local_gradce(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_gradce(dimspec, ce, x, pr, ps)
    end

    function local_gradce!(out_gradce::AbstractMatrix, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradce!(out_gradce, dimspec, ce!, x, pr, ps)
        return nothing
    end


    # Indexed gradce functions
    function local_gradce(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_gradce_i(dimspec, ce_i, x, pr, ps, idxs_ce)
    end

    function local_gradce!(out_gradce_i::AbstractMatrix, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradce_i!(out_gradce_i, dimspec, ce_i!, x, pr, ps, idxs_ce)
        return nothing
    end


    # Full gradci functions
    function local_gradci(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_gradci(dimspec, ci, x, pr, ps)
    end

    function local_gradci!(out_gradci::AbstractMatrix, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradci!(out_gradci, dimspec, ci!, x, pr, ps)
        return nothing
    end


    # Indexed gradci functions
    function local_gradci(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_gradci_i(dimspec, ci_i, x, pr, ps, idxs_ci)
    end

    function local_gradci!(out_gradci_i::AbstractMatrix, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradci_i!(out_gradci_i, dimspec, ci_i!, x, pr, ps, idxs_ci)
        return nothing
    end


 
    # Full gradF functions
    function local_gradF(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_gradF(dimspec, F, x, pr, ps)
    end

    function local_gradF!(out_gradF::AbstractMatrix, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradF!(out_gradF, dimspec, F!, x, pr, ps)
        return nothing
    end

    # Indexed gradF functions
    # NOTE, these yeild a vector in the order that the (l,q) index tuples were in!
    function local_gradF(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}}) where {S <: Real, T <: Real}
        return mm_fd_dn_gradF_i(dimspec, F_i, x, pr, ps, idxs_F)
    end

    function local_gradF(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{Tuple{Int64,Int64}}) where {S <: Real, T <: Real}
        new_idxs_F = map(x -> CartesianIndex{2}(x), idxs_F)
        return mm_fd_dn_gradF_i(dimspec, F_i, x, pr, ps, new_idxs_F)
    end

    function local_gradF!(out_gradF_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradF_i!(out_gradF_i, dimspec, F_i!, x, pr, ps, idxs_F)
        return nothing
    end    

    function local_gradF!(out_gradF_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{Tuple{Int64,Int64}})::Nothing where {S <: Real, T <: Real}
        new_idxs_F = map(x -> CartesianIndex{2}(x), idxs_F)
        mm_fd_dn_gradF_i!(out_gradF_i, dimspec, F_i!, x, pr, ps, new_idxs_F)
        return nothing
    end



    # Full gradFq functions
    function local_gradFq(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_gradFq(dimspec, Fq, x, pr, ps)
    end

    function local_gradFq!(out_gradF::AbstractMatrix, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradFq!(out_gradF, dimspec, Fq!, x, pr, ps)
        return nothing
    end

    # Indexed gradF functions
    function local_gradFq(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_Fq::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_gradFq_i(dimspec, Fq_i, x, pr, ps, idxs_Fq)
    end

    function local_gradFq!(out_gradFq_i::AbstractMatrix, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_Fq::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradFq_i!(out_gradFq_i, dimspec, Fq_i!, x, pr, ps, idxs_Fq)
        return nothing
    end


    # Hessian f functions (we don't do indexed version for this)
    function local_hessf(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_hessf(dimspec, f, x, pr, ps)
    end

    function local_hessf!(out_hessf::AbstractMatrix, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_hessf!(out_hessf, dimspec, f, x, pr, ps)
        return nothing
    end


    # Full hessce functions
        # Vector of length me of 2d Hessian matrices
    function local_hessce(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_hessce(dimspec, ce_i, x, pr, ps)
    end

    # function local_hessce!(out_hessce::AbstractVector{AbstractMatrix{R}}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {R <: Real, S <: Real, T <: Real}
    function local_hessce!(out_hessce::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_hessce!(out_hessce, dimspec, ce_i, x, pr, ps)
        return nothing
    end

    # Indexed hessce functions
    function local_hessce(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_hessce_i(dimspec, ce_i, x, pr, ps, idxs_ce)
    end

    function local_hessce!(out_hessce_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_hessce_i!(out_hessce_i, dimspec, ce_i, x, pr, ps, idxs_ce)
        return nothing
    end




    # Full hessci functions
    # Vector of length me of 2d Hessian matricis
    function local_hessci(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_hessci(dimspec, ci_i, x, pr, ps)
    end

    function local_hessci!(out_hessci::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_hessci!(out_hessci, dimspec, ci_i, x, pr, ps)
        return nothing
    end

    # Indexed hessci functions
    function local_hessci(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_hessci_i(dimspec, ci_i, x, pr, ps, idxs_ci)
    end

    function local_hessci!(out_hessci_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_hessci_i!(out_hessci_i, dimspec, ci_i, x, pr, ps, idxs_ci)
        return nothing
    end
    
    


    # Full hessF functions
    function local_hessF(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}  # ::Matrix{Matrix{S}}
        return mm_fd_dn_hessF(dimspec, F_i, x, pr, ps)
    end

    function local_hessF!(out_hessF::AbstractMatrix, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_hessF!(out_hessF, dimspec, F_i, x, pr, ps)
        return nothing
    end

    # Indexed hessF functions
    # NOTE, these yeild a vector in the order that the (l,q) index tuples were in!
    function local_hessF(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}}) where {S <: Real, T <: Real}
        return mm_fd_dn_hessF_i(dimspec, F_i, x, pr, ps, idxs_F)
    end

    function local_hessF(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{Tuple{Int64,Int64}}) where {S <: Real, T <: Real}
        new_idxs_F = map(x -> CartesianIndex{2}(x), idxs_F)
        return mm_fd_dn_hessF_i(dimspec, F_i, x, pr, ps, new_idxs_F)
    end
    
    function local_hessF!(out_hessF_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_hessF_i!(out_hessF_i, dimspec, F_i, x, pr, ps, idxs_F)
        return nothing
    end

    function local_hessF!(out_hessF_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{Tuple{Int64,Int64}})::Nothing where {S <: Real, T <: Real}
        new_idxs_F = map(x -> CartesianIndex{2}(x), idxs_F)
        mm_fd_dn_hessF_i!(out_hessF_i, dimspec, F_i, x, pr, ps, new_idxs_F)
        return nothing
    end




    # Full hessFq functions
    
    function local_hessFq(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_hessFq(dimspec, Fq_i, x, pr, ps)
    end

    function local_hessFq!(out_hessFq::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        # mm_fd_dn_hessFq!(out_hessFq, dimspec, Fq_i, x, pr, ps)
        error("not implemented")
        return nothing
    end

    # Indexed hessFq functions
    function local_hessFq(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_Fq::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_hessFq_i(dimspec, Fq_i, x, pr, ps, idxs_Fq)
    end

    function local_hessFq!(out_hessFq_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_Fq::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        # mm_fd_dn_hessFq_i!(out_hessFq_i, dimspec, Fq_i, x, pr, ps, idxs_Fq)
        error("not implemented")
        return nothing
    end
    





    # fdp functions
    function local_fdp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}   # ::AbstractVector{T}
        return mm_fd_dn_fdp(dimspec, f, x, pr, ps)
    end

    function local_fdp!(out_fdp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_fdp!(out_fdp, dimspec, f, x, pr, ps)
        return nothing
    end


    # Full cedp functions
    function local_cedp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}  # ::Vector{Vector{S}}
        return mm_fd_dn_cedp(dimspec, ce, x, pr, ps)
    end

    function local_cedp!(out_cedp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_cedp!(out_cedp, dimspec, ce!, x, pr, ps)
        return nothing
    end


    # Indexed cedp functions
    function local_cedp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64}) where {S <: Real, T <: Real}  # Vector{Vector} ret vector of len r, of vectors of len idx
        return mm_fd_dn_cedp_i(dimspec, ce_i, x, pr, ps, idxs_ce)
    end

    function local_cedp!(out_cedp_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_cedp_i!(out_cedp_i, dimspec, ce_i!, x, pr, ps, idxs_ce)
        return nothing
    end


    # Full cidp functions
    function local_cidp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}  # ::Vector{Vector{S}}
        return mm_fd_dn_cidp(dimspec, ci, x, pr, ps)
    end

    function local_cidp!(out_cidp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_cidp!(out_cidp, dimspec, ci!, x, pr, ps)
        return nothing
    end


    # Indexed cidp functions
    function local_cidp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64}) where {S <: Real, T <: Real}  # Vector{Vector} ret vector of len r, of vectors of len idx
        return mm_fd_dn_cidp_i(dimspec, ci_i, x, pr, ps, idxs_ci)
    end

    function local_cidp!(out_cidp_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_cidp_i!(out_cidp_i, dimspec, ci_i!, x, pr, ps, idxs_ci)
        return nothing
    end




    # Full Fdp functions
    function local_Fdp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}   # ::Vector{Matrix{T}}
        return mm_fd_dn_Fdp(dimspec, F, x, pr, ps)
    end

    function local_Fdp!(out_Fdp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_Fdp!(out_Fdp, dimspec, F!, x, pr, ps)
        return nothing
    end

    # Indexed Fdp functions
    # NOTE, these yeild a vector in the order that the (l,q) index tuples were in!
    function local_Fdp(x::AbstractVector{T}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}})::Vector{Vector{T}} where {T <: Real}
        return mm_fd_dn_Fdp_i(dimspec, F_i, x, pr, ps, idxs_F)
    end

    function local_Fdp(x::AbstractVector{T}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{Tuple{Int64,Int64}})::Vector{Vector{T}} where {T <: Real}
        new_idxs_F = map(x -> CartesianIndex{2}(x), idxs_F)
        return mm_fd_dn_Fdp_i(dimspec, F_i, x, pr, ps, new_idxs_F)
    end
    
    function local_Fdp!(out_Fdp_i::AbstractVector, x::AbstractVector{T}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}})::Nothing where {T <: Real}
        mm_fd_dn_Fdp_i!(out_Fdp_i, dimspec, F_i!, x, pr, ps, idxs_F)
        return nothing
    end

    function local_Fdp!(out_Fdp_i::AbstractVector, x::AbstractVector{T}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{Tuple{Int64,Int64}})::Nothing where {T <: Real}
        new_idxs_F = map(x -> CartesianIndex{2}(x), idxs_F)
        mm_fd_dn_Fdp_i!(out_Fdp_i, dimspec, F_i!, x, pr, ps, new_idxs_F)
        return nothing
    end


    # Full gradfdp functions
    function local_gradfdp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}   # ::Vector{Vector{T}}
        return mm_fd_dn_gradfdp(dimspec, f, x, pr, ps)
    end

    function local_gradfdp!(out_Fdp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradfdp!(out_Fdp, dimspec, f, x, pr, ps)
        return nothing
    end



    # Full gradcedp functions
    function local_gradcedp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}  # ::Vector{Matrix{S}}
        return mm_fd_dn_gradcedp(dimspec, ce, x, pr, ps)
    end

    function local_gradcedp!(out_gradcedp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradcedp!(out_gradcedp, dimspec, ce!, x, pr, ps)
        return nothing
    end


    # Indexed gradcedp functions
    function local_gradcedp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64})::Vector{Matrix{T}} where {S <: Real, T <: Real}
        return mm_fd_dn_gradcedp_i(dimspec, ce_i, x, pr, ps, idxs_ce)
    end

    function local_gradcedp!(out_gradcedp_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradcedp_i!(out_gradcedp_i, dimspec, ce_i!, x, pr, ps, idxs_ce)
        return nothing
    end


    # Full gradcidp functions
    function local_gradcidp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_gradcidp(dimspec, ci, x, pr, ps)
    end

    function local_gradcidp!(out_gradcidp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradcidp!(out_gradcidp, dimspec, ci!, x, pr, ps)
        return nothing
    end


    # Indexed gradcidp functions
    function local_gradcidp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64})::Vector{Matrix{S}} where {S <: Real, T <: Real}
        return mm_fd_dn_gradcidp_i(dimspec, ci_i, x, pr, ps, idxs_ci)
    end

    function local_gradcidp!(out_gradcidp_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradcidp_i!(out_gradcidp_i, dimspec, ci_i!, x, pr, ps, idxs_ci)
        return nothing
    end


    # Full gradFdp functions
    function local_gradFdp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
        return mm_fd_dn_gradFdp(dimspec, F_i, x, pr, ps)
    end

    function local_gradFdp!(out_gradFdp::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradFdp!(out_gradFdp, dimspec, F_i, x, pr, ps)
        return nothing
    end

    
    # Indexed gradFdp functions
    function local_gradFdp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}}) where {S <: Real, T <: Real}
        return mm_fd_dn_gradFdp_i(dimspec, F_i, x, pr, ps, idxs_F)
    end

    function local_gradFdp(x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{Tuple{Int64,Int64}}) where {S <: Real, T <: Real}
        new_idxs_F = map(x -> CartesianIndex{2}(x), idxs_F)
        return mm_fd_dn_gradFdp_i(dimspec, F_i, x, pr, ps, new_idxs_F)
    end
    
    function local_gradFdp!(out_gradFdp_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}})::Nothing where {S <: Real, T <: Real}
        mm_fd_dn_gradFdp_i!(out_gradFdp_i, dimspec, F_i, x, pr, ps, idxs_F)
        return nothing
    end

    function local_gradFdp!(out_gradFdp_i::AbstractVector, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{Tuple{Int64,Int64}})::Nothing where {S <: Real, T <: Real}
        new_idxs_F = map(x -> CartesianIndex{2}(x), idxs_F)
        mm_fd_dn_gradFdp_i!(out_gradFdp_i, dimspec, F_i, x, pr, ps, new_idxs_F)
        return nothing
    end



    return MPCCModelDenseForwardDiff(   config,
                        local_f, local_f!,
                        local_ce, local_ce!,
                        local_ci, local_ci!,
                        local_F, local_F!,
                        local_Fq, local_Fq!,
                        local_gradf, local_gradf!,
                        local_gradce, local_gradce!,
                        local_gradci, local_gradci!,
                        local_gradF, local_gradF!,
                        local_gradFq, local_gradFq!,
                        local_hessf, local_hessf!,
                        local_hessce, local_hessce!,
                        local_hessci, local_hessci!,
                        local_hessF, local_hessF!,
                        local_hessFq, local_hessFq!,
                        local_fdp, local_fdp!,
                        local_cedp, local_cedp!,
                        local_cidp, local_cidp!,
                        local_Fdp, local_Fdp!,
                        local_gradfdp, local_gradfdp!,
                        local_gradcedp, local_gradcedp!,
                        local_gradcidp, local_gradcidp!,
                        local_gradFdp, local_gradFdp!
                        # local_do_self_test,
                        # nz_mask 
                    )
end






# Indexed ce functions
function mm_fd_dn_ce_i(dimspec::MPCCDimSpec, ce_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real} 
    P = promote_type(S, T)
    res_ce_i = P[ ce_i[idxs_ce[lp_ce]](x, pr, ps)[1] for lp_ce in eachindex(idxs_ce) ]
    return res_ce_i
end

function mm_fd_dn_ce_i!(out_ce_i::AbstractArray, dimspec::MPCCDimSpec, ce_i!::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    for lp_ce in eachindex(idxs_ce)
        ce_i![lp_ce](@view(out_ce_i[lp_ce:lp_ce]), x, pr, ps)        # As first argument we pass a 1-element vector, so that it can be mutated
    end
    return nothing
end



# Indexed ci functions
function mm_fd_dn_ci_i(dimspec::MPCCDimSpec, ci_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}
    P = promote_type(S, T)
    res_ci_i = P[ ci_i[idxs_ci[lp_ci]](x, pr, ps)[1] for lp_ci in eachindex(idxs_ci) ]
    return res_ci_i
end

function mm_fd_dn_ci_i!(out_ci_i::AbstractArray, dimspec::MPCCDimSpec, ci_i!::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    for lp_ci in eachindex(idxs_ci)
        ci_i![lp_ci](@view(out_ci_i[lp_ci:lp_ci]), x, pr, ps)        # As first argument we pass a 1-element vector, so that it can be mutated
    end
    return nothing
end



# Indexed F functions; NOTE, these yeild a vector in the order that the (l,q) index tuples were in!
function mm_fd_dn_F_i(dimspec::MPCCDimSpec, F_i::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}}) where {F <: Function, S <: Real, T <: Real}
    P = promote_type(S, T)
    res_F_i = P[ F_i[idxs_F[lp_F]](x, pr, ps)[1] for lp_F in eachindex(idxs_F) ]
    return res_F_i
end

function mm_fd_dn_F_i!(out_F_i::AbstractArray, dimspec::MPCCDimSpec, F_i!::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}})::Nothing where {F <: Function, S <: Real, T <: Real}
    for lp_F in eachindex(idxs_F)
        # This is a bawhair away from being a contestant in an "obfuscated Julia" competition...   
        F_i![idxs_F[lp_F]](@view(out_F_i[lp_F:lp_F]), x, pr, ps)
    end
    return nothing
end




# Indexed Fq functions
function mm_fd_dn_Fq_i(dimspec::MPCCDimSpec, Fq_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_Fq::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real} 
    P = promote_type(S, T)
    res_Fq_i = P[ Fq_i[idxs_Fq[lp_Fq]](x, pr, ps)[1] for lp_Fq in eachindex(idxs_Fq) ]
    return res_Fq_i
end

function mm_fd_dn_Fq!(out_Fq_i::AbstractArray, dimspec::MPCCDimSpec, Fq_i!::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_Fq::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    for lp_Fq in eachindex(idxs_Fq)
        Fq_i![lp_Fq](@view(out_Fq_i[lp_Fq:lp_Fq]), x, pr, ps)        # As first argument we pass a 1-element vector, so that it can be mutated
    end
    return nothing
end





# gradf functions

function mm_fd_dn_gradf(dimspec::MPCCDimSpec, f::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    P = promote_type(S, T)
    local_f(z) = f(z, pr, ps)[1]
    return ForwardDiff.gradient(local_f, x)::Vector{P}      # NOTE 20220731: this type annotation is cludgy. fix me.
end


function mm_fd_dn_gradf!(out_gradf::AbstractArray, dimspec::MPCCDimSpec, f::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    local_f(z) = f(z, pr, ps)[1]
    ForwardDiff.gradient!(out_gradf, local_f, x)
    return nothing
end





# Full gradce functions
function mm_fd_dn_gradce(dimspec::MPCCDimSpec, ce::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    P = promote_type(S, T)

    # Shortcut for empty result: ensures type stable (fails assert with FD), and gets dimensions correct
    if 0 == dimspec.me
        return Matrix{P}(undef, 0, dimspec.n)
    end

    # Create closure and do FD calc
    local_ce = (z::AbstractArray) -> ce(z, pr, ps)
    return ForwardDiff.jacobian(local_ce, x)::Matrix{P}
end


function mm_fd_dn_gradce!(out_gradce::AbstractArray, dimspec::MPCCDimSpec, ce!::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    @unpack n, me = dimspec
    local_ce! = (y::AbstractArray, z::AbstractArray) -> ce!(y, z, pr, ps)

    # ForwardDiff asks for a temporary storage area, presumably so that caller can manage memory allocs
    # We do this here for now, perhaps move it upstream later.
    y = zeros(promote_type(S, T), me)
    ForwardDiff.jacobian!(out_gradce, local_ce!, y, x)
    return nothing
end


# Indexed gradce functions
function mm_fd_dn_gradce_i(dimspec::MPCCDimSpec, ce_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}
    P = promote_type(S, T)
    local_ce = (z::AbstractArray) -> mm_fd_dn_ce_i(dimspec, ce_i, z, pr, ps, idxs_ce)
    return ForwardDiff.jacobian(local_ce, x)::Matrix{P}
end


function mm_fd_dn_gradce_i!(out_gradce_i::AbstractArray, dimspec::MPCCDimSpec, ce_i!::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    @unpack n = dimspec    
    len_idxs_ce_i = length(idxs_ce)
    local_ce! = (y::AbstractArray, z::AbstractArray) -> mm_fd_dn_ce_i!(y, dimspec, ce_i!, z, pr, ps, idxs_ce)

    # ForwardDiff asks for a temporary storage area, presumably so that caller can manage memory allocs
    # We do this here for now, perhaps move it upstream later.
    y = zeros(promote_type(S, T), len_idxs_ce_i)
    ForwardDiff.jacobian!(out_gradce_i, local_ce!, y, x)
    return nothing    
end







# Full gradci functions
function mm_fd_dn_gradci(dimspec::MPCCDimSpec, ci::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    P = promote_type(S, T)

    # Shortcut for empty result: ensures type stable (fails assert with FD), and gets dimensions correct
    if 0 == dimspec.mi
        return Matrix{P}(undef, 0, dimspec.n)
    end

    # Create closure and do FD calc    
    local_ci = (z::AbstractArray) -> ci(z, pr, ps)
    return ForwardDiff.jacobian(local_ci, x)::Matrix{P}
end


function mm_fd_dn_gradci!(out_gradci::AbstractArray, dimspec::MPCCDimSpec, ci!::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    @unpack n, mi = dimspec
    local_ci! = (y::AbstractArray, z::AbstractArray) -> ci!(y, z, pr, ps)
    
    # ForwardDiff asks for a temporary storage area, presumably so that caller can manage memory allocs
    # We do this here for now, perhaps move it upstream later.
    y = zeros(promote_type(S, T), mi)
    ForwardDiff.jacobian!(out_gradci, local_ci!, y, x)
    return nothing
end


# Indexed gradci functions
function mm_fd_dn_gradci_i(dimspec::MPCCDimSpec, ci_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}
    P = promote_type(S, T)
    local_ci = (z::AbstractArray) -> mm_fd_dn_ci_i(dimspec, ci_i, z, pr, ps, idxs_ci)
    return ForwardDiff.jacobian(local_ci, x)::Matrix{P}
end


function mm_fd_dn_gradci_i!(out_gradci_i::AbstractArray, dimspec::MPCCDimSpec, ci_i!::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    @unpack n = dimspec    
    len_idxs_ci_i = length(idxs_ci)
    local_ci! = (y::AbstractArray, z::AbstractArray) -> mm_fd_dn_ci_i!(y, dimspec, ci_i!, z, pr, ps, idxs_ci)

    # ForwardDiff asks for a temporary storage area, presumably so that caller can manage memory allocs
    # We do this here for now, perhaps move it upstream later.
    y = zeros(promote_type(S, T), len_idxs_ci_i)
    ForwardDiff.jacobian!(out_gradci_i, local_ci!, y, x)
    return nothing    
end








# Full gradF functions
# # Matrix of dim l x q of n length vectors (the gradients wrt x)

function mm_fd_dn_gradF(dimspec::MPCCDimSpec, F::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}   # ::Matrix{Vector{S}}
    P = promote_type(S, T)
    @unpack n, l, q = dimspec

    # Shortcut for empty result: ensures type stable (fails assert with FD), and gets dimensions correct
    if ( 0 == dimspec.l || 0 == dimspec.q )
        return Matrix{Vector{P}}(undef, dimspec.l, dimspec.q)
    end

    # ForwardDiff.jacobian can work with multidimensional arrays, so we do that
    local_F = (z::AbstractArray) -> F(z, pr, ps)
    gradF_flat = ForwardDiff.jacobian(local_F, x)::Matrix{P}
    gradF = [ gradF_flat[l*(lp_q-1)+lp_l, :] for lp_l in 1:l, lp_q in 1:q ]
    return gradF
end


function mm_fd_dn_gradF!(out_gradF::AbstractArray, dimspec::MPCCDimSpec, F!::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    @unpack n, l, q = dimspec
    # ForwardDiff.jacobian can work with multidimensional arrays, so we do that
    local_F! = (y::AbstractArray, z::AbstractArray) -> F!(y, z, pr, ps)

    # Alloc these locally for now, should really be done by caller to realise efficiency gains
    y = zeros(promote_type(S, T), l, q)
    gradF_flat = zeros(promote_type(S, T), l*q, n)
    ForwardDiff.jacobian!(gradF_flat, local_F!, y, x)
    for lp_q=1:q
        for lp_l=1:l
            out_gradF[lp_l, lp_q] = gradF_flat[l*(lp_q-1)+lp_l, :]      # Surpringly, I think this is the simplest way to do this
        end
    end
    return nothing
end



# Indexed gradF functions
function mm_fd_dn_gradF_i(dimspec::MPCCDimSpec, F_i::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}}) where {F <: Function, S <: Real, T <: Real}
    P = promote_type(S, T)
    @unpack n, l, q = dimspec
    len_idxs_F_i = length(idxs_F)
    local_F = (z::AbstractArray) -> mm_fd_dn_F_i(dimspec, F_i, z, pr, ps, idxs_F)
    gradF_flat = ForwardDiff.jacobian(local_F, x)::Matrix{P}
    gradF = [ gradF_flat[lp_idx, :] for lp_idx=1:len_idxs_F_i ]
    return gradF
end


function mm_fd_dn_gradF_i!(out_gradF_i::AbstractArray, dimspec::MPCCDimSpec, F_i!::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}})::Nothing where {F <: Function, S <: Real, T <: Real}
    @unpack n, l, q = dimspec
    len_idxs_F_i = length(idxs_F)
    local_F! = (y::AbstractArray, z::AbstractArray) -> mm_fd_dn_F_i!(y, dimspec, F_i!, z, pr, ps, idxs_F)    
    # Alloc these locally for now, should really be done by caller to realise efficiency gains
    y = zeros(promote_type(S, T), len_idxs_F_i)
    gradF_flat = zeros(promote_type(S, T), len_idxs_F_i, n)
    ForwardDiff.jacobian!(gradF_flat, local_F!, y, x)
    for lp_idx=1:len_idxs_F_i
        out_gradF_i[lp_idx] = gradF_flat[lp_idx, :]
    end
    return nothing
end




# Full gradFq functions
function mm_fd_dn_gradFq(dimspec::MPCCDimSpec, Fq::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    P = promote_type(S, T)

    # Shortcut for empty result: ensures type stable (fails assert with FD), and gets dimensions correct
    if ( 0 == dimspec.l || 0 == dimspec.q )
        return Matrix{P}(undef, dimspec.q, dimspec.n)
    end
    
    local_Fq = (z::AbstractArray) -> Fq(z, pr, ps) 
    # return convert(Matrix{promote_type(S, T)}, ForwardDiff.jacobian(local_Fq, x))
    return ForwardDiff.jacobian(local_Fq, x)::Matrix{P}
end


function mm_fd_dn_gradFq!(out_gradFq::AbstractArray, dimspec::MPCCDimSpec, Fq!::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    @unpack n, q = dimspec
    local_Fq! = (y::AbstractArray, z::AbstractArray) -> Fq!(y, z, pr, ps)
    
    # ForwardDiff asks for a temporary storage area, presumably so that caller can manage memory allocs
    # We do this here for now, perhaps move it upstream later.
    y = zeros(promote_type(S, T), q)
    ForwardDiff.jacobian!(out_gradFq, local_Fq!, y, x)
    return nothing
end




# Indexed gradFq functions
function mm_fd_dn_gradFq_i(dimspec::MPCCDimSpec, Fq_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_Fq::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real} 
    P = promote_type(S, T)
    local_Fq = (z::AbstractArray) -> mm_fd_dn_Fq_i(dimspec, Fq_i, z, pr, ps, idxs_Fq)
    return ForwardDiff.jacobian(local_Fq, x)::Matrix{P}
end


function mm_fd_dn_gradFq_i!(out_gradFq_i::AbstractArray, dimspec::MPCCDimSpec, Fq_i!::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_Fq::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    @unpack n = dimspec    
    len_idxs_Fq_i = length(idxs_Fq)
    local_Fq! = (y::AbstractArray, z::AbstractArray) -> mm_fd_dn_Fq!(y, dimspec, Fq_i!, z, pr, ps, idxs_Fq)

    # ForwardDiff asks for a temporary storage area, presumably so that caller can manage memory allocs
    # We do this here for now, perhaps move it upstream later.
    y = zeros(promote_type(S, T), len_idxs_Fq_i)
    
    # Unfortuantely, forewarddiff.jacobian returns the transpose of what we want.
    # So, currently, we're stuck allocating again: TODO sort this
    ForwardDiff.jacobian!(out_gradFq_i, local_Fq!, y, x)
    return nothing        
end







# Hessian of f

function mm_fd_dn_hessf(dimspec::MPCCDimSpec, f::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    P = promote_type(S, T)
    local_f(z) = f(z, pr, ps)[1]
    return ForwardDiff.hessian(local_f, x)::Matrix{P}
end


function mm_fd_dn_hessf!(out_hessf::AbstractArray, dimspec::MPCCDimSpec, f::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    local_f(z) = f(z, pr, ps)[1]
    ForwardDiff.hessian!(out_hessf, local_f, x)
    return nothing
end







# Full hessce functions

function mm_fd_dn_hessce(dimspec::MPCCDimSpec, ce_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}
    # Call the indexed version...
    idxs_ce = collect(1:dimspec.me)
    return mm_fd_dn_hessce_i(dimspec, ce_i, x, pr, ps, idxs_ce::AbstractVector{Int64})    
end


function mm_fd_dn_hessce!(out_hessce::AbstractArray, dimspec::MPCCDimSpec, ce_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    # Call the indexed version...
    idxs_ce = collect(1:dimspec.me)
    mm_fd_dn_hessce_i!(out_hessce, dimspec, ce_i, x, pr, ps, idxs_ce)
    return nothing
end


# Indexed hessce functions

function mm_fd_dn_hessce_i(dimspec::MPCCDimSpec, ce_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}
    P = promote_type(S, T)
    @unpack n = dimspec

    len_idxs_ce_i = length(idxs_ce)
    local_ce_i = Vector{Function}(undef, len_idxs_ce_i)
    hessce_i = Vector{Matrix{T}}(undef, len_idxs_ce_i)

    for lp_ce in eachindex(idxs_ce)
        local_ce_i[lp_ce] = (z::AbstractArray) -> ce_i[idxs_ce[lp_ce]](z, pr, ps)[1]
        hessce_i[lp_ce] = ForwardDiff.hessian(local_ce_i[lp_ce], x)::Matrix{P}      
    end
    return hessce_i
end


function mm_fd_dn_hessce_i!(out_hessce_i::AbstractArray, dimspec::MPCCDimSpec, ce_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    len_idxs_ce_i = length(idxs_ce)
    local_ce_i = Vector{Function}(undef, len_idxs_ce_i)

    for lp_ce in eachindex(idxs_ce)
        local_ce_i[lp_ce] = (z::AbstractArray) -> ce_i[idxs_ce[lp_ce]](z, pr, ps)[1]
        ForwardDiff.hessian!(@view(out_hessce_i[lp_ce][:]), local_ce_i[lp_ce], x)
    end
    return nothing
end





# Full hessci functions

function mm_fd_dn_hessci(dimspec::MPCCDimSpec, ci_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}
    # Call the indexed version...
    idxs_ci = collect(1:dimspec.mi)
    return mm_fd_dn_hessci_i(dimspec, ci_i, x, pr, ps, idxs_ci::AbstractVector{Int64})    
end


function mm_fd_dn_hessci!(out_hessci::AbstractArray, dimspec::MPCCDimSpec, ci_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    # Call the indexed version...
    idxs_ci = collect(1:dimspec.mi)
    mm_fd_dn_hessci_i!(out_hessci, dimspec, ci_i, x, pr, ps, idxs_ci)
    return nothing
end


# Indexed hessci functions

function mm_fd_dn_hessci_i(dimspec::MPCCDimSpec, ci_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real} 
    P = promote_type(S, T)
    @unpack n = dimspec

    len_idxs_ci_i = length(idxs_ci)
    local_ci_i = Vector{Function}(undef, len_idxs_ci_i)
    hessci_i = Vector{Matrix{T}}(undef, len_idxs_ci_i)

    for lp_ci in eachindex(idxs_ci)
        local_ci_i[lp_ci] = (z::AbstractArray) -> ci_i[idxs_ci[lp_ci]](z, pr, ps)[1]
        hessci_i[lp_ci] = ForwardDiff.hessian(local_ci_i[lp_ci], x)::Matrix{P}
    end
    return hessci_i
end


function mm_fd_dn_hessci_i!(out_hessci_i::AbstractArray, dimspec::MPCCDimSpec, ci_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    len_idxs_ci_i = length(idxs_ci)
    local_ci_i = Vector{Function}(undef, len_idxs_ci_i)

    for lp_ci in eachindex(idxs_ci)
        local_ci_i[lp_ci] = (z::AbstractArray) -> ci_i[idxs_ci[lp_ci]](z, pr, ps)[1]
        ForwardDiff.hessian!(@view(out_hessci_i[lp_ci][:]), local_ci_i[lp_ci], x)
    end
    return nothing
end








# Full hessF functions
# # Matrix of dim l x q of n length vectors (the hessians wrt x)

function mm_fd_dn_hessF(dimspec::MPCCDimSpec, F_i::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}   # ::Matrix{Matrix{S}}
    P = promote_type(S, T)
    @unpack n, l, q = dimspec

    local_F = Matrix{Function}(undef, l, q)       # TODO do we need to actually store all of them in an array?
    hessF = Matrix{Matrix{T}}(undef, l, q)

    for lp_q=1:q
        for lp_l=1:l
            local_F[lp_l, lp_q] = (z::AbstractArray) -> F_i[lp_l, lp_q](z, pr, ps)[1]
            hessF[lp_l, lp_q] = ForwardDiff.hessian(local_F[lp_l, lp_q], x)::Matrix{P}
        end
    end
    return hessF
end


function mm_fd_dn_hessF!(out_hessF::AbstractArray, dimspec::MPCCDimSpec,  F_i::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    @unpack n, l, q = dimspec

    local_F = Matrix{Function}(undef, l, q)

    for lp_q=1:q
        for lp_l=1:l
            local_F[lp_l, lp_q] = (z::AbstractArray) -> F_i[lp_l, lp_q](z, pr, ps)[1]
            ForwardDiff.hessian!(@view(out_hessF[lp_l, lp_q][:,:]), local_F[lp_l, lp_q], x)
        end
    end
    return nothing
end



# Indexed hessF functions
function mm_fd_dn_hessF_i(dimspec::MPCCDimSpec, F_i::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}}) where {F <: Function, S <: Real, T <: Real}
    P = promote_type(S, T)
    @unpack n, l, q = dimspec
    len_idxs_F_i = length(idxs_F)
    local_F_i = Vector{Function}(undef, len_idxs_F_i)
    hessF_i = Vector{Matrix{T}}(undef, len_idxs_F_i)
    for lp_F in eachindex(idxs_F)
        local_F_i[lp_F] = (z::AbstractArray) -> F_i[idxs_F[lp_F]](z, pr, ps)[1]
        hessF_i[lp_F] = ForwardDiff.hessian(local_F_i[lp_F], x)::Matrix{P}
    end
    return hessF_i
end


function mm_fd_dn_hessF_i!(out_hessF_i::AbstractArray, dimspec::MPCCDimSpec, F_i::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}})::Nothing where {F <: Function, S <: Real, T <: Real}
    @unpack n, l, q = dimspec
    len_idxs_F_i = length(idxs_F)
    local_F_i = Vector{Function}(undef, len_idxs_F_i)       # TODO do we need to actually store all of them in an array?
    for lp_F in eachindex(idxs_F)
        local_F_i[lp_F] = (z::AbstractArray) -> F_i[idxs_F[lp_F]](z, pr, ps)[1]
        ForwardDiff.hessian!(@view(out_hessF_i[lp_F][:,:]), local_F_i[lp_F], x)
    end
    return nothing
end





# Full hessFq functions

function mm_fd_dn_hessFq(dimspec::MPCCDimSpec, Fq_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}
    # Call the indexed version...
    idxs_Fq = collect(1:dimspec.q)
    return mm_fd_dn_hessFq_i(dimspec, Fq_i, x, pr, ps, idxs_Fq::AbstractVector{Int64})    
end


# Indexed hessFq functions

function mm_fd_dn_hessFq_i(dimspec::MPCCDimSpec, Fq_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_Fq::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real} 
    P = promote_type(S, T)
    @unpack n = dimspec

    len_idxs_Fq_i = length(idxs_Fq)
    local_Fq_i = Vector{Function}(undef, len_idxs_Fq_i)
    hessFq_i = Vector{Matrix{T}}(undef, len_idxs_Fq_i)

    for lp_Fq in eachindex(idxs_Fq)
        local_Fq_i[lp_Fq] = (z::AbstractArray) -> Fq_i[idxs_Fq[lp_Fq]](z, pr, ps)[1]
        hessFq_i[lp_Fq] = ForwardDiff.hessian(local_Fq_i[lp_Fq], x)::Matrix{P}
    end
    return hessFq_i
end








# fdp functions

function mm_fd_dn_fdp(dimspec::MPCCDimSpec, f::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}   # ::AbstractVector{T}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    local_f(qr) = f(x, qr, ps)[1]
    return ForwardDiff.gradient(local_f, pr)::Vector{P}
end


function mm_fd_dn_fdp!(out_fdp::AbstractArray, dimspec::MPCCDimSpec, f::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    local_f(z) = f(z, pr, ps)[1]
    ForwardDiff.gradient!(out_fdp, local_f, x)
    return nothing
end






# Full cedp functions
# # Vector of length r of vectors of length me
function mm_fd_dn_cedp(dimspec::MPCCDimSpec, ce::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}    
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack me, r = dimspec

    # Shortcut for empty result: ensures type stable (fails assert with FD), and gets dimensions correct
    if ( 0 == me )
        return Vector{P}[ Vector{P}(undef, 0) for lp_r=1:r ]
    end

    local_ce = (qr::AbstractArray) -> ce(x, qr, ps)
    temp_cedp = ForwardDiff.jacobian(local_ce, pr)::Matrix{P}
    cedp = [ [ temp_cedp[lp_me, lp_r] for lp_me in 1:me ] for lp_r in 1:r ]
    return cedp
end


function mm_fd_dn_cedp!(out_cedp::AbstractArray, dimspec::MPCCDimSpec, ce!::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack me, r = dimspec
    local_ce! = (y::AbstractArray, qr::AbstractArray) -> ce!(y, x, qr, ps)

    # ForwardDiff asks for a temporary storage area, presumably so that caller can manage memory allocs
    # We do this here for now, perhaps move it upstream later.
    y = zeros(promote_type(S, T), me)
    temp_cedp = zeros(promote_type(S, T), me, r)
    ForwardDiff.jacobian!(temp_cedp, local_ce!, y, pr)
    for lp_r=1:r
        for lp_me=1:me
            out_cedp[lp_r][lp_me] = temp_cedp[lp_me, lp_r]
        end
    end
    return nothing
end


# Indexed cedp functions
function mm_fd_dn_cedp_i(dimspec::MPCCDimSpec, ce_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack r = dimspec
    len_idxs_ce_i = length(idxs_ce)
    local_ce = (qr::AbstractArray) -> mm_fd_dn_ce_i(dimspec, ce_i, x, qr, ps, idxs_ce)
    temp_cedp = ForwardDiff.jacobian(local_ce, pr)::Matrix{P}
    cedp = [ [ temp_cedp[lp_me, lp_r] for lp_me in 1:len_idxs_ce_i ] for lp_r in 1:r ]
    return cedp
end

function mm_fd_dn_cedp_i!(out_cedp::AbstractArray, dimspec::MPCCDimSpec, ce_i!::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack me, r = dimspec
    len_idxs_ce_i = length(idxs_ce)
    local_ce! = (y::AbstractArray, qr::AbstractArray) -> mm_fd_dn_ce_i!(y, dimspec, ce_i!, x, qr, ps, idxs_ce)
    # mm_fd_dn_ce_i!(out_ce_i::AbstractVector{R}, dimspec::MPCCDimSpec, ce_i!::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64})::Nothing where {F <: Function, R <: Real, S <: Real, T <: Real}
    # ForwardDiff asks for a temporary storage area, presumably so that caller can manage memory allocs
    # We do this here for now, perhaps move it upstream later.
    y = zeros(promote_type(S, T), len_idxs_ce_i)
    temp_cedp = zeros(promote_type(S, T), len_idxs_ce_i, r)
    ForwardDiff.jacobian!(temp_cedp, local_ce!, y, pr)
    for lp_r=1:r
        for lp_me=1:len_idxs_ce_i
            out_cedp[lp_r][lp_me] = temp_cedp[lp_me, lp_r]
        end
    end
    return nothing
end





# Full cidp functions
# # Vector of length r of vectors of length mi
function mm_fd_dn_cidp(dimspec::MPCCDimSpec, ci::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}    
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack mi, r = dimspec

    # Shortcut for empty result: ensures type stable (fails assert with FD), and gets dimensions correct
    if ( 0 == mi )
        return Vector{P}[ Vector{P}(undef, 0) for lp_r=1:r ]
    end

    local_ci = (qr::AbstractArray) -> ci(x, qr, ps)
    temp_cidp = ForwardDiff.jacobian(local_ci, pr)::Matrix{P}
    cidp = [ [ temp_cidp[lp_mi, lp_r] for lp_mi in 1:mi ] for lp_r in 1:r ]
    return cidp
end


function mm_fd_dn_cidp!(out_cidp::AbstractArray, dimspec::MPCCDimSpec, ci!::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack mi, r = dimspec
    local_ci! = (y::AbstractArray, qr::AbstractArray) -> ci!(y, x, qr, ps)

    # ForwardDiff asks for a temporary storage area, presumably so that caller can manage mimory allocs
    # We do this here for now, perhaps move it upstream later.
    y = zeros(promote_type(S, T), mi)
    temp_cidp = zeros(promote_type(S, T), mi, r)
    ForwardDiff.jacobian!(temp_cidp, local_ci!, y, pr)
    for lp_r=1:r
        for lp_mi=1:mi
            out_cidp[lp_r][lp_mi] = temp_cidp[lp_mi, lp_r]
        end
    end
    return nothing
end


# Indexed cidp functions
function mm_fd_dn_cidp_i(dimspec::MPCCDimSpec, ci_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack r = dimspec
    len_idxs_ci_i = length(idxs_ci)
    local_ci = (qr::AbstractArray) -> mm_fd_dn_ci_i(dimspec, ci_i, x, qr, ps, idxs_ci)
    temp_cidp = ForwardDiff.jacobian(local_ci, pr)::Matrix{P}
    cidp = [ [ temp_cidp[lp_mi, lp_r] for lp_mi in 1:len_idxs_ci_i ] for lp_r in 1:r ]
    return cidp
end

function mm_fd_dn_cidp_i!(out_cidp::AbstractArray, dimspec::MPCCDimSpec, ci_i!::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack mi, r = dimspec
    len_idxs_ci_i = length(idxs_ci)
    local_ci! = (y::AbstractArray, qr::AbstractArray) -> mm_fd_dn_ci_i!(y, dimspec, ci_i!, x, qr, ps, idxs_ci)
    # mm_fd_dn_ci_i!(out_ci_i::AbstractVector{R}, dimspec::MPCCDimSpec, ci_i!::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64})::Nothing where {F <: Function, R <: Real, S <: Real, T <: Real}
    # ForwardDiff asks for a temporary storage area, presumably so that caller can manage mimory allocs
    # We do this here for now, perhaps move it upstream later.
    y = zeros(promote_type(S, T), len_idxs_ci_i)
    temp_cidp = zeros(promote_type(S, T), len_idxs_ci_i, r)
    ForwardDiff.jacobian!(temp_cidp, local_ci!, y, pr)
    for lp_r=1:r
        for lp_mi=1:len_idxs_ci_i
            out_cidp[lp_r][lp_mi] = temp_cidp[lp_mi, lp_r]
        end
    end
    return nothing
end







# Full Fdp functions

function mm_fd_dn_Fdp(dimspec::MPCCDimSpec, F::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack n, l, q, r = dimspec

    # Shortcut for empty result: ensures type stable (fails assert with FD), and gets dimensions correct
    if ( 0 == l || 0 == q )
        return Vector{P}[ Matrix{P}(undef, l, q) for lp_r=1:r ]
    end    

    # ForwardDiff.jacobian can work with multidimensional arrays, so we do that
    local_F = (qr::AbstractArray) -> F(x, qr, ps)
    Fdp = Vector{Matrix{T}}(undef, r)
    Fdp_flat = ForwardDiff.jacobian(local_F, pr)::Matrix{P}
    for lp_pr=1:r
        Fdp[lp_pr] = reshape(Fdp_flat[:, lp_pr], (l, q))
    end
    return Fdp
end



function mm_fd_dn_Fdp!(out_Fdp::AbstractArray, dimspec::MPCCDimSpec, F!::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack n, l, q, r = dimspec

    # ForwardDiff.jacobian can work with multidimensional arrays, so we do that
    local_F! = (z::AbstractArray, qr::AbstractArray) -> F!(z, x, qr, ps)

    # Alloc these locally for now, should really be done by caller to realise efficiency gains
    y = zeros(promote_type(S, T), l, q)
    Fdp_flat = zeros(promote_type(S, T), l*q, r)
    ForwardDiff.jacobian!(Fdp_flat, local_F!, y, pr)
    for lp_pr=1:r
        out_Fdp[lp_pr] = reshape(Fdp_flat[:, lp_pr], (l, q))
    end
    return nothing
end



# Indexed Fdp
function mm_fd_dn_Fdp_i(dimspec::MPCCDimSpec, F_i::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}}) where {F <: Function, S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack n, l, q, r = dimspec

    len_idxs_F_i = length(idxs_F)
    local_F = (qr::AbstractArray) -> mm_fd_dn_F_i(dimspec, F_i, x, qr, ps, idxs_F)
    Fdp = Vector{Vector{T}}(undef, r)
    Fdp_flat = ForwardDiff.jacobian(local_F, pr)::Matrix{P}
    for lp_pr=1:r
        Fdp[lp_pr] = Fdp_flat[(1:len_idxs_F_i)*lp_pr]
    end
    return Fdp
end


function mm_fd_dn_Fdp_i!(out_Fdp::AbstractArray, dimspec::MPCCDimSpec, F_i!::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}})::Nothing where {F <: Function, S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack n, l, q, r = dimspec

    len_idxs_F_i = length(idxs_F)
    local_F! = (z::AbstractArray, qr::AbstractArray) -> mm_fd_dn_F_i!(z, dimspec, F_i!, x, qr, ps, idxs_F)

    # Alloc these locally for now, should really be done by caller to realise efficiency gains
    y = zeros(promote_type(S, T), len_idxs_F_i)
    Fdp_flat = zeros(promote_type(S, T), len_idxs_F_i * r)
    ForwardDiff.jacobian!(Fdp_flat, local_F!, y, pr)
    for lp_pr=1:r
        for lp_F=1:len_idxs_F_i
            out_Fdp[lp_pr][lp_F] = Fdp_flat[lp_F+((lp_pr-1)*len_idxs_F_i)]
        end
    end
    return nothing
end






# Full gradfdp functions

function mm_fd_dn_gradfdp(dimspec::MPCCDimSpec, f::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack n, r = dimspec
    local_gradf = (qr::AbstractArray) -> mm_fd_dn_gradf(dimspec, f, x, qr, ps)
    gradfdp_flat = ForwardDiff.jacobian(local_gradf, pr)::Matrix{P}
    gradfdp = [ gradfdp_flat[:, lp_r] for lp_r in 1:r ]
    return gradfdp
end


function mm_fd_dn_gradfdp!(out_gradfdp::AbstractArray, dimspec::MPCCDimSpec, f::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack n, r = dimspec
    local_gradf! = (z::AbstractArray, qr::AbstractArray) -> mm_fd_dn_gradf!(z, dimspec, f, x, qr, ps)
    y = zeros(promote_type(S, T), n)
    gradfdp_flat = zeros(promote_type(S, T), n, r)
    ForwardDiff.jacobian!(gradfdp_flat, local_gradf!, y, pr)
    for lp_pr=1:r
        out_gradfdp[lp_pr] = gradfdp_flat[:, lp_pr]
    end
    return nothing
end






# Full gradcedp functions
function mm_fd_dn_gradcedp(dimspec::MPCCDimSpec, ce::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack n, me, r = dimspec

    # Shortcut for empty result: ensures type stable (fails assert with FD), and gets dimensions correct
    if ( 0 == me )
        return Matrix{P}[ Matrix{P}(undef, me, n) for lp_r=1:r ]
    end

    local_gradce = (qr::AbstractArray) -> mm_fd_dn_gradce(dimspec, ce, x, qr, ps)
    gradcedp_flat = ForwardDiff.jacobian(local_gradce, pr)::Matrix{P}
    gradcedp = [ reshape(gradcedp_flat[:, lp_pr], (me, n)) for lp_pr in 1:r ]
    return gradcedp
end


function mm_fd_dn_gradcedp!(out_gradcedp::AbstractArray, dimspec::MPCCDimSpec, ce!::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack n, me, r = dimspec
    local_gradce! = (z::AbstractArray, qr::AbstractArray) -> mm_fd_dn_gradce!(z, dimspec, ce!, x, qr, ps)

    # Local alloc, should perhaps do in caller
    y = zeros(promote_type(S, T), me, n)
    gradcedp_flat = zeros(promote_type(S, T), n*me, r)
    ForwardDiff.jacobian!(gradcedp_flat, local_gradce!, y, pr)
    for lp_pr=1:r
        out_gradcedp[lp_pr] = reshape(gradcedp_flat[:, lp_pr], (me, n))
    end
    return nothing
end


# Indexed gradcedp
function mm_fd_dn_gradcedp_i(dimspec::MPCCDimSpec, ce_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack n, me, r = dimspec
    len_idxs_ce_i = length(idxs_ce)
    local_gradce = (qr::AbstractArray) -> mm_fd_dn_gradce_i(dimspec, ce_i, x, qr, ps, idxs_ce)
    gradcedp_flat = ForwardDiff.jacobian(local_gradce, pr)::Matrix{P}
    gradcedp = [ reshape(gradcedp_flat[:, lp_pr], (len_idxs_ce_i, n)) for lp_pr in 1:r ]
    return gradcedp
end



function mm_fd_dn_gradcedp_i!(out_gradcedp::AbstractArray, dimspec::MPCCDimSpec, ce_i!::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ce::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack n, me, r = dimspec
    len_idxs_ce_i = length(idxs_ce)
    local_gradce! = (z::AbstractArray, qr::AbstractArray) -> mm_fd_dn_gradce_i!(z, dimspec, ce_i!, x, qr, ps, idxs_ce)

    # Local alloc, should perhaps do in caller
    y = zeros(promote_type(S, T), len_idxs_ce_i, n)
    gradcedp_flat = zeros(promote_type(S, T), n*me, r)
    ForwardDiff.jacobian!(gradcedp_flat, local_gradce!, y, pr)
    for lp_pr=1:r
        out_gradcedp[lp_pr] = reshape(gradcedp_flat[:, lp_pr], (len_idxs_ce_i, n))
    end
    return nothing
end







# Full gradcidp functions
function mm_fd_dn_gradcidp(dimspec::MPCCDimSpec, ci::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack n, mi, r = dimspec

    # Shortcut for empty result: ensures type stable (fails assert with FD), and gets dimensions correct
    if ( 0 == mi )
        return Matrix{P}[ Matrix{P}(undef, mi, n) for lp_r=1:r ]
    end

    local_gradci = (qr::AbstractArray) -> mm_fd_dn_gradci(dimspec, ci, x, qr, ps)
    gradcidp_flat = ForwardDiff.jacobian(local_gradci, pr)::Matrix{P}
    gradcidp = [ reshape(gradcidp_flat[:, lp_pr], (mi, n)) for lp_pr in 1:r ]
    return gradcidp
end


function mm_fd_dn_gradcidp!(out_gradcidp::AbstractArray, dimspec::MPCCDimSpec, ci!::Function, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64})::Nothing where {S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack n, mi, r = dimspec
    local_gradci! = (z::AbstractArray, qr::AbstractArray) -> mm_fd_dn_gradci!(z, dimspec, ci!, x, qr, ps)

    # Local alloc, should perhaps do in caller
    y = zeros(promote_type(S, T), mi, n)
    gradcidp_flat = zeros(promote_type(S, T), n*mi, r)
    ForwardDiff.jacobian!(gradcidp_flat, local_gradci!, y, pr)
    for lp_pr=1:r
        out_gradcidp[lp_pr] = reshape(gradcidp_flat[:, lp_pr], (mi, n))
    end
    return nothing
end


# Indexed gradcidp
function mm_fd_dn_gradcidp_i(dimspec::MPCCDimSpec, ci_i::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack n, mi, r = dimspec
    len_idxs_ci_i = length(idxs_ci)
    local_gradci = (qr::AbstractArray) -> mm_fd_dn_gradci_i(dimspec, ci_i, x, qr, ps, idxs_ci)
    gradcidp_flat = ForwardDiff.jacobian(local_gradci, pr)::Matrix{P}
    gradcidp = [ reshape(gradcidp_flat[:, lp_pr], (len_idxs_ci_i, n)) for lp_pr in 1:r ]
    return gradcidp
end


function mm_fd_dn_gradcidp_i!(out_gradcidp::AbstractArray, dimspec::MPCCDimSpec, ci_i!::AbstractVector{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_ci::AbstractVector{Int64})::Nothing where {F <: Function, S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack n, mi, r = dimspec
    len_idxs_ci_i = length(idxs_ci)
    local_gradci! = (z::AbstractArray, qr::AbstractArray) -> mm_fd_dn_gradci_i!(z, dimspec, ci_i!, x, qr, ps, idxs_ci)

    # Local alloc, should perhaps do in caller
    y = zeros(promote_type(S, T), len_idxs_ci_i, n)
    gradcidp_flat = zeros(promote_type(S, T), n*mi, r)
    ForwardDiff.jacobian!(gradcidp_flat, local_gradci!, y, pr)
    for lp_pr=1:r
        out_gradcidp[lp_pr] = reshape(gradcidp_flat[:, lp_pr], (len_idxs_ci_i, n))
    end
    return nothing
end







# Full gradFdp functions
function mm_fd_dn_gradFdp(dimspec::MPCCDimSpec, F_i::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}   # ::Vector{Matrix{Vector{S}}}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack n, l, q, r = dimspec
    local_gradF_i = Matrix{Function}(undef, l, q)
    gradFdp = [ [ zeros(promote_type(S, T), n) for lp_l in 1:l, lp_q in 1:q ] for lp_r in 1:r ]
    for lp_q=1:q
        for lp_l=1:l
            local_gradF_i[lp_l, lp_q] = (qr::AbstractArray) -> mm_fd_dn_gradF_i(dimspec, F_i, x, qr, ps, [CartesianIndex(lp_l, lp_q)])[1]
            temp_graddp_flat_i = ForwardDiff.jacobian(local_gradF_i[lp_l, lp_q], pr)::Matrix{P}
            for lp_r=1:r
                gradFdp[lp_r][lp_l, lp_q][:] = temp_graddp_flat_i[:, lp_r]
            end
        end
    end
    return gradFdp
end



function mm_fd_dn_gradFdp!(out_gradFdp::AbstractArray, dimspec::MPCCDimSpec, F_i::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}) where {F <: Function, S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack n, l, q, r = dimspec
    local_gradF_i = Matrix{Function}(undef, l, q)
    # Local alloc, should perhaps do in caller
    y = zeros(promote_type(S, T), n)
    gradFdp_flat_i = zeros(promote_type(S, T), n, r)
    for lp_q=1:q
        for lp_l=1:l
            # Can't get it working fully mutating because it wants to mutate into a vector of vectors,
            # and forwarddiff needs just a plain vector.  So doesn't work.

            # local_gradF_i![lp_l, lp_q] = function (z::AbstractArray, qr::AbstractArray)
            #     vec_of_z = [ z ]
            #     mm_fd_dn_gradF_i!(vec_of_z, dimspec, F_i!, x, qr, ps, [CartesianIndex(lp_l, lp_q)])
            #     println("vec of z: ", vec_of_z)
            # end
            # ForwardDiff.jacobian!(gradFdp_flat_i, local_gradF_i![lp_l, lp_q], y, pr)
            local_gradF_i[lp_l, lp_q] = (qr::AbstractArray) -> mm_fd_dn_gradF_i(dimspec, F_i, x, qr, ps, [CartesianIndex(lp_l, lp_q)])[1]
            ForwardDiff.jacobian!(gradFdp_flat_i, local_gradF_i[lp_l, lp_q], pr)
            for lp_r=1:r
                out_gradFdp[lp_r][lp_l, lp_q][:] = gradFdp_flat_i[:, lp_r]
            end
        end
    end
    return nothing
end



# Indexed gradFdp functions
function mm_fd_dn_gradFdp_i(dimspec::MPCCDimSpec, F_i::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}}) where {F <: Function, S <: Real, T <: Real}   # ::Vector{Matrix{Vector{S}}}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    P = promote_type(S, T)
    @unpack n, l, q, r = dimspec
    len_idxs_F_i = length(idxs_F)
    local_gradF_i = Vector{Function}(undef, len_idxs_F_i)

    # This time we need to do some construction work...
    gradFdp = [ [ zeros(promote_type(S, T), n) for lp_idx in 1:len_idxs_F_i ] for lp_r in 1:r ]

    # Can probably do this in one go rather than loop TODO
    for lp_idx=1:len_idxs_F_i
        local_gradF_i[lp_idx] = (qr::AbstractArray) -> mm_fd_dn_gradF_i(dimspec, F_i, x, qr, ps, [ idxs_F[lp_idx] ])[1]
        temp_graddp_flat_i = ForwardDiff.jacobian(local_gradF_i[lp_idx], pr)::Matrix{P}
        for lp_r=1:r
            gradFdp[lp_r][lp_idx][:] = temp_graddp_flat_i[:, lp_r]
        end
    end
    return gradFdp
end


function mm_fd_dn_gradFdp_i!(out_gradFdp::AbstractArray, dimspec::MPCCDimSpec, F_i::AbstractMatrix{F}, x::AbstractVector{S}, pr::AbstractVector{T}, ps::AbstractVector{Int64}, idxs_F::Vector{CartesianIndex{2}}) where {F <: Function, S <: Real, T <: Real}
    @assert dimspec.r > 0 "r must be positive for parametric calls"
    @unpack n, l, q, r = dimspec
    len_idxs_F_i = length(idxs_F)
    local_gradF_i = Vector{Function}(undef, len_idxs_F_i)
    gradFdp_flat_i = zeros(promote_type(S, T), n, r)

    # Can probably do this in one go rather than loop TODO
    for lp_idx=1:len_idxs_F_i
        local_gradF_i[lp_idx] = (qr::AbstractArray) -> mm_fd_dn_gradF_i(dimspec, F_i, x, qr, ps, [ idxs_F[lp_idx] ])[1]
        ForwardDiff.jacobian!(gradFdp_flat_i, local_gradF_i[lp_idx], pr)
        for lp_r=1:r
            out_gradFdp[lp_r][lp_idx][:] = gradFdp_flat_i[:, lp_r]
        end

    end
    return nothing
end
