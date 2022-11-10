





#--------------------------------------------------
# Convert between F Cartesian and linear indexing
#--------------------------------------------------


"""
    F_cart_to_lin(dimspec::MPCCDimSpec, cart_idx::CartesianIndex{2})

Systematic conversion of 2-tuple CartesianIndex to a linear index. Only depends
on problem dimesions not a specific constraint set.
"""
function F_cart_to_lin(dimspec::MPCCDimSpec, cart_idx::CartesianIndex{2})
    @unpack n, q, l, me, mi = dimspec
    @assert cart_idx[1] <= l && cart_idx[2] <= q
    return (cart_idx[2] - 1) * l + cart_idx[1]
end



"""
    F_lin_to_cart(dimspec::MPCCDimSpec, lin_idx::Int64)

Systematic conversion of linear index to 2-tuple CartesianIndex
"""
function F_lin_to_cart(dimspec::MPCCDimSpec, lin_idx::Int64)
    @unpack n, q, l, me, mi = dimspec
    il = mod(lin_idx - 1, l) + 1
    iq = div(lin_idx - 1, l) + 1
    @assert il <= l && iq <= q
    return CartesianIndex(il, iq)
end





#--------------------------------------------------
# Convert between base index and cpair index
#--------------------------------------------------

"""
    cs_cpair_to_bi(dimspec::MPCCDimSpec, ctype::MPCCActivesetCnstrType, el::MPCCIdxElement)

Convert from a 'cpair' of constraint type (CE, CI, F) and element index, to a
base index in the defined ordering. Only depends on the problem dimensions, not
an actual constraint set.
"""
function cs_cpair_to_bi(dimspec::MPCCDimSpec, ctype::MPCCActivesetCnstrType, el::MPCCIdxElement)
    @unpack n, q, l, me, mi = dimspec

    if (CNSTR_TYPE_CE == ctype)
        @assert typeof(el) <: Integer
        @assert el >= 0 && el <= me
        return el
    elseif (CNSTR_TYPE_CI == ctype)
        @assert typeof(el) <: Integer
        @assert el >= 0 && el <= mi
        return el + me
    elseif (CNSTR_TYPE_F == ctype)
        if (typeof(el) <: CartesianIndex)
            @assert el[1] >= 1 && el[1] <= l
            @assert el[2] >= 1 && el[2] <= q
            return me + mi + F_cart_to_lin(dimspec, el)
        elseif (typeof(el) <: Integer)
            @assert el <= l * q
            return me + mi + el
        end
    else
        # Informative.
        error("Unknown constraint type or other error.")
    end
end



"""
    cs_bi_to_cpair(dimspec::MPCCDimSpec, el::Int64; b_F_cart=true)

Convert a base index, in the defined ordering, to a cpair of constraint type
(CE, CI, F) and element index. Only depends on problem dimensions, not on an
actual constraint set.
"""
function cs_bi_to_cpair(dimspec::MPCCDimSpec, el::Int64; b_F_cart=true)
    @unpack n, q, l, me, mi = dimspec

    @assert el >= 1 && el <= (me + mi + l * q)
    if (el <= me)
        return (CNSTR_TYPE_CE, el)
    elseif (el > me && el <= (me + mi))
        return (CNSTR_TYPE_CI, el - me)
    elseif (el > (me + mi))
        el_F = el - (me + mi)
        if b_F_cart
            return (CNSTR_TYPE_F, F_lin_to_cart(dimspec, el_F))
        else
            return (CNSTR_TYPE_F, el_F)
        end
    end
end







#--------------------------------------------------
# Convert between constraint set index and BI
#--------------------------------------------------


"""
    cs_cnstr_el_to_bi(cnstr_idxs::MPCCConstraintSet, el::Int64)

Convert element of a constraint set to the base index. This DOES depend on the
constraint set.
"""
function cs_cnstr_el_to_bi(cnstr_idxs::MPCCConstraintSet, el::Int64)
    return cnstr_idxs.to_bi[el]
end



"""
    cs_bi_to_cnstr_el(cnstr_idxs::MPCCConstraintSet, bi_idx::Int64)

Convert a base index to the constraint set index.  This may return missing since
the specified BI might not be in the constraint set (not active).  This DOES
depend on the constraint set.
"""
function cs_bi_to_cnstr_el(cnstr_idxs::MPCCConstraintSet, bi_idx::Int64)
    search_uidx = searchsorted(cnstr_idxs.to_bi, bi_idx)
    if (length(search_uidx) > 0)
        return search_uidx.start
    else
        return missing
    end
end




#--------------------------------------------------
# Convert to and from cpair and constraint set element
#--------------------------------------------------

"""
    cs_cnstr_el_to_cpair(cnstr_idxs::MPCCConstraintSet, el::Int64; b_F_cart=true)

Convert element el of constraint set to the cpair format.
"""
function cs_cnstr_el_to_cpair(cnstr_idxs::MPCCConstraintSet, el::Int64; b_F_cart=true)
    return cs_bi_to_cpair(cnstr_idxs.config.dimspec, cnstr_idxs.to_bi[el]; b_F_cart)
end



"""
    cs_cpair_to_cnstr_el(cnstr_idxs::MPCCConstraintSet, ctype::MPCCActivesetCnstrType, el::MPCCIdxElement)

Convert cpair element of a constraint set to its element in the constraint set
(at which place is it in the bi index)
"""
function cs_cpair_to_cnstr_el(cnstr_idxs::MPCCConstraintSet, ctype::MPCCActivesetCnstrType, el::MPCCIdxElement)
    bi_idx = cs_cpair_to_bi(cnstr_idxs.config.dimspec, ctype, el)
    return cs_bi_to_cnstr_el(cnstr_idxs, bi_idx)
end







#--------------------------------------------------
# Manipulation of a constraint set + bitmask stuff
#--------------------------------------------------


"""
    cs_cnstr_idxs_add_by_bi!(cnstr_idxs::MPCCConstraintSet, bi_idx::Int64)

Add element into constraint set by its base index.
"""
function cs_cnstr_idxs_add_by_bi!(cnstr_idxs::MPCCConstraintSet, bi_idx::Int64)
    # Check element doesn't already exist
    search_uidx = searchsorted(cnstr_idxs.to_bi, bi_idx)
    @assert(0 == length(search_uidx), "Cannot add: element already exists in the index set!")

    insert!(cnstr_idxs.to_bi, search_uidx.start, bi_idx)
end



"""
    cs_cnstr_idxs_add_by_cpair!(cnstr_idxs::MPCCConstraintSet, ctype::MPCCActivesetCnstrType, el::MPCCIdxElement)

Add element into the constraint set by its cpair.
"""
function cs_cnstr_idxs_add_by_cpair!(cnstr_idxs::MPCCConstraintSet, ctype::MPCCActivesetCnstrType, el::MPCCIdxElement)
    bi_idx = cs_cpair_to_bi(cnstr_idxs.config.dimspec, ctype, el)
    cs_cnstr_idxs_add_by_bi!(cnstr_idxs, bi_idx)
end



"""
    cs_cnstr_idxs_del_by_cnstr_el!(cnstr_idxs::MPCCConstraintSet, el::Int64)

Delete element of constraint set by its element number.
"""
function cs_cnstr_idxs_del_by_cnstr_el!(cnstr_idxs::MPCCConstraintSet, el::Int64)
    deleteat!(cnstr_idxs.to_bi, el)
end



"""
    cs_cnstr_idxs_del_by_bi!(cnstr_idxs::MPCCConstraintSet, bi_idx::Int64)

Delete element of a constraint set by its base index.
"""
function cs_cnstr_idxs_del_by_bi!(cnstr_idxs::MPCCConstraintSet, bi_idx::Int64)
    el = cs_bi_to_cnstr_el(cnstr_idxs, bi_idx)
    @assert(el !== missing, "Cannot delete: there is no element with that bi index!")
    deleteat!(cnstr_idxs.to_bi, el)
end



"""
    cs_cnstr_idxs_del_by_cpair!(cnstr_idxs::MPCCConstraintSet, ctype::MPCCActivesetCnstrType, el::MPCCIdxElement)

Delete constraint set element by its cpair.
"""
function cs_cnstr_idxs_del_by_cpair!(cnstr_idxs::MPCCConstraintSet, ctype::MPCCActivesetCnstrType, el::MPCCIdxElement)
    el = cs_cpair_to_cnstr_el(cnstr_idxs, ctype, el)
    deleteat!(cnstr_idxs.to_bi, el)
end



"""
    cs_cnstr_get_as_bitmask(cnstr_idxs::MPCCConstraintSet)

Get the constraint set expressed as three bitmasks: BitVector for CE and CI;
BitMatrix for F. Returned as NamedTuple.
"""
function cs_cnstr_get_as_bitmask(cnstr_idxs::MPCCConstraintSet)
    @unpack n, q, l, me, mi = cnstr_idxs.config.dimspec

    # Initialise with zero vectors or matrices
    ce_mask = falses(me)
    ci_mask = falses(mi)
    F_mask = falses(l, q)

    # Iterate to fill in any true values in CE
    for lp_i=1:length(cnstr_idxs)
        (cur_cnstr_type, cur_cnstr_el) = cs_bi_to_cpair(cnstr_idxs.config.dimspec, lp_i)

        if ( cur_cnstr_type == CNSTR_TYPE_CE )
            ce_mask[cur_cnstr_el] = true
        elseif ( cur_cnstr_type == CNSTR_TYPE_CI )
            ci_mask[cur_cnstr_el] = true
        elseif ( cur_cnstr_type == CNSTR_TYPE_F )
            F_mask[cur_cnstr_el] = true
        else
            @error "not a valid constraint type"
        end
    end

    return NamedTuple{(:ce_mask, :ci_mask, :F_mask)}((ce_mask, ci_mask, F_mask))
end



"""
    cs_cnstr_build_from_bitmask(dimspec::MPCCDimSpec, ce_actv::BitVector, ci_actv::BitVector, F_actv::BitMatrix)

Build whole constraint set from bitmasks of the contraints for each type.
"""
function cs_cnstr_build_from_bitmask(config::MPCCModelConfig, ce_actv::BitVector, ci_actv::BitVector, F_actv::BitMatrix)
    @unpack n, q, l, me, mi = config.dimspec

    # TODO we can do this faster if done manually, but can't be bothered at the moment
    cnstr_idxs = MPCCConstraintSet(config, Int64[])

    # ce
    for lp_ce = 1:me
        # Altered 20220814 to allow selection of CE constraints, as might actually need to do this
        if (ce_actv[lp_ce])
            cs_cnstr_idxs_add_by_cpair!(cnstr_idxs, CNSTR_TYPE_CE, lp_ce)
        end
    end

    # ci
    for lp_ci = 1:mi
        if (ci_actv[lp_ci])
            cs_cnstr_idxs_add_by_cpair!(cnstr_idxs, CNSTR_TYPE_CI, lp_ci)
        end
    end

    # F
    for lp_q = 1:q
        for lp_l = 1:l
            if (F_actv[lp_l, lp_q])
                cs_cnstr_idxs_add_by_cpair!(cnstr_idxs, CNSTR_TYPE_F, CartesianIndex(lp_l, lp_q))
            end
        end
    end

    return cnstr_idxs
end



"""
    cs_cnstr_build_all_active(dimspec::MPCCDimSpec)

Build a constraint set with all constraints included (active).
"""
function cs_cnstr_build_all_active(config::MPCCModelConfig)
    @unpack n, q, l, me, mi = config.dimspec

    # TODO we can do this faster is done manually, but can't be bothered at the moment
    cnstr_idxs = MPCCConstraintSet(config, Int64[])

    # ce
    for lp_ce = 1:me
        cs_cnstr_idxs_add_by_cpair!(cnstr_idxs, CNSTR_TYPE_CE, lp_ce)
    end

    # ci
    for lp_ci = 1:mi
        cs_cnstr_idxs_add_by_cpair!(cnstr_idxs, CNSTR_TYPE_CI, lp_ci)
    end

    # F
    for lp_q = 1:q
        for lp_l = 1:l
            cs_cnstr_idxs_add_by_cpair!(cnstr_idxs, CNSTR_TYPE_F, CartesianIndex(lp_l, lp_q))
        end
    end

    return cnstr_idxs
end




#--------------------------------------------------
# Fq (columnwise complementarity) type functions
#--------------------------------------------------


"""
    cs_cnstr_get_Fq_count(cnstr_idxs::MPCCConstraintSet)

Return the number of active constraints in each column q of the F constraints.
"""
function cs_cnstr_get_Fq_count(cnstr_idxs::MPCCConstraintSet)
    @unpack n, q, l, me, mi = cnstr_idxs.config.dimspec

    Fq_cnt = zeros(Int64, q)
    for lp_q = 1:q
        actv_ctr = 0
        for lp_l = 1:l
            cur_idx = cs_cpair_to_cnstr_el(cnstr_idxs, CNSTR_TYPE_F, CartesianIndex(lp_l, lp_q))
            if (cur_idx !== missing)
                actv_ctr += 1
            end
        end
        Fq_cnt[lp_q] = actv_ctr
    end

    return Fq_cnt
end



"""
    cs_cnstr_get_all_active_bi(cnstr_idxs::MPCCConstraintSet)

Returns BI entries for all contraints in `cnstr_idxs`, separated into type in a named tuple.
"""
function cs_cnstr_get_all_active_bi(cnstr_idxs::MPCCConstraintSet)
    @unpack n, q, l, me, mi = cnstr_idxs.config.dimspec

    ce_actv = Vector{Int64}([])
    ci_actv = Vector{Int64}([])
    F_actv = Vector{CartesianIndex{2}}([])
    for lp_i = 1:length(cnstr_idxs)
        cur_cpair = cs_bi_to_cpair(dimspec, cnstr_idxs.to_bi[lp_i])
        if (cur_cpair[1] == CNSTR_TYPE_CE)
            push!(ce_actv, cur_cpair[2])
        elseif (cur_cpair[1] == CNSTR_TYPE_CI)
            push!(ci_actv, cur_cpair[2])
        elseif (cur_cpair[1] == CNSTR_TYPE_F)
            push!(F_actv, cur_cpair[2])
        else
            error("This shouldn't happen: not a valid constriant type.")
        end
    end

    return NamedTuple{(:ce_actv, :ci_actv, :F_actv)}((ce_actv, ci_actv, F_actv))
end



