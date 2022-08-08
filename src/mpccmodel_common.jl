


# TODO 20220708: maybe change function names gradxx to jacxx

# TODO Change Vector{Num} to either Symbolics.Arr{Num, 1} or something more general. May require some finessing


export  MPCCDimSpec,
        MPCCDefinition,
        MPCCFunctions,
        MPCCPointEvalReq,
        MPCCPointEval,
        MPCCModelTestVector,
        MPCCParameterisationDefn,
        MPCCParameterisationFunctions,
        MPCCParameterisations,
        MPCCModelConfig,
        AbstractMPCCModel,
        AbstractMPCCModelDense,
        AbstractMPCCModelSparse,
        MPCCModelDenseForwardDiff,
        MPCCModelSparseSym,
        MPCCJumpFixedFunctions,
        MPCCNewtonPenaltyFunctions

# Trick from https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/15
const Opt{T} = Union{Missing,T}


# Get type parameters from instance of SparseMetrixCSC
sp_getparamtypes(mtrx::SparseMatrixCSC{TT, SS}) where {TT, SS} = NamedTuple{(:T, :S)}((TT, SS))



"""
MPCCDimSpec

Determines the dimensions of a program.

* n: number of spatial dimensions
* q: columns in F (number of complementarity contraints)
* l: row in F (number of expressions in each complementarity contraint)
* me: numebr of equality constraints
* me: numebr of inequality constraints
* r: number of continuous parameters
* s: number of integer parameters

"""
struct MPCCDimSpec
    n::Int64		# Dimension of spatial variables
    q::Int64		# Columns in F
    l::Int64        # Number of F functions per complementarity variable (usually 2)
    me::Int64		# Number of equality constraints
    mi::Int64		# Number of inequality constraints
	r::Int64		# Number of continuous parameters for PMPCC
	s::Int64		# Number of discrete parameters for PMPCC
end
function MPCCDimSpec(dict::Dict)
    return MPCCDimSpec(
        dict["n"],
        dict["q"],
        dict["l"],
        dict["me"],
        dict["mi"],
        dict["r"],
        dict["s"]
    )
end


function show(io::IO, ::MIME"text/plain", dimspec::MPCCDimSpec)
    # print("DimSpec[n=$(dimspec.n), me=$(dimspec.me), me=$(dimspec.mi), q=$(dimspec.q), l=$(dimspec.l), r=$(dimspec.r), s=$(dimspec.s)]")
    printstyled("DimSpec: n="; color=:magenta)
    printstyled(dimspec.n; color=:light_yellow)
    printstyled(", me="; color=:magenta)
    printstyled(dimspec.me; color=:light_yellow)
    printstyled(", mi="; color=:magenta)
    printstyled(dimspec.mi; color=:light_yellow)
    printstyled(", q="; color=:magenta)
    printstyled(dimspec.q; color=:light_yellow)
    printstyled(", r="; color=:magenta)
    printstyled(dimspec.r; color=:light_yellow)
    printstyled(", s="; color=:magenta)
    printstyled(dimspec.s, "\n"; color=:light_yellow)
end





"""
struct MPCCDefinition

Definition of the MPCC/NLP containing Symbolics Num expressions for f, ce, ci,
and F.

* f: objective function
* ce: vector of equality constraints ( == 0 )
* ci: vector of inequality constraints ( >= 0 )
* F: matrix l-by-q of complementarity constraints ( >= 0 )

"""
struct MPCCDefinition
    f::Num
    ce::Vector{Num}
    ci::Vector{Num}
    F::Matrix{Num}
    label_model::Opt{String}
    label_f::Opt{String}
    labels_ce::Vector{Opt{String}}
    labels_ci::Vector{Opt{String}}
    labels_F::Vector{Opt{String}}
end


function show(io::IO, ::MIME"text/plain", defn::MPCCDefinition)
    # Can't bring in an associated dimspec here, so have to count things ourselves.
    me = length(defn.ce)
    mi = length(defn.ci)
    (l, q) = size(defn.F)

    printstyled("---[ MPCCDefinition ]---\n"; color=:magenta)
    printstyled("Label:\t"; color=:magenta)
    printstyled(defn.label_model, "\n"; color=:light_cyan)
    printstyled("Objective function: "; color=:magenta)
    printstyled(defn.f, ";\t"; color=:light_yellow)
    printstyled(defn.label_f, "\n"; color=:light_cyan)
    for lp_me=1:me
        printstyled("ce[$lp_me]: "; color=:magenta)
        printstyled(defn.ce[lp_me], ";\t"; color=:light_yellow)
        printstyled(defn.labels_ce[lp_me], "\n"; color=:light_cyan)
    end
    for lp_mi=1:mi
        printstyled("ci[$lp_mi]: "; color=:magenta)
        printstyled(defn.ci[lp_mi], ";\t"; color=:light_yellow)
        printstyled(defn.labels_ci[lp_mi], "\n"; color=:light_cyan)
    end
    for lp_q=1:q
        printstyled("F[:, $lp_q]: "; color=:magenta)
        printstyled(defn.labels_F[lp_q], "\n"; color=:light_cyan)
        for lp_l=1:l
            printstyled("F[$lp_l, $lp_q]: "; color=:magenta)
            printstyled(defn.F[lp_l, lp_q], "\n"; color=:light_yellow)
        end
    end
end




"""
MPCCFunctions

A MPCCDefinition compiled into functions. This contains both standard and
mutating functions, and individual elements encoded as a function within a
vector or matrix.
"""
struct MPCCFunctions
    f::Function
    f!::Function
    ce::Function                   # Returns vector of length me
    ce!::Function
    ce_i::AbstractVector           # Vector so we can write scalar values in to mutating function argument
    ce_i!::AbstractVector
    ci::Function                   # Returns vector of length mi
    ci!::Function
    ci_i::AbstractVector
    ci_i!::AbstractVector
    F::Function                    # Returns vatrix of dim l x q
    F!::Function
    F_i::AbstractMatrix
    F_i!::AbstractMatrix
    Fq::Function                    # Returns vector of length q
    Fq!::Function
    Fq_i::AbstractVector
    Fq_i!::AbstractVector
end



"""
MPCCPointEvalReq

When making a request to evaluate several things in one request, this bitmask
determines what will be evaluated.
"""
struct MPCCPointEvalReq
    f::Bool
    ce::Bool
    ci::Bool
    F::Bool
    gradf::Bool
    gradce::Bool
    gradci::Bool
    gradF::Bool
    hessf::Bool
    hessce::Bool
    hessci::Bool
    hessF::Bool
    fdp::Bool
    cedp::Bool
    cidp::Bool
    Fdp::Bool
    gradfdp::Bool
    gradcedp::Bool
    gradcidp::Bool
    gradFdp::Bool                           
end
function MPCCPointEvalReq()
    return MPCCPointEvalReq(    true, true, true, true,
                                true, true, true, true,
                                true, true, true, true,
                                true, true, true, true,
                                true, true, true, true
                            )
end




"""
MPCCPointEval

Contains results of a point evaluation.
"""
struct MPCCPointEval{T <: AbstractFloat}
    f::Opt{T}
    ce::Opt{Vector{T}}                          # Vector of length me
    ci::Opt{Vector{T}}                          # Vector of length mi
    F::Opt{Matrix{T}}                           # Matrix of dim l x q
    gradf::Opt{Vector{T}}                       # Vector of length n
    gradce::Opt{Matrix{T}}                      # Changed 20211029: now matrix of me by n
    gradci::Opt{Matrix{T}}                      # Changed 20211029: now matrix of mi by n
    gradF::Opt{Matrix{Vector{T}}}               # Matrix of dim l x q of n length vectors (the gradients wrt x)
    hessf::Opt{Matrix{T}}                       # Matrix dim n x n
    hessce::Opt{Vector{Matrix{T}}}              # Vector of length me of 2d Hessian matrices
    hessci::Opt{Vector{Matrix{T}}}              # Vector of length mi of 2d Hessian matrices
    hessF::Opt{Matrix{Matrix{T}}}               # Matrix of dim l x q of 2d Hessian matrices
    fdp::Opt{Vector{T}}                         # Vector of length r
    cedp::Opt{Vector{Vector{T}}}                # Vector of length r of vectors of length me
    cidp::Opt{Vector{Vector{T}}}                # Vector of length r of vectors of length mi
    Fdp::Opt{Vector{Matrix{T}}}                 # Vector of length r of matrix of dim l x q
    gradfdp::Opt{Vector{Vector{T}}}             # Vector of length r of gradients of f/dp (vectors of length n)
    gradcedp::Opt{Vector{Matrix{T}}}            # Vector of length r of gradce, each /dp
    gradcidp::Opt{Vector{Matrix{T}}}            # Vector of length r of gradci, each /dp
    gradFdp::Opt{Vector{Matrix{Vector{T}}}}     # Vector of length r of gradF, each /dp
end
function MPCCPointEval()  # Quick way to create empty/zero struct
	return MPCCPointEval(
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing,
        missing
    )
end
Base.:(==)(x::MPCCPointEval, y::MPCCPointEval) = (
	x.f == y.f &&
    x.ce == y.ce &&
    x.ci == y.ci &&
	x.F == y.F &&
	x.gradf == y.gradf &&
    x.gradce == y.gradce &&
    x.gradci == y.gradci &&
    x.gradF == y.gradF &&
    x.hessf == y.hessf &&
    x.hessce == y.hessce &&
    x.hessci == y.hessci &&
    x.hessF == y.hessF &&
    x.fdp == y.fdp &&
    x.cedp == y.cedp &&
    x.cidp == y.cidp &&
    x.Fdp == y.Fdp &&
    x.gradfdp == y.gradfdp &&
    x.gradcedp == y.gradcedp &&
    x.gradcidp == y.gradcidp &&
    x.gradFdp == y.gradFdp )
Base.:(≈)(x::MPCCPointEval, y::MPCCPointEval) = (
    x.f ≈ y.f &&
    x.ce ≈ y.ce &&
    x.ci ≈ y.ci &&
    x.F ≈ y.F &&
    x.gradf ≈ y.gradf &&
    x.gradce ≈ y.gradce &&
    x.gradci ≈ y.gradci &&
    x.gradF ≈ y.gradF &&
    x.hessf ≈ y.hessf &&
    x.hessce ≈ y.hessce &&
    x.hessci ≈ y.hessci &&
    x.hessF ≈ y.hessF &&
    x.fdp ≈ y.fdp &&
    x.cedp ≈ y.cedp &&
    x.cidp ≈ y.cidp &&
    x.Fdp ≈ y.Fdp &&
    x.gradfdp ≈ y.gradfdp &&
    x.gradcedp ≈ y.gradcedp &&
    x.gradcidp ≈ y.gradcidp &&
    x.gradFdp ≈ y.gradFdp )






struct MPCCModelTestVector{R <: Real, S <: Real, T <: Real}
    x_val::Vector{S}
    pr_val::Vector{T}
    ps_val::Vector{Int64}
    eval_test_pt::MPCCPointEval{R}
end


"""
MPCCParameterisationDefn

The definition for a parameterisation of a model, i.e. how the `pr` depends on `t`.
"""
struct MPCCParameterisationDefn{R <: Real}
    pr::Vector{Num}
    # prdt::Num
    tspan::Tuple{R, R}
    descr::String
end


function show(io::IO, ::MIME"text/plain", pdefn::MPCCParameterisationDefn{R}) where {R <: Real}
    printstyled("---[ MPCCParameterisationDefn ]---\n"; color=:magenta)
    printstyled("Label: "; color=:magenta)
    printstyled(pdefn.descr, "\n"; color=:light_cyan)
    printstyled("tspan: "; color=:magenta)
    printstyled(pdefn.tspan, "\n"; color=:light_yellow)
    printstyled("pr: "; color=:magenta)
    printstyled(pdefn.pr, "\n"; color=:light_yellow)
end



"""
MPCCParameterisationFunctions

Compilation of `MPCCParameterisationFunctions` into Julia functions, including derivative.
"""
struct MPCCParameterisationFunctions
    pr::Function            # t -> Vector pr
    pr!::Function
    prdt::Function          # t -> Vector of d(pr) / dt
    prdt!::Function
end


"""
MPCCParameterisations

Collection of Num `t` along with definitions of parameterisations and their
compilation.
"""
struct MPCCParameterisations
    t::Num  # Symbolics variable
    defns::Vector{MPCCParameterisationDefn}
    fns::Vector{MPCCParameterisationFunctions}
end



"""
MPCCModelConfig

Model config contains the Num variables for `x`, `pr`, `ps` along with the
`MPCCDimSpec`, definition, compiled definition, any test vectors, any known
solutions, and parameterisations. Basically everything that specifies a model.

- parameterisations: if using a static model with r=0, then this should be set to missing
"""
struct MPCCModelConfig
    x::Vector{Num}      # Symbolics variables
    pr::Vector{Num}     # Symbolics variables
    ps::Vector{Num}     # Symbolics variables
    dimspec::MPCCDimSpec
    defn::MPCCDefinition
    fns::MPCCFunctions
    parameterisations::Opt{MPCCParameterisations}
    testvectors::Vector{MPCCModelTestVector}
    knownsols::Vector{Function}
end





abstract type AbstractMPCCModel end
abstract type AbstractMPCCModelDense <: AbstractMPCCModel end
abstract type AbstractMPCCModelSparse <: AbstractMPCCModel end


"""
MPCCModelDenseForwardDiff

A `MPCCModelConfig` compiled to produce the derivatives we want.
This version uses ForwardDiff.jl and produces dense matrices.
"""
struct MPCCModelDenseForwardDiff <: AbstractMPCCModelDense
    config::MPCCModelConfig
    f::Function
    f!::Opt{Function}
    ce::Function
    ce!::Opt{Function}
    ci::Function
    ci!::Opt{Function}
    F::Function
    F!::Opt{Function}
    Fq::Function
    Fq!::Opt{Function}
    gradf::Function
    gradf!::Opt{Function}
    gradce::Function
    gradce!::Opt{Function}
    gradci::Function
    gradci!::Opt{Function}
    gradF::Function
    gradF!::Opt{Function}
    gradFq::Function
    gradFq!::Opt{Function}
    hessf::Function
    hessf!::Opt{Function}
    hessce::Function
    hessce!::Opt{Function}
    hessci::Function
    hessci!::Opt{Function}
    hessF::Function
    hessF!::Opt{Function}
    hessFq::Function
    hessFq!::Opt{Function}
    fdp::Function
    fdp!::Opt{Function}
    cedp::Function
    cedp!::Opt{Function}
    cidp::Function
    cidp!::Opt{Function}
    Fdp::Function
    Fdp!::Opt{Function}
    gradfdp::Function
    gradfdp!::Opt{Function}
    gradcedp::Function
    gradcedp!::Opt{Function}
    gradcidp::Function
    gradcidp!::Opt{Function}
    gradFdp::Function
    gradFdp!::Opt{Function}
end




# gradf::Opt{Vector{T}}                       # Vector of length n
# gradce::Opt{Matrix{T}}                      # Changed 20211029: now matrix of me by n
# gradci::Opt{Matrix{T}}                      # Changed 20211029: now matrix of mi by n
# gradF::Opt{Matrix{Vector{T}}}               # Matrix of dim l x q of n length vectors (the gradients wrt x)
# hessf::Opt{Matrix{T}}                       # Matrix dim n x n
# hessce::Opt{Vector{Matrix{T}}}              # Vector of length me of 2d Hessian matrices
# hessci::Opt{Vector{Matrix{T}}}              # Vector of length mi of 2d Hessian matrices
# hessF::Opt{Matrix{Matrix{T}}}               # Matrix of dim l x q of 2d Hessian matrices
# fdp::Opt{Vector{T}}                         # Vector of length r
# cedp::Opt{Vector{Vector{T}}}                # Vector of length r of vectors of length me
# cidp::Opt{Vector{Vector{T}}}                # Vector of length r of vectors of length mi
# Fdp::Opt{Vector{Matrix{T}}}                 # Vector of length r of matrix of dim l x q
# gradfdp::Opt{Vector{Vector{T}}}             # Vector of length r of gradients of f/dp (vectors of length n)
# gradcedp::Opt{Vector{Matrix{T}}}            # Vector of length r of gradce, each /dp
# gradcidp::Opt{Vector{Matrix{T}}}            # Vector of length r of gradci, each /dp
# gradFdp::Opt{Vector{Matrix{Vector{T}}}}     # Vector of length r of gradF, each /dp


struct MPCCModelSparseSym <: AbstractMPCCModelSparse
    config::MPCCModelConfig
    sparsity_gradf::SparseMatrixCSC     # {Bool, Int64}
    sparsity_gradce::SparseMatrixCSC
    sparsity_gradci::SparseMatrixCSC
    sparsity_gradF::Matrix{SparseMatrixCSC}
    sparsity_hessf::SparseMatrixCSC
    sparsity_hessce::Vector{SparseMatrixCSC}
    sparsity_hessci::Vector{SparseMatrixCSC}
    sparsity_hessF::Matrix{SparseMatrixCSC}
    sparsity_gradfdp::Vector{SparseMatrixCSC}
    sparsity_gradcedp::Vector{SparseMatrixCSC}
    sparsity_gradcidp::Vector{SparseMatrixCSC}
    sparsity_gradFdp::Vector{Matrix{SparseMatrixCSC}}
    gradf_sparse_num::SparseMatrixCSC
    gradce_sparse_num::SparseMatrixCSC
    gradci_sparse_num::SparseMatrixCSC
    gradF_sparse_num::Matrix{SparseMatrixCSC}
    hessf_sparse_num::SparseMatrixCSC
    hessce_sparse_num::Vector{SparseMatrixCSC}
    hessci_sparse_num::Vector{SparseMatrixCSC}
    hessF_sparse_num::Matrix{SparseMatrixCSC}
    fdp_sparse_num::Vector
    cedp_sparse_num::Vector
    cidp_sparse_num::Vector
    Fdp_sparse_num::Vector{Matrix}
    gradfdp_sparse_num::Vector{SparseMatrixCSC}
    gradcedp_sparse_num::Vector{SparseMatrixCSC}
    gradcidp_sparse_num::Vector{SparseMatrixCSC}
    gradFdp_sparse_num::Vector{Matrix{SparseMatrixCSC}}
    f::Function
    f!::Opt{Function}
    ce::Function
    ce!::Opt{Function}
    ci::Function
    ci!::Opt{Function}
    F::Function
    F!::Opt{Function}
    # Fq::Function
    # Fq!::Opt{Function}
    gradf::Function
    gradf!::Opt{Function}
    gradce::Function
    gradce!::Opt{Function}
    gradci::Function
    gradci!::Opt{Function}
    gradF::Function
    gradF!::Opt{Function}
    # gradFq::Function
    # gradFq!::Opt{Function}
    hessf::Function
    hessf!::Opt{Function}
    hessce::Function
    hessce!::Opt{Function}
    hessci::Function
    hessci!::Opt{Function}
    hessF::Function
    hessF!::Opt{Function}
    # hessFq::Function
    # hessFq!::Opt{Function}
    fdp::Function
    fdp!::Opt{Function}
    cedp::Function
    cedp!::Opt{Function}
    cidp::Function
    cidp!::Opt{Function}
    Fdp::Function
    Fdp!::Opt{Function}
    gradfdp::Function
    gradfdp!::Opt{Function}
    gradcedp::Function
    gradcedp!::Opt{Function}
    gradcidp::Function
    gradcidp!::Opt{Function}
    gradFdp::Function
    gradFdp!::Opt{Function}
end






struct MPCCJumpFixedFunctions
    jump_f::Function
    jump_f_pen::Function
    jump_ce::Vector{Function}
    jump_ci::Vector{Function}
    jump_F::Matrix{Function}
end


# Only std return functions for the moment
struct MPCCNewtonPenaltyFunctions
    ϕ::Function
    gradϕ::Function
    hessϕ::Function
end
