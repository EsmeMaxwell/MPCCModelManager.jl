module MPCCModelManager

# TODO perhaps make more generic by defining all the operations in a custom sense.  Perhaps a 2.0 type of thing.


using SparseArrays, ForwardDiff, LinearAlgebra, UnPack, RuntimeGeneratedFunctions, SymbolicUtils, Symbolics
# StaticArrays

const MPCCMM_ROOT_DIR = @__DIR__

RuntimeGeneratedFunctions.init(@__MODULE__)

@warn "Changed gradients from column to row vectors"


# Shit we might provide a method for using our custom types
import Base.show



include("mpccmodel_common.jl")
include("mpccmodel_proc_model.jl")
include("mpccmodel_calc_forwarddiff_dense.jl")
include("mpccmodel_calc_symdiff_sparse.jl")
include("mpccmodel_calc_forwarddiff_newton_penalty.jl")
include("mpccmodel_pointeval.jl")
include("mpccmodel_test.jl")


export  mpccmodel_build_sym_nums,
        mpccmodel_build_fn_from_defn,
        mpccmodel_construct_config,
        mpccmodel_load_defn_from_file,
        mpccmodel_build_fixed_jump_fns,
        mpccmodel_setup_forwarddiff_dense,
        mpccmodel_setup_symdiff_sparse,
        mpccmodel_setup_newton_penalty




# More recent note: having pr and ps, then gradp, etc, is a bit werid.  cannot take derivative wrt ps, so perhaps change notation

# TODO
# - remove line numbers from generated functions, inline what need inlined, remove bounds checks where appropriate
# - deal with transpose for ce, ci.
# - remove f() in the knownsolutions eval, unncessary and creates a headache


 # - for forwarddiff processing, cast results properly

# NOTE(i): We leave the return types fairly vague and allow caller to specify whether they want sparse, dense, or even static arrays.


# NOTE reminders to self:
# - Type stability: return types should not depend on values but can depend on input types
# - Use argument types for multiple dispatch, it does not improve performance in general
# - Explicit return type annotation seems to degrate performance

# TODO should probably inline the small functions in each model definition
# TODO check ForwardDiff cfg
# TODO saving of lower order derivatives?






end