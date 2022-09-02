module MPCCModelManager

using SparseArrays, ForwardDiff, LinearAlgebra, UnPack, RuntimeGeneratedFunctions, SymbolicUtils, Symbolics

const MPCCMM_ROOT_DIR = @__DIR__

RuntimeGeneratedFunctions.init(@__MODULE__)

# Comment that should only be on test-1 branch

# Shit we might provide a method for using our custom types
import Base.show



include("mpccmodel_common.jl")
include("mpccmodel_proc_model.jl")
include("mpccmodel_calc_forwarddiff_dense.jl")
include("mpccmodel_calc_symdiff_sparse.jl")
include("mpccmodel_calc_forwarddiff_newton_penalty.jl")
include("mpccmodel_constraint_sets.jl")
include("mpccmodel_pointeval.jl")
include("mpccmodel_inc_samples.jl")
include("mpccmodel_test.jl")


# mpccmodel_proc_model.jl
export  mpccmodel_build_sym_nums,
        mpccmodel_build_fn_from_defn,
        mpccmodel_build_parameterisation,
        mpccmodel_construct_config,
        mpccmodel_load_defn_from_file,
        mpccmodel_build_fixed_jump_fns


# mpccmodel_calc_forwarddiff_dense.jl
export  mpccmodel_setup_forwarddiff_dense

# mpccmodel_calc_symdiff_sparse.jl
export  mpccmodel_setup_symdiff_sparse

# mpccmodel_setup_newton_penalty.jl
export  mpccmodel_setup_newton_penalty

# mpccmodel_constraint_sets.jl
export  cs_cpair_to_bi,
        cs_bi_to_cpair,
        cs_cnstr_el_to_bi,
        cs_bi_to_cnstr_el,
        cs_cnstr_el_to_cpair,
        cs_cpair_to_cnstr_el,
        cs_cnstr_idxs_add_by_bi!,
        cs_cnstr_idxs_add_by_cpair!,
        cs_cnstr_idxs_del_by_cnstr_el!,
        cs_cnstr_idxs_del_by_bi!,
        cs_cnstr_idxs_del_by_cpair!,
        cs_cnstr_get_as_bitmask,
        cs_cnstr_build_from_bitmask,
        cs_cnstr_build_all_active,
        cs_cnstr_get_Fq_count,
        cs_cnstr_get_all_active_bi



# TODO
# - remove line numbers from generated functions, inline what need inlined, remove bounds checks where appropriate

# - change notation for jacce, etc, to jacce and then write proper jacce

# TODO perhaps make more generic by defining all the operations in a custom sense.  Perhaps a 2.0 type of thing.


# NOTE reminders to self:
# - Type stability: return types should not depend on values but can depend on input types
# - Use argument types for multiple dispatch, it does not improve performance in general
# - Explicit return type annotation seems to degrate performance

# TODO should probably inline the small functions in each model definition
# TODO check ForwardDiff cfg
# TODO saving of lower order derivatives?



end