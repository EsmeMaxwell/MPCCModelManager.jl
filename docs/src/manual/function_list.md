# Functions and structs

## From ``mpccmodel_proc_model.jl``

These functions are used to process a model definition with end point producing
the model config struct which includes closures to construct `f`, `ce`, `ci`,
`F` functions, parametrisations, etc.

```@docs
MPCCModelMAnager.mpccmodel_build_sym_nums
```

```@docs
MPCCModelManager.mpccmodel_build_fn_from_defn
```

```@docs
MPCCModelManager.mpccmodel_build_parameterisation
```

```@docs
MPCCModelManager.mpccmodel_construct_config
```


## From ``mpccmodel_calc_forwarddiff_dense.jl``

This construsts the full gamut of functions as closures using ForwardDiff.jl.

```@docs
MPCCModelManager.mpccmodel_setup_forwarddiff_dense
```


## From ``mpccmodel_constraint_sets.jl``

This code is used when storing the active set description. It allow manipulation of an active set description in a storage efficient manner.

```@docs
MPCCModelManager.cs_cpair_to_bi
```

```@docs
MPCCModelManager.cs_bi_to_cpair
```

```@docs
MPCCModelManager.cs_cnstr_el_to_bi
```

```@docs
MPCCModelManager.cs_bi_to_cnstr_el
```

```@docs
MPCCModelManager.cs_cnstr_el_to_cpair
```

```@docs
MPCCModelManager.cs_cpair_to_cnstr_el
```

```@docs
MPCCModelManager.cs_cnstr_idxs_add_by_bi!
```

```@docs
MPCCModelManager.cs_cnstr_idxs_add_by_cpair!
```

```@docs
MPCCModelManager.cs_cnstr_idxs_del_by_cnstr_el!
```

```@docs
MPCCModelManager.cs_cnstr_idxs_del_by_bi!
```

```@docs
MPCCModelManager.cs_cnstr_idxs_del_by_cpair!
```

```@docs
MPCCModelManager.cs_cnstr_get_as_bitmask
```

```@docs
MPCCModelManager.cs_cnstr_build_from_bitmask
```

```@docs
MPCCModelManager.cs_cnstr_build_all_active
```

```@docs
MPCCModelManager.cs_cnstr_get_Fq_count
```

```@docs
MPCCModelManager.cs_cnstr_get_all_active_bi
```
