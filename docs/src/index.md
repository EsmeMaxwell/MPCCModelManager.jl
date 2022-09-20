# MPCCModelManager.jl

!!! **DEVELOPMENT VERSION** -- There will be breaking changes to come.

MPCCModelManager.jl is a package for managing parametric MPCC or NLP models. It
takes a symbolic definition and produces a struct containing Julia functions
that can compute standard derivatives (Jacobians, Hessians) including
derivatives with respect to the parameteres. A subset of constraints can be
specified in each function call.

## Installation

The package is not registered yet, so this must be done manually at the moment.

## Usage Description

Each model has a dimension specification `MPCCDimSpec`, which details:

* `n`: The number of spatial coordinates.
* `l`: The number of expressions in each complementarity constraint.
* `q`: The numebr of complementarity constraints.
* `me`: The number of equality constraints.
* `mi`: The number of inequality constraints.
* `r`: The number of continuous parameters, which are assumed implicitly dependent on some real scalar, say "t".
* `s`: The number of integer parameters, not used yet.

Instantiate the model configuration. There are some samples included in the
package (this functionality can be used at runtime though, without world-age
issues)

`model_cfg = model_pmpcc_kj6_build();`

This struct contains, amongst other things, the following fields:

* `dimspec`: Dimensions of the model.
* `x`: Symbolics Num variables vector for `x` (spatial).
* `pr`: Symbolics Num variables vector for parametric variables (actually functions dependant on some supplied "t".
* `ps`: Symbolics Num variables vector for parametric integer variables (not really used).
* `defn`: Symbolic expressions for the model including `f`, `ce`, `ci`, `F`.
* `fns`: Julia native compiled functions mapping directly to `defn`.
* `knownsols`: For test problems, there may be a closed expression for the parametric solution.
* `parameterisations`: Defines the functional relationship of `pr` to some real scalar (this is important when `r > 1`).

Now, from the model config, we'd like to instantiate the functions that will
actually do the work for us. In the future, there will be a choice here, e.g.
there'll be the option of symbolic differentiation producing sparse output.  For
now though, it's only automatic differentiation implemented against
ForwardDiff.jl.

`model_fd = mpccmodel_setup_forwarddiff_dense(model_cfg);`

This returns a struct which includes the model configuration, and also the following functions (standard and mutating)

* `f`, `f!`: Objective function (scalar, I think).
* `ce`, `ce!`: Equality constraints (vector).
* `ci`, `ci!`: Inqquality constraints (vector).
* `F,` `F!`: Comlementarity constraints (matrix).
* `Fq`, `Fq!`: Product of complementarity contraints (vector), for penalty method.
* `gradf`, `gradf!`: Row vector gradient of `f` (vector).
* `jacce`, `jacce!`: Jacobian of `ce`, (`me` by `n` matrix).
* `jacci`, `jacci!`: Jacobian of `ci`, (`mi` by `n` matrix).
* `gradF`, `gradF!`:  Gradients of `F`, (`l` by `q` matrix of length `n` vectors).
* `gradFq`, `gradFq!`: 
* `hessf`, `hessf!`: Hessian of `f`.
* `hessce`, `hessce!`: Hessians of `ce` (length `me` vector of `n` by `n` matrices)
* `hessci`, `hessci!`: Hessians of `ci` (length `mi` vector of `n` by `n` matrices)
* `hessF`, `hessF!`: Hessians of `F` (`l` by `q` matrix of `n` by `n` matrices)
* `hessFq`, `hessFq!`: 
* `fdp`, `fdp!`: Gradient of `f` wrt `pr` (vector of length `r`)
* `cedp`, `cedp!`: Gradient of `ce` wrt `pr` (vector of length `r` of vector of length `me`)
* `cidp`, `cidp!`: Gradient of `ci` wrt `pr` (vector of length `r` of vector of length `mi`)
* `Fdp`, `Fdp!`: Gradient of `F` wrt `pr` (vector of length `r` of matrix of size `l` by `q`)
* `gradfdp`, `gradfdp!`: Gradient wrt `pr` of gradient wrt `x` (vector of length `r` of vector of length `n`)
* `jaccedp`, `jaccedp!`: Gradient wrt `pr` of Jacobian of `ce` (vector of length `r` of matrix of size `me` by `n`)
* `jaccidp`, `jaccidp!`: Gradient wrt `pr` of Jacobian of `ci` (vector of length `r` of matrix of size `mi` by `n`)
* `gradFdp`, `gradFdp!`: Gradient wrt `pr` of spatial gradients of `F` (vector of length `r` or matrix of size `l` by `q` of vector of length `n`)


