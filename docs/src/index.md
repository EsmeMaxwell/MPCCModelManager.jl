# MPCCModelManager.jl

!!! **DEVELOPMENT VERSION** -- There will be breaking changes to come.


MPCCModelManager.jl is a package for managing parametric MPCC or NLP models of
the form

```math
\begin{equation}
\begin{aligned}
\min_{x(p(t))} \quad & f(p(t), x(p(t))) \\
\textrm{subject to} \quad & c_{\mathcal{E},i}(p(t), x(p(t))) = 0 \quad i = 1, \dots, m_{\mathcal{E}}\\
& c_{\mathcal{I},i}(p(t), x(p(t))) \ge 0 \quad i = 1, \dots, m_{\mathcal{I}} \\
& F_{i,j}(p(t), x(p(t))) \ge 0 \quad i = 1, \dots l, \; j = 1, \dots, q \\
& \prod_{i=1}^{l} F_{i,j} = 0 \quad j = 1, \dots, q.
\end{aligned}
\end{equation}
```


The idea is that one can supply a model definition in various formats such as a
native Julia definition, JSON, etc.

The library will then generate Julia native functions that will calculate the
Jacobians, Hessians, and parameteric derivatives thereof. There is choice in how
the derivatives are computed, e.g. dense matrix outputs via automatic
differentiation (ForwardDiff.jl for now), sparse symbolic differentiation, etc.



## Installation

The package is not registered yet, so this must be done manually at the moment.

## Basic Usage Description

Each model has a dimension specification `MPCCDimSpec`, which details:

* `n`: The number of spatial coordinates.
* `l`: The number of expressions in each complementarity constraint.
* `q`: The numebr of complementarity constraints.
* `me`: The number of equality constraints.
* `mi`: The number of inequality constraints.
* `r`: The number of continuous parameters, which are assumed implicitly dependent on some real scalar, say "t".
* `s`: The number of integer parameters, not used yet.

There are some sample models under `./src/model_samples`. In the future, this
will be in a different package and there'll be functionality to import from
JSON, etc. For the moment, when defining you own models, just follow the same
pattern.

Instantiate the model configuration, e.g. using the model kj3:

`model_cfg = model_pmpcc_kj3_build();`

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

This returns a struct which includes the model configuration, and also the following functions (standard and mutating).

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

Each of the functions accepts a primal position `x`, continuous parameters `pr` ($p$ in the model), and discrete parameters `ps`.


