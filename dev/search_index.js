var documenterSearchIndex = {"docs":
[{"location":"manual/function_list/#Functions-and-structs","page":"Functions and structs","title":"Functions and structs","text":"","category":"section"},{"location":"manual/function_list/","page":"Functions and structs","title":"Functions and structs","text":"MPCCModelManager.mpccmodel_build_fn_from_defn","category":"page"},{"location":"manual/function_list/#MPCCModelManager.mpccmodel_build_fn_from_defn","page":"Functions and structs","title":"MPCCModelManager.mpccmodel_build_fn_from_defn","text":"mpccmodel_build_fn_from_defn(dimspec, defnnum, x, pr, ps)\n\nAssemble f, ce, ci, and F Julia native functions from Symbolic definition and returns both non-indexed and indexed variants. Generally no need to call directly as this is called from mpccmodel_load_defn_from_file().\n\nNOTE: To understand how this works, one should also play with a few examples ofbuild_function() in Symbolics.jl.  The output from build_function() is a 2-element vector whereby the first element is the standard function adn the second is a mutating version of the same function. Hence the formulation in the code here.\n\nArguments\n\ndimspec::MPCCDimSpec: The dimension specifications for the proram.\ndefnnum::MPCCDefinition: The Symbolics.jl Num definitions of the expressions.\nx, pr, ps: Symbolics.jl Num variables for x, pr, ps.\n\n\n\n\n\n","category":"function"},{"location":"#MPCCModelManager.jl","page":"Home","title":"MPCCModelManager.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"!!! DEVELOPMENT VERSION – There will be breaking changes to come.","category":"page"},{"location":"","page":"Home","title":"Home","text":"MPCCModelManager.jl is a package for managing parametric MPCC or NLP models. It takes a symbolic definition and produces a struct containing Julia functions that can compute standard derivatives (Jacobians, Hessians) including derivatives with respect to the parameteres. A subset of constraints can be specified in each function call.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package is not registered yet, so this must be done manually at the moment.","category":"page"},{"location":"#Usage-Description","page":"Home","title":"Usage Description","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Each model has a dimension specification MPCCDimSpec, which details:","category":"page"},{"location":"","page":"Home","title":"Home","text":"n: The number of spatial coordinates.\nl: The number of expressions in each complementarity constraint.\nq: The numebr of complementarity constraints.\nme: The number of equality constraints.\nmi: The number of inequality constraints.\nr: The number of continuous parameters, which are assumed implicitly dependent on some real scalar, say \"t\".\ns: The number of integer parameters, not used yet.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Instantiate the model configuration.  The string parameter determines which file to load.","category":"page"},{"location":"","page":"Home","title":"Home","text":"model_cfg = mpccmodel_load_defn_from_file(\"kj6\");","category":"page"},{"location":"","page":"Home","title":"Home","text":"This struct contains, amongst other things, the following fields:","category":"page"},{"location":"","page":"Home","title":"Home","text":"dimspec: Dimensions of the model.\nx: Symbolics Num variables vector for x (spatial).\npr: Symbolics Num variables vector for parametric variables (actually functions dependant on some supplied \"t\".\nps: Symbolics Num variables vector for parametric integer variables (not really used).\ndefn: Symbolic expressions for the model including f, ce, ci, F.\nfns: Julia native compiled functions mapping directly to defn.\nknownsols: For test problems, there may be a closed expression for the parametric solution.\nparameterisations: Defines the functional relationship of pr to some real scalar (this is important when r > 1).","category":"page"},{"location":"","page":"Home","title":"Home","text":"Now, from the model config, we'd like to instantiate the functions that will actually do the work for us. In the future, there will be a choice here, e.g. there'll be the option of symbolic differentiation producing sparse output.  For now though, it's only automatic differentiation implemented against ForwardDiff.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"model_fd = mpccmodel_setup_forwarddiff_dense(model_cfg);","category":"page"},{"location":"","page":"Home","title":"Home","text":"This returns a struct which includes the model configuration, and also the following functions (standard and mutating)","category":"page"},{"location":"","page":"Home","title":"Home","text":"f, f!: Objective function (scalar, I think).\nce, ce!: Equality constraints (vector).\nci, ci!: Inqquality constraints (vector).\nF, F!: Comlementarity constraints (matrix).\nFq, Fq!: Product of complementarity contraints (vector), for penalty method.\ngradf, gradf!: Row vector gradient of f (vector).\ngradce, gradce!: Jacobian of ce, (me by n matrix).\ngradci, gradci!: Jacobian of ci, (mi by n matrix).\ngradF, gradF!:  Gradients of F, (l by q matrix of length n vectors).\ngradFq, gradFq!: \nhessf, hessf!: Hessian of f.\nhessce, hessce!: Hessians of ce (length me vector of n by n matrices)\nhessci, hessci!: Hessians of ci (length mi vector of n by n matrices)\nhessF, hessF!: Hessians of F (l by q matrix of n by n matrices)\nhessFq, hessFq!: \nfdp, fdp!: Gradient of f wrt pr (vector of length r)\ncedp, cedp!: Gradient of ce wrt pr (vector of length r of vector of length me)\ncidp, cidp!: Gradient of ci wrt pr (vector of length r of vector of length mi)\nFdp, Fdp!: Gradient of F wrt pr (vector of length r of matrix of size l by q)\ngradfdp, gradfdp!: Gradient wrt pr of gradient wrt x (vector of length r of vector of length n)\ngradcedp, gradcedp!: Gradient wrt pr of Jacobian of ce (vector of length r of matrix of size me by n)\ngradcidp, gradcidp!: Gradient wrt pr of Jacobian of ci (vector of length r of matrix of size mi by n)\ngradFdp, gradFdp!: Gradient wrt pr of spatial gradients of F (vector of length r or matrix of size l by q of vector of length n)","category":"page"}]
}
