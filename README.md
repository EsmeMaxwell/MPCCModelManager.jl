# MPCCModelManager

!!! **DEVELOPMENT VERSION** -- Only partially implemented and there will be
breaking changes to come.

MPCCModelManager.jl is a package for managing parametric MPCC or NLP models of
the form

$$
\begin{equation}
\begin{aligned}
\min_{x(p(t))} \quad & f(p(t), x(p(t))) \\
\textrm{subject to} \quad & c_{\mathcal{E},i}(p(t), x(p(t))) = 0 \quad i = 1, \dots, m_{\mathcal{E}}\\
& c_{\mathcal{I},i}(p(t), x(p(t))) \ge 0 \quad i = 1, \dots, m_{\mathcal{I}} \\
& F_{i,j}(p(t), x(p(t))) \ge 0 \quad i = 1, \dots l, \; j = 1, \dots, q \\
& \prod_{i=1}^{l} F_{i,j} = 0 \quad j = 1, \dots, q.
\end{aligned}
\end{equation}
$$


The idea is that one can supply a model definition in various formats such as a
native Julia definition, JSON, etc.

The library will then generate Julia native functions that will calculate the
Jacobians, Hessians, and parameteric derivatives thereof. There is choice in how
the derivatives are computed, e.g. dense matrix outputs via automatic
differentiation (ForwardDiff.jl for now), sparse symbolic differentiation, etc.

As of November 2022, there is functionality to import from Julia defined models and compile using FowardDiff.


There's some vague initial documentation here: https://petermaxwell.github.io/MPCCModelManager.jl/dev/



## Installation

The package is not registered yet, so this must be done manually at the moment.


