

# function mpccmodel_test_testvectors_evalpt_basic(model_cfg::MPCCModelConfig, model::AbstractMPCCModel)

#     ptevalreq = MPCCPointEvalReq()

#     nt = length(model_cfg.testvectors)
#     test_results = Vector{Bool}(undef, nt)

#     for lp_t=1:nt

#         # Get the reference point and parameters for this test
#         x0 = model_cfg.testvectors[lp_t].x_val
#         pr0 = model_cfg.testvectors[lp_t].pr_val
#         ps0 = model_cfg.testvectors[lp_t].ps_val

#         # Calculate point using model
#         evalpt = mpccmodel_pointeval_basic(model, ptevalreq, x0, pr0, ps0)

#         # Do the test
#         test_results[lp_t] = model_cfg.testvectors[lp_t].eval_test_pt ≈ evalpt

#     end

#     return test_results
# end