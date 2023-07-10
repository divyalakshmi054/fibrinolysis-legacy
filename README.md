# S-system-fibrinolysis-thrombin-interpolated

To run simulations: Edit visit number and tPA level in `sample_ensemble.jl`, then `include("sample_ensemble.jl")`

To sort real TEG data:   `include("sort_real_teg_data.jl")`

To plot: Make sure you have run `sample_ensemble.jl` and `sort_real_teg_data.jl' or have the resulting data files in your directory. Edit the visit number and tPA level in `generate_plots.jl`, then `include("generate_plots.jl")`

To run Sobol sensitivity analysis: Edit number of samples, bootstrap runs, or confidence interval in `sens_sobol.jl` then `include("sens_sobol.jl")`

To plot Sobol indices: `include("plot_sobol.jl")`; you may need to update the CSV filename or the path to file containing the Sobol results
