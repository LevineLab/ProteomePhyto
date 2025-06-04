# ProteomePhyto

## A proteome allocation model for a phytoplankton cell

### The model

Our coarse-grained model is built on a constrained optimization problem that aims to maximize the steady-state growth rate of a generic phytoplankton cell considering temperature, irradiance levels, and external dissolved inorganic nitrogen and carbon concentrations. 

The model is written in Julia (Version 1.5.3) and contains three files: 
1) par.jl defines the parameter values
2) eqn.jl defines the model equations
3) save.jl runs the model over environmental conditions
4) reopt.jl identifies local minima and helps the optimizer to find the global optima by providing more accurate initial guesses
5) run.jl calls all files in order to run the model and save output.

Julia installation guidelines can be found here (https://docs.julialang.org/en/v1/manual/installation/) and the required packages to run the code above are: JuMP, Ipopt, CSV, DataFrames, Plots

Code developed by Suzana Goncalves Leles and Lara Breithaup.

### Parameterization of the model against the dataset compiled by Anderson et al. (2021)

We parameterized the model to represent two phytoplankton functional types, cyanobacteria and diatoms. To confirm that our parameter choices mimic the general physiology of cyanobacteria, cold-adapted diatoms, and warm-adapted diatoms, we compared the emergent maximum growth rates and growth (left) slopes of the thermal curves from our model against data from 167 species previously compiled by Anderson et al. (2021). The code for this analysis was developed by Arianna I Krinos and is given as the file estimating_growth_slopes_from_anderson_dataset.R. Here we rpovide a sample data to illustrate how this code works (sample_data.csv).

Anderson, S., Barton, A., Clayton, S., Dutkiewicz, S., and Rynearson, T. (2021). Marine phytoplankton functional types exhibit diverse responses to thermal change. Nature communications, 12(1):6413.404
