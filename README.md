# ProteomePhyto

## A proteome allocation model for a phytoplankton cell

### The model

Our coarse-grained model is built on a constrained optimization problem that aims to maximize the steady-state growth rate of a generic phytoplankton cell considering temperature, irradiance levels, and external dissolved inorganic nitrogen and carbon concentrations. 

The model is written in Julia and contains three files: 
1) par.jl defines the parameter values
2) eqn.jl defines the model equations
3) save.jl runs the model over environmental conditions
4) reopt.jl identifies sub-optimal solutions and helps the optimizer by providing more accurate initial guesses
5) run.jl calls all files in order to run the model and save output.

Julia packages required to run the code above: JuMP, Ipopt, CSV, DataFrames, Plots

Code developed by Suzana Goncalves Leles and Lara Breithaup.
