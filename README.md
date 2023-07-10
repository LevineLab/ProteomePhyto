# ProteomePhyto

## A proteome allocation model for a phytoplankton cell

### The model

Our coarse-grained model is built on a constrained optimization problem that aims to maximize the steady-state growth rate of a generic phytoplankton cell considering temperature, irradiance levels, and external dissolved inorganic nitrogen and carbon concentrations. 

The model is written in Julia and contains three files: 
1) par.jl defines the parameter values
2) eqn.jl defines the model equations
3) save.jl runs the model over environmental conditions and saves model output
4) plot.jl plots model output and saves it as a pdf file

Julia packages required to run the code above: JuMP, Ipopt, CSV, DataFrames, Plots
