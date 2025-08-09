# ProteomePhyto

## A proteome allocation model for a phytoplankton cell

### The model

Our coarse-grained model is built on a constrained optimization problem that aims to maximize the steady-state growth rate of a generic phytoplankton cell considering temperature, irradiance levels, and external dissolved inorganic nitrogen and carbon concentrations. 

The model is written in Julia (Version 1.5.3) and contains three files: 
1) par.jl defines the parameter values
2) eqn.jl defines the model equations
3) loop.jl iterates the model over different environmental conditions/parameters
4) reopt.jl identifies local minima and helps the optimizer to find the global optima by providing more accurate initial guesses
5) run.jl calls all files in order to run the model and save output.

Julia installation guidelines can be found here (https://docs.julialang.org/en/v1/manual/installation/) and the required packages to run the code above are: JuMP, Ipopt, CSV, DataFrames

Code developed by Suzana Goncalves Leles and Lara Breithaup.

### Step-by-step to run the model
Open run.jl and follow the steps listed in the file. Briefly, you will need to define:
- type: takes one of 4: "pro" or "syn" or "warm_diatom" or "cold_diatom"
- run_mode: takes one of 2: "thermal_curve" or "latitudinal"
- time_period: takes one of two: "present" or "future"

### Plot model output and generate figures
Model output is organized in different folders:
- output_thermal_curves: contains model output for the thermal curves for each modeled functional type under nutrient-replete and nutrient-deplete conditions. Run the R script "script_thermal_curves.R" to plot all figures in the manuscript associated with this output.
- output_latitudinal: contains model output for the latitudinal simulations for each modeled functional type under present and future climate scenarios. Run the R script "script_lat.R" to plot all figures in the manuscript associated with this output.
  
### Parameterization of the model against the dataset compiled by Anderson et al. (2021)

We parameterized the model to represent four phytoplankton functional types: Prochlorococcus, Synechococcus, warm-adapted diatom, and cold-adapted diatom. To confirm that our parameter choices mimic the general physiology of these groups, we compared the emergent maximum growth rates and growth (left) slopes of the thermal curves from our model against data from 167 species previously compiled by Anderson et al. (2021). The code for this analysis was developed by Arianna I Krinos and is given as the file estimating_growth_slopes_from_anderson_dataset.R. Here we also apply this code to estimate the slope from our modeled thermal growth curves.

Anderson, S., Barton, A., Clayton, S., Dutkiewicz, S., and Rynearson, T. (2021). Marine phytoplankton functional types exhibit diverse responses to thermal change. Nature communications, 12(1):6413.404


