#-----------------------------------------------------------------------#
# Iterate over different environmental conditions and save model output #
#-----------------------------------------------------------------------#

# January 2025

# Importing output from Earth System model
using CSV
using DataFrames

# Optimal solutions looping over different environmental conditions.
env = "Temperature (C°)" # should be one of the following: Light, Temperature (C°) or [DIN]

# Set the environmenal parameters (we will run the model for all combinations for the values set below):
change_DIN = [20.0,0.01]
change_I = 20.0
change_T = 273.0 .+ [5.0:1.0:36;] 

# Set the physiological parameters
change_Ea = collect(range(0.5, stop=1.0, length=10)) 
change_ctag_max = collect(range(2e6, stop = 2e7, length = 10))

# Save variable list and names
variable_list = [logμ, β, ptr, pri, plb, pdp, pp, pru, pld, pgl, pre, ϕtr, ϕri, ϕlb, ϕp, ϕru, ϕld, ϕgl, ϕre, cin, clm, cic, αld, ctag, αlm, αtag, Qcell, Qc, Qn, vru, vres, vp, vgl, vd, vre, vlb, vld, vri, vtr, vol, tot_p, D, γTa, γTd, γTad, γTm, cli, cgu, ecost, ηaan, ϵ, T, DIN, I, Vnu, ctag_max]
variable_strings = ["μ", "β", "ptr", "pri", "plb", "pdp", "pp", "pru", "pld", "pgl", "pre", "ϕtr", "ϕri", "ϕlb", "ϕp", "ϕru", "ϕld", "ϕgl", "ϕre", "cin", "clm", "cic", "αld", "ctag", "αlm", "αtag", "Qcell", "Qc", "Qn", "vru", "vres", "vp", "vgl", "vd", "vre", "vlb", "vld", "vri", "vtr", "vol", "tot_p", "D", "γTa", "γTd", "γTad", "γTm", "cli", "cgu", "ecost", "ηaan", "ϵ",  "T", "DIN", "I", "Vnu", "ctag_max"]

# Define empty dictionary to save values
s = Dict()

# Store model output in empty vectors in dictionary
for i in 1:length(change_DIN)
    set_value(DIN, change_DIN[i])
    for j in 1:length(change_T)
            set_value(T, change_T[j])
            for k in 1:length(change_Ea)
                set_value(Ea, change_Ea[k])
                for l in 1:length(change_ctag_max)
                    set_value(ctag_max, change_ctag_max[l])
                    optimize!(model)
                    # set value in variable list when value calculation changes
                    change_val = findall( x -> x == "Qcell", variable_strings)[1]
                    # loop through variables and save the values in dictionary
                    for n in 1:length(variable_list)
                            model_val = 0
                            if n == 1
                                model_val =  exp.(JuMP.value.(variable_list[n]))
                            elseif n < change_val
                                model_val = JuMP.value.(variable_list[n])
                            else
                                model_val = value(variable_list[n])
                            end
                            # add value to dictionary
                            if !(variable_strings[n] in keys(s))
                                    s[variable_strings[n]] = Array{Any,4}(undef,length(change_T),length(change_DIN),length(change_Ea),length(change_ctag_max))
                            end
                            s[variable_strings[n]][j,i,k,l] = model_val
                    end
                end
            end
    end
end
