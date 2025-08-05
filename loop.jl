#-----------------------------------------------------------------------#
# Iterate over different environmental conditions and parameter values  #
#-----------------------------------------------------------------------#

# July 28, 2025; by Suzana Goncalves Leles

# Save variable list and names
variable_list = [logμ, β, ptr, pri, plb, pdp, pp, pru, pld, pgl, pre, ϕtr, ϕri, ϕlb, ϕp, ϕru, ϕld, ϕgl, ϕre, cin, clm, cic, αld, ctag, αlm, αtag, Qcell, Qc, Qn, vru, vres, vp, vgl, vd, vre, vlb, vld, vri, vtr, vol, tot_p, D, γTa, γTd, γTad, γTm, cli, cgu, ecost, ηaan, ϵ, T, DIN, I, Ea, ctag_max]
variable_strings = ["μ", "β", "ptr", "pri", "plb", "pdp", "pp", "pru", "pld", "pgl", "pre", "ϕtr", "ϕri", "ϕlb", "ϕp", "ϕru", "ϕld", "ϕgl", "ϕre", "cin", "clm", "cic", "αld", "ctag", "αlm", "αtag", "Qcell", "Qc", "Qn", "vru", "vres", "vp", "vgl", "vd", "vre", "vlb", "vld", "vri", "vtr", "vol", "tot_p", "D", "γTa", "γTd", "γTad", "γTm", "cli", "cgu", "ecost", "ηaan", "ϵ",  "T", "DIN", "I", "Ea", "ctag_max"]

# Define empty dictionary to save values
s = Dict()

# Define loop function based on run mode ("thermal_curve" or "latitudinal")

if run_mode == "thermal_curve"
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

elseif run_mode == "latitudinal"
    for i in 1:length(change_DIN)
    set_value(DIN, change_DIN[i])
    #for j in 1:length(change_T)
            #set_value(T, change_T[j])
            set_value(T, change_T[i])
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
                                s[variable_strings[n]] = Array{Any,4}(undef,length(change_DIN),1,length(change_Ea),length(change_ctag_max))
                                #s[variable_strings[n]] = Array{Any,4}(undef,length(change_T),length(change_DIN),length(change_Ea),length(change_ctag_max))
                            end
                            #s[variable_strings[n]][j,i,k,l] = model_val
                            s[variable_strings[n]][i,1,k,l] = model_val
                    end
                end
            end
    #end
end

else
    error("Unrecognized mode: $mode")
end
