#-------------------------------------------------------------------------#
# Iterate over 'bad' model values and repoptimize to correct model output #
#-------------------------------------------------------------------------#

# April 2024; by Lara Breithaupt
using LinearAlgebra
using Statistics

function first_good(growth, size, qcell)
    for temp in 1:length(growth)
        if ( ((growth[temp] > .0001) && (growth[temp] < .001) ) && (size[temp] > 0.1/6) && (size[temp] < 100/6) &&  ((qcell[temp] > 0.005) && (qcell[temp] < 0.5)) )
            return temp
        end
    end
    return 1
end

function bad_vals(x, y, z, w)
    # define empty arrays for new data
    bad_indexes = Array{Float64, 1}(undef, 0)
    # store first two values for reference
    # calculate slope (or use standard difference)
    last_slope = 1.213e-5 
    #last_slope = abs(y[2] - y[1])
    # define the last point that is a 'good' data point
    last_good = first_good(y, z, w)
    # since_last_good = 0
    # loop through data to determine if the next point is 'good' data
    for i in (last_good + 1):(length(y))
        curr_slope = abs(y[i] - y[last_good]) / (x[i] - x[last_good])
        if ( ( curr_slope < 8*last_slope ) && ( (z[i] > 0.1/6) && (z[i] < 100/6) ) && ( (w[i] > 0.005) && (w[i] < 0.5) ) ) #&& curr_slope > .125*last_slope
            last_good = i
            last_slope = curr_slope
        else
            push!(bad_indexes, i)
        end
    end
    return bad_indexes
end

function find_highest_good_index(index, bad_vals)
    # for second optimization, find last good value for start values
    valid_indices = [i for i in 1:index-1 if i ∉ bad_vals]
    if isempty(valid_indices)
        return [index - 1]
    else
        return valid_indices #maximum(valid_indices)
    end
end

function find_lowest_good_index(index, bad_vals)
    # for third optimization, find next good value for start values
    valid_indices = [i for i in index:length(change_T) if i ∉ bad_vals]
    if isempty(valid_indices)
        return [index]
    else
        return valid_indices #minimum(valid_indices)
    end
end

function find_comp(r, lower_indecies, higher_indecies)
    # First run closest below
    if r == 1
        comp = maximum(lower_indecies)
    # Second run closest above
    elseif r == 2
        comp = minimum(higher_indecies)
    # Third run second closest below
    elseif r == 3 && length(lower_indecies) > 2
        comp = lower_indecies[length(lower_indecies) - 1]
    elseif r == 3 && length(lower_indecies) < 2
        comp = maximum(lower_indecies) #CHANGE?
    # Fourth run second closest above
    elseif r == 4 && length(higher_indecies) > 2
        comp = higher_indecies[2]
    elseif r == 4 && length(higher_indecies) < 2
        comp = maximum(higher_indecies) #CHANGE?
    # Fifth run third closest above or third closest 
    elseif r == 5 && length(higher_indecies) > 3#&& length(bad_indecies) > 5
        comp = higher_indecies[3]
    elseif r == 5 && length(lower_indecies) > 3
        comp = lower_indecies[length(lower_indecies) - 2]
    elseif r == 5
        comp = minimum(lower_indeces)
        # if length(higher_indecies) > 3
        #     comp = higher_indecies[3]
        # elseif length(lower_indecies) > 3
        #     comp = lower_indecies[length(lower_indecies) - 2]
        # else
        #     comp = maximum(higher_indecies) #CHANGE?
        # end
    else
        comp = maximum(higher_indecies) #CHANGE?
    end

    return comp
end

# Define variable list, v, and variable strings, vs
vars = [logμ, β, ptr, pri, plb, pdp, pp, pru, pld, pgl, pre, ϕtr, ϕri, ϕlb, ϕp, ϕru, ϕld, ϕgl, ϕre, cin, clm, cic, αld, ctag, αlm, αtag, Qcell, Qc, Qn, vru, vres, vp, vgl, vd, vre, vlb, vld, vri, vtr, vol, tot_p, D, γTa, γTd, γTad, γTm, cli, cgu, ecost, ηaan, ϵ, T, DIN, I, Ea, ctag_max]
vs = ["μ", "β", "ptr", "pri", "plb", "pdp", "pp", "pru", "pld", "pgl", "pre", "ϕtr", "ϕri", "ϕlb", "ϕp", "ϕru", "ϕld", "ϕgl", "ϕre", "cin", "clm", "cic", "αld", "ctag", "αlm", "αtag", "Qcell", "Qc", "Qn", "vru", "vres", "vp", "vgl", "vd", "vre", "vlb", "vld", "vri", "vtr", "vol", "tot_p", "D", "γTa", "γTd", "γTad", "γTm", "cli", "cgu", "ecost", "ηaan", "ϵ",  "T", "DIN", "I", "Ea", "ctag_max"]

# Define index in v array where values change from variables to NL expressions
change_val = findall(x -> x == "Qcell", vs)[1]

reopt = copy(s)

for r in 1:5
    # Loop through the DIN, Ea, and ctag_max indecies and reoptimize the model to correct for the values
    for n in 1:length(change_DIN)
        for v in 1:length(change_Ea)
            for c in 1:length(change_ctag_max)
                # Define the bad indices for the first value in variable list (v) (will extrapolate to all)
                bad_indecies = bad_vals(change_T, reopt["μ"][:, n, v, c], reopt["β"][:, n, v, c], reopt["Qcell"][:, n, v, c])
                # Loop through the 'bad' solutions and reoptimize them
                for b in 1:length(bad_indecies)
                    index = floor(Int, bad_indecies[b])
                    # Reoptimize for each variable
                    # Set start value for all state variables (i.e. not NL expressions)
                    higher_indecies = find_lowest_good_index(index, bad_vals(change_T, reopt["μ"][:, n, v, c], reopt["β"][:, n, v, c], reopt["Qcell"][:, n, v, c]))
                    lower_indecies = find_highest_good_index(index, bad_vals(change_T, reopt["μ"][:, n, v, c], reopt["β"][:, n, v, c], reopt["Qcell"][:, n, v, c]))
                    # Set the comparision value based on the run
                    comp = find_comp(r, lower_indecies, higher_indecies)
                    # Set the start value for each varaiable
                    for var in 1:(change_val - 1)
                        # Set start value for each variable as previous + delta
                        set_start_value.(vars[var], s[vs[var]][comp, n, v, c])
                    end
                    # Set Temp start value and nut value
                    set_value(T, change_T[index])
                    set_value(DIN, change_DIN[n])
                    set_value(Ea, change_Ea[v])
                    set_value(ctag_max, change_ctag_max[c])
                    # Optimize the model for the set value
                    optimize!(model)
                    # Set variable to store the model value for each iteration
                    model_val = 0
                    # Calculate the corrected model value based on variable/NL expression identity
                    for var in 1:length(vars)
                        if var == 1
                            model_val =  exp.(JuMP.value.(vars[var]))
                        elseif var < change_val
                            model_val = JuMP.value.(vars[var])
                        else
                            model_val = value(vars[var])
                        end
                        # Replace the old model value with the new reoptimized value
                        reopt[vs[var]][index, n, v, c] = model_val
                    end
                end
            end
        end
    end
end
