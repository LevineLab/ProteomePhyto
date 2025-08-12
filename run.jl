#-----------------------------------------------------# 
# Set parameters and environment and run the model    #
#-----------------------------------------------------#

# July 28, 2025; by Suzana Goncalves Leles

# Import packages
using CSV
using DataFrames

# STEP 1: Define the phytoplankton functional type (or sensitivity analysis "sens")
type = "pro"                        # takes one of 5: "pro" or "syn" or "warm_diatom" or "cold_diatom" or "sens"
include("par_$(type).jl")              

# STEP 2: Define run mode
run_mode = "thermal_curve"          # takes one of 2: "thermal_curve" or "latitudinal"

# STEP 3: if running in latitudinal mode, define also the time period
time_period = "future"              # takes one of 2: "present" or "future"

# Set the environmenal parameters, run, and save model ouput
if run_mode == "thermal_curve"

    change_DIN = [20.0, 0.1]                  # uM, dissolved inorganic nitrogen
    change_I = 300.0                          # μmol photon m-2 s-1, irradiance

    if type == "sens"
        change_T = 273.0 .+ [18.0:0.25:35;]   # K, temperature   
    else
        change_T = 273.0 .+ [5.0:1.0:36;]     # K, temperature
    end

    output_filename = "output_$(type)_$(run_mode).csv"
 
    include("eqn.jl")                         # defines model equations
    include("loop.jl")                        # iterates over different conditions
    include("reopt.jl")                       # corrects for non-optimal solutions

    # Reshape variables and save model output as a csv file
    vars = ["Qcell", "vru", "vres", "tot_p", "ϕld", "ϕgl", "ϕp", "ϕre", "ϕtr", "ϕlb", "ϕri", "ϕru", "clm", "ptr", "pp", "pdp", "pri", "cli", "cgu", "ctag", "vtr", "Qc", "Qn", "cic", "cin", "vgl", "DIN", "I", "vri", "vlb", "vd", "vre", "D", "ϵ", "pgl", "plb", "pld", "pre", "Ea", "ctag_max", "αtag"]
    re = Dict()
    
    a = length(change_T)
    b = length(change_DIN)
    c = length(change_Ea)
    d = length(change_ctag_max)

    for i in 1:length(vars)
        re[vars[i]] = reshape(s[vars[i]], (a * b * c * d, 1))
    end

    re["T"] = reshape(s["T"] .- 273.0, (a * b * c * d, 1))
    re["μ"] = reshape(s["μ"] *60*24, (a * b * c * d, 1))
    re["β"] = reshape(s["β"]*3*2, (a * b * c * d, 1))
    re["chlvol"] = reshape(s["pp"] + s["pdp"], (a * b * c * d, 1))
    re["ctot"] = reshape(s["cli"] * ηlic .+ s["cgu"] * ηguc .+ s["ctag"] * ηlic, (a * b * c * d, 1))

    for var in vars
        re[var] = re[var][:]
    end

    output =  DataFrame("temp" => re["T"][:], "mu" => re["μ"][:], "size" => re["β"][:], "NC" => re["Qcell"][:], "Ea" => re["Ea"], "ctag_max" => re["ctag_max"], "phot" => re["vru"][:], "resp" => re["vres"][:], "totp" => re["tot_p"][:], 
                    "chlvol" => re["chlvol"][:], "ilpd" => re["ϕld"][:], "icbd" => re["ϕgl"][:], "ichl" => re["ϕp"][:], "icha" => re["ϕre"][:], "itr" => re["ϕtr"][:],
                    "ilpb" => re["ϕlb"][:], "irib" => re["ϕri"][:], "irub" => re["ϕru"][:], "clip" => re["clm"][:], "ctr" => re["ptr"][:], "cp" => re["pp"][:], "cup" => re["pdp"][:], 
                    "cri" => re["pri"][:], "cli" => re["cli"][:], "cgu" => re["cgu"][:], "ctag" => re["ctag"][:], "ctot" => re["ctot"][:], "vtr" => re["vtr"][:], 
                    "Qc" => re["Qc"][:], "Qn" => re["Qn"][:], "cic" => re["cic"][:], "vgl" => re["vgl"][:], "DIN" => re["DIN"][:], "light" => re["I"][:], "vlb" => re["vlb"][:], "D" => re["D"][:], 
                    "esp" => re["ϵ"][:], "cin" => re["cin"][:], "pgl" => re["pgl"][:], "pld" => re["pld"][:], "plb" => re["plb"][:], "pre" => re["pre"][:], "vri" =>re["vri"][:], "vd" => re["vd"][:], "itag" => re["αtag"][:])

    CSV.write(output_filename, output; header=true)


elseif run_mode == "latitudinal"

    cesm = CSV.read("model_CESM.csv", DataFrame)                        # import here the climate projections from Earth System Model
    cesm = dropmissing(cesm)
    output_filename = "output_$(type)_$(run_mode)_$(time_period).csv"

    if time_period == "present"
        change_T = cesm[!,"20m_temperature_1980_2000"] .+ 273           # K, temperature 
        change_DIN = cesm[!,"20m_nitrate_1980_2000"] .* 1000            # uM, dissolved inorganic nitrogen
    elseif time_period == "future"
        change_T = cesm[!,"20m_temperature_2080_2100"] .+ 273           # K, temperature
        change_DIN = cesm[!,"20m_nitrate_2080_2100"] .* 1000            # uM, dissolved inorganic nitrogen
    else
        error("Unrecognized mode: $mode")
    end

    include("eqn.jl")                   # defines model equations
    include("loop.jl")                  # iterates over different conditions

    # Reshape variables and save model output as a csv file
    vars = ["Qcell", "vru", "vres", "tot_p", "ϕld", "ϕgl", "ϕp", "ϕre", "ϕtr", "ϕlb", "ϕri", "ϕru", "clm", "ptr", "pp", "pdp", "pri", "cli", "cgu", "ctag", "vtr", "Qc", "Qn", "cic", "cin", "vgl", "DIN", "I", "vri", "vlb", "vd", "vre", "D", "ϵ", "pgl", "plb", "pld", "pre", "Ea", "ctag_max", "αtag"]
    re = Dict()

    a = length(change_T)
    b = length(change_DIN)
    c = length(change_Ea)
    d = length(change_ctag_max)


    for i in 1:length(vars)
        re[vars[i]] = reshape(s[vars[i]], (a * c * d, 1))
    end

    re["T"] = reshape(s["T"] .- 273.0, (a * c * d, 1))
    re["μ"] = reshape(s["μ"] *60*24, (a * c * d, 1))
    re["β"] = reshape(s["β"]*3*2, (a * c * d, 1))
    re["chlvol"] = reshape(s["pp"] + s["pdp"], (a * c * d, 1))
    re["ctot"] = reshape(s["cli"] * ηlic .+ s["cgu"] * ηguc .+ s["ctag"] * ηlic, (a * c * d, 1))

    for var in vars
        re[var] = re[var][:]
    end
    
    output =  DataFrame("temp" => re["T"][:], "mu" => re["μ"][:], "size" => re["β"][:], "NC" => re["Qcell"][:], "Ea" => re["Ea"], "ctag_max" => re["ctag_max"], "phot" => re["vru"][:], "resp" => re["vres"][:], "totp" => re["tot_p"][:], 
                    "chlvol" => re["chlvol"][:], "ilpd" => re["ϕld"][:], "icbd" => re["ϕgl"][:], "ichl" => re["ϕp"][:], "icha" => re["ϕre"][:], "itr" => re["ϕtr"][:],
                    "ilpb" => re["ϕlb"][:], "irib" => re["ϕri"][:], "irub" => re["ϕru"][:], "clip" => re["clm"][:], "ctr" => re["ptr"][:], "cp" => re["pp"][:], "cup" => re["pdp"][:], 
                    "cri" => re["pri"][:], "cli" => re["cli"][:], "cgu" => re["cgu"][:], "ctag" => re["ctag"][:], "ctot" => re["ctot"][:], "vtr" => re["vtr"][:], 
                    "Qc" => re["Qc"][:], "Qn" => re["Qn"][:], "cic" => re["cic"][:], "vgl" => re["vgl"][:], "DIN" => re["DIN"][:], "light" => re["I"][:], "vlb" => re["vlb"][:], "D" => re["D"][:], 
                    "esp" => re["ϵ"][:], "cin" => re["cin"][:], "pgl" => re["pgl"][:], "pld" => re["pld"][:], "plb" => re["plb"][:], "pre" => re["pre"][:], "vri" =>re["vri"][:], "vd" => re["vd"][:], "itag" => re["αtag"][:])
                    
    CSV.write(output_filename, output; header=true)

else 
    error("Unrecognized mode: $mode")
end
