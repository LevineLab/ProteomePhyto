
include("par.jl")
include("eqn.jl")
include("save.jl")
include("reopt.jl")

# Reshape Variables for Plotting and Save to CSV
using CSV
using DataFrames

vars = ["Qcell", "vru", "vres", "tot_p", "ϕld", "ϕgl", "ϕp", "ϕre", "ϕtr", "ϕlb", "ϕri", "ϕru", "clm", "ptr", "pp", "pdp", "pri", "cli", "cgu", "ctag", "vtr", "Qc", "Qn", "cic", "cin", "vgl", "DIN", "I", "vri", "vlb", "vd", "vre", "D", "ϵ", "pgl", "plb", "pld", "pre", "Vnu", "ctag_max", "αtag"]
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

output =  DataFrame("temp" => re["T"][:], "mu" => re["μ"][:], "size" => re["β"][:], "NC" => re["Qcell"][:], "Vnu" => re["Vnu"], "ctag_max" => re["ctag_max"], "phot" => re["vru"][:], "resp" => re["vres"][:], "totp" => re["tot_p"][:], 
                    "chlvol" => re["chlvol"][:], "ilpd" => re["ϕld"][:], "icbd" => re["ϕgl"][:], "ichl" => re["ϕp"][:], "icha" => re["ϕre"][:], "itr" => re["ϕtr"][:],
                    "ilpb" => re["ϕlb"][:], "irib" => re["ϕri"][:], "irub" => re["ϕru"][:], "clip" => re["clm"][:], "ctr" => re["ptr"][:], "cp" => re["pp"][:], "cup" => re["pdp"][:], 
                    "cri" => re["pri"][:], "cli" => re["cli"][:], "cgu" => re["cgu"][:], "ctag" => re["ctag"][:], "ctot" => re["ctot"][:], "vtr" => re["vtr"][:], 
                    "Qc" => re["Qc"][:], "Qn" => re["Qn"][:], "cic" => re["cic"][:], "vgl" => re["vgl"][:], "DIN" => re["DIN"][:], "light" => re["I"][:], "vlb" => re["vlb"][:], "D" => re["D"][:], 
                    "esp" => re["ϵ"][:], "cin" => re["cin"][:], "pgl" => re["pgl"][:], "pld" => re["pld"][:], "plb" => re["plb"][:], "pre" => re["pre"][:], "vri" =>re["vri"][:], "vd" => re["vd"][:], "itag" => re["αtag"][:])

CSV.write("output.csv", output, header = true)
