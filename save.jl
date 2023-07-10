#-----------------------------------------------------------------------#
# Iterate over different environmental conditions and save model output #
#-----------------------------------------------------------------------#

# March 2023; by Suzana Goncalves Leles

# Optimal solutions looping over different environmental conditions.
env = "Temperature (C°)" # should be one of the following: Light, Temperature (C°) or [DIN]

if env == "[DIN]"
    # Run varying [DIN].
    change_DIN = 10.0 .^ [-7.0:0.18:log10(20.0);]
    change_T = 20.0 + 273.0
    change_I = 20.0
elseif env == "Temperature (C°)"
    # Run varying Temperature.
    change_DIN = 20.0
    change_I = 20.0
    change_T  = 273.0 .+ [3.0:1.0:45.0;]
elseif env == "Light"
    # Run varying Light.
    change_DIN = 20.0
    change_T = 20.0 + 273.0
    change_I = 10.0 .^ [-7.0:0.18:log10(20.0);]
end

# Define empty vectors
s_μ     = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_β     = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ptr   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pri   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_plb   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pdp   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pp    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pru   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pld   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pgl   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_pre   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕtr   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕri   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕlb   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕp    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕru   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕld   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕgl   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϕre   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cin   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_clm   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cic   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_Qcell = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_Qc    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_Qn    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vru   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vres  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vp    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vgl   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vd    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vre    = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vlb   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vld   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vri   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vtr   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_vol   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_totp  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_αld   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_D     = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTa   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTd   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTad  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_γTm   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cli   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_cgu   = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ecost = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ηaan  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ϵ     = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_ctag = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_αlm  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_αtag  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_T  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_DIN  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))
s_I  = Array{Float64,3}(undef,length(change_T),length(change_DIN),length(change_I))

# Store model output in empty vectors
for i in 1:length(change_DIN)
        set_value(DIN, change_DIN[i])
        for j in 1:length(change_T)
                set_value(T, change_T[j])
                for k in 1:length(change_I)
                        set_value(I, change_I[k])
                        optimize!(model)
                        s_μ[j,i,k]   = exp.(JuMP.value.(logμ))
                        s_β[j,i,k]   = JuMP.value.(β)
                        s_ptr[j,i,k] = JuMP.value.(ptr)
                        s_pri[j,i,k] = JuMP.value.(pri)
                        s_plb[j,i,k] = JuMP.value.(plb)
                        s_pdp[j,i,k] = JuMP.value.(pdp)
                        s_pp[j,i,k]  = JuMP.value.(pp)
                        s_pru[j,i,k] = JuMP.value.(pru)
                        s_pld[j,i,k] = JuMP.value.(pld)
                        s_pgl[j,i,k] = JuMP.value.(pgl)
                        s_pre[j,i,k] = JuMP.value.(pre)
                        s_ϕtr[j,i,k] = JuMP.value.(ϕtr)
                        s_ϕri[j,i,k] = JuMP.value.(ϕri)
                        s_ϕlb[j,i,k] = JuMP.value.(ϕlb)
                        s_ϕp[j,i,k]  = JuMP.value.(ϕp)
                        s_ϕru[j,i,k] = JuMP.value.(ϕru)
                        s_ϕld[j,i,k] = JuMP.value.(ϕld)
                        s_ϕgl[j,i,k] = JuMP.value.(ϕgl)
                        s_ϕre[j,i,k] = JuMP.value.(ϕre)
                        s_cin[j,i,k] = JuMP.value.(cin)
                        s_clm[j,i,k] = JuMP.value.(clm)
                        s_cic[j,i,k]  = JuMP.value.(cic)
                        s_Qcell[j,i,k] = value(Qcell)
                        s_Qc[j,i,k] = value(Qc)
                        s_Qn[j,i,k] = value(Qn)
                        s_vru[j,i,k]  = value(vru)
                        s_vres[j,i,k] = value(vres)
                        s_vp[j,i,k]   = value(vp)
                        s_vgl[j,i,k]  = value(vgl)
                        s_vd[j,i,k]   = value(vd)
                        s_vre[j,i,k]   = value(vre)
                        s_vlb[j,i,k]  = value(vlb)
                        s_vld[j,i,k]  = value(vld)
                        s_vri[j,i,k]  = value(vri)
                        s_vtr[j,i,k]  = value(vtr)
                        s_vol[j,i,k]  = value(vol)
                        s_totp[j,i,k] = value(tot_p)
                        s_αld[j,i,k]  = JuMP.value.(αld)
                        s_D[j,i,k]    = value(D)
                        s_γTa[j,i,k]  = value(γTa)
                        s_γTd[j,i,k]  = value(γTd)
                        s_γTad[j,i,k] = value(γTad)
                        s_γTm[j,i,k]  = value(γTm)
                        s_cli[j,i,k]  = value(cli)
                        s_cgu[j,i,k]  = value(cgu)
                        s_ecost[j,i,k] = value(ecost)
                        s_ηaan[j,i,k]  = value(ηaan)
                        s_ϵ[j,i,k]     = value(ϵ)
                        s_ctag[j,i,k]  = JuMP.value.(ctag)
                        s_αlm[j,i,k]  = JuMP.value.(αlm)
                        s_αtag[j,i,k]  = JuMP.value.(αtag)
                        s_T[j,i,k] = value(T)
                        s_DIN[j,i,k] = value(DIN)
                        s_I[j,i,k] = value(I)
                end
        end
end

using CSV
using DataFrames

a = length(change_T)
b = length(change_DIN)
c = length(change_I)

# Reshape the array into a 1-dimensional array
temp = reshape(s_T .- 273.0, (a * b * c, 1))
mu = reshape(s_μ *60*24, (a * b * c, 1))
size = reshape(s_β*3*2, (a * b * c, 1))
NC = reshape(s_Qcell, (a * b * c, 1))
phot = reshape(s_vru, (a * b * c, 1))
resp = reshape(s_vres, (a * b * c, 1))
totp = reshape(s_totp, (a * b * c, 1))
chlvol = reshape(s_pp + s_pdp, (a * b * c, 1))
ilpd = reshape(s_ϕld, (a * b * c, 1))
icbd = reshape(s_ϕgl, (a * b * c, 1))
ichl = reshape(s_ϕp, (a * b * c, 1))
icha = reshape(s_ϕre, (a * b * c, 1))
itr = reshape(s_ϕtr, (a * b * c, 1))
ilpb = reshape(s_ϕlb, (a * b * c, 1))
irib = reshape(s_ϕri, (a * b * c, 1))
irub = reshape(s_ϕru, (a * b * c, 1))
clip = reshape(s_clm, (a * b * c, 1))
ctr = reshape(s_ptr, (a * b * c, 1))
cp = reshape(s_pp, (a * b * c, 1))
cup = reshape(s_pdp, (a * b * c, 1))
cri = reshape(s_pri, (a * b * c, 1))
cli = reshape(s_cli, (a * b * c, 1))
cgu = reshape(s_cgu, (a * b * c, 1))
ctag = reshape(s_ctag, (a * b * c, 1))
ctot = reshape(s_cli * ηlic .+ s_cgu * ηguc .+ s_ctag * ηlic, (a * b * c, 1))
vtr = reshape(s_vtr, (a * b * c, 1))
Qc = reshape(s_Qc, (a * b * c, 1))
Qn = reshape(s_Qn, (a * b * c, 1))
cic = reshape(s_cic, (a * b * c, 1))
vgl = reshape(s_vgl, (a * b * c, 1))
din = reshape(s_DIN, (a * b * c, 1))
light = reshape(s_I, (a * b * c, 1))

output =  DataFrame("temp" => temp[:], "mu" => mu[:], "size" => size[:], "NC" => NC[:], "phot" => phot[:], "resp" => resp[:], "totp" => totp[:], 
                    "chlvol" => chlvol[:], "ilpd" => ilpd[:], "icbd" => icbd[:], "ichl" => ichl[:], "icha" => icha[:], "itr" => itr[:],
                    "ilpb" => ilpb[:], "irib" => irib[:], "irub" => irub[:], "clip" => clip[:], "ctr" => ctr[:], "cp" => cp[:], "cup" => cup[:], 
                    "cri" => cri[:], "cli" => cli[:], "cgu" => cgu[:], "ctag" => ctag[:], "ctot" => ctot[:], "vtr" => vtr[:], 
                    "Qc" => Qc[:], "Qn" => Qn[:], "cic" => cic[:], "vgl" => vgl[:], "DIN" => din[:], "light" => light[:])

CSV.write("model_output.csv", output, header = true)
