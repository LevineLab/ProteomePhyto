#-------------------#
# Plot model output #
#-------------------#

# March 2023; by Suzana Goncalves Leles

using Plots

gr()

# Convert temperature to degree Celsius
change_T = change_T .- 273.0

a = length(change_T)
b = length(change_DIN)
c = length(change_I)

if env == "Temperature (C°)"
        variable = change_T
        a = 1:length(change_T)
elseif env == "[DIN]"
        variable = change_DIN
        b = 1:length(change_DIN)
elseif env == "Light"
        variable = change_I
        c = 1:length(change_I)
end

plot(layout = (4,4), guidefontsize = 9)

plot!(variable, s_β[a,b,c]*3*2, color = "black", label = "β", ylabel = "ESD (μm)", subplot = 1)
scatter!(variable, s_β[a,b,c]*3*2, color = "black", label = "", ylims = (0.0, 100), subplot = 1, legend=:topleft)

plot!(variable, s_μ[a,b,c]*60*24, color = "black", label = "μ", subplot = 2)
scatter!(variable, s_μ[a,b,c]*60*24, color = "black", label = "", ylabel = "growth rate (1/day)", subplot = 2,
        ylims = (0.0, 1.8), legend=:topleft)

plot!(variable,  e_ru*s_vru[a,b,c]*60*24,  color = "olive", label = "vru", ylabel = "moles/vol/day", subplot = 3)
scatter!(variable,  e_ru*s_vru[a,b,c]*60*24, color = "olive", label = "", subplot = 3)
plot!(variable,  e_gl*s_vgl[a,b,c]*60*24 + e_ld*s_vld[a,b,c]*60*24,  color = "black", label = "vres", subplot = 3)
scatter!(variable,  e_gl*s_vgl[a,b,c]*60*24 + e_ld*s_vld[a,b,c]*60*24, color = "black", label = "", subplot = 3)
plot!(variable,  e_p*s_vp[a,b,c]*60*24,  color = "green", label = "vp", subplot = 3)
scatter!(variable,  e_p*s_vp[a,b,c]*60*24, color = "green", label = "", subplot = 3, legend=:topleft)

plot!(variable, s_cin[a,b,c], color = "blue", label = "cin", ylabel = "moles/vol", subplot = 4)
scatter!(variable, s_cin[a,b,c], color = "blue", label = "", subplot = 4)
plot!(variable, s_cic[a,b,c], color = "green", label = "cic", subplot = 4)
scatter!(variable, s_cic[a,b,c], color = "green", label = "", subplot = 4)

plot!(variable, (s_pp[a,b,c]), color = "green",  label = "pp", ylabel = "moles/vol", subplot = 5)
scatter!(variable,  (s_pp[a,b,c]), color = "green", label = "", subplot = 5, legend=:topleft)
plot!(variable, s_pru[a,b,c], color = "olive",  label = "pru", subplot = 5)
scatter!(variable,  s_pru[a,b,c], color = "olive", label = "", subplot = 5)
plot!(variable, s_pld[a,b,c], color = "deeppink",  label = "pld", subplot = 5)
scatter!(variable,  s_pld[a,b,c], color = "deeppink", label = "", subplot = 5)
plot!(variable, s_pgl[a,b,c], color = "purple",  label = "pgl", subplot = 5)
scatter!(variable,  s_pgl[a,b,c], color = "purple", label = "", subplot = 5)
plot!(variable, s_pre[a,b,c], color = "yellow",  label = "pre", subplot = 5)
scatter!(variable,  s_pre[a,b,c], color = "yellow", label = "", subplot = 5)
plot!(variable, s_pri[a,b,c], color = "red",  label = "pri", subplot = 5)
scatter!(variable,  s_pri[a,b,c], color = "red", label = "", subplot = 5)

plot!(variable, s_plb[a,b,c], color = "orange", label = "plb", subplot = 6, legend=:topleft)
scatter!(variable,  s_plb[a,b,c], color = "orange", label = "", ylabel = "moles/vol", subplot = 6)

plot!(variable,  s_clm[a,b,c],  color = "orange", label = "mem_lip", subplot = 7)
scatter!(variable,  s_clm[a,b,c], color = "orange", label = "", ylabel = "lipids", subplot = 7)
plot!(variable, s_ptr[a,b,c], color = "blue", label = "ptr", subplot = 7, legend=:topright)
scatter!(variable, s_ptr[a,b,c], color = "blue", label = "", subplot = 7)

plot!(variable,  s_ϕp[a,b,c],  color = "green", label = "ϕp", ylabel = "allocation factors", subplot = 8)
scatter!(variable,  s_ϕp[a,b,c], color = "green", label = "", xlabel = env, subplot = 8)
plot!(variable, s_ϕru[a,b,c], color = "olive",   label = "ϕru", subplot = 8, legend=:topleft)
scatter!(variable,  s_ϕru[a,b,c], color = "olive", label = "", subplot = 8, ylims = (0.0, 1.0))
plot!(variable, s_ϕre[a,b,c], color = "yellow",   label = "ϕre", subplot = 8)
scatter!(variable,  s_ϕre[a,b,c], color = "yellow", label = "", subplot = 8)
plot!(variable, s_ϕri[a,b,c], color = "red",  label = "ϕri", subplot = 8)
scatter!(variable,  s_ϕri[a,b,c], color = "red", label = "", subplot = 8)

plot!(variable,  s_ctag[a,b,c],  color = "blue", label = "ctag", ylabel = "allocation factors", subplot = 9)
scatter!(variable,  s_ctag[a,b,c], color = "blue", ylims = (0.0, 20), label = "", xlabel = env, subplot = 9)

plot!(variable,  s_ϵ[a,b,c],  color = "black", label = "ϵ", ylabel = "Vother (%)", subplot = 10)
scatter!(variable,  s_ϵ[a,b,c], color = "black", label = "", ylims = (0.0, 1), xlabel = env, subplot = 10)

plot!(variable,  s_ϕld[a,b,c],  color = "deeppink", label = "ϕld", subplot = 11)
scatter!(variable,  s_ϕld[a,b,c], color = "deeppink", label = "", ylabel = "relative investment", subplot = 11)
plot!(variable, s_ϕgl[a,b,c], color = "purple", label = "ϕgl", subplot = 11)
scatter!(variable, s_ϕgl[a,b,c], color = "purple", label = "", subplot = 11, legend=:topright)
plot!(variable,  s_αld[a,b,c],  color = "orange", label = "αld", ylims = (0.0, 1.0), subplot = 11)
scatter!(variable,  s_αld[a,b,c], color = "orange", label = "", ylabel = "allocation factor", subplot = 11)

plot!(variable, s_Qcell[a,b,c], color = "pink",  label = "NC", ylims = (0.005,0.35), subplot = 12)
scatter!(variable, s_Qcell[a,b,c], color = "pink", label = "", subplot = 12, legend=:topright)

plot!(variable, (s_pp[a,b,c]), color = "green",  label = "folded chl", ylabel = "moles/vol", subplot = 13)
scatter!(variable,  s_pp[a,b,c], color = "green", label = "", subplot = 13, legend=:topleft)
plot!(variable, (s_pdp[a,b,c]), color = "red",  label = "unfolded chl", subplot = 13)
scatter!(variable,  s_pdp[a,b,c], color = "red", label = "", subplot = 13)

plot!(variable, s_cgu[a,b,c], color = "purple", label = "starch", subplot = 14)
scatter!(variable, s_cgu[a,b,c], color = "purple", label = "", subplot = 14)
plot!(variable, s_cli[a,b,c], color = "deeppink", label = "TAGs_resp", subplot = 14)
scatter!(variable, s_cli[a,b,c], color = "deeppink", label = "", subplot = 14)
plot!(variable, s_ctag[a,b,c], color = "pink", label = "TAGs_ROS", subplot = 14)
scatter!(variable, s_ctag[a,b,c], color = "pink", label = "", subplot = 14)

nATP = s_vri[a,b,c]*e_ri + s_vtr[a,b,c]*e_tr + s_vlb[a,b,c]*e_lb
cue = 1 .- s_vres[a,b,c]./s_vru[a,b,c]

plot!(variable, cue, color = "black", label = "CUE", ylabel = "carbon use efficiency", subplot = 15)
scatter!(variable, cue, color = "black", label = "", ylims = (0,1), subplot = 15)

plot!(variable, s_αld[a,b,c], color = "deeppink",  label = "αld", ylabel = "fraction vlb", subplot = 16)
scatter!(variable, s_αld[a,b,c], color = "deeppink", label = "", subplot = 16, legend=:topleft)
plot!(variable, s_αlm[a,b,c], color = "blue",  label = "αlm", ylabel = "", subplot = 16)
scatter!(variable, s_αlm[a,b,c], color = "blue", label = "", subplot = 16)
plot!(variable, s_αtag[a,b,c], color = "pink",  label = "αtag", ylabel = "", subplot = 16)
scatter!(variable, s_αtag[a,b,c], color = "pink", ylims = (0,1), label = "", subplot = 16)

plot!(size=(1000,900))

if env == "Temperature (C°)"
        savefig("./model_output.pdf")
elseif env == "[DIN]"
        savefig("./din.pdf")
elseif env == "Light"
        savefig("./light.pdf")
end
