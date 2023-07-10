#---------------------------------------------------------------------------#
# Equations defining the proteome allocation model for a phytoplankton cell #
#---------------------------------------------------------------------------#

# March 2023; by Suzana Goncalves Leles

using JuMP, Ipopt

# Create a JuMP model, using Ipopt as the solver.
model = Model(Ipopt.Optimizer)

# Set the environment:
@NLparameter(model, T   == 20.0 + 273.0)   # K,  temperature
@NLparameter(model, DIN == 20.0)           # uM, external dissolved inorganic nitrogen
@NLparameter(model, DIC == 100.0)          # uM, external dissolved inorganic carbon
@NLparameter(model, I   == 120.0)          # μmol photon m-2 s-1, irradiance

# Define the unknown variables that we will optimize for.
@variables(model, begin
        logμ   , (start = 0.00047)  # 1/minute, log growth rate
        ptr ≥ 0, (start = 369)      # protein: number of  transporters / um3
        pri ≥ 0, (start = 703)      # protein: number of ribosomes / um3
        plb ≥ 0, (start = 337)      # protein: number of lipid biosynthesis pathway proteins / um3
        pp  ≥ 0, (start = 1585)     # protein: number of photosystems / um3
        pdp ≥ 0, (start = 228)      # protein: number of damaged photosystems / um3
        pru ≥ 0, (start = 4241)     # protein: number of rubisco / um3
        pld ≥ 0, (start = 0)        # protein: number of lipid degradation pathway proteins / um3
        pgl ≥ 0, (start = 1840)     # protein: number of glycolysis pathway prteins / um3
        pre ≥ 0, (start = 358)      # protein: number of repair proteins / um3
        cin ≥ 0, (start = 2e6)      # other macromolecule: number of internal nitrogen molecules / um3
        clm ≥ 0, (start = 1.17e6)   # other macromolecule: number of membrane lipids / um3
        cic ≥ 0, (start = 2.4e6)    # other macromolecule: number of internal carbon molecules / um3
        ctag≥ 0, (start = 4.6e6)    # other macromolecule: number of lipid droplets / um3
        ϕtr ≥ 0, (start = 0.002)    # relative proteome investment in: transporters (unitless)
        ϕri ≥ 0, (start = 0.07)     # relative proteome investment in: ribosomes (unitless)
        ϕlb ≥ 0, (start = 0.02)     # relative proteome investment in: lipid biosynthesis (unitless)
        ϕp  ≥ 0, (start = 0.09)     # relative proteome investment in: photosystems (unitless)
        ϕru ≥ 0, (start = 0.14)     # relative proteome investment in: rubisco (unitless)
        ϕld ≥ 0, (start = 0.00)     # relative proteome investment in: lipid degradation (unitless)
        ϕgl ≥ 0, (start = 0.15)     # relative proteome investment in: glycolysis (unitless)
        ϕre ≥ 0, (start = 0.01)     # relative proteome investment in: repair (unitless)
        β   ≥ 0, (start = 1.7)      # um, cellular volume to surface area ratio
        αlm ≥ 0, (start = 0.2)      # fraction of lipid synthesis flux used to: synthesize membrane lipids (unitless)
        αld ≥ 0, (start = 0.0)      # fraction of lipid synthesis flux used to: fuel dark respiration (unitless)
        αtag≥ 0, (start = 0.6)      # fraction of lipid synthesis flux used to: synthesize lipid storage (unitless)
end)

# Set the objective: maximize the log growth rate (logμ).
@objective(model, Max, logμ)

# Define any intermediate calculations.
@NLexpressions(model, begin
        γTa,   exp(Ea*(1.0/(R*Tref) - 1.0/(R*T)))                                              # unitless
        γTd,   exp(Ed*(1.0/(R*Td) - 1.0/(R*T)))                                                # unitless
        γTad,  exp(Etr*(1.0/(R*Tref) - 1.0/(R*T)))                                             # unitless
        γTm,   (Tmax - T)/(Tmax - Tmin)                                                        # unitless
        vp,  ( kref["p"] ) * pp * I / (K["p"] + I)                                             # molecules of photons/μm3/min
        vru, ( kref["ru"] * γTa ) * pru * DIC / ( K["ru"] + DIC )                              # molecules of carbon/μm3/min
        vtr, ( kref["tr"] * γTad ) * ptr * DIN/ ( K["tr"] + DIN)                               # molecules of nitrogen/μm3/min
        vri, ( kref["ri"] * γTa ) * pri * (cic / ( K["ri"] + cic )) * (cin / ( K["ri"] + cin)) # molecules of amino acids/μm3/min
        vlb, ( kref["lb"] * γTa ) * plb * cic / ( K["lb"] + cic )                              # molecules of carbon/μm3/min
        vld, ( kref["ld"] * γTa ) * pld                                                        # molecules of carbon/μm3/min
        vd,  ( kref["d"]  * γTd ) * pp * (ctag_max/ctag)                                       # molecules of photosystems/μm3/min
        vre,  ( kref["re"]  * γTa ) * pre * (pdp/(pdp + pp))                                   # molecules of photosystems/μm3/min
        #vre, ( kref["re"] * γTa ) * pre * (pdp/(pdp + K["f"]))                                 # molecules of photosystems/μm3/min; alternative approach
        vgl,   ( kref["gl"] * γTa ) * pgl                                                      # molecules of carbon/μm3/min
        vres,  vld + vgl                                                                       # molecules of carbon/μm3/min
        cgu,  vgl / ηguc * td                                                                  # molecules of glucose/μm3 (required to fuel dark respiration)
        cli,  vld / ηlic * td                                                                  # molecules of lipid storage/μm3 (required to fuel dark respiration)
        ηaan,  ηaac * ((ϕtr + ϕru + ϕlb + ϕld + ϕgl + ϕre)*Qpt + ϕri*Qri + ϕp*Qp)              # molecules of nitrogen / molecule of amino acid
        Qn,   (ptr*ηtr*Qpt + pru*ηru*Qpt + plb*ηlb*Qpt + pld*ηld*Qpt + pgl*ηgl*Qpt + 
               pre*ηre*Qpt + pri*ηri*Qri + (pp + pdp)*ηp*Qp + cot*ηot*Qpt)*ηaac/ηaa + cin      # total number of molecules of nitrogen/um3
        Qc,   (ptr*ηtr + pru*ηru + plb*ηlb + pld*ηld + pgl*ηgl + pre*ηre + (pp + pdp)*ηp + 
              pri*ηri + cot*ηot)*ηaac/ηaa + cic + clm*ηlic + ctag*ηlic + cli*ηlic + cgu*ηguc   # total number of molecules of carbon
        Qcell, Qn/Qc                                                                           # nitrogen to carbon quota of the cell
        vol,     4/3*π*(β*3)^3                                                                 # um3; cell biovolume
        ϵ,       1.0 - (pp + pdp) * Vp - cli * Vli - ctag * Vli                                # unitless; fraction of the biovolume not occupied by photosystems nor lipid storage  
        tot_p,   (pp + pdp)*ηp + pru*ηru + ptr*ηtr + pri*ηri + plb*ηlb + pld*ηld + pgl*ηgl + 
                  pre*ηre + cot*ηot                                                            # Da/μm3, total protein content
        ecost, vtr*e_tr + vri*e_ri + vlb*e_lb + vd*e_u + vre*e_re                              # molecules of energy/μm3/min
        D, ηru*pru + ηri*pri + ηlb*plb + ηld*pld +
           ηgl*pgl + ηre*pre + ηin*cin + ηic*cic + cgu*ηgu + cot*ηot                           # Da/μm3, concentration of all macromolecules except for photosystems, lipid storage and membrane components
end)

# Define the system of equations and constraints.
@constraints(model, begin
        ϕp + ϕru + ϕtr + ϕri + ϕlb + ϕld + ϕgl + ϕre ==  0.5
        αlm + αld + αtag == 1.0
        ϕld == 0
end)

@NLconstraints(model, begin
        ptr / clm ≤ Mmin + (Mmax * γTm)                                          # constraint on the trasnporter to lipid ratio in the membrane
        β * (s_tr * ptr + s_lm * clm) == 1.0                                     # constraint on the volume to surface area ratio of the cell
        pdp / (pdp + pp) == ( ( kref["d"]  * γTd ) * pp * (ctag_max/ctag) ) / ( kref["re"] * γTa * pre ) # damage (vd) and repair (vre) dynamics assuming vd = vre
        ϕp * vri * 1/(ηp/ηaa) - exp(logμ) * (pp + pdp) == 0.0                    # protein: molecules of photosystems/μm3
        # alternative approach to model folding/unfolding of photosystems:
        #ϕp * vri * 1/(ηp/ηaa) + vre - vd - exp(logμ) * pp == 0.0                  # protein: functional photosystems/μm3
        #vd - vre - exp(logμ) * pdp == 0.0                                         # protein: damaged photosystems/μm3
        ϕru * vri * 1/(ηru/ηaa) - exp(logμ) * pru == 0.0                         # protein: rubisco/μm3
        ϕtr  * vri * 1/(ηtr/ηaa)  - exp(logμ) * ptr  == 0.0                      # protein: transporters/μm3
        ϕri * vri * 1/(ηri/ηaa) - exp(logμ) * pri == 0.0                         # protein: ribosomes/μm3
        ϕlb * vri * 1/(ηlb/ηaa) - exp(logμ) * plb == 0.0                         # protein: plb/μm3
        ϕld * vri * 1/(ηld/ηaa) - exp(logμ) * pld == 0.0                         # protein: pld/μm3
        ϕgl * vri * 1/(ηgl/ηaa) - exp(logμ) * pgl == 0.0                         # protein: pgl/μm3
        ϕre * vri * 1/(ηre/ηaa) - exp(logμ) * pre == 0.0                         # protein: pre/μm3
        vtr - vri*ηaan - exp(logμ) * cin == 0.0                                  # internal pool of nitrogen/μm3
        vru - vri*ηaac - vlb - vgl - exp(logμ) * cic == 0.0                      # internal pool of carbon/μm3
        αlm * vlb/ηlic - exp(logμ) * clm == 0.0                                  # membrane lipids/μm3 (evolved)
        αtag * vlb/ηlic - exp(logμ) * ctag == 0.0                                # lipid storage/μm3 (used to regulate heat-damage; not used in dark respiration)
        ctag ≤ ctag_min + ctag_max                                               # constraint on the maximum concentration of ctag; lipid storage/μm3
        vld == αld * vlb                                                         # lipid degradation rate; molecules of carbon/μm3/min
        # energetic requirements:
        vp * e_p ≥ vru * e_ru                                                    # rubisco costs are paid by photosystems; molecules of ATP/μm3/min
        ecost ≤ (vp*e_p - vru*e_ru)                                              # total ATP cost; molecules of ATP/μm3/min
        ecost*fdr ≤ vgl*e_gl + vld*e_ld                                          # dark respiration; molecules of ATP/μm3/min
        1/ϵ * D ≤ Dmax                                                           # density constraint; Da/μm3
end)

# Solve the optimization problem.
status = optimize!(model)
