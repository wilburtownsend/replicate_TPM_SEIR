#==============================================================================
This script replicates the TPM SEIR vaccine model.

Tested on Julia 1.6.3.
Dependencies:
  CSV v0.9.5
  DataFrames v1.2.2
  XLSX v0.7.8

==============================================================================#
cd("/Users/wtownsend/Documents/GitHub/replicate_TPM_SEIR/")
using LinearAlgebra, XLSX, DataFrames, CSV


#==============================================================================
Set various parameters that will be constant across simulations.
==============================================================================#

# The size of the full population.
N_all = 5000000
# Age-specific populations
age_dist = [5.98, 6.39, 6.56, 6.17, 6.59, 7.40,
            7.44, 6.62, 6.08, 6.41, 6.43, 6.38,
            5.77, 4.90, 4.24, 6.64]./100
N = age_dist.*N_all
# Contact matrix -- from Prem et al
C_raw = XLSX.readxlsx("MUestimates_all_locations_2.xlsx")["New Zealand"][:]
# this makes it balanced.
C = [(C_raw[i,j] + (N[j]/N[i])*C_raw[j,i])/2 for i ∈ 1:16, j ∈ 1:16]
# Latent period
tE = 2.55
# Infectious period
tI = 5
# Average hospital length of stay
tH = 8
# The relative infectiousness of subclinical individuals
τ = 0.5
# Proportion of infections causing clinical disease by age group
pclin = [0.544, 0.555, 0.577, 0.5985, 0.6195, 0.6395, 0.6585,
         0.6770, 0.6950, 0.7118, 0.7273, 0.7418, 0.7553, 0.768,
         0.78, 0.8008]
# Hospitalisation probability by age group (note I slightly inflate the first
# entry to avoid 0/0) (the 2.26 factor is just from the latest paper)
phosp = [1.0e-10, 0.0001, 0.0003, 0.0029, 0.0079, 0.0164, 0.0283, 0.0364,
         0.0405, 0.0523, 0.0718, 0.0907, 0.1089, 0.13, 0.154, 0.178]*2.26
# Hospitalisation probability by age group
pdeath = [0, 0.00003, 0.00006, 0.0001, 0.0002, 0.0004, 0.0007, 0.001,
          0.0014, 0.0027, 0.0049, 0.0093, 0.016, 0.0252, 0.0369, 0.0664]*2.26
# Relative susceptibility by age group (compared to 60-64 year-olds)
u = [0.462, 0.457, 0.445, 0.558, 0.795, 0.934, 0.974, 0.977,
     0.942, 0.931, 0.942, 0.965, 1.00, 0.977, 0.896, 0.856]
# Imported cases per age group per day.
m = 1.0 .* age_dist


#==============================================================================
Define a function to simulate an SEIR model.

Note we store the state in a dictionary with keys (θ,v,i), where 
   θ ∈ ["S", "E", "I", "A", "R", "Imm", "H", "F", "Disch"] indexes compartments;
   v ∈ [0, 1] indexes vax status; and
   i indexes ages.

This function returns a dict of states, indexed by day.
==============================================================================#

function simulate_SEIR(;
    num_days = 730, # Number of days to simulate.
    eI = 0.7, # The efficacy of the vaccine against infection.
    eT = 0.5, # The efficacy of the vaccine against transmission.
    eD = 0.8, # The efficacy of the vaccine against severe disease.
    R0 = 6.0*(1-0.17)*(1-0.20), # The basic reproduction number.
    v = [[0., 0., 0.9*(3/5)]; 0.9*ones(13)] # The vaccination rate for each age
                                            # group.
    )

    # Form the next generation matrix. This is only used to find U, given some R0,
    # so we condition on it on U.
    NGM(U) = [U*u[i]*tI*C[j,i]*(pclin[j] + τ*(1-pclin[j])) for i=1:16, j=1:16] 
    # We set U to match the basic reproduction number (searching over a
    # log-equidistant space).
    R0_given_U(U) = maximum(eigvals(NGM(U)))
    U_list =  exp.(range(log(0.01), log(1.), length=100))
    U = U_list[argmin([(R0 - R0_given_U(U))^2 for U ∈ U_list])]

    # Define infection force λ.
    λ(i, state) = U*(u[i]/N[i]) * sum([(           (state[("I",0,j)] + τ*state[("A",0,j)]) 
                                        + (1 - eT)*(state[("I",1,j)] + τ*state[("A",1,j)])
                                        + tI*m[j]
                                        )*C[j,i] for j ∈ 1:16])

    # Define the differential equations:
        # susceptible
    dS(v,i,state) = -λ(i, state)*state[("S",v,i)]
        # exposed
    dE(v,i,state) =  λ(i, state)*state[("S",v,i)] - state[("E",v,i)]/tE
        # clinical infectious
    dI(v,i,state) =  pclin[i]*state[("E",v,i)]/tE - state[("I",v,i)]/tI
        # subclinical infectious
    dA(v,i,state) =  (1-pclin[i])*state[("E",v,i)]/tE - state[("A",v,i)]/tI
        # recovered
    dR(v,i,state) = (state[("I",v,i)] + state[("A",v,i)])/tI
        # hospitalised
    dH(v,i,state) = (phosp[i]/tE)*state[("E",v,i)]*((1-eD)^v) - state[("H",v,i)]/tH
        # deaths
    dF(v,i,state) = (pdeath[i]/phosp[i])*state[("H",v,i)]/tH
        # discharge
    dDisch(v,i,state) = (1 - (pdeath[i]/phosp[i]))*state[("H",v,i)]/tH


    # This function iterates the differential equations.
    function update(state)
        # The vaccinated, immune population is unchanged.
        new_state = Dict(("Imm",1,i) => state[("Imm",1,i)] for i ∈ 1:16)
        # Update other groups.
        for v ∈ 0:1
            for i ∈ 1:16
                new_state[("S",v,i)] = state[("S",v,i)] + dS(v,i,state)
                new_state[("E",v,i)] = state[("E",v,i)] + dE(v,i,state)
                new_state[("I",v,i)] = state[("I",v,i)] + dI(v,i,state)
                new_state[("A",v,i)] = state[("A",v,i)] + dA(v,i,state)
                new_state[("R",v,i)] = state[("R",v,i)] + dR(v,i,state)
                new_state[("H",v,i)] = state[("H",v,i)] + dH(v,i,state)
                new_state[("F",v,i)] = state[("F",v,i)] + dF(v,i,state)
                new_state[("Disch",v,i)] = state[("Disch",v,i)] + dDisch(v,i,state)
            end
         end
         return new_state
    end

    # Set the initial condition.
    state_initial = Dict()
    for v ∈ 0:1
        for i ∈ 1:16
            state_initial[("E",v,i)] = 0
            state_initial[("I",v,i)] = 0
            state_initial[("A",v,i)] = 0
            state_initial[("R",v,i)] = 0
            state_initial[("H",v,i)] = 0
            state_initial[("F",v,i)] = 0
            state_initial[("Disch",v,i)] = 0
        end
     end
    for i ∈ 1:16
        state_initial[("S",0,i)] = N[i]*(1-v[i])
        state_initial[("S",1,i)] = N[i]*v[i]*(1-eI)
        state_initial[("Imm",1,i)] = N[i]*v[i]*eI
    end

    # Now loop over all days.
    states = Dict(1 => state_initial)
    for day = 2:num_days
        states[day] = update(states[day-1])
    end

    # Return the dictionary of states.
    return states

end


#==============================================================================
Try to replicate existing results.
==============================================================================#

# Run simulation with default specs.
states = simulate_SEIR()
# Count infections + recovered after one year in all groups.
total_infections = sum([states[365][("R",v,i)] + states[365][("I",v,i)]
                        for v ∈ 0:1, i ∈ 1:16])
# Count deaths.
total_deaths = sum([states[365][("F",v,i)] for v ∈ 0:1, i ∈ 1:16])

# Compare to
#    Table "Baseline public health measures and full TTIQ"
#       in the more recent paper.
#    Cent VE, 90% vaxed over 12
println("Total infections were $(round(total_infections; digits=1))")
println("They should be 461893.")
println("Total deaths were $(round(total_deaths; digits=1))")
println("They should be 1557.")

# Check that in each state, the number of people across all possibilities add
# to the total pop -- excluding H, D, Disch, which double count people in other
# categories.
function total_pop(state)
    groups = collect(keys(state))
    groups_no_double_count = groups[[g[1] ∉ ["H", "F", "Disch"] for g in groups]]
    states_no_double_count = [state[g] for g in groups_no_double_count]
    return sum(states_no_double_count)
end
@assert all([total_pop(s) for s in values(states)] .≈ N_all)


#==============================================================================
Calculate P[death | vax] for various vax rates.
==============================================================================#

# This function calculates death probability, by vaxed and unvaxed status.
# Note we exclude the three age groups with partial or no vax coverage.
function death_prob(vaxrate12p,   # vaccination rate for those aged 12+
                    ttiq,         # TTIQ system ∈ ["limited", "full"]
                    vax_efficacy; # vaccines effectiveness ∈ ["H", "C", "L"]
                    num_days=1200 # Number of days to run simulation for
                    )
    @assert vax_efficacy ∈ ["H", "C", "L"]
    if vax_efficacy == "H"
        eI = 0.9
        eT = 0.5
        eD = 0.8
    elseif  vax_efficacy == "C"
        eI = 0.7
        eT = 0.5
        eD = 0.8
    else
        eI = 0.5
        eT = 0.4
        eD = 0.8
    end
    @assert ttiq ∈ ["limited", "full"]
    if ttiq == "limited"
        R0 = 6.0*(1-0.17)*(1-0.10)
    else
        R0 = 6.0*(1-0.17)*(1-0.20)
    end
    states = simulate_SEIR(num_days = num_days,
                           v = [[0., 0., vaxrate12p*(3/5)]; vaxrate12p*ones(13)],
                           R0=R0, eI = eI, eT = eT, eD = eD)
    vaxed_pop       = sum([states[1][(c,1,i)]   for c ∈ ["S", "Imm"], i ∈ 4:16])
    unvaxed_pop     = sum([states[1][("S",0,i)] for i ∈ 4:16])
    vaxed_deaths    = sum([states[num_days][("F",1,i)] for i ∈ 4:16])
    unvaxed_deaths  = sum([states[num_days][("F",0,i)] for i ∈ 4:16])
    overall_prob = (vaxed_deaths + unvaxed_deaths)/(vaxed_pop + unvaxed_pop)
    vaxed_prob   = vaxed_deaths/vaxed_pop
    unvaxed_prob = unvaxed_deaths/unvaxed_pop
    return (overall_prob, vaxed_prob, unvaxed_prob)
end

# Calculate probabilities for various specifications, and export.
df = DataFrame(vaxrate12p = [], ttiq = [], vax_efficacy = [], overall_prob = [],
               vaxed_prob = [], unvaxed_prob = [])
for vaxrate12p ∈ 0.8:0.01:1
    for vax_efficacy ∈ ["H", "C", "L"]
        for ttiq ∈ ["limited", "full"]
            (overall_prob, vaxed_prob, unvaxed_prob) = death_prob(vaxrate12p,
                                                        ttiq, vax_efficacy)
            push!(df, (vaxrate12p, ttiq, vax_efficacy, overall_prob, 
                       vaxed_prob, unvaxed_prob))
        end
    end
end
CSV.write("death_probs.csv", df)







