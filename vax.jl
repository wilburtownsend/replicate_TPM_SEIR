#==============================================================================
This script replicates the TPM SEIR vaccine model.

Sources:
    https://cpb-ap-se2.wpmucdn.com/blogs.auckland.ac.nz/dist/d/75/files/2017/01/a-covid-19-vaccination-model-for-aotearoa.pdf
    https://cpb-ap-se2.wpmucdn.com/blogs.auckland.ac.nz/dist/d/75/files/2017/01/a-covid-19-vaccination-model-for-aotearoa-new-zealand-supplementary.pdf
    https://cpb-ap-se2.wpmucdn.com/blogs.auckland.ac.nz/dist/d/75/files/2017/01/modelling-to-support-a-future-covid-19-strategy-for-aotearoa-new-zealand.pdf

==============================================================================#
cd("/Users/wtownsend/Documents/GitHub/replicate_TPM_SEIR/")
using LinearAlgebra, XLSX


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
# Contact matrix -- from Prem et al.
C_raw = XLSX.readxlsx("MUestimates_all_locations_2.xlsx")["New Zealand"][:]
# this makes it balanced.
C = [(C_raw[i,j] + (N[j]/N[i])*C_raw[j,i])/2 for i ∈ 1:16, j ∈ 1:16]
# Latent period
tE = 2.55
# Infectious period
tI = 5
# The relative infectiousness of subclinical individuals
τ = 0.5
# Proportion of infections causing clinical disease by age group
pclin = [0.544, 0.555, 0.577, 0.5985, 0.6195, 0.6395, 0.6585,
         0.6770, 0.6950, 0.7118, 0.7273, 0.7418, 0.7553, 0.768,
         0.78, 0.8008]
# Relative susceptibility by age group (compared to 60-64 year-olds)
u = [0.462, 0.457, 0.445, 0.558, 0.795, 0.934, 0.974, 0.977,
     0.942, 0.931, 0.942, 0.965, 1.00, 0.977, 0.896, 0.856]
# Imported cases per age group per day.
m = age_dist


#==============================================================================
Define a function to simulate an SEIR model.

Note we store the state in a dictionary with keys (θ,v,i), where 
   θ ∈ ["S", "E", "I", "A", "R", "Imm"] indexes compartments;
   v ∈ [0, 1] indexes vax status; and
   i indexes ages.

This function returns a dict of states, indexed by day.
==============================================================================#

function simulate_SEIR(;
    num_days = 730, # Number of days to simulate.
    eT = 0.5, # The efficacy of the vaccine against transmission.
    eI = 0.7, # The efficacy of the vaccine against infection.
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
Try different simulations.
==============================================================================#

states = simulate_SEIR()
# Count infections + recovered after one year in all groups.
sum([states[365][("R",v,i)] + states[365][("I",v,i)]
    for v ∈ 0:1, i ∈ 1:16])
# Compare to
#    Table "Baseline public health measures and full TTIQ"
#       in the more recent paper.
#    Infections, Cent VE, 90% vaxed over 12
#    Should be 461893

