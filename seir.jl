#==============================================================================
Contains parameters and functions for replicating the TPM SEIR vaccine model.

Tested on Julia 1.6.3.
Dependencies:
  CSV v0.9.5
  DataFrames v1.2.2
  XLSX v0.7.8

==============================================================================#
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
C_filepath = joinpath(@__DIR__, "MUestimates_all_locations_2.xlsx")
C_raw = XLSX.readxlsx(C_filepath)["New Zealand"][:]
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


"""Returns the vaccine efficacy parameters under a given level, `"H"`, `"C"` or
`"L"` (high, central or low, as defined in the TPM report)."""
function vax_efficacy_params(level)
    @assert level ∈ ["H", "C", "L"]
    if level == "H"
        eI = 0.9
        eT = 0.5
        eD = 0.8
    elseif level == "C"
        eI = 0.7
        eT = 0.5
        eD = 0.8
    else
        eI = 0.5
        eT = 0.4
        eD = 0.8
    end
    return (eI, eT, eD)
end


"""Returns the Reff value under a given TTIQ assumption, `"limited"` or `"full"`
as defined in the TPM report."""
function R0_under_ttiq(ttiq)
    @assert ttiq ∈ ["none", "limited", "full"]
    if ttiq == "limited"
        R0 = 6.0*(1-0.17)*(1-0.10)
    elseif ttiq == "full"
        R0 = 6.0*(1-0.17)*(1-0.20)
    else
        R0 = 6.0
    end
    return R0
end


"""Returns a vector containing the proportions of each age group that would be
vaccinated, given an age eligibility settings, either `"5+"` or `"12+"`."""
function age_group_eligibility(eligible)
    @assert eligible ∈ ["5+", "12+"]
    if eligible == "5+"
        return [[0.]; ones(15)]
    elseif eligible == "12+"
        return [[0., 0., 3/5]; ones(13)]
    end
end



#==============================================================================
Define a function to simulate an SEIR model.

Note we store the state in a dictionary with keys (θ,v,i), where
   θ ∈ ["S", "E", "I", "A", "R", "Imm", "H", "F", "Disch"] indexes compartments;
   v ∈ [0, 1] indexes vax status; and
   i indexes ages.

This function returns a dict of states, indexed by day.
==============================================================================#


"""Runs a bisection search to find the value of `x` that would make `f(x)`
within `tol` of `target`, using `lower` and `upper` as initial search
boundaries. Shows a warning if the result is either `lower` or `upper`.
(`f` is assumed to be increasing; behaviour is undefined if it's not.)
"""
function bisection_search(f, target, lower, upper, tol=1e-7)
    if target > f(upper)
        @warn "bisection_search: target > f(upper), returning upper"
        return upper
    elseif target < f(lower)
        @warn "bisection_search: target < f(lower), returning lower"
        return lower
    end

    x = (lower + upper) / 2
    while abs(f(x) - target) > tol
        value = f(x)
        if target < value
            upper = x
        else
            lower = x
        end
        x = (lower + upper) / 2
    end

    return x
end


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
    # We set U to match the basic reproduction number (within a tolerance)
    R0_given_U(U) = maximum(eigvals(NGM(U)))
    U = bisection_search(R0_given_U, R0, 0.01, 1)

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