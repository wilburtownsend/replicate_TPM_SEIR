#==============================================================================
Calculate P[death | vax] for various vax rates.
==============================================================================#

include("seir.jl")

"""Calculates the probabilities of being in a given `compartment` after
`num_days`, by vaxed and unvaxed status. Note we exclude the three ages groups
with partial or no vax coverage."""
function compartment_prob(vaxrate12p,   # vaccination rate for those aged 12+
                          ttiq,         # TTIQ system ∈ ["limited", "full"]
                          vax_efficacy, # vaccines effectiveness ∈ ["H", "C", "L"]
                          compartments, # compartment, e.g. "F" for fatalities
                          num_days=365  # Number of days to run simulation for
                          )
    eI, eT, eD = vax_efficacy_params(vax_efficacy)
    R0 = R0_under_ttiq(ttiq)
    v = vaxrate12p * age_group_eligibility("12+")

    states = simulate_SEIR(num_days = num_days, v = v,
                           R0 = R0, eI = eI, eT = eT, eD = eD)
    vaxed_pop        = sum([states[1][(c,1,i)]   for c ∈ ["S", "Imm"], i ∈ 4:16])
    unvaxed_pop      = sum([states[1][("S",0,i)] for i ∈ 4:16])
    vaxed_affected   = sum([states[num_days][(c,1,i)] for c ∈ compartments, i ∈ 4:16])
    unvaxed_affected = sum([states[num_days][(c,0,i)] for c ∈ compartments, i ∈ 4:16])
    overall_prob = (vaxed_affected + unvaxed_affected)/(vaxed_pop + unvaxed_pop)
    vaxed_prob   = vaxed_affected/vaxed_pop
    unvaxed_prob = unvaxed_affected/unvaxed_pop
    return (overall_prob, vaxed_prob, unvaxed_prob)
end

compartments = ["F"] # deaths
# compartments = ["H", "F", "Disch"] # hospitalisations

# Calculate probabilities for various specifications, and export.
df = DataFrame(vaxrate12p = [], ttiq = [], vax_efficacy = [], overall_prob = [],
               vaxed_prob = [], unvaxed_prob = [])
for vaxrate12p ∈ 0.8:0.01:1, vax_efficacy ∈ ["H", "C", "L"], ttiq ∈ ["limited", "full"]
    print("vaxrate12p $vaxrate12p, vax_efficacy $vax_efficacy, ttiq $ttiq...\r")
    (overall_prob, vaxed_prob, unvaxed_prob) = compartment_prob(vaxrate12p, ttiq,
                                                                vax_efficacy, compartments)
    push!(df, (vaxrate12p, ttiq, vax_efficacy, overall_prob,
               vaxed_prob, unvaxed_prob))
end

CSV.write("probs_$(join(compartments)).csv", df)
