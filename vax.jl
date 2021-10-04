#==============================================================================
Calculate P[death | vax] for various vax rates.
==============================================================================#

include("seir.jl")

# This function calculates death probability, by vaxed and unvaxed status.
# Note we exclude the three age groups with partial or no vax coverage.
function death_prob(vaxrate12p,   # vaccination rate for those aged 12+
                    ttiq,         # TTIQ system ∈ ["limited", "full"]
                    vax_efficacy; # vaccines effectiveness ∈ ["H", "C", "L"]
                    num_days=365  # Number of days to run simulation for
                    )
    eI, eT, eD = vax_efficacy_params(vax_efficacy)
    R0 = R0_under_ttiq(ttiq)

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







