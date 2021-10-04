#==============================================================================
This script replicates the TPM SEIR vaccine model.

Tested on Julia 1.6.3.
Dependencies:
  CSV v0.9.5
  DataFrames v1.2.2
  XLSX v0.7.8

==============================================================================#

include("seir.jl")

pVaxEligs = 0.7:0.05:0.95
ttiqs = ["limited", "full"]
vax_efficacies = ["H", "C", "L"]
eligibilities = ["12+", "5+"]

function age_group_eligibility(eligible)
    @assert eligible ∈ ["5+", "12+"]
    if eligible == "5+"
        return [[0.]; ones(15)]
    elseif eligible == "12+"
        return [[0., 0., 3/5]; ones(13)]
    end
end

# lay out in same format as in TPM report
df = DataFrame(
    infections_H=Float64[], infections_C=Float64[], infections_L=Float64[],
    hospital_H=Float64[], hospital_C=Float64[], hospital_L=Float64[],
    deaths_H=Float64[], deaths_C=Float64[], deaths_L=Float64[],
    peak_beds_H=Float64[], peak_beds_C=Float64[], peak_beds_L=Float64[],
    )

row = 1

for ttiq in ttiqs, eligible in eligibilities, pVaxElig in pVaxEligs

    print("running $ttiq, $eligible, $pVaxElig...\r")
    push!(df, zeros(ncol(df)))

    for vax_efficacy in vax_efficacies

        # The efficacy of the vaccine against infection, transmission, severe disease
        eI, eT, eD = vax_efficacy_params(vax_efficacy)

        # The basic reproduction number.
        R0 = R0_under_ttiq(ttiq)

        # The vaccination rate for each age group.
        v = pVaxElig * age_group_eligibility(eligible)

        # Run simulation with default specs.
        states = simulate_SEIR(eI=eI, eT=eT, eD=eD, R0=R0, v=v)

        # Count infections + recovered after one year in all groups.
        total_infections = sum([states[365][("R",v,i)] + states[365][("I",v,i)]
                                for v ∈ 0:1, i ∈ 1:16])
        # Count hospitalisations.
        total_hospitalisations = sum([states[365][("H",v,i)] + states[365][("Disch",v,i)]
                                     for v ∈ 0:1, i ∈ 1:16])
        # Count deaths.
        total_deaths = sum([states[365][("F",v,i)] for v ∈ 0:1, i ∈ 1:16])
        # Find peak beds occupied.
        peak_beds_occupied = maximum([sum([states[day][("H",v,i)] for v ∈ 0:1, i ∈ 1:16]) for day ∈ 1:365])

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

        df[row, "infections_$vax_efficacy"] = total_infections
        df[row, "hospital_$vax_efficacy"] = total_hospitalisations
        df[row, "deaths_$vax_efficacy"] = total_deaths
        df[row, "peak_beds_$vax_efficacy"] = peak_beds_occupied

    end

    global row += 1
end

CSV.write("replication.csv", df)
