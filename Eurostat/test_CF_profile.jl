clearconsole()
cd(Base.source_dir())

include("EmissionCategorizer.jl")
using .EmissionCategorizer
ec = EmissionCategorizer
include("../GIS/QgisStyleExporter.jl")
using .QgisStyleExporter
qse = QgisStyleExporter

println("[Process]")

nation = "Eurostat"
year = 2015
years = [2015]
nutsLv = 1
onlyNutsInHbs = true
Qtable = "I_CHG_CO2"

# Qtable = "_PRIMAP"

ceIntegrateMode = "cf"      # "ie" (only indirect CE), "de" (only direct CE), or "cf" (integrage direct and indirect CEs)

catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]

cpi_scaling = true; base_year = 2010

substMode = true; if substMode; substTag = "_subst" else substTag = "" end
scaleMode = true; if scaleMode; scaleTag = "Scaled_" else scaleTag = "" end

eqvalMode = false   # [true]: apply square root of household size for equivalance scale
ntWeighMode = true  # [true]: apply NUTS population based weight, [false]:apply HBS weight


popweight = true
grid_pop = true
expenditureMode = false

filePath = Base.source_dir() * "/data/"
indexPath = filePath * "index/"
extrPath = filePath * "extracted/"
emissPath = filePath * "emission/" * string(year) * "/"
indexFile = indexPath * "Eurostat_Index_ver5.0.xlsx"
hhsfile = extrPath * string(year) * "_Households.csv"

ExpenditureFile = extrPath * scaleTag * "Expenditure_matrix_4th" * substTag * ".csv"

incomePeriod = "annual"  # Period: "annual"(default), "monthly", or "daily"

normTag = ["perCapNorm", "perHhNorm"]
categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
subcat=""

print(" Data reading: ")
print("category"); ec.readCategoryData(indexFile, year, nutsLv, except=["None"], subCategory=subcat)
if length(subcat)==0; ec.setCategory(categories); else ec.setCategory(foodCategories); subcat *= "_" end
print(", household"); ec.readHouseholdData(hhsfile, period = incomePeriod, remove = onlyNutsInHbs, alter=true)
print(", population"); ec.readPopulation(year, indexFile, nuts_lv = nutsLv)
print(", gridded population"); ec.readPopGridded(year, indexFile, nuts_lv = [nutsLv], adjust = true)
print(", emission")
if !expenditureMode
    IE_files = []; DE_files = []; ie_nations = []; de_nations = []
    ie_file_tag = "_hhs_"*scaleTag*"IE"*Qtable*".txt"
    if cpi_scaling && base_year != year; ie_file_tag = replace(ie_file_tag, ".txt" => "_converted_" * string(base_year) * ".txt") end
    de_file_tag = "_hhs_"*scaleTag*"DE.txt"
    for f in readdir(emissPath)
        if startswith(f, string(year)) && endswith(f, ie_file_tag); push!(IE_files, emissPath*f); push!(ie_nations, f[6:7])
        elseif startswith(f, string(year)) && endswith(f, de_file_tag); push!(DE_files, emissPath*f); push!(de_nations, f[6:7])
        end
    end
    print("_IE"); ec.readEmissionData(year, ie_nations, IE_files, mode = "ie")
    print("_DE"); ec.readEmissionData(year, de_nations, DE_files, mode = "de")
    print("_CF"); ec.integrateCarbonFootprint(year, mode=ceIntegrateMode)
end
println(" ... complete")

print(" Weights calculating: ")
ec.calculateNutsPopulationWeight(year = year, pop_dens = grid_pop, adjust = true)
println(" ... complete")

print(" Profiling:")
print(" indirect emission")
profileFile = emissPath * string(year) * "_EU_hhs_"*scaleTag*subcat*"_CF_profile.txt"
ec.profilingEmissionByIncome(year, [], profileFile, n_top = 10, mode = "cf", period="annual", sort_mode="income_pc", adjust=true, popWgh=popweight, ntweigh=ntWeighMode)

print(", expenditure")
profileFile = extrPath * string(year) * "_EU_hhs_"*scaleTag*subcat*"_Exp_profile.txt"
expFile = extrPath * string(year) * "_" * scaleTag*"Expenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
ec.profilingExpenditureByIncome(year, [], expFile, profileFile, n_top = 10, period="annual", sort_mode="income_pc", adjust=true, popWgh=popweight, ntweigh=ntWeighMode)


println("\n[Done]")
