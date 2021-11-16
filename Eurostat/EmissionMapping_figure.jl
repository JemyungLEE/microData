# Developed date: 11. Nov. 2021
# Last modified date: 11. Nov. 2021
# Subject: Categorized emission mapping
# Description: Mapping emission through households emissions data, categorizing by district, income-level, and etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())

include("EmissionCategorizer.jl")
using .EmissionCategorizer
ec = EmissionCategorizer

println("[Process]")

nation = "Eurostat"
year = 2015
nutsLv = 1
onlyNutsInHbs = true
Qtable = "_PRIMAP"
ceIntegrateMode = "cf"      # "ie" (only indirect CE), "de" (only direct CE), or "cf" (integrage direct and indirect CEs)
ceProcessMode = ["ie", "de", "cf"]

cpi_scaling = true; base_year = 2010

substMode = true; if substMode; substTag = "_subst" else substTag = "" end
scaleMode = true; if scaleMode; scaleTag = "Scaled_" else scaleTag = "" end

ntWeighMode = true  # [true]: apply NUTS population based weight, [false]:apply HBS weight

popweight = true
grid_pop = true

filePath = Base.source_dir() * "/data/"
indexPath = filePath * "index/"
extrPath = filePath * "extracted/"
emissPath = filePath * "emission/" * string(year) * "/"
indexFile = indexPath * "Eurostat_Index_ver4.6.xlsx"
hhsfile = extrPath * string(year) * "_Households.csv"

ExpenditureFile = extrPath * scaleTag * "Expenditure_matrix_4th" * substTag * ".csv"

incomePeriod = "daily"  # Period: "annual", "monthly"(default), or "daily"

normTag = ["perCapNorm", "perHhNorm"]
categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
subcat=""
# subcat="Food"
foodCategories=["Grain","Vegetable","Fruit","Dairy","Beef","Pork","Poultry","Other meat","Fish",
                "Alcohol","Other beverage","Confectionery","Restaurant","Other food","Food"]

print(" Data reading: ")
print("category"); ec.readCategoryData(indexFile, year, nutsLv, except=["None"], subCategory=subcat)
if length(subcat)==0; ec.setCategory(categories); else ec.setCategory(foodCategories); subcat *= "_" end
print(", household"); ec.readHouseholdData(hhsfile, period = incomePeriod, remove = onlyNutsInHbs, alter=true)
print(", population"); ec.readPopulation(year, indexFile, nuts_lv = nutsLv)
print(", gridded population"); ec.readPopGridded(year, indexFile, nuts_lv = [nutsLv], adjust = true)
print(", emission")
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
println(" ... complete")

print(" Weights calculating: ")
ec.calculateNutsPopulationWeight(year = year, pop_dens = grid_pop, adjust = true)
println(" ... complete")

print(" Categorizing:")
print(" category")
for m in ceProcessMode
    print("_",uppercase(m))
    ec.categorizeHouseholdEmission(year, mode=m, output="", hhsinfo=false, nutsLv=1)
end
println(" ... complete")

print(" Sorting: ")
    fileTag = extrPath * "/sorted/" * "YEAR_NATION_sortedHHsByCF.txt"
    ec.sortHHsByCF(year, [], fileTag, mode = "cf")
println(" ... complete")

println("\n[Done]")