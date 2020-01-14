# Developed date: 27. Dec. 2019
# Last modified date: 14. Jan. 2019
# Subject: Emission mapping
# Description: Mapping emission through households emissions data
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("EmissionCategorizer.jl")

using .EmissionCategorizer

ec = EmissionCategorizer

println("[Process]")

nation = "IND"
year = 2011
emissionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission.txt"
householdFile = Base.source_dir() * "/data/extracted/Households.txt"
sectorFile = "../Eora/data/Eora_HS_match.xlsx"

print(" Data reading: ")
print("household")
ec.readHouseholdData(year, householdFile)
print(", category")
ec.readCategoryData(nation, sectorFile)
print(", emission")
ec.readEmission(year, emissionFile)
println(" ... complete")

categorizedFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_cat_popHhW.txt"
weightMode = 3
print(" Categorizing: ")
ec.categorizeEmission(year, weightMode)
ec.printCategorizedEmission(year, categorizedFile)
println("complete")

#=
print("Emission files comparison: ")
emissionFile = []
push!(emissionFile, path * "2011_IND_hhs_emission.txt")
push!(emissionFile, path * "2011_IND_HH2_emission.txt")
push!(emissionFile, path * "2011_IND_hh3(true)_emission.txt")
push!(emissionFile, path * "2011_IND_hh3(false)_emission.txt")
ec.compareTables(2011, emissionFile)
=#
