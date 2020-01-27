# Developed date: 27. Dec. 2019
# Last modified date: 20. Jan. 2019
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

categorizedFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_cat_perCap.txt"
weightMode = 6  # [0]non-weight, [1]population weighted, [2]household weighted, [3]both population and household weighted
                # ([4],[5]: normalization) [4]per capita, [5]per household
                # (basic information) [6]population and households, ec[1,:]population, ex[2,:]households
eqvalMode = false # [true]apply square root oh household size for equivalance scale
print(" Categorizing: ")
ec.categorizeEmission(year, weightMode, eqvalMode)
ec.printCategorizedEmission(year, categorizedFile)
println("complete")

#exportFile = Base.source_dir() * "/data/emission/2011_IND_hhs_GIS_emission_cat_perCap_.csv"
exportFile = Base.source_dir() * "/data/emission/2011_IND_hhs_GIS_demography.csv"
print(" Exporting: ")
ec.exportEmissionTable(year, "GID_2", exportFile)
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
