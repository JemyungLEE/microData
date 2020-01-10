# Developed date: 27. Dec. 2019
# Last modified date: 27. Dec. 2019
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

print(" Data reading: emission")
ec.readEmission(year, emissionFile)
print(", household")
ec.readHousehold(year, householdFile)
print(", sector")
ec.readSectors(nation, sectorFile)
println(" ... complete")

print(" Categorizing: ")
ec.categorizeEmission
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
