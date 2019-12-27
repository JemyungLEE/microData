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

print("Emission files comparison: ")
path = Base.source_dir()*"/data/emission/"
emissionFile = []
push!(emissionFile, path * "2011_IND_hhs_emission.txt")
push!(emissionFile, path * "2011_IND_HH2_emission.txt")
push!(emissionFile, path * "2011_IND_hh3(true)_emission.txt")
push!(emissionFile, path * "2011_IND_hh3(false)_emission.txt")
ec.compareTables(2011, emissionFile)
println("complete")
