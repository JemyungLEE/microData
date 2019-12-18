# Developed date: 3. Dec. 2019
# Last modified date: 9. Dec. 2019
# Subject: Eora data analyzer
# Description: Proceed data anlysis through Eora MRIO tables
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())

include("EmissionCalculator.jl")
using .EmissionCalculator
ec = EmissionCalculator

year = 2011
path = "/Users/leejimac/github/microData/Eora/data/"
indexFile = path * "Eora_HS_match.xlsx"
tTableFile = path * string(year) * "_eora_t.csv"
vTableFile = path * string(year) *  "_eora_v.csv"
yTableFile = path * string(year) *  "_eora_y.csv"
qTableFile = path * string(year) *  "_eora_q.csv"

println("[Process]")
print(" Index reading: ")
ec.readIndexXlsx(indexFile)
println("complete")

print(" MRIO table reading: ")
ec.readIOTables(year, tTableFile, vTableFile, yTableFile, qTableFile)
println("complete")

print(" Emission calculating: ")
ec.rearrangeIndex()
#ec.getMatchingIndex()
ec.rearrangeTables(year)
ec.calculateEmission(year)
ec.extractHouseholdEmission(year, true)
println("complete")

print(" Emission results printing: ")
nation = "IND"
ec.printEmissions(year, path * string(year) * "_emission.txt")
ec.getNationEmission(year, nation, true)
#ec.getEmissionDataset(year, nation)
ec.printEmissions(year, path * string(year) * "_"*nation*"_emission.txt")
println("complete")

#=
print(" MRIO table printing: ")
ec.printTables(year, path * "Test_print.txt")
println("complete")
=#
