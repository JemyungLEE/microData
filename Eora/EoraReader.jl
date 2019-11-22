# Developed date: 20. Nov. 2019
# Last modified date: 20. Nov. 2019
# Subject: Eora data reader
# Description: Read data from Eora tables
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())

include("FinalDemandReader.jl")
using .FinalDemandReader
fdr = FinalDemandReader

indexFile = Base.source_dir()*"/data/index.xlsx"
finalDemandFile = Base.source_dir()*"/data/2011_eora_y.csv"

println("[Process]")
print(" Index reading: ")
fdr.readIndexData(indexFile)
println("complete")

print(" Final demand reading: ")
fdr.readFinalDemand(finalDemandFile, ["India"])
println("complete")

print(" Test: ")
fdr.test(Base.source_dir()*"/data/test_output.txt")
println("complete")
