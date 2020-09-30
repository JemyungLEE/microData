# Developed date: 20. Aug. 2020
# Last modified date: 30. Sep. 2020
# Subject: EU Household Budget Survey (HBS) microdata integrity check
# Description: Verify the integrity og HBS microdata
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
include("MicroDataReader.jl")
using .MicroDataReader
mdr = MicroDataReader

filePath = Base.source_dir() * "/data/"
categoryFile = filePath * "index/Eurostat_Index_ver0.9.xlsx"

codeSubst = true
fixedComp = true

year = 2010
depthTag = ["1st", "2nd", "3rd", "4th"]
if codeSubst; substTag = "_subst" else substTag = "" end
if fixedComp; fixedTag = "_fixed" else fixedTag = "" end

startDepth = 1
endDepth = 4

expfile = [filePath * "extracted/Expenditure_matrix_"*depthTag[i]*".csv" for i = startDepth:endDepth]
sbstfile = [filePath * "extracted/SubstituteCodes_"*depthTag[i]*".csv" for i = startDepth:endDepth]

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(categoryFile, depth=endDepth)
println("completed")

print(" Depth integrity check: ")
integrityFile = [filePath * "check/Integiry_"*depthTag[i]*substTag*fixedTag*".csv" for i = startDepth:endDepth-1]
mdr.checkDepthIntegrity(year, expfile, integrityFile, startDepth=startDepth, subst=codeSubst, fixed=fixedComp)
println("completed")

println("[done]")
