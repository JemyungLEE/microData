
# Developed date: 5. Aug. 2019
# Last modified date: 5. Aug. 2019
# Subject: commodity classification analyzer, KR
# Description: activate a commodity classification analysis process
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
include("CommClassReader.jl")
using .CommClassReader
ccr = CommClassReader

inputFile = "Comm_Class_KR.txt"
inputFile = Base.source_dir()*"/"*inputFile

codeList = ccr.readCommCode(inputFile)

for c in codeList
    println(c)
end
