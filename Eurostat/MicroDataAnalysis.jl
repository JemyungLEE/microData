# Developed date: 11. Jun. 2020
# Last modified date: 13. Jun. 2020
# Subject: EU Household Budget Survey (HBS) microdata analysis
# Description: proceed data analysis process for EU HBS microdata
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
include("MicroDataReader.jl")
using .MicroDataReader
mdr = MicroDataReader
filePath = Base.source_dir() * "/data/"
categoryFile = filePath * "index/Eurostat_Index_ver0.5.xlsx"
microDataPath = filePath * "microdata/"

year = 2010

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(categoryFile, depth=4)
println("completed")

print(" Micro-data reading:")
mdr.readHouseholdData(year, microDataPath, visible=true)
mdr.readMemberData(year, microDataPath, visible=true)
mdr.buildExpenditureMatrix(year, filePath * "extracted/Expenditure_matrix.csv")
mdr.makeStatistics(year, filePath * "extracted/MicroData_Summary.txt")
println("completed")

print(" Extracted data printing:")
mdr.printHouseholdData(year, filePath * "extracted/Households.csv")
mdr.printMemberData(year, filePath * "extracted/Members.csv")
println("completed")

println("[done]")
