# Developed date: 11. Jun. 2020
# Last modified date: 18. Jun. 2020
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

readDataFromXLSX = false
readDataFromCSV = true
CurrencyConv = false; erfile = filePath * "index/EUR_USD_ExchangeRates.txt"
PPPConv = false; pppfile = filePath * "index/PPP_ConvertingRates.txt"

printData = false

year = 2010

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(categoryFile, depth=4)
println("completed")

hhsfile = filePath * "extracted/Households.csv"
mmsfile = filePath * "extracted/Members.csv"
expfile = filePath * "extracted/Expenditure_matrix.csv"
sttfile = filePath * "extracted/MicroData_Statistics.txt"

if readDataFromXLSX
    print(" Micro-data reading:")
    mdr.readHouseholdData(year, microDataPath, visible=true)
    mdr.readMemberData(year, microDataPath, visible=true)
    mdr.buildExpenditureMatrix(year, expfile)
    mdr.makeStatistics(year, sttfile)
    println("completed")
end

if readDataFromCSV
    mdr.readPrintedHouseholdData(hhsfile)
    mdr.readPrintedMemberData(mmsfile)
    mdr.readPrintedExpenditureData(expfile, buildTable=true)
    mdr.makeStatistics(year, sttfile)
end

if CurrencyConv; mdr.exchangeExpCurrency(erfile) end
if PPPConv; mdr.convertToPPP(pppfile) end

if printData
    print(" Extracted data printing:")
    mdr.printHouseholdData(year, hhsfile)
    mdr.printMemberData(year, mmsfile)
    println("completed")
end

println("[done]")
