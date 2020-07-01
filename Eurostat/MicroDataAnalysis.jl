# Developed date: 11. Jun. 2020
# Last modified date: 01. Jul. 2020
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
categoryFile = filePath * "index/Eurostat_Index_ver0.6.xlsx"
microDataPath = filePath * "microdata/"

readDataFromXLSX = true
readDataFromCSV = false
CurrencyConv = true; erfile = filePath * "index/EUR_USD_ExchangeRates.txt"
PPPConv = true; pppfile = filePath * "index/PPP_ConvertingRates.txt"

printData = true

year = 2010
catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(categoryFile, depth=catDepth)
println("completed")

hhsfile = filePath * "extracted/Households.csv"
mmsfile = filePath * "extracted/Members.csv"
expfile = filePath * "extracted/Expenditure_matrix_"*depthTag[catDepth]*".csv"
sttfile = filePath * "extracted/MicroData_Statistics_"*depthTag[catDepth]*".csv"

if readDataFromXLSX
    print(" Micro-data reading: XLSX")
    mdr.readHouseholdData(year, microDataPath, visible=true)
    mdr.readMemberData(year, microDataPath, visible=true)
    mdr.buildExpenditureMatrix(year, expfile)
    mdr.makeStatistics(year, sttfile)
    println(" completed")
elseif readDataFromCSV
    print(" Micro-data reading: CSV")
    mdr.readPrintedHouseholdData(hhsfile)
    mdr.readPrintedMemberData(mmsfile)
    mdr.readPrintedExpenditureData(expfile, buildTable=true)
    mdr.makeStatistics(year, sttfile)
    println(" completed")
end

if CurrencyConv; print(" Currency exchanging: "); mdr.exchangeExpCurrency(erfile); println("complete") end
if PPPConv; print(" PPP converting: ");  mdr.convertToPPP(pppfile); println("complete") end
if CurrencyConv || PPPConv; mdr.makeStatistics(year, replace(sttfile,".csv"=>"_USD.csv")) end

if printData
    print(" Extracted data printing:")
    mdr.printHouseholdData(year, replace(hhsfile, ".csv"=>"_test.csv"))
    mdr.printMemberData(year, replace(mmsfile, ".csv"=>"_test.csv"))
    println("completed")
end

println("[done]")
