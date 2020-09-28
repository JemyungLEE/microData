# Developed date: 11. Jun. 2020
# Last modified date: 28. Sep. 2020
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
categoryFile = filePath * "index/Eurostat_Index_ver0.9.xlsx"
microDataPath = filePath * "microdata/"

readDataFromXLSX = true
readDataFromCSV = false
CurrencyConv = false; erfile = filePath * "index/EUR_USD_ExchangeRates.txt"
PPPConv = false; pppfile = filePath * "index/PPP_ConvertingRates.txt"

codeSubst = true

printData = true

year = 2010
catDepth = 1
depthTag = ["1st", "2nd", "3rd", "4th"]
# microDataPath = [microDataPath*"BE", microDataPath*"SE"]

ctgfile = filePath * "extracted/Category_"*depthTag[catDepth]*".csv"
hhsfile = filePath * "extracted/Households.csv"
mmsfile = filePath * "extracted/Members.csv"
expfile = filePath * "extracted/Expenditure_matrix_"*depthTag[catDepth]*".csv"
sttfile = filePath * "extracted/MicroData_Statistics_"*depthTag[catDepth]*".csv"
sbstfile = filePath * "extracted/SubstituteCodes_"*depthTag[catDepth]*".csv"

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(categoryFile, depth=catDepth, catFile= ctgfile)
println("completed")

if readDataFromXLSX
    print(" Micro-data reading: XLSX")
    mdr.readHouseholdData(year, microDataPath, visible=true, substitute=codeSubst)
    mdr.readMemberData(year, microDataPath, visible=true)
    mdr.buildExpenditureMatrix(year, expfile, substitute=codeSubst)
    mdr.makeStatistics(year, sttfile, substitute=codeSubst)
    println(" completed")
elseif readDataFromCSV
    print(" Micro-data reading: CSV")
    mdr.readPrintedHouseholdData(hhsfile)
    mdr.readPrintedMemberData(mmsfile)
    mdr.readPrintedExpenditureData(expfile, buildTable=true)
    mdr.makeStatistics(year, sttfile)
    println(" completed")
end

if CurrencyConv; print(" Currency exchanging: ")
    mdr.exchangeExpCurrency(erfile)
    mdr.buildExpenditureMatrix(year, replace(expfile, ".csv"=>"_USD.csv"))
    println("complete")
end
if PPPConv; print(" PPP converting: ")
    mdr.convertToPPP(pppfile)
    hhsfile = replace(hhsfile,".csv"=>"_PPP.csv")
    mmsfile = replace(mmsfile,".csv"=>"_PPP.csv")
    println("complete")
end
if CurrencyConv || PPPConv; mdr.makeStatistics(year, replace(sttfile,".csv"=>"_USD.csv")) end

if printData
    print(" Extracted data printing:")
    mdr.printCategory(year, ctgfile, substitute=codeSubst)
    mdr.printHouseholdData(year, hhsfile)
    mdr.printMemberData(year, mmsfile)
    println("completed")
end

println("[done]")
