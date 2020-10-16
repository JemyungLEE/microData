# Developed date: 11. Jun. 2020
# Last modified date: 15. Oct. 2020
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
categoryFile = filePath * "index/Eurostat_Index_ver1.1.xlsx"
microDataPath = filePath * "microdata/"

readDataFromXLSX = false
readDataFromCSV = true
CurrencyConv = false; erfile = filePath * "index/EUR_USD_ExchangeRates.txt"
PPPConv = false; pppfile = filePath * "index/PPP_ConvertingRates.txt"

codeSubst = false        # recommend 'false' for depth '1st' as there is nothing to substitute
perCap = true

gapMitigation = true    # filling gaps between national account and HBS expenditures
if perCap; eustatsFile = filePath * "index/EU_ConsExp_perCap_COICOP.tsv"
else eustatsFile = filePath * "index/EU_ConsExp_COICOP.tsv"
end

printData = false

year = 2010
catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]
if codeSubst; substTag = "_subst" else substTag = "" end
# microDataPath = [microDataPath*"HU"]
# microDataPath = [microDataPath*"BE", microDataPath*"SE"]

ctgfile = filePath * "extracted/Category_"*depthTag[catDepth]*".csv"
hhsfile = filePath * "extracted/Households.csv"
mmsfile = filePath * "extracted/Members.csv"
expfile = filePath * "extracted/Expenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
sttfile = filePath * "extracted/MicroData_Statistics_"*depthTag[catDepth]*substTag*".csv"
sbstfile = filePath * "extracted/SubstituteCodes_"*depthTag[catDepth]*".csv"

scexpfile = filePath * "extracted/ScaledExpenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
scstatsfile = filePath * "extracted/HBS_COICOP_stats_"*depthTag[catDepth]*substTag*".csv"

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(categoryFile, depth=catDepth, catFile=ctgfile, coicop=gapMitigation)
println("completed")

print(" Micro-data reading: ")
if readDataFromXLSX; print("XLSX")
    mdr.readHouseholdData(year, microDataPath, visible=true, substitute=codeSubst)
    mdr.readMemberData(year, microDataPath, visible=true)
    mdr.buildExpenditureMatrix(year, expfile, substitute=codeSubst)
    mdr.makeStatistics(year, sttfile, substitute=codeSubst)
elseif readDataFromCSV; print("CSV")
    mdr.readPrintedHouseholdData(hhsfile)
    mdr.readPrintedMemberData(mmsfile)
    # mdr.readSubstCodesCSV(sbstfile)
    mdr.readPrintedExpenditureData(expfile, substitute=codeSubst, buildHhsExp=true)
    # mdr.makeStatistics(year, sttfile, substitute=codeSubst)
end
println(" completed")

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

if gapMitigation; print(" HBS-COICOP gap mitigating: ")
    mdr.mitigateExpGap(year, eustatsFile, scexpfile, scstatsfile; percap=perCap, subst=codeSubst, checking=true)
    println("completed")
end

if printData; print(" Extracted data printing:")
    mdr.printCategory(year, ctgfile, substitute=codeSubst)
    # mdr.printHouseholdData(year, hhsfile)
    # mdr.printMemberData(year, mmsfile)
    println("completed")
end

println("[done]")
