# Developed date: 11. Jun. 2020
# Last modified date: 18. Feb. 2021
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
categoryFile = filePath * "index/Eurostat_Index_ver4.0.xlsx"
microDataPath = filePath * "microdata/"

readDataFromXLSX = false
readDataFromCSV = true
CurrencyConv = false; erfile = filePath * "index/EUR_USD_ExchangeRates.txt"
PPPConv = false; pppfile = filePath * "index/PPP_ConvertingRates.txt"

codeSubst = true        # recommend 'false' for depth '1st' as there is nothing to substitute
perCap = true

gapMitigation = true    # filling gaps between national account and HBS expenditures
if perCap; eustatsFile = filePath * "index/EU_ConsExp_perCap_COICOP.tsv"
else eustatsFile = filePath * "index/EU_ConsExp_COICOP.tsv"
end

cpiScaling = true; cpi_std_year = 2010
cpi_file = filePath * "index/EU_hicp.tsv"

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

scexpfile = filePath * "extracted/Scaled_Expenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
scstatsfile = filePath * "extracted/HBS_COICOP_stats_"*depthTag[catDepth]*substTag*".csv"

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(year, categoryFile, depth=catDepth, catFile=ctgfile, coicop=gapMitigation)
println("completed")

print(" Micro-data reading: ")
if readDataFromXLSX; print("XLSX")
    mdr.readHouseholdData(year, microDataPath, visible=true, substitute=codeSubst)
    mdr.readMemberData(year, microDataPath, visible=true)
    mdr.buildExpenditureMatrix(year, expfile, substitute=codeSubst)
elseif readDataFromCSV; print("CSV")
    mdr.readPrintedHouseholdData(hhsfile)
    mdr.readPrintedMemberData(mmsfile)
    if codeSubst; mdr.readSubstCodesCSV(sbstfile) end
    mdr.readPrintedExpenditureData(expfile, substitute=codeSubst, buildHhsExp=true)
end
mdr.makeStatistics(year, sttfile, substitute=codeSubst)
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
    mdr.mitigateExpGap(year, eustatsFile, scexpfile, scstatsfile, percap=perCap, subst=codeSubst, cdrepl=true, alter=true)
    println("completed")
end

if cpiScaling; print(" CPI scaling: ")
    mdr.readCPIs([2010, 2015], cpi_file, idx_sep = ',', freq="A", unit="INX_A_AVG", topLev = "EU")
    mdr.scalingByCPI(year, cpi_std_year, codeDepth=0, topLev = "EU", subst = codeSubst)
    println("completed")
end

if printData; print(" Extracted data printing:")
    mdr.printCategory(year, ctgfile, substitute=codeSubst)
    mdr.printHouseholdData(year, hhsfile)
    mdr.printMemberData(year, mmsfile)
    println("completed")
end

println("[done]")
