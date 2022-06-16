# Developed date: 11. Jun. 2020
# Last modified date: 15. Jun. 2022
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
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"

categoryFile = indexFilePath * "Eurostat_Index_ver5.0.xlsx"
eustatsFile = indexFilePath * "EU_exp_COICOP.tsv"
cpi_file = indexFilePath * "EU_hicp.tsv"

readDataFromXLSX = true; readDataFromCSV = !readDataFromXLSX

CurrencyConv = false; erfile = indexFilePath * "EUR_USD_ExchangeRates.txt"
PPPConv = false; pppfile = indexFilePath * "PPP_ConvertingRates.txt"
ConstConv = true; const_year = 2010

codeSubst = true        # recommend 'false' for depth '1st' as there is nothing to substitute
perCap = true

gapMitigation = true    # filling gaps between national account and HBS expenditures

cpiScaling = false

printData = true

year = 2015
catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]
if codeSubst; substTag = "_subst" else substTag = "" end

microDataPath *= string(year) * "/"
# microDataPath = [microDataPath*"HU",microDataPath*"SE"]

ctgfile = extractedPath * string(year) * "_Category_"*depthTag[catDepth]*".csv"
hhsfile = extractedPath * string(year) * "_Households.csv"
mmsfile = extractedPath * string(year) * "_Members.csv"
expfile = extractedPath * string(year) * "_Expenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
sttfile = extractedPath * string(year) * "_MicroData_Statistics_"*depthTag[catDepth]*substTag*".csv"
sbcdsfile = extractedPath * string(year) * "_SubstituteCodes_" * depthTag[catDepth] * ".csv"
sbctgfile = extractedPath * string(year) * "_Category_" * depthTag[catDepth] * "_subst.csv"

scexpfile = extractedPath * string(year) * "_Scaled_Expenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
scstatsfile = extractedPath * string(year) * "_HBS_COICOP_stats_"*depthTag[catDepth]*substTag*".csv"

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(year, categoryFile, depth=catDepth, catFile=ctgfile, coicop=gapMitigation)
println("completed")

print(" Micro-data reading: ")
if readDataFromXLSX; println("XLSX")
    println("households"); mdr.readHouseholdData(year, microDataPath, visible=true, substitute=codeSubst, per=0.05, pea=1.0)
    println("members"); mdr.readMemberData(year, microDataPath, visible=true)
    print(", matrix"); mdr.buildExpenditureMatrix(year, substitute=codeSubst)
elseif readDataFromCSV; print("CSV")
    sttfile = replace(sttfile, ".csv"=>"_CSV.csv")
    print(" (hhs"); mdr.readPrintedHouseholdData(hhsfile)
    print(", mms"); mdr.readPrintedMemberData(mmsfile)
    if codeSubst; print(", subst"); mdr.readSubstCodesCSV(year, sbctgfile, sbcdsfile) end
    print(", exp)"); mdr.readPrintedExpenditureData(gapMitigation ? scexpfile : expfile, substitute=codeSubst, buildHhsExp=true)
end
print(", statistics"); mdr.makeStatistics(year, sttfile, substitute=codeSubst)
println(" ... complete")

if gapMitigation && readDataFromXLSX; print(" HBS-COICOP gap mitigating: ")
    mdr.mitigateExpGap(year, eustatsFile, percap=perCap, subst=codeSubst, cdrepl=true, alter=true, dist_mode = false)
    println("completed")
end

if CurrencyConv; print(" Currency exchanging: ")
    expfile = replace(expfile, ".csv"=>"_USD.csv")
    scexpfile = replace(scexpfile, ".csv"=>"_USD.csv")
    print(" exchange"); mdr.exchangeExpCurrency(erfile)
    print(", matrix"); mdr.buildExpenditureMatrix(year, substitute=codeSubst)
    println(" ... complete")
end
if PPPConv; print(" PPP converting: ")
    hhsfile = replace(hhsfile,".csv"=>"_PPP.csv")
    mmsfile = replace(mmsfile,".csv"=>"_PPP.csv")
    mdr.convertToPPP(pppfile)
    println("completed")
end
if CurrencyConv || PPPConv; mdr.makeStatistics(year, replace(sttfile,".csv"=>"_USD.csv")) end

if ConstConv || cpiScaling; print(" CPI reading: ")
    mdr.readCPIs([year, const_year], cpi_file, idx_sep = ',', freq="A", unit="INX_A_AVG", topLev = "EU")
    println("completed")
end

if ConstConv && year != const_year; print(" Price converting: to ", const_year," constant year price")
    hhsfile = replace(hhsfile, ".csv" => "_" * string(const_year) * "_constant.csv")
    mmsfile = replace(mmsfile, ".csv" => "_" * string(const_year) * "_constant.csv")
    mdr.convertToConstPrice(year, const_year)
    println(" ... completed")
end

if cpiScaling; print(" CPI scaling: ")
    mdr.scalingByCPI(year, const_year, codeDepth=0, topLev = "EU", subst = codeSubst)
    println("completed")
end

if printData; print(" Extracted data printing:")
    if readDataFromCSV; ctgfile, hhsfile, mmsfile, expfile = [replace(x, ".csv"=>"_CSV.csv") for x in [ctgfile, hhsfile, mmsfile, expfile]] end
    mdr.printCategory(year, ctgfile, substitute=codeSubst)
    mdr.printHouseholdData(year, hhsfile)
    mdr.printMemberData(year, mmsfile)
    if !gapMitigation;
        mdr.printExpenditureMatrix(year, expfile, substitute=codeSubst)
    elseif gapMitigation
        mdr.printExpTable(year, scexpfile, scaled=true, subst=codeSubst)
        mdr.printExpStats(year, scstatsfile, scaled=true)
    end
    println(" completed")
end

println("[done]")
