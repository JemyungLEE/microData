# SPDX-FileCopyrightText: © 2021 Jemyung Lee <jemyung81@gmail.com>
# SPDX-License-Identifier: GPL-3.0

# Developed date: 31. Mar. 2021
# Last modified date: 15. Sep. 2022
# Subject: Household consumption expenditure survey microdata analysis
# Description: proceed microdata analysis process
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
include("MicroDataReader.jl")
using .MicroDataReader
mdr = MicroDataReader

currDict = Dict("IDN"=>"IDR", "IND"=>"INR", "VNM"=>"VND")

year = 2018; exchYear = year
eoraYear = 2015     # eoraYear = year
nation = "IDN"; natCurr = currDict[nation]
quantityMode = true
regionModify = false    # modify region code
skipTitle = true        # skip the first line (= title) of micro-data
dateType = "standard"     # "standard"=YYYYMMDD, "YYMMDD", "DDMMYY"

# year = 2011; exchYear = year
# eoraYear = year
# nation = "IND"; natCurr = currDict[nation]
# quantityMode = false
# regionModify = true     # modify region code
# skipTitle = false        # skip the first line (= title) of micro-data
# dateType = "DDMMYY"     # "standard"=YYYYMMDD, "YYMMDD", "DDMMYY"

filePath = Base.source_dir() * "/data/" * nation * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"

curConv = false; curr_target = "USD"
pppConv = false; pppfile = filePath * "PPP_ConvertingRates.txt"

gapMitigation = false    # filling gaps between national account and HBS expenditures
fitEoraYear = false      # scaling micro-data's expenditure to fit the Eora target year

# input files
erfile = filePath * nation * "_MD_ExchangeRate.txt"
# erfile = indexFilePath * "CurrencyExchangeRates.txt"

popfile = indexFilePath * nation * "_" * string(year) * "_Population.txt"
regfile = indexFilePath * nation * "_" * string(year) * "_Regions.txt"
hidxfile = indexFilePath * nation * "_" * string(year) * "_Households_microdata_attribute.txt"
eidxfile = indexFilePath * nation * "_" * string(year) * "_Expenditures_microdata_attribute.txt"
itemfile = indexFilePath * nation * "_" * string(year) * "_Commodity_items.txt"
# catfile = indexFilePath * nation * "_" * string(year) * "_Commodity_category.txt"
refrevfile = regionModify ? indexFilePath * nation * "_" * string(year) * "_Regions_modified.txt" : ""

# output files
cmmfile = extractedPath * nation * "_" * string(year) * "_MD_Commodities.txt"
hhsfile = extractedPath * nation * "_" * string(year) * "_MD_Households_"*natCurr*".txt"
mmsfile = extractedPath * nation * "_" * string(year) * "_MD_Members.txt"
exdfile = extractedPath * nation * "_" * string(year) * "_MD_Expenditure_"*natCurr*".txt"
regInfoFile = extractedPath * nation * "_" * string(year) * "_MD_RegionInfo.txt"

exmfile = extractedPath * nation * "_" * string(year) * "_MD_ExpenditureMatrix_"*natCurr*".txt"
sttfile = extractedPath * nation * "_" * string(year) * "_MD_Statistics.txt"

scexpfile = extractedPath * nation * "_" * string(year) * "_Scaled_ExpenditureMatrix_"*natCurr*".txt"

println("[Process]")
print(" Region information reading: ")
print("population"); mdr.readPopulation(year, nation, popfile)
print(", region"); mdr.readRegion(year, nation, regfile, region_revised_file = refrevfile)
println(" ... completed")

print(" Micro-data reading: ")
print("microdata"); mdr.readMicroData(year, nation, microDataPath, hidxfile, "", itemfile, eidxfile, hhid_sec = "hhid",
                                        skip_title = skipTitle, periodFiltering=true, ignoreException=true,
                                        region_modify = regionModify, visible = true, date_type = dateType)

if fitEoraYear && eoraYear != nothing && eoraYear != year; print(" Expenditure scaling: from $year to $eoraYear")
    exchYear = eoraYear
    cpiSecFile = indexFilePath * "CPI/CPI_" * nation * "_sectors.txt"
    statFile = indexFilePath * "CPI/CPI_" * nation * "_values.txt"
    linkFile = indexFilePath * "CPI/CPI_" * nation * "_link.txt"
    exmfile = replace(exmfile, ".txt"=>"_scaledTo"*string(eoraYear)*".txt")
    hhsfile = replace(hhsfile, ".txt"=>"_scaledTo"*string(eoraYear)*".txt")
    exdfile = replace(exdfile, ".txt"=>"_scaledTo"*string(eoraYear)*".txt")
    print(", scaling"); mdr.scalingExpByCPI(year, nation, cpiSecFile, statFile, linkFile, eoraYear, period="year", region="district", revHH=true, revMat=false)
    println(" ... completed")
end
if curConv
    print(", exchange"); mdr.exchangeExpCurrency(year, exchYear, nation, natCurr, erfile, target_curr=curr_target, hhs_info = true)
    exdfile = replace(exdfile, "_"*natCurr => "_"*curr_target)
    exmfile = replace(exmfile, "_"*natCurr => "_"*curr_target)
    scexpfile = replace(scexpfile, "_"*natCurr => "_"*curr_target)
    hhsfile = replace(hhsfile, "_"*natCurr => "_"*curr_target)
end
print(", matrix"); mes = mdr.buildExpenditureMatrix(year, nation, period = 365, quantity = quantityMode)
print(", weight"); mdr.calculatePopWeight(year, nation, district=true, province=false)
println(" ... completed")

if gapMitigation; print(" GDP-Survey gap mitigating: ")
    println(" ... completed")
end

print(" Extracted data printing:")
mdr.printCommoditySectors(year, nation, cmmfile)
mdr.printRegionData(year, nation, regInfoFile, region = "district", ur = false)
mdr.printHouseholdData(year, nation, hhsfile, hh_wgh = true, ur_dist = true, surv_date = true)
mdr.printExpenditureData(year, nation, exdfile, quantity = quantityMode)
mdr.printExpenditureMatrix(year, nation, exmfile, quantity = quantityMode, rowErr = mes[4], colErr = mes[5])
println(" ... completed")

println("[done]")
