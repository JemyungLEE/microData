# SPDX-FileCopyrightText: © 2022 Jemyung Lee <jemyung81@gmail.com>
# SPDX-License-Identifier: GPL-3.0

# Developed date: 11. Nov. 2022
# Last modified date: 6. Sep. 2023
# Subject: Categorized emission mapping
# Description: Mapping emission through households emissions data, categorizing by region, living-level, etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())

include("MicroDataReader.jl")
include("EmissionCategorizer.jl")
include("../GIS/QgisStyleExporter.jl")
using .MicroDataReader
using .EmissionCategorizer
using .QgisStyleExporter
mdr = MicroDataReader
ec = EmissionCategorizer
qse = QgisStyleExporter

# year = 2016; exchYear = year
# nation = "Vietnam"
# natA3 = "VNM"
# natCurr = "VND"
# readMembers = false     # read member data
# buildMatrix = false     # read expenditure data and build expenditure matrix
# keyDistMode = true      # set district code as key region code

year = 2018; exchYear = year
nation = "Indonesia"
natA3 = "IDN"
natCurr = "IDR"
readMembers = false     # read member data
buildMatrix = true      # read expenditure data and build expenditure matrix
keyDistMode = true      # set district code as key region code

# year = 2011; exchYear = year
# nation = "India"
# natA3 = "IND"
# natCurr = "INR"
# readMembers = false     # read member data
# buildMatrix = true      # read expenditure data and build expenditure matrix
# keyDistMode = true      # set district code as key region code

# year = 2021; exchYear = year
# nation = "Chinese Taipei"
# natA3 = "TWN"
# natCurr = "TWD"
# readMembers = false     # read member data
# buildMatrix = true      # read expenditure data and build expenditure matrix
# keyDistMode = true      # set district code as key region code

filePath = Base.source_dir() * "/data/" * natA3 * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissionPath = filePath * "emission/" * string(year) * "/"
commonIndexPath = Base.source_dir() * "/data/Common/"
gisIndexPath = commonIndexPath * "gis/"

curConv = true; curr_target = "USD"
pppConv = false; pppfile = filePath * "PPP_ConvertingRates.txt"

Qtable = "I_CHG_CO2"
# Qtable = "PRIMAP"

scaleMode = false
quantMode = false

exceptCategory = ["None", "Taxes"]

if scaleMode; scaleTag = "_Scaled" else scaleTag = "" end
natFileTag = "source/" * string(year) * "/" * natA3 * "_" * string(year)
regInfoFile = filePath * natFileTag * "_MD_RegionInfo.txt"

regInfoFile = filePath * natFileTag * "_MD_RegionInfo.txt"
cmmfile = filePath * natFileTag * "_MD_Commodities.txt"
hhsfile = filePath * natFileTag * "_MD_Households.txt"
mmsfile = filePath * natFileTag * "_MD_Members.txt"
exmfile = filePath * natFileTag * "_MD_Expenditure.txt"
erfile = filePath * natFileTag * "_MD_ExchangeRate.txt"
if !isfile(erfile); erfile = commonIndexPath * "CurrencyExchangeRates.txt" end

expfile = filePath * natFileTag * "_MD_Expenditure_"*natCurr*".txt"

exp_reg_file = extractedPath * natA3 * "_" * string(year) * "_exp_reg.txt"
exp_reg_cat_file = extractedPath * natA3 * "_" * string(year) * "_exp_reg_cat.txt"

println("[Process]")

print(" Micro-data reading:")
print(" "); mdr.readExtractedRegionData(year, natA3, regInfoFile, key_district = keyDistMode, merged_key = true, legacy_mode = true, ignore = false, remove_empty = true)
print(", "); mdr.readExtractedHouseholdData(year, natA3, hhsfile, merged_key = true, skip_empty = true, legacy_mode = true)
print(", "); mdr.readExtractedSectorData(year, natA3, cmmfile)
if readMembers; print(", "); mdr.readExtractedMemberData(year, natA3, mmsfile) end
print(", population weight"); mdr.calculatePopWeight(year, natA3, "", ur_wgh = false, district=true, province=false, hhs_wgh = true)
print(", expenditure"); mdr.readExtractedExpenditureMatrix(year, natA3, exmfile, quantity = quantMode)
print(", filtering"); mdr.filterData(year, natA3, group=false, region="district", quantity=quantMode)
if curConv; print(", currency exchange"); mdr.exchangeExpCurrency(year,exchYear,natA3,natCurr,erfile,target_curr=curr_target, exp_mat=true) end
if pppConv; print(", ppp converting"); mdr.convertToPPP(year, natA3, pppfile); println("complete") end
print(", reshaping"); mdr.reshapeCommoditySectors(year, natA3, except = exceptCategory)
print(", find lost"); mdr.findLostRegion(year,natA3)
println(" ... completed")

print(" Emission categorizing:")

print(" micro-data"); ec.importMicroData(mdr)
print(", category"); ec.setCategory(year, natA3, subgroup = "", except = exceptCategory)
ec.estimateRegionalExpenditure(year, natA3, cat_mode = true, item_mode = true, period="year", popwgh=true, region = "district", qnt_mode = false)
print(", printing")
ec.printRegionalExpenditure(year, natA3, exp_reg_file, region = "district", mode = "item", popwgh=true, ur=false)
ec.printRegionalExpenditure(year, natA3, exp_reg_cat_file, region = "district", mode = "category", popwgh=true, ur=false)
println(" ... completed")

println("[all complete]")
