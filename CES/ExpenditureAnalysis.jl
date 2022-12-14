# Developed date: 11. Nov. 2022
# Last modified date: 15. Nov. 2022
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

filePath = Base.source_dir() * "/data/" * natA3 * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissionPath = filePath * "emission/" * string(year) * "/"
commonIndexPath = Base.source_dir() * "/data/Common/"
gisIndexPath = commonIndexPath * "gis/"

curConv = false; curr_target = "USD"
pppConv = false; pppfile = filePath * "PPP_ConvertingRates.txt"

Qtable = "I_CHG_CO2"
# Qtable = "PRIMAP"

scaleMode = false
quantMode = false

exceptCategory = ["None", "Taxes"]

if scaleMode; scaleTag = "_Scaled" else scaleTag = "" end

natFileTag = natA3 * "_" * string(year)
regInfoFile = filePath * natFileTag * "_MD_RegionInfo.txt"
cmmfile = filePath * natFileTag * "_MD_Commodities_48.txt"
hhsfile = filePath * natFileTag * "_MD_Households_"*natCurr*".txt"
mmsfile = filePath * natFileTag * "_MD_Members.txt"
exmfile = filePath * natFileTag * "_MD_ExpenditureMatrix_"*natCurr*".txt"
erfile = filePath * natFileTag * "_MD_ExchangeRate.txt"
if !isfile(erfile); erfile = commonIndexPath * "CurrencyExchangeRates.txt" end

expfile = filePath * natFileTag * "_MD_Expenditure_"*natCurr*".txt"

exp_reg_file = extractedPath * natFileTag * "_exp_reg.txt"
exp_reg_cat_file = extractedPath * natFileTag * "_exp_reg_cat.txt"

println("[Process]")

print(" Micro-data reading:")
print(" regions"); mdr.readPrintedRegionData(year, natA3, regInfoFile, key_district = keyDistMode)
print(", households"); mdr.readPrintedHouseholdData(year, natA3, hhsfile)
print(", filtering"); mdr.filterRegionData(year, natA3)
if readMembers; print(", members"); mdr.readPrintedMemberData(year, natA3, mmsfile) end
print(", population weight"); mdr.calculatePopWeight(year, natA3, "", ur_wgh = false, district=true, province=false, hhs_wgh = true)
print(", sectors"); mdr.readPrintedSectorData(year, natA3, cmmfile)
if buildMatrix
    print(", expenditures"); mdr.readPrintedExpenditureData(year, natA3, expfile, quantity=quantMode)
    print(", matrix building"); mdr.buildExpenditureMatrix(year, natA3, period = 365, quantity = quantMode)
else print(", expenditure matrix"); mdr.readPrintedExpenditureMatrix(year, natA3, exmfile)
end
if curConv; print(", currency exchange"); mdr.exchangeExpCurrency(year,exchYear,natA3,natCurr,erfile,target_curr=curr_target, exp_mat=true) end
if pppConv; print(", ppp converting"); mdr.convertToPPP(year, natA3, pppfile); println("complete") end
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
