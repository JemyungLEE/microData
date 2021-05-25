# Developed date: 21. May. 2021
# Last modified date: 25. May. 2021
# Subject: Categorized emission mapping
# Description: Mapping emission through households emissions data, categorizing by region, living-level, etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)


### Note: HH currency, EXP currency, period setting parts needed


cd(Base.source_dir())

include("MicroDataReader.jl")
include("EmissionCategorizer.jl")

using .MicroDataReader
using .EmissionCategorizer

mdr = MicroDataReader
ec = EmissionCategorizer

year = 2018; exchYear = year
# nation = "Indonesia"
natA3 = "IDN"
natCurr = "IDR"
readMembers = false     # read member data
readExpends = false      # read expenditure data
buildExpMat = false      # build expenditure matrix

filePath = Base.source_dir() * "/data/" * natA3 * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissionPath = filePath * "emission/"

curConv = true; curr_target = "USD"; erfile = indexFilePath * "CurrencyExchangeRates.txt"
pppConv = false; pppfile = filePath * "PPP_ConvertingRates.txt"

# Qtable = "I_CHG_CO2"
Qtable = "PRIMAP"

quantMode = false
printMode = false

regInfoFile = extractedPath * natA3 * "_" * string(year) * "_RegionInfo.txt"
cmmfile = extractedPath * natA3 * "_" * string(year) * "_Commodities.txt"
hhsfile = extractedPath * natA3 * "_" * string(year) * "_Households_"*natCurr*".txt"
mmsfile = extractedPath * natA3 * "_" * string(year) * "_Members.txt"
itemfile = indexFilePath * natA3 * "_" * string(year) * "_Commodity_items.txt"
expfile = extractedPath * natA3 * "_" * string(year) * "_Expenditure_"*natCurr*".txt"
exmfile = extractedPath * natA3 * "_" * string(year) * "_Expenditure_matrix_"*natCurr*".txt"

deFile = emissionPath * string(year) * "_" * natA3 * "_hhs_"*"DE_org.txt"
ieFile = emissionPath * string(year) * "_" * natA3 * "_hhs_"*"IE_"*Qtable*"_org.txt"

println("[Process]")

print(" Micro-data reading:")
print(" regions"); mdr.readPrintedRegionData(year, natA3, regInfoFile)
print(", households"); mdr.readPrintedHouseholdData(year, natA3, hhsfile)
if readMembers; print(", members"); mdr.readPrintedMemberData(year, natA3, mmsfile) end
print(", sectors"); mdr.readPrintedSectorData(year, natA3, cmmfile)
if readExpends; print(", expenditures"); mdr.readPrintedExpenditureData(year, natA3, expfile, quantity=quantMode) end
if curConv; print(", exchange"); mdr.exchangeExpCurrency(year,exchYear,natA3,natCurr,erfile,target_curr=curr_target, hhs_exp=false, hhs_info=true) end
if pppConv; print(", ppp"); mdr.convertToPPP(year, natA3, pppfile); println("complete") end
if buildExpMat; print(", matrix"); mes = mdr.buildExpenditureMatrix(year, natA3, period = 365, quantity = quantMode) end
println(" ... completed")

if printMode
    print_tag = "_test"
    mdr.printCommoditySectors(year, natA3, replace(cmmfile, ".txt"=>print_tag*".txt") )
    mdr.printRegionData(year, natA3, replace(regInfoFile, ".txt"=>print_tag*".txt") , region = "district", ur = false)
    mdr.printHouseholdData(year, natA3, replace(hhsfile, ".txt"=>print_tag*".txt") , prov_wgh=false, dist_wgh=true, ur_dist=false, surv_date = false)
    if readExpends; mdr.printExpenditureData(year, natA3, replace(expfile, ".txt"=>print_tag*".txt") , quantity = true) end
    if buildExpMat; mdr.printExpenditureMatrix(year, natA3, replace(exmfile, ".txt"=>print_tag*".txt") , rowErr = mes[4], colErr = mes[5]) end
end

print(" Emission categorizing:")
catMode = ["de", "ie", "cf"]
rgCatFile = emissionPath * string(year) * "_" * natA3 * "_region_categorized.txt"

print(" micro-data"); ec.importMicroData(mdr)
print(", DE"); ec.readEmissionData(year, natA3, deFile, mode = "de")
print(", IE"); ec.readEmissionData(year, natA3, ieFile, mode = "ie")
print(", CF"); ec.integrateCarbonFootprint(mode="cf")
print(", category"); ec.setCategory(year, natA3, subgroup = "", except=["None"])
for cm in catMode
    hhCatFile = emissionPath * string(year) * "_" * natA3 * "_hhs_"*uppercase(cm)*"_categorized.txt"
    print(", HHs_"*cm); ec.categorizeHouseholdEmission(year, natA3, mode=cm, output=hhCatFile, hhsinfo=true)
end
for cm in catMode
    print(", Reg_"*cm); ec.categorizeRegionalEmission(year, natA3, mode=cm, period="year", popwgh=true, region="district", ur=false, religion=false)
end
print(", printing"); ec.printRegionalEmission(year, natA3, rgCatFile, region="district", mode=["cf","de","ie"], popwgh=true, ur=false, religion=false)
println(" ... completed")

println("[all complete]")
