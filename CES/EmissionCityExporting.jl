# Developed date: 10. Mar. 2022
# Last modified date: 28. Mar. 2022
# Subject: Exporting City CF and CI web-files
# Description: Export CF and CI data by category for each city through analysis of
#               Customer Expenditure Survey (CES) or Household Budget Survey (HBS) micro-data.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
include("MicroDataReader.jl")
include("ConcMatBuilder.jl")
include("EmissionEstimator.jl")
include("EmissionCategorizer.jl")
include("EmissionCI.jl")
using .MicroDataReader
using .ConcMatBuilder
using .EmissionEstimator
using .EmissionCategorizer
using .EmissionCI
mdr = MicroDataReader
cmb = ConcMatBuilder
ee = EmissionEstimator
ec = EmissionCategorizer
ci = EmissionCI

currDict = Dict("IDN" => "IDR", "IND" => "INR")

# cesYear = 2018; exchYear = cesYear
# years = [cesYear]
# eoraYear = 2015     # eoraYear = cesYear
# natA3 = "IDN"; natCurr = currDict[natA3]
# quantMode = true

cesYear = 2011; exchYear = cesYear
years = [cesYear]
eoraYear = cesYear
natA3 = "IND"; natCurr = currDict[natA3]
quantMode = false

filePath = Base.source_dir() * "/data/" * natA3 * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissionPath = filePath * "emission/" * string(cesYear) * "/"
webIndexPath = indexFilePath * "web/"
gisPath = indexFilePath * "gis/"

scaleMode = false
curConv = false; curr_target = "USD"; erfile = indexFilePath * "CurrencyExchangeRates.txt"
pppConv = false; pppfile = filePath * "PPP_ConvertingRates.txt"

Qtable = "PRIMAP"

if curConv; currTag = curr_target else currTag = natCurr end
if scaleMode; scaleTag = "_Scaled" else scaleTag = "" end

regInfoFile = extractedPath * natA3 * "_" * string(cesYear) * "_RegionInfo.txt"
cmmfile = extractedPath * natA3 * "_" * string(cesYear) * "_Commodities.txt"
hhsfile = extractedPath * natA3 * "_" * string(cesYear) * "_Households_"*currTag*".txt"
mmsfile = extractedPath * natA3 * "_" * string(cesYear) * "_Members.txt"
itemfile = indexFilePath * natA3 * "_" * string(cesYear) * "_Commodity_items.txt"
expfile = extractedPath * natA3 * "_" * string(cesYear) * "_Expenditure_"*currTag*".txt"
exmfile = extractedPath * natA3 * "_" * string(cesYear) * scaleTag * "_Expenditure_matrix_"*currTag*".txt"

gisRegFile = gisPath * "regions.txt"
gisCatFile = gisPath * "cat_labels.txt"
gisConcFile = gisPath * "region_concordance.txt"

deFile = emissionPath * string(cesYear) * "_" * natA3 * "_hhs_" * scaleTag * "DE.txt"
ieFile = emissionPath * string(cesYear) * "_" * natA3 * "_hhs_" * scaleTag * "IE_" * Qtable * ".txt"

webIndexFile = webIndexPath * "keycode_index.txt"

Qtable = "PRIMAP"

categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
web_categories = ["FOOD", "ELECTRICITY", "GAS", "ENERGY", "PUBLIC_TRANS", "PRIVATE_TRANS", "MEDICAL",
                "EDUCATION", "CONSUMABLE", "DURABLE", "SERVICES", "ALL"]
subcat=""

cfav_file, cfac_file = Dict{Int, String}(), Dict{Int, String}()
for y in years
    cfav_file[y] = emissionPath * string(y) * "_" * natA3 * "_gis_" * subcat * "emission_cat_overall_CF_gr.csv"
    cfac_file[y] = emissionPath * string(y) * "_" * natA3 * "_gis_" * subcat * "emission_cat_dr_percap_CF_gr.csv"
end

incomePeriod = "annual"

ce_intgr_mode = "cf"

# expModes = ["ie", "de", "cf"]
# catMode = ["ie", "de", "cf"]
expModes = ["cf"]
catMode = ["cf"]

exceptCategory = ["None", "Taxes"]

city_file_sector = Array{Tuple{String, String}, 1}()
f = open(webIndexFile)
for l in eachline(f)
    s = string.(split(l, '\t'))
    push!(city_file_sector, (s[1], s[2]))
end
close(f)

println("[",cesYear,"]")

print(" Micro-data reading:")
print(" regions"); mdr.readPrintedRegionData(cesYear, natA3, regInfoFile)
print(", households"); mdr.readPrintedHouseholdData(cesYear, natA3, hhsfile)
print(", filtering"); mdr.filterRegionData(cesYear, natA3)
print(", sectors"); mdr.readPrintedSectorData(cesYear, natA3, cmmfile)
print(", expenditures"); mdr.readPrintedExpenditureData(cesYear, natA3, expfile, quantity=quantMode)
if curConv; print(", exchange"); mdr.exchangeExpCurrency(cesYear,exchYear,natA3,natCurr,erfile,target_curr=curr_target, hhs_exp=false, hhs_info=true) end
if pppConv; print(", ppp"); mdr.convertToPPP(cesYear, natA3, pppfile); println("complete") end
print(", matrix"); mes = mdr.buildExpenditureMatrix(cesYear, natA3, period = 365, quantity = quantMode)
println(" ... completed")

print(" Emission-data reading:")
print(" micro-data"); ec.importMicroData(mdr)
print(", DE"); ec.readEmissionData(cesYear, natA3, deFile, mode = "de")
print(", IE"); ec.readEmissionData(cesYear, natA3, ieFile, mode = "ie")
print(", CF"); ec.integrateCarbonFootprint()
print(", category"); ec.setCategory(cesYear, natA3, subgroup = "", except = exceptCategory)
print(", HHs"); for cm in catMode; print("_" * cm); ec.categorizeHouseholdEmission(cesYear, natA3, mode=cm, output="", hhsinfo=true) end
print(", Reg"); for cm in catMode; print("_" * cm); ec.categorizeRegionalEmission(cesYear, natA3, mode=cm, period="year", popwgh=true, region="district", ur=false, religion=false) end
print(", GIS_ID"); ec.readGISinfo(cesYear, natA3, gisRegFile, gisCatFile, id = true)
print(", GIS_conc"); ec.buildGISconc(cesYear, natA3, gisConcFile, region = "district", remove = true)
print(", filtering"); ec.filterRegion(cesYear, natA3, region = "district")
println(" ... completed")

web_city_path = emissionPath * "web/" * "footprint/"
web_center_path = emissionPath * "web/" * "centers/"

ci_rste = 0.95
n_iter = 10000

print(" Bootstrap process:")
print(" data import"); ci.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, cat_filter = true)
print(", CI calculation"); ci.estimateConfidenceIntervals(cesYear, natA3, iter = n_iter, ci_rate = ci_rste, resample_size = 0, replacement = true, boundary="district")
println(" ... completed")

print(" Web-file exporting:")
print(" center"); ec.exportCentersFile(cesYear, natA3, web_center_path)
print(", city"); ci.exportWebsiteCityFiles(cesYear, natA3, web_city_path, web_categories, city_file_sector, cfav_file, cfac_file, boundary="district")
println(" ... completed")
