# SPDX-FileCopyrightText: © 2022 Jemyung Lee <jemyung81@gmail.com>
# SPDX-License-Identifier: GPL-3.0

# Developed date: 10. Mar. 2022
# Last modified date: 15. Jun. 2023
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

currDict = Dict("IDN" => "IDR", "IND" => "INR", "VNM" => "VND", "JPN" => "JPY", "USA" => "USD")

# cesYear = 2016; exchYear = cesYear
# years = [cesYear]
# eoraYear = cesYear
# natA3 = "VNM"; natCurr = currDict[natA3]
# quantMode = false
# readMatrix = true   # [true]: read expenditure matrix, [false]: read expenditure data and build expenditure matrix
# keyDistMode = true  # set district code as key region code
# keyMergMode = false     # set district code as "province_district"
# groupMode = false        # seperate households by survey group

# cesYear = 2018; exchYear = cesYear
# years = [cesYear]
# eoraYear = 2015     # eoraYear = cesYear
# natA3 = "IDN"; natCurr = currDict[natA3]
# quantMode = false
# readMatrix = false
# keyDistMode = true  # set district code as key region code
# keyMergMode = false     # set district code as "province_district"
# groupMode = false        # seperate households by survey group

# cesYear = 2011; exchYear = cesYear
# years = [cesYear]
# eoraYear = cesYear
# natA3 = "IND"; natCurr = currDict[natA3]
# quantMode = false
# readMatrix = false
# keyDistMode = true  # set district code as key region code
# keyMergMode = false     # set district code as "province_district"
# groupMode = false        # seperate households by survey group

# cesYear = 2014; exchYear = cesYear
# years = [cesYear]
# eoraYear = cesYear
# natA3 = "JPN"; natCurr = currDict[natA3]
# quantMode = false
# readMatrix = true
# keyDistMode = true      # set district code as key region code
# keyMergMode = true     # set district code as "province_district"
# groupMode = false        # seperate households by survey group

cesYear = 2008; exchYear = cesYear
years = [cesYear]
eoraYear = cesYear
natA3 = "USA"; natCurr = currDict[natA3]
quantMode = false
readMatrix = true
keyDistMode = true      # set district code as key region code
keyMergMode = true     # set district code as "province_district"
groupMode = true        # seperate households by survey group
groupSplit = false

commPath = Base.source_dir() * "/data/Common/"
filePath = Base.source_dir() * "/data/" * natA3 * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
emissionPath = filePath * "emission/" * string(cesYear) * "/"
webIndexPath = commPath * "web/"
gisPath = commPath * "gis/"

gisLabMode = true     # [true] use "GIS_name" ([false] use "City_name") in "GIS_RegionConc" for map city labeling
scaleMode = false
curConv = false; curr_target = "USD"; erfile = indexFilePath * "CurrencyExchangeRates.txt"
pppConv = false; pppfile = filePath * "PPP_ConvertingRates.txt"

skipNullHhs = true      # [true] exclude household that does not have district code

Qtable = "I_CHG_CO2"; q_tag = "_i_chg_co2"
# Qtable = "PRIMAP"; q_tag = "_primap"

minSamples = 5  # minimum number of sample houses (include the value, >=)
filterMode = true      # exclude regions that have fewere samples than 'minSamples'

if curConv; currTag = curr_target else currTag = natCurr end
if scaleMode; scaleTag = "_Scaled" else scaleTag = "" end

natFileTag = "source/" * string(cesYear) * "/" * natA3 * "_" * string(cesYear)
# natFileTag = natA3 * "_" * string(cesYear)
regInfoFile = filePath * natFileTag * "_MD_RegionInfo.txt"
cmmfile = filePath * natFileTag * "_MD_Commodities.txt"
hhsfile = filePath * natFileTag * "_MD_Households_"*currTag*".txt"
mmsfile = filePath * natFileTag * "_MD_Members.txt"
expfile = filePath * natFileTag * "_MD_Expenditure_"*currTag*".txt"
exmfile = filePath * natFileTag * scaleTag * "_MD_ExpenditureMatrix_"*currTag*".txt"
if !isfile(hhsfile); hhsfile = filePath * natFileTag * "_MD_Households.txt" end
if !isfile(exmfile); exmfile = filePath * natFileTag * "_MD_Expenditure.txt" end
if !isfile(erfile); erfile = filePath * natA3 * "_MD_ExchangeRate.txt" end
if !isfile(erfile); erfile = commPath * "CurrencyExchangeRates.txt" end

web_city_path = filePath * "web/" * "footprint/"
web_center_path = filePath * "web/" * "centers/"

ci_rste = 0.95
n_iter = 10000

gisCatFile = gisPath * "category_labels.txt"
gisRegFile = filePath * natFileTag * "_GIS_RegionInfo.txt"
gisConcFile = filePath * natFileTag * "_GIS_RegionConc.txt"

deFile = emissionPath * string(cesYear) * "_" * natA3 * "_hhs_" * scaleTag * "DE.txt"
ieFile = emissionPath * string(cesYear) * "_" * natA3 * "_hhs_" * scaleTag * "IE" * q_tag * ".txt"

webIndexFile = webIndexPath * "keycode_index.txt"
# webIndexFile = webIndexPath * "keycode_index_all.txt"

ces_categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
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

println("[",cesYear,"]")

print(" Micro-data reading:")
print(" regions"); mdr.readExtractedRegionData(cesYear, natA3, regInfoFile, key_district = keyDistMode, merged_key = keyMergMode, legacy_mode = true)
print(", sectors"); mdr.readExtractedSectorData(cesYear, natA3, cmmfile)
print(", households"); mdr.readExtractedHouseholdData(cesYear, natA3, hhsfile, merged_key = keyMergMode, skip_empty = skipNullHhs, legacy_mode = true)
print(", find lost"); mdr.findLostRegion(cesYear, natA3)
if readMatrix
    print(", expenditure matrix"); mdr.readExtractedExpenditureMatrix(cesYear, natA3, exmfile)
else
    print(", expenditures"); mdr.readExtractedExpenditureData(cesYear, natA3, expfile, quantity=quantMode)
    print(", matrix building"); mdr.buildExpenditureMatrix(cesYear, natA3, period = 365, quantity = quantMode)
end
if filterMode; print(", filtering"); mdr.filterData(cesYear, natA3, group=groupMode, region="district", quantity=quantMode) end
print(", population weight"); mdr.calculatePopWeight(cesYear, natA3, "", ur_wgh = false, district=true, province=false, hhs_wgh = true, gr_wgh = groupMode)
if curConv; print(", currency exchange"); mdr.exchangeExpCurrency(cesYear,exchYear,natA3,natCurr,erfile,target_curr=curr_target, exp_mat=true, hhs_info=true) end
if pppConv; print(", ppp converting"); mdr.convertToPPP(cesYear, natA3, pppfile); println("complete") end
print(", reshape commodities"); mdr.reshapeCommoditySectors(cesYear, natA3, except = exceptCategory, hhs_reshape = !readMatrix)
println(" ... completed")

print(" Emission-data reading:")
print(" micro-data"); ec.importMicroData(mdr)
print(", DE"); ec.readEmissionData(cesYear, natA3, deFile, mode = "de", revise = true)
print(", IE"); ec.readEmissionData(cesYear, natA3, ieFile, mode = "ie", revise = true)
if groupMode
    if groupSplit; print(", split groups"); ec.splitHouseholdGroup(cesYear, natA3, mode = ["ie","de"]) end
    print(", group"); ec.filterGroupEmission(cesYear, natA3, mode = ["ie","de"])
end
print(", CF"); ec.integrateCarbonFootprint()
print(", category"); ec.setCategory(cesYear, natA3, categories = ces_categories, subgroup = "", except = exceptCategory)
print(", HHs"); for cm in catMode; print("_" * cm); ec.categorizeHouseholdEmission(cesYear, natA3, mode=cm, output="", hhsinfo=true, group = groupMode) end
print(", Reg"); for cm in catMode; print("_" * cm); ec.categorizeRegionalEmission(cesYear, natA3, mode=cm, period="year", popwgh=true, region="district", ur=false, religion=false, group=groupMode) end
print(", GIS_ID"); ec.readGISinfo(cesYear, natA3, gisRegFile, gisCatFile, id = true)
print(", GIS_conc"); ec.buildGISconc(cesYear, natA3, gisConcFile, region = "district", remove = true, merged_key = keyMergMode, gis_label_mode = gisLabMode)
print(", filtering"); ec.filterRegion(cesYear, natA3, region = "district", limit = minSamples, group = groupMode)
println(" ... completed")

print(" Bootstrap process:")
print(" data import"); ci.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, cat_filter = true)
print(", CI calculation"); ci.estimateConfidenceIntervals(cesYear, natA3, iter = n_iter, ci_rate = ci_rste, resample_size = 0, replacement = true, boundary="district", group = groupMode)
println(" ... completed")

print(" Web-file exporting:")
print(" center"); ec.exportCentersFile(cesYear, natA3, web_center_path)
print(", web index"); ci.readCityFileSector(webIndexFile)
print(", city"); ci.exportWebsiteCityFiles(cesYear, natA3, web_city_path, web_categories, cfav_file, cfac_file, boundary="district")
println(" ... completed")

println("[all complete]")
