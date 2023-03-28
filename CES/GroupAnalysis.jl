# Developed date: 17. Oct. 2022
# Last modified date: 28. Mar. 2023
# Subject: Analysis group information
# Description: Calculate statistics, CF, or related figures of grouped information
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()

cd(Base.source_dir())

include("MicroDataReader.jl")
include("EmissionCategorizer.jl")
include("GroupEstimator.jl")

using .MicroDataReader
using .EmissionCategorizer
using .GroupEstimator

mdr = MicroDataReader
ec = EmissionCategorizer
gre = GroupEstimator

cesYear = 2018; exchYear = cesYear
eoraYear = 2015
nation = "Indonesia"
natA3 = "IDN"
natCurr = "IDR"
hh_label = "URUT"
member_mode = false
keyMergMode = false    # set district code as "province_district"
resfile = "kor18rt_diseminasi_modified.csv"
qt_label = ["R1401","R1402","R1403","R1404","R1405","R1406","R1407","R1408","R1507","R1508","R1509",
            "R1510A","R1510B","R1511A","R1511B","R1512A"
            # ,"R1512B","R1512C"
            ,"R1513","R1516A","R1518A","R1519","R1520"]
qt_respo = ["1","1","1","1","1","1","1","1",["3","4","5","6","7","8"],["4","5","6","7"],["3","7","8","9"],
            ["3","4","5"],["3","4"],["6","7","8","9","10","11"],["2","8"],"2","1",["2","8"],"4",
            ["1","8","9","10","11","12"],"1"]


curConv = true; curr_target = "USD"
catMode = "cf"
exceptCategory = ["None", "Taxes"]

subcat=""
Qtable = "I_CHG_CO2"

scaleMode = false
quantMode = false

gisLabMode = true   # [true] use "GIS_name" ([false] use "City_name") in "GIS_RegionConc" for map city labeling
minSamples = 1      # minimum number of sample houses (include the value, >=)

incThres = [1.9, 3.0, 5.0]
incLabel = ["pov", "low", "mid", "high"]
incUnit = "USD/day/cap"

if scaleMode; scaleTag = "_Scaled" else scaleTag = "" end
filePath = Base.source_dir() * "/data/" * natA3 * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissionPath = filePath * "emission/" * string(cesYear) * "/"
commonIndexPath = Base.source_dir() * "/data/Common/"
gisIndexPath = commonIndexPath * "gis/"

natFileTag = natA3 * "_" * string(cesYear)
regInfoFile = filePath * natFileTag * "_MD_RegionInfo.txt"
cmmfile = filePath * natFileTag * "_MD_Commodities.txt"
hhsfile = filePath * natFileTag * "_MD_Households_"*natCurr*".txt"
mmsfile = filePath * natFileTag * "_MD_Members.txt"
erfile = filePath * natFileTag * "_MD_ExchangeRate.txt"
if !isfile(hhsfile); hhsfile = filePath * natFileTag * "_MD_Households.txt" end
if !isfile(erfile); erfile = filePath * natA3 * "_MD_ExchangeRate.txt" end
if !isfile(erfile); erfile = commonIndexPath * "CurrencyExchangeRates.txt" end

gisRegFile = filePath * natFileTag * "_GIS_RegionInfo.txt"
gisConcFile = filePath * natFileTag * "_GIS_RegionConc.txt"
gisCatFile = gisIndexPath * "category_labels.txt"

deFile = emissionPath * string(cesYear) * "_" * natA3 * "_hhs_" * scaleTag * "DE.txt"
ieFile = emissionPath * string(cesYear) * "_" * natA3 * "_hhs_" * scaleTag * "IE_" * Qtable * ".txt"

resfile = microDataPath * resfile
nat_grfile = extractedPath * natA3 * "_" * string(cesYear) * "_group_status.txt"
nat_resfile = extractedPath * natA3 * "_" * string(cesYear) * "_group_response.txt"
reg_sttfile = extractedPath * natA3 * "_" * string(cesYear) * "_region_status.txt"
reg_grfile = extractedPath * natA3 * "_" * string(cesYear) * "_region_group_status.txt"
reg_resfile = extractedPath * natA3 * "_" * string(cesYear) * "_region_group_response.txt"
reg_cmpfile = extractedPath * natA3 * "_" * string(cesYear) * "_region_group_compared.txt"

println("[Process]")

print(" Micro-data reading:")
print(" regions"); mdr.readExtractedRegionData(cesYear, natA3, regInfoFile, legacy_mode = true)
print(", households"); mdr.readExtractedHouseholdData(cesYear, natA3, hhsfile, legacy_mode = true)
if member_mode; print(", members"); mdr.readExtractedMemberData(cesYear, natA3, mmsfile) end
print(", population weight"); mdr.calculatePopWeight(cesYear, natA3, "", ur_wgh = false, district=true, province=false, hhs_wgh = true)
print(", sectors"); mdr.readExtractedSectorData(cesYear, natA3, cmmfile)
if curConv; print(", currency exchange"); mdr.exchangeExpCurrency(cesYear,exchYear,natA3,natCurr,erfile,target_curr=curr_target, exp_mat=false) end
println(" ... completed")

print(" Emission categorizing:")
rgCatFile = emissionPath * string(cesYear) * "_" * natA3 * "_region_categorized.txt"

print(" micro-data"); ec.importMicroData(mdr)
print(", DE"); ec.readEmissionData(cesYear, natA3, deFile, mode = "de")
print(", IE"); ec.readEmissionData(cesYear, natA3, ieFile, mode = "ie")
print(", CF"); ec.integrateCarbonFootprint()
print(", category"); ec.setCategory(cesYear, natA3, subgroup = "", except = exceptCategory)
print(", HHs_"*catMode); ec.categorizeHouseholdEmission(cesYear, natA3, mode=catMode, hhsinfo=true)

print(", GIS-info")
ec.readGISinfo(cesYear, natA3, gisRegFile, gisCatFile, id = true)
ec.buildGISconc(cesYear, natA3, gisConcFile, region = "district", remove = true, merged_key = keyMergMode, gis_label_mode = gisLabMode)
ec.filterRegion(cesYear, natA3; region = "district", limit = minSamples)
println(" ... completed")

print( " Region status estimation:")
print(" import"); gre.importData(mdr, ec, mode = "cf")
print(", threshold")
gre.setGroupThresholds(cesYear, natA3, incThres, incLabel, unit = incUnit)
gre.convertGroupThresholds(cesYear, natA3, conv_unit = "USD/year/cap")
print(", response"); gre.readResponseData(cesYear, natA3, resfile, qst_label = qt_label, qst_res = qt_respo, hhid_label = hh_label, weight_mode = true)
print(", estimation")
gre.estimateNationState(cesYear, natA3, group_mode = true, weight_mode = true, income_mode = false, category = "total")
gre.estimateRegionState(cesYear, natA3, group_mode = true, weight_mode = true, region_mode = "district", income_mode = false, category = "total")
print(", print")
gre.printNationalGroupStatus(cesYear, natA3, nat_grfile)
gre.printNationalGroupStatusByResponse(cesYear, natA3, nat_resfile)
gre.printRegionStatus(cesYear, natA3, reg_sttfile)
gre.printRegionalGroupStatus(cesYear, natA3, reg_grfile)
gre.printRegionalGroupStatusByResponse(cesYear, natA3, reg_resfile)
gre.printRegionalGroupCompared(cesYear, natA3, reg_cmpfile)
println(" ... completed")

println("[all complete]")
