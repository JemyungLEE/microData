# Developed date: 21. May. 2021
# Last modified date: 28. Mar. 2023
# Subject: Categorized emission mapping
# Description: Mapping emission through households emissions data, categorizing by region, living-level, etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())

include("MicroDataReader.jl")
include("EmissionCategorizer.jl")
include("MapGenerator.jl")
include("../GIS/QgisStyleExporter.jl")
using .MicroDataReader
using .EmissionCategorizer
using .MapGenerator
using .QgisStyleExporter
mdr = MicroDataReader
ec = EmissionCategorizer
mg = MapGenerator
qse = QgisStyleExporter

# year = 2016; exchYear = year
# nation = "Vietnam"
# natA3 = "VNM"
# natCurr = "VND"
# readMembers = false     # read member data
# buildMatrix = false     # read expenditure data and build expenditure matrix
# keyDistMode = false      # set district code as key region code
# keyMergMode = false     # set district code as "province_district"
# groupMode = false     # seperate households by survey group

# year = 2018; exchYear = year
# nation = "Indonesia"
# natA3 = "IDN"
# natCurr = "IDR"
# readMembers = false     # read member data
# buildMatrix = true      # read expenditure data and build expenditure matrix
# keyDistMode = true      # set district code as key region code
# keyMergMode = false     # set district code as "province_district"
# groupMode = false     # seperate households by survey group

# year = 2011; exchYear = year
# nation = "India"
# natA3 = "IND"
# natCurr = "INR"
# readMembers = false     # read member data
# buildMatrix = true      # read expenditure data and build expenditure matrix
# keyDistMode = true      # set district code as key region code
# keyMergMode = false     # set district code as "province_district"# labelConvMode = false   # convert GeoJSON map's label from GIS_ID to GIS_label
# groupMode = false     # seperate households by survey group

# year = 2014; exchYear = year
# nation = "Japan"
# natA3 = "JPN"
# natCurr = "JPY"
# readMembers = false   # read member data
# buildMatrix = false   # read expenditure data and build expenditure matrix
# keyDistMode = true    # set district code as key region code
# keyMergMode = true    # set district code as "province_district"# labelConvMode = true  # convert GeoJSON map's label from GIS_ID to GIS_label
# groupMode = false     # seperate households by survey group

year = 2015; exchYear = year
nation = "United States"
natA3 = "USA"
natCurr = "USD"
readMembers = false     # read member data
buildMatrix = false     # read expenditure data and build expenditure matrix
keyDistMode = true      # set district code as key region code
keyMergMode = true      # set district code as "province_district"
groupMode = false        # seperate households by survey group

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

gisLabMode = true       # [true] use "GIS_name" ([false] use "City_name") in "GIS_RegionConc" for map city labeling
minSamples = 5          # minimum number of sample houses (include the value, >=)
labelConvMode = true    # convert GeoJSON map's label from GIS_ID to GIS_label
filterMode = true      # exclude regions that have fewere samples than 'minSamples'
skipNullHhs = true      # [true] exclude household that does not have district code

boundary_dict = Dict("IND" => [[[0,20000000]], []], "IDN" =>[[[0, 6000000]], []], "VNM" => [[[0,3000000]], []],
                    "JPN" => [[[0,7000000]], []], "USA" => [[], []])

exportMode = true; minmaxv = boundary_dict[natA3] # {{overall CF min., max.}, {CF per capita min., max.}
exportWebMode = false; unifiedIdMode = true
mapStyleMode = false; colormapReversePerCap=false; labeRevPerCap=true; colormapReverse=false; labeRev=false
mapGenMode = true   # generate GeoJSON maps

# expModes = ["ie", "de", "cf"]
# catMode = ["ie", "de", "cf"]
expModes = ["cf"]
catMode = ["cf"]

exceptCategory = ["None", "Taxes"]

subcat=""
# subcat="Food"

if scaleMode; scaleTag = "_Scaled" else scaleTag = "" end

natFileTag = "source/" * string(year) * "/" * natA3 * "_" * string(year)
# natFileTag = natA3 * "_" * string(year)
regInfoFile = filePath * natFileTag * "_MD_RegionInfo.txt"
cmmfile = filePath * natFileTag * "_MD_Commodities.txt"
hhsfile = filePath * natFileTag * "_MD_Households_"*natCurr*".txt"
mmsfile = filePath * natFileTag * "_MD_Members.txt"
exmfile = filePath * natFileTag * "_MD_ExpenditureMatrix_"*natCurr*".txt"
erfile = filePath * natFileTag * "_MD_ExchangeRate.txt"
if !isfile(hhsfile); hhsfile = filePath * natFileTag * "_MD_Households.txt" end
if !isfile(exmfile); exmfile = filePath * natFileTag * "_MD_Expenditure.txt" end
if !isfile(erfile); erfile = filePath * natA3 * "_MD_ExchangeRate.txt" end
if !isfile(erfile); erfile = commonIndexPath * "CurrencyExchangeRates.txt" end

expfile = filePath * natFileTag * "_MD_Expenditure_"*natCurr*".txt"

gisRegFile = filePath * natFileTag * "_GIS_RegionInfo.txt"
gisConcFile = filePath * natFileTag * "_GIS_RegionConc.txt"
gisCatFile = gisIndexPath * "category_labels.txt"

deFile = emissionPath * string(year) * "_" * natA3 * "_hhs_" * scaleTag * "DE.txt"
ieFile = emissionPath * string(year) * "_" * natA3 * "_hhs_" * scaleTag * "IE_" * Qtable * ".txt"

basemapFile = filePath * natFileTag * ".geojson"
mapListFile = gisIndexPath * "Map_filenames.txt"
mapFilePath = filePath * "maps/" * string(year) * "/"
rgbfile_pc = gisIndexPath * "MPL_RdBu.rgb"
rgbfile_ov = gisIndexPath * "MPL_YlGnBu.rgb"

println("[Process]")

print(" Micro-data reading:")
print(" regions"); mdr.readExtractedRegionData(year, natA3, regInfoFile, key_district = keyDistMode, merged_key = keyMergMode, legacy_mode = true)
print(", households"); mdr.readExtractedHouseholdData(year, natA3, hhsfile, merged_key = keyMergMode, skip_empty = skipNullHhs, legacy_mode = true)
if filterMode; print(", filtering"); mdr.filterRegionData(year, natA3) end
print(", find lost"); mdr.findLostRegion(year,natA3)
if readMembers; print(", members"); mdr.readExtractedMemberData(year, natA3, mmsfile) end
print(", population weight"); mdr.calculatePopWeight(year, natA3, "", ur_wgh = false, district=true, province=false, hhs_wgh = true)
print(", sectors"); mdr.readExtractedSectorData(year, natA3, cmmfile)
if buildMatrix
    print(", expenditures"); mdr.readExtractedExpenditureData(year, natA3, expfile, quantity=quantMode)
    print(", matrix building"); mdr.buildExpenditureMatrix(year, natA3, period = 365, quantity = quantMode)
else print(", expenditure matrix"); mdr.readExtractedExpenditureMatrix(year, natA3, exmfile)
end
if curConv; print(", currency exchange"); mdr.exchangeExpCurrency(year,exchYear,natA3,natCurr,erfile,target_curr=curr_target, exp_mat=true) end
if pppConv; print(", ppp converting"); mdr.convertToPPP(year, natA3, pppfile); println("complete") end
println(" ... completed")

print(" Emission categorizing:")
rgCatFile = emissionPath * string(year) * "_" * natA3 * "_region_categorized.txt"
if groupMode; rgCatGrFile = emissionPath * string(year) * "_" * natA3 * "_region_categorized_grouped.txt" end

print(" micro-data"); ec.importMicroData(mdr)
print(", DE"); ec.readEmissionData(year, natA3, deFile, mode = "de")
print(", IE"); ec.readEmissionData(year, natA3, ieFile, mode = "ie")
print(", CF"); ec.integrateCarbonFootprint()
print(", category"); ec.setCategory(year, natA3, subgroup = "", except = exceptCategory)
for cm in catMode
    hhCatFile = emissionPath * string(year) * "_" * natA3 * "_hhs_"*uppercase(cm)*"_categorized.txt"
    print(", HHs_"*cm); ec.categorizeHouseholdEmission(year, natA3, mode=cm, output=hhCatFile, hhsinfo=true)
    print(", Reg_"*cm); ec.categorizeRegionalEmission(year, natA3, mode=cm, period="year", popwgh=true, region="district", ur=false, religion=false, group=groupMode)
end
print(", printing"); ec.printRegionalEmission(year, natA3, rgCatFile, region="district", mode=catMode, popwgh=true, ur=false, religion=false)
if groupMode; ec.printRegionalGroupEmission(year, natnatA3ion, rgCatGrFile, region="district", mode=catMode, popwgh=true, ur=false, gr=groupMode, religion=false) end
println(" ... completed")

print(" Exporting: ")
if exportMode || exportWebMode || mapStyleMode || mapGenMode;
    print(" GIS-info")
    ec.readGISinfo(year, natA3, gisRegFile, gisCatFile, id = unifiedIdMode)
    ec.buildGISconc(year, natA3, gisConcFile, region = "district", remove = true, merged_key = keyMergMode, gis_label_mode = gisLabMode)
    ec.filterRegion(year, natA3; region = "district", limit = minSamples)

    print(", GIS-exporting")
    gisTag = "District"
    exportFile = emissionPath * "YEAR_"*natA3*"_gis_"*subcat*"emission_cat_OvPcTag.csv"
    exportRateFile = emissionPath * "YEAR_"*natA3*"_gis_"*subcat*"emission_cat_dr_OvPcTag.csv"
    labelList, labelListPerCap = ec.exportRegionalEmission(year, natA3, gisTag, exportFile, region="district", mode=expModes,  nspan=128, minmax=minmaxv, descend=true, empty=false, logarithm=false)
    spanVals, spanValsPerCap = ec.exportEmissionDevRate(year, natA3, gisTag, exportRateFile, mode=expModes, maxr=0.5, minr=-0.5, nspan=128, descend=true, empty=false)
end
if exportWebMode; print(", web-files")
    exportPath = filePath * "webfile/" * string(year) * "/"
    mkpath(exportPath)
    ec.exportWebsiteFiles(year, natA3, exportPath, mode=expModes, rank=true, empty=false)
end
if mapStyleMode; print(", map-style file generating")
    rgbPath = indexFilePath * "gis/rgb/"
    mkpath(rgbPath)
    rgbFile = "../GIS/data/MPL_RdBu.rgb"
    for i=1:length(ec.cat_list)
        qse.readColorMap(rgbFile, reverse=colormapReversePerCap)
        qmlFile = rgbPath * replace(rsplit(rgbFile, '/', limit=2)[end], ".rgb"=>"_percap_"*ec.cat_list[i]*".qml")
        attr = string(year)*"_"*natA3*"_gis_"*subcat*"emission_cat_dr_percap_CF_gr_"*ec.cat_list[i]
        qse.makeQML(qmlFile, attr, empty=false, values=spanValsPerCap[year][natA3][:,i], indexValue=true, labelReverse=labeRevPerCap)
    end

    rgbFile = "../GIS/data/MPL_YlGnBu.rgb"
    for i=1:length(ec.cat_list)
        qse.readColorMap(rgbFile, reverse=colormapReverse)
        qmlFile = rgbPath * replace(rsplit(rgbFile, '/', limit=2)[end], ".rgb"=>"_overall_"*ec.cat_list[i]*".qml")
        attr = string(year)*"_"*natA3*"_gis_"*subcat*"emission_cat_overall_CF_gr_"*ec.cat_list[i]
        qse.makeQML(qmlFile, attr, empty=false, labels=labelList[year][natA3][:,i], indexValue=true, labelReverse=labeRev)
    end
end
if mapGenMode; print(", map-generation")
    mg.importEmissionData(ec, emission = "cf", pc_dev = true, ov_dev = false)
    mg.readBaseMap(year, natA3, basemapFile, remove_reg = false, alter = true, label_conv = labelConvMode)
    mg.readMapInfo(mapListFile)
    mg.convertRgbToHex(mg.readColorMap(rgbfile_ov, reverse=false), mode = "overall")
    mg.convertRgbToHex(mg.readColorMap(rgbfile_pc, reverse=false), mode = "percap")
    mg.mapRegionCF(year, natA3, label_conv = true, blank_color = "#EDEDED", value_mode = true)
    mg.printMapFiles(year, natA3, mapFilePath)
end
println(" ... completed")

println("[all complete]")
