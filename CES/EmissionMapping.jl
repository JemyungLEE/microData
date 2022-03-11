# Developed date: 21. May. 2021
# Last modified date: 11. Mar. 2022
# Subject: Categorized emission mapping
# Description: Mapping emission through households emissions data, categorizing by region, living-level, etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

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

year = 2018; exchYear = year
# nation = "Indonesia"
natA3 = "IDN"
natCurr = "IDR"
readMembers = false     # read member data
readExpends = true      # read expenditure data
buildExpMat = false      # build expenditure matrix

filePath = Base.source_dir() * "/data/" * natA3 * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissionPath = filePath * "emission/" * string(year) * "/"

curConv = true; curr_target = "USD"; erfile = indexFilePath * "CurrencyExchangeRates.txt"
pppConv = false; pppfile = filePath * "PPP_ConvertingRates.txt"

# Qtable = "I_CHG_CO2"
Qtable = "PRIMAP"

scaleMode = false
quantMode = true
printMode = false
exportMode = true; minmaxv = [[[5000, 6000000]], []] # {{overall CF min., max.}, {CF per capita min., max.}
exportWebMode = true; unifiedIdMode = true
mapStyleMode = true; colormapReversePerCap=false; labeRevPerCap=true; colormapReverse=false; labeRev=false

# expModes = ["ie", "de", "cf"]
# catMode = ["ie", "de", "cf"]
expModes = ["cf"]
catMode = ["cf"]

exceptCategory = ["None", "Taxes"]

subcat=""
# subcat="Food"

if scaleMode; scaleTag = "_Scaled" else scaleTag = "" end

regInfoFile = extractedPath * natA3 * "_" * string(year) * "_RegionInfo.txt"
cmmfile = extractedPath * natA3 * "_" * string(year) * "_Commodities.txt"
hhsfile = extractedPath * natA3 * "_" * string(year) * "_Households_"*natCurr*".txt"
mmsfile = extractedPath * natA3 * "_" * string(year) * "_Members.txt"
itemfile = indexFilePath * natA3 * "_" * string(year) * "_Commodity_items.txt"
expfile = extractedPath * natA3 * "_" * string(year) * "_Expenditure_"*natCurr*".txt"
exmfile = extractedPath * natA3 * "_" * string(year) * scaleTag * "_Expenditure_matrix_"*natCurr*".txt"

deFile = emissionPath * string(year) * "_" * natA3 * "_hhs_" * scaleTag * "DE.txt"
ieFile = emissionPath * string(year) * "_" * natA3 * "_hhs_" * scaleTag * "IE_" * Qtable * ".txt"

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
rgCatFile = emissionPath * string(year) * "_" * natA3 * "_region_categorized.txt"

print(" micro-data"); ec.importMicroData(mdr)
print(", DE"); ec.readEmissionData(year, natA3, deFile, mode = "de")
print(", IE"); ec.readEmissionData(year, natA3, ieFile, mode = "ie")
print(", CF"); ec.integrateCarbonFootprint()
print(", category"); ec.setCategory(year, natA3, subgroup = "", except = exceptCategory)
for cm in catMode
    hhCatFile = emissionPath * string(year) * "_" * natA3 * "_hhs_"*uppercase(cm)*"_categorized.txt"
    print(", HHs_"*cm); ec.categorizeHouseholdEmission(year, natA3, mode=cm, output=hhCatFile, hhsinfo=true)
end
for cm in catMode
    print(", Reg_"*cm); ec.categorizeRegionalEmission(year, natA3, mode=cm, period="year", popwgh=true, region="district", ur=false, religion=false)
end
print(", printing"); ec.printRegionalEmission(year, natA3, rgCatFile, region="district", mode=catMode, popwgh=true, ur=false, religion=false)
println(" ... completed")

print(" Exporting: ")
if exportMode || exportWebMode || mapStyleMode;
    print(" GIS-info")
    gisPath = indexFilePath * "gis/"
    regionFile = gisPath * "regions.txt"
    gisCatFile = gisPath * "cat_labels.txt"
    gisConcFile = gisPath * "region_concordance.txt"
    ec.readGISinfo(year, natA3, regionFile, gisCatFile, id = unifiedIdMode)
    ec.buildGISconc(year, natA3, gisConcFile, region = "district", remove = true)

    print(", GIS-exporting")
    # gisTag = "IDN_adm2"
    gisTag = "District"
    exportFile = emissionPath * "YEAR_"*natA3*"_gis_"*subcat*"emission_cat_OvPcTag.csv"
    exportRateFile = emissionPath * "YEAR_"*natA3*"_gis_"*subcat*"emission_cat_dr_OvPcTag.csv"
    labelList, labelListPerCap = ec.exportRegionalEmission(year, natA3, gisTag, exportFile, region="district", mode=expModes,  nspan=128, minmax=minmaxv, descend=true, empty=false, logarithm=false)
    spanVals, spanValsPerCap = ec.exportEmissionDevRate(year, natA3, gisTag, exportRateFile, mode=expModes, maxr=0.5, minr=-0.5, nspan=128, descend=true, empty=false)
end
if exportWebMode; print(", web-files")
    exportPath = emissionPath*"webfile/"
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
println(" ... completed")

println("[all complete]")
