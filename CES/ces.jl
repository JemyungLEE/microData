# Developed date: 11. Apr. 2023
# Last modified date: 12. Apr. 2023
# Subject: Carbon Estimation System
# Description: Read household consumption data, estimate household carbon footprint,
#              categorize CF into eleven categories, and map regional CFs.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

include("MicroDataReader.jl")
include("ConcMatBuilder.jl")
include("EmissionEstimator.jl")
include("EmissionCategorizer.jl")
include("EmissionCI.jl")
include("MapGenerator.jl")
include("Preprocess.jl")
include("../GIS/QgisStyleExporter.jl")
using .MicroDataReader
using .ConcMatBuilder
using .EmissionEstimator
using .EmissionCategorizer
using .EmissionCI
using .MapGenerator
using .Preprocess
using .QgisStyleExporter
mdr = MicroDataReader
cmb = ConcMatBuilder
ee = EmissionEstimator
ec = EmissionCategorizer
ci = EmissionCI
mg = MapGenerator
pp = Preprocess
qse = QgisStyleExporter

cd(Base.source_dir())

year = 0
exchYear = 0
eoraYear = 0
nation = ""
natA3 = ""
natCurr = ""
curr_unit= "USD"
emiss_unit = "tCO2"
keyDistrict = true
keyMerging = true
fitEoraYear = false     # scaling micro-data's expenditure to fit the Eora target year
readMembers = false     # read member data
Conc_float_mode = false
quantMode = false

labelConvMode = false    # convert GeoJSON map's label from GIS_ID to GIS_label
groupMode = false        # seperate households by survey group

curConv = true
pppConv = false

skipNullHhs = true      # [true] exclude household that does not have district code

IE_mode = true             # indirect carbon emission estimation
DE_mode = true              # direct carbon emission estimation
DE_factor_estimate = true   # [true] estimate DE factors from IEA datasets, [false] read DE factors

Qtable = "I_CHG_CO2"

scaleMode = false

gisLabMode = true       # [true] use "GIS_name" ([false] use "City_name") in "GIS_RegionConc" for map city labeling
minSamples = 5          # minimum number of sample houses (include the value, >=)
filterMode = true      # exclude regions that have fewere samples than 'minSamples'

nationDict = Dict("IND" =>"India", "IDN" => "Indonesia", "VNM" => "Viet Nam", "JPN" => "Japan", "USA" => "United States")
currDict = Dict("IDN" => "IDR", "IND" => "INR", "VNM" => "VND", "JPN" => "JPY", "USA" => "USD")
boundary_dict = Dict("IND" => [[[0,20000000]], []], "IDN" =>[[[0, 6000000]], []], "VNM" => [[[0,3000000]], []],
                    "JPN" => [[[0,7000000]], []], "USA" => [[[0, 600000000]], []])

exportMode = true   # export GIS CF resutls
mapGenMode = true   # generate GeoJSON maps

expModes = ["cf"]
catMode = ["cf"]

ces_categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
web_categories = ["FOOD", "ELECTRICITY", "GAS", "ENERGY", "PUBLIC_TRANS", "PRIVATE_TRANS", "MEDICAL",
                "EDUCATION", "CONSUMABLE", "DURABLE", "SERVICES", "ALL"]

exceptCategory = ["None", "Taxes"]

subcat=""

ci_rste = 0.95      # confidence interval: 95% (default)
n_iter = 10000      # bootstrap iteration

clearconsole()

println("[Pre-process]")
print("set conditions")
condition_file = ARGS[1]
cnds = Dict{String, String}()
cfmax, cfmin = 0, 0

f = open(condition_file)
for l in eachline(f)
    s = string.(filter.(x-> !isspace(x), split(lowercase(l), '\t')))
    cnds[s[1]] = s[2]
end
close(f)

lk = "year"; if haskey(cnds, lk) && isa(cnds[lk], Int); ces_vars.year = parse(Int, cnds[lk]) end
lk = "exchangeyear"; if haskey(cnds, lk) && isa(cnds[lk], Int); ces_vars.exchYear = parse(Int, cnds[lk]) end
lk = "eorayear"; if haskey(cnds, lk) && isa(cnds[lk], Int); ces_vars.eoraYear = parse(Int, cnds[lk]) end
lk = "nation"; if haskey(cnds, lk); ces_vars.nation = cnds[lk] end
lk = "nata3"; if haskey(cnds, lk); ces_vars.natA3 = cnds[lk] end
lk = "localcurrency"; if haskey(cnds, lk); ces_vars.natCurr = cnds[lk] end
lk = "globalcurrency"; if haskey(cnds, lk); ces_vars.curr_unit = cnds[lk] end
lk = "emissionunit"; if haskey(cnds, lk); ces_vars.emiss_unit = cnds[lk] end
lk = "keydistrict"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.keyDistrict = parse(Bool, cnds[lk]) end
lk = "keymerging"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.keyMerging = parse(Bool, cnds[lk]) end
lk = "fitEoraYear"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.fitEoraYear = parse(Bool, cnds[lk]) end
lk = "readMembers"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.readMembers = parse(Bool, cnds[lk]) end
lk = "Conc_float_mode"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.Conc_float_mode = parse(Bool, cnds[lk]) end
lk = "quantMode"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.quantMode = parse(Bool, cnds[lk]) end
lk = "labelConvMode"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.labelConvMode = parse(Bool, cnds[lk]) end
lk = "groupMode"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.groupMode = parse(Bool, cnds[lk]) end
lk = "curConv"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.curConv = parse(Bool, cnds[lk]) end
lk = "pppConv"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.pppConv = parse(Bool, cnds[lk]) end
lk = "skipNullHhs"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.skipNullHhs = parse(Bool, cnds[lk]) end
lk = "IE_mode"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.IE_mode = parse(Bool, cnds[lk]) end
lk = "DE_mode"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.DE_mode = parse(Bool, cnds[lk]) end
lk = "DE_factor_estimate"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.DE_factor_estimate = parse(Bool, cnds[lk]) end
lk = "Qtable"; if haskey(cnds, lk); ces_vars.Qtable = cnds[lk] end
lk = "scaleMode"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.scaleMode = parse(Bool, cnds[lk]) end
lk = "gisLabMode"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.gisLabMode = parse(Bool, cnds[lk]) end
lk = "minSamples"; if haskey(cnds, lk) && isa(cnds[lk], Int); ces_vars.minSamples = parse(Int, cnds[lk]) end
lk = "filterMode"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.filterMode = parse(Bool, cnds[lk]) end
lk = "exportMode"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.exportMode = parse(Bool, cnds[lk]) end
lk = "mapGenMode"; if haskey(cnds, lk) && isa(cnds[lk], Bool); ces_vars.mapGenMode = parse(Bool, cnds[lk]) end
lk = "subcat"; if haskey(cnds, lk); ces_vars.subcat = cnds[lk] end
lk = "ci_rste"; if haskey(cnds, lk) && isa(cnds[lk], Number); ces_vars.ci_rste = parse(Float64, cnds[lk]) end
lk = "n_iter"; if haskey(cnds, lk) && isa(cnds[lk], Int); ces_vars.n_iter = parse(Int, cnds[lk]) end
lk = "cfmax"; if haskey(cnds, lk) && isa(cnds[lk], Number); cfmax = parse(Float64, cnds[lk]) end
lk = "cfmin"; if haskey(cnds, lk) && isa(cnds[lk], Number); cfmin = parse(Float64, cnds[lk]) end

ces_vars.nationDict[natA3] = nation
ces_vars.currDict[natA3] = natCurr
if cfmax > 0; ces_vars.boundary_dict[natA3] = [[[cfmin, cfmax]], []] end
println(" ... complete")
println("[prepared]")

if exchYear == 0; exchYear = year end
if eoraYear == 0; eoraYear = year end

filePath = Base.source_dir() * "/data/" * natA3 * "/"
mrioPath = "../Eora/data/"
commonIndexPath = Base.source_dir() * "/data/Common/"

indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissionPath = filePath * "emission/" * string(year) * "/"
concordancePath = filePath * "concordance/"
eoraIndexPath = commonIndexPath * "Eora/"
deDataPath = commonIndexPath* "DE/"
gisIndexPath = commonIndexPath * "gis/"
webIndexPath = commonIndexPath * "web/"

if IE_mode && !DE_mode; quantMode = false end   # quantity mode is available for the direct emission mode only
if curConv; currTag = curr_unit else currTag = natCurr end
if scaleMode; scaleTag = "_Scaled" else scaleTag = "" end
pppfile = filePath * "PPP_ConvertingRates.txt"

q_tag = "_" * lowercase(Qtable)

minmaxv = boundary_dict[natA3]  # {{overall CF min., max.}, {CF per capita min., max.}

natFileTag = "source/" * string(year) * "/" * natA3 * "_" * string(year)
regInfoFile = filePath * natFileTag * "_MD_RegionInfo.txt"
cmmfile = filePath * natFileTag * "_MD_Commodities.txt"
hhsfile = filePath * natFileTag * "_MD_Households_"*natCurr*".txt"
mmsfile = filePath * natFileTag * "_MD_Members.txt"
exmfile = filePath * natFileTag * "_MD_ExpenditureMatrix_"*natCurr*".txt"
erfile = filePath * natFileTag * "_MD_ExchangeRate.txt"
if !isfile(hhsfile); hhsfile = filePath * natFileTag * "_MD_Households.txt" end
if !isfile(exmfile); exmfile = filePath * natFileTag * "_MD_Expenditure.txt" end
if !isfile(erfile); erfile = filePath * "source/" * natA3 * "_MD_ExchangeRate.txt" end
if !isfile(erfile); erfile = commonIndexPath * "CurrencyExchangeRates.txt" end

conmatEoraFile = filePath * natFileTag * "_IOT_ConcMatEora.txt"
conmatDeFile = filePath * natFileTag * "_IOT_ConcMatDe.txt"

expfile = filePath * natFileTag * "_MD_Expenditure_"*natCurr*".txt"

de_sec_file = deDataPath * (!quantMode ? "DE_sectors.txt" : "Emission_sectors.txt")
de_conv_file = commonIndexPath * "Emission_converting_rate.txt"

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

web_city_path = filePath * "web/" * "footprint/"
web_center_path = filePath * "web/" * "centers/"

webIndexFile = webIndexPath * "keycode_index.txt"

cfav_file = emissionPath * string(year) * "_" * natA3 * "_gis_" * subcat * "emission_cat_overall_CF_gr.csv"
cfac_file = emissionPath * string(year) * "_" * natA3 * "_gis_" * subcat * "emission_cat_dr_percap_CF_gr.csv"


println("[CF estimation]")

print(" Micro-data reading:")
print(" regions"); mdr.readExtractedRegionData(year, natA3, regInfoFile, key_district = keyDistrict, merged_key = keyMerging, legacy_mode = true, ignore = false)
print(", households"); mdr.readExtractedHouseholdData(year, natA3, hhsfile, merged_key = keyMerging, skip_empty = skipNullHhs, legacy_mode = true)
if readMembers; print(", members"); mdr.readExtractedMemberData(year, natA3, mmsfile) end
print(", sectors"); mdr.readExtractedSectorData(year, natA3, cmmfile)
print(", expenditure matrix"); mdr.readExtractedExpenditureMatrix(year, natA3, exmfile, quantity = quantMode)
if fitEoraYear && eoraYear != nothing && eoraYear != year
    print(", scaling from $year to $eoraYear")
    exchYear = eoraYear
    cpiSecFile = indexFilePath * "CPI/CPI_" * natA3 * "_sectors.txt"
    statFile = indexFilePath * "CPI/CPI_" * natA3 * "_values.txt"
    linkFile = indexFilePath * "CPI/CPI_" * natA3 * "_link.txt"
    mdr.scalingExpByCPI(year, natA3, cpiSecFile, statFile, linkFile, eoraYear, period="year", region="district", revHH=false, revMat=true)
end
if curConv; print(", currency exchange"); mdr.exchangeExpCurrency(year, exchYear, natA3, natCurr, erfile, target_curr=curr_unit, exp_mat=true) end
if pppConv; print(", ppp converting"); mdr.convertAvgExpToPPP(eoraYear, natA3, pppfile); println("complete") end
println(" ... completed")

print(" Concordance matrix building:")
print(" commodity_code"); cmb.getCommodityCodes(mdr.sc_list[year][natA3])
if IE_mode
    nation_file = eoraIndexPath * "Eora_nations.txt"
    sector_file = eoraIndexPath * "Eora_sectors.txt"
    ie_conc_file = concordancePath * natA3 * "_" * string(year) * "_LinkedSectors_IE.txt"

    print(", IE data reading"); cmb.readIeSectors(nation_file, sector_file)
    print(", IE matrix reading"); cmb.readExtractedIeConMat(conmatEoraFile, float_mode = Conc_float_mode)
    print(", normalization"); cmn_ie = cmb.normConMat()   # {a3, conMat}
end
if DE_mode
    print(", DE data reading")
    ee.setNationDict(nationDict)
    ee.readDirectEmissionData(year, natA3, deDataPath, output_path = filePath * "de/", output_tag = natA3, integrate = true, cpi_scaling = false, cpi_base = 0, cpi_vals = [])
    if DE_factor_estimate
        print(", estimation")
        price_file = filePath * "de/" * "Price_" * natA3 * "_" * string(year) * curr_unit * ".txt"
        emi_intens_file = filePath * "de/" * "Emission_intensity_" * natA3 * "_" * string(year) * "_tCO2_per_" * curr_unit* ".txt"
        if curr_unit != "USD"; ee.exchangeEmCurrency(year, erfile, target = curr_unit, origin = "USD", output = price_file) end
        ee.calculateEmissionRates(year, output = emi_intens_file, currency = curr_unit, quantity = quantMode)
        ee.printEmissionConvRates(year, de_conv_file, emit_unit = emiss_unit, curr_unit = curr_unit, quantity = quantMode, qnt_unit = "kg")
    end
    print(", intensity")
    if quantMode; de_conv_file = replace(de_conv_file, ".txt" => "_qnt.txt") end
    if natA3 == "VNM"; de_conv_file = replace(de_conv_file, ".txt" => "_VNMrevised.txt") end
    ee.readEmissionIntensity(year, natA3, de_sec_file, de_conv_file, emit_unit = emiss_unit, curr_unit = curr_unit, quantity = quantMode)
end
println(" ... complete")

if IE_mode
    print(" MRIO table reading:")
    eora_index = mrioPath * "index/"
    path = mrioPath * string(eoraYear) * "/" * string(eoraYear)
    print(" index"); ee.readIndex(eora_index)
    print(", table"); ee.readIOTables(eoraYear, path*"_eora_t.csv", path*"_eora_v.csv", path*"_eora_y.csv", path*"_eora_q.csv")
    print(", rearranging"); ee.rearrangeMRIOtables(eoraYear, qmode=Qtable)
    print(", Leontief matrix"); ee.calculateLeontief(eoraYear)
    println(" ... complete")
end

print(" Emission calculation: ")
print("data"); ee.getDomesticData(year, natA3, mdr.hh_list[year][natA3], mdr.sc_list[year][natA3], mdr.expMatrix[year][natA3], (quantMode ? mdr.qntMatrix[year][natA3] : []), cmmUnit = (quantMode ? mdr.exportCommodityUnit(year, natA3) : []))
if DE_mode
    de_conc_file = concordancePath * natA3 * "_" * string(year) * "_LinkedSectors_DE.txt"
    de_conc_mat_file = concordancePath * natA3 * "_" * string(year) * "_ConcMat_DE.txt"
    deFile = emissionPath * string(year) * "_" * natA3 * "_hhs_"*scaleTag*"DE.txt"
    if quantMode
        de_conc_file = replace(de_conc_file, ".txt" => "_qnt.txt")
        conmatDeFile = replace(conmatDeFile, ".txt" => "_qnt.txt")
    end
    print(", concordance_DE")
    ee.readDeConcMat(year, natA3, conmatDeFile, norm = true, output = de_conc_mat_file, energy_wgh = true, float_mode = Conc_float_mode)
    if quantMode; print(", convert_DE")
        ee.calculateQuantityConvRate(year, natA3, de_conc_file, qnt_unit = "kg")
    end
    print(", estimate_DE"); ee.calculateDirectEmission(year, natA3, quantity = quantMode, full = true)
    print(", print_DE"); ee.printEmissions(year, natA3, deFile, mode = "de")
end
if IE_mode
    ieFile = emissionPath * string(year) * "_" * natA3 * "_hhs_"*scaleTag*"IE"*q_tag*".txt"
    print(", concordance_IE"); ee.buildWeightedConcMat(year, eoraYear, natA3, con_mat = cmn_ie, output = "")
    print(", estimate_IE"); ee.calculateIndirectEmission(year, eoraYear, natA3, full = true, elapChk=1)
    print(", print_IE"); ee.printEmissions(year, natA3, ieFile, mode = "ie")
end
println(" ... complete")

println(" ... Estimation complete")


println("[CF mapping]")

print(" Micro-data:")
if filterMode; print(", filtering"); mdr.filterRegionData(year, natA3) end
print(", find lost"); mdr.findLostRegion(year,natA3)
if readMembers; print(", members"); mdr.readExtractedMemberData(year, natA3, mmsfile) end
print(", population weight"); mdr.calculatePopWeight(year, natA3, "", ur_wgh = false, district=true, province=false, hhs_wgh = true, gr_wgh = groupMode)
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
    print(", HHs_"*cm); ec.categorizeHouseholdEmission(year, natA3, mode=cm, output=hhCatFile, hhsinfo=true, group = groupMode)
    print(", Reg_"*cm); ec.categorizeRegionalEmission(year, natA3, mode=cm, period="year", popwgh=true, region="district", ur=false, religion=false, group=groupMode)
end
print(", printing"); ec.printRegionalEmission(year, natA3, rgCatFile, region="district", mode=catMode, popwgh=true, ur=false, religion=false)
if groupMode; ec.printRegionalGroupEmission(year, natA3, rgCatGrFile, region="district", mode=catMode, popwgh=true, ur=false, gr=groupMode, religion=false) end
println(" ... completed")

print(" Exporting: ")
if exportMode || mapGenMode;
    print(" GIS-info")
    ec.readGISinfo(year, natA3, gisRegFile, gisCatFile, id = unifiedIdMode)
    ec.buildGISconc(year, natA3, gisConcFile, region = "district", remove = true, merged_key = keyMerging, gis_label_mode = gisLabMode)
    ec.filterRegion(year, natA3; region = "district", limit = minSamples)

    print(", GIS-exporting")
    gisTag = "District"
    exportFile = emissionPath * "YEAR_"*natA3*"_gis_"*subcat*"emission_cat_OvPcTag.csv"
    exportRateFile = emissionPath * "YEAR_"*natA3*"_gis_"*subcat*"emission_cat_dr_OvPcTag.csv"
    labelList, labelListPerCap = ec.exportRegionalEmission(year, natA3, gisTag, exportFile, region="district", mode=expModes,  nspan=128, minmax=minmaxv, descend=true, empty=false, logarithm=false)
    spanVals, spanValsPerCap = ec.exportEmissionDevRate(year, natA3, gisTag, exportRateFile, mode=expModes, maxr=0.5, minr=-0.5, nspan=128, descend=true, empty=false)
end
if mapGenMode; print(", map-generation")
    mg.importEmissionData(ec, emission = "cf", pc_dev = true, ov_dev = false)
    mg.readBaseMap(year, natA3, basemapFile, remove_reg = false, alter = true, label_conv = labelConvMode)
    mg.readMapInfo(mapListFile)
    mg.convertRgbToHex(mg.readColorMap(rgbfile_ov, reverse=false), mode = "overall")
    mg.convertRgbToHex(mg.readColorMap(rgbfile_pc, reverse=false), mode = "percap")
    mg.mapRegionCF(year, natA3, label_conv = labelConvMode, blank_color = "#A9A9A9", value_mode = true)
    mg.printMapFiles(year, natA3, mapFilePath)
end
println(" ... completed")

println(" ... Mapping complete")


println("[City CF exporting]")

print(" Bootstrap process:")
print(" data import"); ci.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, cat_filter = true)
print(", CI calculation"); ci.estimateConfidenceIntervals(year, natA3, iter = n_iter, ci_rate = ci_rste, resample_size = 0, replacement = true, boundary="district")
println(" ... completed")

print(" Web-file exporting:")
print("set category"); ec.setCategory(year, natA3, categories = ces_categories, subgroup = "", except = exceptCategory)
print(", center"); ec.exportCentersFile(year, natA3, web_center_path)
print(", web index"); ci.readCityFileSector(webIndexFile)
print(", city"); ci.exportWebsiteCityFiles(year, natA3, web_city_path, web_categories, city_file_sector, cfav_file, cfac_file, boundary="district")
println(" ... completed")

println(" ... Exporting complete")

println("[", year, " ", nation, "(",natA3, ")", " all complete]")
