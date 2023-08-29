# SPDX-FileCopyrightText: © 2023 Jemyung Lee <jemyung81@gmail.com>
# SPDX-License-Identifier: GPL-3.0

# Developed date: 11. Apr. 2023
# Last modified date: 29. Aug. 2023
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
include("../GIS/QgisStyleExporter.jl")
using .MicroDataReader
using .ConcMatBuilder
using .EmissionEstimator
using .EmissionCategorizer
using .EmissionCI
using .MapGenerator
using .QgisStyleExporter
mdr = MicroDataReader
cmb = ConcMatBuilder
ee = EmissionEstimator
ec = EmissionCategorizer
ci = EmissionCI
mg = MapGenerator
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

labelConvMode = true    # convert GeoJSON map's label from GIS_ID to GIS_label
groupMode = false       # seperate households by survey group
gr_all_label = "Mixed"  # label that indicate "include all groups"
groupSplit = false      # split household data, which is included in multiple groups, tagging group for HHID (ex."ID0001"->"SURVEY_ID00001")

curConv = true
pppConv = false

skipNullHhs = true      # [true] exclude household that does not have district code
skipNullReg = true      # [true] exclude regions that does not have population

IE_mode = true          # indirect carbon emission estimation
DE_mode = true          # direct carbon emission estimation
DE_factor_estimate = true   # [true] estimate DE factors from IEA datasets, [false] read DE factors
IE_elap = 0             # [n] display IE calculation remained time every 'n' iteration, [0] not display

print_hhs_ie = false    # print household level indirect emission by sector
print_hhs_de = false    # print household level direct emission by sector
print_hhs_cf = false    # print household level CF by sector
print_reg_cf = false    # print region level CF by sector

Qtable = "I_CHG_CO2"

scaleMode = false

gisLabMode = true       # [true] use "GIS_name" ([false] use "City_name") in "GIS_RegionConc" for map city labeling
minSamples = 5          # minimum number of sample houses (include the value, >=)
filterMode = true       # exclude regions that have fewere samples than 'minSamples'
emptyRegRemove = false  # remove empty region from the map

nationDict = Dict{String, String}()
currDict = Dict{String, String}()
boundary_dict = Dict{String, Any}()

# nationDict = Dict("IND" =>"India", "IDN" => "Indonesia", "VNM" => "Viet Nam", "JPN" => "Japan", "USA" => "United States")
# currDict = Dict("IDN" => "IDR", "IND" => "INR", "VNM" => "VND", "JPN" => "JPY", "USA" => "USD")
# boundary_dict = Dict("IND" => [[[0,20000000]], []], "IDN" =>[[[0, 6000000]], []], "VNM" => [[[0,3000000]], []],
#                     "JPN" => [[[0,7000000]], []], "USA" => [[[0, 600000000]], []])

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

# clearconsole()

println("[Pre-process]")

if length(ARGS) > 0
    print(" set conditions")
    condition_file = Base.source_dir() * "/init/" * ARGS[1]
    if condition_file[end-3] != '.'; condition_file *= ".txt"
    elseif !endswith(condition_file, "txt"); println("Note: please enter TXT file name.")
    end
    cnds = Dict{String, String}()
    cfmax, cfmin = 0, 0

    f = open(condition_file)
    for l in eachline(f)
        s = string.(split(l, '\t'))
        cnds[string(filter(x-> !isspace(x), lowercase(s[1])))] = s[2]
    end
    close(f)

    lk = "year"; if haskey(cnds, lk) && tryparse(Int, cnds[lk]) > 0; global year = parse(Int, cnds[lk]) end
    lk = "exchangeyear"; if haskey(cnds, lk) && tryparse(Int, cnds[lk]) > 0; global exchYear = parse(Int, cnds[lk]) end
    lk = "eorayear"; if haskey(cnds, lk) && tryparse(Int, cnds[lk]) > 0; global eoraYear = parse(Int, cnds[lk]) end
    lk = "nation"; if haskey(cnds, lk); global nation = cnds[lk] end
    lk = "nata3"; if haskey(cnds, lk); global natA3 = cnds[lk] end
    lk = "localcurrency"; if haskey(cnds, lk); global natCurr = cnds[lk] end
    lk = "globalcurrency"; if haskey(cnds, lk); global curr_unit = cnds[lk] end
    lk = "emissionunit"; if haskey(cnds, lk); global emiss_unit = cnds[lk] end
    lk = "keydistrict"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global keyDistrict = parse(Bool, cnds[lk]) end
    lk = "keymerging"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global keyMerging = parse(Bool, cnds[lk]) end
    lk = "fiteorayear"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global fitEoraYear = parse(Bool, cnds[lk]) end
    lk = "readmembers"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global readMembers = parse(Bool, cnds[lk]) end
    lk = "conc_float_mode"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global Conc_float_mode = parse(Bool, cnds[lk]) end
    lk = "quantmode"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global quantMode = parse(Bool, cnds[lk]) end
    lk = "labelconvmode"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global labelConvMode = parse(Bool, cnds[lk]) end
    lk = "groupmode"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global groupMode = parse(Bool, cnds[lk]) end
    lk = "groupsplit"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global groupSplit = parse(Bool, cnds[lk]) end
    lk = "curconv"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global curConv = parse(Bool, cnds[lk]) end
    lk = "pppconv"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global pppConv = parse(Bool, cnds[lk]) end
    lk = "skipnullhhs"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global skipNullHhs = parse(Bool, cnds[lk]) end
    lk = "skipnullreg"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global skipNullReg = parse(Bool, cnds[lk]) end
    lk = "ie_mode"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global IE_mode = parse(Bool, cnds[lk]) end
    lk = "de_mode"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global DE_mode = parse(Bool, cnds[lk]) end
    lk = "de_factor_estimate"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global DE_factor_estimate = parse(Bool, cnds[lk]) end
    lk = "ie_elap"; if haskey(cnds, lk) && tryparse(Int, cnds[lk]) > 0; global IE_elap = parse(Int, cnds[lk]) end
    lk = "qtable"; if haskey(cnds, lk); global Qtable = cnds[lk] end
    lk = "scalemode"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global scaleMode = parse(Bool, cnds[lk]) end
    lk = "gisLabmode"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global gisLabMode = parse(Bool, cnds[lk]) end
    lk = "minsamples"; if haskey(cnds, lk) && tryparse(Int, cnds[lk]) != nothing; global minSamples = parse(Int, cnds[lk]) end
    lk = "filtermode"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global filterMode = parse(Bool, cnds[lk]) end
    lk = "exportmode"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global exportMode = parse(Bool, cnds[lk]) end
    lk = "mapgenmode"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global mapGenMode = parse(Bool, cnds[lk]) end
    lk = "subcat"; if haskey(cnds, lk); global subcat = cnds[lk] end
    lk = "ci_rste"; if haskey(cnds, lk) && tryparse(Float64, cnds[lk]) != nothing; global ci_rste = parse(Float64, cnds[lk]) end
    lk = "n_iter"; if haskey(cnds, lk) && tryparse(Int, cnds[lk]) != nothing; global n_iter = parse(Int, cnds[lk]) end
    lk = "cfmax"; if haskey(cnds, lk) && tryparse(Float64, cnds[lk]) != nothing; cfmax = parse(Float64, cnds[lk]) end
    lk = "cfmin"; if haskey(cnds, lk) && tryparse(Float64, cnds[lk]) != nothing; cfmin = parse(Float64, cnds[lk]) end
    lk = "rmemptymap"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global emptyRegRemove = parse(Bool, cnds[lk]) end

    lk = "printhousholdiebysector"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global print_hhs_ie = parse(Bool, cnds[lk]) end
    lk = "printhousholddebysector"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global print_hhs_de = parse(Bool, cnds[lk]) end
    lk = "printhousholdcfbysector"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global print_hhs_cf = parse(Bool, cnds[lk]) end
    lk = "printregioncfbysector"; if haskey(cnds, lk) && tryparse(Bool, cnds[lk]) != nothing; global print_reg_cf = parse(Bool, cnds[lk]) end
else println("\nERROR: Init file is not defined")
end

if length(ARGS) > 1 && tryparse(Int, ARGS[2]) != nothing; global year = parse(Int, ARGS[2])
elseif year == 0; println("\nERROR: Year value is not defined.")
end

global nationDict[natA3] = nation
global currDict[natA3] = natCurr
if cfmax > 0; global boundary_dict[natA3] = [[[cfmin, cfmax]], []] else boundary_dict[natA3] = [[], []] end
println(" ... complete")

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
hhs_cf_file = emissionPath * string(year) * "_" * natA3 * "_hhs_" * scaleTag * "CF.txt"
reg_cf_file = emissionPath * string(year) * "_" * natA3 * "_region_" * scaleTag * "CF.txt"

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
print(" "); mdr.readExtractedRegionData(year, natA3, regInfoFile, key_district = keyDistrict, merged_key = keyMerging, legacy_mode = true, ignore = false, remove_empty = skipNullReg)
print(", "); hh_ngr = mdr.readExtractedHouseholdData(year, natA3, hhsfile, merged_key = keyMerging, skip_empty = skipNullHhs, legacy_mode = true, group_split = groupSplit)
print(", "); sc_ngr = mdr.readExtractedSectorData(year, natA3, cmmfile)
if hh_ngr == sc_ngr > 0; groupMode = true
elseif hh_ngr != sc_ngr; print("Warning: group numbers in houehold and sector files are different.")
end
if readMembers; print(", "); mdr.readExtractedMemberData(year, natA3, mmsfile) end
print(", expenditure"); mdr.readExtractedExpenditureMatrix(year, natA3, exmfile, quantity = quantMode, group_split = groupSplit)
if groupMode; print(", group filtering"); mdr.filterGroupExpenditure(year, natA3, all_label = gr_all_label) end
if fitEoraYear && eoraYear != nothing && eoraYear != year
    print(", scaling from $year to $eoraYear")
    exchYear = eoraYear
    cpiSecFile = indexFilePath * "CPI/CPI_" * natA3 * "_sectors.txt"
    statFile = indexFilePath * "CPI/CPI_" * natA3 * "_values.txt"
    linkFile = indexFilePath * "CPI/CPI_" * natA3 * "_link.txt"
    mdr.scalingExpByCPI(year, natA3, cpiSecFile, statFile, linkFile, eoraYear, period="year", region="district", revHH=false, revMat=true)
end
if curConv; print(", currency"); mdr.exchangeExpCurrency(year, exchYear, natA3, natCurr, erfile, target_curr=curr_unit, exp_mat=true) end
if pppConv; print(", ppp converting"); mdr.convertAvgExpToPPP(eoraYear, natA3, pppfile); println("complete") end
print(", reshaping"); mdr.reshapeCommoditySectors(year, natA3, except = exceptCategory)
print(", find lost"); mdr.findLostRegion(year,natA3)
print(", data")
ee.getDomesticData(year, natA3, mdr.hh_list[year][natA3], mdr.sc_list[year][natA3], mdr.expMatrix[year][natA3], (quantMode ? mdr.qntMatrix[year][natA3] : []), cmmUnit = (quantMode ? mdr.exportCommodityUnit(year, natA3) : []))
cmb.getCommodityCodes(mdr.sc_list[year][natA3])
println(" ... complete")

if DE_mode
    de_conc_file = concordancePath * natA3 * "_" * string(year) * "_LinkedSectors_DE.txt"
    deFile = emissionPath * string(year) * "_" * natA3 * "_hhs_"*scaleTag*"DE.txt"
    if quantMode
        de_conc_file = replace(de_conc_file, ".txt" => "_qnt.txt")
        conmatDeFile = replace(conmatDeFile, ".txt" => "_qnt.txt")
        de_conv_file = replace(de_conv_file, ".txt" => "_qnt.txt")
        if natA3 == "VNM"; de_conv_file = replace(de_conv_file, ".txt" => "_VNMrevised.txt") end
    end

    print(" DE:")
    print(" data reading")
    ee.setNationDict(nationDict)
    ee.readDirectEmissionData(year, natA3, deDataPath, output_path = "", output_tag = "", integrate = true, cpi_scaling = false, cpi_base = 0, cpi_vals = [])
    if DE_factor_estimate
        print(", converting rate")
        price_file = filePath * "de/" * "Price_" * natA3 * "_" * string(year) * curr_unit * ".txt"
        if curr_unit != "USD"; ee.exchangeEmCurrency(year, erfile, target = curr_unit, origin = "USD", output = price_file) end
        ee.calculateEmissionRates(year, output = "", currency = curr_unit, quantity = quantMode)
        ee.printEmissionConvRates(year, de_conv_file, emit_unit = emiss_unit, curr_unit = curr_unit, quantity = quantMode, qnt_unit = "kg")
    end
    print(", intensity")
    ee.readEmissionIntensity(year, natA3, de_sec_file, de_conv_file, emit_unit = emiss_unit, curr_unit = curr_unit, quantity = quantMode)

    print(", concordance")
    ee.readDeConcMat(year, natA3, conmatDeFile, norm = true, output = "", energy_wgh = true, float_mode = Conc_float_mode)
    if quantMode; print(", convert_DE")
        ee.calculateQuantityConvRate(year, natA3, de_conc_file, qnt_unit = "kg")
    end
    print(", estimation"); ee.calculateDirectEmission(year, natA3, quantity = quantMode, full = true)
    if print_hhs_de; print(", printing"); ee.printEmissions(year, natA3, deFile, mode = "de") end
    println(" ... complete")
end

if IE_mode
    nation_file = eoraIndexPath * "Eora_nations.txt"
    sector_file = eoraIndexPath * "Eora_sectors.txt"
    ie_conc_file = concordancePath * natA3 * "_" * string(year) * "_LinkedSectors_IE.txt"
    ieFile = emissionPath * string(year) * "_" * natA3 * "_hhs_"*scaleTag*"IE"*q_tag*".txt"
    eora_index = mrioPath * "index/"
    path = mrioPath * string(eoraYear) * "/" * string(eoraYear)

    print(" IE:")
    print(" concordance")
    cmb.readIeSectors(nation_file, sector_file)
    cmb.readExtractedIeConMat(conmatEoraFile, float_mode = Conc_float_mode)
    cmn_ie = cmb.normConMat()   # {a3, conMat}

    print(", MRIO table")
    ee.readIndex(eora_index)
    ee.readIOTables(eoraYear, path*"_eora_t.csv", path*"_eora_v.csv", path*"_eora_y.csv", path*"_eora_q.csv")
    ee.rearrangeMRIOtables(eoraYear, qmode=Qtable)
    print(", Leontief matrix"); ee.calculateLeontief(eoraYear)

    print(", estimation")
    ee.buildWeightedConcMat(year, eoraYear, natA3, con_mat = cmn_ie, output = "")
    ee.calculateIndirectEmission(year, eoraYear, natA3, full = true, elapChk = IE_elap)
    if print_hhs_ie; print(", printing"); ee.printEmissions(year, natA3, ieFile, mode = "ie") end
    println(" ... complete")
end


println("[CF mapping]")

print(" Micro-data:")
if filterMode; print(" filtering"); mdr.filterData(year, natA3, group=groupMode, region="district", quantity=quantMode) end
print(", population weight"); mdr.calculatePopWeight(year, natA3, "", ur_wgh = false, district=true, province=false, hhs_wgh = true, gr_wgh = groupMode)
print(", import data"); ec.importMicroData(mdr)
println(" ... complete")

print(" Emission categorizing:")
rgCatFile = emissionPath * string(year) * "_" * natA3 * "_region_categorized.txt"
em_mode = (IE_mode ? (DE_mode ? ["ie","de"] : ["ie"]) : (DE_mode ? ["de"] : []))
if groupMode; rgCatGrFile = emissionPath * string(year) * "_" * natA3 * "_region_categorized_grouped.txt" end
if IE_mode || DE_mode; print("import"); ec.importEmissionData(year, natA3, ee, mode = em_mode, revise = true) end
if !DE_mode; print(", DE"); ec.readEmissionData(year, natA3, deFile, mode = "de", revise = true) end
if !IE_mode; print(", IE"); ec.readEmissionData(year, natA3, ieFile, mode = "ie", revise = true) end
if groupMode
    if groupSplit; print(", split groups"); ec.splitHouseholdGroup(year, natA3, mode = ["ie","de"], all_gr = gr_all_label) end
    print(", group filtering"); ec.filterGroupEmission(year, natA3, mode = ["ie","de"], all_gr = gr_all_label)
end
print(", CF"); ec.integrateCarbonFootprint()
print(", category"); ec.setCategory(year, natA3, categories = ces_categories, subgroup = "", except = exceptCategory)
for cm in catMode
    hhCatFile = emissionPath * string(year) * "_" * natA3 * "_hhs_"*uppercase(cm)*"_categorized.txt"
    print(", HHs_"*cm); ec.categorizeHouseholdEmission(year, natA3, mode=cm, output=hhCatFile, hhsinfo=true, group = groupMode, all_gr = gr_all_label)
    print(", Reg_"*cm); ec.categorizeRegionalEmission(year, natA3, mode=cm, period="year", popwgh=true, region="district", ur=false, religion=false, group=groupMode)
end
print(", printing"); ec.printRegionalEmission(year, natA3, rgCatFile, region="district", mode=catMode, popwgh=true, ur=false, religion=false)
if print_hhs_cf; print("_hhs"); ec.printHouseholdEmissionsBySector(year, natA3, hhs_cf_file, mode=["cf"], em_mode = ["household"]) end
if print_reg_cf; print("_region"); ec.printRegionalEmissionBySector(year, natA3, reg_cf_file, region = "district", mode=["cf"], em_mode = ["overall", "percap"]) end
if groupMode; ec.printRegionalGroupEmission(year, natA3, rgCatGrFile, region="district", mode=catMode, popwgh=true, ur=false, gr=groupMode, religion=false) end
println(" ... complete")

print(" Exporting:")
if exportMode || mapGenMode;
    print(" GIS-info")
    ec.readGISinfo(year, natA3, gisRegFile, gisCatFile, id = true)
    ec.buildGISconc(year, natA3, gisConcFile, region = "district", remove = true, merged_key = keyMerging, gis_label_mode = gisLabMode)
    ec.filterRegion(year, natA3; region = "district", limit = minSamples, group = groupMode)

    print(", GIS-exporting")
    gisTag = "District"
    exportFile = emissionPath * "YEAR_"*natA3*"_gis_"*subcat*"emission_cat_OvPcTag.csv"
    exportRateFile = emissionPath * "YEAR_"*natA3*"_gis_"*subcat*"emission_cat_dr_OvPcTag.csv"
    labelList, labelListPerCap = ec.exportRegionalEmission(year, natA3, gisTag, exportFile, region="district", mode=expModes,  nspan=128, minmax=minmaxv, descend=true, empty=false, logarithm=false)
    spanVals, spanValsPerCap = ec.exportEmissionDevRate(year, natA3, gisTag, exportRateFile, mode=expModes, maxr=0.5, minr=-0.5, nspan=128, descend=true, empty=false)
end
if mapGenMode; print(", map-generation")
    mg.importEmissionData(ec, emission = "cf", pc_dev = true, ov_dev = false)
    mg.readBaseMap(year, natA3, basemapFile, remove_reg = emptyRegRemove, alter = true, label_conv = labelConvMode)
    mg.readMapInfo(mapListFile)
    mg.convertRgbToHex(mg.readColorMap(rgbfile_ov, reverse=false), mode = "overall")
    mg.convertRgbToHex(mg.readColorMap(rgbfile_pc, reverse=false), mode = "percap")
    mg.mapRegionCF(year, natA3, label_conv = labelConvMode, blank_color = "#A9A9A9", value_mode = true)
    mg.printMapFiles(year, natA3, mapFilePath)
end
println(" ... complete")


println("[City CF exporting]")

print(" Bootstrap:")
print(" data import"); ci.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, cat_filter = true)
print(", CI calculation"); ci.estimateConfidenceIntervals(year, natA3, iter = n_iter, ci_rate = ci_rste, resample_size = 0, replacement = true, boundary="district", group = groupMode)
println(" ... complete")

print(" Web-file exporting:")
print(", center"); ec.exportCentersFile(year, natA3, web_center_path)
print(", web index"); ci.readCityFileSector(webIndexFile)
print(", city"); ci.exportWebsiteCityFiles(year, natA3, web_city_path, web_categories, cfav_file, cfac_file, boundary="district")
println(" ... complete")

println("[", year, " ", nation, "(",natA3, ")", " all complete]")
