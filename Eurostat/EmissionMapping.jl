# Developed date: 5. Aug. 2020
# Last modified date: 21. Oct. 2021
# Subject: Categorized emission mapping
# Description: Mapping emission through households emissions data, categorizing by district, income-level, and etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("EmissionCategorizer.jl")
using .EmissionCategorizer
ec = EmissionCategorizer
include("../GIS/QgisStyleExporter.jl")
using .QgisStyleExporter
qse = QgisStyleExporter

println("[Process]")

nation = "Eurostat"
year = 2015
years = [2015]
nutsLv = 1
onlyNutsInHbs = true
# Qtable = "_I_CHG_CO2"
Qtable = "_PRIMAP"
ceIntegrateMode = "cf"      # "ie" (only indirect CE), "de" (only direct CE), or "cf" (integrage direct and indirect CEs)
ceProcessMode = ["ie", "de", "cf"]
cePrintMode = ["ie", "de", "cf"]
ceExportMode = "cf"         # "ie" (only indirect CE), "de" (only direct CE), or "cf" (integrage direct and indirect CEs

cpi_scaling = true; base_year = 2010

substMode = true; if substMode; substTag = "_subst" else substTag = "" end
scaleMode = true; if scaleMode; scaleTag = "Scaled_" else scaleTag = "" end

eqvalMode = false   # [true]: apply square root of household size for equivalance scale
ntWeighMode = true  # [true]: apply NUTS population based weight, [false]:apply HBS weight

exportMode = true
minmaxv = [[[0,1.2*10^9]], []] # {{overall CF min., max.}, {CF per capita min., max.}
expNtMode = "hbs"   # ; expNtMode = "gis"
exportWebMode = true
# buildWebFolder = false
mapStyleMode = true; colormapReversePerCap=false; labeRevPerCap=true; colormapReverse=false; labeRev=false

popweight = true
grid_pop = true
expenditureMode = false

# incomeMode = false; relativeMode=false
# religionMode = false
# incomeByReligionMode = false
# expenditureRangeMode = false
# emissionLevelMode = false
#
# costEstimationMode = false
# costEstimationByThresholdMode = false
# costEstimationByReligionMode = false

filePath = Base.source_dir() * "/data/"
indexPath = filePath * "index/"
extrPath = filePath * "extracted/"
emissPath = filePath * "emission/" * string(year) * "/"
indexFile = indexPath * "Eurostat_Index_ver4.6.xlsx"
hhsfile = extrPath * string(year) * "_Households.csv"

ExpenditureFile = extrPath * scaleTag * "Expenditure_matrix_4th" * substTag * ".csv"

incomePeriod = "daily"  # Period: "annual", "monthly"(default), or "daily"

normTag = ["perCapNorm", "perHhNorm"]
categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
subcat=""
# subcat="Food"
foodCategories=["Grain","Vegetable","Fruit","Dairy","Beef","Pork","Poultry","Other meat","Fish",
                "Alcohol","Other beverage","Confectionery","Restaurant","Other food","Food"]

print(" Data reading: ")
print("category"); ec.readCategoryData(indexFile, year, nutsLv, except=["None"], subCategory=subcat)
if length(subcat)==0; ec.setCategory(categories); else ec.setCategory(foodCategories); subcat *= "_" end
print(", household"); ec.readHouseholdData(hhsfile, period = incomePeriod, remove = onlyNutsInHbs, alter=true)
print(", population"); ec.readPopulation(year, indexFile, nuts_lv = nutsLv)
print(", gridded population"); ec.readPopGridded(year, indexFile, nuts_lv = [nutsLv], adjust = true)
print(", emission")
if !expenditureMode
    IE_files = []; DE_files = []; ie_nations = []; de_nations = []
    ie_file_tag = "_hhs_"*scaleTag*"IE"*Qtable*".txt"
    if cpi_scaling && base_year != year; ie_file_tag = replace(ie_file_tag, ".txt" => "_converted_" * string(base_year) * ".txt") end
    de_file_tag = "_hhs_"*scaleTag*"DE.txt"
    for f in readdir(emissPath)
        if startswith(f, string(year)) && endswith(f, ie_file_tag); push!(IE_files, emissPath*f); push!(ie_nations, f[6:7])
        elseif startswith(f, string(year)) && endswith(f, de_file_tag); push!(DE_files, emissPath*f); push!(de_nations, f[6:7])
        end
    end
    print("_IE"); ec.readEmissionData(year, ie_nations, IE_files, mode = "ie")
    print("_DE"); ec.readEmissionData(year, de_nations, DE_files, mode = "de")
    print("_CF"); ec.integrateCarbonFootprint(year, mode=ceIntegrateMode)
# elseif expenditureMode ec.readExpenditure(year, nations, ExpenditureFile)
end
println(" ... complete")

print(" Weights calculating: ")
ec.calculateNutsPopulationWeight(year = year, pop_dens = grid_pop, adjust = true)
println(" ... complete")

print(" Categorizing:")
if expenditureMode; tag = "_exp" else tag = "" end
print(" category")
# hhsDeFile = Base.source_dir() * "/data/emission/2010_EU_hhs_"*scaleTag*subcat*"de_cat.csv"
NutsEmissionFile = emissPath * string(year) * "_EU_nuts_"*scaleTag*subcat*"emission_cat.csv"
for m in ceProcessMode
    print("_",uppercase(m))
    hhsEmissionFile = emissPath * string(year) * "_EU_hhs_"*scaleTag*subcat*uppercase(m)*"_cat.csv"
    ec.categorizeHouseholdEmission(year, mode=m, output=hhsEmissionFile, hhsinfo=false, nutsLv=1)
    ec.categorizeRegionalEmission(year, mode=m, nutsLv=1, period=incomePeriod, adjust=true, religion=false, popWgh=popweight, ntweigh=ntWeighMode)
end
ec.printRegionalEmission(year, NutsEmissionFile, mode=cePrintMode, totm=true, expm=true, popm=true, relm=false, wghm=true, povm=false, ntweigh=ntWeighMode)

print(" National abstract: ")
ec.makeNationalSummary(year, emissPath * string(year) * "_National_summary_"*scaleTag*Qtable*".txt", nuts_mode=true)
println(" ... complete")

# if DEmode; print("_DE");ec.categorizeDirectEmission(years; output=hhsDeFile, hhsinfo=false, nutsLv=1) end
# ec.calculateDistrictPoverty(year, povline=1.9, popWgh=popweight)

if exportMode || exportWebMode || mapStyleMode; print(", GIS-exporting")
    gisTag = "NUTS"
    exportFile = emissPath * "YEAR_EU_NUTS_gis_"*scaleTag*subcat*"emission_cat_OvPcTag.csv"
    exportRateFile = emissPath * "YEAR_EU_NUTS_gis_"*scaleTag*subcat*"emission_cat_dr_OvPcTag.csv"
    labelList, labelListPerCap = ec.exportRegionalEmission(years, gisTag, exportFile, mode=ceExportMode, nutsmode=expNtMode, nspan=128, minmax=minmaxv, descend=true, empty=false, logarithm=false)
    spanVals, spanValsPerCap = ec.exportEmissionDiffRate(years, gisTag, exportRateFile, 0.5, -0.5, 128, descend=true, empty=false)
end
if exportWebMode; print(", web-files")
    exportPath = Base.source_dir() * "/data/emission/"*scaleTag*"webfile/"
    ec.exportWebsiteFiles(years, exportPath, nutsmode=expNtMode, rank=true, empty=false, major=true)
end
# if buildWebFolder; print(", web-folders")
#     centerpath = Base.source_dir() * "/data/index/gis/"
#     outputpath = Base.source_dir() * "/data/emission/"*scaleTag*"webfolder/"
#     ec.buildWebsiteFolder(years, centerpath, outputpath, nutsmode=expNtMode, rank=true)
# end
if mapStyleMode; print(", map-style file generating")
    rgbFile = "../GIS/data/EU/MPL_RdBu.rgb"
    for i=1:length(ec.catList)
        qse.readColorMap(rgbFile, reverse=colormapReversePerCap)
        qmlFile = replace(rgbFile, ".rgb"=>"_percap_"*scaleTag*ec.catList[i]*".qml")
        attr = string(years[1])*"_EU_NUTS_gis_"*subcat*"emission_cat_dr_percap_gr_"*ec.catList[i]
        qse.makeQML(qmlFile, attr, empty=false, values=spanValsPerCap[years[1]][:,i], indexValue=true, labelReverse=labeRevPerCap)
    end

    rgbFile = "../GIS/data/EU/MPL_YlGnBu.rgb"
    for i=1:length(ec.catList)
        qse.readColorMap(rgbFile, reverse=colormapReverse)
        qmlFile = replace(rgbFile, ".rgb"=>"_overall_"*scaleTag*ec.catList[i]*".qml")
        attr = string(years[1])*"_EU_NUTS_gis_"*subcat*"emission_cat_overall_gr_"*ec.catList[i]
        qse.makeQML(qmlFile, attr, empty=false, labels=labelList[years[1]][:,i], indexValue=true, labelReverse=labeRev)
    end
end

# if incomeMode
#     print(" income")
#     incomeFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_inc_"*tag*".csv"
#     intervals = [0.2,0.4,0.6,0.8,1.0]; absint=false; descendig=false
#     # intervals = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]; absint=false; descendig=false
#     # intervals = [150,30]; absint = true
#     eData = ec.categorizeHouseholdByIncome(year,intervals,normMode,popWgh=popweight,sqrRt=eqvalMode,absIntv=absint,perCap=percapita,desOrd=descendig,wghmode=popwghmode)
#     ec.printEmissionByIncome(year, incomeFile, intervals, eData[4], eData[5], eData[6], absIntv=absint, relative=relativeMode)
# end
# if religionMode
#     print(" religion")
#     religionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_rel_"*tag*".csv"
#     eData = ec.categorizeHouseholdByReligion(year, normMode, sqrRt=eqvalMode, popWgh=popweight, wghmode=popwghmode)
#     ec.printEmissionByReligion(year, religionFile, eData[4], eData[5], eData[6])
# end
# if incomeByReligionMode
#     print(" income-religion")
#     intervals = [0.2,0.4,0.6,0.8,1.0]; absint=false; descendig=false
#     incomeReligionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_incByRel_"*tag*".csv"
#     eData = ec.categorizeHouseholdByIncomeByReligion(year,intervals,normMode,popWgh=popweight,sqrRt=eqvalMode,absIntv=absint,perCap=percapita,desOrd=descendig,wghmode=popwghmode)
#     ec.printEmissionByIncomeByReligion(year,incomeReligionFile,intervals,eData[4],eData[5],eData[6],absIntv=absint,desOrd=descendig)
# end
# if expenditureRangeMode
#     print(" expenditure-range")
#     ranges = [1.25, 1.9, 3.0, 4.0, 5.0, 10.0]; absint=true; absSpan=true
#     expRngFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_rng_"*tag*".csv"
#     eData = ec.categorizeHouseholdByExpRange(year,ranges,normMode,over=0.1,less=0.1,absRng=absint,absSpn=absSpan,perCap=percapita,popWgh=popweight,wghmode=popwghmode)
#     ec.printEmissionByRange(year,expRngFile,eData[6], eData[3], eData[4], eData[5], eData[7])
# end
# if emissionLevelMode
#     print(" district")
#     intervals = [0.1, 0.8, 0.1]; absint = false
#     districlevelFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_dis_"*tag*".csv"
#     ec.categorizeDistrictByEmissionLevel(year, normMode, intervals)
#     ec.categorizeHouseholdByEmissionLevel(year, intervals, normMode, squareRoot=eqvalMode, absintv=absint)
#     ec.printEmissionByDistEmLev(year, districlevelFile, intervals)
# end
# println(" ... complete")
#
#
# if costEstimationMode||costEstimationByThresholdMode||costEstimationByReligionMode; print(" Cost estimating:") end
# if costEstimationMode
#     print(" district")
#     intervals = [0.2,0.4,0.6,0.8,1.0]; absint=false; descendig=false; nameMode=true
#     costDistrictFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cost_dis.csv"
#     expDistrictFile = [Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cost_dis_gis.csv", "Censuscode"]
#     ec.estimateEmissionCostByDistrict(year, intervals, normMode, perCap=percapita, popWgh=popweight, desOrd=descendig, name=nameMode, output=costDistrictFile, exportFile=expDistrictFile, wghmode=popwghmode)
# end
# if costEstimationByThresholdMode
#     print(" threshold")
#     thresholds = [1.9, 3.0, 5.0]; overRange=0.1;lessRange=0.1; absint=false; absSpan=true; descendig=false; nameMode=true; stackedMode=false
#     costThresholdFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cost_threshold.csv"
#     expThresholdFile = [Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cost_threshold_gis.csv", "Censuscode"]
#     ec.estimateEmissionCostByDistrictForThreshold(year,thresholds,normMode,perCap=percapita,stacked=stackedMode,absSpn=absSpan,over=overRange,less=lessRange, popWgh=popweight, desOrd=descendig, name=nameMode, output=costThresholdFile, exportFile=expThresholdFile, wghmode=popwghmode)
# end
# if costEstimationByReligionMode
#     print(" religion")
#     intervals = [0.2,0.4,0.6,0.8,1.0]; absint=false; descendig=false; nameMode=true
#     costReligionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cost_rel.csv"
#     expReligionFile = [Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cost_rel_gis.csv", "Censuscode"]
#     ec.estimateEmissionCostByDistrictByReligion(year, intervals, normMode, absIntv=absint, perCap=percapita, popWgh=popweight, desOrd=descendig, name=nameMode, output=costReligionFile, exportFile=expReligionFile, wghmode=popwghmode)
# end
# if costEstimationMode||costEstimationByThresholdMode||costEstimationByReligionMode; println(" ... complete") end

println("\n[Done]")
