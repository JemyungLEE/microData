# Developed date: 5. Aug. 2020
# Last modified date: 9. Feb. 2021
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
year = 2010
years = [2010]
nutsLv = 1
# Qtable = "_I_CHG_CO2"
Qtable = "_PRIMAP"
CFmode = true
DEmode = true
substMode = true
scaleMode = true
if substMode; substTag = "_subst" else substTag = "" end
if scaleMode; scaleTag = "Scaled" else scaleTag = "" end

EmissionFilePath = Base.source_dir() * "/data/emission/"
ExpenditureFilePath = Base.source_dir()*"/data/extracted/"*scaleTag*"Expenditure_matrix_4th"*substTag*".csv"
householdFile = Base.source_dir() * "/data/extracted/Households.csv"
indexFile = Base.source_dir() *"/data/index/Eurostat_Index_ver3.1.xlsx"

eqvalMode = false   # [true]: apply square root of household size for equivalance scale
ntWeighMode = true  # [true]: apply NUTS population based weight, [false]:apply HBS weight

exportMode = true; minmaxv = [[[0,2.0*10^8]], []] # {{overall CF min., max.}, {CF per capita min., max.}
expNtMode = "hbs"   # ; expNtMode = "gis"
exportWebMode = false
buildWebFolder = false
mapStyleMode = true; colormapReversePerCap=false; labeRevPerCap=true; colormapReverse=false; labeRev=false

popweight = true
expenditureMode = false

incomeMode = false; relativeMode=false
religionMode = false
incomeByReligionMode = false
expenditureRangeMode = false
emissionLevelMode = false

costEstimationMode = false
costEstimationByThresholdMode = false
costEstimationByReligionMode = false

incomePeriod = "daily"  # Period: "annual", "monthly"(default), or "daily"

normTag = ["perCapNorm", "perHhNorm"]
categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
subcat=""
# subcat="Food"
foodCategories=["Grain","Vegetable","Fruit","Dairy","Beef","Pork","Poultry","Other meat","Fish",
                "Alcohol","Other beverage","Confectionery","Restaurant","Other food","Food"]

print(" Data reading: ")
print("category"); ec.readCategoryData(indexFile, years, nutsLv, except=["None"], subCategory=subcat, nuts3pop=true)
if length(subcat)==0; ec.setCategory(categories); else ec.setCategory(foodCategories); subcat *= "_" end
print(", household"); ec.readHouseholdData(householdFile, period=incomePeriod)
print(", emission")
if !expenditureMode
    CF_files = []; DE_files = []; nations = []
    for f in readdir(EmissionFilePath)
        if endswith(f, "_hhs_emission"*Qtable*".txt"); push!(CF_files, EmissionFilePath*f); push!(nations, f[6:7])
        elseif endswith(f, "_hhs_DE.txt"); push!(DE_files, EmissionFilePath*f)
        end
    end
    print("_CF"); ec.readCarbonFootprint(year, nations, CF_files)
    print("_DE"); ec.readDirectEmission(year, nations, DE_files)
# elseif expenditureMode ec.readExpenditure(year, nations, ExpenditureFilePath)
end
println(" ... complete")

# print(" National abstract: ")
# ec.makeNationalSummary(year, EmissionFilePath * "National_summary"*Qtable*".txt")
# println(" ... complete")

print(" Weights calculating: ")
ec.calculateNutsPopulationWeight()
println(" ... complete")

print(" Categorizing:")
if expenditureMode; tag = "_exp" else tag = "" end
print(" category")
hhsEmissionFile = Base.source_dir() * "/data/emission/2010_EU_hhs_"*subcat*"emission_cat.csv"
hhsDeFile = Base.source_dir() * "/data/emission/2010_EU_hhs_"*subcat*"de_cat.csv"
NutsEmissionFile = Base.source_dir() * "/data/emission/2010_EU_nuts_"*subcat*"emission_cat.csv"
if CFmode; print("_CF"); ec.categorizeHouseholdEmission(years, output=hhsEmissionFile, hhsinfo=false, nutsLv=1) end
if DEmode; print("_DE");ec.categorizeDirectEmission(years; output=hhsDeFile, hhsinfo=false, nutsLv=1) end
# ec.calculateDistrictPoverty(year, povline=1.9, popWgh=popweight)
ec.categorizeRegionalEmission(years, cf=CFmode, de=DEmode, nutsLv=1, period=incomePeriod, religion=false, popWgh=popweight, ntweigh=ntWeighMode)
ec.printRegionalEmission(years, cf=CFmode, de=DEmode, NutsEmissionFile, totm=true, expm=true, popm=true, relm=false, wghm=true, povm=false, ntweigh=ntWeighMode)

if exportMode || exportWebMode || mapStyleMode; print(", GIS-exporting")
    gisTag = "NUTS"
    exportFile = Base.source_dir() * "/data/emission/YEAR_EU_NUTS_gis_"*subcat*"emission_cat_OvPcTag.csv"
    exportRateFile = Base.source_dir() * "/data/emission/YEAR_EU_NUTS_gis_"*subcat*"emission_cat_dr_OvPcTag.csv"
    label_list = ec.exportRegionalEmission(years, gisTag, exportFile, cf=CFmode, de=DEmode, nutsmode=expNtMode, nspan=128, minmax=minmaxv, descend=true, empty=false, logarithm=false)
    # labelList, labelListPerCap = ec.exportRegionalEmission(years, gisTag, exportFile, cf=CFmode, de=DEmode, nutsmode=expNtMode, nspan=128, minmax=minmaxv, descend=true, empty=false, logarithm=false)
    spanVals, spanValsPerCap = ec.exportEmissionDiffRate(years, gisTag, exportRateFile, 0.5, -0.5, 128, descend=true, empty=false)
    if CFmode; labelList, labelListPerCap = label_list[1], label_list[2]
        if DEmode; labelListDe, labelListPerCapDe = label_list[3], label_list[4] end
    elseif DEmode; labelListDe, labelListPerCapDe = label_list[1], label_list[2]
    end
end
if exportWebMode; print(", web-files")
    exportPath = Base.source_dir() * "/data/emission/webfile/"
    ec.exportWebsiteFiles(years, exportPath, nutsmode=expNtMode, rank=true, empty=false, major=true)
end
if buildWebFolder; print(", web-folders")
    centerpath = Base.source_dir() * "/data/index/gis/"
    outputpath = Base.source_dir() * "/data/emission/webfolder/"
    ec.buildWebsiteFolder(years, centerpath, outputpath, nutsmode=expNtMode, rank=true)
end
if mapStyleMode; print(", map-style file generating")
    rgbFile = "../GIS/data/EU/MPL_RdBu.rgb"
    for i=1:length(ec.catList)
        qse.readColorMap(rgbFile, reverse=colormapReversePerCap)
        qmlFile = replace(rgbFile, ".rgb"=>"_percap_"*ec.catList[i]*".qml")
        attr = string(years[1])*"_EU_NUTS_gis_"*subcat*"emission_cat_dr_percap_gr_"*ec.catList[i]
        qse.makeQML(qmlFile, attr, empty=false, values=spanValsPerCap[years[1]][:,i], indexValue=true, labelReverse=labeRevPerCap)
    end

    rgbFile = "../GIS/data/EU/MPL_YlGnBu.rgb"
    for i=1:length(ec.catList)
        qse.readColorMap(rgbFile, reverse=colormapReverse)
        qmlFile = replace(rgbFile, ".rgb"=>"_overall_"*ec.catList[i]*".qml")
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
