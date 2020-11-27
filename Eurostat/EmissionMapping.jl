# Developed date: 5. Aug. 2020
# Last modified date: 27. Nov. 2020
# Subject: Categorized emission mapping
# Description: Mapping emission through households emissions data, categorizing by district, income-level, and etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("EmissionCategorizer.jl")
using .EmissionCategorizer
ec = EmissionCategorizer

println("[Process]")

nation = "Eurostat"
year = 2010
years = [2010]
nutsLv = 1
# Qtable = "_I_CHG_CO2"
Qtable = "_PRIMAP"
substMode = true
scaleMode = true
if substMode; substTag = "_subst" else substTag = "" end
if scaleMode; scaleTag = "Scaled" else scaleTag = "" end

EmissionFilePath = Base.source_dir() * "/data/emission/"
ExpenditureFilePath = Base.source_dir()*"/data/extracted/"*scaleTag*"Expenditure_matrix_4th"*substTag*".csv"
householdFile = Base.source_dir() * "/data/extracted/Households.csv"
indexFile = Base.source_dir() *"/data/index/Eurostat_Index_ver1.9.xlsx"

weightMode = 1  # [0]non-weight, [1]per capita, [2]per household
normMode = 1    # [0]non-weight, [1]per capita, [2]per household
eqvalMode = false   # [true]apply square root of household size for equivalance scale
ntWeighMode = true  # [true]:apply NUTS population based weight, [false]:apply HBS weight

exportMode = false; if weightMode==1; minmaxv = [[0,20000000]] elseif weightMode==4; minmaxv = [] end
exportWebMode = false
mapStyleMode = false; colormapReverse = false; labeRev = false

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

incomePeriod = "daily"

weightTag = ["perCap", "perHH"]
normTag = ["perCap", "perHH"]
categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
subcat=""
# subcat="Food"
foodCategories=["Grain","Vegetable","Fruit","Dairy","Beef","Pork","Poultry","Other meat","Fish",
                "Alcohol","Other beverage","Confectionery","Restaurant","Other food","Food"]

print(" Data reading: ")
print("category"); ec.readCategoryData(indexFile, years, nutsLv, except=["None"], subCategory=subcat)
if length(subcat)==0; ec.setCategory(categories); else ec.setCategory(foodCategories); subcat *= "_" end
print(", household"); ec.readHouseholdData(householdFile, period=incomePeriod)
print(", emission")
if !expenditureMode
    CF_files = []; CE_files = []; nations = []
    for f in readdir(EmissionFilePath)
        if endswith(f, "_hhs_emission"*Qtable*".txt"); push!(CF_files, EmissionFilePath*f); push!(nations, f[6:7])
        elseif endswith(f, "_hhs_CE.txt"); push!(CE_files, EmissionFilePath*f)
        end
    end
    ec.readCarbonFootprint(year, nations, CF_files)
    ec.readDirectEmission(year, nations, CE_files)
# elseif expenditureMode ec.readExpenditure(year, nations, ExpenditureFilePath)
end
println(" ... complete")

print(" National abstract: ")
# ec.makeNationalSummary(year, EmissionFilePath * "National_summary"*Qtable*".txt")
println(" ... complete")

print(" Weights calculating: ")
ec.calculateNutsPopulationWeight()
println(" ... complete")

print(" Categorizing:")
print(" category")
if weightMode==0; tag="non" else tag = weightTag[weightMode] end
if expenditureMode; tag *= "_exp" end
hhsEmissionFile = Base.source_dir() * "/data/emission/2011_EU_hhs_"*subcat*"emission_cat.csv"
NutsEmissionFile = Base.source_dir() * "/data/emission/2011_EU_nuts_"*subcat*"emission_cat_"*tag*".csv"
ec.categorizeHouseholdEmission(years, output=hhsEmissionFile, hhsinfo=false, nutsLv=1)
# ec.calculateDistrictPoverty(year, povline=1.9, popWgh=popweight)
ec.categorizeRegionalEmission(years, weightMode, nutsLv=1, period="daily", religion=false, popWgh=popweight, ntweigh=ntWeighMode)
    # Period for MPCE: "annual", "monthly"(default), or "daily"
if weightMode == 1; overallMode = true else overallMode = false end
ec.printRegionalEmission(years, NutsEmissionFile, totm=overallMode, expm=true, popm=true, relm=false, wghm=true, povm=false, ntweigh=ntWeighMode)
println()

#
# if exportMode || exportWebMode || mapStyleMode
#     print(", exporting")
#     #gidTag = "GID_2"
#     gidTag = "Censuscode"
#     exportFile = Base.source_dir() * "/data/emission/2011_IND_dist_GIS_"*subcat*"emission_cat_"*tag*".csv"
#     exportRateFile = Base.source_dir() * "/data/emission/2011_IND_dist_GIS_"*subcat*"emission_cat_dr_"*tag*".csv"
#     exportRankFile = Base.source_dir() * "/data/emission/2011_IND_dist_GIS_"*subcat*"emission_cat_rnk_"*tag*".csv"
#     exportOrderFile = Base.source_dir() * "/data/emission/2011_IND_dist_GIS_"*subcat*"emission_cat_ord_"*tag*"_gr.csv"
#     gData = ec.exportDistrictEmission(year, gidTag, exportFile, weightMode, logarithm=false, descend=true, empty=true, minmax=minmaxv)
#     drData = ec.exportEmissionDiffRate(year, gidTag, exportRateFile, 0.5, -0.5, 128, descend=true, empty=true)
# #    ec.exportEmissionValGroup(year, gidTag, exportRankFile, 128, descend=true, logscl=false)
# #    ec.exportEmissionRankGroup(year, gidTag, exportOrderFile, 128, descend=true)
# end
# if exportWebMode
#     print(", web-file exporting")
#     exportPath = Base.source_dir() * "/data/emission/webfile/"
#     ec.exportWebsiteFiles(year,exportPath,weightMode,gData[2],gData[5],gData[6],gData[3],gData[7],rank=true,empty=true)
# end
# if mapStyleMode
#     print(", map-style file generating")
#     if weightMode==1||weightMode==2; rgbFile = "../GIS/data/MPL_YlGnBu.rgb"
#     elseif weightMode==4||weightMode==5; rgbFile = "../GIS/data/MPL_RdBu.rgb"
#     end
#     for i=1:length(ec.catList)
#         qse.readColorMap(rgbFile, reverse=colormapReverse)
#         qmlFile = replace(rgbFile, ".rgb"=>"_"*tag*"_"*ec.catList[i]*".qml")
#         if weightMode==1||weightMode==2
#             attr = "2011_IND_dist_GIS_emission_cat_"*tag*"_gr_"*ec.catList[i]
#             qse.makeQML(qmlFile, attr, empty=true, labels=gData[9][:,i], indexValue=true, labelReverse=labeRev)
#         elseif weightMode==4||weightMode==5
#             attr = "2011_IND_dist_GIS_emission_cat_dr_"*tag*"_gr_"*ec.catList[i]
#             qse.makeQML(qmlFile, attr, empty=true, values=drData[4][:,i], indexValue=true, labelReverse=labeRev)
#         end
#     end
# end
#
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


println("[Done]")
