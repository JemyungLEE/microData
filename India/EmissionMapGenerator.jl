# Developed date: 27. Dec. 2019
# Last modified date: 3. Sep. 2020
# Subject: Emission mapping
# Description: Mapping emission through households emissions data
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("EmissionCategorizer.jl")
using .EmissionCategorizer
ec = EmissionCategorizer
include("SamplingError.jl")
using .SamplingError
se = SamplingError
include("../Plot/EmissionPlots.jl")
using .EmissionPlots
ep = EmissionPlots
include("../GIS/QgisStyleExporter.jl")
using .QgisStyleExporter
qse = QgisStyleExporter

println("[Process]")

nation = "IND"
year = 2011
emissionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission.txt"
householdFile = Base.source_dir() * "/data/extracted/Households.txt"
sectorFile = Base.source_dir() *"/data/index/IND_index_match_v1.3.xlsx"

mergingMode = true # true: proceed district merging, default=false

weightMode = 4  # [0]non-weight, [1]population weighted, [2]household weighted, [3]both population and household weighted
                # ([4],[5]: normalization) [4]per capita, [5]per household
                # (basic information) [6]population and households, [1,:]population, [2,:]households
normMode = 1    # [0]non-weight, [1]per capita, [2]per houehold,
                # (basic information) [3]population and households by religions, [1,:]population, [2,:]households
eqvalMode = false   # [true]apply square root of household size for equivalance scale

exportMode = false; if weightMode==1; minmaxv = [[0,20000000]] elseif weightMode==4; minmaxv = [] end
exportWebMode = false
mapStyleMode = false; colormapReverse = false; labeRev = true

percapita = true; popweight = true; popwghmode="district"
expenditureMode = false

incomeMode = true; relativeMode=true
religionMode = true
incomeByReligionMode = true
expenditureRangeMode = false
emissionLevelMode = false

costEstimationMode = false
costEstimationByThresholdMode = false
costEstimationByReligionMode = false

bootstrapMode = false
violinPlotting = false
stackedBarMode = false
bubbleChartMode = false
emissionByExp_plotting = false

mpcePeriod = "daily"

weightTag = ["popW", "hhW", "popWhhW", "perCap", "perHH", "demography"]
normTag = ["perCap", "perHH", "demography"]
categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
subcat=""
# subcat="Food"
foodCategories=["Grain","Vegetable","Fruit","Dairy","Beef","Pork","Poultry","Other meat","Fish",
                "Alcohol","Other beverage","Confectionery","Restaurant","Other food","Food"]

print(" Data reading: ")
print("category")
ec.readCategoryData(nation, sectorFile, except=["None"],subCategory=subcat)
if length(subcat)==0; ec.setCategory(categories)
else ec.setCategory(foodCategories); subcat*="_"
end
print(", household")
ec.readHouseholdData(year, householdFile, mergingMode, period=mpcePeriod)
print(", emission")
if !expenditureMode; ec.readEmission(year, emissionFile)
elseif expenditureMode ec.readExpenditure(year, Base.source_dir()*"/data/extracted/Expend_Matrix.txt")
end
println(" ... complete")

print(" Weights calculating: ")
print("state")
statePopulationFile = Base.source_dir()*"/data/statistics/StatePopulation.csv"
ec.calculateStatePopulationWeight(statePopulationFile)
print(", district")
districtPopulationFile = Base.source_dir()*"/data/statistics/DistrictPopulation.csv"
ec.calculateDistrictPopulationWeight(districtPopulationFile, sectorFile)
println(" ... complete")

print(" Categorizing:")
print(" category")
if weightMode>0; tag = weightTag[weightMode]; else tag="non" end
if expenditureMode; tag *= "_exp" end

hhsEmissionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cat.csv"
DistEmissionFile = Base.source_dir() * "/data/emission/2011_IND_dist_"*subcat*"emission_cat_"*tag*".csv"
ec.categorizeHouseholdEmission(year, output=hhsEmissionFile, hhsinfo=true, wghmode=popwghmode)
ec.calculateDistrictPoverty(year, povline=1.9, popWgh=popweight, period=mpcePeriod)
eData = ec.categorizeDistrictEmission(year, weightMode, sqrRoot=eqvalMode, period=mpcePeriod, religion=true, popWgh=popweight)
    # Period for MPCE: "annual", "monthly"(default), or "daily"
if weightMode == 4; overallMode = true else overallMode = false end
ec.printEmissionByDistrict(year, DistEmissionFile,eData[7],eData[6],name=true,totm=overallMode,expm=true,popm=true,hhsm=false,relm=true,wghm=true,denm=true,povm=true)

#ec.plotHHsEmission(year)

if exportMode || exportWebMode || mapStyleMode
    print(", exporting")
    #gidTag = "GID_2"
    gidTag = "Censuscode"
    exportFile = Base.source_dir() * "/data/emission/2011_IND_dist_GIS_"*subcat*"emission_cat_"*tag*".csv"
    exportRateFile = Base.source_dir() * "/data/emission/2011_IND_dist_GIS_"*subcat*"emission_cat_dr_"*tag*".csv"
    exportRankFile = Base.source_dir() * "/data/emission/2011_IND_dist_GIS_"*subcat*"emission_cat_rnk_"*tag*".csv"
    exportOrderFile = Base.source_dir() * "/data/emission/2011_IND_dist_GIS_"*subcat*"emission_cat_ord_"*tag*"_gr.csv"
    gData = ec.exportDistrictEmission(year, gidTag, exportFile, weightMode, logarithm=false, descend=true, empty=true, minmax=minmaxv)
    drData = ec.exportEmissionDiffRate(year, gidTag, exportRateFile, 0.5, -0.5, 128, descend=true, empty=true)
#    ec.exportEmissionValGroup(year, gidTag, exportRankFile, 128, descend=true, logscl=false)
#    ec.exportEmissionRankGroup(year, gidTag, exportOrderFile, 128, descend=true)
end
if exportWebMode
    print(", web-file exporting")
    exportPath = Base.source_dir() * "/data/emission/webfile/"
    ec.exportWebsiteFiles(year,exportPath,weightMode,gData[2],gData[5],gData[6],gData[3],gData[7],rank=true,empty=true)
end
if mapStyleMode
    print(", map-style file generating")
    if weightMode==1||weightMode==2; rgbFile = "../GIS/data/MPL_YlGnBu.rgb"
    elseif weightMode==4||weightMode==5; rgbFile = "../GIS/data/MPL_RdBu.rgb"
    end
    for i=1:length(ec.catList)
        qse.readColorMap(rgbFile, reverse=colormapReverse)
        qmlFile = replace(rgbFile, ".rgb"=>"_"*tag*"_"*ec.catList[i]*".qml")
        if weightMode==1||weightMode==2
            attr = "2011_IND_dist_GIS_emission_cat_"*tag*"_gr_"*ec.catList[i]
            qse.makeQML(qmlFile, attr, empty=true, labels=gData[9][:,i], indexValue=true, labelReverse=labeRev)
        elseif weightMode==4||weightMode==5
            attr = "2011_IND_dist_GIS_emission_cat_dr_"*tag*"_gr_"*ec.catList[i]
            qse.makeQML(qmlFile, attr, empty=true, values=drData[4][:,i], indexValue=true, labelReverse=labeRev)
        end
    end
end

if incomeMode
    print(" income")
    incomeFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_inc_"*tag*".csv"
    intervals = [0.2,0.4,0.6,0.8,1.0]; absint=false; descendig=false
    # intervals = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]; absint=false; descendig=false
    # intervals = [150,30]; absint = true
    eData = ec.categorizeHouseholdByIncome(year,intervals,normMode,popWgh=popweight,sqrRt=eqvalMode,absIntv=absint,perCap=percapita,desOrd=descendig,wghmode=popwghmode)
    ec.printEmissionByIncome(year, incomeFile, intervals, eData[4], eData[5], eData[6], absIntv=absint, relative=relativeMode)
end
if religionMode
    print(" religion")
    religionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_rel_"*tag*".csv"
    eData = ec.categorizeHouseholdByReligion(year, normMode, sqrRt=eqvalMode, popWgh=popweight, wghmode=popwghmode)
    ec.printEmissionByReligion(year, religionFile, eData[4], eData[5], eData[6])
end
if incomeByReligionMode
    print(" income-religion")
    intervals = [0.2,0.4,0.6,0.8,1.0]; absint=false; descendig=false
    incomeReligionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_incByRel_"*tag*".csv"
    eData = ec.categorizeHouseholdByIncomeByReligion(year,intervals,normMode,popWgh=popweight,sqrRt=eqvalMode,absIntv=absint,perCap=percapita,desOrd=descendig,wghmode=popwghmode)
    ec.printEmissionByIncomeByReligion(year,incomeReligionFile,intervals,eData[4],eData[5],eData[6],absIntv=absint,desOrd=descendig)
end
if expenditureRangeMode
    print(" expenditure-range")
    ranges = [1.25, 1.9, 3.0, 4.0, 5.0, 10.0]; absint=true; absSpan=true
    expRngFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_rng_"*tag*".csv"
    eData = ec.categorizeHouseholdByExpRange(year,ranges,normMode,over=0.1,less=0.1,absRng=absint,absSpn=absSpan,perCap=percapita,popWgh=popweight,wghmode=popwghmode)
    ec.printEmissionByRange(year,expRngFile,eData[6], eData[3], eData[4], eData[5], eData[7])
end
if emissionLevelMode
    print(" district")
    intervals = [0.1, 0.8, 0.1]; absint = false
    districlevelFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_dis_"*tag*".csv"
    ec.categorizeDistrictByEmissionLevel(year, normMode, intervals)
    ec.categorizeHouseholdByEmissionLevel(year, intervals, normMode, squareRoot=eqvalMode, absintv=absint)
    ec.printEmissionByDistEmLev(year, districlevelFile, intervals)
end
println(" ... complete")


if costEstimationMode||costEstimationByThresholdMode||costEstimationByReligionMode; print(" Cost estimating:") end
if costEstimationMode
    print(" district")
    intervals = [0.2,0.4,0.6,0.8,1.0]; absint=false; descendig=false; nameMode=true
    costDistrictFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cost_dis.csv"
    expDistrictFile = [Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cost_dis_gis.csv", "Censuscode"]
    ec.estimateEmissionCostByDistrict(year, intervals, normMode, perCap=percapita, popWgh=popweight, desOrd=descendig, name=nameMode, output=costDistrictFile, exportFile=expDistrictFile, wghmode=popwghmode)
end
if costEstimationByThresholdMode
    print(" threshold")
    thresholds = [1.9, 3.0, 5.0]; overRange=0.1;lessRange=0.1; absint=false; absSpan=true; descendig=false; nameMode=true; stackedMode=false
    costThresholdFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cost_threshold.csv"
    expThresholdFile = [Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cost_threshold_gis.csv", "Censuscode"]
    ec.estimateEmissionCostByDistrictForThreshold(year,thresholds,normMode,perCap=percapita,stacked=stackedMode,absSpn=absSpan,over=overRange,less=lessRange, popWgh=popweight, desOrd=descendig, name=nameMode, output=costThresholdFile, exportFile=expThresholdFile, wghmode=popwghmode)
end
if costEstimationByReligionMode
    print(" religion")
    intervals = [0.2,0.4,0.6,0.8,1.0]; absint=false; descendig=false; nameMode=true
    costReligionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cost_rel.csv"
    expReligionFile = [Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cost_rel_gis.csv", "Censuscode"]
    ec.estimateEmissionCostByDistrictByReligion(year, intervals, normMode, absIntv=absint, perCap=percapita, popWgh=popweight, desOrd=descendig, name=nameMode, output=costReligionFile, exportFile=expReligionFile, wghmode=popwghmode)
end
if costEstimationMode||costEstimationByThresholdMode||costEstimationByReligionMode; println(" ... complete") end

if bootstrapMode
    print(" Bootstrap proceeding: ")
    intervals = [0.2,0.4,0.6,0.8,1.0]
    cd(Base.source_dir())
    plotFile = "../Plot/chart/Emission_Bootstrap.png"
    plotFile = "../Plot/chart/Emission_District_Bootstrap.png"
    se.migrateData(year, ec)
    se.proceedExpenditureBootstrap(year, eData, intervals; perCap=false, disp=true, output=plotFile)
    #se.proceedDistrictBootstrap(year, eData; perCap=false, disp=true, output=plotFile)
    println(" ... complete")

end

if violinPlotting
    print(" Emission plotting: ")
    intervals = [0.2,0.4,0.6,0.8,1.0]
    cd(Base.source_dir())
    plotFile = "../Plot/chart/Emission_ViolinPlot.png"
    ep.migrateData(year, ec)
    ep.plotExpCatViolin(year, eData, intervals; perCap=true, boxplot=true, disp=true, output=plotFile)
    println(" ... complete")
end

if stackedBarMode
    print(" Stacked bar plotting: ")
    cd(Base.source_dir())
    emissionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_inc_perCap.csv"
    hhsFile = Base.source_dir() * "/data/emission/2011_IND_hhs_"*subcat*"emission_cat.csv"
    chartFile = "../Plot/chart/Stacked_bar_chart.png"
    paletteFile = "../Plot/Table color palettes.txt"
    ep.migrateData(year, ec)
    ep.readColorPalette(paletteFile, rev=true, tran=true)
    ep.plotExpStackedBarChart(emissionFile, chartFile, disp=true, hhsCF=hhsFile)
    println("completed")
end

if bubbleChartMode
    print(" Bubble chart plotting: ")
    cd(Base.source_dir())
    chartFile = "../Plot/chart/Bubble_chart.png"
    dataFile = "../Plot/chart/Bubble_chart.csv"
    ep.migrateData(year, ec)
    ep.plotCfBubbleChart(year, chartFile, disp=true, dataoutput=dataFile, povline=1.9)
    println("completed")
end

if emissionByExp_plotting
    print(" Plotting: ")
    print("emission by expenditure")
    outputFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_cost.csv"
    efc.printEmissionByExp(year, outputFile, period="daily", percap=false, plot=false, dispmode=false, guimode=false)
    println(" ... complete")
end

println("[Done]")
