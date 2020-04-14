# Developed date: 18. Feb. 2019
# Last modified date: 14. Apr. 2020
# Subject: Food carbon emission analysis
# Description: Calculate food sections' househld carbon emissions
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())


include("EmissionCategorizer.jl")
include("EmissionFoodCategorizer.jl")

using .EmissionCategorizer
using .EmissionFoodCategorizer

ec = EmissionCategorizer
efc = EmissionFoodCategorizer

println("[Process]")

nation = "IND"
year = 2011
emissionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission.txt"
householdFile = Base.source_dir() * "/data/extracted/Households.txt"
sectorFile = Base.source_dir() * "/data/index/IND_index_match_v1.3.xlsx"

mergingMode = true # true: proceed district merging, default=false

normMode = 1    # [0]non-weight, [1]per capita, [2]per houehold,
                # (basic information) [3]population and households by religions, [1,:]population, [2,:]households
eqvalMode = false   # [true]apply square root of household size for equivalance scale

incomeMode = true; percapita = true; popweight = true
religionMode = false
incomeByReligionMode = false
emissionLevelMode = false

emissionByExp_plotting = false

normTag = ["perCap", "perHH", "demography"]

print(" Data reading: ")
print("category")
ec.readCategoryData(nation, sectorFile, subCategory="Food")
print(", household")
ec.readHouseholdData(year, householdFile, mergingMode)
print(", emission")
ec.readEmission(year, emissionFile)
#ec.readExpenditure(year, Base.source_dir()*"/data/extracted/Expend_Matrix.txt")
println(" ... complete")

print(" Data migrating: ")
efc.migrateData(year, ec)
println("complete")

print(" Categorizing:")
tag = normTag[normMode]
distEmissionFile = Base.source_dir() * "/data/emission/2011_IND_dist_food_emission_"*tag*".csv"
print(" households")
efc.categorizeHouseholdEmission(year)
print(" districts")
eData = efc.categorizeDistrictEmission(year, sqrRoot=eqvalMode, period="daily", religion=true)
efc.printEmissionByDistrict(year,distEmissionFile,eData[7],eData[6],name=true,expm=true,popm=true,hhsm=false,relm=true)
#efc.printEmissionHHs(year, outputFile*"hhs_"*tag*".csv")
#efc.printHHsEmissionData(year, Base.source_dir() * "/data/emission/2011_IND_HHsEmissions.csv", sorting=true)

if incomeMode
    print(" income")
    incomeFile = Base.source_dir() * "/data/emission/2011_IND_hhs_food_emission_inc_"*tag*".csv"
#    intervals = [0.032208036488199,0.1,0.194936331961496,0.5,0.9,1.0]; absint=false; descendig=false
#    intervals = [0.25,0.5,0.75,1.00]; absint=false; descendig=false
    intervals = [0.2,0.4,0.6,0.8,1.0]; absint=false; descendig=false
#    intervals = [1/6,2/6,3/6,4/6,5/6,1.0]; absint=false; descendig=false
#    intervals = [0.2,0.8,1.0]; absint=false; descendig=false   # poverty line $1.9
    #intervals = [150,30]; absint = true
    eData = efc.categorizeHouseholdByIncome(year,intervals,normMode,popWgh=popweight,sqrRt=eqvalMode,absIntv=absint,perCap=percapita,desOrd=descendig)
    efc.printEmissionByIncome(year, incomeFile, intervals, eData[4], eData[5], eData[6], absIntv=absint)
end
if religionMode
    print(" religion")
    religionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_food_emission_rel_"*tag*".csv"
    eData = efc.categorizeHouseholdByReligion(year, normMode, sqrRt=eqvalMode, popWgh=popweight)
    efc.printEmissionByReligion(year, religionFile, eData[4], eData[5], eData[6])
end
if incomeByReligionMode
    print(" income-religion")
    intervals = [0.2,0.4,0.6,0.8,1.0]; absint=false; descendig=false
    incomeReligionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_food_emission_incByRel_"*tag*".csv"
    eData = efc.categorizeHouseholdByIncomeByReligion(year,intervals,normMode,popWgh=popweight,sqrRt=eqvalMode,absIntv=absint,perCap=percapita,desOrd=descendig)
    efc.printEmissionByIncomeByReligion(year,incomeReligionFile,intervals,eData[4],eData[5],eData[6],absIntv=absint,desOrd=descendig)
end

println(" ... complete")

if emissionByExp_plotting
    print(" Plotting: ")
    print("emission by expenditure")
    outputFile = Base.source_dir() * "/data/emission/2011_IND_food_emission_byExpenditure.txt"
    efc.printEmissionByExp(year, outputFile, period="daily", percap=false, plot=false, dispmode=false, guimode=false)
    println(" ... complete")
end

println("[Done]")
