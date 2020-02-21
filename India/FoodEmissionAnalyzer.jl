# Developed date: 18. Feb. 2019
# Last modified date: 20. Feb. 2020
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
sectorFile = "../Eora/data/Eora_HS_match.xlsx"

mergingMode = true # true: proceed district merging, default=false

normMode = 2    # [0]non-weight, [1]per capita, [2]per houehold,
                # (basic information) [3]population and households by religions, [1,:]population, [2,:]households
eqvalMode = true   # [true]apply square root of household size for equivalance scale

single_categorizing = false
multi_categorizing = false
emissionByExp_plotting = true

normTag = ["perCap", "perHH", "demography"]

print(" Data reading: ")
print("category")
ec.readCategoryData(nation, sectorFile)
print(", household")
ec.readHouseholdData(year, householdFile, mergingMode)
print(", emission")
ec.readEmission(year, emissionFile)
println(" ... complete")

print(" Data migrating: ")
efc.migrateData(year, ec)
efc.readCategoryData(nation, sectorFile)
println("complete")

print(" Categorizing:")
tag = normTag[normMode]
outputFile = Base.source_dir() * "/data/emission/2011_IND_food_emission_"
print(" households")
efc.categorizeEmissionHouseholds(year)
efc.printEmissionHHs(year, outputFile*"hhs_"*tag*".csv")

if single_categorizing
    print(", religion")
    efc.categorizeEmissionReligion(year, eqvalMode)
    efc.printEmissionRel(year, outputFile*"rel_"*tag*".csv")
    print(", district")
    efc.categorizeEmissionDistrict(year, eqvalMode)
    efc.printEmissionDis(year, outputFile*"dis_"*tag*".csv", false)
    intervals = [0.1, 0.2, 0.4, 0.2, 0.1]
    print(", income")
    efc.categorizeEmissionIncome(year, intervals, normMode, eqvalMode)
    efc.printEmissionInc(year, outputFile*"inc_"*tag*".csv", intervals)
    print(", emission")
    efc.categorizeEmissionLevel(year, intervals, normMode, eqvalMode)
    efc.printEmissionLev(year, outputFile*"lev_"*tag*".csv", intervals)
end
if multi_categorizing
    intervals = [0.1, 0.2, 0.4, 0.2, 0.1]
    print(", religion-income")
    efc.categorizeEmissionReligionIncome(year, intervals, normMode, eqvalMode)
    efc.printEmissionRelInc(year, outputFile*"RelInc_"*tag*".csv", intervals)
    print(", religion-emission_level")
    efc.categorizeEmissionReligionLevel(year, intervals, normMode, eqvalMode)
    efc.printEmissionRelLev(year, outputFile*"RelLev_"*tag*".csv", intervals)
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
