# Developed date: 27. Dec. 2019
# Last modified date: 6. Mar. 2020
# Subject: Emission mapping
# Description: Mapping emission through households emissions data
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("EmissionCategorizer.jl")
using .EmissionCategorizer
ec = EmissionCategorizer

println("[Process]")

nation = "IND"
year = 2011
emissionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission.txt"
householdFile = Base.source_dir() * "/data/extracted/Households.txt"
sectorFile = "../Eora/data/Eora_HS_match.xlsx"

mergingMode = true # true: proceed district merging, default=false

weightMode = 4  # [0]non-weight, [1]population weighted, [2]household weighted, [3]both population and household weighted
                # ([4],[5]: normalization) [4]per capita, [5]per household
                # (basic information) [6]population and households, [1,:]population, [2,:]households
normMode = 1    # [0]non-weight, [1]per capita, [2]per houehold,
                # (basic information) [3]population and households by religions, [1,:]population, [2,:]households
eqvalMode = false   # [true]apply square root of household size for equivalance scale

categorizeMode = true
exportMode = false
exportWebMode = false

districtMode = false
religionMode = false
incomeMode = true
emissionByExp_plotting = false

weightTag = ["popW", "hhW", "popWhhW", "perCap", "perHH", "demography"]
normTag = ["perCap", "perHH", "demography"]

print(" Data reading: ")
print("category")
ec.readCategoryData(nation, sectorFile)
print(", household")
ec.readHouseholdData(year, householdFile, mergingMode)
print(", emission")
ec.readEmission(year, emissionFile)
println(" ... complete")

print(" Categorizing:")
if categorizeMode
    print(" category")
    if weightMode>0; tag = weightTag[weightMode]; else tag="non" end
    hhsEmissionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_cat.csv"
    categorizedFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_cat_"*tag*".txt"
    ec.categorizeHouseholdsEmission(year, output = hhsEmissionFile, hhsinfo=true)
    ec.categorizeEmission(year, weightMode, eqvalMode)
    ec.printCategorizedEmission(year, categorizedFile, name=true)

    if exportMode || exportWebMode
        exportFile = Base.source_dir() * "/data/emission/2011_IND_hhs_GIS_emission_cat_"*tag*".csv"
        exportRateFile = Base.source_dir() * "/data/emission/2011_IND_hhs_GIS_emission_cat_dr_"*tag*".csv"
        exportRankFile = Base.source_dir() * "/data/emission/2011_IND_hhs_GIS_emission_cat_rnk_"*tag*".csv"
        exportOrderFile = Base.source_dir() * "/data/emission/2011_IND_hhs_GIS_emission_cat_ord_"*tag*"_gr.csv"
        print(", exporting")
        gData = ec.exportEmissionTable(year, "GID_2", exportFile, weightMode, false)
        ec.exportEmissionDiffRate(year, "GID_2", exportRateFile, false, 0.5, -0.5, 128, true)
        ec.exportEmissionValGroup(year, "GID_2", exportRankFile, false, 128, descend=true, logscl=false)
        ec.exportEmissionRankGroup(year, "GID_2", exportOrderFile, false, 128, descend=true)
    end
    if exportWebMode
        exportPath = Base.source_dir() * "/data/emission/webfile/"
        print(", web-file exporting")
        ec.exportWebsiteFiles(year,exportPath,weightMode,gData[2],gData[5],gData[6],gData[3],gData[7],rank=true,empty=true)
    end
end
if districtMode
    print(" district")
    intervals = [0.1, 0.8, 0.1]
    districlevelFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_dis_"*tag*".txt"
    eData = ec.categorizeEmission(year, 0, eqvalMode)
    ec.categorizeDistrict(year, eData[4], normMode, eData[6], eData[7], intervals)
    ec.printCategorizedDistrict(year, districlevelFile,intervals)
end
if religionMode
    print(" religion")
    religionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_rel_"*tag*".txt"
    ec.categorizeReligion(year, normMode, eqvalMode)
    ec.printCategorizedReligion(year, religionFile)
end
if incomeMode
    print(" income")
    incomeFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_inc_"*tag*".txt"
    intervals = [0.1, 0.8, 0.1]; absint = false
    #intervals = [150,30]; absint = true
    eData = ec.categorizeIncome(year, intervals, normMode, eqvalMode, absintv=absint, percap=true)
    ec.printCategorizedIncome(year, incomeFile, intervals, eData[4], eData[5], absintv=absint)
end
if emissionByExp_plotting
    print(" Plotting: ")
    print("emission by expenditure")
    outputFile = Base.source_dir() * "/data/emission/2011_IND_food_emission_byExpenditure.txt"
    efc.printEmissionByExp(year, outputFile, period="daily", percap=false, plot=false, dispmode=false, guimode=false)
    println(" ... complete")
end

println(" ... complete")
