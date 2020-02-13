# Developed date: 27. Dec. 2019
# Last modified date: 13. Feb. 2020
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
exportMode = true
exportWebMode = true
districtMode = false
religionMode = false
incomeMode = false

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
    categorizedFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_cat_"*tag*".txt"
    ec.categorizeEmission(year, weightMode, eqvalMode)
    ec.printCategorizedEmission(year, categorizedFile, true)

    if exportMode || exportWebMode
        exportFile = Base.source_dir() * "/data/emission/2011_IND_hhs_GIS_emission_cat_"*tag*".csv"
        exportRateFile = Base.source_dir() * "/data/emission/2011_IND_hhs_GIS_emission_cat_dr_"*tag*".csv"
        print(", exporting")
        gData = ec.exportEmissionTable(year, "GID_2", exportFile, weightMode, false)
        ec.exportEmissionDiffRate(year, "GID_2", exportRateFile, false, 0.5, -0.5, 128, true)
    end
    if exportWebMode
        exportPath = Base.source_dir() * "/data/emission/webfile/"
        print(", web-file exporting")
        ec.exportWebsiteFiles(year, exportPath, weightMode, gData[2], gData[5], gData[6], gData[3], gData[7])
    end
end
if districtMode
    print(" district")
    intervals = [0.1, 0.8, 0.1]
    districlevelFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_dis_perCap.txt"
    eData = ec.categorizeEmission(year, 0, eqvalMode)
    ec.categorizeDistrict(year, eData[4], normMode, eData[6], eData[7], intervals)
    ec.printCategorizedDistrict(year, districlevelFile,intervals)
end
if religionMode
    print(" religion")
    religionFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_rel_perCap.txt"
    ec.categorizeReligion(year, normMode, eqvalMode)
    ec.printCategorizedReligion(year, religionFile)
end
if incomeMode
    print(" income")
    incomeFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_inc_perCap.txt"
    intervals = [0.1, 0.2, 0.4, 0.2, 0.1]
    ec.categorizeIncome(year, intervals, normMode, eqvalMode)
    ec.printCategorizedIncome(year, incomeFile, intervals,)
end

println(" ... complete")
