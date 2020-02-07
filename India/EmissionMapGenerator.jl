# Developed date: 27. Dec. 2019
# Last modified date: 7. Feb. 2019
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

print(" Data reading: ")
print("household")
ec.readHouseholdData(year, householdFile)
print(", category")
ec.readCategoryData(nation, sectorFile)
print(", emission")
ec.readEmission(year, emissionFile)
println(" ... complete")

weightMode = 4  # [0]non-weight, [1]population weighted, [2]household weighted, [3]both population and household weighted
                # ([4],[5]: normalization) [4]per capita, [5]per household
                # (basic information) [6]population and households, [1,:]population, [2,:]households
normMode = 1    # [0]non-weight, [1]per capita, [2]per houehold,
                # (basic information) [3]population and households by religions, [1,:]population, [2,:]households
eqvalMode = false # [true]apply square root of household size for equivalance scale
categorizeMode = true
exportMode = false
districtMode = false
religionMode = false
incomeMode = false

weightTag = ["popW", "hhW", "popWhhW", "perCap", "perHH", "demography"]
normTag = ["perCap", "perHH", "demography"]

print(" Categorizing:")
if categorizeMode
    print(" category")
    if weightMode>0; tag = weightTag[weightMode]; else tag="non" end
    categorizedFile = Base.source_dir() * "/data/emission/2011_IND_hhs_emission_cat_"*tag*".txt"
    ec.categorizeEmission(year, weightMode, eqvalMode)
    ec.printCategorizedEmission(year, categorizedFile, true)

    if exportMode
        exportFile = Base.source_dir() * "/data/emission/2011_IND_hhs_GIS_emission_cat_"*tag*".csv"
        print(" exporting")
        ec.exportEmissionTable(year, "GID_2", exportFile, true, nation, sectorFile)
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
