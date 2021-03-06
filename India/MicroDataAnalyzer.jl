# Developed date: 21. Oct. 2019
# Last modified date: 12. Feb. 2021
# Subject: India microdata analyzer
# Description: proceed data analysis process for India household consumption microdata
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
include("MicroDataReader.jl")
include("PovertyMap.jl")
using .MicroDataReader
using .PovertyMap
mdr = MicroDataReader
pm = PovertyMap
catFile = Base.source_dir()*"/data/index/ProductCategory.txt"

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(catFile)
println("completed")

#for type 1
tag = "T1_"
path = Base.source_dir()*"/data/type_1/"
hhdata = []
push!(hhdata, [path*"test_lv1.txt", [36, 24, 2, 37, 38, 6, [7, 8, 9, 10, 14]]])   # level 1: household data
push!(hhdata, [path*"test_lv2.txt", [40, 20, 24]])             # level 2
push!(hhdata, [path*"test_lv3.txt", [35, 29]])                 # level 3, Please divide MPCE by 100 to get figures in ( Rs. 0.00)
push!(hhdata, [path*"test_lv4.txt", [38, 23, 22, 24, 25, 21]])     # level 4
microdata = []
push!(microdata, [path*"test_lv5.txt", [31, 20, 24, 23, 30, 25, 22, 21]])  # level 5, Rs., quantity for 30 and 365 days (0.000)
push!(microdata, [path*"test_lv6.txt", [30, 20, 24, 23, 365]])  # level 6, Rs., Please divide Last_365days_Quantity by 1000 to get figures in (0.000) for consumption of clothing, bedding etc. and consumption of footwear during last 365 days is in no. of pairs.
push!(microdata, [path*"test_lv7.txt", [28, 20, 22, -1, 365]])  # level 7
push!(microdata, [path*"test_lv8.txt", [27, 20, 21, -1, 30]])   # level 8, Rs.
push!(microdata, [path*"test_lv9.txt", [40, 20, 34, -1, 365]])  # level 9, Rs.
print(" Household data reading: $tag")
mdr.readHouseholdData(hhdata, tag)
println("completed")
print(" Expenditure data reading: $tag")
mdr.readMicroData(microdata, tag)
println("completed")

#for type 2
tag = "T2_"
path = Base.source_dir()*"/data/type_2/"
hhdata = []
push!(hhdata, [path*"test_lv1.txt", [33, 23, 2, 34, 35, 6, [7, 8, 9, 10, 14]]])   # level 1: household data
push!(hhdata, [path*"test_lv2.txt", [37, 19, 23]])             # level 2
push!(hhdata, [path*"test_lv3.txt", [31, 27]])                 # level 3, Please divide MPCE by 100 to get figures in ( Rs. 0.00)
push!(hhdata, [path*"test_lv4.txt", [35, 22, 21, 23, 24, 20]])     # level 4
microdata = []
push!(microdata, [path*"test_lv5.txt", [28, 19, 23, 22, 30, 24, 21, 20, [7,180,329]]])  # level 5, Rs., quantity 0.000, Reference period for last 7 days: Edible oil; egg, fish & meat; vegetables, fruits, spices, beverages and processed foods; pan, tobacco & intoxicants
push!(microdata, [path*"test_lv6.txt", [25, 19, 21, 20, 365]])  # level 6, Rs., Please divide Last_365days_Quantity by 1000 to get figures in (0.000) for consumption of clothing, bedding etc. and consumption of footwear during last 365 days is in no. of pairs.
push!(microdata, [path*"test_lv7.txt", [24, 19, 20, -1, 365]])  # level 7, Rs.
push!(microdata, [path*"test_lv8.txt", [24, 19, 20, -1, 30]])   # level 8, Rs.
push!(microdata, [path*"test_lv9.txt", [32, 20, 28, -1, 365]])  # level 9, Rs.
print(" Household data reading: $tag")
mdr.readHouseholdData(hhdata, tag)
println("completed")
print(" Expenditure data reading: $tag")
mdr.readMicroData(microdata, tag)
println("completed")


exchCurr = false
pppConv = true
povApply = false
weightMode = false

expMat = true
printMat = true

exchRateFile = Base.source_dir()*"/data/index/CurrencyExchangeRates.txt"
pppFile = Base.source_dir()*"/data/index/PPPs.txt"          # befor PPP revision at May, 2020
# pppFile = Base.source_dir()*"/data/index/PPPs_revised.txt"  # after PPP revision at May, 2020
#exchangeRate = 46.6226  # 2011 average exchange rate, USD to Indian Rupee
#ppp = 15.109            # 2011, India/USD

if pppConv; tag = "_ppp" else tag = "" end
povertyLineFile = Base.source_dir()*"/data/index/PovertyLine"*tag*".csv"
povertyOutputFile = Base.source_dir()*"/data/extracted/PovertyStatistics"*tag*".txt"
stateIndexFile = Base.source_dir()*"/data/index/IND_gis_index.xlsx"
povertyGIS_file = Base.source_dir()*"/data/extracted/IND_povert_st"*tag*".csv"

print(" Currency converting:")
if exchCurr; print(" expenditure"); mdr.exchangeExpCurrency(exchRateFile, inverse=true) end
if pppConv; print(" PPP"); mdr.convertMpceToPPP(pppFile) end
println("... complete")

if povApply
    print(" Poverty line applying: ")
    pm.migrateData(mdr)
    pm.applyPovertyLine(povertyLineFile, povertyOutputFile)
    pm.calculatePopulationWeight()
    pm.exportPovertyMap(stateIndexFile, povertyGIS_file)
    println("completed" )
end

if weightMode
    print(" Population weight calculating: ")
    print("state")
    statePopulationFile = Base.source_dir()*"/data/statistics/StatePopulation.csv"
    districtPopulationFile = Base.source_dir()*"/data/statistics/DistrictPopulation.csv"
    indexFile = Base.source_dir()*"/data/index/IND_index_match_v1.3.xlsx"
    mdr.calculateStatePopulationWeight(statePopulationFile)
    mdr.calculateDistrictPopulationWeight(statePopulationFile, indexFile)
    println(", completed" )
end

if expMat
    print(" Expenditure matrix building: ")
    expenditureMatrixFile = Base.source_dir()*"/data/extracted/Expend_Matrix.txt"
    householdDataFrameFile = Base.source_dir()*"/data/extracted/Household_DataFrame.txt"
    mdr.makeExpenditureMatrix(expenditureMatrixFile)
    println("completed")
    #print(" Household DataFrame building: ")
    #mdr.convertHouseholdData(householdDataFrameFile)
    #println("completed")
end

if printMat
    print(" Extracted data printing: ")
    householdsFile = Base.source_dir()*"/data/extracted/Households.txt"
    memberFile = Base.source_dir()*"/data/extracted/Members.txt"
    expenditureFile = Base.source_dir()*"/data/extracted/Expenditures.txt"
    mdr.printHouseholdData(householdsFile, addPov=povApply, addWgh=weightMode)
    mdr.printMemberData(memberFile)
    if expMat; mdr.printMicroData(expenditureFile) end
    println("completed")
end

println("[done]")
