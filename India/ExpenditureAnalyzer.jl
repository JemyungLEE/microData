# Developed date: 22. Jan. 2020
# Last modified date: 4. Feb. 2020
# Subject: Analyze household consumer expenditure
# Description: Calculate household expenditures by hh size and by categorizes
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("MicroDataReader.jl")
include("ExpenditureCategorizer.jl")

using .MicroDataReader
using .ExpenditureCategorizer

mdr = MicroDataReader
ec = ExpenditureCategorizer

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(Base.source_dir()*"/data/index/ProductCategory.txt")
println("complete")

# India Household Expenditure micro-data reading process

tag = "T1_"
path = Base.source_dir()*"/data/type_1/"
hhdata = []
push!(hhdata, [path*"test_lv1.txt", [36, 24, 2, 37, 38, 6]])   # level 1: household data
push!(hhdata, [path*"test_lv2.txt", [40, 20, 24]])             # level 2
push!(hhdata, [path*"test_lv3.txt", [35, 29]])                 # level 3, Please divide MPCE by 100 to get figures in ( Rs. 0.00)
push!(hhdata, [path*"test_lv4.txt", [38, 23, 22, 24, 25]])     # level 4
microdata = []
push!(microdata, [path*"test_lv5.txt", [31, 20, 24, 23, 30, 25, 22, 21]])  # level 5, Rs., quantity for 30 and 365 days (0.000)
push!(microdata, [path*"test_lv6.txt", [30, 20, 24, 23, 365]])  # level 6, Rs., Please divide Last_365days_Quantity by 1000 to get figures in (0.000) for consumption of clothing, bedding etc. and consumption of footwear during last 365 days is in no. of pairs.
push!(microdata, [path*"test_lv7.txt", [28, 20, 22, -1, 365]])  # level 7
push!(microdata, [path*"test_lv8.txt", [27, 20, 21, -1, 30]])   # level 8, Rs.
push!(microdata, [path*"test_lv9.txt", [40, 20, 34, -1, 365]])  # level 9, Rs.
print(" Household data reading: $tag")
mdr.readHouseholdData(hhdata, tag)
println("complete")
print(" Expenditure data reading: $tag")
mdr.readMicroData(microdata, tag)
println("complete")

tag = "T2_"
path = Base.source_dir()*"/data/type_2/"
hhdata = []
push!(hhdata, [path*"test_lv1.txt", [33, 23, 2, 34, 35, 6]])   # level 1: household data
push!(hhdata, [path*"test_lv2.txt", [37, 19, 23]])             # level 2
push!(hhdata, [path*"test_lv3.txt", [31, 27]])                 # level 3, Please divide MPCE by 100 to get figures in ( Rs. 0.00)
push!(hhdata, [path*"test_lv4.txt", [35, 22, 21, 23, 24]])     # level 4
microdata = []
push!(microdata, [path*"test_lv5.txt", [28, 19, 23, 22, 30, 24, 21, 20, [7,180,329]]])  # level 5, Rs., quantity 0.000, Reference period for last 7 days: Edible oil; egg, fish & meat; vegetables, fruits, spices, beverages and processed foods; pan, tobacco & intoxicants
push!(microdata, [path*"test_lv6.txt", [25, 19, 21, 20, 365]])  # level 6, Rs., Please divide Last_365days_Quantity by 1000 to get figures in (0.000) for consumption of clothing, bedding etc. and consumption of footwear during last 365 days is in no. of pairs.
push!(microdata, [path*"test_lv7.txt", [24, 19, 20, -1, 365]])  # level 7, Rs.
push!(microdata, [path*"test_lv8.txt", [24, 19, 20, -1, 30]])   # level 8, Rs.
push!(microdata, [path*"test_lv9.txt", [32, 20, 28, -1, 365]])  # level 9, Rs.
print(" Household data reading: $tag")
mdr.readHouseholdData(hhdata, tag)
println("complete")
print(" Expenditure data reading: $tag")
mdr.readMicroData(microdata, tag)
println("complete")

print(" Currency exchanging: ")
exchangeRate = 0.01888      # 2011-12-26, Indian Rupee to USD
mdr.currencyExchange(exchangeRate)
println("complete")

year = 2011
nat = "IND"

print(" Expenditure analyzing: ")
expData = mdr.makeExpenditureMatrix()   # [1]:expenditure matrix(hhid, sec), [2]:hhid, [3]: India sectors
ec.getExpenditureData(year, expData)
ec.getHouseholdData(year, mdr.households)
println("complete")

print(" Expenditure categorizing: ")
eoraIndexFile = "../Eora/data/Eora_HS_match.xlsx"
ec.readCategoryData(nat, eoraIndexFile)
ec.categorizeExpenditure(year)
println("complete")

print(" Households by expenditure counting: ")
countingFile = Base.source_dir()*"/data/expenditure/"*string(year)*"_"*nat*"_count.txt"
maxexp = [800,800,400,100,3200,100,500,2,200,800,300,300,8000]
print("counting"); expData = ec.countByExpenditure(year, 20, maxexp, [], 20)
print(", non-parametric regression "); ec.nonparreg(year, expData[5], false)
ec.printCountedResult(year, countingFile, expData[4], expData[5], expData[6], expData[7], expData[8])
println("complete")

print(" Expenditure heatmap printing: ")
heatmapFile = Base.source_dir()*"/data/expenditure/"*string(year)*"_"*nat*"_heatmap"
plt = ec.plotHeatmap(year, expData[4], expData[5], true, false, false)
println("complete")


print(" Results printing: ")
path = Base.source_dir()*"/data/expenditure/"
ec.printCategorizedExpenditureByHHsize(path*string(year)*"_"*nat*"_categorized_expenditure.txt")
ec.printCategorizedExpenditureByHH(year,path*string(year)*"_"*nat*"_categorized_expenditure_HH.txt")
println("complete")
