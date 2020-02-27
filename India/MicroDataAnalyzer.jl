# Developed date: 21. Oct. 2019
# Last modified date: 26. Feb. 2020
# Subject: India microdata analyzer
# Description: proceed data analysis process for India household consumption microdata
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
include("MicroDataReader.jl")
using .MicroDataReader
mdr = MicroDataReader
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
push!(hhdata, [path*"test_lv4.txt", [38, 23, 22, 24, 25]])     # level 4

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
push!(hhdata, [path*"test_lv4.txt", [35, 22, 21, 23, 24]])     # level 4

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

print(" Currency exchanging: ")
exchangeRate = 1/46.6226        # 2011 average exchange rate, Indian Rupee to USD
ppp = 15.109                # 2011, India/USD
mdr.currencyExchange(exchangeRate, ppp)
println("complete")

path = Base.source_dir()*"/data/extracted/"
expenditureMatrixFile = path*"Expend_Matrix.txt"
householdDataFrameFile = path*"Household_DataFrame.txt"

print(" Expenditure matrix building: ")
mdr.makeExpenditureMatrix(expenditureMatrixFile)
println("completed")
print(" Household DataFrame building: ")
mdr.convertHouseholdData(householdDataFrameFile)
println("completed")


householdsFile = path*"Households.txt"
memberFile = path*"Members.txt"
expenditureFile = path*"Expenditures.txt"

mdr.printHouseholdData(householdsFile)
mdr.printMemberData(memberFile)
mdr.printMicroData(expenditureFile)


println("[completed]")
