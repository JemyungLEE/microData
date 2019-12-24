# Developed date: 5. Dec. 2019
# Last modified date: 24. Dec. 2019
# Subject: Analyze carbon emissions by final demands of Eora and Comtrade data
# Description: Calculate carbon emissions by utilizing Eora T, V, Y, and Q tables.
#              Commodity sectors: Eora or Comtrade data, Service sectors: Eora data
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

#clearconsole()
cd(Base.source_dir())

include("MicroDataReader.jl")
include("EmissionEstimator.jl")
include("ImportTransformer.jl")
include("../Comtrade/HsDataReader.jl")
include("../Comtrade/TradeMatrixBuilder.jl")
include("../converting/HsConcMatBuilder.jl")
include("../converting/XLSXextractor.jl")

using .MicroDataReader
using .EmissionEstimator
using .ImportTransformer
using .HsDataReader
using .TradeMatrixBuilder
using .HsConcMatBuilder
using .XLSXextractor

mdr = MicroDataReader
ee = EmissionEstimator
it = ImportTransformer
hdr = HsDataReader
tmb = TradeMatrixBuilder
hb = HsConcMatBuilder
xls = XLSXextractor

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(Base.source_dir()*"/index/ProductCategory.txt")
println("complete")

# India Household Expenditure micro-data reading process

tag = "T1_"
path = Base.source_dir()*"/type_1/"
hhdata = []
push!(hhdata, [path*"test_lv1.txt", [36, 24, 2, 37, 38, 6]])   # level 1: household data
push!(hhdata, [path*"test_lv2.txt", [40, 20]])                 # level 2
push!(hhdata, [path*"test_lv3.txt", [35, 28, 29]])             # level 3
push!(hhdata, [path*"test_lv4.txt", [38, 23, 22, 24, 25]])     # level 4
microdata = []
push!(microdata, [path*"test_lv5.txt", [31, 20, 24, 23, 30, 25, 22, 21]])  # level 5
push!(microdata, [path*"test_lv6.txt", [30, 20, 24, 23, 365]])  # level 6, Please divide Last_365days_Quantity by 1000 to get figures in (0.000) for consumption of clothing, bedding etc. and consumption of footwear during last 365 days is in no. of pairs.
push!(microdata, [path*"test_lv7.txt", [28, 20, 22, -1, 365]])  # level 7
push!(microdata, [path*"test_lv8.txt", [27, 20, 21, -1, 30]])   # level 8
push!(microdata, [path*"test_lv9.txt", [40, 20, 34, -1, 365]])  # level 9

print(" Household data reading: $tag")
mdr.readHouseholdData(hhdata, tag)
println("complete")
print(" Expenditure data reading: $tag")
mdr.readMicroData(microdata, tag)
println("complete")

tag = "T2_"
path = Base.source_dir()*"/type_2/"
hhdata = []
push!(hhdata, [path*"test_lv1.txt", [33, 23, 2, 34, 35, 6]])   # level 1: household data
push!(hhdata, [path*"test_lv2.txt", [37, 19]])                 # level 2
push!(hhdata, [path*"test_lv3.txt", [31, 27, 27]])             # level 3, Please divide MPCE by 100 to get figures in ( Rs. 0.00)
push!(hhdata, [path*"test_lv4.txt", [35, 22, 21, 23, 24]])     # level 4
microdata = []
push!(microdata, [path*"test_lv5.txt", [28, 19, 23, 22, 30, 24, 21, 20]])  # level 5
push!(microdata, [path*"test_lv6.txt", [25, 19, 21, 20, 365]])  # level 6, Please divide Last_365days_Quantity by 1000 to get figures in (0.000) for consumption of clothing, bedding etc. and consumption of footwear during last 365 days is in no. of pairs.
push!(microdata, [path*"test_lv7.txt", [24, 19, 20, -1, 365]])  # level 7
push!(microdata, [path*"test_lv8.txt", [24, 19, 20, -1, 30]])   # level 8
push!(microdata, [path*"test_lv9.txt", [32, 20, 28, -1, 365]])  # level 9
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

print(" Expenditure matrix building: ")
expData = mdr.makeExpenditureMatrix()   # [1]:expenditure matrix(hhid, sec), [2]:hhid, [3]: Indi sectors
ee.getDomesticData(transpose(expData[1]), expData[2], expData[3])
println("complete")

print(" Memory securing: ")
mdr.initVars()
println("complete")

# Converting process of Eora final demand data to India micro-data format
path = "../converting/data/"
concordanceFile = path * "India(STAT) vs EORA_Ver1.2.xlsx"

print(" Concordance matrix building: Eora ...")
xls.readXlsxData(concordanceFile, "India")
xls.buildConMat()
cmn = xls.normConMat()   # {a3, {conMat, sumEora, sumNat}}
println("complete")

# Eora household's final-demand import sector data reading process
year = 2011
nation = "IND"
path = "../Eora/data/" * string(year)
eoraIndexFile = "../Eora/data/Eora_HS_match.xlsx"

print(" Eora index reading: ")
ee.readIndexXlsx(eoraIndexFile)
println("complete")

print(" MRIO table reading: ")
ee.readIOTables(year, path*"_eora_t.csv", path*"_eora_v.csv", path*"_eora_y.csv", path*"_eora_q.csv")
ee.rearrangeIndex()
ee.rearrangeTables(year)
println("complete")

print(" Weighted concordance matrix building: ")
ee.buildWeightedConcMat(year, nation, cmn)
println("complete")

print(" Emission calculating: ")
path = Base.source_dir()*"/data/emission/"
emissionFile = path * string(year) * "_" * nation * "_hh_emission.txt"
ee.calculateEmission(year, false, 1)
ee.printEmissions(year, emissionFile)
println("complete")

println(" ... all complete")

#=
# Import accounts analysis part
print(" Expenditure analysis: ")
it.analyzeExpenditures(expData[1], expData[2], expData[3])
println("complete")

path = Base.source_dir()*"/data/emission/"
print(" Emission data trasform: Eora ...")
it.transformEORAtoIND(ceData[3], ceData[2], ceData[1], cmn, ceData[4])
it.printTransfImport(path*"Emission_tranformed_Eora_to_IND.txt")
println("complete")

# Households' share calculation process
print(" Household share calculation: ")
hhShareFile = path*"ImportHH_tranformed_fromEor_toIND.txt"
it.calculateHHtotal(hhShareFile, false, true)
=#
