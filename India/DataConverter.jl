# Developed date: 15. Nov. 2019
# Last modified date: 22. Nov. 2019
# Subject: Import data converter
# Description: Convert import data sectors of Comtrade and Eora to those of India micro-data.
#              Comtrade sectors for product accounts and Eora sectors for both product and sercive accounts.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())

include("MicroDataReader.jl")
include("ImportTransformer.jl")
include("../Comtrade/HsDataReader.jl")
include("../Comtrade/TradeMatrixBuilder.jl")
include("../Eora/FinalDemandReader.jl")
include("../converting/HsConcMatBuilder.jl")
include("../converting/XLSXextractor.jl")

using .MicroDataReader
using .ImportTransformer
using .HsDataReader
using .TradeMatrixBuilder
using .FinalDemandReader
using .HsConcMatBuilder
using .XLSXextractor

mdr = MicroDataReader
it = ImportTransformer
hdr = HsDataReader
tmb = TradeMatrixBuilder
fdr = FinalDemandReader
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

print(" Expenditure matrix building: ")
expData = mdr.makeExpenditureMatrix()
println("complete")

print(" Expenditure analysis: ")
it.analyzeExpenditures(expData[1], expData[2], expData[3])
println("complete")

mdr.initVars()
#=
# UN Comtrade nation-by-nation trade data reading process

path = "/Users/leejimac/github/microData/Comtrade/data/"

print(" Trade data reading: ")
hdr.readTradeData(path * "hs-2011b.csv", "H3")
println("completed")
print(" Trade matrix builing: ")
tmb.readHsCategories(path*"H3_classification.txt", [6])
tdData = tmb.buildTradeMatrix(hdr.trades, ["India"], "Import")
println("completed")

hdr.initVars()

# Converting process of Comtrade data to India micro-data format

path = "/Users/leejimac/github/microData/converting/data/"
inputHS = path*"India-HS converting table_Ver1.3.xlsx"

print(" Concordance matrix building: ")
hb.readXlsxData(inputHS, "India", "India_HS_Lv6")
conMatNorm = hb.normConMat(hb.buildConMat())
println("complete")

path = Base.source_dir()*"/data/transformed/"
print(" Import data trasform: ")
it.transformHStoIND(tdData[2], tdData[3], tdData[1]["India"], conMatNorm.conMat)
it.printTransfImport(path*"Import_tranformed_fromHS_toIND.txt")
println("complete")
=#
# Eora household's final-demand import sector data reading process

path = "/Users/leejimac/github/microData/Eora/data/"
eoraIndexFile = path * "index.xlsx"
finalDemandFile = path * "2011_eora_y.csv"

print(" EORA Final demand reading: ")
fdr.readIndexData(eoraIndexFile)
fdr.readFinalDemand(finalDemandFile, ["India"])
fdData = fdr.getFinalDemand("India")    # fdData = [fdMat, sec, abb]
println("complete")

# Converting process of Eora final demand data to India micro-data format

path = "/Users/leejimac/github/microData/converting/data/"
concordanceFile = path * "India(STAT) vs EORA_Ver1.2.xlsx"

print(" Concordance matrix building: ")
xls.readXlsxData(concordanceFile, "India")
xls.buildConMat()
cmn = xls.normConMat()   # {a3, {ConMat, sumEora, sumNat}}
println("complete")

path = Base.source_dir()*"/data/transformed/"
print(" Import data trasform: ")
it.transformEORAtoIND(fdData[3], fdData[2], fdData[1], cmn)
it.printTransfImport(path*"Import_tranformed_fromEora_toIND.txt")
println("complete")

# Households' share calculation process

print(" Household share calculation: ")

#it.calculateHHshare(path*"ImportHH_tranformed_fromHS_toIND.txt", false)
it.calculateHHshare(path*"ImportHH_tranformed_fromEora_toIND.txt", false)
#it.calculateHHshare("", false)
println("complete")
