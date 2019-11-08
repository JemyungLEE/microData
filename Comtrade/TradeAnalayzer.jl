# Developed date: 21. Oct. 2019
# Last modified date: 8. Nov. 2019
# Subject: Harmonized System (HS) UN comtrade data analyzer
# Description: analysis UN comtrade trading data
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())

using DataFrames
include("HsDataReader.jl")
include("TradeMatrixBuilder.jl")
using .HsDataReader
using .TradeMatrixBuilder
hdr = HsDataReader
tmb = TradeMatrixBuilder

path = Base.source_dir() * "/data/"
inputData = path * "hs-2011b.csv"
outputFile = path * "HS_dataframe.txt"

hsCodeLevel = "H3"
nation = ["India"]

print("Trade data reading: ")
hdr.readTradeData(inputData, hsCodeLevel)
println("completed")
#hdr.analyzeStatus()
#hdr.analyzeStatus(nation)
#hdr.analyzeStatus([], nation)
#hdr.analyzeStatus(nation, nation, "Import")
#hdr.readTradeDataCSV(inputData)
#hdr.exportDataFrames()
#hdr.printDataFrames(outputFile)

#=
matchResultFile = path*"HS_"*nation[1]*"_matching.txt"
print("Trade data matching test: ")
hdr.matchingTest(nation[1], matchResultFile)
println("completed")
=#
flow = "Import"
#flow = "Export"
#flow = "Net"
tradeMatrixFile = path * "TradeMatrix_$flow.txt"

print("Trade matrix builing: ")
tmb.readHsCategories(path*"H3_classification.txt", [6])
tmb.buildTradeMatrix(hdr.trades, nation, flow)
tmb.printTradeMatrix(tradeMatrixFile)
println("completed")
