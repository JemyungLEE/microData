# Developed date: 21. Oct. 2019
# Last modified date: 29. Oct. 2019
# Subject: Harmonized System (HS) UN comtrade data analyzer
# Description: analysis UN comtrade trading data
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())

using DataFrames
include("HsDataReader.jl")
using .HsDataReader
hdr = HsDataReader

inputData = "hs-2011b.csv"

inputData = Base.source_dir() * "/data/" * inputData

hdr.readTradeData(inputData)
#hdr.readTradeDataCSV(inputData)
df = hdr.exportDataFrames()

println(describe(df))
