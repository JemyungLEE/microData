# Developed date: 31. Mar. 2021
# Last modified date: 31. Mar. 2021
# Subject: Household consumption expenditure survey microdata analysis
# Description: proceed microdata analysis process
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
include("MicroDataReader.jl")
using .MicroDataReader
mdr = MicroDataReader

year = 2018
nation = "IDN"

filePath = Base.source_dir() * "/data/" * nation * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath =filePath * "extracted/"

curConv = true; erfile = indexFilePath * "CurrencyExchangeRates.txt"
pppConv = false; pppfile = filePath * "PPP_ConvertingRates.txt"

gapMitigation = false    # filling gaps between national account and HBS expenditures

printData = true

# input files
popfile = indexFilePath * nation * "_" * string(year) * "_Population.txt"
regfile = indexFilePath * nation * "_" * string(year) * "_Regions.txt"
hidxfile = indexFilePath * nation * "_" * string(year) * "_Households_microdata_attribute.txt"
eidxfile = indexFilePath * nation * "_" * string(year) * "_Expenditures_microdata_attribute.txt"
itemfile = indexFilePath * nation * "_" * string(year) * "_Commodity_items.txt"
catfile = indexFilePath * nation * "_" * string(year) * "_Commodity_category.txt"

# output files
hhsfile = extractedPath * nation * "_" * string(year) * "_Households.txt"
mmsfile = extractedPath * nation * "_" * string(year) * "_Members.txt"
exdfile = extractedPath * nation * "_" * string(year) * "_Expenditure.txt"
wghfile = extractedPath * nation * "_" * string(year) * "_Weight.txt"

exmfile = extractedPath * nation * "_" * string(year) * "_Expenditure_matrix.txt"
sttfile = extractedPath * nation * "_" * string(year) * "_MicroData_statistics.txt"

scexpfile = extractedPath * nation * "_" * string(year) * "_Scaled_Expenditure_matrix.txt"

println("[Process]")
print(" Region information reading: ")
print("population"); mdr.readPopulation(year, nation, popfile)
print(", region"); mdr.readRegion(year, nation, regfile)
println(" ... completed")

print(" Micro-data reading: ")


print("microdata"); mdr.readMicroData(year, nation, microDataPath, hidxfile, "", itemfile, eidxfile)
if curConv; print(", exchange"); mdr.exchangeExpCurrency(year, nation, erfile, inverse=true) end
print(", matrix"); mdr.buildExpenditureMatrix(year, nation, exmfile, print_err=true)
print(", weight"); mdr.calculatePopWeight(year, nation, wghfile, district=true, province=false)
println(" ... completed")

if gapMitigation; print(" GDP-Survey gap mitigating: ")

    println("completed")
end

if printData; print(" Extracted data printing:")
    mdr.printHouseholdData(year, nation, hhsfile, prov_wgh=false, dist_wgh=true, ur_dist=false)
    mdr.printExpenditureData(year, nation, exdfile)
    println(" ... completed")
end

println("[done]")
