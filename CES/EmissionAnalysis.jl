# Developed date: 13. Apr. 2021
# Last modified date: 5. May. 2021
# Subject: Estimate carbon footprint by household consumptions
# Description: Calculate direct and indirect carbon emissions
#              by linking household consumptions and global supply chain,
#              utilizing Eora T, V, Y, and Q tables.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("MicroDataReader.jl")
include("ConcMatBuilder.jl")
include("EmissionEstimator.jl")

using .MicroDataReader
using .ConcMatBuilder
using .EmissionEstimator

mdr = MicroDataReader
cmb = ConcMatBuilder
ee = EmissionEstimator

year = 2018
eoraYear = 2015     # eoraYear = year
nation = "Indonesia"
natA3 = "IDN"

filePath = Base.source_dir() * "/data/" * natA3 * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissionPath = filePath * "emission/"

curConv = true; erfile = indexFilePath * "CurrencyExchangeRates.txt"
pppConv = false; pppfile = filePath * "PPP_ConvertingRates.txt"

DE_conv = filePath * "index/EmissionCovertingRate.txt"

IE_mode = false      # indirect carbon emission estimation
DE_mode = false      # direct carbon emission estimation

# Qtable = "I_CHG_CO2"
Qtable = "PRIMAP"
scaleMode = true

if scaleMode; scaleTag = "_Scaled" else scaleTag = "" end

regInfoFile = extractedPath * nation * "_" * string(year) * "_RegionInfo.txt"
cmmfile = extractedPath * nation * "_" * string(year) * "_Commodities.txt"
hhsfile = extractedPath * nation * "_" * string(year) * "_Households.txt"
mmsfile = extractedPath * nation * "_" * string(year) * "_Members.txt"
itemfile = indexFilePath * nation * "_" * string(year) * "_Commodity_items.txt"
exdfile = extractedPath * nation * "_" * string(year) * "_Expenditure.txt"
exmfile = extractedPath * nation * "_" * string(year) * scaleTag * "_Expenditure_matrix.txt"

println("[Process]")
# print(" Category codes reading: ")
# mdr.readCategory(categoryFile, depth=catDepth, catFile=ctgfile, inclAbr=abrExpMode)
# println("completed")

print(" Micro-data reading: ")
print(", regions"); mdr.readPrintedRegionData(year, nation, regInfoFile)
print(", households"); mdr.readPrintedHouseholdData(year, nation, hhsfile)
# print(", members"); mdr.readPrintedMemberData(year, nation, mmsfile)
print(", sectors"); mdr.readPrintedSectorData(year, nation, cmmfile)
print(", expenditures(", expfile, ")"); mdr.readPrintedExpenditureData(readPrintedHouseholdDataexpfile)
if curConv; print(", exchange"); mdr.exchangeExpCurrency(year, nation, erfile, inverse=true) end
if PPPConv; print(", ppp"); mdr.convertToPPP(pppfile); println("complete") end
print(", build matrix"); mdr.buildExpenditureMatrix(year, nation, expfile)
# print(", matrix"); mdr.readPrintedExpenditureMatrix(year, nation, exmfile)
println(" ... completed")

if IE_mode
    # Converting process of Eora final demand data to India micro-data format
    concordanceFile = filePath * "/" * natA3 * "/index/"* natA3 *"_EORA_Conc_ver0.9.xlsx"
    conMatFile = filePath * "/" * natA3 * "/index/concordance/ConcMat.txt"
    conSumMatFile = filePath * "/" * natA3 * "/index/concordance/ConcSumMat.txt"
    print(" Concordance matrix building:")
    # print(" xlsx reading"); cmb.readXlsxData(concordanceFile, nation, weight=false)
    nationFile = filePath * "/" * natA3 * "/index/concordance/Eora_nations.txt"
    sectorFile = filePath * "/" * natA3 * "/index/concordance/Eora_sectors.txt"
    concFile = filePath * "/" * natA3 * "/index/concordance/LinkedSectors_IE.txt"
    print(" data reading"); cmb.readConcMatFile(nationFile, sectorFile, concFile, nation, weight=false)
    # linkSecFile = filePath * "/" * natA3 * "/index/concordance/LinkedSectors_IE.txt"
    # print(", linkage printing"); cmb.exportLinkedSectors(linkSecFile, natA3, mrio="Eora")
    print(", matrix builing"); cmb.buildConMat()
    print(", normalization"); cmn = cmb.normConMat()   # {a3, conMat}

    cmb.printConMat(conMatFile, nation, norm = true, categ = true)
    # cmb.printSumNat(conSumMatFile, nation, norm = true)
    println(" complete")

    # Eora household's final-demand import sector data reading process
    print(" MRIO table reading:")
    path = "../Eora/data/" * string(eoraYear) * "/" * string(eoraYear)
    print(" index"); ee.readIndex("../Eora/data/index/revised/")
    print(", table"); ee.readIOTables(eoraYear, path*"_eora_t.csv", path*"_eora_v.csv", path*"_eora_y.csv", path*"_eora_q.csv")
    print(", rearrange"); ee.rearrangeIndex(qmode=Qtable); ee.rearrangeTables(year, qmode=Qtable)
    print(", Leontief matrix"); ee.calculateLeontief(eoraYear)
    println(" complete")
end
# if DE_mode
#     print(" Direct emission converting indices reading:")
#     ee.readEmissionRates(year, categoryFile, DE_conv)
#     println(" complete")
# end

println(" Emission calculation: ")
ieFile = emissionPath * string(year) * "_" * natA3 * "_hhs_"*scaleTag*"IE_"*Qtable*".txt"
deFile = emissionPath * string(year) * "_" * natA3 * "_hhs_"*scaleTag*"DE.txt"
print(" data"); ee.getDomesticData(year, mdr.expTable[year], mdr.hhsList[year])
# if DE_mode; print(", DE")
#     ee.calculateDirectEmission(year, n)
#     ee.printDirectEmissions(year, deFile)
# end
if IE_mode
    print(", concordance"); ee.buildWeightedConcMat(year, ee.abb[mdr.nationNames[n]], cmn, output = filePath*"index/concordance/"*scaleTag*"Eurostat_Eora_weighted_concordance_table.csv")
    print(", IE"); ee.calculateIndirectEmission(year)
    ee.printIndirectEmissions(year, emissionFile)
end

println("complete")

println(" ... all complete")
