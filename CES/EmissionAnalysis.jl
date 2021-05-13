# Developed date: 13. Apr. 2021
# Last modified date: 7. May. 2021
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

cesYear = 2018; exchYear = cesYear
eoraYear = 2015
nation = "Indonesia"
natA3 = "IDN"
natCurr = "IDR"
fitEoraYear = true      # scaling micro-data's expenditure to fit the Eora target year
readMembers = false     # read member data

filePath = Base.source_dir() * "/data/" * natA3 * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissionPath = filePath * "emission/"
concordancePath = filePath * "index/concordance/"

curConv = true; curr_target = "USD"; erfile = indexFilePath * "CurrencyExchangeRates.txt"
pppConv = false; pppfile = filePath * "PPP_ConvertingRates.txt"

DE_conv = filePath * "index/EmissionCovertingRate.txt"

IE_mode = true      # indirect carbon emission estimation
DE_mode = false      # direct carbon emission estimation

# Qtable = "I_CHG_CO2"
Qtable = "PRIMAP"
scaleMode = false
sparseMode = true
quantMode = true

if scaleMode; scaleTag = "_Scaled" else scaleTag = "" end

regInfoFile = extractedPath * natA3 * "_" * string(cesYear) * "_RegionInfo.txt"
cmmfile = extractedPath * natA3 * "_" * string(cesYear) * "_Commodities.txt"
hhsfile = extractedPath * natA3 * "_" * string(cesYear) * "_Households.txt"
mmsfile = extractedPath * natA3 * "_" * string(cesYear) * "_Members.txt"
itemfile = indexFilePath * natA3 * "_" * string(cesYear) * "_Commodity_items.txt"
expfile = extractedPath * natA3 * "_" * string(cesYear) * "_Expenditure_"*natCurr*".txt"
exmfile = extractedPath * natA3 * "_" * string(cesYear) * scaleTag * "_Expenditure_matrix_"*natCurr*".txt"

println("[Process]")

print(" Micro-data reading: ")
print(", regions"); mdr.readPrintedRegionData(cesYear, natA3, regInfoFile)
print(", households"); mdr.readPrintedHouseholdData(cesYear, natA3, hhsfile)
if readMembers; print(", members"); mdr.readPrintedMemberData(cesYear, natA3, mmsfile) end
print(", sectors"); mdr.readPrintedSectorData(cesYear, natA3, cmmfile)
print(", expenditures(", rsplit(expfile, '/', limit=2)[2], ")"); mdr.readPrintedExpenditureData(cesYear, natA3, expfile, quantity=quantMode)
if fitEoraYear && eoraYear != nothing && eoraYear != cesYear; print(" Expenditure scaling: from $cesYear to $eoraYear")
    exchYear = eoraYear
    cpiSecFile = indexFilePath * "CPI/CPI_"*natA3*"_sectors.txt"
    statFile = indexFilePath * "CPI/CPI_"*natA3*"_values.txt"
    linkFile = indexFilePath * "CPI/CPI_"*natA3*"_link.txt"
    print(", scaling"); mdr.scalingExpByCPI(cesYear, natA3, cpiSecFile, statFile, linkFile, eoraYear, period="year", region="district", revHH=true, revMat=false)
    println(" ... completed")
end
if curConv; print(", exchange"); mdr.exchangeExpCurrency(cesYear, exchYear, natA3, natCurr, erfile, target_curr=curr_target) end
if pppConv; print(", ppp"); mdr.convertToPPP(eoraYear, natA3, pppfile); println("complete") end
print(", matrix"); mes = mdr.buildExpenditureMatrix(cesYear, natA3, period = 365, quantity = quantMode)
# print(", matrix"); mdr.readPrintedExpenditureMatrix(cesYear, natA3, exmfile)
println(" ... completed")

if IE_mode
    conMatFile = concordancePath * "ConcMat.txt"
    conSumMatFile = concordancePath * "ConcSumMat.txt"
    nationFile = concordancePath * "Eora_nations.txt"
    sectorFile = concordancePath * "Eora_sectors.txt"
    concFile = concordancePath * "LinkedSectors_IE.txt"

    print(" Concordance matrix building:")
    print(" data reading"); cmb.readIeConcMatFile(nationFile, sectorFile, concFile, weight=false)
    print(", matrix builing"); cmb.buildIeConMat()
    print(", normalization"); cmn = cmb.normConMat()   # {a3, conMat}
    print(", print matrix"); cmb.printConMat(conMatFile, natA3, norm = true, categ = true)
    # cmb.printSumNat(conSumMatFile, natA3, norm = true)
    println(" complete")

    print(" MRIO table reading:")
    path = "../Eora/data/" * string(eoraYear) * "/" * string(eoraYear)
    print(" index"); ee.readIndex("../Eora/data/index/revised/")
    print(", table"); ee.readIOTables(eoraYear, path*"_eora_t.csv", path*"_eora_v.csv", path*"_eora_y.csv", path*"_eora_q.csv")
    print(", rearrange"); ee.rearrangeIndex(qmode=Qtable); ee.rearrangeTables(eoraYear, qmode=Qtable)
    print(", Leontief matrix"); ee.calculateLeontief(eoraYear)
    println(" complete")
end
if DE_mode
    print(" Direct emission preparing:")
    print(" data reading"); ee.readEmissionIntensity(cesYear, natA3, sectorFile, intensityFile)
    print(", concordance matrix building"); cmb.buildDeConcMat(natA3, deCodeFile, concFile, norm = true, output = "")
    println(" complete")
end

println(" Emission calculation: ")
if DE_mode || IE_mode
    ieFile = emissionPath * string(cesYear) * "_" * natA3 * "_hhs_"*scaleTag*"IE_"*Qtable*".txt"
    deFile = emissionPath * string(cesYear) * "_" * natA3 * "_hhs_"*scaleTag*"DE.txt"
    print(" data"); ee.getDomesticData(cesYear, natA3, mdr.hh_list[cesYear][natA3], mdr.sc_list[cesYear][natA3], mdr.expMatrix[cesYear][natA3])
end
if DE_mode; print(", DE")
    ee.calculateDirectEmission(cesYear, natA3, commodity = quantMode, sparseMat = sparseMode)
    ee.printDirectEmissions(cesYear, natA3, outputFile)
end
if IE_mode
    conWghMatFile = filePath*"index/concordance/"*scaleTag*"Eurostat_Eora_weighted_concordance_table.csv"
    print(", concordance"); ee.buildWeightedConcMat(cesYear, eoraYear, natA3, cmn, output = conWghMatFile)
    print(", estimate_IE"); ee.calculateIndirectEmission(cesYear, eoraYear, natA3, sparseMat = sparseMode, elapChk = 1)
    print(", print_IE"); ee.printIndirectEmissions(cesYear, natA3, ieFile)
end
println(" ... complete")

println("[all complete]")
