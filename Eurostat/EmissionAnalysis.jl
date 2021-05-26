# Developed date: 28. Jul. 2020
# Last modified date: 26. May. 2021
# Subject: Estimate carbon footprint by final demands of Eora
# Description: Calculate carbon emissions by utilizing Eora T, V, Y, and Q tables.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("MicroDataReader.jl")
include("EmissionEstimator.jl")
include("../converting/XLSXextractor.jl")

using .MicroDataReader
using .XLSXextractor
using .EmissionEstimator

mdr = MicroDataReader
xls = XLSXextractor
ee = EmissionEstimator

filePath = Base.source_dir() * "/data/"
categoryFile = filePath * "index/Eurostat_Index_ver3.2.xlsx"
microDataPath = filePath * "microdata/"

CSV_reading = true     # reading micro-data from extracted CSV files
XLSX_reading = false     # reading micro-data from original XLSX files

CurrencyConv = true; erfile = filePath * "index/EUR_USD_ExchangeRates.txt"
PPPConv = false; pppfile = filePath * "index/PPP_ConvertingRates.txt"

DE_conv = filePath * "index/EmissionCovertingRate.txt"

IE_mode = true      # indirect carbon emission estimation
DE_mode = true      # direct carbon emission estimation

nation = "Eurostat"
year = 2010
catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]

# Qtable = "I_CHG_CO2"
Qtable = "PRIMAP"

abrExpMode = false
substMode = true
eoraRevised = true
scaleMode = true

if substMode; substTag = "_subst" else substTag = "" end
if scaleMode; scaleTag = "Scaled_" else scaleTag = "" end

hhsfile = filePath * "extracted/Households.csv"
mmsfile = filePath * "extracted/Members.csv"
expfile = filePath * "extracted/"*scaleTag*"Expenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
ctgfile = filePath * "extracted/Category_" * depthTag[catDepth] * ".csv"
sttfile = filePath * "extracted/MicroData_Statistics_"*scaleTag*depthTag[catDepth]*substTag*".csv"
sbstfile = filePath * "extracted/SubstituteCodes_" * depthTag[catDepth] * ".csv"

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(categoryFile, depth=catDepth, catFile=ctgfile, inclAbr=abrExpMode)
println("completed")

print(" Micro-data reading: ")
if CSV_reading
    print("CSV")
    if substMode; print(" substitutes"); mdr.readSubstCodesCSV(sbstfile) end
    print(", households"); mdr.readPrintedHouseholdData(hhsfile)
    print(", members"); mdr.readPrintedMemberData(mmsfile)
    print(", expenditures(", expfile, ")"); mdr.readPrintedExpenditureData(expfile, substitute=substMode, buildHhsExp=true)
elseif XLSX_reading
    println("XLSX")
    println(", households"); mdr.readHouseholdData(year, microDataPath, visible=true, substitute=substMode)
    println(", members"); mdr.readMemberData(year, microDataPath, visible=true)
    print(", expenditures"); mdr.buildExpenditureMatrix(year, expfile, substitute=substMode)
    print(", statistics");mdr.makeStatistics(year, sttfile, substitute=substMode)
else println()
end
print(", sector data"); ee.getSectorData(year, vcat(mdr.heCodes, mdr.heSubst))
println(" ... completed")

if CurrencyConv; print(" Currency exchanging:")
    print(" USD transform"); mdr.exchangeExpCurrency(erfile)
    print(" rebuild expenditure matrix"); mdr.buildExpenditureMatrix(year, replace(expfile, ".csv"=>"_USD.csv"), substitute=substMode)
    println(" complete")
end
if PPPConv; print(" PPP converting: "); mdr.convertToPPP(pppfile); println("complete") end

if IE_mode
    # Converting process of Eora final demand data to India micro-data format
    concordanceFile = filePath * "index/EU_EORA_Conc_ver1.2.xlsx"
    print(" Concordance matrix building:")
    print(" xlsx reading"); xls.readXlsxData(concordanceFile, nation)
    print(", matrix builing"); xls.buildConMat()
    print(", substitution appending"); xls.addSubstSec(mdr.heSubst, mdr.heRplCd, mdr.heCats)
    print(", normalization"); cmn = xls.normConMat()   # {a3, conMat}
    conMatFile = filePath * "index/concordance/ConcMat.txt"
    conSumMatFile = filePath * "index/concordance/ConcSumMat.txt"
    xls.printConMat(conMatFile, nation, norm = true, categ = true)
    xls.printSumNat(conSumMatFile, nation, norm = true)
    println(" complete")

    # Eora household's final-demand import sector data reading process
    print(" Eora index reading: ")
    if eoraRevised; ee.readIndexXlsx("../Eora/data/index/revised/", revised=true)
    else ee.readIndexXlsx("../Eora/data/index/Eora_index.xlsx")
    end
    println("complete")

    print(" MRIO table reading:")
    path = "../Eora/data/" * string(year) * "/" * string(year)
    print(" IO table"); ee.readIOTables(year, path*"_eora_t.csv", path*"_eora_v.csv", path*"_eora_y.csv", path*"_eora_q.csv")
    print(", rearrange"); ee.rearrangeIndex(qmode=Qtable); ee.rearrangeTables(year, qmode=Qtable)
    print(", Leontief matrix"); ee.calculateLeontief(year)
    println(" complete")
end
if DE_mode
    print(" Direct emission converting indices reading:")
    ee.readEmissionRates(year, categoryFile, DE_conv)
    println(" complete")
end

println(" Emission calculation: ")
path = Base.source_dir()*"/data/emission/"
ns = length(mdr.nations)

# read consuming time
tf = open(filePath * "index/"*Qtable*"_time.txt")
tp = Dict{String, Float64}()
readline(tf)
for l in eachline(tf); s = split(l, '\t'); tp[s[1]] = parse(Float64, s[2]) end
close(tf)

ttp = 0
st = time()    # check start time
for i = 1:ns
    n = mdr.nations[i]
    emissionFile = path * string(year) * "_" * n * "_hhs_"*scaleTag*"IE_"*Qtable*".txt"
    deFile = path * string(year) * "_" * n * "_hhs_"*scaleTag*"DE.txt"
    print("\t", n, ":")
    print(" data"); ee.getDomesticData(year, mdr.expTable[year][n], mdr.hhsList[year][n])
    if DE_mode; print(", DE")
        ee.calculateDirectEmission(year, n)
        ee.printDirectEmissions(year, deFile)
    end
    if IE_mode
        print(", concordance"); ee.buildWeightedConcMat(year, ee.abb[mdr.nationNames[n]], cmn, output = filePath*"index/concordance/"*scaleTag*"Eurostat_Eora_weighted_concordance_table.csv")
        print(", IE"); ee.calculateIndirectEmission(year, false, 0, reuseLti=true)
        ee.printIndirectEmissions(year, emissionFile)
    end

    global ttp += tp[n]
    elap = floor(Int, time() - st)
    (eMin, eSec) = fldmod(elap, 60)
    (eHr, eMin) = fldmod(eMin, 60)
    (rMin, rSec) = fldmod(floor(Int, elap/ttp*(1-ttp)), 60)
    (rHr, rMin) = fldmod(rMin, 60)
    println(", ",n," (",i,"/",ns,") ",eHr,":",eMin,":",eSec," elapsed, total ",rHr,":",rMin,":",rSec," remained")
end
println("complete")

println(" ... all complete")
