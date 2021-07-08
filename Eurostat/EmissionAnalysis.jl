# Developed date: 28. Jul. 2020
# Last modified date: 7. Jul. 2021
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
indexPath = filePath * "index/"
extrPath = filePath * "extracted/"
concPath = indexPath * "concordance/"
categoryFile = indexPath * "Eurostat_Index_ver4.1.xlsx"
concFiles = Dict(2010 => indexPath*"2010_EU_EORA_Conc_ver1.5.xlsx", 2015 => indexPath*"2015_EU_EORA_Conc_ver1.1.xlsx")
natLabels = Dict(2010 => "Eurostat", 2015 => "Eurostat_2015")
conMatFile = concPath * "ConcMat.txt"
conSumMatFile = concPath * "ConcSumMat.txt"

CSV_reading = true     # reading micro-data from extracted CSV files
XLSX_reading = false     # reading micro-data from original XLSX files

CurrencyConv = true; erfile = indexPath * "EUR_USD_ExchangeRates.txt"
PPPConv = false; pppfile = indexPath * "PPP_ConvertingRates.txt"

DE_conv = indexPath * "EmissionCovertingRate.txt"

IE_mode = true      # indirect carbon emission estimation
DE_mode = false      # direct carbon emission estimation

nation = "Eurostat"
year = 2015
catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]

microDataPath = filePath * "microdata/" * string(year) * "/"

# Qtable = "I_CHG_CO2"
Qtable = "PRIMAP"

abrExpMode = false
substMode = true
eoraRevised = true
scaleMode = true
testMode = true; test_nats = ["HU", "SE"]

if substMode; substTag = "_subst" else substTag = "" end
if scaleMode; scaleTag = "Scaled_" else scaleTag = "" end

hhsfile = extrPath * string(year) * "_Households.csv"
mmsfile = extrPath * string(year) * "_Members.csv"
expfile = extrPath * string(year) * "_" * scaleTag*"Expenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
ctgfile = extrPath * string(year) * "_Category_" * depthTag[catDepth] * ".csv"
sttfile = extrPath * string(year) * "_MicroData_Statistics_"*scaleTag*depthTag[catDepth]*substTag*".csv"
sbstfile = extrPath * string(year) * "_SubstituteCodes_" * depthTag[catDepth] * ".csv"

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(year, categoryFile, depth=catDepth, catFile=ctgfile, coicop=scaleMode, inclAbr=abrExpMode)
println("completed")

print(" Micro-data reading: ")
if CSV_reading
    print("CSV")
    print(", households"); mdr.readPrintedHouseholdData(hhsfile)
    print(", members"); mdr.readPrintedMemberData(mmsfile)
    if substMode; print(" substitutes"); mdr.readSubstCodesCSV(sbstfile) end
    print(", expenditures(", expfile, ")"); mdr.readPrintedExpenditureData(expfile, substitute=substMode, buildHhsExp=true)
elseif XLSX_reading
    println("XLSX")
    println(", households"); mdr.readHouseholdData(year, microDataPath, visible=true, substitute=substMode)
    println(", members"); mdr.readMemberData(year, microDataPath, visible=true)
    print(", expenditures"); mdr.buildExpenditureMatrix(year, substitute=substMode)
    print(", statistics");mdr.makeStatistics(year, sttfile, substitute=substMode)
else println()
end
print(", sector data"); ee.getSectorData(year, mdr.heCodes, mdr.heSubst)
println(" ... completed")

if CurrencyConv; print(" Currency exchanging:")
    print(" USD transform"); mdr.exchangeExpCurrency(erfile)
    print(" rebuild expenditure matrix"); mdr.buildExpenditureMatrix(year, substitute=substMode)
    println(" complete")
end
if PPPConv; print(" PPP converting: "); mdr.convertToPPP(pppfile); println("complete") end

if IE_mode
    # Converting process of Eora final demand data to India micro-data format
    print(" Concordance matrix building:")
    print(" concordance"); xls.readXlsxData(concFiles[year], nation, nat_label = natLabels[year])
    print(", matrix"); xls.buildConMat()
    print(", substitution"); xls.addSubstSec(year, mdr.heSubst, mdr.heRplCd, mdr.heCats)
    print(", normalization"); cmn = xls.normConMat()   # {a3, conMat}
    print(", printing"); xls.printConMat(conMatFile, nation, norm = true, categ = true)
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
if !testMode; nat_list = mdr.nations else nat_list = test_nats end
ns = length(nat_list)

# read consuming time
tf = open(indexPath * Qtable * "_time.txt")
tp = Dict{String, Float64}()
readline(tf)
for l in eachline(tf); s = split(l, '\t'); tp[s[1]] = parse(Float64, s[2]) end
close(tf)

ttp = 0
st = time()    # check start time
for i = 1:ns
    n = nat_list[i]
    emissionFile = path * string(year) * "_" * n * "_hhs_"*scaleTag*"IE_"*Qtable*".txt"
    deFile = path * string(year) * "_" * n * "_hhs_"*scaleTag*"DE.txt"
    print("\t", n, ":")
    print(" data"); ee.getDomesticData(year, n, mdr.expTable, mdr.hhsList)
    if DE_mode; print(", DE")
        ee.calculateDirectEmission(year, n)
        ee.printDirectEmissions(year, deFile)
    end
    if IE_mode
        wgh_conc_file = concPath * scaleTag * "Eurostat_Eora_weighted_concordance_table.csv"
        print(", concordance"); ee.buildWeightedConcMat(year, ee.abb[mdr.nationNames[n]], cmn, output = wgh_conc_file)
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
