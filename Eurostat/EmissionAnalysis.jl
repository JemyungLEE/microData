# Developed date: 28. Jul. 2020
# Last modified date: 4. Aug. 2020
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
categoryFile = filePath * "index/Eurostat_Index_ver0.7.xlsx"
microDataPath = filePath * "microdata/"

CurrencyConv = true; erfile = filePath * "index/EUR_USD_ExchangeRates.txt"
PPPConv = true; pppfile = filePath * "index/PPP_ConvertingRates.txt"

nation = "Eurostat"
year = 2010
catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]

abrExpMode = false
substMode = true
hhsExpMode = true

hhsfile = filePath * "extracted/Households.csv"
mmsfile = filePath * "extracted/Members.csv"
expfile = filePath * "extracted/Expenditure_matrix_"*depthTag[catDepth]*".csv"
ctgfile = filePath * "extracted/Category_" * depthTag[catDepth] * ".csv"
sbstfile = filePath * "extracted/SubstituteCodes_" * depthTag[catDepth] * ".csv"

println("[Process]")
print(" Category codes reading: ")
mdr.readCategory(categoryFile, depth=catDepth, catFile=ctgfile, inclAbr=abrExpMode)
println("completed")

print(" Micro-data CSV reading:")
if substMode; print(" substitutes"); mdr.readSubstCodesCSV(sbstfile) end
print(", households"); mdr.readPrintedHouseholdData(hhsfile)
print(", members"); mdr.readPrintedMemberData(mmsfile)
print(", expenditures"); mdr.readPrintedExpenditureData(expfile, substitute=substMode, buildHhsExp=true)
println(" completed")

if CurrencyConv; print(" Currency exchanging: ")
    print(" USD transform"); mdr.exchangeExpCurrency(erfile)
    print(" rebuild expenditure matrix"); mdr.buildExpenditureMatrix(year, replace(expfile, ".csv"=>"_USD.csv"), substitute=substMode)
    println(" complete")
end
if PPPConv; print(" PPP converting: "); mdr.convertToPPP(pppfile); println("complete") end

# Converting process of Eora final demand data to India micro-data format
concordanceFile = filePath * "index/EU_EORA_Conc_ver1.0.xlsx"
print(" Concordance matrix building:")
print(" xlsx reading"); xls.readXlsxData(concordanceFile, nation)
print(", matrix builing"); xls.buildConMat()
print(", substitution appending"); xls.addSubstSec(mdr.heSubst, mdr.heRplCd, mdr.heCats)
print(", normalization"); cmn = xls.normConMat()   # {a3, conMat}
# conMatFile = filePath * "index/concordance/ConcMat.txt"
# conSumMatFile = filePath * "index/concordance/ConcSumMat.txt"
# xls.printConMat(conMatFile, nation, norm = true, categ = true)
# xls.printSumNat(conSumMatFile, nation, norm = true)
println(" complete")

# Eora household's final-demand import sector data reading process
eoraIndexFile = "../Eora/data/index/Eora_index.xlsx"
print(" Eora index reading: ")
ee.readIndexXlsx(eoraIndexFile)
println("complete")

print(" MRIO table reading:")
path = "../Eora/data/" * string(year) * "/" * string(year)
print("IO table"); ee.readIOTables(year, path*"_eora_t.csv", path*"_eora_v.csv", path*"_eora_y.csv", path*"_eora_q.csv")
print(", rearrange"); ee.rearrangeIndex()
print(""); ee.rearrangeTables(year)
print(", Leontief matrix"); ee.calculateLeontief(year)
println(" complete")

println(" CF calculation: ")
path = Base.source_dir()*"/data/emission/"
ns = length(mdr.nations)
nhhs = [length(mdr.hhsList[year][x]) for x in mdr.nations]
nrh = sum(nhhs); nch = 0
st = time()    # check start time
for i = 1:ns
    n = mdr.nations[i]
    emissionFile = path * string(year) * "_" * n * "_hhs_emission.txt"
    print("\t", n, ":")
    print(" data"); ee.getDomesticData(mdr.expTable[year][n], mdr.hhsList[year][n], vcat(mdr.heCodes, mdr.heSubst))
    print(", concordance"); ee.buildWeightedConcMat(year, ee.abb[mdr.nationNames[n]], cmn, output = filePath*"index/concordance/Eurostat_Eora_weighted_concordance_table.csv")
    print(", CF"); ee.calculateEmission(year, false, 0, reuseLti=true)
    ee.printEmissions(year, emissionFile)

    global nch += nhhs[i]; global nrh -= nhhs[i]
    elap = floor(Int, time() - st)
    (eMin, eSec) = fldmod(elap, 60)
    (eHr, eMin) = fldmod(eMin, 60)
    (rMin, rSec) = fldmod(floor(Int, elap/nch*nrh), 60)
    (rHr, rMin) = fldmod(rMin, 60)
    println(", ",n," (",i,"/",ns,") ",eHr,":",eMin,":",eSec," elapsed, total ",rHr,":",rMin,":",rSec," remained")
end
println("complete")

println(" ... all complete")
