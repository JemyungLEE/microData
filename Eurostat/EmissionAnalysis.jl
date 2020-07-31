# Developed date: 28. Jul. 2020
# Last modified date: 31. Jul. 2020
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
    println("complete")
end
if PPPConv; print(" PPP converting: "); mdr.convertToPPP(pppfile); println("complete") end

# Converting process of Eora final demand data to India micro-data format
concordanceFile = filePath * "index/EU_EORA_Conc_ver1.0.xlsx"
print(" Concordance matrix building:")
print(" xlsx reading"); xls.readXlsxData(concordanceFile, nation)
print(", matrix builing"); xls.buildConMat()
print(", substitutin adding"); xls.addSubstSec(mdr.heSubst, mdr.heRplCd, mdr.heCats)
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

print(" MRIO table reading: ")
path = "../Eora/data/" * string(year) * "/" * string(year)
ee.readIOTables(year, path*"_eora_t.csv", path*"_eora_v.csv", path*"_eora_y.csv", path*"_eora_q.csv")
ee.rearrangeIndex()
ee.rearrangeTables(year)
println("complete")

println("CF calculation: ")
path = Base.source_dir()*"/data/emission/"
st = time()     # check start time
et = time()
i = 1
ns = length(mdr.nations)
for n in mdr.nations
    emissionFile = path * string(year) * "_" * nation * "_hhs_emission.txt"
    print(n, ":")
    print(" domestic data")
    ee.getDomesticData(mdr.expTable[year][n], mdr.hhsList[year][n], vcat(mdr.heCodes, mdr.heSubst))
    print(", weighted concordance matrix")
    ee.buildWeightedConcMat(year, ee.abb[mdr.nationNames[n]], cmn, output = filePath*"index/concordance/Eurostat_Eora_weighted_concordance_table.csv")
    print(", carbon footprint")
    ee.calculateEmission(year, false, 0)
    ee.printEmissions(year, emissionFile)

    (eMin, eSec) = fldmod(floor(Int, time() - et), 60)
    (eHr, eMin) = fldmod(eMin, 60)
    et = time()
    (rMin, rSec) = fldmod(floor(Int, ((time()-st)/i)*(ns-i)), 60)
    (rHr, rMin) = fldmod(rMin, 60)
    println(", ",n," (",i,"/",ns,") ",eHr,":",eMin,":",eSec," elapsed, total ",rHr,":",rMin,":",rSec," remained")
    i += 1
end
println("complete")

println(" ... all complete")
