# Developed date: 28. Jul. 2020
# Last modified date: 30. Aug. 2021
# Subject: Estimate carbon footprint by final demands of Eora
# Description: Calculate carbon emissions by utilizing Eora T, V, Y, and Q tables.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("MicroDataReader.jl")
include("EmissionEstimator.jl")
include("../converting/XLSXextractor.jl")

using Statistics
using .MicroDataReader
using .XLSXextractor
using .EmissionEstimator

mdr = MicroDataReader
xls = XLSXextractor
ee = EmissionEstimator

CSV_reading = true     # reading micro-data from extracted CSV files
XLSX_reading = false     # reading micro-data from original XLSX files

# DE_conv = indexPath * "EmissionCovertingRate.txt"

IE_mode = false             # indirect carbon emission estimation
DE_mode = true              # direct carbon emission estimation
DE_factor_estimate = true   # [true] estimate DE factors from IEA datasets, [false] read DE factors

nation = "Eurostat"
year = 2015
catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]
emission_unit = "tCO2"
currency_unit = "EUR"
reg_group = "EU"

filePath = Base.source_dir() * "/data/"
indexPath = filePath * "index/"
extrPath = filePath * "extracted/"
concPath = indexPath * "concordance/"
emit_path = indexPath * "de/"
microDataPath = filePath * "microdata/" * string(year) * "/"

categoryFile = indexPath * "Eurostat_Index_ver4.2.xlsx"
CurrencyConv = false; erfile = indexPath * "EUR_USD_ExchangeRates.txt"
PPPConv = false; pppfile = indexPath * "PPP_ConvertingRates.txt"

concFiles = Dict(2010 => indexPath*"2010_EU_EORA_Conc_ver1.5.xlsx", 2015 => indexPath*"2015_EU_EORA_Conc_ver1.1.xlsx")
natLabels = Dict(2010 => "Eurostat", 2015 => "Eurostat_2015")
conMatFile = concPath * "ConcMat.txt"
conSumMatFile = concPath * "ConcSumMat.txt"
deSectorFile = emit_path * "DE_sectors.txt"
deHbsLinkFile = emit_path * "DE_linked_sectors_" * string(year) * ".txt"
deIntensityFile = emit_path * "Emission_converting_rate_"*string(year)*"_EU.txt"

# Qtable = "I_CHG_CO2"
Qtable = "PRIMAP"

abrExpMode = false
substMode = true
eoraRevised = true
scaleMode = true
testMode = false; test_nats = ["HU", "SE"]

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
    print(", assembe Conc_mat"); ee.assembleConcMat(year, cmn)
    println(" ... complete")
end
if DE_mode
    print(" Direct emission converting indices reading:")
    # ee.readEmissionRates(year, categoryFile, DE_conv)
    ee.readEmissionData(year, mdr.nationNames, emit_path, output_path = emit_path * "EU/", output_tag = reg_group, integrate = true)
    if DE_factor_estimate
        price_file = emit_path * "EU/" * "Price_" * string(year) * "_EU_" * currency_unit * ".txt"
        emi_intens_file = emit_path * "EU/" * "Emission_intensity_" * string(year) * "_EU_tCO2_per_" * currency_unit* ".txt"
        de_conv_file = emit_path * "Emission_converting_rate_" * string(year) * "_EU.txt"
        ee.exchangeEmCurrency(year, erfile, target = currency_unit, origin = "USD", output = price_file)
        ee.calculateEmissionRates(year, output = emi_intens_file, currency = currency_unit)
        ee.printEmissionConvRates(year, de_conv_file, emit_unit = emission_unit, curr_unit = currency_unit)
    end
    ee.readEmissionIntensity(year, mdr.nations, deSectorFile, deIntensityFile, emit_unit = emission_unit, curr_unit = currency_unit)
    println(" complete")
end

println(" Emission calculation: ")
path = Base.source_dir()*"/data/emission/"
if !testMode; nat_list = mdr.nations else nat_list = test_nats end
ns = length(nat_list)

# read consuming time
tf = open(indexPath * string(year) * "_" * Qtable * "_time.txt")
tp = Dict{String, Float64}()
readline(tf)
for l in eachline(tf); s = split(l, '\t'); tp[s[1]] = parse(Float64, s[2]) end
close(tf)
tp_avg = mean([tp[n] for n in filter(n->n in nat_list, collect(keys(tp)))])
tp_sum = sum([haskey(tp, n) ? tp[n] : tp_avg for n in nat_list])

ttp = 0
st = time()    # check start time
for i = 1:ns
    n = nat_list[i]
    emissionFile = path * string(year) * "_" * n * "_hhs_"*scaleTag*"IE_"*Qtable*".txt"
    deFile = path * string(year) * "_" * n * "_hhs_"*scaleTag*"DE.txt"
    print("\t", n, ":")
    print(" data"); ee.getDomesticData(year, n, mdr.expTable, mdr.hhsList)
    if DE_mode; print(", DE")
        # ee.calculateDirectEmission(year, n)
        # ee.printDirectEmissions(year, deFile)
        de_conc_file = concPath  * string(year) * "_"* n * "_DE_conc_mat.txt"
        ee.buildDeConcMat(year, n, deSectorFile, deHbsLinkFile, norm = true, energy_wgh = true, output = de_conc_file, group = reg_group)
        ee.calculateDirectEmission(year, n, sparseMat = false, enhance = false, full = true)
        ee.printEmissions(year, deFile, mode = "de")
    end
    if IE_mode
        wgh_conc_file = concPath * scaleTag * n * "_Eora_weighted_concordance_table.csv"
        print(", concordance"); ee.buildWeightedConcMat(year, ee.abb[mdr.nationNames[n]], output = wgh_conc_file)
        print(", IE"); ee.calculateIndirectEmission(year, false, 0)
        ee.printIndirectEmissions(year, emissionFile)
    end

    if haskey(tp, n); global ttp += tp[n] else global ttp += tp_avg end
    elap = floor(Int, time() - st)
    (eMin, eSec) = fldmod(elap, 60)
    (eHr, eMin) = fldmod(eMin, 60)
    (rMin, rSec) = fldmod(floor(Int, elap/ttp*(tp_sum-ttp)), 60)
    (rHr, rMin) = fldmod(rMin, 60)
    println(", ",n," (",i,"/",ns,") ",eHr,":",eMin,":",eSec," elapsed, total ",rHr,":",rMin,":",rSec," remained")
end
println("complete")

println(" ... all complete")
