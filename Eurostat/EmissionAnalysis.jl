# Developed date: 28. Jul. 2020
# Last modified date: 6. Jun. 2022
# Subject: Estimate carbon footprint by final demands of Eora
# Description: Calculate carbon emissions by utilizing Eora T, V, Y, and Q tables.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("MicroDataReader.jl")
include("ConcMatBuilder.jl")
include("EmissionEstimator.jl")
include("EmissionCategorizer.jl")
include("EmissionDecomposer.jl")

using Statistics
using .MicroDataReader
using .ConcMatBuilder
using .EmissionEstimator
using .EmissionCategorizer
using .EmissionDecomposer

mdr = MicroDataReader
cmb = ConcMatBuilder
ee = EmissionEstimator
ec = EmissionCategorizer
ed = EmissionDecomposer

# opr_mode = "pc"
opr_mode = "server"

if opr_mode == "pc"
    clearconsole()
    filePath = Base.source_dir() * "/data/"
    mrioPath = "../Eora/data/"
elseif opr_mode == "server"
    filePath = "/import/mary/lee/Eurostat/data/"
    mrioPath = "/import/mary/lee/Eora/data/"
end

IE_mode = true             # indirect carbon emission estimation
DE_mode = false              # direct carbon emission estimation
DE_factor_estimate = true   # [true] estimate DE factors from IEA datasets, [false] read DE factors

nation = "Eurostat"
year = 2010
catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]
emission_unit = "tCO2"
currency_unit = "EUR"
reg_group = "EU"

indexPath = filePath * "index/"
extrPath = filePath * "extracted/"
concPath = indexPath * "concordance/"
emit_path = indexPath * "de/"
microDataPath = filePath * "microdata/" * string(year) * "/"
cePath = filePath * "emission/"

categoryFile = indexPath * "Eurostat_Index_ver5.0.xlsx"
erfile = indexPath * "EUR_USD_ExchangeRates.txt"

PPPConv = false; pppfile = indexPath * "PPP_ConvertingRates.txt"

concFiles = Dict(2010 => indexPath*"2010_EU_EORA_Conc_ver1.5.xlsx", 2015 => indexPath*"2015_EU_EORA_Conc_ver1.1.xlsx")
natLabels = Dict(2010 => "Eurostat", 2015 => "Eurostat_2015")
conMatFile = concPath * string(year) * "_ConcMat.txt"
conSumMatFile = concPath * string(year) * "_ConcSumMat.txt"
deSectorFile = emit_path * "DE_sectors.txt"
deHbsLinkFile = emit_path * "DE_linked_sectors_" * string(year) * ".txt"
deIntensityFile = emit_path * "Emission_converting_rate_"*string(year)*"_EU.txt"
cpi_file = indexPath * "EU_hicp.tsv"
domfile = indexPath * string(year) * "_domestic_sectors.csv"

Qtable = "I_CHG_CO2"
# Qtable = "PRIMAP"

abrExpMode = false
substMode = true
scaleMode = true
cpiScaling = true; base_year = 2010
CurrencyConv = true

all_wgh_mode = true    # apply all related sub-sectors for calculating substitution codes' concordance table
adjustConc = false
domestic_mode = false

testMode = false; test_nats = nats = ["BE", "BG"]

if IE_mode && !DE_mode; CurrencyConv = true
elseif !IE_mode && DE_mode; CurrencyConv = false
elseif IE_mode && DE_mode; println("Be careful selecting emission estimation mode.")
end
if substMode; substTag = "_subst" else substTag = "" end
if scaleMode; scaleTag = "Scaled_" else scaleTag = "" end

hhsfile = extrPath * string(year) * "_Households.csv"
mmsfile = extrPath * string(year) * "_Members.csv"
expfile = extrPath * string(year) * "_" * scaleTag*"Expenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
ctgfile = extrPath * string(year) * "_Category_" * depthTag[catDepth] * ".csv"
sttfile = extrPath * string(year) * "_MicroData_Statistics_"*scaleTag*depthTag[catDepth]*substTag*".csv"
sbcdsfile = extrPath * string(year) * "_SubstituteCodes_" * depthTag[catDepth] * ".csv"
sbctgfile = extrPath * string(year) * "_Category_" * depthTag[catDepth] * "_subst.csv"

println("[Process] year = ", year)
print(" Category codes reading: ")
mdr.readCategory(year, categoryFile, depth=catDepth, catFile=ctgfile, coicop=scaleMode, inclAbr=abrExpMode)
println("completed")

print(" Micro-data reading:")
print(" households"); mdr.readPrintedHouseholdData(hhsfile)
# print(", members"); mdr.readPrintedMemberData(mmsfile)
if substMode; print(" substitutes"); mdr.readSubstCodesCSV(year, sbctgfile, sbcdsfile) end
print(", expenditures(", rsplit(expfile, "/",limit=2)[end], ")")
mdr.readPrintedExpenditureData(expfile, substitute=substMode, buildHhsExp=true)
print(", sector data"); ee.getSectorData(year, mdr.heCodes, mdr.heSubst)
println(" ... completed")

if CurrencyConv; print(" Currency exchanging:")
    print(" USD transform"); mdr.exchangeExpCurrency(erfile)
    # print(", rebuild matrix"); mdr.buildExpenditureMatrix(year, substitute=substMode)
    println(" complete")
end
if PPPConv; print(" PPP converting: "); mdr.convertToPPP(pppfile); println("complete") end

if cpiScaling && year != base_year; print(" CPI scaling:")
    print(" read CPIs"); mdr.readCPIs([base_year, year], cpi_file, idx_sep = ',', freq="A", unit="INX_A_AVG", topLev = "EU")
    if IE_mode
        print(", scaling"); mdr.scalingByCPI(year, base_year, codeDepth=0, topLev = "EU", subst = substMode)
        print(", rebuild matrix"); mdr.buildExpenditureMatrix(year, substitute=substMode)
    end
    println(" ... complete")
end

if IE_mode
    # Converting process of Eora final demand data to India micro-data format
    print(" Concordance matrix building:")
    print(" concordance"); cmb.readXlsxData(year, concFiles[year], nation, nat_label = natLabels[year], domestic_codes = [])
    print(", matrix"); cmb.buildConMat(year)
    print(", substitution"); cmb.addSubstSec(year, mdr.heSubst, all_wgh_mode ? mdr.heSubHrr : mdr.heRplCd, mdr.heCats, exp_table = mdr.expTable, norm = true, wgh_all = all_wgh_mode)
    print(", normalization"); cmn = cmb.normConMat(year, domestic_nat = "")   # {a3, conMat}
    print(", printing"); cmb.printConMat(year, conMatFile, nation, norm = true, categ = true)
    # cmb.printSumNat(year, conSumMatFile, nation, norm = true)
    println(" complete")

    # Eora household's final-demand import sector data reading process
    print(" MRIO table reading:")
    eora_index = mrioPath * "index/"
    path = mrioPath * string(year) * "/" * string(year)
    print(" index"); ee.readIOindex(eora_index)
    print(", table"); ee.readIOTables(year, path*"_eora_t.csv", path*"_eora_v.csv", path*"_eora_y.csv", path*"_eora_q.csv")
    print(", rearrange"); ee.rearrangeMRIOtables(year, qmode=Qtable)
    if !domestic_mode; print(", assembe Conc_mat"); ee.assembleConcMat(year, cmn, dom_nat = "")
    else print(", read domestic sectors"); ee.readDomesticSectors(year, domfile)
    end
    if !cpiScaling || year == base_year; print(", Leontief matrix"); ee.calculateLeontief(year) end
    println(" ... complete")
end
if DE_mode
    print(" Direct emission converting indices reading:")
    print(" source data")
    ee.readEmissionData(year, mdr.nationNames, emit_path, output_path = emit_path * "EU/", output_tag = reg_group, integrate = true, cpi_scaling = true, cpi_base = base_year, cpi_vals = mdr.cpis)
    if DE_factor_estimate
        print(", estimation")
        price_file = emit_path * "EU/" * "Price_" * string(year) * "_EU_" * currency_unit * ".txt"
        emi_intens_file = emit_path * "EU/" * "Emission_intensity_" * string(year) * "_EU_tCO2_per_" * currency_unit* ".txt"
        de_conv_file = emit_path * "Emission_converting_rate_" * string(year) * "_EU.txt"
        ee.exchangeEmCurrency(year, erfile, target = currency_unit, origin = "USD", output = price_file)
        ee.calculateEmissionRates(year, output = emi_intens_file, currency = currency_unit)
        ee.printEmissionConvRates(year, de_conv_file, emit_unit = emission_unit, curr_unit = currency_unit)
    end
    print(", intensity")
    ee.readEmissionIntensity(year, mdr.nations, deSectorFile, deIntensityFile, emit_unit = emission_unit, curr_unit = currency_unit)
    println(" ... complete")
end

println(" Emission calculation: ")
if !testMode; nat_list = mdr.nations else nat_list = test_nats end
ns = length(nat_list)

# read consuming time
time_file = indexPath * string(year) * "_" * Qtable * "_time.txt"
if cpiScaling && year != base_year; time_file = replace(time_file, "_time" => "_converted_time") end
tf = open(time_file)
tp = Dict{String, Float64}()
readline(tf)
for l in eachline(tf); s = split(l, '\t'); tp[s[1]] = parse(Float64, s[2]) end
close(tf)
tp_avg = mean([tp[n] for n in filter(n->n in nat_list, collect(keys(tp)))])
tp_sum = sum([haskey(tp, n) ? tp[n] : tp_avg for n in nat_list])

t_bp, t_tax, t_sub, v_bp, y_bp = [], [], [], [], []
if IE_mode && cpiScaling && year != base_year
    println(" read mrio"); t_bp, t_tax, t_sub, v_bp, y_bp = ed.setMrioTables(year, mrioPath)
end

ttp = 0
st = time()    # check start time
for i = 1:ns
    global t_bp, t_tax, t_sub, v_bp, y_bp
    n = nat_list[i]
    emissionFile = cePath * string(year) * "_" * n * "_hhs_"*scaleTag*"IE_"*Qtable*".txt"
    if cpiScaling && year != base_year; emissionFile = replace(emissionFile, ".txt" => "_converted_" * string(base_year) * ".txt") end
    deFile = cePath * string(year) * "_" * n * "_hhs_"*scaleTag*"DE.txt"
    print("\t", n, ":")
    print(" data"); ee.getDomesticData(year, n, mdr.expTable, mdr.hhsList)
    if DE_mode; print(", DE")
        de_conc_file = concPath  * string(year) * "_"* n * "_DE_conc_mat.txt"
        ee.buildDeConcMat(year, n, deSectorFile, deHbsLinkFile, norm = true, energy_wgh = true, output = de_conc_file, group = reg_group)
        ee.calculateDirectEmission(year, n, sparseMat = false, enhance = false, full = true)
        ee.printEmissions(year, deFile, mode = "de")
    end
    if IE_mode
        wgh_conc_file = concPath * scaleTag * n * "_Eora_weighted_concordance_table.csv"
        print(", concordance")
        if domestic_mode; ee.assembleConcMat(year, cmn, dom_nat = n) end
        ee.buildWeightedConcMat(year, ee.abb[mdr.nationNames[n]], output = wgh_conc_file, adjust = adjustConc)
        if cpiScaling && year != base_year; print(", mrio converting")
            ed.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, nations = [])
            ed.convertTable(year, n, base_year, mrioPath, t_bp = t_bp, t_tax = t_tax, t_sub = t_sub, v_bp = v_bp, y_bp = y_bp)
            ee.mTables[year].t = ed.mrio_tabs_conv[year][n].t
            ee.mTables[year].v = ed.mrio_tabs_conv[year][n].v
            ee.mTables[year].y = ed.mrio_tabs_conv[year][n].y
            ee.calculateLeontief(year)
        end
        print(", IE"); ee.calculateIndirectEmission(year, false, 0)
        ee.printIndirectEmissions(year, emissionFile)
        if year != base_year; ed.clearFactors(nation = n) end
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
