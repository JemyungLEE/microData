# Developed date: 16. Nov. 2021
# Last modified date: 20. Aug. 2022
# Subject: Structual Decomposition Analysis (server version)
# Description: Process for Input-Output Structural Decomposition Analysis
#              reading and decomposing multi-year micro-data
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())
include("MicroDataReader.jl")
include("ConcMatBuilder.jl")
include("EmissionEstimator.jl")
include("EmissionCategorizer.jl")
include("EmissionDecomposer.jl")
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

years = [2010, 2015]
base_year = 2010

# opr_mode = "pc"
opr_mode = "server"

if opr_mode == "pc"
    test_mode = true
    mem_clear_mode = true
    # clearconsole()
    filePath = Base.source_dir() * "/data/"
    mrioPath = "../Eora/data/"
elseif opr_mode == "server"
    test_mode = true
    mem_clear_mode = false
    filePath = "/import/mary/lee/Eurostat/data/"
    mrioPath = "/import/mary/lee/Eora/data/"
end
nats_test = []

indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissDataPath = filePath* "emission/"

Qtable = "I_CHG_CO2"
# Qtable = "PRIMAP"

scaleMode = true; if scaleMode; scaleTag = "Scaled_" else scaleTag = "" end

nation = "Eurostat"
nutsLv = 1

categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
subcat=""

categoryFile = indexFilePath * "Eurostat_Index_ver5.0.xlsx"
eustatsFile = indexFilePath * "EU_exp_COICOP.tsv"
cpi_file = indexFilePath * "EU_hicp.tsv"

concFiles = Dict(2010 => indexFilePath*"2010_EU_EORA_Conc_ver1.5.xlsx", 2015 => indexFilePath*"2015_EU_EORA_Conc_ver1.1.xlsx")
natLabels = Dict(2010 => "Eurostat", 2015 => "Eurostat_2015")

CurrencyConv = true; erfile = indexFilePath * "EUR_USD_ExchangeRates.txt"
PPPConv = false; pppfile = indexFilePath * "PPP_ConvertingRates.txt"
ConstConv = true        # convert household income and expenditure to constant year (=base_year) price

codeSubst = true        # recommend 'false' for depth '1st' as there is nothing to substitute
perCap = true
grid_pop = true

removeNTZ = true
adjustNTZ = removeNTZ ? false : true
all_wgh_mode = true    # apply all related sub-sectors for calculating substitution codes' concordance table

catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]
if codeSubst; substTag = "_subst" else substTag = "" end

for year in years

    global filePath, indexFilePath, microDataPath, extractedPath, emissDataPath
    global Qtable, scaleMode, scaleTag, nation, nutsLv, categories, subcat
    global categoryFile, eustatsFile, cpi_file, concFiles, natLabels
    global CurrencyConv, erfile, PPPConv, pppfile, codeSubst, perCap
    global indexFilePath, catDepth, depthTag, codeSubst, substTag, grid_pop, mrioPath

    const_tag = ConstConv && year != base_year ? "_" * string(base_year) * "_constant" : ""

    println("[",year,"]")
    microDataPath *= string(year) * "/"
    ctgfile = extractedPath * string(year) * "_Category_"*depthTag[catDepth]*".csv"
    hhsfile = extractedPath * string(year) * "_Households" * const_tag * ".csv"
    mmsfile = extractedPath * string(year) * "_Members" * const_tag * ".csv"
    expfile = extractedPath * string(year) * "_" * scaleTag * "Expenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
    sbcdsfile = extractedPath * string(year) * "_SubstituteCodes_" * depthTag[catDepth] * ".csv"
    sbctgfile = extractedPath * string(year) * "_Category_" * depthTag[catDepth] * "_subst.csv"
# if year == 2010; hhsfile = replace(hhsfile, ".csv" => "_NT0.csv") end

    print(" Category codes reading:")
    mdr.readCategory(year, categoryFile, depth=catDepth, catFile=ctgfile, coicop=scaleMode)
    ec.readCategoryData(categoryFile, year, nutsLv, except=["None"], subCategory=subcat)
    ec.setCategory(categories)
    println(" ... complete")

    print(" Micro-data reading:")
    print(" hhs")
    ec.readHouseholdData(hhsfile, period = "annual", alter=true, remove_nt = true, remove_z = removeNTZ)
    mdr.readPrintedHouseholdData(hhsfile, regions = ec.regList[year])
    # print(", mms"); mdr.readPrintedMemberData(mmsfile)
    if codeSubst; print(", subst"); mdr.readSubstCodesCSV(year, sbctgfile, sbcdsfile) end
    print(", exp"); mdr.readPrintedExpenditureData(expfile, substitute=codeSubst, buildHhsExp=true)
    println(" ... complete")

    if CurrencyConv; print(" Currency exchanging: ")
        print(" exchange"); mdr.exchangeExpCurrency(erfile, year = year)
        # print(" rebuild matrix"); mdr.buildExpenditureMatrix(year, substitute=codeSubst)
        println(" ... complete")
    end
    if PPPConv; print(" PPP converting:")
        mdr.convertToPPP(pppfile)
        println(" ... complete")
    end

    if year != base_year; print(" CPI scaling:")
        print(" read CPIs"); mdr.readCPIs(years, cpi_file, idx_sep = ',', freq="A", unit="INX_A_AVG", topLev = "EU")
        print(", scaling"); mdr.scalingByCPI(year, base_year, codeDepth=0, topLev = "EU", subst = codeSubst)
        print(", rebuild matrix"); mdr.buildExpenditureMatrix(year, substitute=codeSubst)
        println(" ... complete")
    end

    print(" Concordance matrix building:")
    print(" concordance"); cmb.readXlsxData(year, concFiles[year], nation, nat_label = natLabels[year])
    print(", matrix"); cmb.buildConMat(year)
    print(", substitution"); cmb.addSubstSec(year, mdr.heSubst, all_wgh_mode ? mdr.heSubHrr : mdr.heRplCd, mdr.heCats, exp_table = mdr.expTable, norm = true, wgh_all = all_wgh_mode)
    print(", normalization"); cmn = cmb.normConMat(year)   # {a3, conMat}
    print(", memory clear"); cmb.initVars(year = year)
    println(" ... complete")

    print(" MRIO table reading:")
    eora_index = mrioPath * "index/"
    m_path = mrioPath * string(year) * "/" * string(year)
    print(" index"); ee.readIOindex(eora_index)
    print(", IO table"); ee.readIOTables(year, m_path*"_eora_t.csv", m_path*"_eora_v.csv", m_path*"_eora_y.csv", m_path*"_eora_q.csv")
    print(", rearrange"); ee.rearrangeMRIOtables(year, qmode=Qtable)
    println(" ... complete")

    print(" Data import:")
    print(" sector"); ee.getSectorData(year, mdr.heCodes, mdr.heSubst)
    print(", assemble"); ee.assembleConcMat(year, cmn)
    print(", population"); ec.readPopulation(year, categoryFile, nuts_lv = nutsLv)
    print(", gridded population"); ec.readPopGridded(year, categoryFile, nuts_lv = [nutsLv], adjust = true)
    print(", nuts weight"); ec.calculateNutsPopulationWeight(year = year, pop_dens = grid_pop, adjust = adjustNTZ)

    de_file_path = emissDataPath * string(year) * "/"
    DE_files = filter(f->startswith(f, string(year)) && endswith(f, "_hhs_"*scaleTag*"DE.txt"), readdir(de_file_path))
    de_nations = [f[6:7] for f in DE_files]
    DE_files = de_file_path .* DE_files

    print(", DE"); ec.readEmissionData(year, de_nations, DE_files, mode = "de")
    print(", importing"); ed.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, nations = [])
    print(", convert NUTS"); ed.convertNUTS(year = year)
    print(", detect NUTS"); ed.storeNUTS(year, cat_data = ec)
    print(", nuts weight"); ed.storeNutsWeight(year = year)
    println(" ... completed")
end

reuse_mem = true
# sda_mode = "penta"
# sda_mode = "hexa"
sda_mode = "categorized"

sda_path = emissDataPath * "SDA/"
factorPath = sda_path * "factors/"

target_year = 2015
println("[SDA process]")
pop_dens = 0        # [1] Densely populated, [2] Intermediate, [3] Sparsely populated
pop_label = Dict(0 => "", 1 => "_dense", 2 => "_inter", 3 => "_sparse")
const_tag = ConstConv ? "_" * string(base_year) * "_constant" : ""
nats = ed.filterNations()

if length(ARGS) > 0; nats = map(x -> string(x), ARGS)
elseif test_mode; nats = nats_test
end

sda_file_tag = "_"* nats[1] * (length(nats)>1 ? "to" * nats[end] : "")
delta_file = sda_path * string(target_year) * "_" * string(base_year) * "_deltas_" * sda_mode * pop_label[pop_dens] * const_tag * sda_file_tag* ".txt"

if pop_dens in [1, 2, 3]
    print(" Population density:")
    nat_flt = ed.filterNonPopDens(years, nats, pop_dens = pop_dens)
    print(", filter nation"); filter!(x -> x in nat_flt, nats)
    println(" ... complete")
end
print(" Integrate NUTS: from ", target_year, " to ", base_year)
ed.integrateNUTS(target_year, base_year, categoryFile, modify = true, pop_dens = true)
println(" ... complete")

st = time()
println(" ", nats)
ed.printDeltaTitle(delta_file, mode = sda_mode, cf_print = true, st_print = true)
for n in nats
    print(n, ":")
    print(" decomposing")
    for y in years
        conc_mat_wgh = ee.buildWeightedConcMat(y, ee.abb[mdr.nationNames[n]], adjust = false)[1]
        ed.storeConcMat(y, n, conc_mat_wgh)
        ed.decomposeFactors(y, base_year, n, mrioPath, visible = false, pop_dens = pop_dens, mode = sda_mode)
    end

    print(", factors"); ed.prepareDeltaFactors(target_year, base_year, nation = n, mode = sda_mode, reuse = reuse_mem)
    print(", sda"); ed.structuralAnalysis(target_year, base_year, n, mode = sda_mode, reuse = reuse_mem)
    print(", printing"); ed.printDeltaValues(delta_file, n, mode = sda_mode, cf_print = true, st_print = true)

    if mem_clear_mode
        ec_clear = (n == nats[end])
        mdr.initVars(year = years, nation = n)
        ec.initVars(year = years, nation = n, clear_all = ec_clear)
        ed.clearFactors(nation = n)
    end

    elap = floor(Int, time() - st); (eMin, eSec) = fldmod(elap, 60); (eHr, eMin) = fldmod(eMin, 60)
    println(",\t", eHr, ":", eMin, ":", eSec, " elapsed")
end

println(" ... completed")
