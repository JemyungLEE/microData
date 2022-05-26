# Developed date: 31. Aug. 2021
# Last modified date: 26. May. 2022
# Subject: Structual Decomposition Analysis
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

base_year, target_year = 2010, 2015
years = [base_year, target_year]

# opr_mode = "pc"
opr_mode = "server"

if opr_mode == "pc"
    mem_clear_mode = true
    # clearconsole()
    filePath = Base.source_dir() * "/data/"
    mrioPath = "../Eora/data/"
elseif opr_mode == "server"
    mem_clear_mode = false
    filePath = "/import/mary/lee/Eurostat/data/"
    mrioPath = "/import/mary/lee/Eora/data/"
end

indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissDataPath = filePath* "emission/"

Qtable = "PRIMAP"
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

codeSubst = true        # recommend 'false' for depth '1st' as there is nothing to substitute
perCap = true
grid_pop = true

all_wgh_mode = true    # apply all related sub-sectors for calculating substitution codes' concordance table
adjustConc = false
domestic_mode = false

removeNTZ = true
adjustNTZ = removeNTZ ? false : true

catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]
if codeSubst; substTag = "_subst" else substTag = "" end

SDA_test = false; test_nats = ["FR"];
if SDA_test; test_tag = "_test" else test_tag = "" end

reuse_mem = true
sda_mode = "penta"
# sda_mode = "hexa"
# sda_mode = "categorized"

sda_path = emissDataPath * "SDA/"

conc_mat = Dict{Int, Dict{String, Array{Float64,2}}}()

for year in years

    global filePath, indexFilePath, microDataPath, extractedPath, emissDataPath, mrioPath
    global Qtable, scaleMode, scaleTag, nation, nutsLv, categories, subcat, adjustConc, domestic_mode
    global categoryFile, eustatsFile, cpi_file, concFiles, natLabels
    global CurrencyConv, erfile, PPPConv, pppfile, codeSubst, perCap, conc_mat
    global catDepth, depthTag, codeSubst, substTag, grid_pop

    println("[",year,"]")
    microDataPath *= string(year) * "/"
    domfile = indexFilePath * string(year) * "_domestic_sectors.csv"
    ctgfile = extractedPath * string(year) * "_Category_"*depthTag[catDepth]*".csv"
    hhsfile = extractedPath * string(year) * "_Households.csv"
    mmsfile = extractedPath * string(year) * "_Members.csv"
    expfile = extractedPath * string(year) * "_" * scaleTag * "Expenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
    sbcdsfile = extractedPath * string(year) * "_SubstituteCodes_" * depthTag[catDepth] * ".csv"
    sbctgfile = extractedPath * string(year) * "_Category_" * depthTag[catDepth] * "_subst.csv"

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
    print(", normalization"); conc_mat[year] = cmb.normConMat(year)   # {a3, conMat}
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
    if !domestic_mode; print(", assembe Conc_mat"); ee.assembleConcMat(year, conc_mat[year], dom_nat = "")
    else print(", read domestic sectors"); ee.readDomesticSectors(year, domfile)
    end
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

println("[SDA process]")
nt_lv0_mode = false       # nation level (NUTS lv0) SDA mode
pop_dens = 0        # [1] Densely populated, [2] Intermediate, [3] Sparsely populated
pop_label = Dict(0 => "", 1 => "_dense", 2 => "_inter", 3 => "_sparse")
delta_file = sda_path * string(target_year) * "_" * string(base_year) * "_deltas_" * sda_mode * pop_label[pop_dens] * test_tag* ".txt"
nats = ed.filterNations()
if SDA_test; nats = test_nats end

if pop_dens in [1, 2, 3]
    print(" Population density:")
    nat_flt = ed.filterNonPopDens(years, nats, pop_dens = pop_dens)
    print(", filter nation"); filter!(x -> x in nat_flt, nats)
    println(" ... complete")
end
print(" Integrate NUTS: from ", target_year, " to ", base_year)
ed.integrateNUTS(target_year, base_year, categoryFile, modify = true, pop_dens = false, nt0_mode = nt_lv0_mode)
println(" ... complete")

st = time()
println(" ", nats)
ed.printDeltaTitle(delta_file, mode = sda_mode, cf_print = true, st_print = true)
for n in nats
    print(n, ":")
    print(" decomposing")
    for y in years
        print("_", y)
        if domestic_mode
            ee.getDomesticData(y, n, mdr.expTable, mdr.hhsList)
            conc_mat_org = ee.assembleConcMat(y, conc_mat[y], dom_nat = n)[1]
        else conc_mat_org = []
        end
        conc_mat_wgh = ee.buildWeightedConcMat(y, ee.abb[mdr.nationNames[n]], adjust = adjustConc)[1]
        ed.storeConcMat(y, n, conc_mat_wgh, conc_mat_nw = conc_mat_org)
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
