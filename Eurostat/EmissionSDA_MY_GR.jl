# Developed date: 11. Dec. 2021
# Last modified date: 27. Feb. 2023
# Subject: Structual Decomposition Analysis (grouped)
# Description: Process for Input-Output Structural Decomposition Analysis
#              reading and decomposing multi-year micro-data
#              by population density, CF level, and income level
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

all_wgh_mode = true    # apply all related sub-sectors for calculating substitution codes' concordance table
adjustConc = false
domestic_mode = false

removeNTZ = true
adjustNTZ = removeNTZ ? false : true

catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]
if codeSubst; substTag = "_subst" else substTag = "" end

SDA_test = false; sda_test_nats = ["RO", "SE", "SK"];
if SDA_test; test_tag = "_test" else test_tag = "" end

test_tag = ""
if length(ARGS) > 0
    sda_test_nats = map(x -> string(x), ARGS)
    test_tag = "_"* sda_test_nats[1] * (length(sda_test_nats)>1 ? "to" * sda_test_nats[end] : "")
elseif SDA_test
    test_tag = "_test"
end

reuse_mem = true
sda_mode = "penta"
# sda_mode = "hexa"
# sda_mode = "categorized"

sda_path = emissDataPath * "SDA/"

nt_lv0_mode = true          # nation level (NUTS lv0) SDA mode
pd_mode = false              # grouping by population density
cf_group = false             # grouping by CF per capita, stacked proportion
inc_group = false            # grouping by income per capita, stacked proportion
cf_boundary = true          # grouping by CF per capita, boundary
inc_boundary = true         # grouping by income per capita, boundary
ce_intgr_mode = "cf"        # "ie" (only indirect CE), "de" (only direct CE), or "cf" (integrage direct and indirect CEs)
HHs_sorting_mode = "hhs"    # houseolg sorting mode for grouping: "hhs" household account, "percap" individual account

pop_dens = pd_mode ? [1,2,3] : []   # [1] Densely populated, [2] Intermediate, [3] Sparsely populated
cf_gr = cf_group ? [0.1, 0.9, 1.0] : []
inc_gr = inc_group ? [0.1, 0.9, 1.0] : []
cf_bnd = cf_boundary ? [0, 8, 50] : []
inc_bnd = inc_boundary ? [0, 7000, 70000] : []

conc_mat = Dict{Int, Dict{String, Array{Float64,2}}}()
pos_cf = Dict{Int, Dict{String, Dict{String, Float64}}}()
pos_inc = Dict{Int, Dict{String, Dict{String, Float64}}}()

nt_tag = "NT" * (nt_lv0_mode ? "0" : string(nutsLv)) * "_"

for year in years

    global filePath, indexFilePath, microDataPath, extractedPath, emissDataPath, mrioPath, base_year
    global Qtable, scaleMode, scaleTag, nation, nutsLv, categories, subcat, adjustConc, domestic_mode
    global categoryFile, eustatsFile, cpi_file, concFiles, natLabels
    global CurrencyConv, erfile, PPPConv, pppfile, codeSubst, perCap, conc_mat
    global indexFilePath, catDepth, depthTag, codeSubst, substTag, grid_pop, pos_cf, pos_inc

    const_tag = ConstConv && year != base_year ? "_" * string(base_year) * "_constant" : ""

    println("[",year,"]")
    microDataPath *= string(year) * "/"
    domfile = indexFilePath * string(year) * "_domestic_sectors.csv"
    ctgfile = extractedPath * string(year) * "_Category_"*depthTag[catDepth]*".csv"
    hhsfile = extractedPath * string(year) * "_Households" * const_tag * ".csv"
    mmsfile = extractedPath * string(year) * "_Members" * const_tag * ".csv"
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

    ce_file_path = emissDataPath * string(year) * "/"
    de_file_tag = "_hhs_"*scaleTag*"DE.txt"
    DE_files = filter(f->startswith(f, string(year)) && endswith(f, de_file_tag), readdir(ce_file_path))
    de_nations = [f[6:7] for f in DE_files]
    DE_files = ce_file_path .* DE_files
    print(", DE"); ec.readEmissionData(year, de_nations, DE_files, mode = "de")

    if cf_group || inc_group || cf_boundary || inc_boundary
        ie_file_tag = "_hhs_"*scaleTag*"IE_"*Qtable*".txt"
        if year != base_year; ie_file_tag = replace(ie_file_tag, ".txt" => "_converted_" * string(base_year) * ".txt") end
        IE_files = filter(f->startswith(f, string(year)) && endswith(f, ie_file_tag), readdir(ce_file_path))
        ie_nations = [f[6:7] for f in IE_files]
        IE_files = ce_file_path .* IE_files
        print(", IE"); ec.readEmissionData(year, ie_nations, IE_files, mode = "ie")
        print(", CF"); ec.integrateCarbonFootprint(year, mode = ce_intgr_mode)
        print(", categorizing"); ec.categorizeHouseholdEmission(year, mode = ce_intgr_mode, output="", hhsinfo=false, nutsLv=1)

        if cf_group && HHs_sorting_mode == "hhs"
            pos_cf[year] = ec.sortHHsByStatus(year, ie_nations, mode = "cf", sort_mode="cf")[year]
        elseif cf_group && HHs_sorting_mode == "percap"
            pos_cf[year] = ec.sortHHsByStatus(year, ie_nations, mode = "cf", sort_mode="cfpc")[year]
        end
        if inc_group && HHs_sorting_mode == "hhs"
            pos_inc[year] = ec.sortHHsByStatus(year, ie_nations, mode = "cf", sort_mode="income")[year]
        elseif inc_group && HHs_sorting_mode == "percap"
            pos_inc[year] = ec.sortHHsByStatus(year, ie_nations, mode = "cf", sort_mode="income_pc")[year]
        end
    end
    print(", importing"); ed.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, nations = [])
    print(", convert NUTS"); ed.convertNUTS(year = year)
    print(", detect NUTS"); ed.storeNUTS(year, cat_data = ec)
    print(", nuts weight"); ed.storeNutsWeight(year = year)
    println(" ... completed")
end

println("[SDA process]")

pop_label = Dict(true => "_byGroup", false => "")
const_tag = ConstConv ? "_" * string(base_year) * "_constant" : ""
pl_chk = pd_mode || cf_group || inc_group || cf_boundary || inc_boundary
delta_file = sda_path * string(target_year) * "_" * string(base_year) * "_deltas_" * nt_tag * sda_mode * pop_label[pl_chk] * const_tag * test_tag * ".txt"
nats = ed.filterNations()
if length(ARGS) > 0 || SDA_test; nats = sda_test_nats end

if pd_mode
    print(" Filtering nation:")
    print(" by populatoin density"); nat_flt = ed.getNonPopDens(years, nats, pop_dens = pop_dens, remove = false)
    println(" ... complete")
end
print(" Integrate NUTS: from ", target_year, " to ", base_year)
ed.integrateNUTS(target_year, base_year, categoryFile, modify = true, pop_dens = pd_mode, nt0_mode = nt_lv0_mode)
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
        ed.decomposeFactorsByGroup(y, base_year, n, mrioPath, mode = sda_mode, pop_dens = pop_dens, visible = false,
                                    cf_intv = cf_gr, inc_intv = inc_gr, hpos_cf = pos_cf, hpos_inc = pos_inc,
                                    cf_bndr = cf_bnd, inc_bndr = inc_bnd, bndr_mode = HHs_sorting_mode)
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
