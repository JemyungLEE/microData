clearconsole()
cd(Base.source_dir())
include("MicroDataReader.jl")
include("ConcMatBuilder.jl")
include("EmissionEstimator.jl")
include("EmissionCategorizer.jl")
include("EmissionDecomposer_hexa.jl")
using .MicroDataReader
using .ConcMatBuilder
using .EmissionEstimator
using .EmissionCategorizer
using .EmissionDecomposer_hexa
mdr = MicroDataReader
cmb = ConcMatBuilder
ee = EmissionEstimator
ec = EmissionCategorizer
ed = EmissionDecomposer_hexa

years = [2010, 2015]
base_year = 2010

for year in years

    filePath = Base.source_dir() * "/data/"
    indexFilePath = filePath * "index/"
    microDataPath = filePath * "microdata/"
    extractedPath = filePath * "extracted/"
    emissDataPath = filePath* "emission/"

    # Qtable = "I_CHG_CO2"
    Qtable = "PRIMAP"
    scaleMode = true; if scaleMode; scaleTag = "Scaled_" else scaleTag = "" end

    nation = "Eurostat"
    nutsLv = 1

    categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                    "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
    subcat=""

    categoryFile = indexFilePath * "Eurostat_Index_ver4.3_NT0.xlsx"
    eustatsFile = indexFilePath * "EU_exp_COICOP.tsv"
    cpi_file = indexFilePath * "EU_hicp.tsv"

    concFiles = Dict(2010 => indexFilePath*"2010_EU_EORA_Conc_ver1.5.xlsx", 2015 => indexFilePath*"2015_EU_EORA_Conc_ver1.1.xlsx")
    natLabels = Dict(2010 => "Eurostat", 2015 => "Eurostat_2015")

    CurrencyConv = true; erfile = indexFilePath * "EUR_USD_ExchangeRates.txt"
    PPPConv = false; pppfile = indexFilePath * "PPP_ConvertingRates.txt"

    codeSubst = true        # recommend 'false' for depth '1st' as there is nothing to substitute
    perCap = true

    eoraRevised = true

    catDepth = 4
    depthTag = ["1st", "2nd", "3rd", "4th"]
    if codeSubst; substTag = "_subst" else substTag = "" end

    println("[",year,"]")
    microDataPath *= string(year) * "/"
    ctgfile = extractedPath * string(year) * "_Category_"*depthTag[catDepth]*".csv"
    hhsfile = extractedPath * string(year) * "_Households.csv"
    mmsfile = extractedPath * string(year) * "_Members.csv"
    expfile = extractedPath * string(year) * "_" * scaleTag * "Expenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
    sbstfile = extractedPath * string(year) * "_SubstituteCodes_"*depthTag[catDepth]*".csv"
    if year == 2010; hhsfile = replace(hhsfile, ".csv" => "_NT0.csv") end

    print(" Category codes reading:")
    mdr.readCategory(year, categoryFile, depth=catDepth, catFile=ctgfile, coicop=scaleMode)
    println(" ... complete")

    print(" Micro-data reading:")
    print(" hhs"); mdr.readPrintedHouseholdData(hhsfile)
    print(", mms"); mdr.readPrintedMemberData(mmsfile)
    if codeSubst; print(", subst"); mdr.readSubstCodesCSV(sbstfile) end
    print(", exp"); mdr.readPrintedExpenditureData(expfile, substitute=codeSubst, buildHhsExp=true)
    println(" ... complete")

    if CurrencyConv; print(" Currency exchanging: ")
        print(" exchange"); mdr.exchangeExpCurrency(erfile, year = year)
        print(" rebuild matrix"); mdr.buildExpenditureMatrix(year, substitute=codeSubst)
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
    print(", substitution"); cmb.addSubstSec(year, mdr.heSubst, mdr.heRplCd, mdr.heCats, exp_table = [])
    print(", normalization"); cmn = cmb.normConMat(year)   # {a3, conMat}
    println(" ... complete")

    print(" MRIO table reading:")
    if eoraRevised; eora_index = "../Eora/data/index/revised/" else eora_index = "../Eora/data/index/Eora_index.xlsx" end
    path = "../Eora/data/" * string(year) * "/" * string(year)
    print(" index"); ee.readIndexXlsx("../Eora/data/index/revised/", revised = eoraRevised, initiate = true)
    print(", IO table"); ee.readIOTables(year, path*"_eora_t.csv", path*"_eora_v.csv", path*"_eora_y.csv", path*"_eora_q.csv")
    print(", rearrange"); ee.rearrangeIndex(qmode=Qtable); ee.rearrangeTables(year, qmode=Qtable)
    println(" ... complete")

    print(" Data import:")
    print(" sector"); ee.getSectorData(year, mdr.heCodes, mdr.heSubst)
    print(", assemble"); ee.assembleConcMat(year, cmn)
    print(", category"); ec.readCategoryData(categoryFile, year, nutsLv, except=["None"], subCategory=subcat, nuts3pop=true)
    ec.setCategory(categories)

    print(", household"); ec.readHouseholdData(hhsfile, period = "daily", remove = true)
    print(", nuts weight"); ec.calculateNutsPopulationWeight(year = year)

    de_file_path = emissDataPath * string(year) * "/"
    DE_files = filter(f->startswith(f, string(year)) && endswith(f, "_hhs_"*scaleTag*"DE.txt"), readdir(de_file_path))
    de_nations = [f[6:7] for f in DE_files]
    DE_files = de_file_path .* DE_files

    print(", DE"); ec.readEmissionData(year, de_nations, DE_files, mode = "de")
    print(", importing"); ed.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, nations = [])

    print(", detect NUTS"); ed.storeNUTS(year, cat_data = ec)
    print(", nuts weight"); ed.storeNutsWeight(year = year)
    println(" ... completed")
end

SDA_test = false; test_nats = ["BE", "BG"]
factorPrintMode = false
mem_clear_mode = true

filePath = Base.source_dir() * "/data/"
indexFilePath = filePath * "index/"
emissDataPath = filePath* "emission/"
sda_path = emissDataPath * "SDA/"
factorPath = sda_path * "factors/"

target_year = 2015
println(" SDA process:")
n_factor = 6
pop_dens = 0        # [1] Densely populated, [2] Intermediate, [3] Sparsely populated
pop_label = Dict(0 => "", 1 => "_dense", 2 => "_inter", 3 => "_sparse")
delta_file = sda_path * string(target_year) * "_" * string(base_year) * "_deltas_hexa" * pop_label[pop_dens] *".txt"
if SDA_test; nats = test_nats else nats = ed.filterNations() end

if mem_clear_mode && !(pop_dens in [1, 2, 3]); mdr.initVars() end
if pop_dens in [1, 2, 3]
    nuts_lv = 3
    pop_dens_file = indexFilePath * "Eurostat_population_density.tsv"
    print(" Population density:")
    for year in years
        print(" $year read density"); ed.readPopDensity(year, pop_dens_file)
        print(", filter population"); ed.filterPopByDensity(year, nuts_lv = nuts_lv)
        print(" and nation"); filter!(x -> x in ed.filterNonPopDens(year, pop_dens = pop_dens), nats)
        # print(", print results"); ed.printPopByDens(year, indexFilePath * string(year) * "_Population_by_density_NUTS_Lv" * string(nuts_lv) * ".txt")
    end
    println(" ... complete")
end

println(" ", nats)
for n in nats
    print(n, ":")
    print(" decomposing")
    for y in years
        conc_mat_wgh = ee.buildWeightedConcMat(y, ee.abb[mdr.nationNames[n]])[1][:,:]
        ed.storeConcMat(y, n, conc_mat_wgh)
        ed.decomposeFactors(y, base_year, n, visible = false, pop_nuts3 = false, pop_dens = pop_dens)
        if factorPrintMode
            if y != base_year; ed.printLmatrix(y, factorPath, nation = n, base_year = base_year) end
            ed.printFactors(factorPath, year = y, nation = n)
        end
    end

    print(", factors")
    ed.prepareDeltaFactors(target_year, base_year, nation = n)
    print(", sda")
    ed.structuralAnalysis(target_year, base_year, n, n_factor)
    println()
    if mem_clear_mode; ed.clearFactors(nation = n) end
end
if factorPrintMode; print(", printing")
    for y in years
        ed.printNUTS(y, factorPath)
        if y == base_year; ed.printLmatrix(y, factorPath) end
    end
end
println(" ... completed")

print(" SDA results printing: "); ed.printDelta(delta_file)
println(" ... completed")
