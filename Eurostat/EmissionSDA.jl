# Developed date: 30. Jul. 2021
# Last modified date: 16. Aug. 2021
# Subject: Structual Decomposition Analysis
# Description: Process for Input-Output Structural Decomposition Analysis
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
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

filePath = Base.source_dir() * "/data/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissDataPath = filePath* "emission/"
sda_path = emissDataPath * "SDA/"
factorPath = sda_path * "factors/"

# Qtable = "I_CHG_CO2"
Qtable = "PRIMAP"
scaleMode = true; if scaleMode; scaleTag = "Scaled_" else scaleTag = "" end

nation = "Eurostat"
nutsLv = 1

categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
subcat=""
# subcat="Food"
foodCategories=["Grain","Vegetable","Fruit","Dairy","Beef","Pork","Poultry","Other meat","Fish",
                "Alcohol","Other beverage","Confectionery","Restaurant","Other food","Food"]

categoryFile = indexFilePath * "Eurostat_Index_ver4.2_NT0.xlsx"
eustatsFile = indexFilePath * "EU_exp_COICOP.tsv"
cpi_file = indexFilePath * "EU_hicp.tsv"

concFiles = Dict(2010 => indexFilePath*"2010_EU_EORA_Conc_ver1.5.xlsx", 2015 => indexFilePath*"2015_EU_EORA_Conc_ver1.1.xlsx")
natLabels = Dict(2010 => "Eurostat", 2015 => "Eurostat_2015")

CurrencyConv = true; erfile = indexFilePath * "EUR_USD_ExchangeRates.txt"
PPPConv = false; pppfile = indexFilePath * "PPP_ConvertingRates.txt"

codeSubst = true        # recommend 'false' for depth '1st' as there is nothing to substitute
perCap = true

eoraRevised = true

factorEstimateMode = false
factorPrintMode = false

SDA_mode = true

year = 2015
base_year = 2010
catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]
if codeSubst; substTag = "_subst" else substTag = "" end

microDataPath *= string(year) * "/"
# microDataPath = [microDataPath*"BE", microDataPath*"HU"]

ctgfile = extractedPath * string(year) * "_Category_"*depthTag[catDepth]*".csv"
hhsfile = extractedPath * string(year) * "_Households.csv"
mmsfile = extractedPath * string(year) * "_Members.csv"
expfile = extractedPath * string(year) * "_" * scaleTag * "Expenditure_matrix_"*depthTag[catDepth]*substTag*".csv"
sbstfile = extractedPath * string(year) * "_SubstituteCodes_"*depthTag[catDepth]*".csv"
if year == 2010; hhsfile = replace(hhsfile, ".csv" => "_NT0.csv") end

println("[Process]")

if factorEstimateMode
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
        print(" exchange"); mdr.exchangeExpCurrency(erfile)
        print(" rebuild matrix"); mdr.buildExpenditureMatrix(year, substitute=codeSubst)
        println(" ... complete")
    end
    if PPPConv; print(" PPP converting:")
        mdr.convertToPPP(pppfile)
        println(" ... complete")
    end

    if year != base_year; print(" CPI scaling:")
        print(" read CPIs"); mdr.readCPIs([2010, 2015], cpi_file, idx_sep = ',', freq="A", unit="INX_A_AVG", topLev = "EU")
        print(", scaling"); mdr.scalingByCPI(year, base_year, codeDepth=0, topLev = "EU", subst = codeSubst)
        print(", rebuild matrix"); mdr.buildExpenditureMatrix(year, substitute=codeSubst)
        println(" ... complete")
    end

    print(" Concordance matrix building:")
    print(" concordance"); cmb.readXlsxData(concFiles[year], nation, nat_label = natLabels[year])
    print(", matrix"); cmb.buildConMat()
    print(", substitution"); cmb.addSubstSec(year, mdr.heSubst, mdr.heRplCd, mdr.heCats, exp_table = [])
    print(", normalization"); cmn = cmb.normConMat()   # {a3, conMat}
    println(" ... complete")

    print(" MRIO table reading:")
    if eoraRevised; eora_index = "../Eora/data/index/revised/" else eora_index = "../Eora/data/index/Eora_index.xlsx" end
    path = "../Eora/data/" * string(year) * "/" * string(year)
    print(" index"); ee.readIndexXlsx("../Eora/data/index/revised/", revised = eoraRevised)
    print(", IO table"); ee.readIOTables(year, path*"_eora_t.csv", path*"_eora_v.csv", path*"_eora_y.csv", path*"_eora_q.csv")
    print(", rearrange"); ee.rearrangeIndex(qmode=Qtable); ee.rearrangeTables(year, qmode=Qtable)
    println(" ... complete")

    print(" Data import:")
    print(" sector"); ee.getSectorData(year, mdr.heCodes, mdr.heSubst)
    print(", assemble"); ee.assembleConcMat(year, cmn)
    print(", category"); ec.readCategoryData(categoryFile, year, nutsLv, except=["None"], subCategory=subcat, nuts3pop=true)
    if length(subcat)==0; ec.setCategory(categories); else ec.setCategory(foodCategories); subcat *= "_" end

    print(", household"); ec.readHouseholdData(hhsfile, period = "daily", remove = true)
    print(", nuts weight"); ec.calculateNutsPopulationWeight()

    DE_files = filter(f->startswith(f, string(year)) && endswith(f, "_hhs_"*scaleTag*"DE.txt"), readdir(emissDataPath))
    de_nations = [f[6:7] for f in DE_files]
    DE_files = emissDataPath .* DE_files
    print(", DE"); ec.readEmissionData(year, de_nations, DE_files, mode = "de")
    print(", importing"); ed.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, nations = [])
    print(", detect NUTS"); ed.storeNUTS(year, cat_data = ec)
    print(", nuts weight"); ed.storeNutsWeight()
    println(" ... completed")

    print(" SDA decomposition:")
    print(" concordance")
    for n in ed.nat_list
        conc_mat_wgh = ee.buildWeightedConcMat(year, ee.abb[mdr.nationNames[n]])[1][:,:]
        ed.storeConcMat(year, n, conc_mat_wgh)
    end
    print(", decomposing"); ed.decomposeFactors(year, base_year, visible = true)
    if factorPrintMode; print(", factor printing"); ed.printFactors(factorPath) end
    println(" ... completed")
end

if SDA_mode
    print(" Factors reading:")
    print(" category"); ec.readCategoryData(categoryFile, year, nutsLv, except=["None"], subCategory=subcat, nuts3pop=true)
    print(", import"); ed.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, nations = [])
    print(", detect nation"); ed.detectNations(factorPath, year, base_year, factor_file_tag = "_factors.txt")
    print(", detect NUTS"); ed.storeNUTS(year, cat_data = ec)
    print(", ", year); ed.readFactors(year, factorPath, visible = false)
    print(", ", base_year); ed.readFactors(base_year, factorPath, visible = false)
    println(" ... completed")

    n_factor = 5
    delta_file = sda_path * "deltas.txt"
    print(" SDA process:")
    print(" preparing delta, "); ed.prepareDeltaFactors(year, base_year)
    print(" structural analysis"); for n in ed.nat_list; ed.structuralAnalysis(year, base_year, n, n_factor) end
    print(", print delta"); ed.printDelta(delta_file)
    println(" ... completed")
end
