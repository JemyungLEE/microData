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

target_year = 2015
base_year = 2010
years = [base_year, target_year]

filePath = Base.source_dir() * "/data/"
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

categoryFile = indexFilePath * "Eurostat_Index_ver4.3_NT0.xlsx"
eustatsFile = indexFilePath * "EU_exp_COICOP.tsv"
cpi_file = indexFilePath * "EU_hicp.tsv"

concFiles = Dict(2010 => indexFilePath*"2010_EU_EORA_Conc_ver1.5.xlsx", 2015 => indexFilePath*"2015_EU_EORA_Conc_ver1.1.xlsx")
natLabels = Dict(2010 => "Eurostat", 2015 => "Eurostat_2015")

CurrencyConv = true; erfile = indexFilePath * "EUR_USD_ExchangeRates.txt"
PPPConv = false; pppfile = indexFilePath * "PPP_ConvertingRates.txt"

convertToBaseYear = true

codeSubst = true        # recommend 'false' for depth '1st' as there is nothing to substitute
perCap = true

eoraRevised = true

catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]
if codeSubst; substTag = "_subst" else substTag = "" end

ie_file_tag = "_hhs_"*scaleTag*"IE_"*Qtable*".txt"
de_file_tag = "_hhs_"*scaleTag*"DE.txt"

for year in years

    global filePath, indexFilePath, microDataPath, extractedPath, emissDataPath
    global Qtable, scaleMode, scaleTag, nation, nutsLv, categories, subcat
    global categoryFile, eustatsFile, cpi_file, concFiles, natLabels
    global CurrencyConv, erfile, PPPConv, pppfile, codeSubst, perCap, eoraRevised
    global catDepth, depthTag, codeSubst, substTag, ie_file_tag, de_file_tag

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
    println(" ... complete")

    print(" Data import:")
    print(" sector"); ee.getSectorData(year, mdr.heCodes, mdr.heSubst)
    print(", category"); ec.readCategoryData(categoryFile, year, nutsLv, except=["None"], subCategory=subcat, nuts3pop=true)
    ec.setCategory(categories)
    print(", household"); ec.readHouseholdData(hhsfile, period = "daily", remove = true)
    print(", nuts weight"); ec.calculateNutsPopulationWeight(year = year)

    cf_file_path = emissDataPath * string(year) * "/"
    if convertToBaseYear && year != base_year; ie_ftag = replace(ie_file_tag, ".txt" => "_converted_" * string(base_year) * "_price.txt")
    else ie_ftag = ie_file_tag
    end
    IE_files = filter(f->startswith(f, string(year)) && endswith(f, ie_ftag), readdir(cf_file_path))
    DE_files = filter(f->startswith(f, string(year)) && endswith(f, de_file_tag), readdir(cf_file_path))
    ie_nations = [f[6:7] for f in IE_files]
    de_nations = [f[6:7] for f in DE_files]
    IE_files = cf_file_path .* IE_files
    DE_files = cf_file_path .* DE_files

    print(", IE"); ec.readEmissionData(year, ie_nations, IE_files, mode = "ie")
    print(", DE"); ec.readEmissionData(year, de_nations, DE_files, mode = "de")
    print(", importing"); ed.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, nations = [])

    print(", detect NUTS"); ed.storeNUTS(year, cat_data = ec)
    print(", nuts weight"); ed.storeNutsWeight(year = year)
    println(" ... completed")
end

sda_path = emissDataPath * "SDA/"
pop_label = Dict(0 => "", 1 => "dense", 2 => "inter", 3 => "sparse")
pop_dens = 0        # [1] Densely populated, [2] Intermediate, [3] Sparsely populated
ci_rste = 0.95
n_iter = 10000

CI_test = false; test_nats = ["BE"];
if CI_test; test_tag = "_test" else test_tag = "" end
ci_file = sda_path * string(target_year) * "_" * string(base_year) * "_CI_" * pop_label[pop_dens] * test_tag * ".txt"
if convertToBaseYear; ci_file = replace(ci_file, ".txt" => "_converted_" * string(base_year) * "_price.txt") end

println(" Bootstrap process:")

if CI_test; nats = test_nats else nats = ed.filterNations() end

if pop_dens in [1, 2, 3]
    print(" read gridded population"); ed.getPopGridded(years, categoryFile)
    nat_flt = ed.filterNonPopDens(years, nats, pop_dens = pop_dens)
    print(", filter nation"); filter!(x -> x in nat_flt, nats)
end

print(", bootstrap")
for y in years
    ed.estimateConfidenceIntervals(y, nats, iter = n_iter, ci_rate = ci_rste, resample_size = 0, replacement = true, pop_dens = pop_dens, visible = true)
end
print(", printing resutls"); ed.printConfidenceIntervals(years, ci_file, nats, pop_dens = pop_dens, ci_rate = ci_rste)

println(" ... completed")
