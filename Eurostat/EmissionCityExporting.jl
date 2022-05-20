# Developed date: 22. Dec. 2021
# Last modified date: 20. May. 2022
# Subject: Exporting City CF and CI web-files
# Description: Export CF and CI data by category for each city
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

base_year, target_year = 2010, 2015
years = [base_year, target_year]

filePath = Base.source_dir() * "/data/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissDataPath = filePath* "emission/"
webFilePath = indexFilePath * "web/"
mrioPath = "../Eora/data/"

Qtable = "PRIMAP"
scaleMode = true; if scaleMode; scaleTag = "Scaled_" else scaleTag = "" end

nation = "Eurostat"
nutsLv = 1

categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
web_categories = ["FOOD", "ELECTRICITY", "GAS", "ENERGY", "PUBLIC_TRANS", "PRIVATE_TRANS", "MEDICAL",
                "EDUCATION", "CONSUMABLE", "DURABLE", "SERVICES", "ALL"]
subcat=""

categoryFile = indexFilePath * "Eurostat_Index_ver5.0.xlsx"
eustatsFile = indexFilePath * "EU_exp_COICOP.tsv"
cpi_file = indexFilePath * "EU_hicp.tsv"

web_index_file = webFilePath * "keycode_index.txt"

concFiles = Dict(2010 => indexFilePath*"2010_EU_EORA_Conc_ver1.5.xlsx", 2015 => indexFilePath*"2015_EU_EORA_Conc_ver1.1.xlsx")
natLabels = Dict(2010 => "Eurostat", 2015 => "Eurostat_2015")

CurrencyConv = true; erfile = indexFilePath * "EUR_USD_ExchangeRates.txt"
PPPConv = false; pppfile = indexFilePath * "PPP_ConvertingRates.txt"

cfav_file, cfac_file = Dict{Int, String}(), Dict{Int, String}()
for y in years
    cfav_file[y] = emissDataPath * string(y) * "/" * string(y) *  "_EU_NUTS_gis_Scaled_emission_cat_overall_CF_gr.csv"
    cfac_file[y] = emissDataPath * string(y) * "/" * string(y) *  "_EU_NUTS_gis_Scaled_emission_cat_dr_percap_gr.csv"
end

codeSubst = true        # recommend 'false' for depth '1st' as there is nothing to substitute
perCap = true
grid_pop = true

adjustConc = false
domestic_mode = false

removeNTZ = true
adjustNTZ = false

incomePeriod = "annual"

catDepth = 4
depthTag = ["1st", "2nd", "3rd", "4th"]
if codeSubst; substTag = "_subst" else substTag = "" end

ce_intgr_mode = "cf"

ie_file_tag = "_hhs_"*scaleTag*"IE_"*Qtable*".txt"
de_file_tag = "_hhs_"*scaleTag*"DE.txt"

city_file_sector = Array{Tuple{String, String}, 1}()
f = open(web_index_file)
for l in eachline(f)
    s = string.(split(l, '\t'))
    push!(city_file_sector, (s[1], s[2]))
end
close(f)

for year in years

    global filePath, indexFilePath, microDataPath, extractedPath, emissDataPath
    global Qtable, scaleMode, scaleTag, nation, nutsLv, categories, subcat, adjustConc, domestic_mode
    global categoryFile, eustatsFile, cpi_file, concFiles, natLabels
    global CurrencyConv, erfile, PPPConv, pppfile, codeSubst, perCap
    global catDepth, depthTag, codeSubst, substTag, grid_pop

    println("[",year,"]")
    microDataPath *= string(year) * "/"
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
    print(" household")
    ec.readHouseholdData(hhsfile, period = incomePeriod, alter=true, remove_nt = true, remove_z = removeNTZ)
    mdr.readPrintedHouseholdData(hhsfile, regions = ec.regList[year])
    if codeSubst; print(", subst code"); mdr.readSubstCodesCSV(year, sbctgfile, sbcdsfile) end
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

    print(" Data import:")
    print(" sector"); ee.getSectorData(year, mdr.heCodes, mdr.heSubst)
    print(", population"); ec.readPopulation(year, categoryFile, nuts_lv = nutsLv)
    print(", gridded population"); ec.readPopGridded(year, categoryFile, nuts_lv = [nutsLv], adjust = true)

    cf_file_path = emissDataPath * string(year) * "/"
    if year != base_year; ie_ftag = replace(ie_file_tag, ".txt" => "_converted_" * string(base_year) * ".txt")
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
    print(", CF"); ec.integrateCarbonFootprint(year, mode = ce_intgr_mode)
    print(", nuts weight"); ec.calculateNutsPopulationWeight(year = year, pop_dens = true, adjust = adjustNTZ)
    print(", categorizing hhs"); ec.categorizeHouseholdEmission(year, mode = ce_intgr_mode, output="", hhsinfo=false, nutsLv=1)
    print(" reg"); ec.categorizeRegionalEmission(year, mode=ce_intgr_mode, nutsLv=1, period=incomePeriod, adjust=adjustNTZ, religion=false, popWgh=true, ntweigh=true)
    print(", importing"); ed.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, nations = [])
    print(", convert NUTS"); ed.convertNUTS(year = year)
    print(", detect NUTS"); ed.storeNUTS(year, cat_data = ec)
    print(", store weight"); ed.storeNutsWeight(year = year)
    println(" ... completed")
end

web_path = emissDataPath * "web/"
ci_rste = 0.95
n_iter = 10000

CI_test = false; test_nats = ["BE","ES","FR"]

println(" Bootstrap process:")

if CI_test; nats = test_nats; web_path *= "test/" else nats = [] end


print(", bootstrap")
for y in years
    ed.estimateConfidenceIntervals(y, nats, iter = n_iter, ci_rate = ci_rste, resample_size = 0, replacement = true, visible = true, adjust = adjustNTZ)
end
print(", web-file exporting")
ed.exportWebsiteCityFiles(years, nats, web_path, web_categories, city_file_sector, cfav_file, cfac_file)

println(" ... completed")
