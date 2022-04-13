# Developed date: 9. Nov. 2021
# Last modified date: 12. Nov. 2021
# Subject: Categorized emission gap mapping
# Description: Mapping emission gaps through households emissions data, categorizing by district, income-level, and etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())

include("MicroDataReader.jl")
include("EmissionEstimator.jl")
include("EmissionCategorizer.jl")
include("EmissionDecomposer.jl")
include("../GIS/QgisStyleExporter.jl")
using .MicroDataReader
using .EmissionEstimator
using .EmissionCategorizer
using .EmissionDecomposer
using .QgisStyleExporter
mdr = MicroDataReader
ee = EmissionEstimator
ec = EmissionCategorizer
ed = EmissionDecomposer
qse = QgisStyleExporter

println("[Process]")

nation = "Eurostat"
target_year, base_year = 2015, 2010
int_year = [target_year + 10000 * base_year]
nutsLv = 1
onlyNutsInHbs = true
Qtable = "_PRIMAP"
ceIntegrateMode = "cf"      # "ie" (only indirect CE), "de" (only direct CE), or "cf" (integrage direct and indirect CEs)
ceProcessMode = ["ie", "de", "cf"]
cePrintMode = ["ie", "de", "cf"]
ceExportMode = "cf"         # "ie" (only indirect CE), "de" (only direct CE), or "cf" (integrage direct and indirect CEs

cpi_scaling = true

substMode = true; if substMode; substTag = "_subst" else substTag = "" end
scaleMode = true; if scaleMode; scaleTag = "Scaled_" else scaleTag = "" end

eqvalMode = false   # [true]: apply square root of household size for equivalance scale
ntWeighMode = true  # [true]: apply NUTS population based weight, [false]:apply HBS weight

minmaxv = [[[-1.5*10^8,1.5*10^8]], [[-8.0,8.0]]] # {{overall CF min., max.}, {CF per capita min., max.}
expNtMode = "hbs"
exportWebMode = true
mapStyleMode = true; colormapReversePerCap=false; labeRevPerCap=false; colormapReverse=false; labeRev=false

popweight = true
grid_pop = true

filePath = Base.source_dir() * "/data/"
indexPath = filePath * "index/"
extrPath = filePath * "extracted/"
indexFile = indexPath * "Eurostat_Index_ver5.0.xlsx"

ExpenditureFile = extrPath * scaleTag * "Expenditure_matrix_4th" * substTag * ".csv"

incomePeriod = "annual"  # Period: "annual"(default), "monthly", or "daily"

normTag = ["perCapNorm", "perHhNorm"]
categories = ["Food", "Electricity", "Gas", "Other energy", "Public transport", "Private transport", "Medical care",
                "Education", "Consumable goods", "Durable goods", "Other services", "Total"]
subcat=""

for year in [target_year, base_year]
    global emissPath, filePath, extrPath, indexFile, ExpenditureFile, nutsLv, subcat, base_year
    global onlyNutsInHbs, Qtable, ceIntegrateMode, ceProcessMode, cePrintMode, ceExportMode, cpi_scaling
    global substMode, substTag, scaleMode, scaleTag, eqvalMode, ntWeighMode
    global minmaxv, expNtMode, exportWebMode, mapStyleMode, colormapReversePerCap, labeRevPerCap, colormapReverse, labeRev
    global popweight, grid_pop, incomePeriod, normTag, categories

    println(year)
    emissPath = filePath * "emission/" * string(year) * "/"
    hhsfile = extrPath * string(year) * "_Households.csv"
    print(" Data reading: ")
    print(" hhs micro-data"); mdr.readPrintedHouseholdData(hhsfile)
    print(", category"); ec.readCategoryData(indexFile, year, nutsLv, except=["None"], subCategory=subcat)
    ec.setCategory(categories)
    print(", household"); ec.readHouseholdData(hhsfile, period = incomePeriod, remove = onlyNutsInHbs, alter=true)
    print(", population"); ec.readPopulation(year, indexFile, nuts_lv = nutsLv)
    print(", gridded population"); ec.readPopGridded(year, indexFile, nuts_lv = [nutsLv], adjust = true)
    print(", emission")
    IE_files = []; DE_files = []; ie_nations = []; de_nations = []
    ie_file_tag = "_hhs_"*scaleTag*"IE"*Qtable*".txt"
    if cpi_scaling && base_year != year; ie_file_tag = replace(ie_file_tag, ".txt" => "_converted_" * string(base_year) * ".txt") end
    de_file_tag = "_hhs_"*scaleTag*"DE.txt"
    for f in readdir(emissPath)
        if startswith(f, string(year)) && endswith(f, ie_file_tag); push!(IE_files, emissPath*f); push!(ie_nations, f[6:7])
        elseif startswith(f, string(year)) && endswith(f, de_file_tag); push!(DE_files, emissPath*f); push!(de_nations, f[6:7])
        end
    end
    print("_IE"); ec.readEmissionData(year, ie_nations, IE_files, mode = "ie")
    print("_DE"); ec.readEmissionData(year, de_nations, DE_files, mode = "de")
    print("_CF"); ec.integrateCarbonFootprint(year, mode=ceIntegrateMode)
    println(" ... complete")

    print(" Weights calculating: ")
    ec.calculateNutsPopulationWeight(year = year, pop_dens = grid_pop, adjust = true)
    println(" ... complete")

    print(" Categorizing:")
    print(" category")
    for m in ceProcessMode
        print("_",uppercase(m))
        ec.categorizeHouseholdEmission(year, mode=m, output="", hhsinfo=false, nutsLv=1)
        ec.categorizeRegionalEmission(year, mode=m, nutsLv=1, period=incomePeriod, adjust=true, religion=false, popWgh=popweight, ntweigh=ntWeighMode)
    end

    print(", importing"); ed.importData(hh_data = mdr, mrio_data = ee, cat_data = ec, nations = [], cat_filter = false)
    print(", convert NUTS"); ed.convertNUTS(year = year)
    print(", detect NUTS"); ed.storeNUTS(year, cat_data = ec)
    println(" ... complete")
end

# println("Print results: ")
# NutsEmissionFile = emissPath * string(year) * "_EU_nuts_"*scaleTag*subcat*"emission_cat.csv"
# for iy in int_year
#     print(iy, ":")
#     print(" National abstract")
#     ec.makeNationalSummary(iy, emissPath * string(iy) * "_National_summary_"*scaleTag*Qtable*".txt", nuts_mode=true)
#     print(", regional emission")
#     ec.printRegionalEmission(iy, NutsEmissionFile, mode=cePrintMode, totm=true, expm=true, popm=true, relm=false, wghm=true, povm=false, ntweigh=ntWeighMode)
#     println(" ... complete")
# end

println(target_year, "-", base_year)
print(" Gap exporting:")
nats = ed.filterNations()
emissPath = filePath * "emission/" * string(int_year[1]) * "/"
mkpath(emissPath)
print(" NUTS integration"); ed.integrateNUTS(target_year, base_year, indexFile, modify = true, pop_dens = true)
print(", importing"); ec.importIntegratedNUTS(ed.nuts_intg, ed.nuts_intg_list)
gisTag = "NUTS"
exportFile = emissPath * "YEAR_EU_NUTS_gis_"*scaleTag*subcat*"emission_cat_OvPcTag.csv"
exportRateFile = emissPath * "YEAR_EU_NUTS_gis_"*scaleTag*subcat*"emission_cat_dr_OvPcTag.csv"
print(", region"); ec.exportRegionalEmission([target_year, base_year], gisTag, "", mode=ceExportMode, nutsmode=expNtMode, nspan=128, minmax=minmaxv, descend=true, empty=false, logarithm=false)
print(", gap"); labelList, labelListPerCap = ec.calculateTemporalGap(target_year, base_year, exportFile, nats, mode="cf", nspan=128, minmax=minmaxv, descend=true, logarithm=false, tag="NUTS")
print(", diff"); spanVals, spanValsPerCap = ec.exportEmissionDiffRate(int_year, gisTag, exportRateFile, 0.5, -0.5, 128, descend=false, empty=false)

if exportWebMode; print(", web-files")
    exportPath = Base.source_dir() * "/data/emission/"*scaleTag*"webfile/"
    ec.exportWebsiteFiles(int_year, exportPath, nutsmode=expNtMode, rank=true, empty=false, major=true)
end

if mapStyleMode; print(", map-style file generating")
    rgbFile = "../GIS/data/EU/MPL_RdBu.rgb"
    for i=1:length(ec.catList)
        qse.readColorMap(rgbFile, reverse=colormapReversePerCap)
        qmlFile = replace(rgbFile, ".rgb"=>"_percap_"*scaleTag*ec.catList[i]*".qml")
        attr = string(int_year[1])*"_EU_NUTS_gis_"*subcat*"emission_cat_dr_percap_gr_"*ec.catList[i]
        qse.makeQML(qmlFile, attr, empty=false, values=spanValsPerCap[int_year[1]][:,i], indexValue=true, labelReverse=labeRevPerCap)
    end

    rgbFile = "../GIS/data/EU/MPL_YlGnBu.rgb"
    for i=1:length(ec.catList)
        qse.readColorMap(rgbFile, reverse=colormapReverse)
        qmlFile = replace(rgbFile, ".rgb"=>"_overall_"*scaleTag*ec.catList[i]*".qml")
        attr = string(int_year[1])*"_EU_NUTS_gis_"*subcat*"emission_cat_overall_gr_"*ec.catList[i]
        qse.makeQML(qmlFile, attr, empty=false, labels=labelList[int_year[1]][:,i], indexValue=true, labelReverse=labeRev)
    end
end
println(" ... complete")

println("\n[Done]")
