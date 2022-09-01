# Developed date: 13. Apr. 2021
# Last modified date: 1. Sep. 2022
# Subject: Estimate carbon footprint by household consumptions
# Description: Calculate direct and indirect carbon emissions
#              by linking household consumptions and global supply chain,
#              utilizing Eora T, V, Y, and Q tables.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()

cd(Base.source_dir())

include("MicroDataReader.jl")
include("ConcMatBuilder.jl")
include("EmissionEstimator.jl")

using .MicroDataReader
using .ConcMatBuilder
using .EmissionEstimator

mdr = MicroDataReader
cmb = ConcMatBuilder
ee = EmissionEstimator

cesYear = 2016; exchYear = cesYear
eoraYear = cesYear
nation = "Vietnam"
natA3 = "VNM"
natCurr = "VND"
curr_unit= "USD"
emiss_unit = "tCO2"
fitEoraYear = false     # scaling micro-data's expenditure to fit the Eora target year
readMembers = false     # read member data
buildMatrix = false     # read expenditure data and build matrix
buildIeConc = false       # build Eora-CES concordance matrix
buildDeConc = false     # build direct emission concordance matrix

# cesYear = 2011; exchYear = cesYear
# eoraYear = cesYear
# nation = "India"
# natA3 = "IND"
# natCurr = "INR"
# curr_unit= "USD"
# emiss_unit = "tCO2"
# fitEoraYear = false     # scaling micro-data's expenditure to fit the Eora target year
# readMembers = false     # read member data
# buildMatrix = true      # read expenditure data and build matrix
# buildIeConc = true      # build Eora-CES concordance matrix
# buildDeConc = true      # build direct emission concordance matrix


# cesYear = 2018; exchYear = cesYear
# eoraYear = 2015
# nation = "Indonesia"
# natA3 = "IDN"
# natCurr = "IDR"
# curr_unit= "USD"
# emiss_unit = "tCO2"
# fitEoraYear = true      # scaling micro-data's expenditure to fit the Eora target year
# readMembers = false     # read member data
# buildMatrix = true      # read expenditure data and build matrix
# buildIeConc = true      # build Eora-CES concordance matrix
# buildDeConc = true      # build direct emission concordance matrix

filePath = Base.source_dir() * "/data/" * natA3 * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
emissionPath = filePath * "emission/"
commonIndexPath = Base.source_dir() * "/data/Common/"
concordancePath = filePath * "concordance/"
eoraIndexPath = commonIndexPath * "Eora/"
deDataPath = commonIndexPath* "DE/"

curConv = true; curr_target = "USD"; erfile = commonIndexPath * "CurrencyExchangeRates.txt"
pppConv = false; pppfile = filePath * "PPP_ConvertingRates.txt"

IE_mode = true             # indirect carbon emission estimation
DE_mode = false              # direct carbon emission estimation
DE_factor_estimate = false   # [true] estimate DE factors from IEA datasets, [false] read DE factors

Qtable = "I_CHG_CO2"; q_tag = "_i_chg_co2"
# Qtable = "PRIMAP"l q_tag = "_primap"

scaleMode = false
sparseMode = false
enhanceMode = false
fullMode = true
quantMode = false

memorySecure = false

if scaleMode; scaleTag = "_Scaled" else scaleTag = "" end

natFileTag = natA3 * "_" * string(cesYear)
regInfoFile = filePath * natFileTag * "_MD_RegionInfo.txt"
cmmfile = filePath * natFileTag * "_MD_Commodities.txt"
hhsfile = filePath * natFileTag * "_MD_Households_"*natCurr*".txt"
mmsfile = filePath * natFileTag * "_MD_Members.txt"
exmfile = filePath * natFileTag * "_MD_ExpenditureMatrix_"*natCurr*".txt"
erfile = filePath * natFileTag * "_MD_ExchangeRate.txt"

conmatEoraFile = filePath * natFileTag * "_IOT_ConcMatEora.txt"
conmatDeFile = filePath * natFileTag * "_IOT_ConcMatDe.txt"

expfile = filePath * natFileTag * "_MD_Expenditure_"*natCurr*".txt"

de_sec_file = deDataPath * "DE_sectors.txt"
de_conv_file = commonIndexPath * "Emission_converting_rate.txt"

# regInfoFile = extractedPath * natA3 * "_" * string(cesYear) * "_RegionInfo.txt"
# cmmfile = extractedPath * natA3 * "_" * string(cesYear) * "_Commodities.txt"
# hhsfile = extractedPath * natA3 * "_" * string(cesYear) * "_Households_"*natCurr*".txt"
# mmsfile = extractedPath * natA3 * "_" * string(cesYear) * "_Members.txt"
# itemfile = indexFilePath * natA3 * "_" * string(cesYear) * "_Commodity_items.txt"
# expfile = extractedPath * natA3 * "_" * string(cesYear) * "_Expenditure_"*natCurr*".txt"
# exmfile = extractedPath * natA3 * "_" * string(cesYear) * scaleTag * "_Expenditure_matrix_"*natCurr*".txt"

println("[Process]")

print(" Micro-data reading:")
print(" regions"); mdr.readPrintedRegionData(cesYear, natA3, regInfoFile)
print(", households"); mdr.readPrintedHouseholdData(cesYear, natA3, hhsfile)
if readMembers; print(", members"); mdr.readPrintedMemberData(cesYear, natA3, mmsfile) end
print(", sectors"); mdr.readPrintedSectorData(cesYear, natA3, cmmfile)
if buildMatrix
    print(", expenditures"); mdr.readPrintedExpenditureData(cesYear, natA3, expfile, quantity=quantMode)
    print(", matrix building"); mdr.buildExpenditureMatrix(cesYear, natA3, period = 365, quantity = quantMode)
else print(", expenditure matrix"); mdr.readPrintedExpenditureMatrix(cesYear, natA3, exmfile)
end
if fitEoraYear && eoraYear != nothing && eoraYear != cesYear
    print(", scaling from $cesYear to $eoraYear")
    exchYear = eoraYear
    cpiSecFile = indexFilePath * "CPI/CPI_"*natA3*"_sectors.txt"
    statFile = indexFilePath * "CPI/CPI_"*natA3*"_values.txt"
    linkFile = indexFilePath * "CPI/CPI_"*natA3*"_link.txt"
    mdr.scalingExpByCPI(cesYear, natA3, cpiSecFile, statFile, linkFile, eoraYear, period="year", region="district", revHH=buildMatrix, revMat=!buildMatrix)
end
if curConv; print(", currency exchange"); mdr.exchangeExpCurrency(cesYear, exchYear, natA3, natCurr, erfile, target_curr=curr_target, exp_mat=true) end
if pppConv; print(", ppp converting"); mdr.convertAvgExpToPPP(eoraYear, natA3, pppfile); println("complete") end
println(" ... completed")

print(" Concordance matrix building:")
print(" commodity_code"); cmb.getCommodityCodes(mdr.sc_list[cesYear][natA3])
if IE_mode
    nationFile = eoraIndexPath * "Eora_nations.txt"
    sectorFile = eoraIndexPath * "Eora_sectors.txt"
    concFile = concordancePath * natA3 * "_" * string(cesYear) * "_LinkedSectors_IE.txt"
    conMatFile = concordancePath * "ConcMat.txt"
    conSumMatFile = concordancePath * "ConcSumMat.txt"

    print(", IE data reading"); cmb.readIeSectors(nationFile, sectorFile)
    if buildIeConc
        print(", IE links reading"); cmb.readIeConcMatFile(concFile, weight=false)
        print(", IE matrix builing"); cmb.buildIeConMat()
    else
        print(", IE matrix reading"); cmb.readPrintedIeConMat(conmatEoraFile)
    end
    print(", normalization"); cmn_ie = cmb.normConMat()   # {a3, conMat}
    print(", print matrix"); cmb.printConMat(conMatFile, natA3, norm = true, categ = true)
    cmb.printSumNat(conSumMatFile, natA3, norm = true)
end
if DE_mode
    print(", DE data reading")
    ee.setNationDict(Dict("IND" =>"India", "IDN" => "Indonesia", "VNM" => "Viet Nam"))
    ee.readDirectEmissionData(cesYear, natA3, deDataPath, output_path = filePath * "de/", output_tag = natA3, integrate = true, cpi_scaling = false, cpi_base = 0, cpi_vals = [])
    if DE_factor_estimate
        print(", estimation")
        price_file = filePath * "de/" * "Price_" * natA3 * "_" * string(cesYear) * curr_unit * ".txt"
        emi_intens_file = filePath * "de/" * "Emission_intensity_" * natA3 * "_" * string(cesYear) * "_tCO2_per_" * curr_unit* ".txt"
        if curr_unit != "USD"; ee.exchangeEmCurrency(cesYear, erfile, target = curr_unit, origin = "USD", output = price_file) end
        ee.calculateEmissionRates(cesYear, output = emi_intens_file, currency = curr_unit)
        ee.printEmissionConvRates(cesYear, de_conv_file, emit_unit = emiss_unit, curr_unit = curr_unit)
    end
    print(", intensity")
    if natA3 == "VNM"; de_conv_file = replace(de_conv_file, ".txt" => "_VNMrevised.txt") end
    ee.readEmissionIntensity(cesYear, natA3, de_sec_file, de_conv_file, emit_unit = emiss_unit, curr_unit = curr_unit)
end
println(" ... complete")

if IE_mode
    print(" MRIO table reading:")
    eora_index = "../Eora/data/index/"
    path = "../Eora/data/" * string(eoraYear) * "/" * string(eoraYear)
    print(" index"); ee.readIndex(eora_index)
    print(", table"); ee.readIOTables(eoraYear, path*"_eora_t.csv", path*"_eora_v.csv", path*"_eora_y.csv", path*"_eora_q.csv")
    print(", rearranging"); ee.rearrangeMRIOtables(eoraYear, qmode=Qtable)
    print(", Leontief matrix"); ee.calculateLeontief(eoraYear)
    println(" ... complete")
end

print(" Emission calculation: ")
print("data"); ee.getDomesticData(cesYear, natA3, mdr.hh_list[cesYear][natA3], mdr.sc_list[cesYear][natA3], mdr.expMatrix[cesYear][natA3], (quantMode ? mdr.qntMatrix[cesYear][natA3] : []))
if memorySecure; print(", clear"); mdr.initVars() end
if DE_mode
    de_conc_file = concordancePath * natA3 * "_" * string(cesYear) * "_LinkedSectors_DE.txt"
    de_conc_mat_file = concordancePath * natA3 * "_" * string(cesYear) * "_ConcMat_DE.txt"
    deFile = emissionPath * "/" * string(cesYear) * "/" * string(cesYear) * "_" * natA3 * "_hhs_"*scaleTag*"DE.txt"
    print(", concordance_DE")
    if buildDeConc
        ee.buildDeConcMat(cesYear, natA3, de_conc_file, norm = true, output = de_conc_mat_file, energy_wgh = true)
    else
        ee.readDeConcMat(cesYear, natA3, conmatDeFile, norm = true, output = de_conc_mat_file, energy_wgh = true)
    end
    print(", estimate_DE"); ee.calculateDirectEmission(cesYear, natA3, quantity=quantMode, sparseMat=sparseMode, enhance=enhanceMode, full=fullMode)
    print(", print_DE"); ee.printEmissions(cesYear, natA3, deFile, mode = "de")
end
if IE_mode
    ieFile = emissionPath * "/" * string(cesYear) * "/" * string(cesYear) * "_" * natA3 * "_hhs_"*scaleTag*"IE"*q_tag*".txt"
    conWghMatFile = ""
    conSumMatFile = filePath * "concordance/" * natFileTag * "_ConcSumMat.txt"
    # conWghMatFile = filePath*"index/concordance/"*natA3*"_Eora"*scaleTag*"_weighted_concordance_table.csv"
    print(", concordance_IE"); ee.buildWeightedConcMat(cesYear, eoraYear, natA3, con_mat = cmn_ie, output = conWghMatFile)
    print(", estimate_IE"); ee.calculateIndirectEmission(cesYear, eoraYear, natA3, sparseMat=sparseMode, enhance=enhanceMode, full=fullMode, elapChk=1)
    print(", print_IE"); ee.printEmissions(cesYear, natA3, ieFile, mode = "ie")
end
println(" ... complete")

println("[all complete]")
