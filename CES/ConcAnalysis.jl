# Developed date: 12. May. 2021
# Last modified date: 18. Mar. 2022
# Subject: Concordance data analysis
# Description: read concordance data and extract linked sectors, and build concordance matrix
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("ConcMatBuilder.jl")
using .ConcMatBuilder
cmb = ConcMatBuilder

nation = "Indonesia"
natA3 = "IDN"
natLabel = ""
year = 2018

# nation = "India"
# natA3 = "IND"
# natLabel = ""
# year = 2011

# nation = "Eurostat"
# natA3 = "EU"
# natLabel = ""
# year = 2010

# year = 2015
# nation = "Eurostat"
# natLabel = "Eurostat_" * string(year)
# natA3 = "EU"

xlsxMode = true
textMode = false

filePath = Base.source_dir() * "/data/" * natA3 * "/"
concordancePath = filePath * "index/concordance/"
conMatFile = concordancePath * string(year) * "_ConcMat.txt"
conSumMatFile = concordancePath * string(year) * "_ConcSumMat.txt"
nationFile = concordancePath * string(year) * "_Eora_nations.txt"
sectorFile = concordancePath * string(year) * "_Eora_sectors.txt"
concFile = concordancePath * string(year) * "_LinkedSectors_IE.txt"

if (year, nation) == (2018, "Indonesia"); xlsxFile = filePath * "index/IDN_EORA_Conc_ver1.1.xlsx"
elseif (year, nation) == (2011, "India"); xlsxFile = filePath * "index/IND_EORA_Conc_ver1.2.xlsx"
elseif (year, nation) == (2010, "Eurostat"); xlsxFile = "/Users/leejmacbook/github/microData/Eurostat/data/index/" * "2010_EU_EORA_Conc_ver1.5.xlsx"
elseif (year, nation) == (2015, "Eurostat"); xlsxFile = "/Users/leejmacbook/github/microData/Eurostat/data/index/" * "2015_EU_EORA_Conc_ver1.1.xlsx"
end

if xlsxMode
    print("Xlsx file mode: ")
    print(" reading"); cmb.readXlsxData(xlsxFile, nation, weight=false, nat_label = natLabel)
    print(", matrix"); cmb.buildIeConMat()
    print(", normalization"); cmb.normConMat()
    print(", exporting"); cmb.exportLinkedSectors(concFile, nation, mrio="Eora")
    print(", printint"); cmb.printConMat(conMatFile, natA3, norm = true, categ = true)
    cmb.printSumNat(conSumMatFile, natA3, norm = true)
    println(" complete")
end

if textMode
    print(" Concordance matrix building:")
    print(" data reading"); cmb.readIeConcMatFile(nationFile, sectorFile, concFile, weight=false)
    print(", matrix builing"); cmb.buildIeConMat()
    print(", normalization"); cmn = cmb.normConMat()   # {a3, conMat}
    print(", print matrix"); cmb.printConMat(conMatFile, natA3, norm = true, categ = true)
    # cmb.printSumNat(conSumMatFile, natA3, norm = true)
    println(" complete")
end
