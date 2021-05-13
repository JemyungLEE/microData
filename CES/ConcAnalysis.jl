# Developed date: 12. May. 2021
# Last modified date: 13. May. 2021
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
xlsxMode = true
textMode = false

filePath = Base.source_dir() * "/data/" * natA3 * "/"
concordancePath = filePath * "index/concordance/"
conMatFile = concordancePath * "ConcMat.txt"
conSumMatFile = concordancePath * "ConcSumMat.txt"
nationFile = concordancePath * "Eora_nations.txt"
sectorFile = concordancePath * "Eora_sectors.txt"
concFile = concordancePath * "LinkedSectors_IE.txt"
xlsxFile = filePath * "index/IDN_EORA_Conc_ver1.0.xlsx"

if xlsxMode
    print("Xlsx file mode: ")
    print(" reading"); cmb.readXlsxData(xlsxFile, nation, weight=false)
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
