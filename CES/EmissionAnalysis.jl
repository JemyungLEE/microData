# Developed date: 13. Apr. 2021
# Last modified date: 13. Apr. 2021
# Subject: Estimate carbon footprint by household consumptions
# Description: Calculate direct and indirect carbon emissions
#              by linking household consumptions and global supply chain,
#              utilizing Eora T, V, Y, and Q tables.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)


cd(Base.source_dir())

include("MicroDataReader.jl")
include("../converting/XLSXextractor.jl")

using .MicroDataReader
using .XLSXextractor

mdr = MicroDataReader
xls = XLSXextractor

filePath = Base.source_dir() * "/data/"

nation = "Indonesia"
natA3 = "IDN"

# Converting process of Eora final demand data to India micro-data format
concordanceFile = filePath * "/" * natA3 * "/index/"* natA3 *"_EORA_Conc_ver0.9.xlsx"
linkSecFile = filePath * "/" * natA3 * "/index/concordance/LinkedSectors.txt"
conMatFile = filePath * "/" * natA3 * "/index/concordance/ConcMat.txt"
conSumMatFile = filePath * "/" * natA3 * "/index/concordance/ConcSumMat.txt"
print(" Concordance matrix building:")
print(" xlsx reading"); xls.readXlsxData(concordanceFile, nation)
print(", linkage printing"); xls.exportLinkedSectors(linkSecFile, natA3)
print(", matrix builing"); xls.buildConMat()
print(", normalization"); cmn = xls.normConMat()   # {a3, conMat}

xls.printConMat(conMatFile, nation, norm = true, categ = true)
xls.printSumNat(conSumMatFile, nation, norm = true)
println(" complete")
