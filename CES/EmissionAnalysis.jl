# Developed date: 13. Apr. 2021
# Last modified date: 14. Apr. 2021
# Subject: Estimate carbon footprint by household consumptions
# Description: Calculate direct and indirect carbon emissions
#              by linking household consumptions and global supply chain,
#              utilizing Eora T, V, Y, and Q tables.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)


cd(Base.source_dir())

include("MicroDataReader.jl")
include("ConcMatBuilder.jl")

using .MicroDataReader
using .ConcMatBuilder

mdr = MicroDataReader
cmb = ConcMatBuilder

filePath = Base.source_dir() * "/data/"

nation = "Indonesia"
natA3 = "IDN"

# Converting process of Eora final demand data to India micro-data format
concordanceFile = filePath * "/" * natA3 * "/index/"* natA3 *"_EORA_Conc_ver0.9.xlsx"
linkSecFile = filePath * "/" * natA3 * "/index/concordance/LinkedSectors_IE.txt"
conMatFile = filePath * "/" * natA3 * "/index/concordance/ConcMat.txt"
conSumMatFile = filePath * "/" * natA3 * "/index/concordance/ConcSumMat.txt"
print(" Concordance matrix building:")
# print(" xlsx reading"); cmb.readXlsxData(concordanceFile, nation, weight=false)

nationFile = filePath * "/" * natA3 * "/index/concordance/Eora_nations.txt"
sectorFile = filePath * "/" * natA3 * "/index/concordance/Eora_sectors.txt"
concFile = filePath * "/" * natA3 * "/index/concordance/LinkedSectors_IE.txt"
print(" concordance data reading"); cmb.readConcMatFile(nationFile, sectorFile, concFile, nation, weight=false)

print(", linkage printing"); cmb.exportLinkedSectors(linkSecFile, natA3, mrio="Eora")
print(", matrix builing"); cmb.buildConMat()
print(", normalization"); cmn = cmb.normConMat()   # {a3, conMat}

cmb.printConMat(conMatFile, nation, norm = true, categ = true)
cmb.printSumNat(conSumMatFile, nation, norm = true)
println(" complete")
