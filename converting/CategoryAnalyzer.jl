# Developed date: 27. Aug. 2019
# Last modified date: 17. Oct. 2019
# Subject: classification category analyzer, India-Eora
# Description: match categories of India commodity classificartion and Eora industry classification
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
#include("ConcMatBuilder.jl")
#include("XLSXextractor.jl")
#using .ConcMatBuilder
#cmb = ConcMatBuilder
#using .XLSXextractor
#xe = XLSXextractor

nation = "India"    # converting nation
#=
inputFile = "India_Eora_singleLink.txt"
#inputFile = "India_Eora_multipleLinks.txt"
#inputFile = "India(Eora)_Basic(Eora).txt"
inputXLSX = "India(STAT) vs EORA_Ver1.2.xlsx"
outputXLSX = "India_EORA_concordance.txt"
outputXLSXsum = "India_EORA_sumIndSec.txt"
outputXLSXnorm = "India_EORA_concordance_norm.txt"
outputXLSXsumNorm = "India_EORA_sumIndSec_norm.txt"

conMatFile = "ConcordanceMatrix_"*inputFile
inputFile = Base.source_dir()*"/data/"*inputFile
conMatFile = Base.source_dir()*"/data/"*conMatFile
inputXLSX = Base.source_dir()*"/data/"*inputXLSX

xe.readXlsxData(inputXLSX, nation)
xe.buildConMat()
xe.normConMat()
xe.printConMat(outputXLSX, nation)
xe.printSumNat(outputXLSXsum, nation)
xe.printConMat(outputXLSXnorm, nation, true)
xe.printSumNat(outputXLSXsumNorm, nation, true)
=#
#cmb.readXlsxFile(inputXLSX)
#cmb.readClassCodes(inputFile)
#cmb.readClassCodes(inputFile, "Basic", "India")
#cmb.makeConMat()
#cmb.printConMat(conMatFile)

include("HsConcMatBuilder.jl")
using .HsConcMatBuilder
hb = HsConcMatBuilder

inputHS = "India-HS converting table_Ver1.3.xlsx"
outputHS = "India_HS_concordance.txt"
outputHSnorm = "India_HS_concordance_norm.txt"

inputHS = Base.source_dir()*"/data/"*inputHS
outputHS = Base.source_dir()*"/data/"*outputHS
outputHSnorm = Base.source_dir()*"/data/"*outputHSnorm

hb.readXlsxData(inputHS, nation, "India_HS_Lv6")
ct = hb.buildConMat()
#hb.printConMat(outputHS, ct, nation)
ctn = hb.normConMat(ct)
hb.printConMat(outputHSnorm, ctn, nation)

println("All processes completed.")
