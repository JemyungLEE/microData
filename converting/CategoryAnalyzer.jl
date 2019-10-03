# Developed date: 27. Aug. 2019
# Last modified date: 3. Oct. 2019
# Subject: classification category analyzer, India-Eora
# Description: match categories of India commodity classificartion and Eora industry classification
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
#include("ConcMatBuilder.jl")
include("XLSXextractor.jl")
#using .ConcMatBuilder
#cmb = ConcMatBuilder
using .XLSXextractor
xe = XLSXextractor

nation = "India"    # converting nation

inputFile = "India_Eora_singleLink.txt"
#inputFile = "India_Eora_multipleLinks.txt"
#inputFile = "India(Eora)_Basic(Eora).txt"
inputXLSX = "India(STAT) vs EORA_Ver1.2.xlsx"
outputXLSX = "India_EORA_concordance.txt"
outputXLSXsum = "India_EORA_sumIndSec.txt"
outputXLSXnorm = "India_EORA_concordance_norm.txt"
outputXLSXsumNorm = "India_EORA_sumIndSec_norm.txt"

conMatFile = "ConcordanceMatrix_"*inputFile
inputFile = Base.source_dir()*"/"*inputFile
conMatFile = Base.source_dir()*"/"*conMatFile
inputXLSX = Base.source_dir()*"/"*inputXLSX

xe.readXlsxData(inputXLSX, nation)
xe.buildConMat()
xe.normConMat()
xe.printConMat(outputXLSX, nation)
xe.printSumNat(outputXLSXsum, nation)
xe.printConMat(outputXLSXnorm, nation, true)
xe.printSumNat(outputXLSXsumNorm, nation, true)

#cmb.readXlsxFile(inputXLSX)
#cmb.readClassCodes(inputFile)
#cmb.readClassCodes(inputFile, "Basic", "India")
#cmb.makeConMat()
#cmb.printConMat(conMatFile)

println("All processes completed.")
