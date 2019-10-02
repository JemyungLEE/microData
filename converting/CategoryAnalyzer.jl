# Developed date: 27. Aug. 2019
# Last modified date: 1. Oct. 2019
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

nation = "India"

inputFile = "India_Eora_singleLink.txt"
#inputFile = "India_Eora_multipleLinks.txt"
#inputFile = "India(Eora)_Basic(Eora).txt"
inputXLSX = "India(STAT) vs EORA_Ver1.1.xlsx"

conMatFile = "ConcordanceMatrix_"*inputFile
inputFile = Base.source_dir()*"/"*inputFile
conMatFile = Base.source_dir()*"/"*conMatFile
inputXLSX = Base.source_dir()*"/"*inputXLSX

xe.readXlsxFile(inputXLSX, nation)
#cmb.readXlsxFile(inputXLSX)
#cmb.readClassCodes(inputFile)
#cmb.readClassCodes(inputFile, "Basic", "India")
#cmb.makeConMat()
#cmb.printConMat(conMatFile)
