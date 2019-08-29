# Developed date: 27. Aug. 2019
# Last modified date: 27. Aug. 2019
# Subject: classification category analyzer, India-Eora
# Description: match categories of India commodity classificartion and Eora industry classification
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
include("ConcMatBuilder.jl")
using .ConcMatBuilder
cmb = ConcMatBuilder

inputFile = "India_Eora_singleLink.txt"
#inputFile = "India_Eora_multipleLinks.txt"
conMatFile = "ConcordanceMatrix_"*inputFile
inputFile = Base.source_dir()*"/"*inputFile
conMatFile = Base.source_dir()*"/"*conMatFile

cmb.readClassCodes(inputFile)
cmb.makeConMat()
cmb.printConMat(conMatFile)
