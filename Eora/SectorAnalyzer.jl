# Developed date: 2. Sep. 2019
# Last modified date: 4. Sep. 2019
# Subject: Analysis of sector classification
# Description: classify sector categories of all nations in Eora dataset.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
include("SectorClassifier.jl")
using .SectorClassifier
sc = SectorClassifier

inputFile = "indexSector.txt"

inputFile = Base.source_dir()*"/"*inputFile
abstractFile = replace(inputFile, ".txt" => "_abstract.txt")
comparedFile = replace(inputFile, ".txt" => "_compared.txt")
containedFile = replace(inputFile, ".txt" => "_contained.txt")

sc.readSectorData(inputFile)
sc.compareSectors(comparedFile, containedFile)
#sc.printNationData(abstractFile)
