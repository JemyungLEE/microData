# Developed date: 13. Feb. 2020
# Last modified date: 11. Mar. 2019
# Subject: GIS map style file generator
# Description: Make map style files including QML
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("QgisStyleExporter.jl")
using .QgisStyleExporter
qse = QgisStyleExporter

println("[Process]")
print(" RGB file reading: ")
#rgbFile = Base.source_dir()*"/data/MPL_RdBu.rgb"
rgbFile = Base.source_dir()*"/data/MPL_YlGnBu.rgb"
qse.readColorMap(rgbFile)
println("complete")

print(" QML file exporting: ")
attr = "2011_IND_hhs_GIS_emission_cat_dr_perCap_gr_Total"
qmlFile = replace(rgbFile, ".rgb"=>".qml")
qse.makeQML(qmlFile, attr, empty=true)
println("complete")

println("[Complete]")
