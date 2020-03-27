# Developed date: 27. Mar. 2020
# Last modified date: 27. Mar. 2020
# Subject: Plotting charts
# Description: Plot charts
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("BoxChart.jl")
using .BoxChart
bc = BoxChart

println("[Plotting]")

emissionData = "../India/data/emission/2011_IND_hhs_emission_rng_perCap.csv"
expendirtureData = "../India/data/emission/2011_IND_hhs_emission_rng_perCap_exp.csv"

xtitle = "Daily expenditure (USD/day/capita)"
ytitle = "Annual carbon emission (tCO2/year/capita)"

println(" Values reading")
bc.readValues(expendirtureData, emissionData)
println(" Plot printing")
bc.plotBoxChart(xtitle, ytitle, disp=true, dash=true, output=Base.source_dir()*"/chart/BoxChart.png")

println("[Done]")
