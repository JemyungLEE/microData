# Developed date: 27. Mar. 2020
# Last modified date: 27. Mar. 2020
# Subject: Plotting charts
# Description: Plot box and violin charts
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

cd(Base.source_dir())

include("BoxChart.jl")
using .BoxChart
bc = BoxChart
include("EmissionPlots.jl")
using .EmissionPlots
ep = EmissionPlots

boxPlotMode = false
stackedBarMode = true

println("[Plotting]")

if stackedBarMode
    print(" Stacked bar chart plotting: ")
    emissionFile = "../India/data/emission/2011_IND_hhs_emission_inc_perCap.csv"
    hhsFile = "../India/data/emission/2011_IND_hhs_emission_cat.csv"
    chartFile = Base.source_dir()*"/chart/Stacked_bar_chart.png"
    ep.readColorPalette("Table color palettes.txt", rev=true, tran=true)
    ep.plotExpStackedBarChart(emissionFile, chartFile, disp=true, hhsCF=hhsFile)
    println("completed")
end

if boxPlotMode
    emissionData = "../India/data/emission/2011_IND_hhs_emission_rng_perCap.csv"
    expendirtureData = "../India/data/emission/2011_IND_hhs_emission_rng_perCap_exp.csv"

    xtitle = "Daily expenditure (USD/day/capita)"
    ytitle = "Annual carbon emission (tCO2/year/capita)"

    println(" Values reading")
    bc.readValues(expendirtureData, emissionData)
    println(" Plot printing")
    bc.plotBoxChart(xtitle, ytitle, disp=true, dash=true, output=Base.source_dir()*"/chart/BoxChart.png")
end

println("[Done]")
