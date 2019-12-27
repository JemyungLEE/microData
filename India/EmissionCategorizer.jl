module EmissionCategorizer

# Developed date: 20. Dec. 2019
# Last modified date: 27. Dec. 2019
# Subject: Categorize India households carbon emissions
# Description: Categorize emissions by districts (province, city, etc) and by expenditure categories
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX

sec = Array{String, 1}()    # India products or services sectors
hhid = Array{String, 1}()   # Household ID

cat = Dict{String, String}()    # category dictionary: {sector, category}
dis = Dict{String, String}()    # hhid's district: {hhid, district}

emissions = Dict{Int16, Array{Float64, 2}}()
emissionsCat = Dict{Int16, Array{Float64, 2}}()     # categozied emission

function readEmission(year, inputFile)

    global sec, hhid

    f = open(inputFile)

    hhid = deleteat!(split(readline(f), '\t'), 1)   # hhid list
    e = zeros(Float64, 0, length(hhid))             # emission matrix: {sector, hhid}
    for l in eachline(f)
        l = split(l, '\t')
        push!(sec, l[1])                            # sec list
        e = vcat(e, map(x->parse(Float64,x),deleteat!(l, 1))')
    end

    global emissions[year] = e
    close(f)

    return e
end

function readSectors(nat, inputFile)

    global cat
    xf = XLSX.readxlsx(inputFile)

    sh = xf[nat]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; cat[r[2]] = r[4] end end
    sh = xf[nat*"_dis"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; dis[r[2]] = r[5] end end

    close(xf)
end

function categorizeEmission(nat, inputFile, outputFile)

    global sec, hhid, cat, dis
    global emissions, emissionsCat

    ns = length(sec)
    nh = length(hhid)
    nc = length(cat)
    nd = length(dis)

    

end

function compareTables(year, inputFiles)

    e = Array{Array{Float64, 2}, 1}()

    for f in inputFiles; push!(e, readEmission(year, f)) end

    println()
    for i = 1:length(inputFiles); print("\t$i") end
    println()
    for i = 1:length(inputFiles)
        print(i)
        for j = 1:i; print("\t") end
        for j = (i+1):length(inputFiles)
            print("\t", isapprox(e[i],e[j]))
        end
        println()
    end
end

end
