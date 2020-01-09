module EmissionCategorizer

# Developed date: 20. Dec. 2019
# Last modified date: 9. Jan. 2020
# Subject: Categorize India households carbon emissions
# Description: Categorize emissions by districts (province, city, etc) and by expenditure categories
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX

sec = Array{String, 1}()    # India products or services sectors
hhid = Array{String, 1}()   # Household ID

cat = Dict{String, String}()    # category dictionary: {sector, category}
dis = Dict{String, String}()    # hhid's district: {hhid, district}
siz = Dict{String, Int}()       # hhid's family size: {hhid, number of members}

emissions = Dict{Int16, Array{Float64, 2}}()        # {year, table}
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

function readHousehold(year, inputFile)

    global siz
    f = open(inputFile)

    readline(f)
    for l in eachline(f)
        l = split(l, '\t')
        siz[l[1]] = l[7]
    end

    close(f)
end

function readSectors(nat, inputFile)

    global cat, dis
    xf = XLSX.readxlsx(inputFile)

    sh = xf[nat]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; cat[r[2]] = r[4] end end
    sh = xf[nat*"_dis"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; dis[r[2]] = r[5] end end

    close(xf)

    return cat, dis
end

function categorizeEmission(nat, inputFile, outputFile)

    global sec, hhid, cat, dis, siz
    global emissions, emissionsCat

    ns = length(sec)
    nh = length(hhid)

    cl = sort(unique(values(cat)))      # category list
    dl = sort(unique(values(dis)))      # district list

    indCat = Dict{String, String}()     # index dictionary of category
    indDis = Dict{String, String}()     # index dictionary of district

    # make index dictionaries
    for k in collect(keys(cat)); indCat[k] = findfirst(x->x==cat[k], cl) end
    for k in collect(keys(dis)); indDis[k] = findfirst(x->x==dis[k], dl) end

    # summerize households and members by districts
    thbd = zeros(Int, length(dl))   # total households by district
    tpbd = zeros(Int, length(dl))   # total members of households by district
    for i = 1:nh
        thbd[indDis[hhid[j]]] += 1
        tpbd[indDis[hhid[j]]] += siz[hhid[j]]
    end

    # categorize emission data
    for y in collect(keys(emissions))
        e = emissions[y]
        ec = zeros(Float64, length(cl), length(dl))

        # categorizing
        for i=1:ns; for j=1:nh; ec[indCat[sec[i]],indDis[hhid[j]]] += e[i,j] end end
        # weighting
        for i=1:length(cl); for j=1:length(dl); ec[i,j] *= POP / tpbd[j] end end
        for i=1:length(cl); for j=1:length(dl); ec[i,j] *= HOU / thbd[j] end end
        for i=1:length(cl); for j=1:length(dl); ec[i,j] *= 0.5 * (POP / tpbd[j] + HOU / thbd[j]) end end

        emissionsCat[y] = ec
    end

    return emissionsCat
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
