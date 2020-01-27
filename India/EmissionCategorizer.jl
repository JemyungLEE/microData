module EmissionCategorizer

# Developed date: 20. Dec. 2019
# Last modified date: 20. Jan. 2020
# Subject: Categorize India households carbon emissions
# Description: Categorize emissions by districts (province, city, etc) and by expenditure categories
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX

sec = Array{String, 1}()    # India products or services sectors
hhid = Array{String, 1}()   # Household ID

cat = Dict{String, String}()    # category dictionary: {sector code, category}
dis = Dict{String, String}()    # hhid's district: {hhid, district code}
siz = Dict{String, Int}()       # hhid's family size: {hhid, number of members}
pop = Dict{String, Tuple{Int,Int}}() # population by district: {district code, (population, number of households)}
hhs = Dict{String, Int}()       # number of households by district: {district code}

gid = Dict{String, String}()    # districts' gis_codes: {district code, gis id}

emissions = Dict{Int16, Array{Float64, 2}}()        # {year, table}

catList = Array{String, 1}()    # category list
disList = Array{String, 1}()    # district list
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

function readHouseholdData(year, inputFile)

    global siz, dis
    f = open(inputFile)

    readline(f)
    for l in eachline(f)
        l = split(l, '\t')
        siz[l[1]] = parse(Int,l[7])
        dis[l[1]] = l[5]
    end

    close(f)
end

function readCategoryData(nat, inputFile)

    global cat, dis, gid
    xf = XLSX.readxlsx(inputFile)

    sh = xf[nat*"_sec"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; cat[string(r[1])] = r[4] end end
    sh = xf[nat*"_dist"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; gid[string(r[1])] = r[3] end end
    sh = xf[nat*"_pop"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1 && !ismissing(r[3]); pop[string(r[3])] = (r[9], r[8]) end end
    close(xf)

    return cat, dis, gid
end

function categorizeEmission(year, weightMode = 0, squareRoot = false)
    # weightMode: [0]non-weight, [1]population weighted, [2]household weighted, [3]both population and household weighted
    #             ([4],[5]: normalization) [4]per capita, [5]per household
    # squareRoot: [true]apply square root of household size for an equivalance scale

    global sec, hhid, cat, dis, siz
    global emissions, emissionsCat

    global catList = sort(unique(values(cat)))      # category list
    push!(catList, "Total")
    global disList = sort(unique(values(dis)))      # district list
    nc = length(catList)
    nd = length(disList)

    indCat = Dict{String, Int}()     # index dictionary of category
    indDis = Dict{String, Int}()     # index dictionary of district

    # make index dictionaries
    for s in sec; indCat[s] = findfirst(x->x==cat[s], catList) end
    for h in hhid; indDis[h] = findfirst(x->x==dis[h], disList) end

    # sum households and members by districts
    thbd = zeros(Float64, length(disList))   # total households by district
    tpbd = zeros(Float64, length(disList))   # total members of households by district
    for h in hhid
        thbd[indDis[h]] += 1
        if squareRoot; tpbd[indDis[h]] += sqrt(siz[h])
        else tpbd[indDis[h]] += siz[h]
        end
    end

    # categorize emission data
    e = emissions[year]
    ec = zeros(Float64, nc, nd)
    # categorizing
    for i=1:length(sec); for j=1:length(hhid); ec[indCat[sec[i]],indDis[hhid[j]]] += e[i,j] end end
    # summing
    for i = 1:nd; for j = 1:nc-1; ec[nc, i] += ec[j,i] end end

    # weighting
    if weightMode == 1
        for i=1:nc; for j=1:nd
            if haskey(pop, disList[j]); ec[i,j] *= pop[disList[j]][1]/tpbd[j]
            else ec[i,j] = 0
            end
        end end
    elseif weightMode == 2
        for i=1:nc; for j=1:nd
            if haskey(pop, disList[j]); ec[i,j] *= pop[disList[j]][2]/thbd[j]
            else ec[i,j] = 0
            end
        end end
    elseif weightMode == 3
        for i=1:nc; for j=1:nd
            if haskey(pop, disList[j]); ec[i,j] *= 0.5*(pop[disList[j]][1]/tpbd[j]+pop[disList[j]][2]/thbd[j])
            else ec[i,j] = 0
            end
        end end
    # normalizing
    elseif weightMode == 4
        for i=1:nc; for j=1:nd; ec[i,j] /= tpbd[j] end end
    elseif weightMode == 5
        for i=1:nc; for j=1:nd; ec[i,j] /= thbd[j] end end
    # basic information
    elseif weightMode == 6
        for i=1:nd;
            ec[1,i] = tpbd[i]
            ec[2,i] = thbd[i]
        end
    end

    emissionsCat[year] = ec

    return ec, catList, disList
end

function exportEmissionTable(year, tag, outputFile)

    global gid, catList, disList, emissionsCat
    ec = emissionsCat[year]

    # preparing exporting table
    gidList = sort(unique(values(gid)))
    tb = zeros(Float64, length(gidList), length(catList))
    for i=1:length(disList)
        idx = findfirst(x->x==gid[disList[i]],gidList)
        for j=1:length(catList)
            tb[idx, j] += ec[j, i]
        end
    end

    # exporting table
    f = open(outputFile, "w")
    print(f, tag)
    for c in catList; print(f,",",c) end
    println(f)
    for i = 1:size(tb, 1)
        print(f, gidList[i])
        for j = 1:size(tb, 2); print(f, ",", tb[i,j]) end
        println(f)
    end

    close(f)

end

function printCategorizedEmission(year, outputFile)

    global catList, disList, emissionsCat
    ec = emissionsCat[year]

    f = open(outputFile, "w")

    for d in disList; print(f, "\t", d) end
    println(f)
    for i = 1:length(catList)
        print(f, catList[i])
        for j = 1:length(disList); print(f, "\t", ec[i,j]) end
        println(f)
    end

    close(f)
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
