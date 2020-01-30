module ExpenditureCategorizer

# Developed date: 22. Jan. 2020
# Last modified date: 30. Jan. 2020
# Subject: Categorize India household consumer expenditures
# Description: Categorize expenditures by districts (province, city, etc) and by expenditure categories
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX
using Plots; plotly()
using ORCA

sec = Array{String, 1}()    # India products or services sectors
hhid = Array{String, 1}()   # Household ID

cat = Dict{String, String}()    # category dictionary: {sector code, category}

dis = Dict{String, String}()    # hhid's district: {hhid, district code}
siz = Dict{String, Int}()       # hhid's family size: {hhid, number of members}
rel = Dict{String, Int}()       # hhid's religion: {hhid, religion code}
inc = Dict{String, Float64}()   # hhid's income, monthly per capita expenditure: {hhid, mcpe}

pop = Dict{String, Int}()       # population by district: {district code, population}
hhs = Dict{String, Int}()       # number of households by district: {district code, number of households}

gid = Dict{String, String}()    # districts' gis_codes: {district code, gis id}
catlist = Array{String, 1}()    # category list
dislist = Array{String, 1}()    # district list

exp = Dict{Int16, Array{Float64, 2}}()      # expenditure: {year, {households, India sectors}}
expcat = Dict{Int16, Array{Float64, 2}}()   # categozied expenditure: {year, {households, categories}}
expcnt = Dict{Int16, Array{Array{Int64, 2},1}}()    # households count by expenditure: {year, {category, {expenditure range, hh size}}}

hhsize = Dict{Int16, Dict{Int, Int}}()      # hh frequancy by size: {year, {size, number of households}}
hhexp = Dict{Int16, Dict{Int, Float64}}()   # hh total expenditure by size: {year, {size, total expenditure}}
hhavg = Dict{Int16, Dict{Int, Float64}}()   # hh average expenditure by size+: {year, {category, {size}}}

function getExpenditureData(year, expData)

    global exp[year] = expData[1]
    global hhid = expData[2]
    global sec = expData[3]
end

function getHouseholdData(year, households)

    global dis, siz, rel, inc

    for h in collect(keys(households))
        dis[h] = households[h].district
        siz[h] = households[h].size
        rel[h] = households[h].religion
        inc[h] = households[h].mpceMrp
    end
end

function readCategoryData(nat, inputFile)

    global cat, gid, pop, hhs
    xf = XLSX.readxlsx(inputFile)

    sh = xf[nat*"_sec"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; cat[string(r[1])] = r[4] end end
    sh = xf[nat*"_dist"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; gid[string(r[1])] = r[3] end end
    sh = xf[nat*"_pop"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1 && !ismissing(r[3])
        pop[string(r[3])] = r[9]
        hhs[string(r[3])] = r[8]
    end end

    close(xf)
end

function categorizeExpenditure(year)

    global sec, hhid, cat, dis, siz
    global exp, expcat
    global hhsize, hhexp

    global catlist = sort(unique(values(cat)))      # category list
    push!(catlist, "Total")
    global dislist = sort(unique(values(dis)))      # district list
    nc = length(catlist)
    nd = length(dislist)

    indCat = Dict{String, Int}()     # index dictionary of category
    indDis = Dict{String, Int}()     # index dictionary of district

    # make index dictionaries
    for s in sec; indCat[s] = findfirst(x->x==cat[s], catlist) end
    for h in hhid; indDis[h] = findfirst(x->x==dis[h], dislist) end

    # categorize expenditure data by ten sectors
    e = exp[year]
    ec = zeros(Float64, length(hhid), nc)
    for i=1:length(hhid); for j=1:length(sec); ; ec[i, indCat[sec[j]]] += e[i,j] end end
    # summing
    for i=1:length(hhid); for j=1:nc-1; ec[i, nc] += ec[i,j] end end

    # categorize expenditure data by household size
    for i=1:length(hhid)
        if haskey(hhsize, siz[hhid[i]])
            hhsize[siz[hhid[i]]] += 1
            hhexp[siz[hhid[i]]] += ec[i,:]
        else
            hhsize[siz[hhid[i]]] = 1
            hhexp[siz[hhid[i]]] = ec[i,:]
        end
    end
    for s in collect(keys(hhsize)); hhexp[s] /= hhsize[s] end

    expcat[year] = ec

    return ec, catlist, dislist, hhexp, hhsize
end

function setIntervals(year, nrow = 20, rmax = 1, rmin = 0)  # number of row sections, ratios of over-max and under-min.
    global expcat, catlist
    ec = expcat[year]
    max = zeros(length(catlist))
    min = zeros(length(catlist))



    return max, min
end

function countByExpenditure(year, nrow = 20, maxexp=[], minexp=[], maxsiz = 20)

    global expcnt, expcat, expavg
    global hhid, catlist
    ec = expcat[year]

    # prepare counting
    expdic = Dict{String, Array{Float64, 1}}()
    for i = 1:length(hhid); expdic[hhid[i]] = ec[i,:] end
    maxhhsiz = maximum(collect(keys(hhsize)))
    maxhhexp = collect(maximum(ec[:,i]) for i=1:length(catlist))
    minhhexp = collect(minimum(ec[:,i]) for i=1:length(catlist))

    # counting process
    cntcat = []     # categorized counts
    avgcat = []     # categorized average
    col = [1:maxsiz;]
    rowlist = []
    for i = 1:length(catlist)
        tmpmin = 0.0
        if length(minexp) > 0; tmpmin = minexp[i] end
        if length(maxexp) > 0; row = collect((j*(maxexp[i]-tmpmin)/nrow + tmpmin) for j=0:nrow)
        else row = collect((j*(maxhhexp[i]-tmpmin)/nrow + tmpmin) for j=0:nrow)
        end

        cnt = zeros(Int, length(row), length(col)+1)
        for j=1:nrow
            for k=1:maxsiz; cnt[j,k]=count(x->((col[k]==siz[x])&&(row[j]<=expdic[x][i]<row[j+1])), hhid) end
            cnt[j,end] = count(x->((maxsiz<siz[x])&&(row[j]<=expdic[x][i]<row[j+1])), hhid)
        end
        for k=1:maxsiz; cnt[end,k]=count(x->((col[k]==siz[x])&&(row[end]<=expdic[x][i])), hhid) end
        cnt[end,end] = count(x->((maxsiz<siz[x])&&(row[end]<=expdic[x][i])), hhid)

        avg = zeros(Float64, length(col)+1)
        for j=1:maxsiz

        end

        push!(cntcat, cnt)
        push!(avgcat, avg)
        push!(rowlist, row)
    end

    expcnt[year] = cntcat

    return cntcat, rowlist, col, maxhhsiz, maxhhexp, minhhexp
end

function printCountedResult(year, outputFile, rowlist=[], col=[], maxhhsiz=0, maxhhexp=[], minhhexp=[])

    global expcnt, catlist
    cntcat = expcnt[year]

    f = open(outputFile, "w")

    for i=1:length(catlist)
        println(f, "[",catlist[i],"]")
        cnt = cntcat[i]

        if maxhhsiz>0; println(f,"Max. HH size: ", maxhhsiz) end
        if length(maxhhexp)>0; println(f,"Max. HH expenditure: ", maxhhexp[i]) end
        if length(minhhexp)>0; println(f,"Min. HH expenditure: ", minhhexp[i]) end
        print(f, "Exp./Size")
        if length(col)>0
            for j in col; print(f,"\t",j) end
            println(f,"\t<",col[end])
        end

        for j=1:size(cnt,1)
            if length(rowlist)>0; print(f, "<", trunc(Int, rowlist[i][j])) end
            for k=1:size(cnt,2); print(f,"\t", cnt[j,k]) end
            println(f)
        end
        println(f)
    end

    close(f)
end

function plotHeatmap(year, rowlist=[], col=[], dispmode =false, filename="")
    global expcnt, catlist
    cntcat = expcnt[year]

    plotlist = []
    for i=1:length(catlist)
        if length(rowlist)>0 && length(col)>0
            push!(plotlist, Plots.heatmap(col, rowlist[i], cntcat[i], title = catlist[i]))
        else push!(plotlist, Plots.heatmap(cntcat[i], title = catlist[i]))
        end
    end

    if dispmode; for p in plotlist; display(p) end end
    if length(filename)>0; for i=1:length(catlist); ORCA.savefig(plotlist[i], filename*"_"*catlist[i]*".png") end end

    return plotlist
end

function printCategorizedExpenditureByHHsize(outputFile)

    global hhexp, hhsize, catlist

    f = open(outputFile, "w")

    print(f, "HH_size\tNumber")
    for c in catlist; print(f, "\t", c) end
    println(f)
    for s in sort(collect(keys(hhexp)))
        print(f, s, "\t", hhsize[s])
        for i = 1:length(catlist)
            print(f, "\t", hhexp[s][i])
        end
        println(f)
    end

    close(f)
end

function printCategorizedExpenditureByHH(year, outputFile)

    global expcat, hhid, catlist
    ec = expcat[year]

    f = open(outputFile, "w")

    print(f, "HHID\tSize")
    for c in catlist; print(f, "\t", c) end
    println(f)
    for i = 1:length(hhid)
        print(f, hhid[i], "\t", siz[hhid[i]])
        for j = 1:length(catlist)
            print(f, "\t", ec[i,j])
        end
        println(f)
    end

    close(f)
end

end
