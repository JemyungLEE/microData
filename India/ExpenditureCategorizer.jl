module ExpenditureCategorizer

# Developed date: 22. Jan. 2020
# Last modified date: 4. Feb. 2020
# Subject: Categorize India household consumer expenditures
# Description: Categorize expenditures by districts (province, city, etc) and by expenditure categories
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using XLSX
using Plots
using Statistics
using KernelEstimator

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
expcnt = Dict{Int16, Array{Array{Int64, 2},1}}()    # household frequancy: {year, {category, {expenditure range, hh size}}}
expavg = Dict{Int16, Array{Array{Float64, 1},1}}()  # average expenditure: {year, {category, {hh size}}}
expsrs = Dict{Int16, Array{Array{Float64, 1},1}}()  # hh size square root scale expenditure: {year, {category, {hh size}}}
expavgreg = Dict{Int16, Array{Array{Float64, 1},1}}()   # regression estimation of average expenditure: {year, {category, {hh size}}}
expsrsreg = Dict{Int16, Array{Array{Float64, 1},1}}()   # regression estimation of hh size square root scale expenditure: {year, {category, {hh size}}}

hhsize = Dict{Int16, Dict{Int, Int}}()      # hh frequancy by hh size: {year, {size, number of households}}
hhexp = Dict{Int16, Dict{Int, Array{Float64, 1}}}() # average hh expenditure by hh size: {year, {size, {category}}}

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

    indCat = Dict{String, Int}()    # index dictionary of category
    indDis = Dict{String, Int}()    # index dictionary of district
    hs = Dict{Int, Int}()           # hh frequancy by hh size
    he = Dict{Int, Array{Float64, 1}}() # average hh total expenditure by hh size

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
        if haskey(hs, siz[hhid[i]])
            hs[siz[hhid[i]]] += 1
            he[siz[hhid[i]]] += ec[i,:]
        else
            hs[siz[hhid[i]]] = 1
            he[siz[hhid[i]]] = ec[i,:]
        end
    end
    for s in collect(keys(hs)); he[s] /= hs[s] end

    expcat[year] = ec
    hhexp[year] = he
    hhsize[year] = hs

    return ec, catlist, dislist, he, hs
end

function setIntervals(year, nrow = 20, rmax = 1, rmin = 0, logscale=false)
                            # number of row sections, ratios of over-max and under-min, logarithm scale
    global expcat, catlist
    ec = expcat[year]
    max = zeros(length(catlist))
    min = zeros(length(catlist))



    return max, min
end

function countByExpenditure(year, nrow = 20, maxexp=[], minexp=[], maxsiz = 20)

    global expcnt, expcat, expavg, expsrs
    global hhid, catlist
    ec = expcat[year]

    # Prepare counting
    expdic = Dict{String, Array{Float64, 1}}()      # expenditure: {hhid, {category}}
    for i = 1:length(hhid); expdic[hhid[i]] = ec[i,:] end
    maxhhsiz = maximum(collect(keys(hhsize)))
    maxhhexp = collect(maximum(ec[:,i]) for i=1:length(catlist))
    minhhexp = collect(minimum(ec[:,i]) for i=1:length(catlist))

    # Counting process
    cntcat = []     # categorized counts
    avgcat = []     # categorized average
    srscat = []     # categorized hh size square root scale expenditure average
    col = [1:maxsiz;]
    rowlist = []
    for i = 1:length(catlist)
        tmpmin = 0.0
        if length(minexp) > 0; tmpmin = minexp[i] end
        if length(maxexp) > 0; row = collect((j*(maxexp[i]-tmpmin)/nrow + tmpmin) for j=0:nrow)
        else row = collect((j*(maxhhexp[i]-tmpmin)/nrow + tmpmin) for j=0:nrow)
        end

        # Count frequancy
        cnt = zeros(Int, length(row), maxsiz+1)
        for j=1:nrow
            for k=1:maxsiz; cnt[j,k]=count(x->((col[k]==siz[x])&&(row[j]<=expdic[x][i]<row[j+1])), hhid) end
            cnt[j,end] = count(x->((col[end]<siz[x])&&(row[j]<=expdic[x][i]<row[j+1])), hhid)
        end
        for k=1:maxsiz; cnt[end,k]=count(x->((col[k]==siz[x])&&(row[end]<=expdic[x][i])), hhid) end
        cnt[end,end] = count(x->((col[end]<siz[x])&&(row[end]<=expdic[x][i])), hhid)

        # Calculate average
        avg = zeros(Float64, maxsiz+1)
        for j=1:maxsiz; avg[j] = mean(expdic[k][i] for k in hhid if col[j]==siz[k]) end
        avg[end] = mean(expdic[k][i] for k in hhid if col[end]<siz[k])

        # Calculate square root scale average
        srs = zeros(Float64, maxsiz+1)
        for j=1:maxsiz; srs[j] = avg[j]/sqrt(j) end
        srs[end] = mean(expdic[k][i]/sqrt(siz[k]) for k in hhid if col[end]<siz[k])

        push!(cntcat, cnt)
        push!(avgcat, avg)
        push!(srscat, srs)
        push!(rowlist, row)
    end

    expcnt[year] = cntcat
    expavg[year] = avgcat
    expsrs[year] = srscat

    return cntcat, avgcat, srscat, rowlist, col, maxhhsiz, maxhhexp, minhhexp
end

function nonparreg(year, col, linear = true)  #

    global expavg, expsrs, expavgreg, expsrsreg
    avgcat = expavg[year]
    srscat = expsrs[year]
    avgreg = []
    srsreg = []

    for i = 1:length(catlist)
        avg = avgcat[i]
        srs = srscat[i]

        if linear
            push!(avgreg, npr(col, avg[1:end-1], xeval=col, reg=locallinear))
            push!(srsreg, npr(col, srs[1:end-1], xeval=col, reg=locallinear))
        else
            push!(avgreg, npr(col, avg[1:end-1], xeval=col, reg=localconstant))
            push!(srsreg, npr(col, srs[1:end-1], xeval=col, reg=localconstant))
        end
    end

    expavgreg[year] = avgreg
    expsrsreg[year] = srsreg

    return avgreg, srsreg
end

function printCountedResult(year, outputFile, rowlist=[], col=[], maxhhsiz=0, maxhhexp=[], minhhexp=[])

    global expcnt, expavg, catlist
    cntcat = expcnt[year]
    avgcat = expavg[year]

    f = open(outputFile, "w")

    for i=1:length(catlist)
        println(f, "[",catlist[i],"]")
        cnt = cntcat[i]
        avg = avgcat[i]

        if maxhhsiz>0; println(f,"Max. HH size: ", maxhhsiz) end
        if length(maxhhexp)>0; println(f,"Max. HH expenditure: ", maxhhexp[i]) end
        if length(minhhexp)>0; println(f,"Min. HH expenditure: ", minhhexp[i]) end
        print(f, "Exp./Size")
        if length(col)>0
            for j in col; print(f,"\t",j) end
            println(f,"\t<",col[end])
        else for j = 1:size(cnt,2); print(f,"\t",j) end
        end

        for j=1:size(cnt,1)
            if length(rowlist)>0; print(f, "<", trunc(Int, rowlist[i][j]))
            else print(f, j)
            end
            for k=1:size(cnt,2); print(f,"\t", cnt[j,k]) end
            println(f)
        end

        print(f, "Avg.")
        for j=1:size(cnt,2); print(f,"\t", avg[j]) end
        println(f)

        println(f)
    end

    close(f)
end

function plotHeatmap(year, rowlist=[], col=[], dispmode=false, guimode=false, logmode=false, filename="")

    plotly()
    global expcnt, expavg, expsrs, catlist
    cntcat = expcnt[year]
    avgcat = expavg[year]
    srscat = expsrs[year]
    avgreg = expavgreg[year]
    srsreg = expsrsreg[year]

    if length(col)>0;
        collist = [string(c) for c in col]
        push!(collist, "<"*collist[end])
    end

    plotlist = []
    for i=1:length(catlist)
        if logmode; yx=("Exp.(USD)", :log)
        else yx=("Exp.(USD)")
        end
        if length(rowlist)>0 && length(col)>0
            p = heatmap(collist, rowlist[i], cntcat[i], title=catlist[i], xaxis=("HH size"), yaxis=yx, legend=:outertopright)
            p = plot!(collist, avgcat[i], label = "Avg.", width=3, legend=:inside)
            p = plot!(collist, srscat[i], label = "SqRtSc.", width=3, legend=:inside)
            p = plot!(collist[1:end-1], avgreg[i], label = "Avg.Reg.", width=1, legend=:inside)
            p = plot!(collist[1:end-1], srsreg[i], label = "SqRtSc.Reg.", width=1, legend=:inside)
            push!(plotlist, p)
        else
            p = heatmap(cntcat[i], title = catlist[i], xaxis="HH size", yaxis=yx, legend=:outertopright)
            p = plot!(avgcat[i], label = "Avg.", width=3, legend=:inside)
            p = plot!(srscat[i], label = "SqRtSc.", width=3, legend=:inside)
            p = plot!(avgreg[i], label = "Avg.Reg.", width=1, legend=:inside)
            p = plot!(srsreg[i], label = "SqRtSc.Reg.", width=1, legend=:inside)
            push!(plotlist, p)
        end

        if dispmode; display(p) end
        if guimode; gui() end
        #if length(filename)>0; png(p, filename*"_"*catlist[i]) end
    end

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
