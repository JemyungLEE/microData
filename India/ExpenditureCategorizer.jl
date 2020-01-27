module ExpenditureCategorizer

# Developed date: 22. Jan. 2020
# Last modified date: 27. Jan. 2020
# Subject: Categorize India household consumer expenditures
# Description: Categorize expenditures by districts (province, city, etc) and by expenditure categories
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX

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
expcnt = Dict{Int16, Array{Array{Int64, 2}},1}()     # households count by expenditure: {year, {expenditure range, members range}}

hhsize = Dict{Int, Int}()                   # {size, number of households}
hhexp = Dict{Int, Array{Float64, 1}}()      # {size, total expenditure}

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

    return ec, catList, disList, hhexp, hhsize
end

function countHhByExpenditure(year, max=[], min=[])

    global

end

function printCategorizedExpenditureByHHsize(outputFile)

    global hhexp, hhsize, catlist

    f = open(outputFile, "w")

    print(f, "HH_size\tNumber")
    for c in catList; print(f, "\t", c) end
    println(f)
    for s in sort(collect(keys(hhexp)))
        print(f, s, "\t", hhsize[s])
        for i = 1:length(catList)
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
    for c in catList; print(f, "\t", c) end
    println(f)
    for i = 1:length(hhid)
        print(f, hhid[i], "\t", siz[hhid[i]])
        for j = 1:length(catList)
            print(f, "\t", ec[i,j])
        end
        println(f)
    end

    close(f)
end

end
