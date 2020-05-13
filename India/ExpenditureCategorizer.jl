module ExpenditureCategorizer

# Developed date: 22. Jan. 2020
# Last modified date: 13. May. 2020
# Subject: Categorize India household consumer expenditures
# Description: Categorize expenditures by districts (province, city, etc) and by expenditure categories
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using XLSX
using Plots
using Statistics
using KernelEstimator

sec = Array{String, 1}()    # India products or services sectors
secName = Dict{String, String}()    # sector name dictionary: {sector code, name}
hhid = Array{String, 1}()   # Household ID

cat = Dict{String, String}()    # category dictionary: {sector code, category}

sta = Dict{String, String}()    # hhid's state: {hhid, state code}
dis = Dict{String, String}()    # hhid's district: {hhid, district code}
siz = Dict{String, Int}()       # hhid's family size: {hhid, number of members}
typ = Dict{String, String}()    # hhid's sector type, urban or rural: {hhid, "urban" or "rural"}
rel = Dict{String, Int}()       # hhid's religion: {hhid, religion code}
inc = Dict{String, Float64}()   # hhid's income, monthly per capita expenditure: {hhid, mcpe}

# pop = Dict{String, Int}()       # population by district: {district code, population}
# hhs = Dict{String, Int}()       # number of households by district: {district code, number of households}
pop = Dict{String, Tuple{Int,Int,Float64,Float64}}()    # population by district: {district code, (population,households,area(km^2),density(persons/km^2))}

nam = Dict{String, String}()    # districts' name: {district code, district name}
gid = Dict{String, String}()    # districts' gis_codes: {district code, gis id (GIS_2)}
gidData = Dict{String, Tuple{String, String, String, String}}() # GID code data: {GID_2, {district code, district name, state code, state name}}
merDist = Dict{String, String}()    # list of merged district: {merged district's code, remained district's code}
misDist = Array{String, 1}()    # list of missing district: {GID_2}

catList = Array{String, 1}()    # category list
staList = Array{String, 1}()    # state list
disList = Array{String, 1}()    # district list
incList = Array{Float64, 1}()   # income sector list

wghSta = Dict{String, Float64}()    # hhid's state-population weight: {hhid, weight}
wghDis = Dict{String, Float64}()    # hhid's district-population weight: {hhid, weight}

exp = Dict{Int16, Array{Float64, 2}}()      # expenditure: {year, {households, India sectors}}
expcat = Dict{Int16, Array{Float64, 2}}()   # categozied expenditure: {year, {households, categories}}
expinc = Dict{Int16, Array{Float64, 2}}()   # categozied expenditure: {year, {expenditure-level, categories}}

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
end

function getHouseholdData(year, households, merging=false; period="monthly")

    global sta, dis, siz, rel, inc, staList, disList

    for h in collect(keys(households))
        if merging==true&&haskey(merDist, households[h].district); dis[h] = merDist[households[h].district]
        else dis[h] = households[h].district
        end
        sta[h] = households[h].state
        siz[h] = households[h].size
        typ[h] = households[h].sector
        rel[h] = households[h].religion
        inc[h] = households[h].mpceMrp
    end

    # convert MPCE's period
    if period=="daily"; for h in hhid; inc[h] = inc[h]/30 end
    elseif period=="annual"; mmtoyy = 365/30; for h in hhid; inc[h] = inc[h]*mmtoyy end
    end

    staList = sort(unique(values(sta)))      # state list
    disList = sort(unique(values(dis)))      # district list
end

function readCategoryData(nat, inputFile; subCategory="", except=[])

    global sec, secName, cat, gid, nam, pop, hhs, catList, gidData, merDist, misDist
    xf = XLSX.readxlsx(inputFile)

    sh = xf[nat*"_sec"]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r)>1
            secCode = string(r[1])  # sector code
            push!(sec, secCode)
            if length(subCategory)==0 && !ismissing(r[4]) && !(string(r[4]) in except); cat[secCode] = string(r[4])
            elseif subCategory=="Food" && !ismissing(r[5]); cat[secCode] = string(r[5])
            end
            secName[secCode]=string(r[2])
        end
    end
    sh = xf[nat*"_dist"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1; gid[string(r[1])] = string(r[3]); nam[string(r[1])] = string(r[2]) end end
    sh = xf[nat*"_pop"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1 && !ismissing(r[3]); pop[string(r[3])] = (r[9], r[8], r[12], r[9]/r[12]) end end
    sh = xf["2011GID"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1; gidData[string(r[3])]=(string(r[7]),string(r[4]),string(r[6]),string(r[5])) end end
    #=
    sh = xf[nat*"_gid"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1
        codes = split(r[3],r"[._]")
        gidData[string(r[3])]=(codes[3],r[4],codes[2],r[2])
    end end
    =#
    # Read merging districts
    sh = xf[nat*"_mer"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1; merDist[string(r[3])] = string(r[1]) end end
    close(xf)

    # Search data-missing district(s)
    gidList = collect(values(gid))
    for gid in collect(keys(gidData)) if !(gid in gidList); push!(misDist, gid) end end

    catList = sort(unique(values(cat)))
    if length(subCategory)==0; push!(catList, "Total") # category list
    else push!(catList, subCategory)
    end

end

function calculateStatePopulationWeight(populationFile)

    global staList, wghSta

    stapop = Dict{String, Tuple{Int, Int, Int}}()   # State population, {State code, population{total, rural, urban}}
    stasmp = Dict{String, Array{Int, 1}}()          # State sample size, {State code, sample number{total, rural, urban}}
    stawgh = Dict{String, Array{Float64, 1}}()      # State population weight, {State code, population{total, rural, urban}}

    # read population data
    f = open(populationFile)
    readline(f)
    for l in eachline(f)
        s = split(l, ",")
        push!(staList,string(s[1]))
        stapop[string(s[1])] = (parse(Int, s[4]), parse(Int, s[6]), parse(Int, s[8]))
        stasmp[string(s[1])] = zeros(Int, 3)
    end
    close(f)

    # count sample number
    for h in hhid
        if typ[h] == "rural"; stidx=1; elseif typ[h] == "urban"; stidx=2
        else println("HH sector error: not \"urban\" nor \"rural\"")
        end
        stasmp[sta[h]][1] += siz[h]
        stasmp[sta[h]][stidx+1] += siz[h]
    end

    # calculate weights
    for st in staList
        stawgh[st] = zeros(Float64, 3)
        for i=1:3; stawgh[st][i] = stapop[st][i]/stasmp[st][i] end
    end

    for h in hhid
        if typ[h] == "rural"; wghSta[h] = stawgh[sta[h]][2]
        elseif typ[h] == "urban"; wghSta[h] = stawgh[sta[h]][3]
        else println("Household ",h," sector is wrong")
        end
    end
end

function calculateDistrictPopulationWeight(populationFile, concordanceFile)

    global disList, wghDis

    dispop = Dict{String, Array{Int, 1}}()      # District population, {Survey district code, population{total, rural, urban}}
    dissmp = Dict{String, Array{Int, 1}}()      # District sample size, {Survey district code, sample number{total, rural, urban}}
    diswgh = Dict{String, Array{Float64, 1}}()  # District population weight, {Survey district code, population{total, rural, urban}}

    disConc = Dict{String, String}()     # Population-Survey concordance, {Statistics district code, Survey district code}

    # read concordance data
    xf = XLSX.readxlsx(concordanceFile)
    sh = xf["IND_pop"]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r)>1 && !ismissing(r[3]); disConc[string(r[2])]=string(r[3]) end
    end
    close(xf)

    # read population data
    f = open(populationFile)
    readline(f)
    for l in eachline(f)
        s = string.(split(l, ","))
        discode = s[2]
        if haskey(disConc, discode)
            if s[4]=="Total"; dispop[disConc[discode]] = [parse(Int, s[6]), 0, 0]
            elseif s[4]=="Rural"; dispop[disConc[discode]][2] = parse(Int, s[6])
            elseif s[4]=="Urban"; dispop[disConc[discode]][3] = parse(Int, s[6])
            else println(discode, " does not have \"Total\", \"Urban\" nor \"Rural\" data")
            end
        end
    end
    close(f)
    for dc in disList; dissmp[dc] = zeros(Int, 3) end

    # count sample number
    for h in hhid
        stidx = parse(Int,typ[h]); if !(stidx==1||stidx==2) println("sector index error") end
        dissmp[dis[h]][1] += siz[h]
        dissmp[dis[h]][stidx+1] += siz[h]
    end

    # calculate weights
    for dc in disList
        diswgh[dc] = zeros(Float64, 3)
        diswgh[dc][1] = dispop[dc][1]/dissmp[dc][1]
        if dissmp[dc][2]==0 && dissmp[dc][3]==0; println(dc," does not have samples in both rural and urban areas.")
        elseif dissmp[dc][2]==0; diswgh[dc][2]=0; diswgh[dc][3]=dispop[dc][1]/dissmp[dc][3]
        elseif dissmp[dc][3]==0; diswgh[dc][3]=0; diswgh[dc][2]=dispop[dc][1]/dissmp[dc][2]
        else for i in [2,3]; diswgh[dc][i]=dispop[dc][i]/dissmp[dc][i] end
        end
    end

    for h in hhid
        stidx = parse(Int,typ[h])
        if stidx==1||stidx==2; wghDis[h] = diswgh[dis[h]][stidx+1] else println("sector index error") end

        # if typ[h] == "rural"; wghDis[h] = diswgh[dis[h]][2]
        # elseif typ[h] == "urban"; wghDis[h] = diswgh[dis[h]][3]
        # else println("Household ",h," sector is wrong")
        # end
    end
end

function categorizeExpenditure(year)

    global sec, hhid, cat, dis, siz, catList, disList
    global exp, expcat
    global hhsize, hhexp

    nc = length(catList)
    nd = length(disList)

    indCat = Dict{String, Int}()    # index dictionary of category
    indDis = Dict{String, Int}()    # index dictionary of district
    hs = Dict{Int, Int}()           # hh frequancy by hh size
    he = Dict{Int, Array{Float64, 1}}() # average hh total expenditure by hh size

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

    return ec, catList, disList, he, hs
end

function categorizeExpenditureByIncome(year,intv=[],normMode=0,topExpFile=""; perCap=false,desOrd=false,popWgh=false,wghmode="district")
                                            # intv: proportions between invervals of highest to lowest
                                            # normmode: [1]per capita, [2]per houehold
                                            # perCap: [true] per capita mode, [false] per household mode
                                            # desOrd: [true] descening order of 'intv[]', [false] ascending order
    global wghSta, wghDis; if wghmode=="state"; wgh=wghSta elseif wghmode=="district"; wgh=wghDis end
    global sec, hhid, cat, dis, siz, inc, catList, incList
    global exp, expinc

    if length(intv) == 0; intv = [0.25,0.5,0.75,1.00]
    elseif sort!(intv)[end] != 1; intv /= intv[end]
    end

    nh = length(hhid)
    ns = length(sec)
    nc = length(catList)
    ni = length(intv)

    incArray = [inc[h] for h in hhid]
    incOrder = sortperm(incArray, rev=desOrd)   # desOrd: [true] descening order, [false] ascending order

    # sort India sectors by caterory
    secList = Dict{String, Array{String, 1}}()
    for c in catList; secList[c] = filter(x->cat[x]==c, collect(keys(cat))) end
    secListIdx = Dict{String, Array{Int, 1}}()
    for c in catList; secListIdx[c] = [findfirst(x->x==secList[c][i],sec) for i=1:length(secList[c])] end

    # make index dictionaries
    indCat = Dict{String, Int}()    # index dictionary of category
    for s in sec; indCat[s] = findfirst(x->x==cat[s], catList) end

    # determine sections' starting index and values
    pcidx = []  # current sector's starting index for 'per capita' emissions
    indInc = zeros(Int, nh)
    if perCap   # determine sections if interval ratios are for 'per capita' emissions
        accpop = 0
        idx = 1
        push!(pcidx, idx)
        push!(incList, incArray[incOrder[idx]])
        if popWgh; totpop = sum([siz[h]*wgh[h] for h in hhid])
        else totpop = sum(collect(values(siz)))
        end
        for i=1:nh
            if popWgh; accpop += siz[hhid[incOrder[i]]] * wgh[hhid[incOrder[i]]]
            else accpop += siz[hhid[incOrder[i]]]
            end
            if accpop/totpop > intv[idx] && i < nh
                push!(pcidx, i)
                push!(incList, incArray[incOrder[i]])
                idx += 1
            end
            indInc[incOrder[i]] = idx
        end
    else
        push!(incList, incArray[incOrder[1]])
        for i=1:length(intv); push!(incList, incArray[incOrder[trunc(Int, intv[i]*nh)]]) end
    end

    # set sector index
    if !perCap
        i = 1
        for s = 1:ni
            while i <= trunc(Int, nh*intv[s])
                indInc[incOrder[i]] = s
                i += 1
            end
        end
        indInc[incOrder[nh]] = ni
    end

    # sum households and members by districts
    thbi = zeros(Float64, ni)   # total households by income level
    tpbi = zeros(Float64, ni)   # total members of households by income level
    twpbi = zeros(Float64, ni)  # total state/district-population weighted members of households by income level
    for i= 1:nh; thbi[indInc[i]] += 1 end
    for i= 1:nh; tpbi[indInc[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; twpbi[indInc[i]] += siz[hhid[i]] * wgh[hhid[i]] end end

    # categorize expenditure data
    eh = exp[year]
    ei = zeros(Float64, ni, nc)
    ebc = Dict{String, Array{Float64,2}}()    # overall expenditures by expenditure-level, category, and India sector
    for c in catList; ebc[c]= zeros(Float64, ni, length(secList[c])) end
    if popWgh
        for i=1:nh; for j=1:ns; ; ei[indInc[i], indCat[sec[j]]] += eh[i,j] * wgh[hhid[i]] end end
        for i=1:ni; for j=1:nc-1; ei[i, nc] += eh[i,j] end end
        for c in catList; for i=1:nh; for j=1:length(secList[c]); ebc[c][indInc[i],j] += eh[i,secListIdx[c][j]] * wgh[hhid[i]] end end end
    else for i=1:nh; ei[indInc[i],:] += eh[i,:] end
    end

    # normalizing
    if normMode == 1
        if popWgh
            for i=1:nc; ei[:,i] ./= twpbi end
            for c in catList; for i=1:length(secList[c]); ebc[c][:,i] ./= twpbi end end
        else for i=1:nc; ei[:,i] ./= tpbi end
        end
    elseif normMode == 2; for i=1:nc; ei[:,i] ./= thbi end
    end

    # sort top-expending Inida sectors
    ntc = 10    # number of top-expenditure items per category
    tecn = Array{String,3}(undef, ni, nc, ntc)  # Names of top-expenditure India expending items by category
    tecv = zeros(Float64, ni, nc, ntc)          # Values of top-expenditure India expending items by category
    tecp = zeros(Float64, ni, nc, ntc)          # Values' proportions of top-expenditure India expending items by category
    for i=1:ni
        for j=1:nc-1
            ord = sortperm(ebc[catList[j]][i,:],rev=true)
            if ntc<length(ord); lng=ntc else lng = length(ord) end
            totval = sum(ebc[catList[j]][i,:])
            for k=1:lng
                tecn[i,j,k] = secName[secList[catList[j]][ord[k]]]
                tecv[i,j,k] = ebc[catList[j]][i,ord[k]]
                tecp[i,j,k] = tecv[i,j,k]/totval
            end
        end
    end

    if length(topExpFile)>0
        f = open(topExpFile, "w")
        for i=1:nc-1
            println(f,catList[i])
            for j=1:ni
                print(f,intv[j])
                if ntc<length(secList[catList[i]]); lng=ntc else lng = length(secList[catList[i]]) end
                for k=1:lng; print(f,",\"",tecn[j,i,k],"\",",round(tecp[j,i,k]*100,digits=1),"%") end
                println(f)
            end
            println(f)
        end
        close(f)
    end

    expinc[year] = ei

    return ei, catList, incList, tpbi, thbi, twpbi, indInc
end

function analyzeCategoryComposition(year, output="")
    global sec, secNam, hhid, cat, catList
    global exp, expcat

    nhc = 5 # number of high composition sectors

    nh = length(hhid)
    ns = length(sec)
    nc = length(catList)

    e = exp[year]   # {households, India sectors}}
    ec = expcat[year]   # {households, categories}

    te = [sum(e[:,i]) for i=1:ns]
    tec = [sum(ec[:,i]) for i=1:nc]
    # make index dictionaries
    indCat = [findfirst(x->x==cat[s], catList) for s in sec]

    # analyze composition
    orderSec = Array{Array{String, 1}, 1}()  # high composition sectors' id: {category, {high composition sectors}}
    propSec = Array{Array{Float64, 1}, 1}()  # high composition sectors' proportion: {category, {high composition sectors}}
    for i=1:nc

        catidx = findall(x->x==i, indCat)
        teorder = sortperm([te[idx] for idx in catidx], rev=true)

        nts = length(catidx)
        if nts>nhc; nts = nhc end

        push!(orderSec, [sec[catidx[teorder[j]]] for j=1:nts])
        push!(propSec, [te[catidx[teorder[j]]]/tec[i] for j=1:nts])
    end

    if length(output)>0
        f = open(output, "w")
        print(f, "Category"); for i=1:nhc; print(f, ",Sector_no.",i) end; println(f)
        for i=1:nc
            print(f, catList[i])
            for j=1:length(orderSec[i]); print(f, ",\"",secNam[orderSec[i][j]]," (",round(propSec[i][j],digits=3),")\"") end
            println(f)
        end
        close(f)
    end
end

function setIntervals(year, nrow = 20, rmax = 1, rmin = 0, logscale=false)
                            # number of row sections, ratios of over-max and under-min, logarithm scale
    global expcat, catList
    ec = expcat[year]
    max = zeros(length(catList))
    min = zeros(length(catList))



    return max, min
end

function countByExpenditure(year, nrow = 20, maxexp=[], minexp=[], maxsiz = 20)

    global expcnt, expcat, expavg, expsrs
    global hhid, catList
    ec = expcat[year]

    # Prepare counting
    expdic = Dict{String, Array{Float64, 1}}()      # expenditure: {hhid, {category}}
    for i = 1:length(hhid); expdic[hhid[i]] = ec[i,:] end
    maxhhsiz = maximum(collect(keys(hhsize)))
    maxhhexp = collect(maximum(ec[:,i]) for i=1:length(catList))
    minhhexp = collect(minimum(ec[:,i]) for i=1:length(catList))

    # Counting process
    cntcat = []     # categorized counts
    avgcat = []     # categorized average
    srscat = []     # categorized hh size square root scale expenditure average
    col = [1:maxsiz;]
    rowlist = []
    for i = 1:length(catList)
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

    for i = 1:length(catList)
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

    global expcnt, expavg, catList
    cntcat = expcnt[year]
    avgcat = expavg[year]

    f = open(outputFile, "w")

    for i=1:length(catList)
        println(f, "[",catList[i],"]")
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
    global expcnt, expavg, expsrs, catList
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
    for i=1:length(catList)
        if logmode; yx=("Exp.(USD)", :log)
        else yx=("Exp.(USD)")
        end
        if length(rowlist)>0 && length(col)>0
            p = heatmap(collist, rowlist[i], cntcat[i], title=catList[i], xaxis=("HH size"), yaxis=yx, legend=:outertopright)
            p = plot!(collist, avgcat[i], label = "Avg.", width=3, legend=:inside)
            p = plot!(collist, srscat[i], label = "SqRtSc.", width=3, legend=:inside)
            p = plot!(collist[1:end-1], avgreg[i], label = "Avg.Reg.", width=1, legend=:inside)
            p = plot!(collist[1:end-1], srsreg[i], label = "SqRtSc.Reg.", width=1, legend=:inside)
            push!(plotlist, p)
        else
            p = heatmap(cntcat[i], title = catList[i], xaxis="HH size", yaxis=yx, legend=:outertopright)
            p = plot!(avgcat[i], label = "Avg.", width=3, legend=:inside)
            p = plot!(srscat[i], label = "SqRtSc.", width=3, legend=:inside)
            p = plot!(avgreg[i], label = "Avg.Reg.", width=1, legend=:inside)
            p = plot!(srsreg[i], label = "SqRtSc.Reg.", width=1, legend=:inside)
            push!(plotlist, p)
        end

        if dispmode; display(p) end
        if guimode; gui() end
        #if length(filename)>0; png(p, filename*"_"*catList[i]) end
    end

    return plotlist
end

function printCategorizedExpenditureByHHsize(outputFile)

    global hhexp, hhsize, catList

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

    global expcat, hhid, catList
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
