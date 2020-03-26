module EmissionCategorizer

# Developed date: 20. Dec. 2019
# Last modified date: 26. Mar. 2020
# Subject: Categorize India households carbon emissions
# Description: Categorize emissions by districts (province, city, etc) and by expenditure categories
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using Plots
using XLSX
using Statistics

hhid = Array{String, 1}()   # Household ID
sec = Array{String, 1}()    # India products or services sectors
secName = Dict{String, String}()    # sector name dictionary: {sector code, name}
relName = ["Hinduism","Islam","Christianity","Sikhism","Jainism","Buddhism","Zoroastrianism", "Others", "None"]

cat = Dict{String, String}()    # category dictionary: {sector code, category}
dis = Dict{String, String}()    # hhid's district: {hhid, district code}
sta = Dict{String, String}()    # hhid's state: {hhid, state code}
typ = Dict{String, String}()    # hhid's sector type, urban or rural: {hhid, "urban" or "rural"}
siz = Dict{String, Int}()       # hhid's family size: {hhid, number of members}
inc = Dict{String, Float64}()   # hhid's income: {hhid, monthly per capita expenditure (mixed reference period)}
rel = Dict{String, Int}()       # hhid's religion: {hhid, religion code} [1]Hinduism,[2]Islam,[3]Christianity,[4]Sikhism,[5]Jainism,[6]Buddhism,[7]Zoroastrianism,[9]Others,[0]None
wgh = Dict{String, Float64}()   # hhid's population weight: {hhid, weight}

sam = Dict{String, Tuple{Int,Int}}()    # sample population and households by districct: {district code, (population, number of households)}
pop = Dict{String, Tuple{Int,Int}}()    # population by district: {district code, (population, number of households)}
ave = Dict{String, Float64}()   # average annual expenditure per capita, USD/yr: {district code, mean Avg.Exp./cap/yr}
nam = Dict{String, String}()    # districts' name: {district code, district name}
gid = Dict{String, String}()    # districts' gis_codes: {district code, gis id (GIS_2)}
gidData = Dict{String, Tuple{String, String, String, String}}() # GID code data: {GID_2, {district code, district name, state code, state name}}
merDist = Dict{String, String}()    # list of merged district: {merged district's code, remained district's code}
misDist = Array{String, 1}()    # list of missing district: {GID_2}

emissions = Dict{Int16, Array{Float64, 2}}()        # {year, table}

catList = Array{String, 1}()    # category list
disList = Array{String, 1}()    # district list
relList = Array{String, 1}()    # religion list
incList = Array{Float64, 1}()   # income sector list
levList = Array{Float64, 1}()   # carbon emission level sector list

emissionsHHs = Dict{Int16, Array{Float64, 2}}()     # categozied emission by household: {year, {hhid, category}}
emissionsDis = Dict{Int16, Array{Float64, 2}}()     # categozied emission by district: {year, {district, category}}
emissionsRel = Dict{Int16, Array{Float64, 2}}()     # categozied emission by religion: {year, {religion, category}}
emissionsInc = Dict{Int16, Array{Float64, 2}}()     # categozied emission by incomes: {year, {income level, category}}
emissionsIncRel = Dict{Int16, Array{Float64, 3}}()  # categozied emission by incomes: {year, {religion, income level, category}}
emissionsDisLev = Dict{Int16, Array{Float64, 2}}()  # categozied emission by district emission level: {year, {emission level, category}}
emissionsRng = Dict{Int16, Array{Float64, 2}}()     # categozied emission by expenditure range: {year, {range, category}}
emissionsLev = Dict{Int16, Array{Float64, 2}}()     # categozied emission by emission level: {year, {emission level, category}}


emissionsDisDiff = Dict{Int16, Array{Float64, 2}}() # differences of categozied emission by district: (emission-mean)/mean, {year, {district, category}}

gisDistrictEmission = Dict{Int16, Array{Float64, 2}}()    # GIS version, categozied emission by district: {year, {category, district(GID)}}
gisDistrictEmissionDiff = Dict{Int16, Array{Float64, 2}}() # GIS version, differences of categozied emission by district: (emission-mean)/mean, {year, {category, district(GID)}}
gisDistrictEmissionDiffRank = Dict{Int16, Array{Int, 2}}() # GIS version, difference ranks of categozied emission by district: (emission-mean)/mean, {year, {category, district(GID)}}

function readCategoryData(nat, inputFile; subCategory="")

    global sec, cat, gid, nam, pop, gidData, merDist, misDist
    xf = XLSX.readxlsx(inputFile)

    sh = xf[nat*"_sec"]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r)>1
            secCode = string(r[1])  # sector code
            push!(sec, secCode)
            if length(subCategory)==0; cat[secCode] = string(r[4])
            elseif subCategory=="Food" && !ismissing(r[5]); cat[secCode] = string(r[5])
            end
            secName[secCode]=string(r[2])
        end
    end
    sh = xf[nat*"_dist"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1; gid[string(r[1])] = r[3]; nam[string(r[1])] = r[2] end end
    sh = xf[nat*"_pop"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1 && !ismissing(r[3]); pop[string(r[3])] = (r[9], r[8]) end end
    sh = xf[nat*"_gid"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1
        codes = split(r[3],r"[._]")
        gidData[string(r[3])]=(codes[3],r[4],codes[2],r[2])
    end end
    # Read merging districts
    sh = xf[nat*"_mer"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1; merDist[string(r[3])] = string(r[1]) end end
    close(xf)
    # Search data-missing district(s)
    gidList = collect(values(gid))
    for gid in collect(keys(gidData)) if !(gid in gidList); push!(misDist, gid) end end

    global catList = sort(unique(values(cat)))
    if length(subCategory)==0; push!(catList, "Total") # category list
    else push!(catList, subCategory)
    end
end

function readHouseholdData(year, inputFile, merging=false; period="monthly")
                                                            # period: "monthly"(default), "daily", or "annual"

    global hhid, siz, dis, typ, inc, rel, disList, relList
    f = open(inputFile)

    readline(f)
    for l in eachline(f)
        l = string.(split(l, '\t'))
        push!(hhid, l[1])
        siz[l[1]] = parse(Int,l[7])
        sta[l[1]] = l[4]
        if merging==true&&haskey(merDist, l[5]); dis[l[1]] = merDist[l[5]]
        else dis[l[1]] = l[5]
        end
        typ[l[1]] = l[6]
        inc[l[1]] = parse(Float64,l[8])
        rel[l[1]] = parse(Int,l[9])
        wgh[l[1]] = parse(Float64,l[12])
    end
    # convert MPCE's period
    if period=="daily"; for h in hhid; inc[h] = inc[h]/30 end
    elseif period=="annual"; mmtoyy = 365/30; for h in hhid; inc[h] = inc[h]*mmtoyy end
    end

    close(f)

    global disList = sort(unique(values(dis)))      # district list
    global relList = sort(unique(values(rel)))      # religion list
    if 0 in relList; deleteat!(relList, findall(x->x==0, relList)); push!(relList, 0) end
end

function readEmission(year, inputFile)

    global hhid

    f = open(inputFile)
    readline(f)
    e = zeros(Float64, length(sec), length(hhid))             # emission matrix: {sector, hhid}
    i = 1
    for l in eachline(f)
        l = split(l, '\t')[2:end]
        e[i,:] = map(x->parse(Float64,x),l)
        i += 1
    end
    global emissions[year] = e
    close(f)
end

function readExpenditure(year, inputFile)

    global sec, hhid
    nh = length(hhid)
    f = open(inputFile)
    readline(f)
    e = zeros(Float64, length(sec), nh)      # expenditure matrix: {sector, hhid}
    i = 1
    for l in eachline(f)
        l = split(l, '\t')[2:end-1]
        if i<=nh; e[:,i] = map(x->parse(Float64,x),l) end
        i += 1
    end

    global emissions[year] = e
    close(f)

    return e
end

function categorizeHouseholdEmission(year; output="", hhsinfo=false)
    global sec, hhid, cat, catList
    global emissions, emissionsHHs

    nc = length(catList)
    nh = length(hhid)
    ns = length(sec)

    # make index dictionaries
    indCat = [if haskey(cat, s); findfirst(x->x==cat[s], catList) end for s in sec]

    # categorize emission data
    e = emissions[year]
    ec = zeros(Float64, nh, nc)
    # categorizing
    for i=1:nh; for j=1:ns; if indCat[j]!=nothing; ec[i,indCat[j]] += e[j,i] end end end
    # summing
    for i=1:nc-1; ec[:, nc] += ec[:,i] end

    # save the results
    emissionsHHs[year] = ec

    # print the results
    if length(output)>0
        f = open(output, "w")
        print(f,"HHID"); for c in catList; print(f, ",", c) end
        if hhsinfo; print(f, ",HH_size,MPCE") end; println(f)
        for i = 1:length(hhid)
            print(f, hhid[i]); for j = 1:length(catList); print(f, ",", ec[i,j]) end
            if hhsinfo; print(f, ",",siz[hhid[i]],",",inc[hhid[i]]); println(f) end
        end
        close(f)
    end
end

function analyzeCategoryComposition(year, output="")
    global sec, hhid, cat, catlist
    global emissions, emissionsHHs

    nhc = 5 # number of high composition sectors

    nh = length(hhid)
    ns = length(sec)
    nc = length(catlist)

    e = emissions[year]         # {India sectors, hhid}}
    ec = emissionsHHs[year]     # {hhid, category}

    te = [sum(e[i,:]) for i=1:ns]
    tec = [sum(ec[:,i]) for i=1:nc]
    # make index dictionaries
    indCat = [findfirst(x->x==cat[s], catlist) for s in sec]

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
        print(f, "Category"); for i=1:nts; print(f, "\tSector_no.",i) end; println(f)
        for i=1:nc
            print(f, catlist[i])
            for j=1:length(orderSec[i]); print(f, "\t",secName[orderSec[i][j]]," (",round(propSec[i][j],digits=3),")") end
            println(f)
        end
        close(f)
    end
end

function categorizeDistrictEmission(year, weightMode=0; sqrRoot=false, period="monthly", religion=false)
    # weightMode: [0]non-weight, [1]population weighted, [2]household weighted, [3]both population and household weighted
    #             ([4],[5]: normalization) [4]per capita, [5]per household
    # sqrRoot: [true]apply square root of household size for an equivalance scale
    # period: "monthly", "daily", or "annual"
    # religion: [true] categorize districts' features by religions

    global hhid, cat, dis, siz, inc, sam, ave, rel
    global emissionsHHs, catList, disList, relList
    global emissionsDis, emissionsDisDiff

    nh = length(hhid)
    nc = length(catList)
    nd = length(disList)
    nr = length(relList)

    # make index arrays
    indDis = [findfirst(x->x==dis[hhid[i]], disList) for i=1:nh]
    if religion; indRel = [findfirst(x->x==rel[hhid[i]], relList) for i=1:nh] end

    # sum sample households and members by districts
    thbd = zeros(Float64, nd)   # total households by district
    tpbd = zeros(Float64, nd)   # total members of households by district
    for i=1:nh; thbd[indDis[i]] += 1 end
    if sqrRoot; for i=1:nh; tpbd[indDis[i]] += sqrt(siz[hhid[i]]) end
    else for i=1:nh; tpbd[indDis[i]] += siz[hhid[i]] end
    end
    for i=1:nd; sam[disList[i]] = (tpbd[i], thbd[i]) end

    # sum sample households and members by districts and by religions
    if religion
        thbdr = zeros(Float64, nd, nr)  # total households by district, by religion
        tpbdr = zeros(Float64, nd, nr)  # total members of households by district, by religion
        for i=1:nh
            thbdr[indDis[i],indRel[i]] += 1
            if sqrRoot; tpbdr[indDis[i],indRel[i]] += sqrt(siz[hhid[i]])
            else tpbdr[indDis[i],indRel[i]] += siz[hhid[i]]
            end
        end
    end

    # calculate average monthly expenditure per capita by districts
    totexp = zeros(Float64, nd)     # total expenditures by district
    for i=1:nh; totexp[indDis[i]] += inc[hhid[i]]*siz[hhid[i]] end
    if sqrRoot; for i=1:nd; ave[disList[i]] = totexp[i]/sqrt(tpbd[i]) end
    else for i=1:nd; ave[disList[i]] = totexp[i]/tpbd[i] end
    end
    # convert 'AVEpC' to annual or daily
    if period=="annual"; mmtoyy = 365/30; for i=1:nd; ave[disList[i]] = ave[disList[i]] * mmtoyy end
    elseif period=="daily"; for i=1:nd; ave[disList[i]] = ave[disList[i]] / 30 end
    end

    # categorize emission data
    ec = emissionsHHs[year]
    ed = zeros(Float64, nd, nc)
    if !sqrRoot; for i=1:nh; ed[indDis[i],:] += ec[i,:] end
    elseif sqrRoot && weightMode==5; for i=1:nh; ed[indDis[i],:] += ec[i,:]/sqrt(siz[hhid[i]]) end
    end

    # weighting
    if weightMode == 1
        for i=1:nd
            if haskey(pop, disList[i]); ed[i,:] *= pop[disList[i]][1]/tpbd[i]
            else ed[i,:] .= 0
            end
        end
    elseif weightMode == 2
        for i=1:nd
            if haskey(pop, disList[i]); ed[i,:] *= pop[disList[i]][2]/thbd[i]
            else ed[i,:] .= 0
            end
        end
    elseif weightMode == 3
        for i=1:nd
            if haskey(pop, disList[i]); ed[i,:] *= 0.5*(pop[disList[i]][1]/tpbd[i]+pop[disList[i]][2]/thbd[i])
            else ed[i,:] .= 0
            end
        end
    # normalizing
    elseif weightMode == 4; for i=1:nc; ed[:,i] ./= tpbd end
    elseif weightMode == 5; for i=1:nc; ed[:,i] ./= thbd end
    # basic information
    elseif weightMode == 6; ed[:,1], ed[:,2] = tpbd[:], thbd[:]
    end

    # calculate differences
    avg = mean(ed, dims=1)
    ecd = zeros(size(ed))
    for i=1:size(ed,2); ecd[:,i] = (ed[:,i].-avg[i])/avg[i] end

    # save the results
    emissionsDis[year] = ed
    emissionsDisDiff[year] = ecd

    if religion; return ed, catList, disList, thbd, tpbd, thbdr, tpbdr
    else return ed, catList, disList, thbd, tpbd
    end
end

function categorizeDistrictByEmissionLevel(year, normMode = 0, intv=[])
                                            # intv: proportions between invervals of highest to lowest
                                            # normmode: [1]per capita, [2]per household, [3]basic information
    global hhid, sam, catList, disList
    global emissionsDis
    ed = emissionsDis[year]

    if length(intv) == 0; intv = [0.25,0.25,0.25,0.25]
    elseif sum(intv) != 1; sumi = sum(intv); intv /= sumi
    end

    nh = length(hhid)
    nc = length(catList)
    nd = length(disList)
    ni = length(intv)

    # make index dictionaries
    disOrder = sortperm(ed[:,end], rev=true)    #descending order indexing, [1]highest, [end]lowest values' indexes
    indDis = zeros(Int, nd)     # index dictionary of districts
    i = 1
    for s = 1:ni
        while i <= trunc(Int, nd*sum(intv[1:s]))
            indDis[disOrder[i]] = s
            i += 1
        end
    end
    indDis[disOrder[nd]] = ni   # for the last index

    # categorize emission data
    edl = zeros(Float64, ni, nc)
    tp = zeros(Int, ni)
    th = zeros(Int, ni)
    for i=1:nd
        edl[indDis[i],:] += ed[i,:]
        tp[indDis[i]] += sam[disList[i]][1]
        th[indDis[i]] += sam[disList[i]][2]
    end
    # normalizing
    if normMode == 1; for i=1:nc; edl[:,i] ./= tp end
    elseif normMode == 2 ;for i=1:nc; edl[:,i] ./= th end
    # basic information
    elseif normMode == 3; ed[:,1], ed[:,2] = tp[:], th[:]
    end

    emissionsDisLev[year] = edl

    return edl, catList, disList
end

function categorizeHouseholdByReligion(year, normMode=0; sqrRt=false, popWgh=false)
                                            # normmode: [1]per capita, [2]per houehold, [3]basic information
    global hhid, cat, dis, siz, rel, catList, relList
    global emissionsHHs, emissionsRel

    nh = length(hhid)
    nc = length(catList)
    nr = length(relList)

    # make an index dictionarie array of religion
    indRel = [findfirst(x->x==rel[hhid[i]], relList) for i=1:nh]

    # sum households and members by districts
    thbr = zeros(Float64, nr)   # total households by religion
    tpbr = zeros(Float64, nr)   # total members of households by religion
    twpbr = zeros(Float64, nr)  # total state-population weighted members of households by religion
    for i=1:nh; thbr[indRel[i]] += 1 end
    if !sqrRt; for i=1:nh; tpbr[indRel[i]] += siz[hhid[i]] end
    elseif sqrRt && normMode==2; for i=1:nh; tpbr[indRel[i]] += sqrt(siz[hhid[i]]) end
    end
    if popWgh
        if sqrRt; for i=1:nh; twpbr[indRel[i]] += sqrt(siz[hhid[i]]) * wgh[hhid[i]] end
        else for i=1:nh; twpbr[indRel[i]] += siz[hhid[i]] * wgh[hhid[i]] end
        end
    end

    # categorize emission data
    ec = emissionsHHs[year]
    er = zeros(Float64, nr, nc)
    if !sqrRt
        if popWgh; for i=1:nh; er[indRel[i],:] += ec[i,:] * wgh[hhid[i]] end
        else for i=1:nh; er[indRel[i],:] += ec[i,:] end
        end
    elseif sqrRt && normMode==2; for i=1:nh; er[indRel[i],:] += ec[i,:]/sqrt(siz[hhid[i]]) end
    end

    # normalizing
    if normMode == 1
        if popWgh; for i=1:nc; er[:,i] ./= twpbr end
        else for i=1:nc; er[:,i] ./= tpbr end
        end
    elseif normMode == 2; for i=1:nc; er[:,i] ./= thbr end
    # basic information
    elseif normMode == 3; er[:,1], er[:,2] = tpbr[:], thbr[:]
    end

    emissionsRel[year] = er

    return er, catList, relList, tpbr, thbr, twpbr
end

function categorizeHouseholdByIncome(year,intv=[],normMode=0; sqrRt=false,absIntv=false,perCap=false,desOrd=false,popWgh=false)
                                            # intv: proportions between invervals of highest to lowest
                                            # absIntv: if "true", then intv[] is a list of income values, descending order
                                            # normmode: [1]per capita, [2]per houehold
                                            # perCap: [true] per capita mode, [false] per household mode
                                            # desOrd: [true] descening order of 'intv[]', [false] ascending order
    global sec, hhid, cat, dis, siz, inc, catList, incList
    global emissionsHHs, emissionsInc

    if !absIntv && length(intv) == 0; intv = [0.25,0.5,0.75,1.00]
    elseif sort!(intv)[end] != 1; intv /= intv[end]
    end

    nh = length(hhid)
    nc = length(catList)
    ni = length(intv); if absIntv==true; ni +=1 end

    incArray = [inc[h] for h in hhid]
    incOrder = sortperm(incArray, rev=desOrd)   # desOrd: [true] descening order, [false] ascending order

    # determine sections' starting index and values
    pcidx = []  # current sector's starting index for 'per capita' emissions
    if perCap && !absIntv   # determine sections if interval ratios are for 'per capita' emissions
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
            if accpop/totpop > intv[idx]
                push!(pcidx, i)
                push!(incList, incArray[incOrder[i]])
                idx += 1
            end
        end
        if !(nh in pcidx); push!(pcidx, nh); push!(incList, incArray[incOrder[nh]]) end
    elseif absIntv; incList = copy(intv)
    else
        push!(incList, incArray[incOrder[1]])
        for i=1:length(intv); push!(incList, incArray[incOrder[trunc(Int, intv[i]*nh)]]) end
    end

    # set sector index
    indInc = zeros(Int, nh)
    if absIntv  # indInc: '1' for over intv[1], '3' for below intv[2], [2] for others
        for i=1:nh
            if incArray[i] >= intv[1]; indInc[i] = 1
            elseif incArray[i] < intv[2]; indInc[i] = 3
            else indInc[i] = 2
            end
        end
    elseif perCap
        idx = 1
        for i=1:nh
            if idx<ni && i>=pcidx[idx+1]; idx +=1 end
            indInc[incOrder[i]] = idx
        end
    else
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
    twpbi = zeros(Float64, ni)  # total state-population weighted members of households by income level
    for i= 1:nh; thbi[indInc[i]] += 1 end
    if sqrRt; for i= 1:nh; tpbi[indInc[i]] += sqrt(siz[hhid[i]]) end
    else for i= 1:nh; tpbi[indInc[i]] += siz[hhid[i]] end
    end
    if popWgh
        if sqrRt; for i= 1:nh; twpbi[indInc[i]] += sqrt(siz[hhid[i]]) * wgh[hhid[i]] end
        else for i= 1:nh; twpbi[indInc[i]] += siz[hhid[i]] * wgh[hhid[i]] end
        end
    end

    # categorize emission data
    ec = emissionsHHs[year]
    ei = zeros(Float64, ni, nc)
    if !sqrRt
        if popWgh; for i=1:nh; ei[indInc[i],:] += ec[i,:] * wgh[hhid[i]] end
        else for i=1:nh; ei[indInc[i],:] += ec[i,:] end
        end
    elseif sqrRt && normMode==2; for i=1:nh; ei[indInc[i],:] += ec[i,:]/sqrt(siz[hhid[i]]) end
    end

    # normalizing
    if normMode == 1
        if popWgh; for i=1:nc; ei[:,i] ./= twpbi end
        else for i=1:nc; ei[:,i] ./= tpbi end
        end
    elseif normMode == 2; for i=1:nc; ei[:,i] ./= thbi end
    end

    emissionsInc[year] = ei

    return ei, catList, incList, tpbi, thbi, twpbi
end

function categorizeHouseholdByExpRange(year,rng=[],normMode=0; absRng=false,perCap=false,popWgh=false,over=0.1,less=0.1)
                                            # absRng: [true] apply absolute range, [false] apply population ratio range
                                            # rng: standard values of ranges for grouping
                                            # over/less: range from 'stdv', ratios of samples, househods, or population
                                            # normMode: [1]per capita, [2]per houehold
                                            # perCap: [true] per capita mode, [false] per household mode
    global sec, hhid, cat, dis, siz, inc, catList, incList
    global emissionsHHs, emissionsRng

    sort!(rng)

    nh = length(hhid)
    nc = length(catList)
    nr = length(rng)

    incList = rng
    expArray = [inc[h] for h in hhid]
    expOrder = sortperm(expArray, rev=false)    # rev: [true] descening order, [false] ascending order

    # determine each section's starting and ending indexes and values
    rsidx = zeros(Int, nr, 3)   # sector's standard, starting, and ending index: [standard, start, end]
    if absRng
        for i=1:nr
            rsidx[i,1] = findmin(abs.([expArray[x]-rng[i] for x in expOrder]))[2]     # [value, index]
            rsidx[i,2] = findmin(abs.([expArray[x]-(rng[i]*(1-less)) for x in expOrder]))[2]     # [value, index]
            rsidx[i,3] = findmin(abs.([expArray[x]-(rng[i]*(1+over)) for x in expOrder]))[2]     # [value, index]
        end
    else
        if perCap       # determine sections if ranges are for 'per capita' emissions
            if popWgh; totpop = sum([siz[h]*wgh[h] for h in hhid])
            else totpop = sum(collect(values(siz)))
            end
            for i=1:nr
                # find standard index of i_th range
                rsidx[i,1] = findmin(abs.([expArray[x]-rng[i] for x in expOrder]))   # [value, index]
                # find bottom index of i_th range
                idx = rsidx[i,1]
                lp = totpop * less
                while lp > 0 && idx > 0
                    h = hhid[expOrder[idx]]
                    if popWgh; lp -= siz[h] * wgh[h] else lp -= siz[h] end
                    if lp>0; rsidx[i,2] = idx end
                    idx -= 1
                end
                # find top index of i_th range
                idx = rsidx[i,1]
                lp = totpop * over
                while lp > 0 && idx <= nh
                    h = hhid[expOrder[idx]]
                    if popWgh; lp -= siz[h] * wgh[h] else lp -= siz[h] end
                    if lp>0; rsidx[i,3]= idx end
                    idx += 1
                end
            end
        else
            for i=1:nr
                rsidx[i,1] = findmin(abs.([expArray[x]-rng[i] for x in expOrder]))   # [value, index]
                rsidx[i,2] = rsidx[i,1] - (less * nh); if rsidx[i,2] <= 0; rsidx[i,2] = 1 end
                rsidx[i,3] = rsidx[i,1] + (over * nh); if rsidx[i,2] > nh; rsidx[i,2] = nh end
            end
        end
    end

    # sum households and members by districts
    thber = zeros(Float64, nr)   # total households by expenditure range
    tpber = zeros(Float64, nr)   # total members of households by expenditure range
    twpber = zeros(Float64, nr)  # total state-population weighted members of households by expenditure range
    for i=1:nr; thber[i] = rsidx[i,3] - rsidx[i,2] + 1 end
    for i=1:nr; for j=rsidx[i,2]:rsidx[i,3]; tpber[i] += siz[hhid[expOrder[j]]] end end
    if popWgh;
        for i=1:nr; for j=rsidx[i,2]:rsidx[i,3]; twpber[i] += siz[hhid[expOrder[j]]] * wgh[hhid[expOrder[j]]] end end
    end

    # categorize emission data
    ec = emissionsHHs[year]
    er = zeros(Float64, nr, nc)
    if popWgh; for i=1:nr; for j=rsidx[i,2]:rsidx[i,3]; er[i,:] += ec[expOrder[j],:] * wgh[hhid[expOrder[j]]] end end
    else for i=1:nr; for j=rsidx[i,2]:rsidx[i,3]; er[i,:] += ec[expOrder[j],:] end end
    end

    # normalizing
    if normMode == 1
        if popWgh; for i=1:nc; er[:,i] ./= twpber end
        else for i=1:nc; er[:,i] ./= tpber end
        end
    elseif normMode == 2; for i=1:nc; er[:,i] ./= thber end
    end

    emissionsRng[year] = er

    return er, catList, thber, tpber, twpber, rsidx, expOrder
end

function categorizeHouseholdByEmissionLevel(year, intv=[], normMode = 0; squareRoot = false, absintv=false)
                                                    # intv: proportions between invervals of highest to lowest
                                                    # normmode: [1]per capita, [2]per houehold
    global hhid, sec, catList, levList, cat, dis, siz
    global emissionsHHs, emissionsLev

    if !absintv && length(intv) == 0; intv = [0.25,0.25,0.25,0.25]
    elseif !absintv && sum(intv) != 1; sumi = sum(intv); intv /= sumi
    end

    ec = emissionsHHs[year]
    nh = length(hhid)
    nc = length(catList)
    nl = length(intv); if absintv==true; nl +=1 end

    # make index list
    if squareRoot; levArray = [ec[i,nc]/sqrt(siz[hhid[i]]) for i=1:nh]
    elseif normMode==1; levArray = [ec[i,nc]/siz[hhid[i]] for i=1:nh]
    elseif normMode==2; levArray = copy(ec[:,nc])
    end
    levOrder = sortperm(levArray, rev=true)     #descending order indexing, [1]highest, [end]lowest values' indexes
    for i=1:nl; push!(levList, levArray[levOrder[trunc(Int, sum(intv[1:i])*nh)]]) end

    indLev = zeros(Int, nh)     # index dictionary of income sections
    i = 1
    for s = 1:nl
        while i <= trunc(Int, nh*sum(intv[1:s]))
            indLev[levOrder[i]] = s
            i += 1
        end
    end
    indLev[levOrder[nh]] = nl

    # sum households and members by districts
    thbl = zeros(Int, nl)   # total households by income level
    tpbl = zeros(Int, nl)   # total members of households by income level
    for i=1:nh
        thbl[indLev[i]] += 1
        if !squareRoot; tpbl[indLev[i]] += siz[hhid[i]]
        elseif squareRoot && normMode==2; tpbl[indLev[i]] += sqrt(siz[hhid[i]])
        end
    end

    # categorize emission data
    ec = emissionsHHs[year]
    el = zeros(Float64, nl, nc)
    if !squareRoot; for i=1:nh; el[indLev[i],:] += ec[i,:] end
    elseif squareRoot && normMode==2; for i=1:nh; el[indLev[i],:] += ec[i,:]/sqrt(siz[hhid[i]]) end
    end

    # normalizing
    if normMode == 1; for i=1:nc; el[:,i] ./= tpbl end
    elseif normMode == 2; for i=1:nc; el[:,i] ./= thbl end
    # basic information
    elseif normMode == 3; el[:,1], el[:,2] = tpbl[:], thbl[:]
    end

    emissionsLev[year] = el
end

function categorizeHouseholdByIncomeByReligion(year,intv=[],normMode=0; sqrRt=false,absIntv=false,perCap=false,desOrd=false,popWgh=false)
                                            # intv: proportions between invervals of highest to lowest
                                            # absIntv: if "true", then intv[] is a list of income values, descending order
                                            # normmode: [1]per capita, [2]per houehold, [3]basic information
                                            # perCap: [true] per capita mode, [false] per household mode
                                            # desOrd: [true] descening order of 'intv[]', [false] ascending order
    global sec, hhid, cat, dis, siz, inc, rel, catList, incList, relList
    global emissionsHHs, emissionsIncRel

    if !absIntv && length(intv) == 0; intv = [0.25,0.5,0.75,1.00]
    elseif sort!(intv)[end] != 1; intv /= intv[end]
    end

    nh = length(hhid)
    nc = length(catList)
    nr = length(relList)
    ni = length(intv); if absIntv==true; ni +=1 end

    incArray = [inc[h] for h in hhid]
    incOrder = sortperm(incArray, rev=desOrd)   # desOrd: [true] descening order, [false] ascending order

    # make index dictionarie of religion
    indRel = [findfirst(x->x==rel[hhid[i]], relList) for i=1:nh]

    # determine sections' starting index and values of income index
    pcidx = []  # current sector's starting index for 'per capita' emissions
    if perCap && !absIntv   # determine sections if interval ratios are for 'per capita' emissions
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
            if accpop/totpop > intv[idx]
                push!(pcidx, i)
                push!(incList, incArray[incOrder[i]])
                idx += 1
            end
        end
        if !(nh in pcidx); push!(pcidx, nh); push!(incList, incArray[incOrder[nh]]) end
    elseif absIntv; incList = copy(intv)
    else
        push!(incList, incArray[incOrder[1]])
        for i=1:length(intv); push!(incList, incArray[incOrder[trunc(Int, intv[i]*nh)]]) end
    end

    # set sector index
    indInc = zeros(Int, nh)
    if absIntv  # indInc: '1' for over intv[1], '3' for below intv[2], [2] for others
        for i=1:nh
            if incArray[i] >= intv[1]; indInc[i] = 1
            elseif incArray[i] < intv[2]; indInc[i] = 3
            else indInc[i] = 2
            end
        end
    elseif perCap
        idx = 1
        for i=1:nh
            if idx<ni && i>=pcidx[idx+1]; idx +=1 end
            indInc[incOrder[i]] = idx
        end
    else
        i = 1
        for s = 1:ni
            while i <= trunc(Int, nh*intv[s])
                indInc[incOrder[i]] = s
                i += 1
            end
        end
        indInc[incOrder[nh]] = ni
    end

    # sum households and members by religion and income
    thbir = zeros(Float64, nr, ni)   # total households by income level
    tpbir = zeros(Float64, nr, ni)   # total members of households by income level
    twpbir = zeros(Float64, nr, ni)  # total state-population weighted members of households by income level
    for i= 1:nh; thbir[indRel[i],indInc[i]] += 1 end
    if sqrRt; for i= 1:nh; tpbir[indRel[i],indInc[i]] += sqrt(siz[hhid[i]]) end
    else for i= 1:nh; tpbir[indRel[i],indInc[i]] += siz[hhid[i]] end
    end
    if popWgh
        if sqrRt; for i= 1:nh; twpbir[indRel[i],indInc[i]] += sqrt(siz[hhid[i]]) * wgh[hhid[i]] end
        else for i= 1:nh; twpbir[indRel[i],indInc[i]] += siz[hhid[i]] * wgh[hhid[i]] end
        end
    end

    # categorize emission data
    ec = emissionsHHs[year]
    eir = zeros(Float64, nr, ni, nc)
    if !sqrRt
        if popWgh; for i=1:nh; eir[indRel[i],indInc[i],:] += ec[i,:] * wgh[hhid[i]] end
        else for i=1:nh; eir[indRel[i],indInc[i],:] += ec[i,:] end
        end
    elseif sqrRt && normMode==2; for i=1:nh; eir[indRel[i],indInc[i],:] += ec[i,:]/sqrt(siz[hhid[i]]) end
    end

    # normalizing
    if normMode == 1
        if popWgh; for i=1:nc; eir[:,:,i] ./= twpbir end
        else for i=1:nc; eir[:,:,i] ./= tpbir end
        end
    elseif normMode == 2; for i=1:nc; eir[:,:,i] ./= thbir end
    # basic information
    elseif normMode == 3; eir[:,:,1], eir[:,:,2] = tpbir, thbir
    end

    emissionsIncRel[year] = eir

    return eir, catList, incList, tpbir, thbir, twpbir
end

function printEmissionByDistrict(year,outputFile,tpbdr=[],thbdr=[]; name=false,expm=false,popm=false,hhsm=false,relm=false)
    # expm: print average expenditure per capita
    # popm: print population related figures
    # hhsm: print households related figures
    # relm: print relgion related figures

    global nam, catList, disList, relList, relName, pop, ave, emissionsDis
    ed = emissionsDis[year]

    f = open(outputFile, "w")

    print(f,"District")
    for c in catList; print(f, ",", c) end
    if expm; print(f, ",Exp") end
    if popm; print(f, ",Pop") end
    if hhsm; print(f, ",HHs") end
    if relm && popm; print(f, ",RelByPop")  end
    if relm && hhsm; print(f, ",RelByHHs")  end
    if relm && !popm && !hhsm; println("Religion mode should be operated with 'pop' or 'hhs' mode.") end
    println(f)
    for i = 1:length(disList)
        if name; print(f, nam[disList[i]]) else print(f, disList[i]) end
        for j = 1:length(catList); print(f, ",", ed[i,j]) end
        if expm; print(f, ",", ave[disList[i]]) end
        if popm; print(f, ",", pop[disList[i]][1]) end
        if hhsm; print(f, ",", pop[disList[i]][2]) end
        if relm && popm; print(f, ",", relName[partialsortperm(tpbdr[i,:],1,rev=true)])  end
        if relm && hhsm; print(f, ",", relName[partialsortperm(thbdr[i,:],1,rev=true)])  end
        println(f)
    end

    close(f)
end

function printEmissionByReligion(year, outputFile, tpbr=[], thbr=[], twpbr=[])

    global catList, relList, relName, emissionsRel
    er = emissionsRel[year]

    f = open(outputFile, "w")

    print(f,"Religion")
    for c in catList; print(f, ",", c) end
    if length(tpbr)>0; print(f,",Pop.") end
    if length(thbr)>0; print(f,",HH.") end
    if length(twpbr)>0; print(f,",WghPop.") end
    println(f)
    for i = 1:length(relList)
        print(f, relName[i])
        for j = 1:length(catList); print(f, ",", er[i,j]) end
        if length(tpbr)>0; print(f,",",tpbr[i]) end
        if length(thbr)>0; print(f,",",thbr[i]) end
        if length(twpbr)>0; print(f,",",twpbr[i]) end
        println(f)
    end

    close(f)
end

function printEmissionByIncome(year, outputFile, intv=[], tpbi=[], thbi=[], twpbi=[]; absIntv=false, desOrd=false)

    global catList, incList, emissionsInc
    ei = emissionsInc[year]
    ni = length(intv); if absIntv; ni += 1 end

    f = open(outputFile, "w")

    print(f,"Exp_Lv")
    for c in catList; print(f, ",", c) end
    if length(tpbi)>0; print(f,",Pop.") end
    if length(thbi)>0; print(f,",HH.") end
    if length(twpbi)>0; print(f,",WghPop.") end
    println(f)
    for i = 1:ni
        if absIntv
            if i==1; print(f, "< ",intv[1])
            elseif i==2; print(f, "< ",intv[2])
            elseif i==3; print(f, "> ",intv[2])
            end
        elseif i==1; print(f, "0-", round(intv[1]*100,digits=1),"% (",incList[i+1],")")
        else print(f, round(intv[i-1]*100,digits=1),"-",round(intv[i]*100,digits=1),"% (",incList[i+1],")")
        end

        for j = 1:length(catList); print(f, ",", ei[i,j]) end
        if length(tpbi)>0; print(f,",",tpbi[i]) end
        if length(thbi)>0; print(f,",",thbi[i]) end
        if length(twpbi)>0; print(f,",",twpbi[i]) end
        println(f)
    end

    close(f)
end

function printEmissionByRange(year, outputFile, rsidx=[], thber=[], tpber=[], twpber=[], order=[])

    global catList, emissionsRng
    er = emissionsRng[year]
    nr = size(rsidx,1)

    f = open(outputFile, "w")

    print(f,"Expenditure")
    for c in catList; print(f, ",", c) end
    if length(tpber)>0; print(f,",Pop.") end
    if length(thber)>0; print(f,",HH.") end
    if length(twpber)>0; print(f,",WghPop.") end
    println(f)
    for i = 1:nr
        print(f,round(inc[hhid[order[rsidx[i,2]]]],digits=2),"-",round(inc[hhid[order[rsidx[i,3]]]],digits=2))

        for j = 1:length(catList); print(f, ",", er[i,j]) end
        if length(tpber)>0; print(f,",",tpber[i]) end
        if length(thber)>0; print(f,",",thber[i]) end
        if length(twpber)>0; print(f,",",twpber[i]) end
        println(f)
    end

    close(f)
end

function printEmissionByIncomeByReligion(year, outputFile, intv=[], tpbir=[], thbir=[], twpbir=[]; absIntv=false, desOrd=false)

    global catList, incList, relList, emissionsIncRel
    eir = emissionsIncRel[year]
    nc = length(catList)
    nr = length(relList)
    ni = length(intv); if absIntv; ni += 1 end
    relName = ["Hinduism","Islam","Christianity","Sikhism","Jainism","Buddhism","Zoroastrianism", "Others", "None"]

    f = open(outputFile, "w")

    for i = 1:nr
        println(f, relName[i])
        print(f,"Exp_Lv")
        for c in catList; print(f, ",", c) end
        if length(tpbir)>0; print(f,",Pop.") end
        if length(thbir)>0; print(f,",HH.") end
        if length(twpbir)>0; print(f,",WghPop.") end
        println(f)
        for j = 1:ni
            if absIntv
                if j==1; print(f, "< ",intv[1])
                elseif j==2; print(f, "< ",intv[2])
                elseif j==3; print(f, "> ",intv[2])
                end
            elseif j==1; print(f, "0-", round(intv[1]*100,digits=1),"% (",incList[j+1],")")
            else print(f, round(intv[j-1]*100,digits=1),"-",round(intv[j]*100,digits=1),"% (",incList[j+1],")")
            end

            for k = 1:nc; print(f, ",", eir[i,j,k]) end
            if length(tpbir)>0; print(f,",",tpbir[i,j]) end
            if length(thbir)>0; print(f,",",thbir[i,j]) end
            if length(twpbir)>0; print(f,",",twpbir[i,j]) end
            println(f)
        end
        println(f)
    end
    close(f)
end

#=
function plotHHsEmission(year)

    ec = emissionsHHs[year]

    plotly()

    e = ec[:,end]

    p = violin(e, y="CO2 emission (tCO2)", box=true, points='all')

    if dispmode; display(p) end
    if guimode; gui() end

end

function plotEmissionByIncome(year, intv=[]; desOrd=false)

    plotly()
    e = emissionsInc[year]

end
=#
function printEmissionByDistEmLev(year, outputFile, intv=[])

    global catList, disList, emissionsDisLev
    edl = emissionsDisLev[year]

    f = open(outputFile, "w")

    print(f,"CF_Lv")
    for c in catList; print(f, ",", c) end
    println(f)
    for i = 1:length(intv)
        print(f, "\t<", trunc(Int, sum(intv[1:i])*100),"%")
        for j = 1:length(catList); print(f, ",", edl[i,j]) end
        println(f)
    end

    close(f)
end

function printEmissionByHhsEmLev(year, outputFile, intv=[])

    global catList, disList, emissionsLev
    el = emissionsLev[year]

    f = open(outputFile, "w")

    print(f,"CF_Lv")
    for c in catList; print(f, ",", c) end
    println(f)
    for i = 1:length(intv)
        print(f, "\t<", trunc(Int, sum(intv[1:i])*100),"%")
        for j = 1:length(catList); print(f, ",", el[i,j]) end
        println(f)
    end

    close(f)
end

function exportDistrictEmission(year, tag, outputFile, weightMode::Int; name=false)

    global sam, pop, ave, gid, gidData, catList, disList, emissionsCatNW
    ec = emissionsDis[year]
    nc = length(catList)

    # making exporting table
    gidList = sort(unique(values(gid)))
    ngid = length(gidList)

    tb = zeros(Float64, ngid, nc)
    spo = zeros(Float64, ngid)   # number of sample population by district
    shh = zeros(Float64, ngid)   # number of sample households by district
    tpo = zeros(Float64, ngid)   # total number of population by district
    thh = zeros(Float64, ngid)   # total number of households by district
    aec = zeros(Float64, ngid)   # average expenditure per capita by district
    for i=1:length(disList)
        idx = findfirst(x->x==gid[disList[i]],gidList)
        for j=1:nc; tb[idx,j] += ec[i,j] end
        spo[idx] += sam[disList[i]][1]
        shh[idx] += sam[disList[i]][2]
        tpo[idx] += pop[disList[i]][1]
        thh[idx] += pop[disList[i]][2]
        aec[idx] += ave[disList[i]]*sam[disList[i]][1]
    end
    for i=1:ngid; aec[i] /= spo[i] end
    if weightMode==1; for i=1:ngid; for j=1:nc; tb[i,j] *= tpo[i]/spo[i] end end
    elseif weightMode==2; for i=1:ngid; for j=1:nc; tb[i,j] *= thh[i]/shh[i] end end
    elseif weightMode==4; for i=1:ngid; for j=1:nc; tb[i,j] /= spo[i] end end
    elseif weightMode==5; for i=1:ngid; for j=1:nc; tb[i,j] /= shh[i] end end
    end

    # exporting table
    f = open(outputFile, "w")
    print(f, tag); for c in catList; print(f,",",c) end; println(f)
    for i = 1:size(tb, 1)
        if name; print(f, gidData[gidList[i]][2]) else print(f, gidList[i]) end
        for j = 1:size(tb, 2); print(f, ",", tb[i,j]) end
        println(f)
    end
    close(f)

    gisDistrictEmission[year] = tb

    return tb, gidList, spo, shh, tpo, thh, aec
end

function exportEmissionDiffRate(year,tag,outputFile,maxr=0.5,minr=-0.5,nspan=128; descend=false,name=false,empty=false)

    global gid, gidData, misDist
    gde = gisDistrictEmission[year]
    gidList = sort(unique(values(gid)))

    # calculate difference rates
    avg = mean(gde, dims=1)
    gded = zeros(size(gde))
    for i=1:size(gde,2); gded[:,i] = (gde[:,i].-avg[i])/avg[i] end

    # grouping by ratios; ascending order
    span = [(maxr-minr)*(i-1)/(nspan-2)+minr for i=1:nspan-1]
    rank = zeros(Int, size(gded))
    for j=1:size(gded,2)    # category number
        for i=1:size(gded,1)    # gid district number
            if gded[i,j]>=maxr; rank[i,j] = nspan
            else rank[i,j] = findfirst(x->x>gded[i,j],span)
            end
        end
    end
    # for descending order, if "descend == true".
    if descend; for j=1:size(gded,2); for i=1:size(gded,1); rank[i,j] = nspan - rank[i,j] + 1 end end end

    # exporting difference table
    f = open(outputFile, "w")
    print(f, tag)
    for c in catList; print(f,",",c) end
    println(f)
    for i = 1:size(gded, 1)
        if name; print(f, gidData[gidList[i]][2])
        else print(f, gidList[i])
        end
        for j = 1:size(gded, 2); print(f, ",", gded[i,j]) end
        println(f)
    end
    if empty
        for i=1:length(misDist)
            if name; print(f, gidData[misDist[i]][2])
            else print(f, misDist[i])
            end
            for j = 1:size(gded, 2); print(f, ",0") end
            println(f)
        end
    end

    close(f)

    # exporting difference group table
    f = open(replace(outputFile,".csv"=>"_gr.csv"), "w")
    print(f, tag)
    for c in catList; print(f,",",c) end
    println(f)
    for i = 1:size(rank, 1)
        if name; print(f, gidData[gidList[i]][2])
        else print(f, gidList[i])
        end
        for j = 1:size(rank, 2); print(f, ",", rank[i,j]) end
        println(f)
    end
    if empty
        for i=1:length(misDist)
            if name; print(f, gidData[misDist[i]][2])
            else print(f, misDist[i])
            end
            for j = 1:size(gded, 2); print(f, ",0") end
            println(f)
        end
    end
    close(f)

    gisDistrictEmissionDiff[year] = gded
    gisDistrictEmissionDiffRank[year] = rank

    return gded, avg, rank
end

function exportEmissionValGroup(year, tag, outputFile, nspan=128; max=0, min=0, descend=false, logscl=false, name=false)

    global gid, gidData
    gde = gisDistrictEmission[year]
    gidList = sort(unique(values(gid)))
    if max==min==0; setNode = true; else setNode = false end

    # grouping by emission amount; ascending order
    rank = zeros(Int, size(gde))
    for j=1:size(gde,2)    # category number
        if setNode; max = maximum(gde[:,j]); min = minimum(gde[:,j]) end
        if logscl; max = log(max); if min>0; min = log(min) end end
        span = [(max-min)*i/nspan + min for i=1:nspan]
        if logscl
            for i=1:size(gde,1)    # gid district number
                if log(gde[i,j])==max; rank[i,j] = nspan
                else rank[i,j] = findfirst(x->x>log(gde[i,j]),span)
                end
            end
        else
            for i=1:size(gde,1)    # gid district number
                if gde[i,j]==max; rank[i,j] = nspan
                else rank[i,j] = findfirst(x->x>gde[i,j],span)
                end
            end
        end

    end
    # for descending order, if "descend == true".
    if descend; for j=1:size(gde,2); for i=1:size(gde,1); rank[i,j] = nspan - rank[i,j] + 1 end end end

    # exporting difference table
    f = open(outputFile, "w")
    print(f, tag)
    for c in catList; print(f,",",c) end
    println(f)
    for i = 1:size(gde, 1)
        if name; print(f, gidData[gidList[i]][2])
        else print(f, gidList[i])
        end
        for j = 1:size(gde, 2); print(f, ",", gde[i,j]) end
        println(f)
    end
    close(f)

    # exporting difference group table
    f = open(replace(outputFile,".csv"=>"_gr.csv"), "w")
    print(f, tag)
    for c in catList; print(f,",",c) end
    println(f)
    for i = 1:size(rank, 1)
        if name; print(f, gidData[gidList[i]][2])
        else print(f, gidList[i])
        end
        for j = 1:size(rank, 2); print(f, ",", rank[i,j]) end
        println(f)
    end
    close(f)

    return rank
end

function exportEmissionRankGroup(year, tag, outputFile, nspan = 128; descend=false, name=false)

    global gid, gidData
    gde = gisDistrictEmission[year]
    gidList = sort(unique(values(gid)))
    nd = length(gidList)

    # grouping by emission amount; ascending order
    rank = zeros(Int, size(gde))
    for j=1:size(gde,2)    # category number
        order = sortperm(gde[:,j], rev=descend)
        for i=1:nd; rank[order[i],j] = ceil(Int, i/nd*nspan) end
    end

    # exporting group table
    f = open(outputFile, "w")
    print(f, tag)
    for c in catList; print(f,",",c) end
    println(f)
    for i = 1:size(rank, 1)
        if name; print(f, gidData[gidList[i]][2])
        else print(f, gidList[i])
        end
        for j = 1:size(rank, 2); print(f, ",", rank[i,j]) end
        println(f)
    end
    close(f)

    return rank
end

function exportWebsiteFiles(year, path, weightMode, gidList, totalPop, totalHH, sampPop, avgExp; rank=false, empty=false)

    global nam, gid, gidData, catList, misDist
    global gisEmissionCat, gisEmissionCatDif
    gde = gisDistrictEmission[year]
    gded = gisDistrictEmissionDiff[year]
    gdedr = gisDistrictEmissionDiffRank[year]

    for j=1:length(catList)
        f = open(path*"CFAC_"*catList[j]*"_"*string(year)*".txt","w")
        println(f, "KEY_CODE\tSTATE\tDISTRICT\tSTATE_NAME\tDISTRICT_NAME\tGENHH_ALLPPC")
        for i=1:length(gidList)
            gd = gidData[gidList[i]]
            print(f, gidList[i],"\t",gd[3],"\t",gd[1],"\t",gd[4],"\t",gd[2],"\t")
            if rank; print(f, gdedr[i,j])
            else print(f, gded[i,j])
            end
            println(f)
        end
        if empty
            for g in misDist
                gd = gidData[g]
                print(f, g,"\t",gd[3],"\t",gd[1],"\t",gd[4],"\t",gd[2],"\t")
                if rank; print(f,"0") end
                println(f)
            end
        end
        close(f)

        f = open(path*"CFAV_"*catList[j]*"_"*string(year)*".txt","w")
        print(f, "KEY_CODE\tSTATE\tDISTRICT\tSTATE_NAME\tDISTRICT_NAME\tGENHH_ALLP\tGENHH_APPPC")
        if catList[j]=="Total" || catList[j]=="All"; println(f, "\tANEXPPC\tPOP")
        else println(f)
        end
        for i=1:length(gidList)
            gd = gidData[gidList[i]]
            print(f, gidList[i],"\t",gd[3],"\t",gd[1],"\t",gd[4],"\t",gd[2],"\t")
            if weightMode == 1; print(f,gde[i,j],"\t",gde[i,j]/totalPop[i])
            elseif weightMode == 2; print(f,gde[i,j],"\t",gde[i,j]/totalHH[i])
            elseif weightMode == 4; print(f,gde[i,j]*totalPop[i],"\t",gde[i,j])
            elseif weightMode == 5; print(f,gde[i,j]*totalHH[i],"\t",gde[i,j])
            elseif weightMode == 0; print(f,gde[i,j]*totalHH[i]/sampPop[i],"\t",gde[i,j]/sampPop[i])
            end
            if catList[j]=="Total" || catList[j]=="All"; println(f,"\t",avgExp[i],"\t",convert(Int, totalPop[i]))
            else println(f)
            end
        end
        if empty
            for g in misDist
                gd = gidData[g]
                print(f, g,"\t",gd[3],"\t",gd[1],"\t",gd[4],"\t",gd[2],"\t\t")
                if catList[j]=="Total" || catList[j]=="All"; print(f,"\t\t") end
                println(f)
            end
        end
        close(f)
    end

end

function printEmissionByExp(year, outputFile=""; percap=false, period="monthly", plot=false, dispmode=false, guimode=false)
                                                #xaxis: "daily", "weekly", "monthly", or "annual"
    global hhid, catList, inc, siz, emissionsHHs
    ec = emissionsHHs[year]

    nh= length(hhid)
    nc = length(catList)
    fce = zeros(Float64, nh)    # food carbon emission
    exp = zeros(Float64, nh)    # consumption expenditure
    mms = zeros(Int, nh)        # household members

    for i=1:nh
        fce[i] = ec[i,nc]
        exp[i] = inc[hhid[i]]
        mms[i] = siz[hhid[i]]
    end
    if percap; for i=1:nh; exp[i] /= sqrt(mms[i]) end end
    if period=="daily"; exp /= 30
    elseif period=="weekly"; exp *= 7/30
    elseif period=="annual"; exp *= 365/30
    elseif period!="monthly"; println("temporal axis is wrong")
    end

    order = sortperm(exp, rev=true)

    if length(outputFile)>0
        f = open(outputFile, "w")
        println(f,"HHID\tSize\tExpenditure_",period,"\tEmission_annual")
        for i=1:nh; println(f,hhid[order[i]],"\t",siz[hhid[order[i]]],"\t",exp[order[i]],"\t",fce[order[i]]) end
        close(f)
    end

    if plot
        plotly()
        nspan = 100
        span = [(maximum(exp)-minimum(exp))/(nspan-1)*(i-1)+minimum(exp) for i=1:nspan]
        reg = npr(exp, fce, xeval=span, reg=localconstant)
        #p = scatter(exp, fce, xaxis=("Expenditure("*xaxis*",USD)"), yaxis=("CF (tCO2/yr)"))
        #p = plot!(span, reg, width=2)
        p = plot(span, reg, width=2, xaxis=("Expenditure("*period*",USD)"), yaxis=("CF (tCO2/yr)"))
        if dispmode; display(p) end
        if guimode; gui() end
    end
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
