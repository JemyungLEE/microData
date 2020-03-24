module EmissionFoodCategorizer

# Developed date: 17. Feb. 2020
# Last modified date: 24. Mar. 2020
# Subject: Categorize carbon emissions in the food sectors
# Description: Categorize food emissions by districts (province, city, etc) and by religion,
#              expending the functions of the 'EmissionCategorizer' module.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using XLSX
using Plots
using Statistics
using KernelEstimator

hhid = Array{String, 1}()   # Household ID
sec = Array{String, 1}()    # India products or services sectors
secName = Dict{String, String}()    # sector name dictionary: {sector code, name}
relName = ["Hinduism","Islam","Christianity","Sikhism","Jainism","Buddhism","Zoroastrianism", "Others", "None"]

cat = Dict{String, String}()    # category dictionary: {sector code, sub-food-category}
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

catList = Array{String, 1}()    # category list: sub-category of food sections
disList = Array{String, 1}()    # district list
relList = Array{String, 1}()    # religion list
incList = Array{Float64, 1}()   # income sector list
levList = Array{Float64, 1}()   # carbon emission level sector list

emissionsHHs = Dict{Int16, Array{Float64, 2}}() # categozied emission by household: {year, {hhid, category}}
emissionsDis = Dict{Int16, Array{Float64, 2}}() # categozied emission by district: {year, {district, category}}
emissionsRel = Dict{Int16, Array{Float64, 2}}() # categozied emission by religion: {year, {religion, category}}
emissionsInc = Dict{Int16, Array{Float64, 2}}() # categozied emission by income: {year, {income level, category}}
emissionsIncRel = Dict{Int16, Array{Float64, 3}}()  # categozied emission by incomes: {year, {religion, income level, category}}
emissionsDisLev = Dict{Int16, Array{Float64, 2}}()  # categozied emission by district emission level: {year, {emission level, category}}
emissionsLev = Dict{Int16, Array{Float64, 2}}() # categozied emission by emission level: {year, {emission level, category}}

emissionsRelLev = Dict{Int16, Array{Float64, 3}}() # categozied emission by religion and by emission level: {year, {religion, emission level, category}}

function migrateData(year, ec)
    global sec, hhid = ec.sec, ec.hhid
    global emissions[year] = ec.emissions[year]
    global dis, sta, typ, siz, inc, rel, wgh = ec.dis, ec.sta, ec.typ, ec.siz, ec.inc, ec.rel, ec.wgh
    global gid, nam, pop, gidData, merDist, misDist = ec.gid, ec.nam, ec.pop, ec.gidData, ec.merDist, ec.misDist

    global disList = sort(unique(values(dis)))      # district list
    global relList = sort(unique(values(rel)))      # religion list
    if 0 in relList; deleteat!(relList, findall(x->x==0, relList)); push!(relList, 0) end
end

function readCategoryData(nat, inputFile)
    global cat, catList
    xf = XLSX.readxlsx(inputFile)
    sh = xf[nat*"_sec"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1 && !ismissing(r[5]); cat[string(r[1])] = r[5] end end
    close(xf)

    catList = sort(unique(values(cat)))      # category list
    push!(catList, "Food")
end

function categorizeHouseholdEmission(year)
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
end

function printHHsEmissionData(year, outputFile; sorting=false)
    global hhid, catList, emissionsHHs, siz, inc

    ec = emissionsHHs[year]
    nc = length(catList)
    nh = length(hhid)
    incomes = [inc[h] for h in hhid]

    f = open(outputFile, "w")

    if sorting; incOrder = sortperm(incomes, rev=true) end

    println(f, "HHID,Income,Emission,HH_size")
    for i=1:nh
        if sorting; idx = incOrder[i]
        else idx = i
        end
        print(f, hhid[idx],",",incomes[idx],",",ec[idx,nc],",",siz[hhid[idx]])
        println(f)
    end

    close(f)
end

function categorizeDistrictEmission(year; sqrRoot=false, period="monthly", religion=false)
    # weightMode: [0]non-weight, [1]population weighted, [2]household weighted, [3]both population and household weighted
    #             ([4],[5]: normalization) [4]per capita, [5]per household
    # sqrRoot: [true]apply square root of household size for an equivalance scale
    # period: "monthly", "daily", or "annual"
    # religion: [true] categorize districts' features by religions

    global hhid, cat, dis, siz, inc, sam, ave, rel
    global catList, disList, relList
    global emissionsHHs, emissionsDis

    nh = length(hhid)
    nc = length(catList)
    nd = length(disList)
    nr = length(relList)

    # make index arrays
    indDis = [findfirst(x->x==dis[hhid[i]], disList) for i=1:nh]
    if religion; indRel = [findfirst(x->x==rel[hhid[i]], relList) for i=1:nh] end

    # sum households and members by districts
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
        for i=1:nh; thbdr[indDis[i],indRel[i]] += 1 end
        if sqrRoot; for i=1:nh; tpbdr[indDis[i],indRel[i]] += sqrt(siz[hhid[i]]) end
        else for i=1:nh; tpbdr[indDis[i],indRel[i]] += siz[hhid[i]] end
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
    elseif sqrRoot; for i=1:nh; ed[indDis[i],:] += ec[i,:]/sqrt(siz[hhid[i]]) end
    end

    # normalizing
    if sqrRoot; for i=1:nc; ed[:,i] ./= thbd end
    else for i=1:nc; ed[:,i] ./= tpbd end
    end

    emissionsDis[year] = ed

    if religion; return ed, catList, disList, thbd, tpbd, thbdr, tpbdr
    else return ed, catList, disList, thbd, tpbd
    end
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
                                            # normmode: [1]per capita, [2]per houehold, [3]basic information
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
    # basic information
    elseif normMode == 3; ei[:,1], ei[:,2] = tpbi[:], thbi[:]
    end

    emissionsInc[year] = ei

    return ei, catList, incList, tpbi, thbi, twpbi
end

function categorizeEmissionLevel(year, intv=[], normMode = 0, squareRoot = false)
                                                    # intv: proportions between invervals of highest to lowest
                                                    # normmode: [1]per capita, [2]per houehold
    global hhid, catList, levList, cat, dis, siz
    global emissionsHHs, emissionsLev

    if length(intv) == 0; intv = [0.25,0.25,0.25,0.25]
    elseif sum(intv) != 1; sumi = sum(intv); intv /= sumi end

    ec = emissionsHHs[year]
    nc = length(catList)
    nl = length(intv)
    nh = length(hhid)

    # make index list
    levArray = zeros(Float64, nh)
    if squareRoot; for i=1:nh; levArray[i] = ec[i,nc]/sqrt(siz[hhid[i]]) end
    elseif normMode==1; for i=1:nh; levArray[i] = ec[i,nc]/siz[hhid[i]] end
    elseif normMode==2; for i=1:nh; levArray[i] = ec[i,nc] end
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
        tpbl[indLev[i]] += siz[hhid[i]]
    end

    # categorize emission data
    ec = emissionsHHs[year]
    el = zeros(Float64, nl, nc)
    if squareRoot; for i=1:nh; for j=1:nc; el[indLev[i],j] += ec[i,j]/sqrt(siz[hhid[i]]) end end
    else for i=1:nh; for j=1:nc; el[indLev[i],j] += ec[i,j] end end
    end

    # normalizing
    if squareRoot; for i=1:nc; el[:,i] ./= thbl end
    else for i=1:nc; el[:,i] ./= tpbl end
    end

    emissionsLev[year] = el
end

function categorizeHouseholdByIncomeByReligion(year,intv=[],normMode=0; sqrRt=false,absIntv=false,perCap=false,desOrd=false,popWgh=false)
                                            # intv: proportions between invervals of highest to lowest
                                            # absIntv: if "true", then intv[] is a list of income values, descending order
                                            # normmode: [1]per capita, [2]per houehold, [3]basic information
                                            # perCap: [true] per capita mode, [false] per household mode
                                            # desOrd: [true] descening order of 'intv[]', [false] ascending order
    global sec, hhid, cat, dis, siz, rel, inc, catList, relList, incList
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

function categorizeEmissionReligionLevel(year, intv=[], normMode = 0, squareRoot = false)
                                                    # intv: proportions between invervals of highest to lowest
                                                    # normmode: [1]per capita, [2]per household
    global hhid, catList, relList, levList, cat, dis, siz, rel
    global emissionsHHs, emissionsRelLev

    if length(intv) == 0; intv = [0.25,0.25,0.25,0.25]
    elseif sum(intv) != 1; sumi = sum(intv); intv /= sumi end

    ec = emissionsHHs[year]
    nc = length(catList)
    nr = length(relList)
    nl = length(intv)
    nh = length(hhid)

    # make religion index array
    indRel = zeros(Int, nh)     # index array of religion
    for i=1:nh; indRel[i] = findfirst(x->x==rel[hhid[i]], relList) end
    # make income index array
    levArray = zeros(Float64, nh)
    if squareRoot; for i=1:nh; levArray[i] = ec[i,nc]/sqrt(siz[hhid[i]]) end
    elseif normMode==1; for i=1:nh; levArray[i] = ec[i,nc]/siz[hhid[i]] end
    elseif normMode==2; for i=1:nh; levArray[i] = ec[i,nc] end
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
    if i == nh; indLev[levOrder[i]] = nl end

    # sum households and members by religion by income
    thbrbl = zeros(Int, nr, nl)   # total households by religion by income level
    tpbrbl = zeros(Int, nr, nl)   # total members of households by religion by income level
    for i=1:nh
        thbrbl[indRel[i],indLev[i]] += 1
        tpbrbl[indRel[i],indLev[i]] += siz[hhid[i]]
    end

    # categorize emission data
    erl = zeros(Float64, nr, nl, nc)
    if squareRoot; for i=1:nh; for j=1:nc; erl[indRel[i],indLev[i],j] += ec[i,j]/sqrt(siz[hhid[i]]) end end
    else for i=1:nh; for j=1:nc; erl[indRel[i],indLev[i],j] += ec[i,j] end end
    end

    # normalizing
    if squareRoot || normMode==2; for i=1:nc; erl[:,:,i] ./= thbrbl end
    elseif normMode==1 for i=1:nc; erl[:,:,i] ./= tpbrbl end
    end

    emissionsRelLev[year] = erl
end

function printHouseholdEmission(year, outputFile; hhsinfo=false)

    global catList, emissionsHHs, siz, inc
    eh = emissionsHHs[year]

    f = open(outputFile, "w")

    print(f,"HHID"); for c in catList; print(f, ",", c) end
    if hhsinfo; print(f, ",HH_size,MPCE") end; println(f)
    for i = 1:length(hhid)
        print(f, hhid[i]); for j = 1:length(catList); print(f, ",", eh[i,j]) end
        if hhsinfo; print(f, ",",siz[hhid[i]],",",inc[hhid[i]]); println(f) end
    end
    close(f)
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

function printEmissionLev(year, outputFile, intv=[])

    global catList, emissionsLev
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

function printEmissionRelLev(year, outputFile, intv=[])

    global catList, relList, emissionsRelLev
    erl = emissionsRelLev[year]
    relName = ["Hinduism","Islam","Christianity","Sikhism","Jainism","Buddhism","Zoroastrianism", "Others", "None"]

    f = open(outputFile, "w")

    for i=1:length(relList)
        println(f,"[",relName[i],"]")
        print(f,"CF_Lv")
        for c in catList; print(f, ",", c) end
        println(f)
        for j = 1:length(intv)
            print(f, "\t<", trunc(Int, sum(intv[1:j])*100),"%")
            for k = 1:length(catList); print(f, ",", erl[i,j,k]) end
            println(f)
        end
        println(f)
    end

    close(f)
end

function printEmissionByExp(year, outputFile=""; percap=true, period="monthly", plot=false, dispmode=false, guimode=false)
                                                #xaxis: "daily", "weekly", "monthly", or "annual"
    global hhid, catList, inc, siz, emissionsHHs
    ce = emissionsHHs[year]

    nh= length(hhid)
    nc = length(catList)
    fce = zeros(Float64, nh)    # food carbon emission
    exp = zeros(Float64, nh)    # consumption expenditure
    mms = zeros(Int, nh)        # household members

    for i=1:nh
        fce[i] = ce[i,nc]
        exp[i] = inc[hhid[i]]
        mms[i] = siz[hhid[i]]
    end

    if !percap; for i=1:nh; exp[i] *= mms[i] end end

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

end
