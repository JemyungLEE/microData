module EmissionFoodCategorizer

# Developed date: 17. Feb. 2020
# Last modified date: 20. Feb. 2020
# Subject: Categorize carbon emissions in the food sectors
# Description: Categorize food emissions by districts (province, city, etc) and by religion,
#              expending the functions of the 'EmissionCategorizer' module.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using XLSX
using Plots
using Statistics
using KernelEstimator

sec = Array{String, 1}()    # India products or services sectors
hhid = Array{String, 1}()   # Household ID

cat = Dict{String, String}()    # category dictionary: {sector code, sub-food-category}
dis = Dict{String, String}()    # hhid's district: {hhid, district code}
typ = Dict{String, String}()    # hhid's sector type, urban or rural: {hhid, "urban" or "rural"}
siz = Dict{String, Int}()       # hhid's family size: {hhid, number of members}
inc = Dict{String, Float64}()   # hhid's income: {hhid, monthly per capita expenditure (mixed reference period)}
rel = Dict{String, Int}()       # hhid's religion: {hhid, religion code} [1]Hinduism,[2]Islam,[3]Christianity,[4]Sikhism,[5]Jainism,[6]Buddhism,[7]Zoroastrianism,[9]Others,[0]None

sam = Dict{String, Tuple{Int,Int}}()    # sample population and households by districct: {district code, (population, number of households)}
pop = Dict{String, Tuple{Int,Int}}()    # population by district: {district code, (population, number of households)}
ave = Dict{String, Float64}()   # average annual expenditure per capita, USD/yr: {district code, mean Avg.Exp./cap/yr}
nam = Dict{String, String}()    # districts' name: {district code, district name}
gid = Dict{String, String}()    # districts' gis_codes: {district code, gis id (GIS_2)}
gidData = Dict{String, Tuple{String, String, String, String}}() # GID code data: {GID_2, {district code, district name, state code, state name}}
merDist = Dict{String, String}()    # list of merged district: {merged district's code, remained district's code}

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
emissionsLev = Dict{Int16, Array{Float64, 2}}() # categozied emission by emission level: {year, {emission level, category}}

emissionsRelInc = Dict{Int16, Array{Float64, 3}}() # categozied emission by religion and by income: {year, {religion, income level, category}}
emissionsRelLev = Dict{Int16, Array{Float64, 3}}() # categozied emission by religion and by emission level: {year, {religion, emission level, category}}

function migrateData(year, ec)
    global sec, hhid = ec.sec, ec.hhid
    global emissions[year] = ec.emissions[year]
    global siz, dis, typ, inc, rel = ec.siz, ec.dis, ec.typ, ec.inc, ec.rel
    global gid, nam, pop, gidData, merDist = ec.gid, ec.nam, ec.pop, ec.gidData, ec.merDist

    global disList = sort(unique(values(dis)))      # district list
    global relList = sort(unique(values(rel)))     # religion list
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

function categorizeEmissionHouseholds(year)
    global sec, hhid, catList, cat
    global emissions, emissionsHHs

    nc = length(catList)
    nh = length(hhid)
    ns = length(sec)

    indCat = Dict{String, Int}()     # index dictionary of category

    # make index dictionaries
    for s in sec; if haskey(cat, s); indCat[s] = findfirst(x->x==cat[s], catList) end end

    # categorize emission data
    e = emissions[year]
    ec = zeros(Float64, nh, nc)
    # categorizing
    for i=1:nh; for j=1:ns; if haskey(indCat, sec[j]); ec[i,indCat[sec[j]]] += e[j,i] end end end
    # summing
    for i=1:nh; for j=1:nc-1; ec[i, nc] += ec[i,j] end end

    # save the results
    emissionsHHs[year] = ec
end

function categorizeEmissionReligion(year, squareRoot = false)
    global hhid, catList, relList, cat, dis, siz, rel
    global emissionsHHs, emissionsRel

    nc = length(catList)
    nr = length(relList)
    nh = length(hhid)

    # make index array
    indRel = zeros(Int, nh)     # index array of religion
    for i=1:nh; indRel[i] = findfirst(x->x==rel[hhid[i]], relList) end

    # sum households and members by districts
    thbr = zeros(Float64, nr)   # total households by religion
    tpbr = zeros(Float64, nr)   # total members of households by religion
    for i=1:nh
        thbr[indRel[i]] += 1
        tpbr[indRel[i]] += siz[hhid[i]]
    end

    # categorize emission data
    ec = emissionsHHs[year]
    er = zeros(Float64, nr, nc)
    if squareRoot; for i=1:nh; for j=1:nc; er[indRel[i],j] += ec[i,j]/sqrt(siz[hhid[i]]) end end
    else for i=1:nh; for j=1:nc; er[indRel[i],j] += ec[i,j] end end
    end

    # normalizing
    if squareRoot; for i=1:nc; er[:,i] ./= thbr end
    else for i=1:nc; er[:,i] ./= tpbr end
    end

    emissionsRel[year] = er
end

function categorizeEmissionDistrict(year, squareRoot = false)
    global hhid, catList, disList, cat, dis, siz
    global emissionsHHs, emissionsDis

    nc = length(catList)
    nd = length(disList)
    nh = length(hhid)

    # make index array
    indDis = zeros(Int, nh)     # index array of religion
    for i=1:nh; indDis[i] = findfirst(x->x==dis[hhid[i]], disList) end

    # sum households and members by districts
    thbd = zeros(Float64, nd)   # total households by religion
    tpbd = zeros(Float64, nd)   # total members of households by religion
    for i=1:nh
        thbd[indDis[i]] += 1
        tpbd[indDis[i]] += siz[hhid[i]]
    end
    for i=1:nd; sam[disList[i]] = (tpbd[i], thbd[i]) end

    # categorize emission data
    ec = emissionsHHs[year]
    ed = zeros(Float64, nd, nc)
    if squareRoot; for i=1:nh; for j=1:nc; ed[indDis[i],j] += ec[i,j]/sqrt(siz[hhid[i]]) end end
    else for i=1:nh; for j=1:nc; ed[indDis[i],j] += ec[i,j] end end
    end

    # normalizing
    if squareRoot; for i=1:nc; ed[:,i] ./= thbd end
    else for i=1:nc; ed[:,i] ./= tpbd end
    end

    emissionsDis[year] = ed
end

function categorizeEmissionIncome(year, intv=[], normMode=0, squareRoot=false; absintv=false)
                                    # intv: proportions between invervals of highest to lowest
                                    # absintv: if "true", then intv[] is a list of income values, descending order
                                    # normmode: [1]per capita, [2]per household
    global hhid, catList, incList, cat, dis, siz, inc
    global emissionsHHs, emissionsInc

    if !absintv && length(intv) == 0; intv = [0.25,0.25,0.25,0.25]
    elseif !absintv && sum(intv) != 1; sumi = sum(intv); intv /= sumi end

    nc = length(catList)
    ni = length(intv)
    nh = length(hhid)

    # make index array
    incArray = zeros(Float64, nh)
    if squareRoot; for i=1:nh; incArray[i] = inc[hhid[i]]/sqrt(siz[hhid[i]]) end
    elseif normMode==1; for i=1:nh; incArray[i] = inc[hhid[i]]/siz[hhid[i]] end
    elseif normMode==2; for i=1:nh; incArray[i] = inc[hhid[i]] end
    end
    incOrder = sortperm(incArray, rev=true)  #descending order indexing, [1]highest, [end]lowest values' indexes
    if absintv; incList[:] = intv[:]
    else for i=1:ni; push!(incList, incArray[incOrder[trunc(Int, sum(intv[1:i])*nh)]]) end
    end

    indInc = zeros(Int, nh)     # index array of income sections
    if absintv  # indInc: '1' for over intv[1], '3' for below intv[2], [2] for others
        for i=1:nh
            if incArray[incOrder[i]] >= intv[1]; indInc[incOrder[i]] = 1
            elseif incArray[incOrder[i]] < intv[1]; indInc[incOrder[i]] = 3
            else intv[1]; indInc[incOrder[i]] = 2
            end
        end
    else
        i = 1
        for s = 1:length(intv)
            while i <= trunc(Int, nh*sum(intv[1:s]))
                indInc[incOrder[i]] = s
                i += 1
            end
        end
        if i == nh; indInc[incOrder[i]] = length(intv) end
    end

    # sum households and members by districts
    thbi = zeros(Int, ni)   # total households by income level
    tpbi = zeros(Int, ni)   # total members of households by income level
    for i=1:nh
        thbi[indInc[i]] += 1
        tpbi[indInc[i]] += siz[hhid[i]]
    end

    # categorize emission data
    ec = emissionsHHs[year]
    ei = zeros(Float64, ni, nc)
    if squareRoot; for i=1:nh; for j=1:nc; ei[indInc[i],j] += ec[i,j]/sqrt(siz[hhid[i]]) end end
    else for i=1:nh; for j=1:nc; ei[indInc[i],j] += ec[i,j] end end
    end

    # normalizing
    if squareRoot || normMode==2; for i=1:nc; ei[:,i] ./= thbi end
    elseif normMode==1 for i=1:nc; ei[:,i] ./= tpbi end
    end

    emissionsInc[year] = ei
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
    if i == nh; indLev[levOrder[i]] = nl end

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

function categorizeEmissionReligionIncome(year, intv=[], normMode = 0, squareRoot = false)
                                                    # intv: proportions between invervals of highest to lowest
                                                    # normmode: [1]per capita, [2]per household
    global hhid, catList, relList, incList, cat, dis, siz, rel, inc
    global emissionsHHs, emissionsRelInc

    if length(intv) == 0; intv = [0.25,0.25,0.25,0.25]
    elseif sum(intv) != 1; sumi = sum(intv); intv /= sumi end

    nc = length(catList)
    nr = length(relList)
    ni = length(intv)
    nh = length(hhid)

    # make religion index array
    indRel = zeros(Int, nh)     # index array of religion
    for i=1:nh; indRel[i] = findfirst(x->x==rel[hhid[i]], relList) end
    # make income index array
    incArray = zeros(Float64, nh)
    if squareRoot; for i=1:nh; incArray[i] = inc[hhid[i]]/sqrt(siz[hhid[i]]) end
    elseif normMode==1; for i=1:nh; incArray[i] = inc[hhid[i]]/siz[hhid[i]] end
    elseif normMode==2; for i=1:nh; incArray[i] = inc[hhid[i]] end
    end
    incOrder = sortperm(incArray, rev=true)  #descending order indexing, [1]highest, [end]lowest values' indexes
    for i=1:ni; push!(incList, incArray[incOrder[trunc(Int, sum(intv[1:i])*nh)]]) end

    indInc = zeros(Int, nh)     # index array of income sections
    i = 1
    for s = 1:ni
        while i <= trunc(Int, nh*sum(intv[1:s]))
            indInc[incOrder[i]] = s
            i += 1
        end
    end
    if i == nh; indInc[incOrder[i]] = ni end

    # sum households and members by religion by income
    thbrbi = zeros(Int, nr, ni)   # total households by religion by income level
    tpbrbi = zeros(Int, nr, ni)   # total members of households by religion by income level
    for i=1:nh
        thbrbi[indRel[i],indInc[i]] += 1
        tpbrbi[indRel[i],indInc[i]] += siz[hhid[i]]
    end

    # categorize emission data
    ec = emissionsHHs[year]
    eri = zeros(Float64, nr, ni, nc)
    if squareRoot; for i=1:nh; for j=1:nc; eri[indRel[i],indInc[i],j] += ec[i,j]/sqrt(siz[hhid[i]]) end end
    else for i=1:nh; for j=1:nc; eri[indRel[i],indInc[i],j] += ec[i,j] end end
    end

    # normalizing
    if squareRoot || normMode==2; for i=1:nc; eri[:,:,i] ./= thbrbi end
    elseif normMode==1 for i=1:nc; eri[:,:,i] ./= tpbrbi end
    end

    emissionsRelInc[year] = eri
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

function printEmissionHHs(year, outputFile)

    global catList, emissionsHHs
    eh = emissionsHHs[year]

    f = open(outputFile, "w")

    print(f,"HHID")
    for c in catList; print(f, ",", c) end
    println(f)
    for i = 1:length(hhid)
        print(f, hhid[i])
        for j = 1:length(catList); print(f, ",", eh[i,j]) end
        println(f)
    end

    close(f)
end

function printEmissionDis(year, outputFile, name=false)

    global nam, catList, disList, emissionsDis
    ed = emissionsDis[year]

    f = open(outputFile, "w")

    print(f,"District")
    for c in catList; print(f, ",", c) end
    println(f)
    for i = 1:length(disList)
        if name; print(f, nam[disList[i]]) else print(f, disList[i]) end
        for j = 1:length(catList); print(f, ",", ed[i,j]) end
        println(f)
    end

    close(f)
end

function printEmissionRel(year, outputFile)

    global catList, relList, emissionsRel
    er = emissionsRel[year]
    relName = ["Hinduism","Islam","Christianity","Sikhism","Jainism","Buddhism","Zoroastrianism", "Others", "None"]

    f = open(outputFile, "w")

    print(f,"Religion")
    for c in catList; print(f, ",", c) end
    println(f)
    for i = 1:length(relList)
        print(f, relName[i])
        for j = 1:length(catList); print(f, ",", er[i,j]) end
        println(f)
    end

    close(f)
end

function printEmissionInc(year, outputFile, intv=[]; absintv=false)

    global catList, emissionsInc
    ei = emissionsInc[year]

    f = open(outputFile, "w")

    print(f,"Exp_Lv")
    for c in catList; print(f, ",", c) end
    println(f)
    for i = 1:length(intv)
        if absintv
            if i==1; print(f, "< ",intv[1])
            elseif i==2; print(f, "< ",intv[2])
            elseif i==3; print(f, "> ",intv[2])
            end
        else print(f, "\t<", trunc(Int, sum(intv[1:i])*100),"%")
        end
        for j = 1:length(catList); print(f, ",", ei[i,j]) end
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

function printEmissionRelInc(year, outputFile, intv=[]; plot=false, gui=false)

    global catList, relList, emissionsRelInc
    eri = emissionsRelInc[year]
    relName = ["Hinduism","Islam","Christianity","Sikhism","Jainism","Buddhism","Zoroastrianism", "Others", "None"]

    f = open(outputFile, "w")
    for i=1:length(relList)
        println(f,"[",relName[i],"]")
        print(f,"Exp_Lv")
        for c in catList; print(f, ",", c) end
        println(f)
        for j = 1:length(intv)
            print(f, "\t<", trunc(Int, sum(intv[1:j])*100),"%")
            for k = 1:length(catList); print(f, ",", eri[i,j,k]) end
            println(f)
        end
        println(f)
    end
    close(f)

    if plot
        plotly()

    end
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

function printEmissionByExp(year, outputFile=""; percap=false, period="monthly", plot=false, dispmode=false, guimode=false)
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

end
