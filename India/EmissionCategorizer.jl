module EmissionCategorizer

# Developed date: 20. Dec. 2019
# Last modified date: 12. Feb. 2020
# Subject: Categorize India households carbon emissions
# Description: Categorize emissions by districts (province, city, etc) and by expenditure categories
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX
using Statistics

sec = Array{String, 1}()    # India products or services sectors
hhid = Array{String, 1}()   # Household ID

cat = Dict{String, String}()    # category dictionary: {sector code, category}
dis = Dict{String, String}()    # hhid's district: {hhid, district code}
typ = Dict{String, String}()    # hhid's sector type, urban or rural: {hhid, "urban" or "rural"}
siz = Dict{String, Int}()       # hhid's family size: {hhid, number of members}
inc = Dict{String, Float64}()   # hhid's income: {hhid, monthly per capita expenditure (mixed reference period)}
rel = Dict{String, Int}()       # hhid's religion: {hhid, religion code} [1]Hinduism,[2]Islam,[3]Christianity,[4]Sikhism,[5]Jainism,[6]Buddhism,[7]Zoroastrianism,[9]Others,[0]None

sam = Dict{String, Tuple{Int,Int}}()    # sample population and households by districct: {district code, (population, number of households)}
pop = Dict{String, Tuple{Int,Int}}()    # population by district: {district code, (population, number of households)}
ave = Dict{String, Float64}()    # average annual expenditure per capita, USD/yr: {district code, mean Avg.Exp./cap/yr}
nam = Dict{String, String}()    # districts' name: {district code, district name}
gid = Dict{String, String}()    # districts' gis_codes: {district code, gis id (GIS_2)}
gidData = Dict{String, Tuple{String, String, String, String}}() # GID code data: {GID_2, {district code, district name, state code, state name}}
merDist = Dict{String, String}()    # list of merged district: {merged district's code, remained district's code}

emissions = Dict{Int16, Array{Float64, 2}}()        # {year, table}

catList = Array{String, 1}()    # category list
disList = Array{String, 1}()    # district list
relList = Array{String, 1}()    # religion list
incList = Array{Float64, 1}()   # income sector list

emissionsCatNW = Dict{Int16, Array{Float64, 2}}()   # categozied emission by district, non-weighted: {year, {category, district}}

emissionsCat = Dict{Int16, Array{Float64, 2}}()     # categozied emission by district: {year, {category, district}}
emissionsRel = Dict{Int16, Array{Float64, 2}}()     # categozied emission by religion: {year, {category, religion}}
emissionsInc = Dict{Int16, Array{Float64, 2}}()     # categozied emission by incomes: {year, {category, income level}}
emissionsDis = Dict{Int16, Array{Float64, 2}}()     # categozied emission by district emission level: {year, {category, emission level}}
emissionsCatDif = Dict{Int16, Array{Float64, 2}}()  # differences of categozied emission by district: (emission-mean)/mean, {year, {category, district}}

gisEmissionCat = Dict{Int16, Array{Float64, 2}}()    # GIS version, categozied emission by district: {year, {category, district(GID)}}
gisEmissionCatDif = Dict{Int16, Array{Float64, 2}}() # GIS version, differences of categozied emission by district: (emission-mean)/mean, {year, {category, district(GID)}}

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

function readHouseholdData(year, inputFile, merging=false)

    global siz, dis, typ, inc, rel
    f = open(inputFile)

    readline(f)
    for l in eachline(f)
        l = split(l, '\t')
        siz[l[1]] = parse(Int,l[7])
        if merging==true&&haskey(merDist, l[5]); dis[l[1]] = merDist[l[5]]
        else dis[l[1]] = l[5]
        end
        typ[l[1]] = l[6]
        inc[l[1]] = parse(Float64,l[8])
        rel[l[1]] = parse(Int,l[9])
    end

    close(f)
end

function readCategoryData(nat, inputFile)

    global cat, gid, nam, pop, gidData, merDist
    xf = XLSX.readxlsx(inputFile)

    sh = xf[nat*"_sec"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1; cat[string(r[1])] = r[4] end end
    sh = xf[nat*"_dist"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1; gid[string(r[1])] = r[3]; nam[string(r[1])] = r[2] end end
    sh = xf[nat*"_pop"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1 && !ismissing(r[3]); pop[string(r[3])] = (r[9], r[8]) end end
    sh = xf[nat*"_gid"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1
        codes = split(r[3],r"[._]")
        gidData[string(r[3])]=(codes[3],r[4],codes[2],r[2])
    end end
    sh = xf[nat*"_mer"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1; merDist[string(r[3])] = string(r[1]) end end
    close(xf)
end

function categorizeEmission(year, weightMode=0, squareRoot=false)
    # weightMode: [0]non-weight, [1]population weighted, [2]household weighted, [3]both population and household weighted
    #             ([4],[5]: normalization) [4]per capita, [5]per household
    # squareRoot: [true]apply square root of household size for an equivalance scale

    global sec, hhid, cat, dis, siz, inc, sam, ave
    global emissions, emissionsCat, emissionsCatDif, emissionsCatNW

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

    # sum sample households and members by districts
    thbd = zeros(Float64, length(disList))   # total households by district
    tpbd = zeros(Float64, length(disList))   # total members of households by district
    for h in hhid
        thbd[indDis[h]] += 1
        if squareRoot; tpbd[indDis[h]] += sqrt(siz[h])
        else tpbd[indDis[h]] += siz[h]
        end
    end
    for i=1:nd; sam[disList[i]] = (tpbd[i], thbd[i]) end

    # calculate average annual expenditure per capita by district
    texp = zeros(Float64, length(disList))  # total expenditures by district
    temporalMultiplier = 365/30             # convert mothly values to annual ones
    for h in hhid; texp[indDis[h]] += inc[h]*siz[h] end
    for i=1:nd; ave[disList[i]] = texp[i]/tpbd[i]*temporalMultiplier end

    # categorize emission data
    e = emissions[year]
    ec = zeros(Float64, nc, nd)
    # categorizing
    for i=1:length(sec); for j=1:length(hhid); ec[indCat[sec[i]],indDis[hhid[j]]] += e[i,j] end end
    # summing
    for i = 1:nd; for j = 1:nc-1; ec[nc, i] += ec[j,i] end end

    # backup
    ecnw = deepcopy(ec)

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

    # calculate differences
    avg = mean(ec, dims=2)
    ecd = zeros(size(ec))
    for i=1:size(ec,1); ecd[i,:] = (ec[i,:].-avg[i])/avg[i] end

    # save the results
    emissionsCat[year] = ec
    emissionsCatDif[year] = ecd
    emissionsCatNW[year] = ecnw

    return ec, catList, disList, indCat, indDis, thbd, tpbd
end

function categorizeDistrict(year, indCat, normMode = 0, thbd=[], tpbd=[], intv=[])
                                            #intv: proportions between invervals of highest to lowest
                                            # normmode: [1]per capita, [2]per houehold, [3]basic information
    global sec, hhid, cat, dis, siz
    global catList, disList
    global emissionsCat

    ec = emissionsCat[year]

    if length(intv) == 0; intv = [0.25,0.25,0.25,0.25]
    elseif sum(intv) != 1; sumi = sum(intv); intv /= sumi end

    nh = length(hhid)
    nc = length(catList)
    nd = length(disList)
    ni = length(intv)

    disArray = []
    for i = 1:length(disList); push!(disArray, ec[end,i]) end
    disOrder = sortperm(disArray, rev=true)  #descending order indexing, [1]highest, [end]lowest values' indexes

    # make index dictionaries
    indDis = Dict{String, Int}()        # index dictionary of districts
    i = 1
    for s = 1:length(intv)
        while i <= trunc(Int, nd*sum(intv[1:s]))
            indDis[disList[disOrder[i]]] = s
            i += 1
        end
    end
    if i == nd; indDis[disList[disOrder[i]]] = length(intv) end

    # categorize emission data
    ed = zeros(Float64, nc, ni)
    # categorizing
    for i=1:nc; for j=1:nd; ed[i,indDis[disList[j]]] += ec[i,j] end end

    tp = zeros(Int, ni)
    th = zeros(Int, ni)
    for i=1:nd; tp[indDis[disList[i]]] += tpbd[i] end
    for i=1:nd; tp[indDis[disList[i]]] += thbd[i] end
    # normalizing
    if normMode == 1; for i=1:nc; for j=1:ni; ed[i,j] /= tp[j] end end
    elseif normMode == 2 ;for i=1:nc; for j=1:ni; ed[i,j] /= th[j] end end
    # basic information
    elseif normMode == 3; for i=1:ni; ed[1,i] = tp[i]; ed[2,i] = th[i] end
    end

    emissionsDis[year] = ed

    return ed, catList, disList
end

function categorizeReligion(year, normMode = 0, squareRoot = false, indCat=[])
                                            # normmode: [1]per capita, [2]per houehold, [3]basic information
    global sec, hhid, cat, dis, siz, rel
    global emissions, emissionsRel

    global catList = sort(unique(values(cat)))      # category list
    push!(catList, "Total")
    global relList = sort(unique(values(rel)))     # religion list
    if 0 in relList; deleteat!(relList, findall(x->x==0, relList)); push!(relList, 0) end

    nc = length(catList)
    nr = length(relList)

    # make index dictionaries
    if length(indCat)==0
        indCat = Dict{String, Int}()     # index dictionary of category
        for s in sec; indCat[s] = findfirst(x->x==cat[s], catList) end
    end
    indRel = Dict{String, Int}()     # index dictionary of religion
    for h in hhid; indRel[h] = findfirst(x->x==rel[h], relList) end

    # sum households and members by districts
    thbr = zeros(Float64, length(relList))   # total households by religion
    tpbr = zeros(Float64, length(relList))   # total members of households by religion
    for h in hhid
        thbr[indRel[h]] += 1
        if squareRoot; tpbr[indRel[h]] += sqrt(siz[h])
        else tpbr[indRel[h]] += siz[h]
        end
    end

    # categorize emission data
    e = emissions[year]
    er = zeros(Float64, nc, nr)
    # categorizing
    for i=1:length(sec); for j=1:length(hhid); er[indCat[sec[i]],indRel[hhid[j]]] += e[i,j] end end
    # summing
    for i = 1:nr; for j = 1:nc-1; er[nc, i] += er[j,i] end end

    # normalizing
    if normMode == 1
        for i=1:nc; for j=1:nr; er[i,j] /= tpbr[j] end end
    elseif normMode == 2
        for i=1:nc; for j=1:nr; er[i,j] /= thbr[j] end end
    # basic information
    elseif normMode == 3
        for i=1:nr; er[1,i] = tpbr[i]; er[2,i] = thbr[i] end
    end

    emissionsRel[year] = er

    return er, catList, relList
end

function categorizeIncome(year, intv=[], normMode = 0, squareRoot = false, indCat=[])
                                            #intv: proportions between invervals of highest to lowest
                                            # normmode: [1]per capita, [2]per houehold, [3]basic information
    global sec, hhid, cat, dis, siz, inc
    global catList, incList
    global emissions, emissionsInc

    catList = sort(unique(values(cat)))      # category list
    push!(catList, "Total")

    if length(intv) == 0; intv = [0.25,0.25,0.25,0.25]
    elseif sum(intv) != 1; sumi = sum(intv); intv /= sumi end

    nh = length(hhid)
    nc = length(catList)
    ni = length(intv)

    incArray = []
    for h in hhid; push!(incArray, inc[h]) end
    incOrder = sortperm(incArray, rev=true)  #descending order indexing, [1]highest, [end]lowest values' indexes
    for i=1:length(intv); push!(incList, incArray[incOrder[trunc(Int, sum(intv[1:i])*nh)]]) end

    # make index dictionaries
    if length(indCat)==0
        indCat = Dict{String, Int}()    # index dictionary of category
        for s in sec; indCat[s] = findfirst(x->x==cat[s], catList) end
    end
    indInc = Dict{String, Int}()        # index dictionary of income sections
    i = 1
    for s = 1:length(intv)
        while i <= trunc(Int, nh*sum(intv[1:s]))
            indInc[hhid[incOrder[i]]] = s
            i += 1
        end
    end
    if i == nh; indInc[hhid[incOrder[i]]] = length(intv) end

    # sum households and members by districts
    thbi = zeros(Float64, length(incList))   # total households by religion
    tpbi = zeros(Float64, length(incList))   # total members of households by religion
    for h in hhid
        thbi[indInc[h]] += 1
        if squareRoot; tpbi[indInc[h]] += sqrt(siz[h])
        else tpbi[indInc[h]] += siz[h]
        end
    end

    # categorize emission data
    e = emissions[year]
    ei = zeros(Float64, nc, ni)
    # categorizing
    for i=1:length(sec); for j=1:length(hhid); ei[indCat[sec[i]],indInc[hhid[j]]] += e[i,j] end end
    # summing
    for i = 1:ni; for j = 1:nc-1; ei[nc, i] += ei[j,i] end end

    # normalizing
    if normMode == 1
        for i=1:nc; for j=1:ni; ei[i,j] /= tpbi[j] end end
    elseif normMode == 2
        for i=1:nc; for j=1:ni; ei[i,j] /= thbi[j] end end
    # basic information
    elseif normMode == 3
        for i=1:ni; ei[1,i] = tpbi[i]; ei[2,i] = thbi[i] end
    end

    emissionsInc[year] = ei

    return ei, catList, incList
end

function printCategorizedEmission(year, outputFile, name = false)

    global nam, catList, disList, emissionsCat
    ec = emissionsCat[year]

    f = open(outputFile, "w")

    if name
        for c in catList; print(f, "\t", c) end
        println(f)
        for i = 1:length(disList)
            print(f, nam[disList[i]])
            for j = 1:length(catList); print(f, "\t", ec[j,i]) end
            println(f)
        end
    else
        for d in disList; print(f, "\t", d) end
        println(f)
        for i = 1:length(catList)
            print(f, catList[i])
            for j = 1:length(disList); print(f, "\t", ec[i,j]) end
            println(f)
        end
    end

    close(f)
end

function printCategorizedReligion(year, outputFile)

    global catList, relList, emissionsRel
    er = emissionsRel[year]

    relName = ["Hinduism","Islam","Christianity","Sikhism","Jainism","Buddhism","Zoroastrianism", "Others", "None"]

    f = open(outputFile, "w")

    for r in relName; print(f, "\t", r) end
    println(f)
    for i = 1:length(catList)
        print(f, catList[i])
        for j = 1:length(relList); print(f, "\t", er[i,j]) end
        println(f)
    end

    close(f)
end

function printCategorizedIncome(year, outputFile, intv=[])

    global catList, incList, emissionsInc
    ei = emissionsInc[year]

    f = open(outputFile, "w")

    for i in intv; print(f, "\t<", trunc(Int, i*100),"%") end
    println(f)
    for i = 1:length(catList)
        print(f, catList[i])
        for j = 1:length(incList); print(f, "\t", ei[i,j]) end
        println(f)
    end

    close(f)
end

function printCategorizedDistrict(year, outputFile, intv=[])

    global catList, disList, emissionsDis
    ed = emissionsDis[year]

    f = open(outputFile, "w")

    for i in intv; print(f, "\t<", trunc(Int, i*100),"%") end
    println(f)
    for i = 1:length(catList)
        print(f, catList[i])
        for j = 1:length(intv); print(f, "\t", ed[i,j]) end
        println(f)
    end

    close(f)
end

function exportEmissionTable(year, tag, outputFile, weightMode::Int, name=false)

    global sam, pop, ave, gid, gidData, catList, disList, emissionsCatNW
    ecnw = emissionsCatNW[year]

    # making exporting table
    gidList = sort(unique(values(gid)))
    tb = zeros(Float64, length(gidList), length(catList))
    spo = zeros(Float64, length(gidList))   # number of sample population by district
    shh = zeros(Float64, length(gidList))   # number of sample households by district
    tpo = zeros(Float64, length(gidList))   # total number of population by district
    thh = zeros(Float64, length(gidList))   # total number of households by district
    aec = zeros(Float64, length(gidList))   # average expenditure per capita by district
    for i=1:length(disList)
        idx = findfirst(x->x==gid[disList[i]],gidList)
        for j=1:length(catList); tb[idx, j] += ecnw[j, i] end
        spo[idx] += sam[disList[i]][1]
        shh[idx] += sam[disList[i]][2]
        tpo[idx] += pop[disList[i]][1]
        thh[idx] += pop[disList[i]][2]
        aec[idx] += ave[disList[i]]*sam[disList[i]][1]
    end
    for i=1:length(gidList); aec[i] /= spo[i] end
    if weightMode==1; for i=1:length(gidList); for j=1:length(catList); tb[i,j] *= tpo[i]/spo[i] end end
    elseif weightMode==2; for i=1:length(gidList); for j=1:length(catList); tb[i,j] *= thh[i]/shh[i] end end
    elseif weightMode==4; for i=1:length(gidList); for j=1:length(catList); tb[i,j] /= spo[i] end end
    elseif weightMode==5; for i=1:length(gidList); for j=1:length(catList); tb[i,j] /= shh[i] end end
    end

    # exporting table
    f = open(outputFile, "w")
    print(f, tag)
    for c in catList; print(f,",",c) end
    println(f)
    for i = 1:size(tb, 1)
        if name; print(f, gidData[gidList[i]][2])
        else print(f, gidList[i])
        end
        for j = 1:size(tb, 2); print(f, ",", tb[i,j]) end
        println(f)
    end
    close(f)

    gisEmissionCat[year] = tb

    return tb, gidList, spo, shh, tpo, thh, aec
end

function exportEmissionDiffRate(year, tag, outputFile, name=false)

    global gid, gidData
    gec = gisEmissionCat[year]
    gidList = sort(unique(values(gid)))

    # calculate difference rates
    avg = mean(gec, dims=1)
    gecd = zeros(size(gec))
    for i=1:size(gec,2); gecd[:,i] = (gec[:,i].-avg[i])/avg[i] end

    # exporting table
    f = open(outputFile, "w")
    print(f, tag)
    for c in catList; print(f,",",c) end
    println(f)
    for i = 1:size(gecd, 1)
        if name; print(f, gidData[gidList[i]][2])
        else print(f, gidList[i])
        end
        for j = 1:size(gecd, 2); print(f, ",", gecd[i,j]) end
        println(f)
    end
    close(f)

    gisEmissionCatDif[year] = gecd
end

function exportWebsiteFiles(year, path, weightMode, gidList, totalPop, totalHH, sampPop, avgExp)

    global nam, gid, gidData, catList
    global gisEmissionCat, gisEmissionCatDif
    gec = gisEmissionCat[year]
    gecd = gisEmissionCatDif[year]

    for j=1:length(catList)
        f = open(path*"CFAC_"*catList[j]*"_"*string(year)*".txt","w")
        println(f, "KEY_CODE\tSTATE\tDISTRICT\tSTATE_NAME\tDISTRICT_NAME\tGENHH_ALLPPC")
        for i=1:length(gidList)
            gd = gidData[gidList[i]]
            println(f, gidList[i],"\t",gd[3],"\t",gd[1],"\t",gd[4],"\t",gd[2],"\t",gecd[i,j])
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
            if weightMode == 1; print(f,gec[i,j],"\t",gec[i,j]/totalPop[i])
            elseif weightMode == 2; print(f,gec[i,j],"\t",gec[i,j]/totalHH[i])
            elseif weightMode == 4; print(f,gec[i,j]*totalPop[i],"\t",gec[i,j])
            elseif weightMode == 5; print(f,gec[i,j]*totalHH[i],"\t",gec[i,j])
            elseif weightMode == 0; print(f,gec[i,j]*totalHH[i]/sampPop[i],"\t",gec[i,j]/sampPop[i])
            end
            if catList[j]=="Total" || catList[j]=="All"; println(f,"\t",avgExp[i],"\t",convert(Int, totalPop[i]))
            else println(f)
            end
        end
        close(f)
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
