module EmissionCategorizer

# Developed date: 3. Aug. 2020
# Last modified date: 2. Dec. 2020
# Subject: Categorize EU households' carbon footprints
# Description: Read household-level CFs and them by consumption category, district, expenditure-level, and etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using XLSX
using Statistics

hhsList = Dict{Int, Dict{String, Array{String, 1}}}()   # Household ID: {year, {nation, {hhid}}}
sec = Array{String, 1}()            # Consumption products' or services' sectors
secName = Dict{String, String}()    # sector name dictionary: {sector code, name}
ceSec = Array{String, 1}()          # direct emission causing consumption sectors
ceCodes = Array{Array{String, 1}, 1}()  # direct emission related consumption sectors: {CE category, {expenditure sectors}}

# hhid -> nation(2digit)_hhid: some HHIDs are duplicated across multiple countries
cat = Dict{String, String}()                # category dictionary: {sector code, category}
nat = Dict{Int, Dict{String, String}}()     # hhid's nation: {year, {hhid, nation code}}
reg = Dict{Int, Dict{String, String}}()     # hhid's NUTS: {year, {hhid, NUTS code}}
typ = Dict{Int, Dict{String, String}}()     # hhid's sector type, urban or rural: {year, {hhid, "urban" or "rural"}}
siz = Dict{Int, Dict{String, Int}}()        # hhid's family size: {year, {hhid, number of members}}
eqs = Dict{Int, Dict{String, Float64}}()    # hhid's family equivalent size (OECD scale): {year, {hhid, number of members}}
meqs = Dict{Int, Dict{String, Float64}}()   # hhid's family equivalent size (modified OECD scale): {year, {hhid, number of members}}
inc = Dict{Int, Dict{String, Float64}}()    # hhid's income: {year, {hhid, total income}}
exp = Dict{Int, Dict{String, Float64}}()    # hhid's domestic expenditure: {year, {hhid, total domestic expenditure}}
pds = Dict{Int, Dict{String, Int}}()        # hhid region's population density: {year, {hhid, district's population density}}
rel = Dict{Int, Dict{String, Int}}()        # hhid's religion: {year, {hhid, religion code}}
wgh = Dict{Int, Dict{String, Float64}}()    # hhid's weight: {year, {hhid, weight}}
wghNuts = Dict{Int, Dict{String, Float64}}()    # hhid's NUTS weight: {year, {hhid, weight}}

nutsLv = 0      # NUTS level
nuts = Dict{Int, Dict{String, String}}()       # NUTS: {year, {code, label}}
nutsList = Dict{Int, Dict{String, Array{String, 1}}}()      # NUTS code list: {year, {nation_code, {NUTS_code}}}
pop = Dict{Int, Dict{String, Float64}}()        # Population: {year, {NUTS_code, population}}
popList = Dict{Int, Dict{String, Dict{String, Float64}}}()  # Population list: {year, {nation_code, {NUTS_code, population}}}

emissions = Dict{Int, Dict{String, Array{Float64, 2}}}()    # carbon footprint: {year, {nation, {table}}}
directCE = Dict{Int, Dict{String, Array{Float64, 2}}}()     # direct carbon emission: {year, {nation, {table}}}

yrList = Array{Int, 1}()        # year list
catList = Array{String, 1}()    # category list
natList = Array{String, 1}()    # nation list
natName = Dict{String,String}() # nation's code and full-name: {nation code, full-name}

emissionsHHs = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied emission by household: {year, {nation, {hhid, category}}}
emissionsReg = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied emission by region: {year, {nation, {region, category}}}
emissionsRegDiff = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied emission differences by region: {year, {nation, {region, category}}}

gisNutsList = Dict{Int, Array{String, 1}}()             # GIS version, NUTS list: {year, {region(giscd)}}
gisRegionalEmission = Dict{Int, Array{Float64, 2}}()    # GIS version, categozied emission by district: {year, {category, region(giscd)}}
gisRegionalEmissionRank = Dict{Int, Array{Float64, 2}}()# GIS version, categozied emission rank by district: {year, {category, region(giscd)}}
gisRegionalEmissionDiff = Dict{Int, Array{Float64, 2}}()# GIS version, differences of categozied emission by district: (emission-mean)/mean, {year, {category, district(GID)}}
gisRegionalEmissionDiffRank = Dict{Int, Array{Int, 2}}()# GIS version, difference ranks of categozied emission by district: (emission-mean)/mean, {year, {category, district(GID)}}

gisTotPop = Dict{Int, Array{Float64, 1}}()      # GIS version, total population by NUTS
gisSamPop = Dict{Int, Array{Float64, 1}}()      # GIS version, total sample members by NUTS
gisAvgExp = Dict{Int, Array{Float64, 1}}()      # GIS version, average expenditure by NUTS

regList = Dict{Int, Array{String, 1}}() # district list

sam = Dict{Int, Dict{String, Tuple{Int,Int}}}() # sample population and households by districct: {district code, (population, number of households)}
ave = Dict{Int, Dict{String, Float64}}()        # average annual expenditure per capita, USD/yr: {district code, mean Avg.Exp./cap/yr}

hbscd = Dict{Int, Dict{String, String}}()   # concordance NUTS code: {year, {HBS NUTS code, replaced HBS NUTS code}}
popcd = Dict{Int, Dict{String, String}}()   # concordance NUTS code: {year, {HBS NUTS code, Population NUTS code}}
giscd = Dict{Int, Dict{String, String}}()   # concordance NUTS code: {year, {GIS NUTS code, corresponding NUTS code}}
popcdlist = Dict{Int, Array{String, 1}}()   # Population NUTS code list: {year, Population NUTS code}
giscdlist = Dict{Int, Array{String, 1}}()   # GIS NUTS code list: {year, GIS NUTS code}

###################
typList = Array{String, 1}()    # area type list
relList = Array{String, 1}()    # religion list
incList = Array{Float64, 1}()   # income sector list
levList = Array{Float64, 1}()   # carbon emission level sector list

emissionsRel = Dict{Int, Array{Float64, 2}}()     # categozied emission by religion: {year, {religion, category}}
emissionsInc = Dict{Int, Array{Float64, 2}}()     # categozied emission by incomes: {year, {income level, category}}
emissionsIncRel = Dict{Int, Array{Float64, 3}}()  # categozied emission by incomes by religion: {year, {religion, income level, category}}
emissionsDisLev = Dict{Int, Array{Float64, 2}}()  # categozied emission by district emission level: {year, {emission level, category}}
emissionsRng = Dict{Int16, Array{Float64, 2}}()     # categozied emission by expenditure range: {year, {range, category}}
emissionsLev = Dict{Int16, Array{Float64, 2}}()     # categozied emission by emission level: {year, {emission level, category}}

emissionsDisDiff = Dict{Int16, Array{Float64, 2}}() # differences of categozied emission by district: (emission-mean)/mean, {year, {district, category}}

emissionCostDis = Dict{Int16, Array{Float64, 2}}()     # categozied emission by district: {year, {district, category}}
gisDistrictEmissionCost = Dict{Int16, Array{Float64, 2}}()    # GIS version, total emission cost for poverty alleivation: {year, {category, district(GID)}}

function makeNationalSummary(year, outputFile)

    global hhsList, natList, natName, siz, wgh
    global emissions, directCE

    nn = length(natList)
    natsam = zeros(Int, nn)
    nateqs = zeros(Float64, nn)
    natmeqs = zeros(Float64, nn)
    natwgh = zeros(Float64, nn)
    natcf = zeros(Float64, nn)      # Overall CF
    natcfpc = zeros(Float64, nn)    # CF per capita
    natcfph = zeros(Float64, nn)    # CF per household
    natcfpeqs = zeros(Float64, nn)  # CF per equivalent size
    natcfpmeqs = zeros(Float64, nn) # CF per modified equivalent size
    natce = zeros(Float64, nn)      # Overall CE
    natcepc = zeros(Float64, nn)    # CE per capita

    for i=1:nn
        n = natList[i]
        for j = 1:length(hhsList[year][n])
            h = hhsList[year][n][j]
            cf = sum(emissions[year][n][:,j])
            natsam[i] += siz[year][h]
            nateqs[i] += eqs[year][h]
            natmeqs[i] += meqs[year][h]
            natwgh[i] += wgh[year][h] * siz[[year]h]
            natcf[i] += wgh[year][h] * cf
            natcfpc[i] += cf
            natcfph[i] += cf
            natcfpeqs[i] += cf
            natcfpmeqs[i] += cf
            ce = sum(directCE[year][n][:,j])
            natce[i] += wgh[year][h] * ce
            natcepc[i] += ce
        end
        natcfpc[i] /= natsam[i]
        natcfph[i] /= length(hhsList[year][n])
        natcfpeqs[i] /= nateqs[i]
        natcfpmeqs[i] /= natmeqs[i]
        natcepc[i] /= natsam[i]
    end

    f = open(outputFile, "w")
    println(f, "Nation\tHHs\tMMs\tWeights\tCF_overall\tCF_percapita\tCF_perhh\tCF_pereqs\tCF_permeqs\tCE_overall\tCE_percapita")
    for i=1:nn
        print(f, natList[i],"\t",length(hhsList[year][natList[i]]),"\t",natsam[i],"\t",natwgh[i])
        print(f, "\t",natcf[i],"\t",natcfpc[i],"\t",natcfph[i],"\t",natcfpeqs[i],"\t",natcfpmeqs[i])
        print(f, "\t", natce[i], "\t", natcepc[i])
        println(f)
    end
    close(f)
end

function readCategoryData(inputFile, year=[], ntlv=0; subCategory="", except=[])

    global sec, ceSec, cat, gid, nam, pop, gidData, misDist, ceCodes
    global natList, natName, nuts, nutslList, pop, popList
    global hbscd, popcd, giscd, popcdlist, giscdlist
    global nutsLv = ntlv
    xf = XLSX.readxlsx(inputFile)

    sh = xf["Sector"]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r)>1 && !ismissing(r[1])
            secCode = string(r[1])  # sector code
            push!(sec, secCode)
            if length(subCategory)==0 && !ismissing(r[4]) && !(string(r[4]) in except); cat[secCode] = string(r[4])
            elseif subCategory=="Food" && !ismissing(r[5]); cat[secCode] = string(r[5])
            end
            secName[secCode]=string(r[2])
        end
    end
    sh = xf["CE_sector"]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r)>1 && !ismissing(r[1])
            push!(ceSec, string(r[1]))
            push!(ceCodes, [string(r[i]) for i=2:4 if !ismissing(r[i])])
        end
    end
    sh = xf["Nation"]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r)>1 && !ismissing(r[1])
            push!(natList, string(r[1]))
            natName[natList[end]] = string(r[2])
        end
    end
    for y in year
        nuts[y] = Dict{String, String}()
        nutsList[y] = Dict{String, Array{String, 1}}()
        pop[y] = Dict{String, Float64}()
        popList[y] = Dict{String, Dict{String, Float64}}()
        for n in natList
            nutsList[y][n] = Array{String, 1}()
            popList[y][n] = Dict{String, Float64}()
        end
        popcd[y] = Dict{String, String}()
        hbscd[y] = Dict{String, String}()
        giscd[y] = Dict{String, String}()
        popcdlist[y] = Array{String, 1}()
        giscdlist[y] = Array{String, 1}()
    end
    for y in year
        sh = xf["NUTS"*string(y)]
        for r in XLSX.eachrow(sh)
            if XLSX.row_number(r)>1
                lv = parse(Int, string(r[3]))
                ntcd = string(r[1])
                n = string(r[4])
                if length(ntcd)==lv+2; nuts[y][ntcd] = string(r[2]) else println("NUTS level error: ", ntcd, "\t", lv) end
                if n in natList && lv==ntlv && ntcd[end] != 'Z'; push!(nutsList[y][n], ntcd) end
            end
        end

        # add aggregated regions for region-code absent nations
        for n in ["DE","FR","EL"]
            rg = n
            for i=1:ntlv; rg *= "0" end
            push!(nutsList[y][n], rg)
        end

        sh = xf["Conc"*string(y)]
        for r in XLSX.eachrow(sh)
            if XLSX.row_number(r)>1
                hbscd[y][string(r[1])] = string(r[2])
                giscd[y][string(r[1])] = string(r[3])
                if !(string(r[3]) in giscdlist[y]); push!(giscdlist[y], string(r[3])) end
            end
        end
        sh = xf["PopCd"*string(y)]
        for r in XLSX.eachrow(sh)
            if XLSX.row_number(r)>1 && !ismissing(r[2])
                popcd[y][string(r[1])] = string(r[2])
                if !(string(r[1]) in popcdlist[y]); push!(popcdlist[y], string(r[1])) end
            end
        end
    end

    sh = xf["Pop_mod"]
    yridx = []
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r)==1
            str = [string(r[i]) for i=1:13]
            yridx = [findfirst(x->x==string(y), str) for y in year]
        else
            ntcd = string(split(r[1], ',')[end])
            n = ntcd[1:2]
            if n in natList
                for i=1:length(year)
                    if haskey(popcd[year[i]], ntcd)
                        ncd = popcd[year[i]][ntcd]
                        if ncd in nutsList[year[i]][n]
                            s = strip(replace(string(r[yridx[i]]), ['b','d','e','p']=>""))
                            if tryparse(Float64, s)!==nothing
                                pop[year[i]][ntcd] = parse(Float64, s)
                                if !haskey(popList[year[i]][n], ncd); popList[year[i]][n][ncd] = 0 end
                                popList[year[i]][n][ncd] += pop[year[i]][ntcd]
                            end
                        end
                    end
                end
            end
        end
    end
    close(xf)

    # for y in year, n in natList
    #     for c in sort(collect(keys(popList[y][n])))
    #         for p in popList[y][n][c]
    #             println(y,"\t",n,"\t",c,"\t",p)
    #         end
    #     end
    # end
    # for y in year; println(y,"\t",giscdlist[y]) end
    # for y in year, n in natList; println(y,"\t",n,"\t",nutsList[y][n]) end
    # for y in year, n in natList
    #     for cd in sort(collect(keys(popList[y][n])))
    #         println(y,"\t",n,"\t",cd)
    #     end
    # end

    global catList = sort(unique(values(cat)))
    if length(subCategory)==0; push!(catList, "Total")      # category list
    else push!(catList, subCategory)
    end
end

function setCategory(list::Array{String,1})
    global catList
    if sort(catList)==sort(list); catList = list
    else println("Category items are different.\n",catList,"\n",sort(list))
    end
end

function readHouseholdData(inputFile; period="monthly", sampleCheck=false)  # period: "monthly"(default), "daily", or "annual"
    global yrList, hhsList, natList, regList, relList, disSta, nutsLv
    global nat, reg, siz, eqs, meqs, typ, inc, rel

    year = 0
    nation = ""
    f = open(inputFile)
    readline(f)
    for l in eachline(f)
        s = string.(split(l, ','))
        year = parse(Int, s[1])
        if !(year in yrList)
            push!(yrList, year)
            hhsList[year] = Dict{String, Array{String, 1}}()
            nat[year] = Dict{String, String}()
            reg[year] = Dict{String, String}()
            typ[year] = Dict{String, String}()
            siz[year] = Dict{String, Int}()
            eqs[year] = Dict{String, Float64}()
            meqs[year] = Dict{String, Float64}()
            inc[year] = Dict{String, Float64}()
            exp[year] = Dict{String, Float64}()
            pds[year] = Dict{String, Int}()
            rel[year] = Dict{String, Int}()
            wgh[year] = Dict{String, Float64}()
        end
        if !haskey(hhsList[year], s[2]); hhsList[year][s[2]] = Array{String, 1}() end
        hh = s[2]*"_"*s[3]      # replaced from "hh = s[3]"
        push!(hhsList[year][s[2]], hh)
        siz[year][hh] = parse(Int,s[5])
        eqs[year][hh] = parse(Float64,s[12])
        meqs[year][hh] = parse(Float64,s[13])
        nat[year][hh] = s[2]
        reg[year][hh] = giscd[year][s[4]][1:(nutsLv+2)]
        wgh[year][hh] = parse(Float64,s[6])
        inc[year][hh] = parse(Float64,s[7])
        exp[year][hh] = parse(Float64,s[9])
        pds[year][hh] = parse(Int,s[11])
    end

    # convert household's income and expenditure data period
    if period=="daily"
        for n in collect(keys(hhsList[year]))
            for h in hhsList[year][n]; inc[year][h] = inc[year][h]/30; exp[year][h] = exp[year][h]/30 end
        end
    elseif period=="annual"
        mmtoyy=365/30
        for n in collect(keys(hhsList[year]))
            for h in hhsList[year][n]; inc[year][h] = inc[year][h]*mmtoyy; exp[year][h] = exp[year][h]*mmtoyy end
        end
    end
    close(f)

    regList[year] = sort(unique(values(reg[year])))

    # check sample numbers by district's population density
    if sampleCheck
        samples = zeros(Int, length(regList[year]), 4)    # {total, high_dens, middle_dens, low_dens}
        for n in collect(keys(hhsList[year]))
            for h in hhsList[year][n]
                idx = findfirst(x->x==reg[year][h], regList[year])
                samples[idx,1] += siz[year][h]
                samples[idx,pds[h]+1] += siz[year][h]
            end
        end
        f = open(Base.source_dir() * "/data/extracted/SampleNumber.txt","w")
        println(f,"District\tAll\tHigh_density\tMiddle_density\tLow_density")
        for i=1:length(regList[year]); print(f,regList[year][i]); for sn in samples[i,:]; println(f,"\t",sn) end end
        close(f)
    end
end

function readCarbonFootprint(year, nations, inputFiles)

    global sec, hhsList
    global emissions[year] = Dict{String, Array{Float64, 2}}()

    ns = length(sec)
    nn = length(nations)
    if length(inputFiles) == nn
        for i = 1:nn
            n = nations[i]
            nh = length(hhsList[year][n])
            f = open(inputFiles[i])
            readline(f)
            e = zeros(Float64, ns, nh)
            for l in eachline(f)
                l = split(l, '\t')
                j = findfirst(x->x==string(l[1]), sec)
                l = l[2:end]
                e[j,:] = map(x->parse(Float64,x),l)
            end
            emissions[year][n] = e
            close(f)
        end
    else println("Sizes of nation list (", nn,") and emission files (", length(inputFiles),") do not match")
    end
end

function readDirectEmission(year, nations, inputFiles)

    global ceSec, hhsList
    global directCE[year] = Dict{String, Array{Float64, 2}}()

    ns = length(ceSec)
    nn = length(nations)
    if length(inputFiles) == nn
        for i = 1:nn
            n = nations[i]
            nh = length(hhsList[year][n])
            f = open(inputFiles[i])
            readline(f)
            ce = zeros(Float64, ns, nh)
            for l in eachline(f)
                l = split(l, '\t')
                j = findfirst(x->x==string(l[1]), ceSec)
                l = l[2:end]
                ce[j,:] = map(x->parse(Float64,x),l)
            end
            directCE[year][n] = ce
            close(f)
        end
    else println("Sizes of nation list (", nn,") and emission files (", length(inputFiles),") do not match")
    end
end

function readExpenditure(year, nations, inputFiles)

    # global sec, hhsList, emissions
    #
    # ns = length(sec)
    # nn = length(nations)
    # if length(inputFiles) == nn
    #     for i = 1:nn
    #         n = nations[i]
    #         nh = length(hhsList[year][n])
    #         f = open(inputFiles[i])
    #         readline(f)
    #         e = zeros(Float64, ns, nh)      # expenditure matrix: {sector, hhid}
    #         j = 1
    #         for l in eachline(f)
    #             l = split(l, '\t')[2:end-1]
    #             if j<=nh; e[:,j] = map(x->parse(Float64,x),l) end
    #             j += 1
    #         end
    #
    #         global emissions[year][n] = e
    #         close(f)
    #     end
    # else  println("Sizes of nation list (", nn,") and emission files (", length(inputFiles),") do not match")
    # end
end

function calculateNutsPopulationWeight()

    global yrList, natList, hhsList, nutsList, popList, wghNuts
    global siz, reg, pds

    ntpop = Dict{Int,Dict{String,Dict{String,Array{Int,1}}}}()      # NUTS population:{year,{nation,{NUTS,population{total,dense,mid,sparse,none}}}}
    ntsmp = Dict{Int,Dict{String,Dict{String,Array{Int,1}}}}()      # NUTS sample size:{year,{nation,{NUTS,sample number{total,dense,mid,sparse,none}}}}
    ntwgh = Dict{Int,Dict{String,Dict{String,Array{Float64,1}}}}()  # NUTS population weight:{year,{nation,{NUTS,population{total,dense,mid,sparse,none}}}}

    # count region samples
    typeidx = [2,3,4,0,0,0,0,0,5]
    for y in yrList
        ntsmp[y] = Dict{String, Dict{String, Array{Int,1}}}()
        for n in natList
            ntsmp[y][n] = Dict{String, Array{Int,1}}()
            for nt in nutsList[y][n]; ntsmp[y][n][nt] = zeros(Int, 5) end
            for hh in hhsList[y][n]
                ntsmp[y][n][reg[y][hh]][typeidx[pds[y][hh]]] += siz[y][hh]
                ntsmp[y][n][reg[y][hh]][1] += siz[y][hh]
            end
        end
    end

    # calculate weights
    for y in yrList
        ntwgh[y] = Dict{String, Dict{String, Array{Int,1}}}()
        for n in natList
            ntwgh[y][n] = Dict{String, Array{Int,1}}()
            for nt in nutsList[y][n]
                if nt in giscdlist[y]
                    ntwgh[y][n][nt] = zeros(Float64, 5)
                    ntwgh[y][n][nt][1] = popList[y][n][nt] / ntsmp[y][n][nt][1]
                end
            end
        end
    end

    # allocate household weight
    for y in yrList
        wghNuts[y] = Dict{String, Float64}()
        for n in natList, hh in hhsList[y][n]; wghNuts[y][hh] = ntwgh[y][n][reg[y][hh]][1] end
    end
end

function categorizeHouseholdEmission(years; output="", hhsinfo=false, nutsLv=1)
    global wgh, sec, hhid, cat, siz, inc, catList, natList
    global emissions, emissionsHHs

    nc = length(catList)
    ns = length(sec)

    # make an index dict
    catidx = [if haskey(cat, s); findfirst(x->x==cat[s], catList) end for s in sec]

    # categorize emission data
    if isa(years, Number); years = [years] end
    for y in years
        emissionsHHs[y] = Dict{String, Array{Float64, 2}}()
        for n in natList
            nh = length(hhsList[y][n])
            e = emissions[y][n]
            ec = zeros(Float64, nh, nc)
            # categorizing
            for i=1:nh; for j=1:ns; if catidx[j]!=nothing; ec[i,catidx[j]] += e[j,i] end end end
            # summing
            for i=1:nc-1; ec[:, nc] += ec[:,i] end
            # save the results
            emissionsHHs[y][n] = ec
        end
    end

    # print the results
    if length(output)>0
        f = open(output, "w")
        print(f,"Year,Nation,HHID"); for c in catList; print(f, ",", c) end
        if hhsinfo; print(f, ",HH_size,MPCE,PopWgh") end; println(f)
        for y in years
            for n in natList
                for i = 1:length(hhsList[y][n])
                    hh = hhsList[y][n][i]
                    print(f, y,',',n,',',hh)
                    for j = 1:length(catList); print(f, ",", emissionsHHs[y][n][i,j]) end
                    if hhsinfo; print(f, ",",siz[y][hh],",",inc[y][hh],",",wgh[y][hh]) end
                    println(f)
                end
            end
        end
        close(f)
    end
end

function categorizeRegionalEmission(years=[], weightMode=0; nutsLv=1, period="monthly", religion=false, popWgh=false, ntweigh=false)
    # weightMode: [0]non-weight, [1]per capita, [2]per household
    # period: "monthly", "daily", or "annual"
    # religion: [true] categorize districts' features by religions

    global hhsList, natList, cat, reg, siz, inc, sam, ave, rel, pop, wgh, wghNuts
    global catList, nutsList, relList, popcd
    global emissionsHHs, emissionsReg, emissionsRegDiff

    if ntweigh; pwgh = wghNuts else pwgh = wgh end

    nc = length(catList)
    nr = length(relList)

    for y in years
        emissionsReg[y] = Dict{String, Array{Float64, 2}}()
        emissionsRegDiff[y] = Dict{String, Array{Float64, 2}}()
        sam[y] = Dict{String, Tuple{Int,Int}}()
        ave[y] = Dict{String, Float64}()

        for n in natList
            hhs = hhsList[y][n]
            nts = nutsList[y][n]
            nh = length(hhs)
            nn = length(nts)

            # make index arrays
            ntidx = [findfirst(x->x==reg[y][hhs[i]], nts) for i=1:nh]
            if religion; relidx = [findfirst(x->x==rel[y][hhs[i]], relList) for i=1:nh] end

            # sum sample households and members by regions
            thbd = zeros(Float64, nn)   # total households by region
            tpbd = zeros(Float64, nn)   # total members of households by region
            for i=1:nh; thbd[ntidx[i]] += 1 end
            # if sqrRoot; for i=1:nh; tpbd[ntidx[i]] += sqrt(siz[hhs[i]]) end
            # else
                for i=1:nh; tpbd[ntidx[i]] += siz[y][hhs[i]] end
            # end
            for i=1:nn; sam[y][nts[i]] = (tpbd[i], thbd[i]) end

            # sum sample households and members by regions and by religions
            if religion
                thbdr = zeros(Float64, nn, nr)  # total households by district, by religion
                tpbdr = zeros(Float64, nn, nr)  # total members of households by district, by religion
                for i=1:nh
                    thbdr[ntidx[i],relidx[i]] += 1
                    # if sqrRoot; tpbdr[ntidx[i],relidx[i]] += sqrt(siz[hhs[i]])
                    # else
                        tpbdr[ntidx[i],relidx[i]] += siz[y][hhs[i]]
                    # end
                end
            end

            # calculate average monthly expenditure per capita by region
            totexp = zeros(Float64, nn)     # total expenditures by region
            if popWgh
                for i=1:nh; totexp[ntidx[i]] += inc[y][hhs[i]]*siz[y][hhs[i]]*pwgh[y][hhs[i]] end
                for i=1:nn; if nts[i] in giscdlist[y]; ave[y][nts[i]] = totexp[i]/popList[y][n][nts[i]] end end
            else
                for i=1:nh; totexp[ntidx[i]] += inc[y][hhs[i]]*siz[y][hhs[i]] end
                for i=1:nn; ave[y][nts[i]] = totexp[i]/tpbd[i] end
            end
            # convert 'AVEpC' to annual or daily
            if period=="annual"; mmtoyy = 365/30; for i=1:nn; ave[y][nts[i]] = ave[y][ntd[i]] * mmtoyy end
            elseif period=="daily"; for i=1:nn; if nts[i] in giscdlist[y]; ave[y][nts[i]] = ave[y][nts[i]]/30 end end
            end

            # categorize emission data
            ec = emissionsHHs[y][n]
            en = zeros(Float64, nn, nc)
            # if !sqrRoot
                if popWgh; for i=1:nh; en[ntidx[i],:] += ec[i,:]*pwgh[y][hhs[i]] end
                else for i=1:nh; en[ntidx[i],:] += ec[i,:] end
                end
            # elseif sqrRoot && weightMode==2; for i=1:nh; en[ntidx[i],:] += ec[i,:]/sqrt(siz[hhs[i]]) end
            # end

            # normalizing
            if weightMode == 1
                if popWgh; for i=1:nn; if nts[i] in giscdlist[y]; for j=1:nc; en[i,j] /= popList[y][n][nts[i]] end end end
                else for i=1:nc; en[:,i] ./= tpbd end
                end
            elseif weightMode == 2; for i=1:nc; en[:,i] ./= thbd end
            end

            # calculate differences
            avg = mean(en, dims=1)
            ecn = zeros(size(en))
            for i=1:size(en,2); ecn[:,i] = (en[:,i].-avg[i])/avg[i] end

            # save the results
            emissionsReg[y][n] = en
            emissionsRegDiff[y][n] = ecn
        end
    end
end

function printRegionalEmission(years, outputFile; totm=false,expm=false,popm=false,relm=false,wghm=false,denm=false,povm=false,ntweigh=false)
    # expm: print average expenditure per capita
    # popm: print population related figures
    # hhsm: print households related figures
    # relm: print relgion related figures
    # denm: print population density
    # povm: print poverty rates

    global natList, hhsList, catList, nutsList, relList, relName, popList, ave, emissionsReg
    global wgh, wghNuts

    if ntweigh; pwgh = wghNuts else pwgh = wgh end
    f = open(outputFile, "w")

    nc = length(catList)
    print(f,"Year,Nation,NUTS")
    for c in catList; print(f, ",", c) end
    if totm; print(f, ",Overall_CF") end
    if expm; print(f, ",Income") end
    if popm; print(f, ",Population") end
    if wghm; print(f, ",PopWeight") end
    # if povm; print(f, ",PovertyRatio") end
    println(f)
    for y in years, n in natList
        nts = nutsList[y][n]
        en = emissionsReg[y][n]
        for i in [x for x = 1:length(nts) if !isnan(en[x,end])&&en[x,end]!=0]
            print(f, y,',',n,',',nts[i])
            for j = 1:nc; print(f, ",", en[i,j]) end
            if haskey(popList[y][n], nts[i])
                if totm; print(f, ",", en[i,end]*popList[y][n][nts[i]]) end
                if expm; if nts[i] in giscdlist[y]; print(f, ",", ave[y][nts[i]]) else print(f, ",NA") end end
                if popm; print(f, ",", popList[y][n][nts[i]]) end
                if wghm; print(f, ",", sum([pwgh[y][h]*siz[y][h] for h in filter(x->reg[y][x]==nts[i],hhsList[y][n])])) end
                # if povm; print(f, ",", disPov[regList[i]]) end
            end
            println(f)
        end
    end

    close(f)
end

function exportRegionalEmission(years=[], tag="", outputFile=""; percap=false, nspan=128, minmax=[], descend=false, empty=false, logarithm=false)

    global catList, nuts, nutsList, gisNutsList, giscd, natList, popList, sam, ave, reg
    global emissionsReg, gisRegionalEmission, gisRegionalEmissionRank

    nc = length(catList)

    for y in years
        nts = Dict{String, Array{String, 1}}()
        ntslist = Array{String, 1}()
        ntkeys = sort(collect(keys(nuts[y])))
        for n in natList
            nts[n] = filter(x->x in ntkeys, nutsList[y][n])
            append!(ntslist, nts[n])
        end
        nn = length(ntslist)

        # making exporting table
        tb = zeros(Float64, nn, nc)
        spo = zeros(Float64, nn)   # number of sample population by region
        tpo = zeros(Float64, nn)   # total number of population by region
        aec = zeros(Float64, nn)   # average expenditure per capita by region

        for n in natList
            ec = emissionsReg[y][n]
            for nt in nts[n]
                gnt = giscd[y][nt]
                gidx = findfirst(x->x==nt, ntslist)
                ntidx = findfirst(x->x==gnt, nutsList[y][n])
                tb[gidx,:] += ec[ntidx,:] * popList[y][n][gnt]
                spo[gidx] += sam[y][gnt][1]
                tpo[gidx] += popList[y][n][gnt]
                aec[gidx] += ave[y][gnt]*popList[y][n][gnt]
            end
        end
        # normalizing
        if percap
            for i=1:nn
                aec[i] /= tpo[i]
                for j=1:nc; tb[i,j] /= tpo[i] end
            end
        end

        # find min. and max.
        if length(minmax)==1; maxde = [minmax[1][2] for i=1:nc]; minde = [minmax[1][1] for i=1:nc]
        elseif length(minmax)==nc; maxde = [minmax[i][2] for i=1:nc]; minde = [minmax[i][1] for i=1:nc]
        elseif logarithm; maxde = [log10(maximum(tb[:,i])) for i=1:nc]; minde = [log10(minimum(tb[:,i])) for i=1:nc]
        else maxde = [maximum(tb[:,i]) for i=1:nc]; minde = [minimum(tb[:,i]) for i=1:nc]
        end
        replace!(minde, Inf=>0, -Inf=>0)

        # grouping by ratios; ascending order
        span = zeros(Float64, nspan+1, nc)
        for j=1:nc; span[:,j] = [(maxde[j]-minde[j])*(i-1)/nspan+minde[j] for i=1:nspan+1] end
        if logarithm; for i=1:size(span,1); for j=1:nc; span[i,j] = 10^span[i,j] end end end

        # println(maxde)
        # println(minde)
        # for i=1:size(span)[2]
        #     println(span[:,i])
        # end

        # grouping by rank; ascending order
        rank = zeros(Int, nn, nc)
        for j=1:nc
            for i=1:nn
                if tb[i,j]>=span[end-1,j]; rank[i,j] = nspan
                elseif tb[i,j] <= span[1,j]; rank[i,j] = 1
                else rank[i,j] = findfirst(x->x>=tb[i,j],span[:,j]) - 1
                end
            end
        end
        # for descending order, if "descend == true"
        if descend
            for i=1:nc; span[:,i] = reverse(span[:,i]) end
            for j=1:nc; for i=1:ngid; rank[i,j] = nspan - rank[i,j] + 1 end end
        end
        # prepare labels
        labels = Array{String,2}(undef,nspan,nc)
        for j=1:nc; labels[:,j] = [string(round(span[i,j],digits=0))*"-"*string(round(span[i+1,j],digits=0)) for i=1:nspan] end

        # exporting table
        outputFile = replace(outputFile,"YEAR_"=>string(y)*"_")
        f = open(outputFile, "w")
        print(f, tag); for c in catList; print(f,",",c) end; println(f)
        for i = 1:size(tb, 1)
            print(f, ntslist[i])
            for j = 1:size(tb, 2); print(f, ",", tb[i,j]) end
            println(f)
        end
        if empty
            # for not-covered NUTS data
        end
        close(f)

        # exporting group table
        f = open(replace(outputFile,".csv"=>"_gr.csv"), "w")
        print(f, tag); for c in catList; print(f,",",c) end; println(f)
        for i = 1:nn
            print(f, ntslist[i])
            for j = 1:nc; print(f, ",", rank[i,j]) end
            println(f)
        end
        if empty
            # for not-covered NUTS data
        end
        close(f)

        gisTotPop[y] = tpo
        gisSamPop[y] = spo
        gisAvgExp[y] = aec
        gisNutsList[y] = ntslist
        gisRegionalEmission[y] = tb
        gisRegionalEmissionRank[y] = rank
    end
end

function exportEmissionDiffRate(years=[], tag="", outputFile="", maxr=0.5, minr=-0.5, nspan=128; descend=false, empty=false)

    global catList, nutsList, gisNutsList, gisRegionalEmission, gisRegionalEmissionDiff, gisRegionalEmissionDiffRank

    nc = length(catList)

    for y in years
        ntslist = gisNutsList[y]
        gre = gisRegionalEmission[y]

        # calculate difference rates
        avg = mean(gre, dims=1)
        gred = zeros(size(gre))
        for i=1:size(gre,2); gred[:,i] = (gre[:,i].-avg[i])/avg[i] end

        # grouping by ratios; ascending order
        span = [(maxr-minr)*(i-1)/(nspan-2)+minr for i=1:nspan-1]
        spanval = zeros(Float64, nspan, nc)
        for i=1:nc
            spanval[1:end-1,i] = span[:].*avg[i].+avg[i]
            spanval[end,i] = spanval[end-1,i]
        end

        rank = zeros(Int, size(gred))
        for j=1:size(gred,2)    # category number
            for i=1:size(gred,1)    # gid district number
                if gred[i,j]>=maxr; rank[i,j] = nspan
                else rank[i,j] = findfirst(x->x>gred[i,j],span)
                end
            end
        end
        # for descending order, if "descend == true".
        if descend
            for j=1:size(gred,2); for i=1:size(gred,1); rank[i,j] = nspan - rank[i,j] + 1 end end
        end

        # exporting difference table
        outputFile = replace(outputFile,"YEAR_"=>string(y)*"_")
        f = open(outputFile, "w")
        print(f, tag); for c in catList; print(f,",",c) end; println(f)
        for i = 1:size(gred, 1)
            print(f, ntslist[i])
            for j = 1:size(gred, 2); print(f, ",", gred[i,j]) end
            println(f)
        end
        if empty
            # for not-covered NUTS data
        end
        close(f)

        # exporting difference group table
        f = open(replace(outputFile,".csv"=>"_gr.csv"), "w")
        print(f, tag); for c in catList; print(f,",",c) end; println(f)
        for i = 1:size(rank, 1)
            print(f, ntslist[i])
            for j = 1:size(rank, 2); print(f, ",", rank[i,j]) end
            println(f)
        end
        if empty
            # for not-covered NUTS data
        end
        close(f)

        gisRegionalEmissionDiff[y] = gred
        gisRegionalEmissionDiffRank[y] = rank
    end
end

function exportWebsiteFiles(years, path; percap=false, rank=false, empty=false)

    global natName, nuts, catList, gisNutsList, gisTotPop, gisAvgExp
    global gisRegionalEmission, gisRegionalEmissionDiff, gisRegionalEmissionDiffRank
    global gisEmissionCat, gisEmissionCatDif

    for y in years
        ntslist = gisNutsList[y]
        gre = gisRegionalEmission[y]
        gred = gisRegionalEmissionDiff[y]
        gredr = gisRegionalEmissionDiffRank[y]

        for j=1:length(catList)
            if percap
                f = open(path*string(y)*"/"*"CFAC_"*catList[j]*"_"*string(y)*".txt","w")
                println(f, "KEY_CODE\tCOUNTRY\tNUTS1\tCOUNTRY_NAME\tNUTS_NAME\tGENHH_ALLPPC")
                for i=1:length(ntslist)
                    nt = ntslist[i]
                    print(f, nt,"\t",nt[1:2],"\t",nt,"\t",natName[nt[1:2]],"\t",nuts[y][nt],"\t")
                    if rank; print(f, gredr[i,j]) else print(f, gred[i,j]) end
                    println(f)
                end
                if empty
                    # for not-covered NUTS data
                end
                close(f)
            else
                f = open(path*string(y)*"/"*"CFAV_"*catList[j]*"_"*string(y)*".txt","w")
                print(f, "KEY_CODE\tCOUNTRY\tNUTS1\tCOUNTRY_NAME\tNUTS_NAME\tGENHH_ALLP\tGENHH_APPPC")
                if catList[j]=="Total" || catList[j]=="All"; println(f, "\tANEXPPC\tPOP")
                else println(f)
                end
                for i=1:length(ntslist)
                    nt = ntslist[i]
                    print(f, nt,"\t",nt[1:2],"\t",nt,"\t",natName[nt[1:2]],"\t",nuts[y][nt],"\t")
                    print(f, gre[i,j],"\t",gre[i,j]/gisTotPop[y][i])
                    if catList[j]=="Total" || catList[j]=="All"; println(f,"\t",gisAvgExp[y][i],"\t",convert(Int, gisTotPop[y][i]))
                    else println(f)
                    end
                end
                if empty
                    # for not-covered NUTS data
                end
                close(f)
            end
        end
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
        print(f, "Category"); for i=1:nts; print(f, ",Sector_no.",i) end; println(f)
        for i=1:nc
            print(f, catlist[i])
            for j=1:length(orderSec[i]); print(f, ",",secName[orderSec[i][j]]," (",round(propSec[i][j],digits=3),")") end
            println(f)
        end
        close(f)
    end
end

function calculateDistrictPoverty(year; povline=1.9, popWgh=false)

    global hhid, dis, siz, inc, pop, sam, disPov
    global regList
    nd = length(regList)

    povr = zeros(Float64, nd)
    for h in hhid
        if inc[h]<povline
            idx = findfirst(x->x==dis[h], regList)
            if popWgh; povr[idx] += siz[h]*wghDis[h]
            else povr[idx] += siz[h]
            end
        end
    end
    if popWgh; for i=1:nd; povr[i] /= pop[regList[i]][1] end
    else for i=1:nd; povr[i] /= sam[regList[i]][1] end
    end

    for i=1:nd; disPov[regList[i]] = povr[i] end
end

function categorizeDistrictByEmissionLevel(year, normMode = 0, intv=[])
                                            # intv: proportions between invervals of highest to lowest
                                            # normmode: [1]per capita, [2]per household, [3]basic information
    global hhid, sam, catList, regList
    global emissionsDis
    ed = emissionsDis[year]

    if length(intv) == 0; intv = [0.25,0.25,0.25,0.25]
    elseif sum(intv) != 1; sumi = sum(intv); intv /= sumi
    end

    nh = length(hhid)
    nc = length(catList)
    nd = length(regList)
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
        tp[indDis[i]] += sam[regList[i]][1]
        th[indDis[i]] += sam[regList[i]][2]
    end
    # normalizing
    if normMode == 1; for i=1:nc; edl[:,i] ./= tp end
    elseif normMode == 2 ;for i=1:nc; edl[:,i] ./= th end
    # basic information
    elseif normMode == 3; ed[:,1], ed[:,2] = tp[:], th[:]
    end

    emissionsDisLev[year] = edl

    return edl, catList, regList
end

function categorizeHouseholdByReligion(year, normMode=0; sqrRt=false, popWgh=false, wghmode="district")
                                            # normmode: [1]per capita, [2]per houehold, [3]basic information
    global wghSta, wghDis; if wghmode=="state"; wgh=wghSta elseif wghmode=="district"; wgh=wghDis end
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

function categorizeHouseholdByIncome(year,intv=[],normMode=0; sqrRt=false,absIntv=false,perCap=false,desOrd=false,popWgh=false,wghmode="district")
                                            # intv: proportions between invervals of highest to lowest
                                            # absIntv: if "true", then intv[] is a list of income values, descending order
                                            # normmode: [1]per capita, [2]per houehold
                                            # perCap: [true] per capita mode, [false] per household mode
                                            # desOrd: [true] descening order of 'intv[]', [false] ascending order
    global wghSta, wghDis; if wghmode=="state"; wgh=wghSta elseif wghmode=="district"; wgh=wghDis end
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
    indInc = zeros(Int, nh)
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
            if accpop/totpop > intv[idx] && i < nh
                push!(pcidx, i)
                push!(incList, incArray[incOrder[i]])
                idx += 1
            end
            indInc[incOrder[i]] = idx
        end
    elseif absIntv; incList = copy(intv)
    else
        push!(incList, incArray[incOrder[1]])
        for i=1:length(intv); push!(incList, incArray[incOrder[trunc(Int, intv[i]*nh)]]) end
    end

    # set sector index
    if absIntv  # indInc: '1' for over intv[1], '3' for below intv[2], [2] for others
        for i=1:nh
            if incArray[i] >= intv[1]; indInc[i] = 1
            elseif incArray[i] < intv[2]; indInc[i] = 3
            else indInc[i] = 2
            end
        end
    elseif !perCap
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

    return ei, catList, incList, tpbi, thbi, twpbi, indInc
end

function categorizeHouseholdByExpRange(year,rng=[],normMode=0; perCap=false,popWgh=false,absRng=false,absSpn=false,over=0.1,less=0.1,wghmode="district")
                                            # absRng: [true] apply absolute range, [false] apply population ratio range
                                            # rng: standard values of ranges for grouping
                                            # over/less: range from 'stdv', ratios of samples, househods, or population
                                            # normMode: [1]per capita, [2]per houehold
                                            # perCap: [true] per capita mode, [false] per household mode
    global wghSta, wghDis; if wghmode=="state"; wgh=wghSta elseif wghmode=="district"; wgh=wghDis end
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
        if absSpn
            for i=1:nr
                rsidx[i,1] = findmin(abs.([expArray[x]-rng[i] for x in expOrder]))[2]     # [value, index]
                rsidx[i,2] = findmin(abs.([expArray[x]-(rng[i]-less) for x in expOrder]))[2]     # [value, index]
                rsidx[i,3] = findmin(abs.([expArray[x]-(rng[i]+over) for x in expOrder]))[2]     # [value, index]
            end
        else
            for i=1:nr
                rsidx[i,1] = findmin(abs.([expArray[x]-rng[i] for x in expOrder]))[2]     # [value, index]
                rsidx[i,2] = findmin(abs.([expArray[x]-(rng[i]*(1-less)) for x in expOrder]))[2]     # [value, index]
                rsidx[i,3] = findmin(abs.([expArray[x]-(rng[i]*(1+over)) for x in expOrder]))[2]     # [value, index]
            end
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
    twpber = zeros(Float64, nr)  # total state/district-population weighted members of households by expenditure range
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

function categorizeHouseholdByIncomeByReligion(year,intv=[],normMode=0; sqrRt=false,absIntv=false,perCap=false,desOrd=false,popWgh=false,wghmode="district")
    # intv: proportions between invervals of highest to lowest
    # absIntv: if "true", then intv[] is a list of income values, descending order
    # normmode: [1]per capita, [2]per houehold, [3]basic information
    # perCap: [true] per capita mode, [false] per household mode
    # desOrd: [true] descening order of 'intv[]', [false] ascending order

    global wghSta, wghDis; if wghmode=="state"; wgh=wghSta elseif wghmode=="district"; wgh=wghDis end
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
    indInc = zeros(Int, nh)
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
            if accpop/totpop > intv[idx] && i < nh
                push!(pcidx, i)
                push!(incList, incArray[incOrder[i]])
                idx += 1
            end
            indInc[incOrder[i]] = idx
        end
    elseif absIntv; incList = copy(intv)
    else
        push!(incList, incArray[incOrder[1]])
        for i=1:length(intv); push!(incList, incArray[incOrder[trunc(Int, intv[i]*nh)]]) end
    end

    # set sector index
    if absIntv  # indInc: '1' for over intv[1], '3' for below intv[2], [2] for others
        for i=1:nh
            if incArray[i] >= intv[1]; indInc[i] = 1
            elseif incArray[i] < intv[2]; indInc[i] = 3
            else indInc[i] = 2
            end
        end
    elseif !perCap
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

function estimateEmissionCostByDistrict(year,expIntv=[],normMode=0; perCap=false,popWgh=false,desOrd=false,name=false,output="",exportFile=[],wghmode="district")
    # Estimate poverty alleviation leading CF (CO2 emissions) increase by district
    # considering expenditure-level, and distric consumption patterns
    # intv: proportions between invervals of highest to lowest
    # normmode: [1]per capita, [2]per houehold
    # perCap: [true] per capita mode, [false] per household mode

    global wghSta, wghDis; if wghmode=="state"; wgh=wghSta elseif wghmode=="district"; wgh=wghDis end
    global sec, hhid, siz, cat, dis, typ, inc, rel, pop, sam
    global staList, regList, catList, incList, relList, typList, disSta
    global emissionsHHs, emissionCostDis

    nh = length(hhid)
    ns = length(staList)
    nd = length(regList)
    nt = length(typList)+1
    nc = length(catList)
    ne = length(expIntv)

    # make index arrays
    idxSta = [findfirst(x->x==sta[hhid[i]], staList) for i=1:nh]    # state index
    idxDis = [findfirst(x->x==dis[hhid[i]], regList) for i=1:nh]    # district index
    idxTyp = [findfirst(x->x==typ[hhid[i]], typList) for i=1:nh]    # area type index

    # determine expenditure sections' starting index and values of income index
    incList = []
    pcidx = []  # current sector's starting index for 'per capita' emissions
    idxExp = zeros(Int, nh)
    expArray = [inc[h] for h in hhid]
    expOrder = sortperm(expArray, rev=desOrd)   # desOrd: [true] descening order, [false] ascending order

    if perCap   # determine sections if interval ratios are for 'per capita' emissions
        accpop = 0
        idx = 1
        push!(pcidx, idx)
        push!(incList, expArray[expOrder[idx]])
        if popWgh; totpop = sum([siz[h]*wgh[h] for h in hhid])
        else totpop = sum(collect(values(siz)))
        end
        for i=1:nh
            if popWgh; accpop += siz[hhid[expOrder[i]]] * wgh[hhid[expOrder[i]]]
            else accpop += siz[hhid[expOrder[i]]]
            end
            if accpop/totpop > expIntv[idx] && i < nh
                push!(pcidx, i)
                push!(incList, expArray[expOrder[i]])
                idx += 1
            end
            idxExp[expOrder[i]] = idx
        end
    else
        idx = 1
        push!(incList, expArray[expOrder[idx]])
        for i=1:ne
            push!(incList, expArray[expOrder[trunc(Int, expIntv[i]*nh)]])
            while idx <= trunc(Int, nh*expIntv[i])
                idxExp[idxOrder[idx]] = i
                idx += 1
            end
        end
        idxExp[idxOrder[nh]] = ne
    end

    # sum households and members by district
    hhsd = zeros(Float64, nd, nt, ne)   # total households by income level
    popd = zeros(Float64, nd, nt, ne)   # total members of households by income level
    wpopd = zeros(Float64, nd, nt, ne)  # total state-population weighted members of households by income level
    for i= 1:nh; hhsd[idxDis[i],idxTyp[i],idxExp[i]] += 1 end
    for i= 1:nh; popd[idxDis[i],idxTyp[i],idxExp[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpopd[idxDis[i],idxTyp[i],idxExp[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    for i= 1:nh; hhsd[idxDis[i],nt,idxExp[i]] += 1 end
    for i= 1:nh; popd[idxDis[i],nt,idxExp[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpopd[idxDis[i],nt,idxExp[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    # sum households and members by state
    hhss = zeros(Float64, ns, nt, ne)   # total households by income level
    pops = zeros(Float64, ns, nt, ne)   # total members of households by income level
    wpops = zeros(Float64, ns, nt, ne)  # total state-population weighted members of households by income level
    for i= 1:nh; hhss[idxSta[i],idxTyp[i],idxExp[i]] += 1 end
    for i= 1:nh; pops[idxSta[i],idxTyp[i],idxExp[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpops[idxSta[i],idxTyp[i],idxExp[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    for i= 1:nh; hhss[idxSta[i],nt,idxExp[i]] += 1 end
    for i= 1:nh; pops[idxSta[i],nt,idxExp[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpops[idxSta[i],nt,idxExp[i]] += siz[hhid[i]] * wgh[hhid[i]] end end

    # calculate average emission by the all categories
    eh = emissionsHHs[year]
    avgExpDis = zeros(Float64, nd, nt, ne, nc)     # average emission by district
    avgExpSta = zeros(Float64, ns, nt, ne, nc)     # average emission by state
    if popWgh
        for i=1:nh; avgExpDis[idxDis[i],idxTyp[i],idxExp[i],:] += eh[i,:] * wgh[hhid[i]] end
        for i=1:nh; avgExpDis[idxDis[i],nt,idxExp[i],:] += eh[i,:] * wgh[hhid[i]] end
        for i=1:nh; avgExpSta[idxSta[i],idxTyp[i],idxExp[i],:] += eh[i,:] * wgh[hhid[i]] end
        for i=1:nh; avgExpSta[idxSta[i],nt,idxExp[i],:] += eh[i,:] * wgh[hhid[i]] end
    else
        for i=1:nh; avgExpDis[idxDis[i],idxTyp[i],idxExp[i],:] += eh[i,:] end
        for i=1:nh; avgExpDis[idxDis[i],nt,idxExp[i],:] += eh[i,:] end
        for i=1:nh; avgExpSta[idxSta[i],idxTyp[i],idxExp[i],:] += eh[i,:] end
        for i=1:nh; avgExpSta[idxSta[i],nt,idxExp[i],:] += eh[i,:] end
    end
    if normMode == 1
        if popWgh
            for i=1:nc; avgExpDis[:,:,:,i] ./= wpopd end
            for i=1:nc; avgExpSta[:,:,:,i] ./= wpops end
        else
            for i=1:nc; avgExpDis[:,:,:,i] ./= popd end
            for i=1:nc; avgExpSta[:,:,:,i] ./= pops end
        end
    elseif normMode == 2
        for i=1:nc; avgExpDis[:,:,:,i] ./= hhsd end
        for i=1:nc; avgExpSta[:,:,:,i] ./= hhss end
    end

    # determine the average expenditures
    avg = zeros(Float64, nd, nt, ne, nc)   # determined average emission by all classifications
    for i=1:nd; for j=1:nt-1; for k=1:ne
        st = findfirst(x->x==disSta[regList[i]], staList)
        if avgExpDis[i,j,k,end]>0; avg[i,j,k,:] = avgExpDis[i,j,k,:]          # if there is district rural/urban data
        elseif avgExpDis[i,nt,k,end]>0; avg[i,j,k,:] = avgExpDis[i,nt,k,:]    # if there is no district rural/urban data
        elseif avgExpSta[st,j,k,end]>0; avg[i,j,k,:] = avgExpSta[st,j,k,:]    # if there is no district data but state rural/urban data
        elseif avgExpSta[st,nt,k,end]>0; avg[i,j,k,:] = avgExpSta[st,nt,k,:]  # if there is no state rural/urban data
        else println(regList[i]," ",typList[i]," ",expIntv[k]," does not have matching average emission data.")
        end
    end end end

    # estimate district emission cost per capita
    ecpc = zeros(Float64, ne-1, nd, nc)   # emission cost per capita or hhs, {alleviation target, district, category}
    if !desOrd
        if normMode == 1
            for i=2:ne; for j=1:nh; if idxExp[j]<i
                cost = avg[idxDis[j],idxTyp[j],i,:]*siz[hhid[j]] - eh[j,:]
                if cost[end]>0 ; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= sam[regList[i]][1] end
        else
            for i=2:ne; for j=1:nh; if idxExp[j]<i
                cost = avg[idxDis[j],idxTyp[j],i,:] - eh[j,:]
                if cost>0; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= sam[regList[i]][2] end
        end
    else
        if normMode == 1; for i=1:ne-1; for j=1:nh;
            if idxExp[j]>i
                cost = avg[idxDis[j],idxTyp[j],i,:]*siz[hhid[j]] - eh[j,:]
                if cost[end]>0 ; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= sam[regList[i]][1] end
        else for i=1:ne-1; for j=1:nh;
            if idxExp[j]>i
                cost = avg[idxDis[j],idxTyp[j],i,:] - eh[j,:]
                if cost>0; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= sam[regList[i]][2] end
        end
    end
    # estimate district total emission cost
    ec = zeros(Float64, ne-1, nd, nc)   # district emission cost, {alleviation target, district, category}
    if normMode == 1; for i=1:nd; ec[:,i,:] = ecpc[:,i,:] .* pop[regList[i]][1] end
    else for i=1:nd; ec[:,i,:] = ecpc[:,i,:] .* pop[regList[i]][2] end
    end

    #println("NaN: ", count(x->isnan(x),avg))

    # print CSV files
    if length(output)>0
        target = ne; if target>ne; target = ne end
        # total emission cost
        f = open(replace(output,".csv"=>"_total.csv"), "w")
        for i=1:target-1
            if !desOrd; println(f,"Target: ",expIntv[i+1]*100,"%") else println(f,"Target: ",expIntv[i]*100,"%") end
            print(f,"District"); for j=1:nc; print(f,",",catList[j]) end; println(f)
            for j=1:nd
                if name; print(f, nam[regList[j]]) else print(f, regList[j]) end
                for k=1:nc; print(f,",",ec[i,j,k]) end
                println(f)
            end
            println(f)
        end
        close(f)
        # emission cost per capita
        f = open(replace(output,".csv"=>"_perCap.csv"), "w")
        for i=1:target-1
            if !desOrd; println(f,"Target: ",expIntv[i+1]*100,"%") else println(f,"Target: ",expIntv[i]*100,"%") end
            print(f,"District"); for j=1:nc; print(f,",",catList[j]) end; println(f)
            for j=1:nd
                if name; print(f, nam[regList[j]]) else print(f, regList[j]) end
                for k=1:nc; print(f,",",ecpc[i,j,k]) end
                println(f)
            end
            println(f)
        end
        close(f)
    end

    # export map-making CSV files
    if length(exportFile)>0
        global ave, gid, gidData

        target = 2; target -= 1

        # making exporting table
        gidList = sort(unique(values(gid)))
        ng = length(gidList)

        tc = zeros(Float64, ng, nc)     # total emission cost
        tcpc = zeros(Float64, ng, nc)   # emission cost per capita
        spo = zeros(Float64, ng)   # number of sample population by district
        shh = zeros(Float64, ng)   # number of sample households by district
        tpo = zeros(Float64, ng)   # total number of population by district
        thh = zeros(Float64, ng)   # total number of households by district
        aec = zeros(Float64, ng)   # average expenditure per capita by district
        for i=1:nd
            idx = findfirst(x->x==gid[regList[i]],gidList)
            tc[idx,:] += ec[target,i,:]
            spo[idx] += sam[regList[i]][1]
            shh[idx] += sam[regList[i]][2]
            tpo[idx] += pop[regList[i]][1]
            thh[idx] += pop[regList[i]][2]
            aec[idx] += ave[regList[i]]*sam[regList[i]][1]
        end
        for i=1:ng; aec[i] /= spo[i] end
        if normMode==1; for i=1:ng; for j=1:nc; tcpc[i,j] = tc[i,j]/tpo[i] end end
        elseif normMode==2; for i=1:ng; for j=1:nc; tcpc[i,j] = tc[i,j]/thh[i] end end
        end

        # exporting total emission cost
        f = open(replace(exportFile[1],".csv"=>"_total.csv"), "w")
        print(f, exportFile[2]); for c in catList; print(f,",",c) end; println(f)
        for i = 1:size(tc, 1)
            print(f, gidList[i])
            for j = 1:size(tc, 2); print(f, ",", tc[i,j]) end
            println(f)
        end
        close(f)
        # exporting emission cost per capita
        f = open(replace(exportFile[1],".csv"=>"_perCap.csv"), "w")
        print(f, exportFile[2]); for c in catList; print(f,",",c) end; println(f)
        for i = 1:size(tcpc, 1)
            print(f, gidList[i])
            for j = 1:size(tcpc, 2); print(f, ",", tcpc[i,j]) end
            println(f)
        end
        close(f)

        gisDistrictEmissionCost[year] = tc
    end

    return ecpc, ec, hhsd, popd, wpopd, avgExpDis
end

function estimateEmissionCostByDistrictForThreshold(year,rng=[],normMode=0; stacked=false, absSpn=false,over=0.1,less=0.1,perCap=false,popWgh=false,desOrd=false,name=false,output="",exportFile=[],wghmode="district")
    # Estimate poverty alleviation leading CF (CO2 emissions) increase by district
    # considering expenditure-level, and distric consumption patterns
    # intv: proportions between invervals of highest to lowest
    # normmode: [1]per capita, [2]per houehold
    # perCap: [true] per capita mode, [false] per household mode

    global wghSta, wghDis; if wghmode=="state"; wgh=wghSta elseif wghmode=="district"; wgh=wghDis end
    global sec, hhid, siz, cat, dis, typ, inc, rel, pop, sam
    global staList, regList, catList, incList, relList, typList, disSta
    global emissionsHHs, emissionCostDis

    nh = length(hhid)
    ns = length(staList)
    nd = length(regList)
    nt = length(typList)+1
    nc = length(catList)
    nr = length(rng)
    ng = nr+1   # number of groups

    # make index arrays
    idxSta = [findfirst(x->x==sta[hhid[i]], staList) for i=1:nh]    # state index
    idxDis = [findfirst(x->x==dis[hhid[i]], regList) for i=1:nh]    # district index
    idxTyp = [findfirst(x->x==typ[hhid[i]], typList) for i=1:nh]    # area type index

    # determine each section's starting and ending indexes and values
    idxExp = zeros(Int, nh)

    expArray = [inc[h] for h in hhid]
    expOrder = sortperm(expArray, rev=desOrd)   # desOrd: [true] descening order, [false] ascending order
    incList = rng
    rsidx = zeros(Int, nr, 3)   # sector's standard, starting, and ending index: [standard, start, end]
    if absSpn
        for i=1:nr
            rsidx[i,1] = findmin(abs.([expArray[x]-rng[i] for x in expOrder]))[2]     # [value, index]
            rsidx[i,2] = findmin(abs.([expArray[x]-(rng[i]-less) for x in expOrder]))[2]     # [value, index]
            rsidx[i,3] = findmin(abs.([expArray[x]-(rng[i]+over) for x in expOrder]))[2]     # [value, index]
        end
    else
        for i=1:nr
            rsidx[i,1] = findmin(abs.([expArray[x]-rng[i] for x in expOrder]))[2]     # [value, index]
            rsidx[i,2] = findmin(abs.([expArray[x]-(rng[i]*(1-less)) for x in expOrder]))[2]     # [value, index]
            rsidx[i,3] = findmin(abs.([expArray[x]-(rng[i]*(1+over)) for x in expOrder]))[2]     # [value, index]
        end
    end
    # set each hh's group index
    for i=1:nh
        if expArray[i]>=rng[end]; idxExp[i]=ng
        else idxExp[i] = findfirst(x->x>expArray[i],rng)
        end
    end

    # sum households and members by districts by expenditure range
    hhsdg = zeros(Float64, nd, nt, ng)   # total households by enpenditure group
    popdg = zeros(Float64, nd, nt, ng)   # total members of households by enpenditure group
    wpopdg = zeros(Float64, nd, nt, ng)  # total state-population weighted members of households by enpenditure group
    for i= 1:nh; hhsdg[idxDis[i],idxTyp[i],idxExp[i]] += 1 end
    for i= 1:nh; popdg[idxDis[i],idxTyp[i],idxExp[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpopdg[idxDis[i],idxTyp[i],idxExp[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    for i= 1:nh; hhsdg[idxDis[i],nt,idxExp[i]] += 1 end
    for i= 1:nh; popdg[idxDis[i],nt,idxExp[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpopdg[idxDis[i],nt,idxExp[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    # sum households and members by state
    hhssg = zeros(Float64, ns, nt, ng)   # total households by enpenditure group
    popsg = zeros(Float64, ns, nt, ng)   # total members of households by enpenditure group
    wpopsg = zeros(Float64, ns, nt, ng)  # total state-population weighted members of households by enpenditure group
    for i= 1:nh; hhssg[idxSta[i],idxTyp[i],idxExp[i]] += 1 end
    for i= 1:nh; popsg[idxSta[i],idxTyp[i],idxExp[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpopsg[idxSta[i],idxTyp[i],idxExp[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    for i= 1:nh; hhssg[idxSta[i],nt,idxExp[i]] += 1 end
    for i= 1:nh; popsg[idxSta[i],nt,idxExp[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpopsg[idxSta[i],nt,idxExp[i]] += siz[hhid[i]] * wgh[hhid[i]] end end

    # sum households and members by district by expenditure threshold
    hhsdth = zeros(Float64, nd, nt, nr)   # total households by expenditure threshold by district
    popdth = zeros(Float64, nd, nt, nr)   # total members of households by expenditure threshold by district
    wpopdth = zeros(Float64, nd, nt, nr)  # total district/district-population weighted members of households by expenditure threshold by district
    hhssth = zeros(Float64, ns, nt, nr)   # total households by expenditure threshold by state
    popsth = zeros(Float64, ns, nt, nr)   # total members of households by expenditure threshold by state
    wpopsth = zeros(Float64, ns, nt, nr)  # total state/district-population weighted members of households by expenditure threshold by state
    for i=1:nr
        for j=rsidx[i,2]:rsidx[i,3]
            hidx = expOrder[j]
            hhsdth[idxDis[hidx],idxTyp[hidx],i] += 1
            popdth[idxDis[hidx],idxTyp[hidx],i] += siz[hhid[hidx]]
            if popWgh; wpopdth[idxDis[hidx],idxTyp[hidx],i] += siz[hhid[hidx]] * wgh[hhid[hidx]] end
            hhsdth[idxDis[hidx],nt,i] += 1
            popdth[idxDis[hidx],nt,i] += siz[hhid[hidx]]
            if popWgh; wpopdth[idxDis[hidx],nt,i] += siz[hhid[hidx]] * wgh[hhid[hidx]] end
            hhssth[idxSta[hidx],idxTyp[hidx],i] += 1
            popsth[idxSta[hidx],idxTyp[hidx],i] += siz[hhid[hidx]]
            if popWgh; wpopsth[idxSta[hidx],idxTyp[hidx],i] += siz[hhid[hidx]] * wgh[hhid[hidx]] end
            hhssth[idxSta[hidx],nt,i] += 1
            popsth[idxSta[hidx],nt,i] += siz[hhid[hidx]]
            if popWgh; wpopsth[idxSta[hidx],nt,i] += siz[hhid[hidx]] * wgh[hhid[hidx]] end
        end
    end

    # calculate average emission for each expenditure threshold
    eh = emissionsHHs[year]
    avgExpDis = zeros(Float64, nd, nt, nr, nc)     # average emission by district
    avgExpSta = zeros(Float64, ns, nt, nr, nc)     # average emission by state
    if popWgh
        for i=1:nr
            for j=rsidx[i,2]:rsidx[i,3]
                hidx = expOrder[j]
                avgExpDis[idxDis[hidx],idxTyp[hidx],i,:] += eh[hidx,:] * wgh[hhid[hidx]]
                avgExpDis[idxDis[hidx],nt,i,:] += eh[hidx,:] * wgh[hhid[hidx]]
                avgExpSta[idxSta[hidx],idxTyp[hidx],i,:] += eh[hidx,:] * wgh[hhid[hidx]]
                avgExpSta[idxSta[hidx],nt,i,:] += eh[hidx,:] * wgh[hhid[hidx]]
            end
        end
    else
        for i=1:nr
            for j=rsidx[i,2]:rsidx[i,3]
                hidx = expOrder[j]
                avgExpDis[idxDis[hidx],idxTyp[hidx],i,:] += eh[hidx,:]
                avgExpDis[idxDis[hidx],nt,i,:] += eh[hidx,:]
                avgExpSta[idxSta[hidx],idxTyp[hidx],i,:] += eh[hidx,:]
                avgExpSta[idxSta[hidx],nt,i,:] += eh[hidx,:]
            end
        end
    end
    if normMode == 1
        if popWgh
            for i=1:nc; avgExpDis[:,:,:,i] ./= wpopdth end
            for i=1:nc; avgExpSta[:,:,:,i] ./= wpopsth end
        else
            for i=1:nc; avgExpDis[:,:,:,i] ./= popdth end
            for i=1:nc; avgExpSta[:,:,:,i] ./= popsth end
        end
    elseif normMode == 2
        for i=1:nc; avgExpDis[:,:,:,i] ./= hhsdth end
        for i=1:nc; avgExpSta[:,:,:,i] ./= hhssth end
    end

    # determine the average expenditures
    avg = zeros(Float64, nd, nt, nr, nc)   # determined average emission by all classifications
    for i=1:nd; for j=1:nt-1; for k=1:nr
        st = findfirst(x->x==disSta[regList[i]], staList)
        if avgExpDis[i,j,k,end]>0; avg[i,j,k,:] = avgExpDis[i,j,k,:]          # if there is district rural/urban data
        elseif avgExpDis[i,nt,k,end]>0; avg[i,j,k,:] = avgExpDis[i,nt,k,:]    # if there is no district rural/urban data
        elseif avgExpSta[st,j,k,end]>0; avg[i,j,k,:] = avgExpSta[st,j,k,:]    # if there is no district data but state rural/urban data
        elseif avgExpSta[st,nt,k,end]>0; avg[i,j,k,:] = avgExpSta[st,nt,k,:]  # if there is no state rural/urban data
        else println(regList[i]," ",typList[i]," ",rng[k]," does not have matching average emission data.")
        end
    end end end

    # estimate district emission cost per capita
    ecpc = zeros(Float64, nr, nd, nc)   # emission cost per capita or hhs, {alleviation target, district, category}
    if !desOrd
        if normMode == 1
            if stacked
                for i=1:nr; for j=1:nh; if idxExp[j]<=i
                    cost = avg[idxDis[j],idxTyp[j],i,:]*siz[hhid[j]] - eh[j,:]
                    if cost[end]>0 ; ecpc[i,idxDis[j],:] += cost end
                end end end
            else
                for i=1:nr; for j=1:nh; if idxExp[j]==i
                    cost = avg[idxDis[j],idxTyp[j],i,:]*siz[hhid[j]] - eh[j,:]
                    if cost[end]>0 ; ecpc[i,idxDis[j],:] += cost end
                end end end
            end
            for i=1:nd; ecpc[:,i,:] ./= sam[regList[i]][1] end
        else
            if stacked
                for i=1:nr; for j=1:nh; if idxExp[j]<=i
                    cost = avg[idxDis[j],idxTyp[j],i,:] - eh[j,:]
                    if cost>0; ecpc[i,idxDis[j],:] += cost end
                end end end
            else
                for i=1:nr; for j=1:nh; if idxExp[j]==i
                    cost = avg[idxDis[j],idxTyp[j],i,:] - eh[j,:]
                    if cost>0; ecpc[i,idxDis[j],:] += cost end
                end end end
            end
            for i=1:nd; ecpc[:,i,:] ./= sam[regList[i]][2] end
        end
    else    # incomplete code
        if normMode == 1; for i=1:nr; for j=1:nh;
            if idxExp[j]>i
                cost = avg[idxDis[j],idxTyp[j],i,:]*siz[hhid[j]] - eh[j,:]
                if cost[end]>0 ; ecpc[i,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= sam[regList[i]][1] end
        else for i=1:nr; for j=1:nh;
            if idxExp[j]>i
                cost = avg[idxDis[j],idxTyp[j],i,:] - eh[j,:]
                if cost>0; ecpc[i,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= sam[regList[i]][2] end
        end
    end
    # estimate district total emission cost
    ec = zeros(Float64, nr, nd, nc)   # district emission cost, {alleviation target, district, category}
    if normMode == 1; for i=1:nd; ec[:,i,:] = ecpc[:,i,:] .* pop[regList[i]][1] end
    else for i=1:nd; ec[:,i,:] = ecpc[:,i,:] .* pop[regList[i]][2] end
    end

    #println("NaN: ", count(x->isnan(x),avg))

    # print CSV files
    if length(output)>0
        # total emission cost
        f = open(replace(output,".csv"=>"_total.csv"), "w")
        for i=1:nr
            println(f,"Target: ",rng[i])
            print(f,"District"); for j=1:nc; print(f,",",catList[j]) end
            print(f,",sample_hhs,sample_hhs,weighted_pop")
            println(f)
            for j=1:nd
                if name; print(f, nam[regList[j]]) else print(f, regList[j]) end
                for k=1:nc; print(f,",",ec[i,j,k]) end
                print(f,",",hhsdth[j,nt,i],",",popdth[j,nt,i],",",wpopdth[j,nt,i])
                println(f)
            end
            println(f)
        end
        close(f)
        # emission cost per capita
        f = open(replace(output,".csv"=>"_perCap.csv"), "w")
        for i=1:nr
            println(f,"Target: ",rng[i])
            print(f,"District"); for j=1:nc; print(f,",",catList[j]) end
            print(f,",sample_hhs,sample_hhs,weighted_pop")
            println(f)
            for j=1:nd
                if name; print(f, nam[regList[j]]) else print(f, regList[j]) end
                for k=1:nc; print(f,",",ecpc[i,j,k]) end
                print(f,",",hhsdth[j,nt,i],",",popdth[j,nt,i],",",wpopdth[j,nt,i])
                println(f)
            end
            println(f)
        end
        close(f)
    end

    # export map-making CSV files
    if length(exportFile)>0
        global ave, gid, gidData

        target = 1

        # making exporting table
        gidList = sort(unique(values(gid)))
        ng = length(gidList)

        tc = zeros(Float64, ng, nc)     # total emission cost
        tcpc = zeros(Float64, ng, nc)   # emission cost per capita
        spo = zeros(Float64, ng)   # number of sample population by district
        shh = zeros(Float64, ng)   # number of sample households by district
        tpo = zeros(Float64, ng)   # total number of population by district
        thh = zeros(Float64, ng)   # total number of households by district
        aec = zeros(Float64, ng)   # average expenditure per capita by district
        for i=1:nd
            idx = findfirst(x->x==gid[regList[i]],gidList)
            tc[idx,:] += ec[target,i,:]
            spo[idx] += sam[regList[i]][1]
            shh[idx] += sam[regList[i]][2]
            tpo[idx] += pop[regList[i]][1]
            thh[idx] += pop[regList[i]][2]
            aec[idx] += ave[regList[i]]*sam[regList[i]][1]
        end
        for i=1:ng; aec[i] /= spo[i] end
        if normMode==1; for i=1:ng; for j=1:nc; tcpc[i,j] = tc[i,j]/tpo[i] end end
        elseif normMode==2; for i=1:ng; for j=1:nc; tcpc[i,j] = tc[i,j]/thh[i] end end
        end

        # exporting total emission cost
        f = open(replace(exportFile[1],".csv"=>"_total.csv"), "w")
        print(f, exportFile[2]); for c in catList; print(f,",",c) end; println(f)
        for i = 1:size(tc, 1)
            print(f, gidList[i])
            for j = 1:size(tc, 2); print(f, ",", tc[i,j]) end
            println(f)
        end
        close(f)
        # exporting emission cost per capita
        f = open(replace(exportFile[1],".csv"=>"_perCap.csv"), "w")
        print(f, exportFile[2]); for c in catList; print(f,",",c) end; println(f)
        for i = 1:size(tcpc, 1)
            print(f, gidList[i])
            for j = 1:size(tcpc, 2); print(f, ",", tcpc[i,j]) end
            println(f)
        end
        close(f)

        gisDistrictEmissionCost[year] = tc
    end

    return ecpc, ec
end

function estimateEmissionCostByDistrictByReligion(year, expIntv=[], normMode=0; absIntv=false, perCap=false, popWgh=false, desOrd=false, name=false, output="", exportFile=[], wghmode="district")
    # Estimate poverty alleviation leading CF (CO2 emissions) increase by district
    # considering religion, expenditure-level, and distric consumption patterns
    # intv: proportions between invervals of highest to lowest
    # normmode: [1]per capita, [2]per houehold
    # perCap: [true] per capita mode, [false] per household mode

    global wghSta, wghDis; if wghmode=="state"; wgh=wghSta elseif wghmode=="district"; wgh=wghDis end
    global sec, hhid, siz, cat, dis, typ, inc, rel, pop, sam
    global staList, regList, catList, incList, relList, typList, disSta
    global emissionsHHs, emissionCostDis

    nh = length(hhid)
    ns = length(staList)
    nd = length(regList)
    nt = length(typList)+1
    nc = length(catList)
    ne = length(expIntv)
    nr = length(relList)

    # make index arrays
    idxSta = [findfirst(x->x==sta[hhid[i]], staList) for i=1:nh]    # state index
    idxDis = [findfirst(x->x==dis[hhid[i]], regList) for i=1:nh]    # district index
    idxTyp = [findfirst(x->x==typ[hhid[i]], typList) for i=1:nh]    # area type index
    idxRel = [findfirst(x->x==rel[hhid[i]], relList) for i=1:nh]    # religion index

    # determine expenditure sections' starting index and values of income index
    pcidx = []  # current sector's starting index for 'per capita' emissions
    incList = []
    idxExp = zeros(Int, nh)
    expArray = [inc[h] for h in hhid]
    expOrder = sortperm(expArray, rev=desOrd)   # desOrd: [true] descening order, [false] ascending order
    if perCap   # determine sections if interval ratios are for 'per capita' emissions
        accpop = 0
        idx = 1
        push!(pcidx, idx)
        push!(incList, expArray[expOrder[idx]])
        if popWgh; totpop = sum([siz[h]*wgh[h] for h in hhid])
        else totpop = sum(collect(values(siz)))
        end
        for i=1:nh
            if popWgh; accpop += siz[hhid[expOrder[i]]] * wgh[hhid[expOrder[i]]]
            else accpop += siz[hhid[expOrder[i]]]
            end
            if accpop/totpop > expIntv[idx] && i < nh
                push!(pcidx, i)
                push!(incList, expArray[expOrder[i]])
                idx += 1
            end
            idxExp[expOrder[i]] = idx
        end
    else
        idx = 1
        push!(incList, expArray[expOrder[idx]])
        for i=1:ne
            push!(incList, expArray[expOrder[trunc(Int, expIntv[i]*nh)]])
            while idx <= trunc(Int, nh*expIntv[i])
                idxExp[idxOrder[idx]] = i
                idx += 1
            end
        end
        idxExp[idxOrder[nh]] = ne
    end

    # sum households and members by all classifications
    hhsc = zeros(Float64, nd, nt, ne, nr)   # total households by income level
    popc = zeros(Float64, nd, nt, ne, nr)   # total members of households by income level
    wpopc = zeros(Float64, nd, nt, ne, nr)  # total state-population weighted members of households by income level
    for i= 1:nh; hhsc[idxDis[i],idxTyp[i],idxExp[i],idxRel[i]] += 1 end
    for i= 1:nh; popc[idxDis[i],idxTyp[i],idxExp[i],idxRel[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpopc[idxDis[i],idxTyp[i],idxExp[i],idxRel[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    for i= 1:nh; hhsc[idxDis[i],nt,idxExp[i],idxRel[i]] += 1 end
    for i= 1:nh; popc[idxDis[i],nt,idxExp[i],idxRel[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpopc[idxDis[i],nt,idxExp[i],idxRel[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    # sum households and members by district
    hhsd = zeros(Float64, nd, nt, ne)   # total households by income level
    popd = zeros(Float64, nd, nt, ne)   # total members of households by income level
    wpopd = zeros(Float64, nd, nt, ne)  # total state-population weighted members of households by income level
    for i= 1:nh; hhsd[idxDis[i],idxTyp[i],idxExp[i]] += 1 end
    for i= 1:nh; popd[idxDis[i],idxTyp[i],idxExp[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpopd[idxDis[i],idxTyp[i],idxExp[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    for i= 1:nh; hhsd[idxDis[i],nt,idxExp[i]] += 1 end
    for i= 1:nh; popd[idxDis[i],nt,idxExp[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpopd[idxDis[i],nt,idxExp[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    # sum households and members by state
    hhss = zeros(Float64, ns, nt, ne)   # total households by income level
    pops = zeros(Float64, ns, nt, ne)   # total members of households by income level
    wpops = zeros(Float64, ns, nt, ne)  # total state-population weighted members of households by income level
    for i= 1:nh; hhss[idxSta[i],idxTyp[i],idxExp[i]] += 1 end
    for i= 1:nh; pops[idxSta[i],idxTyp[i],idxExp[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpops[idxSta[i],idxTyp[i],idxExp[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    for i= 1:nh; hhss[idxSta[i],nt,idxExp[i]] += 1 end
    for i= 1:nh; pops[idxSta[i],nt,idxExp[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpops[idxSta[i],nt,idxExp[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    hhssr = zeros(Float64, ns, nt, ne, nr)   # total households by income level
    popsr = zeros(Float64, ns, nt, ne, nr)   # total members of households by income level
    wpopsr = zeros(Float64, ns, nt, ne, nr)  # total state-population weighted members of households by income level
    for i= 1:nh; hhssr[idxSta[i],idxTyp[i],idxExp[i],idxRel[i]] += 1 end
    for i= 1:nh; popsr[idxSta[i],idxTyp[i],idxExp[i],idxRel[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpopsr[idxSta[i],idxTyp[i],idxExp[i],idxRel[i]] += siz[hhid[i]] * wgh[hhid[i]] end end
    for i= 1:nh; hhssr[idxSta[i],nt,idxExp[i],idxRel[i]] += 1 end
    for i= 1:nh; popsr[idxSta[i],nt,idxExp[i],idxRel[i]] += siz[hhid[i]] end
    if popWgh; for i= 1:nh; wpopsr[idxSta[i],nt,idxExp[i],idxRel[i]] += siz[hhid[i]] * wgh[hhid[i]] end end

    # calculate average emission by the all categories
    eh = emissionsHHs[year]
    avgExp = zeros(Float64, nd, nt, ne, nr, nc)   # average emission by all classifications
    avgExpDis = zeros(Float64, nd, nt, ne, nc)     # average emission by district
    avgExpSta = zeros(Float64, ns, nt, ne, nc)     # average emission by state
    avgExpStaRel = zeros(Float64, ns, nt, ne, nr, nc)     # average emission by state by religion
    if popWgh
        for i=1:nh; avgExp[idxDis[i],idxTyp[i],idxExp[i],idxRel[i],:] += eh[i,:] * wgh[hhid[i]] end
        for i=1:nh; avgExp[idxDis[i],nt,idxExp[i],idxRel[i],:] += eh[i,:] * wgh[hhid[i]] end
        for i=1:nh; avgExpDis[idxDis[i],idxTyp[i],idxExp[i],:] += eh[i,:] * wgh[hhid[i]] end
        for i=1:nh; avgExpDis[idxDis[i],nt,idxExp[i],:] += eh[i,:] * wgh[hhid[i]] end
        for i=1:nh; avgExpSta[idxSta[i],idxTyp[i],idxExp[i],:] += eh[i,:] * wgh[hhid[i]] end
        for i=1:nh; avgExpSta[idxSta[i],nt,idxExp[i],:] += eh[i,:] * wgh[hhid[i]] end
        for i=1:nh; avgExpStaRel[idxSta[i],idxTyp[i],idxExp[i],idxRel[i],:] += eh[i,:] * wgh[hhid[i]] end
        for i=1:nh; avgExpStaRel[idxSta[i],nt,idxExp[i],idxRel[i],:] += eh[i,:] * wgh[hhid[i]] end
    else
        for i=1:nh; avgExp[idxDis[i],idxTyp[i],idxExp[i],idxRel[i],:] += eh[i,:] end
        for i=1:nh; avgExp[idxDis[i],nt,idxExp[i],idxRel[i],:] += eh[i,:] end
        for i=1:nh; avgExpDis[idxDis[i],idxTyp[i],idxExp[i],:] += eh[i,:] end
        for i=1:nh; avgExpDis[idxDis[i],nt,idxExp[i],:] += eh[i,:] end
        for i=1:nh; avgExpSta[idxSta[i],idxTyp[i],idxExp[i],:] += eh[i,:] end
        for i=1:nh; avgExpSta[idxSta[i],nt,idxExp[i],:] += eh[i,:] end
        for i=1:nh; avgExpStaRel[idxSta[i],idxTyp[i],idxExp[i],idxRel[i],:] += eh[i,:] end
        for i=1:nh; avgExpStaRel[idxSta[i],nt,idxExp[i],idxRel[i],:] += eh[i,:] end
    end
    if normMode == 1
        if popWgh
            for i=1:nc; avgExp[:,:,:,:,i] ./= wpopc end
            for i=1:nc; avgExpDis[:,:,:,i] ./= wpopd end
            for i=1:nc; avgExpSta[:,:,:,i] ./= wpops end
            for i=1:nc; avgExpStaRel[:,:,:,:,i] ./= wpopsr end
        else
            for i=1:nc; avgExp[:,:,:,:,i] ./= popc end
            for i=1:nc; avgExpDis[:,:,:,i] ./= popd end
            for i=1:nc; avgExpSta[:,:,:,i] ./= pops end
            for i=1:nc; avgExpStaRel[:,:,:,:,i] ./= popsr end
        end
    elseif normMode == 2
        for i=1:nc; avgExp[:,:,:,:,i] ./= hhsc end
        for i=1:nc; avgExpDis[:,:,:,i] ./= hhsd end
        for i=1:nc; avgExpSta[:,:,:,i] ./= hhss end
        for i=1:nc; avgExpStaRel[:,:,:,:,i] ./= hhssr end
    end

    # determine the average expenditures
    avg = zeros(Float64, nd, nt, ne, nr, nc)   # determined average emission by all classifications
    for i=1:nd; for j=1:nt-1; for k=1:ne; for l=1:nr
        st = findfirst(x->x==disSta[regList[i]], staList)
        if avgExp[i,j,k,l,end]>0; avg[i,j,k,l,:] = avgExp[i,j,k,l,:]    # if there is district religion-rural/urban data
        elseif avgExpDis[i,j,k,end]>0; avg[i,j,k,l,:] = avgExpDis[i,j,k,:]  # if there is no district religion, but rural/urban data
        elseif avgExpDis[i,nt,k,end]>0; avg[i,j,k,l,:] = avgExpDis[i,nt,k,:]    # if there is no district religion nor rural/urban data
        elseif avgExpStaRel[st,j,k,l,end]>0; avg[i,j,k,l,:] = avgExpStaRel[st,j,k,l,:]  # if there is state religion-rural/urban data
        elseif avgExpSta[st,j,k,end]>0; avg[i,j,k,l,:] = avgExpSta[st,j,k,:]    # if there is no state religion, but rural/urban data
        elseif avgExpSta[st,nt,k,end]>0; avg[i,j,k,l,:] = avgExpSta[st,nt,k,:]  # if there is no state religion nor rural/urban data
        else println(regList[i]," ",typList[i]," ",expIntv[k]," ",relList[l]," does not have matching average emission data.")
    end end end end end

    # estimate district emission cost per capita
    ecpc = zeros(Float64, ne-1, nd, nc)   # emission cost per capita or hhs, {alleviation target, district, category}
    if !desOrd
        if normMode == 1
            for i=2:ne; for j=1:nh; if idxExp[j]<i
                cost = avg[idxDis[j],idxTyp[j],i,idxRel[j],:]*siz[hhid[j]] - eh[j,:]
                if cost[end]>0 ; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= sam[regList[i]][1] end
        else
            for i=2:ne; for j=1:nh; if idxExp[j]<i
                cost = avg[idxDis[j],idxTyp[j],i,idxRel[j],:] - eh[j,:]
                if cost>0; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= sam[regList[i]][2] end
        end
    else
        if normMode == 1; for i=1:ne-1; for j=1:nh;
            if idxExp[j]>i
                cost = avg[idxDis[j],idxTyp[j],i,idxRel[j],:]*siz[hhid[j]] - eh[j,:]
                if cost[end]>0 ; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= sam[regList[i]][1] end
        else for i=1:ne-1; for j=1:nh;
            if idxExp[j]>i
                cost = avg[idxDis[j],idxTyp[j],i,idxRel[j],:] - eh[j,:]
                if cost>0; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= sam[regList[i]][2] end
        end
    end
    # estimate district total emission cost
    ec = zeros(Float64, ne-1, nd, nc)   # district emission cost, {alleviation target, district, category}
    if normMode == 1; for i=1:nd; ec[:,i,:] = ecpc[:,i,:] .* pop[regList[i]][1] end
    else for i=1:nd; ec[:,i,:] = ecpc[:,i,:] .* pop[regList[i]][2] end
    end

    #println("NaN: ", count(x->isnan(x),avg))

    # print CSV files
    if length(output)>0
        target = ne; if target>ne; target = ne end
        # total emission cost
        f = open(replace(output,".csv"=>"_total.csv"), "w")
        for i=1:target-1
            if !desOrd; println(f,"Target: ",expIntv[i+1]*100,"%") else println(f,"Target: ",expIntv[i]*100,"%") end
            print(f,"District"); for j=1:nc; print(f,",",catList[j]) end; println(f)
            for j=1:nd
                if name; print(f, nam[regList[j]]) else print(f, regList[j]) end
                for k=1:nc; print(f,",",ec[i,j,k]) end
                println(f)
            end
            println(f)
        end
        close(f)
        # emission cost per capita
        f = open(replace(output,".csv"=>"_perCap.csv"), "w")
        for i=1:target-1
            if !desOrd; println(f,"Target: ",expIntv[i+1]*100,"%") else println(f,"Target: ",expIntv[i]*100,"%") end
            print(f,"District"); for j=1:nc; print(f,",",catList[j]) end; println(f)
            for j=1:nd
                if name; print(f, nam[regList[j]]) else print(f, regList[j]) end
                for k=1:nc; print(f,",",ecpc[i,j,k]) end
                println(f)
            end
            println(f)
        end
        close(f)
    end

    # export map-making CSV files
    if length(exportFile)>0
        global ave, gid, gidData

        target = 2; target -= 1

        # making exporting table
        gidList = sort(unique(values(gid)))
        ng = length(gidList)

        tc = zeros(Float64, ng, nc)     # total emission cost
        tcpc = zeros(Float64, ng, nc)   # emission cost per capita
        spo = zeros(Float64, ng)   # number of sample population by district
        shh = zeros(Float64, ng)   # number of sample households by district
        tpo = zeros(Float64, ng)   # total number of population by district
        thh = zeros(Float64, ng)   # total number of households by district
        aec = zeros(Float64, ng)   # average expenditure per capita by district
        for i=1:nd
            idx = findfirst(x->x==gid[regList[i]],gidList)
            tc[idx,:] += ec[target,i,:]
            spo[idx] += sam[regList[i]][1]
            shh[idx] += sam[regList[i]][2]
            tpo[idx] += pop[regList[i]][1]
            thh[idx] += pop[regList[i]][2]
            aec[idx] += ave[regList[i]]*sam[regList[i]][1]
        end
        for i=1:ng; aec[i] /= spo[i] end
        if normMode==1; for i=1:ng; for j=1:nc; tcpc[i,j] = tc[i,j]/tpo[i] end end
        elseif normMode==2; for i=1:ng; for j=1:nc; tcpc[i,j] = tc[i,j]/thh[i] end end
        end

        # exporting total emission cost
        f = open(replace(exportFile[1],".csv"=>"_total.csv"), "w")
        print(f, exportFile[2]); for c in catList; print(f,",",c) end; println(f)
        for i = 1:size(tc, 1)
            print(f, gidList[i])
            for j = 1:size(tc, 2); print(f, ",", tc[i,j]) end
            println(f)
        end
        close(f)
        # exporting emission cost per capita
        f = open(replace(exportFile[1],".csv"=>"_perCap.csv"), "w")
        print(f, exportFile[2]); for c in catList; print(f,",",c) end; println(f)
        for i = 1:size(tcpc, 1)
            print(f, gidList[i])
            for j = 1:size(tcpc, 2); print(f, ",", tcpc[i,j]) end
            println(f)
        end
        close(f)

        gisDistrictEmissionCost[year] = tc
    end

    return ecpc, ec, hhsc, popc, wpopc, avgExp
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

function printEmissionByIncome(year, outputFile, intv=[], tpbi=[], thbi=[], twpbi=[]; absIntv=false, desOrd=false, relative=false)

    global catList, incList, emissionsInc
    ei = emissionsInc[year]
    ni = length(intv); if absIntv; ni += 1 end
    nc = length(catList)

    f = open(outputFile, "w")
    print(f,"Exp_Lv,Value")
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
        else
            if i==1; print(f, "Bottom ", round(intv[i]*100,digits=1),"%,Less ",round(incList[i+1],digits=2))
            elseif i==ni; print(f, "Top ", round((intv[i]-intv[i-1])*100,digits=1),"%,Over ",round(incList[i],digits=2))
            else
                print(f, round(intv[i-1]*100,digits=1),"-",round(intv[i]*100,digits=1),"%")
                print(f,",",round(incList[i],digits=2),"-",round(incList[i+1],digits=2))
            end
        end

        for j = 1:nc; print(f, ",", ei[i,j]) end
        if length(tpbi)>0; print(f,",",tpbi[i]) end
        if length(thbi)>0; print(f,",",thbi[i]) end
        if length(twpbi)>0; print(f,",",twpbi[i]) end
        println(f)
    end
    close(f)

    if relative
        rel = zeros(Float64,ni,nc)
        medidx = Int(round((1+ni)/2,digits=0))
        for i=1:ni; for j=1:nc; rel[i,j] = ei[i,j]/ei[medidx,j]*100 end end
        f = open(replace(outputFile,".csv"=>"_relative.csv"), "w")
        print(f,"Exp_Lv,Value")
        for c in catList; print(f, ",", c) end
        println(f)
        for i = 1:ni
            if i==1; print(f, "Bottom ", round(intv[i]*100,digits=1),"%,Less ",round(incList[i+1],digits=2))
            elseif i==ni; print(f, "Top ", round((intv[i]-intv[i-1])*100,digits=1),"%,Over ",round(incList[i],digits=2))
            else
                print(f, round(intv[i-1]*100,digits=1),"-",round(intv[i]*100,digits=1),"%")
                print(f,",",round(incList[i],digits=2),"-",round(incList[i+1],digits=2))
            end
            for j = 1:nc; print(f, ",", rel[i,j]) end
            println(f)
        end
        close(f)
    end

end

function printEmissionByRange(year, outputFile, rsidx=[], thber=[], tpber=[], twpber=[], order=[])

    global catList, emissionsRng
    er = emissionsRng[year]
    nr = size(rsidx,1)

    f = open(outputFile, "w")

    print(f,"Expenditure,Range")
    for c in catList; print(f, ",", c) end
    if length(tpber)>0; print(f,",Pop.") end
    if length(thber)>0; print(f,",HH.") end
    if length(twpber)>0; print(f,",WghPop.") end
    println(f)
    for i = 1:nr
        print(f,round(inc[hhid[order[rsidx[i,1]]]],digits=2))
        print(f,",",round(inc[hhid[order[rsidx[i,2]]]],digits=2),"-",round(inc[hhid[order[rsidx[i,3]]]],digits=2))
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

function printEmissionByDistEmLev(year, outputFile, intv=[])

    global catList, regList, emissionsDisLev
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

    global catList, regList, emissionsLev
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

end
