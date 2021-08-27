module EmissionCategorizer

# Developed date: 3. Aug. 2020
# Last modified date: 27. Aug. 2021
# Subject: Categorize EU households' carbon footprints
# Description: Read household-level CFs and them by consumption category, district, expenditure-level, and etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using XLSX
using Statistics
using Formatting: printfmt

hhsList = Dict{Int, Dict{String, Array{String, 1}}}()   # Household ID: {year, {nation, {hhid}}}
sec = Array{String, 1}()                # Consumption products' or services' sectors
secName = Dict{String, String}()        # sector name dictionary: {sector code, name}
# deSec = Dict{Int, Array{String, 1}}()                   # direct emission causing consumption sectors: {year, {DE sector}}
# deCodes = Dict{Int, Dict{String, Array{String, 1}}}()   # direct emission related consumption sectors: {year, {DE sector, {HBS sector}}}
# deHbsList = Dict{Int, Array{String , 1}}()              # direct emission corresponding HBS sector list: {year, {HBS sector}}
# deCat = Dict{Int, Dict{String, String}}()               # direct emission category linkages: {year, {DE sector, DE category}}

# hhid -> nation(2digit)_hhid: some HHIDs are duplicated across multiple countries
cat = Dict{Int, Dict{String, String}}()         # category dictionary: {sector code, category}
nat = Dict{Int, Dict{String, String}}()         # hhid's nation: {year, {hhid, nation code}}
reg = Dict{Int, Dict{String, String}}()         # hhid's NUTS: {year, {hhid, NUTS code}}
typ = Dict{Int, Dict{String, String}}()         # hhid's sector type, urban or rural: {year, {hhid, "urban" or "rural"}}
siz = Dict{Int, Dict{String, Int}}()            # hhid's family size: {year, {hhid, number of members}}
eqs = Dict{Int, Dict{String, Float64}}()        # hhid's family equivalent size (OECD scale): {year, {hhid, number of members}}
meqs = Dict{Int, Dict{String, Float64}}()       # hhid's family equivalent size (modified OECD scale): {year, {hhid, number of members}}
inc = Dict{Int, Dict{String, Float64}}()        # hhid's income: {year, {hhid, total income}}
exp = Dict{Int, Dict{String, Float64}}()        # hhid's domestic expenditure: {year, {hhid, total domestic expenditure}}
pds = Dict{Int, Dict{String, Int}}()            # hhid region's population density: {year, {hhid, district's population density}}
rel = Dict{Int, Dict{String, Int}}()            # hhid's religion: {year, {hhid, religion code}}
wgh = Dict{Int, Dict{String, Float64}}()        # hhid's weight: {year, {hhid, weight}}
wghNuts = Dict{Int, Dict{String, Float64}}()    # hhid's NUTS weight: {year, {hhid, weight}}

nutsLv = 0      # NUTS level
nuts = Dict{Int, Dict{String, String}}()       # NUTS: {year, {code, label}}
nutsList = Dict{Int, Dict{String, Array{String, 1}}}()      # NUTS code list: {year, {nation_code, {NUTS_code}}}
pop = Dict{Int, Dict{String, Float64}}()        # Population: {year, {NUTS_code, population}}
popList = Dict{Int, Dict{String, Dict{String, Float64}}}()  # Population list: {year, {nation_code, {NUTS_code, population}}}
poplb = Dict{Int, Dict{String, String}}()       # populaton NUTS label: {year, {NUTS_code, NUTS_label}}

directCE = Dict{Int, Dict{String, Array{Float64, 2}}}()     # direct carbon emission: {year, {nation, {table}}}
indirectCE = Dict{Int, Dict{String, Array{Float64, 2}}}()   # indirect carbon emission: {year, {nation, {table}}}
integratedCF = Dict{Int, Dict{String, Array{Float64, 2}}}() # carbon footprint: {year, {nation, {table}}}

yrList = Array{Int, 1}()        # year list
catList = Array{String, 1}()    # category list
deCatList = Dict{Int, Array{String, 1}}()  # DE category list: {year, {DE category}}
natList = Array{String, 1}()    # nation list
natName = Dict{String,String}() # nation's code and full-name: {nation code, full-name}
natA3 = Dict{String,String}()   # nation's A3: {nation code, A3}

ieHHs = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied indirect emission by household: {year, {nation, {hhid, category}}}
ieReg = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied indirect emission by region: {year, {nation, {region, category}}}
ieRegDiff = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied indirect emission differences by region: {year, {nation, {region, category}}}
deHHs = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied direct emission by household: {year, {nation, {hhid, category}}}
deReg = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied direct emission by region: {year, {nation, {region, category}}}
deRegDiff = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied direct emission differences by region: {year, {nation, {region, category}}}
cfHHs = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied carbon footprint by household: {year, {nation, {hhid, category}}}
cfReg = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied carbon footprint by region: {year, {nation, {region, category}}}
cfRegDiff = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied carbon footprint differences by region: {year, {nation, {region, category}}}

ntpop = Dict{Int,Dict{String,Dict{String,Array{Int,1}}}}()      # NUTS population:{year,{nation,{NUTS,population{total,dense,mid,sparse,none}}}}
ntsmp = Dict{Int,Dict{String,Dict{String,Array{Int,1}}}}()      # NUTS sample size:{year,{nation,{NUTS,sample number{total,dense,mid,sparse,none}}}}
ntwgh = Dict{Int,Dict{String,Dict{String,Array{Float64,1}}}}()  # NUTS population weight:{year,{nation,{NUTS,population{total,dense,mid,sparse,none}}}}


# GIS data
gisNutsList = Dict{Int, Array{String, 1}}()           # NUTS list: {year, {region(hbscd)}}
gisRegionalIe = Dict{Int, Array{Float64, 2}}()        # categozied indirect emission by district: {year, {region(hbscd), category}}
gisRegionalIeRank = Dict{Int, Array{Int, 2}}()        # categozied indirect emission rank by district: {year, {region(hbscd), category}}
gisRegionalIePerCap = Dict{Int, Array{Float64, 2}}()  # categozied indirect emission per capita by district: {year, {region(hbscd), category}}
gisRegionalIeRankPerCap = Dict{Int, Array{Int, 2}}()  # categozied indirect emission per capita rank by district: {year, {region(hbscd), category}}
gisRegionalDe = Dict{Int, Array{Float64, 2}}()        # categozied direct emission by district: {year, {region(hbscd), DE category}}
gisRegionalDePerCap = Dict{Int, Array{Float64, 2}}()  # categozied direct emission per capita by district: {year, {region(hbscd), DE category}}
gisRegionalDeRank = Dict{Int, Array{Int, 2}}()        # categozied direct emission rank by district: {year, {region(hbscd), DE category}}
gisRegionalDeRankPerCap = Dict{Int, Array{Int, 2}}()  # categozied direct emission per capita rank by district: {year, {region(hbscd), DE category}}

gisRegionalCF = Dict{Int, Array{Float64, 2}}()        # categozied carbon footprint by district: {year, {region(hbscd), category}}
gisRegionalCFrank = Dict{Int, Array{Int, 2}}()        # categozied carbon footprint rank by district: {year, {region(hbscd), category}}
gisRegionalCFperCap = Dict{Int, Array{Float64, 2}}()  # categozied carbon footprint per capita by district: {year, {region(hbscd), category}}
gisRegionalCFrankPerCap = Dict{Int, Array{Int, 2}}()  # categozied carbon footprint per capita rank by district: {year, {region(hbscd), category}}
gisRegionalCFdiff = Dict{Int, Array{Float64, 2}}()    # differences of categozied carbon footprint by district: (emission-mean)/mean, {year, {district(GID), category}}
gisRegionalCFdiffRank = Dict{Int, Array{Int, 2}}()    # difference ranks of categozied carbon footprint by district: (emission-mean)/mean, {year, {district(GID), category}}
gisRegionalCFdiffPerCap = Dict{Int, Array{Float64, 2}}()  # differences of categozied carbon footprint per capita by district: (emission-mean)/mean, {year, {district(GID), category}}
gisRegionalCFdiffRankPerCap = Dict{Int, Array{Int, 2}}()  # difference ranks of categozied carbon footprint per capita by district: (emission-mean)/mean, {year, {district(GID), category}}

gisTotPop = Dict{Int, Array{Float64, 1}}()      # GIS version, total population by NUTS
gisSamPop = Dict{Int, Array{Float64, 1}}()      # GIS version, total sample members by NUTS
gisAvgExp = Dict{Int, Array{Float64, 1}}()      # GIS version, average expenditure by NUTS

regList = Dict{Int, Array{String, 1}}() # district list

sam = Dict{Int, Dict{String, Tuple{Int,Int}}}() # sample population and households by districct: {district code, (population, number of households)}
ave = Dict{Int, Dict{String, Float64}}()        # average annual expenditure per capita, USD/yr: {district code, mean Avg.Exp./cap/yr}

popcd = Dict{Int, Dict{String, String}}()       # concordance NUTS code: {year, {NUTS code, replaced population NUTS code}}
pophbscd = Dict{Int, Dict{String, String}}()    # concordance NUTS code: {year, {population NUTS code, HBS NUTS code}}
hbscd = Dict{Int, Dict{String, String}}()       # concordance NUTS code: {year, {NUTS code, HBS NUTS code}}
popgiscd = Dict{Int, Dict{String, String}}()    # concordance NUTS code: {year, {population NUTS code, GIS NUTS code}}
popcdlist = Dict{Int, Array{String, 1}}()       # Population NUTS code list: {year, Population NUTS code}
ntcdlist = Dict{Int, Array{String, 1}}()        # NUTS code list: {year, NUTS code}
gispopcdlist = Dict{Int, Dict{String, Array{String, 1}}}()   # Population NUTS code list: {year, {GIS_NUTS coded, {Population NUTS code}}}
giscdlist = Dict{Int, Array{String, 1}}()       # GIS NUTS code list: {year, GIS NUTS code}
hbspopcdlist = Dict{Int, Dict{String, Array{String, 1}}}()   # Population NUTS code list: {year, {HBS_NUTS coded, {Population NUTS code}}}
hbscdlist = Dict{Int, Array{String, 1}}()       # HBS NUTS code list: {year, HBS NUTS code}
majorCity = Dict{Int, Dict{String, String}}()   # major city of NUTS: {year, {NUTS_upper, major sub-NUTS}}
gisCoord = Dict{Int, Dict{String, Tuple{Float64, Float64}}}()   # GIS coordinates: {year, {GIS_NUTS, {X, Y}}}
giscatlab = Dict{String, String}()              # Web-exporting category matching table: {Category label in program, in web-site files}

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

function makeNationalSummary(year, outputFile; nuts_mode=false)

    global hhsList, natList, natName, siz, wgh, wghNuts, popList
    global indirectCE, directCE, integratedCF
    if nuts_mode; w = wghNuts else w = wgh end

    nn = length(natList)
    natsam = zeros(Int, nn)
    nateqs = zeros(Float64, nn)
    natmeqs = zeros(Float64, nn)
    natwgh = zeros(Float64, nn)
    natie = zeros(Float64, nn)      # Overall IE
    natiepc = zeros(Float64, nn)    # IE per capita
    natde = zeros(Float64, nn)      # Overall DE
    natdepc = zeros(Float64, nn)    # DE per capita
    natcf = zeros(Float64, nn)      # Overall CF
    natcfpc = zeros(Float64, nn)    # CF per capita
    natcfph = zeros(Float64, nn)    # CF per household
    natcfpeqs = zeros(Float64, nn)  # CF per equivalent size
    natcfpmeqs = zeros(Float64, nn) # CF per modified equivalent size

    for i=1:nn
        n = natList[i]
        for j = 1:length(hhsList[year][n])
            h = hhsList[year][n][j]
            ie = sum(indirectCE[year][n][:,j])
            de = sum(directCE[year][n][:,j])
            cf = sum(integratedCF[year][n][:,j])
            natsam[i] += siz[year][h]
            nateqs[i] += eqs[year][h]
            natmeqs[i] += meqs[year][h]
            natwgh[i] += w[year][h] * siz[year][h]
            natcf[i] += w[year][h] * cf
            natcfph[i] += cf
            natcfpeqs[i] += cf
            natcfpmeqs[i] += cf
            natie[i] += w[year][h] * ie
            natde[i] += w[year][h] * de
        end
        natcfpc[i] = natcf[i] / natwgh[i]
        natcfph[i] /= length(hhsList[year][n])
        natcfpeqs[i] /= nateqs[i]
        natcfpmeqs[i] /= natmeqs[i]
        natiepc[i] = natie[i] / natwgh[i]
        natdepc[i] = natde[i] / natwgh[i]
    end

    f = open(outputFile, "w")
    println(f, "Nation\tHHs\tMMs\tPop\tWeights\tCF_overall\tCF_percapita\tCF_perhh\tCF_pereqs\tCF_permeqs\tIE_overall\tIE_percapita\tDE_overall\tDE_percapita")
    for i=1:nn
        print(f, natList[i],"\t",length(hhsList[year][natList[i]]),"\t",natsam[i],"\t", pop[year][natList[i]],"\t",natwgh[i])
        print(f, "\t",natcf[i],"\t",natcfpc[i],"\t",natcfph[i],"\t",natcfpeqs[i],"\t",natcfpmeqs[i])
        print(f, "\t", natie[i], "\t", natiepc[i], "\t", natde[i], "\t", natdepc[i])
        println(f)
    end
    close(f)
end

function readCategoryData(inputFile, year, ntlv=0; subCategory="", except=[], nuts3pop=false)

    global sec, deSec, cat, gid, nam, pop, poplb, gidData, misDist
    global natList, natName, natA3, nuts, nutslList, pop, popList, catList
    global popcd, pophbscd, hbscd, gisCoord, popgiscd, popcdlist, ntcdlist
    global giscdlist, gispopcdlist, giscatlab, hbspopcdlist, hbscdlist
    global nutsLv = ntlv
    nuts_yr = Dict(2010=>"2010", 2015=>"2010")

    xf = XLSX.readxlsx(inputFile)
    if isa(year, Number); year = [year] end
    cdlen = ntlv+2

    sh = xf["Nation"]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r)>1 && !ismissing(r[1]) && !(string(r[1]) in natList)
            push!(natList, string(r[1]))
            natName[natList[end]] = string(r[2])
            natA3[natList[end]] = string(r[3])
        end
    end

    for y in year
        nuts[y], nutsList[y] = Dict{String, String}(), Dict{String, Array{String, 1}}()
        pop[y], poplb[y], popList[y] = Dict{String, Float64}(), Dict{String, String}(), Dict{String, Dict{String, Float64}}()
        for n in natList; nutsList[y][n], popList[y][n] = Array{String, 1}(), Dict{String, Float64}() end
        pophbscd[y], popcd[y], hbscd[y] = Dict{String, String}(), Dict{String, String}(), Dict{String, String}()
        popgiscd[y], popcdlist[y] = Dict{String, String}(), Array{String, 1}()
        ntcdlist[y], giscdlist[y], hbscdlist[y] = Array{String, 1}(), Array{String, 1}(), Array{String, 1}()
        gispopcdlist[y], hbspopcdlist[y] = Dict{String, Array{String, 1}}(), Dict{String, Array{String, 1}}()
        gisCoord[y] = Dict{String, Tuple{Float64, Float64}}()

        cat[y] = Dict{String, String}()
        tb = xf["Sector_" * string(y)][:]
        ci = findfirst(x -> x == (length(subCategory) == 0 ? "Category" : subCategory * "_category"), tb[1,:])
        for i in filter(x -> !ismissing(tb[x,1]), collect(2:size(tb,1)))
            s = string(tb[i,1])
            push!(sec, s)
            secName[s] = string(tb[i,2])
            if !ismissing(tb[i,ci]) && !(string(tb[i,ci]) in except); cat[y][s] = string(tb[i,ci]) end
        end
        # for r in XLSX.eachrow(sh)
        #     if XLSX.row_number(r)>1 && !ismissing(r[1])
        #         secCode = string(r[1])  # sector code
        #         push!(sec, secCode)
        #         if length(subCategory)==0 && !ismissing(r[4]) && !(string(r[4]) in except); cat[y][secCode] = string(r[4])
        #         elseif subCategory=="Food" && !ismissing(r[5]); cat[y][secCode] = string(r[5])
        #         elseif subCategory=="Electricity" && !ismissing(r[6]); cat[y][secCode] = string(r[6])
        #         elseif subCategory=="Transport" && !ismissing(r[7]); cat[y][secCode] = string(r[7])
        #         end
        #         secName[secCode]=string(r[2])
        #     end
        # end

        tb = xf["NUTS" * nuts_yr[y]][:]
        for i in filter(x -> !ismissing(tb[x,1]), collect(2:size(tb,1)))
            lv, ntcd, n = parse(Int, string(tb[i,3])), string(tb[i,1]), string(tb[i,4])
            if length(ntcd)==lv+2; nuts[y][ntcd] = ('/' in tb[i,2] ? strip(split(tb[i,2], '/')[1]) : strip(tb[i,2]))
            else println("NUTS level error: ", ntcd, "\t", lv)
            end
            if n in natList && lv == ntlv && ntcd[end] != 'Z'; push!(nutsList[y][n], ntcd) end
        end
        # sh = xf[nuts_sheet[y]]
        # for r in XLSX.eachrow(sh)
        #     if XLSX.row_number(r)>1
        #         lv, ntcd, n = parse(Int, string(r[3])), string(r[1]), string(r[4])
        #         if length(ntcd)==lv+2
        #             if '/' in r[2]; nuts[y][ntcd] = strip(split(r[2], '/')[1]) else nuts[y][ntcd] = strip(r[2]) end
        #         else println("NUTS level error: ", ntcd, "\t", lv)
        #         end
        #         if n in natList && lv==ntlv && ntcd[end] != 'Z'; push!(nutsList[y][n], ntcd) end
        #     end
        # end

        # add aggregated regions for region-code absent nations
        # agg_reg = ["FR"]
        agg_reg = natList
        for n in filter(x -> x in natList, agg_reg)
            rg = n
            for i=1:ntlv; rg *= "0" end
            if !(rg in nutsList[y][n])
                push!(nutsList[y][n], rg)
                nuts[y][rg] = nuts[y][n]
            end
        end

        tb = xf["Conc" * nuts_yr[y]][:]
        for i in filter(x -> !ismissing(tb[x,1]), collect(2:size(tb,1)))
            popcd[y][string(tb[i,1])], hbscd[y][string(tb[i,1])]  = string(tb[i,2]), string(tb[i,3])
            if !(string(tb[i,3]) in ntcdlist[y]); push!(ntcdlist[y], string(tb[i,3])) end
        end
        # sh = xf["Conc" * nuts_yr[y]]
        # for r in XLSX.eachrow(sh)
        #     if XLSX.row_number(r)>1
        #         popcd[y][string(r[1])], hbscd[y][string(r[1])]  = string(r[2]), string(r[3])
        #         if !(string(r[3]) in ntcdlist[y]); push!(ntcdlist[y], string(r[3])) end
        #     end
        # end

        tb = xf["PopCd" * nuts_yr[y]][:]
        hbs_i, gis_i = [findfirst(x -> x == lb * string(y), tb[1,:]) for lb in ["NUTS_", "GIS_"]]
        for i in filter(x -> !ismissing(tb[x,2]), collect(2:size(tb,1)))
            nt_map, nt_hbs, nt_gis = string(tb[i,1]), string(tb[i,hbs_i]), string(tb[i,gis_i])
            pophbscd[y][nt_map], popgiscd[y][nt_map] = nt_hbs, nt_gis
            if !(nt_map in popcdlist[y]); push!(popcdlist[y], nt_map) end
            if !haskey(hbspopcdlist[y], nt_hbs); hbspopcdlist[y][nt_hbs] = Array{String, 1}() end
            push!(hbspopcdlist[y][nt_hbs], nt_map)
            if !haskey(gispopcdlist[y], nt_gis); gispopcdlist[y][nt_gis] = Array{String, 1}() end
            push!(gispopcdlist[y][nt_gis], nt_map)
        end
        # sh = xf["PopCd" * nuts_yr[y]]
        # for r in XLSX.eachrow(sh)
        #     if XLSX.row_number(r)>1 && !ismissing(r[2])
        #         nt_map, nt_hbs, nt_gis = string(r[1]), string(r[2]), string(r[3])
        #         pophbscd[y][nt_map], popgiscd[y][nt_map] = nt_hbs, nt_gis
        #         if !(nt_map in popcdlist[y]); push!(popcdlist[y], nt_map) end
        #         if !haskey(hbspopcdlist[y], nt_hbs); hbspopcdlist[y][nt_hbs] = Array{String, 1}() end
        #         push!(hbspopcdlist[y][nt_hbs], nt_map)
        #         if !haskey(gispopcdlist[y], nt_gis); gispopcdlist[y][nt_gis] = Array{String, 1}() end
        #         push!(gispopcdlist[y][nt_gis], nt_map)
        #     end
        # end
        hbscdlist[y] = sort(filter(x->length(x)==cdlen, collect(keys(hbspopcdlist[y]))))
        giscdlist[y] = sort(filter(x->length(x)==cdlen, collect(keys(gispopcdlist[y]))))
    end

    if nuts3pop
        sntpop = Dict{Int, Dict{String, Float64}}(); for y in year; sntpop[y] = Dict{String, Float64}() end
        nant = Dict{Int, Array{String, 1}}(); for y in year; nant[y] = Array{String, 1}() end

        tb = xf["Pop_mod_NUTS3"][:]
        yi = [findfirst(x -> x == y, tb[1,:]) for y in year]
        for ri = 2:size(tb,1)
            ntcd, ntlb = strip.(split(replace(tb[ri, 1], ",Total,Total,Number"=>""), " - "))    # code, label
            if '/' in ntlb; ntlb = strip(split(ntlb, '/')[1]) end
            n = ntcd[1:2]
            if n in natList
                for i=1:length(year)
                    y = year[i]
                    s = strip(replace(string(tb[ri, yi[i]]), ['b','d','e','p',',']=>""))
                    if tryparse(Float64, s) != nothing && !haskey(pop[y], ntcd); pop[y][ntcd] = parse(Float64, s) end
                    if haskey(pophbscd[y], ntcd) && pophbscd[y][ntcd] in nutsList[y][n]
                        ncd = pophbscd[y][ntcd]
                        if !haskey(popList[y][n], ncd); popList[y][n][ncd] = 0; sntpop[y][ncd] = 0 end
                        if tryparse(Float64, s) != nothing
                            if length(ntcd) == cdlen; popList[y][n][ncd] += pop[y][ntcd]
                            elseif length(ntcd) == cdlen+1; sntpop[y][ncd] += pop[y][ntcd]
                            end
                        end
                    end
                    if !haskey(poplb[y], ntcd); poplb[y][ntcd] = ntlb end
                end
            end
        end
        for y in year, n in natList, nt in collect(keys(popList[y][n]))
            if popList[y][n][nt] == 0; popList[y][n][nt] = sntpop[y][nt] end
        end
        # sh = xf["Pop_mod_NUTS3"]
        # yridx = []
        # sntpop = Dict{Int, Dict{String, Float64}}(); for y in year; sntpop[y] = Dict{String, Float64}() end
        # nant = Dict{Int, Array{String, 1}}(); for y in year; nant[y] = Array{String, 1}() end
        # for r in XLSX.eachrow(sh)
        #     if XLSX.row_number(r)==1
        #         str = [string(r[i]) for i=1:21]
        #         yridx = [findfirst(x->x==string(y), str) for y in year]
        #     else
        #         ntlb = split(replace(r[1],",Total,Total,Number"=>""), " - ")
        #         ntcd, ntlb = ntlb[1], ntlb[2]   # code, label
        #         if '/' in ntlb; ntlb = strip(split(ntlb, '/')[1]) end
        #         n = ntcd[1:2]
        #         if n in natList
        #             for i=1:length(year)
        #                 y = year[i]
        #                 s = strip(replace(string(r[yridx[i]]), ['b','d','e','p',',']=>""))
        #                 if haskey(pophbscd[y], ntcd)
        #                     ncd = pophbscd[y][ntcd]
        #                     if ncd in nutsList[y][n]
        #                         if !haskey(popList[y][n], ncd); popList[y][n][ncd] = 0; sntpop[y][ncd] = 0 end
        #                         if tryparse(Float64, s)!==nothing
        #                             if !haskey(pop[y], ntcd); pop[y][ntcd] = parse(Float64, s) end
        #                             if length(ntcd) == cdlen; popList[y][n][ncd] += pop[y][ntcd]
        #                             elseif length(ntcd) == cdlen+1; sntpop[y][ncd] += pop[y][ntcd]
        #                             end
        #                         end
        #                     end
        #                 end
        #                 if !haskey(poplb[y], ntcd); poplb[y][ntcd] = ntlb end
        #             end
        #         end
        #     end
        # end
        # for y in year, n in natList, nt in collect(keys(popList[y][n]))
        #     if popList[y][n][nt] == 0; popList[y][n][nt] = sntpop[y][nt] end
        # end
    else
        tb = xf["Pop_mod"][:]
        yi = [findfirst(x -> x == y, tb[1,:]) for y in year]
        for ri = 2:size(tb,1)
            ntcd = string(rsplit(tb[ri, 1], ',', limit=2)[2])
            n = ntcd[1:2]
            for i = 1:length(year)
                if n in natList && haskey(pophbscd[year[i]], ntcd) && pophbscd[year[i]][ntcd] in nutsList[year[i]][n]
                    ncd = pophbscd[year[i]][ntcd]
                    s = strip(replace(string(tb[ri, yi[i]]), ['b','d','e','p']=>""))
                    if tryparse(Float64, s) != nothing
                        pop[year[i]][ntcd] = parse(Float64, s)
                        if !haskey(popList[year[i]][n], ncd); popList[year[i]][n][ncd] = 0 end
                        popList[year[i]][n][ncd] += pop[year[i]][ntcd]
                    end
                end
            end
        end
        # sh = xf["Pop_mod"]
        # yridx = []
        # for r in XLSX.eachrow(sh)
        #     if XLSX.row_number(r)==1
        #         str = [string(r[i]) for i=1:13]
        #         yridx = [findfirst(x->x==string(y), str) for y in year]
        #     else
        #         ntcd = string(split(r[1], ',')[end])
        #         n = ntcd[1:2]
        #         if n in natList
        #             for i=1:length(year)
        #                 if haskey(pophbscd[year[i]], ntcd)
        #                     ncd = pophbscd[year[i]][ntcd]
        #                     if ncd in nutsList[year[i]][n]
        #                         s = strip(replace(string(r[yridx[i]]), ['b','d','e','p']=>""))
        #                         if tryparse(Float64, s)!==nothing
        #                             pop[year[i]][ntcd] = parse(Float64, s)
        #                             if !haskey(popList[year[i]][n], ncd); popList[year[i]][n][ncd] = 0 end
        #                             popList[year[i]][n][ncd] += pop[year[i]][ntcd]
        #                         end
        #                     end
        #                 end
        #             end
        #         end
        #     end
        # end
    end

    sh = xf["GIS_coor"]
    for y in year, r in XLSX.eachrow(sh); if XLSX.row_number(r)>1; gisCoord[y][r[1]] = (r[2], r[3]) end end

    sh = xf["GIS_cat"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1; giscatlab[string(r[1])] = string(r[2]) end end
    close(xf)

    # sh = xf["DE_sector"]
    # for y in year
    #     deSec[y] = Array{String, 1}()
    #     deCodes[y] = Dict{String, Array{String, 1}}()
    #     deHbsList[y] = Array{String, 1}()
    #     deCat[y] = Dict{String, String}()
    #     deCatList[y] = Array{String, 1}()
    # end
    # for r in XLSX.eachrow(sh)
    #     if XLSX.row_number(r)>1 && !ismissing(r[1])
    #         if !(r[2] in deSec[r[1]])
    #             push!(deSec[r[1]], r[2])
    #             deCodes[r[1]][r[2]] = Array{String, 1}()
    #         end
    #         push!(deCodes[r[1]][r[2]], r[3])
    #         if !(r[3] in deHbsList[r[1]]); push!(deHbsList[r[1]], r[3]) end
    #         deCat[r[1]][r[2]] = r[4]
    #     end
    # end
    # for y in year
    #     sort!(deHbsList[y])
    #     deCatList[y] = sort(unique(collect(values(deCat[y]))))
    #     if length(subCategory)==0; push!(deCatList[y], "Total") else push!(deCatList[y], subCategory) end
    # end

    for y in year
        if length(catList) == 0; catList = sort(unique(values(cat[y])))
        elseif sort(catList) != sort(unique(values(cat[y]))); println("Categories are different in ", year)
        end
    end
    if length(subCategory)==0 && !("Total" in catList); push!(catList, "Total")
    elseif !(subCategory in catList); push!(catList, subCategory)
    end
end

function setCategory(list::Array{String,1})
    global catList
    if sort(catList)==sort(list); catList = list
    else println("Category items are different.\n",sort(catList),"\n",sort(list))
    end
end

function readHouseholdData(inputFile; period="monthly", sampleCheck=false, remove=false)  # period: "monthly"(default), "daily", or "annual"

    global yrList, hhsList, natList, regList, relList, disSta, nutsLv
    global pop, nat, reg, siz, eqs, meqs, typ, inc, exp, rel, wgh, pds

    year = 0
    nation = ""
    f = open(inputFile)
    readline(f)
    for l in eachline(f)
        s = string.(split(l, ','))
        year = parse(Int, s[1])
        if !(year in yrList)
            push!(yrList, year)
            hhsList[year], nat[year], reg[year]  = Dict{String, Array{String, 1}}(), Dict{String, String}(), Dict{String, String}()
            typ[year], siz[year]  = Dict{String, String}(), Dict{String, Int}()
            eqs[year], meqs[year]  = Dict{String, Float64}(), Dict{String, Float64}()
            inc[year], exp[year] = Dict{String, Float64}(), Dict{String, Float64}()
            pds[year], rel[year], wgh[year] = Dict{String, Int}(), Dict{String, Int}(), Dict{String, Float64}()
        end
        if !haskey(hhsList[year], s[2]); hhsList[year][s[2]] = Array{String, 1}() end
        hh = s[2]*"_"*s[3]      # replaced from "hh = s[3]"
        push!(hhsList[year][s[2]], hh)
        siz[year][hh] = parse(Int,s[5])
        eqs[year][hh], meqs[year][hh] = parse(Float64,s[12]), parse(Float64,s[13])
        nat[year][hh], reg[year][hh] = s[2], hbscd[year][s[4]][1:(nutsLv+2)]
        wgh[year][hh], pds[year][hh] = parse(Float64,s[6]), parse(Int,s[11])
        inc[year][hh], exp[year][hh] = parse(Float64,s[7]), parse(Float64,s[9])
    end
    filter!(x->x in unique(collect(values(nat[year]))), natList)

    tmp_reg = sort(unique(collect(values(reg[year]))))
    if remove; for y in yrList, n in natList; nutsList[y][n] = filter(x->x in tmp_reg, nutsList[y][n]) end end
    for y in yrList, n in natList
        nt = n; for i=1:nutsLv; nt *= "0" end
        if nt in nutsList[y][n] && !haskey(popList[y][n], nt)
            popList[y][n][nt] = sum([popList[y][n][x] for x in filter(x -> !(x in tmp_reg), collect(keys(popList[y][n])))])
        end
    end

    # convert household's income and expenditure data period
    if period=="daily"; cvr = 1/30 elseif period=="annual" cvr = 365/30 end
    for n in collect(keys(hhsList[year]))
        for h in hhsList[year][n]; inc[year][h] = inc[year][h] * cvr; exp[year][h] = exp[year][h] * cvr end
    end
    # if period=="daily"
    #     for n in collect(keys(hhsList[year]))
    #         for h in hhsList[year][n]; inc[year][h] = inc[year][h]/30; exp[year][h] = exp[year][h]/30 end
    #     end
    # elseif period=="annual"
    #     mmtoyy=365/30
    #     for n in collect(keys(hhsList[year]))
    #         for h in hhsList[year][n]; inc[year][h] = inc[year][h]*mmtoyy; exp[year][h] = exp[year][h]*mmtoyy end
    #     end
    # end
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

function readEmissionData(year, nations, inputFiles; mode = "ie")
    # mode: [ie] indirect emission, [de] direct emission

    global sec, hhsList, indirectCE, directCE
    if lowercase(mode) == "ie"; indirectCE[year] = Dict{String, Array{Float64, 2}}()
    elseif lowercase(mode) == "de"; directCE[year] = Dict{String, Array{Float64, 2}}()
    end

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
                e[findfirst(x->x==string(l[1]), sec),:] = map(x->parse(Float64,x), l[2:end])
            end
            if lowercase(mode) == "ie"; indirectCE[year][n] = e
            elseif lowercase(mode) == "de"; directCE[year][n] = e
            end
            close(f)
        end
    else println("Sizes of nation list (", nn,") and emission files (", length(inputFiles),") do not match")
    end
end

# function readIndirectEmission(year, nations, inputFiles)
#
#     global sec, hhsList
#     global indirectCE[year] = Dict{String, Array{Float64, 2}}()
#
#     ns = length(sec)
#     nn = length(nations)
#     if length(inputFiles) == nn
#         for i = 1:nn
#             n = nations[i]
#             nh = length(hhsList[year][n])
#             f = open(inputFiles[i])
#             readline(f)
#             e = zeros(Float64, ns, nh)
#             for l in eachline(f)
#                 l = split(l, '\t')
#                 e[findfirst(x->x==string(l[1]), sec),:] = map(x->parse(Float64,x), l[2:end])
#             end
#             indirectCE[year][n] = e
#             close(f)
#         end
#     else println("Sizes of nation list (", nn,") and emission files (", length(inputFiles),") do not match")
#     end
# end
#
# function readDirectEmission(year, nations, inputFiles)
#
#     global deSec, hhsList, deHbsList
#     global directCE[year] = Dict{String, Array{Float64, 2}}()
#
#     ns = length(deHbsList[year])
#     nn = length(nations)
#     if length(inputFiles) == nn
#         for i = 1:nn
#             n = nations[i]
#             nh = length(hhsList[year][n])
#             f = open(inputFiles[i]); readline(f)
#             de = zeros(Float64, ns, nh)
#             for l in eachline(f)
#                 l = split(l, '\t')
#                 j = findfirst(x->x==string(l[1]), deHbsList[year])
#                 l = l[2:end]
#                 de[j,:] = map(x->parse(Float64,x), l)
#             end
#             directCE[year][n] = de
#             close(f)
#         end
#     else println("Sizes of nation list (", nn,") and emission files (", length(inputFiles),") do not match")
#     end
#
#     # ns = length(deSec[year])
#     # nn = length(nations)
#     # if length(inputFiles) == nn
#     #     for i = 1:nn
#     #         n = nations[i]
#     #         nh = length(hhsList[year][n])
#     #         f = open(inputFiles[i]); readline(f)
#     #         de = zeros(Float64, ns, nh)
#     #         for l in eachline(f)
#     #             l = split(l, '\t')
#     #             j = findfirst(x->x==string(l[1]), deSec[year])
#     #             l = l[2:end]
#     #             de[j,:] = map(x->parse(Float64,x),l)
#     #         end
#     #         directCE[year][n] = de
#     #         close(f)
#     #     end
#     # else println("Sizes of nation list (", nn,") and emission files (", length(inputFiles),") do not match")
#     # end
# end

function integrateCarbonFootprint(year; mode="cf")

    global natList, hhsList, sec, directCE, indirectCE
    global integratedCF[year] = Dict{String, Array{Float64, 2}}()

    ie, de = indirectCE[year], directCE[year]
    # nds = length(deHbsList[year])
    # dsidx = [findfirst(x->x==deHbsList[year][i], sec) for i=1:nds]

    for n in filter(x -> haskey(hhsList[year], x), natList)
        nh = length(hhsList[year][n])
        if mode == "cf"; integratedCF[year][n] = ie[n] + de[n]
        elseif mode == "ie"; integratedCF[year][n] = ie[n]
        elseif mode == "de"; integratedCF[year][n] = de[n]
        end
    end
end

function readExpenditure(year, nations, inputFiles)

    # global sec, hhsList, indirectCE
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
    #         global indirectCE[year][n] = e
    #         close(f)
    #     end
    # else  println("Sizes of nation list (", nn,") and emission files (", length(inputFiles),") do not match")
    # end
end

function calculateNutsPopulationWeight(;year = 0)

    global yrList, natList, hhsList, nutsList, popList, wghNuts
    global siz, reg, pds, ntpop, ntsmp, ntwgh
    if year > 0; yrs = [year] else yrs = yrList end

    for y in yrs
        ntpop[y] = Dict{String, Dict{String, Array{Int, 1}}}()
        ntsmp[y] = Dict{String, Dict{String, Array{Int, 1}}}()
        ntwgh[y] = Dict{String, Dict{String, Array{Float64, 1}}}()
    end

    # count region samples
    typeidx = [2,3,4,0,0,0,0,0,5]
    for y in yrs
        ntsmp[y] = Dict{String, Dict{String, Array{Int,1}}}()
        for n in filter(x -> haskey(hhsList[y], x), natList)
            ntsmp[y][n] = Dict{String, Array{Int,1}}()
            for nt in nutsList[y][n]; ntsmp[y][n][nt] = zeros(Int, 5) end
            for hh in hhsList[y][n]
                ntsmp[y][n][reg[y][hh]][typeidx[pds[y][hh]]] += siz[y][hh]
                ntsmp[y][n][reg[y][hh]][1] += siz[y][hh]
            end
        end
    end

    # calculate weights
    for y in yrs
        ntwgh[y] = Dict{String, Dict{String, Array{Int,1}}}()
        for n in filter(x -> haskey(hhsList[y], x), natList)
            ntwgh[y][n] = Dict{String, Array{Int,1}}()
            for nt in filter(x -> x in ntcdlist[y], nutsList[y][n])
                ntwgh[y][n][nt] = zeros(Float64, 5)
                ntwgh[y][n][nt][1] = popList[y][n][nt] / ntsmp[y][n][nt][1]
            end
        end
    end

    # allocate household weight
    for y in yrs
        wghNuts[y] = Dict{String, Float64}()
        for n in filter(x -> haskey(hhsList[y], x), natList), hh in hhsList[y][n]; wghNuts[y][hh] = ntwgh[y][n][reg[y][hh]][1] end
    end

    return wghNuts
end

function categorizeHouseholdEmission(year; mode="cf", output="", hhsinfo=false, nutsLv=1)
    # mode: [ie] indirect CE, [de] direct CE, [cf] integrated CF
    global wgh, sec, hhid, cat, siz, inc, catList, natList, deHbsList
    global indirectCE, ieHHs, directCE, deHHs, integratedCF, cfHHs

    nc = length(catList)
    sl = Dict{Int, Array{String, 1}}()
    if isa(year, Number); year = [year] end
    for y in year; sl[y] = sec end

    if mode == "ie"; et = indirectCE; ht = ieHHs
    elseif mode == "de"; et = directCE; ht = deHHs
    elseif mode == "cf"; et = integratedCF; ht = cfHHs
    else println("wrong emission categorizing mode: ", mode)
    end

    # categorize emission data
    for y in year
        ns = length(sl[y])
        ht[y] = Dict{String, Array{Float64, 2}}()
        catidx = [if haskey(cat[y], s); findfirst(x->x==cat[y][s], catList) end for s in sl[y]] # make an index dict
        for n in natList
            nh = length(hhsList[y][n])
            e = et[y][n]
            ec = zeros(Float64, nh, nc)
            for i=1:nh, j=1:ns; if catidx[j]!=nothing; ec[i,catidx[j]] += e[j,i] end end    # categorizing
            for i=1:nc-1; ec[:, nc] += ec[:,i] end      # summing
            ht[y][n] = ec       # save the results
        end
    end

    # print the results
    if length(output)>0
        f = open(output, "w")
        print(f,"Year,Nation,HHID"); for c in catList; print(f, ",", c) end
        if hhsinfo; print(f, ",HH_size,MPCE,PopWgh") end; println(f)
        for y in year, n in natList, i = 1:length(hhsList[y][n])
            hh = hhsList[y][n][i]
            print(f, y,',',n,',',hh)
            for j = 1:length(catList); print(f, ",", ht[y][n][i,j]) end
            if hhsinfo; print(f, ",",siz[y][hh],",",inc[y][hh],",",wgh[y][hh]) end
            println(f)
        end
        close(f)
    end
end

function categorizeRegionalEmission(year=[]; mode = "cf", nutsLv=1, period="monthly", religion=false, popWgh=false, ntweigh=false)
    # mode: "ie" indirect CE, "de" direct CE, "cf" integrated CF
    # period: "monthly", "daily", or "annual"
    # religion: [true] categorize districts' features by religions

    global hhsList, natList, cat, reg, siz, inc, sam, ave, rel, pop, wgh, wghNuts, catList, nutsList, relList, pophbscd
    global ieHHs, ieReg, ieRegDiff, deHHs, deReg, deRegDiff, cfHH, cfReg, cfRegDiff
    global deCat, deCatList

    if isa(year, Number); year = [year] end
    if ntweigh; pwgh = wghNuts else pwgh = wgh end
    nc = length(catList)
    nr = length(relList)

    for y in year
        er, erd = Dict{String, Array{Float64, 2}}(), Dict{String, Array{Float64, 2}}()
        sam[y], ave[y] = Dict{String, Tuple{Int,Int}}(), Dict{String, Float64}()

        for n in natList
            hhs, nts = hhsList[y][n], nutsList[y][n]
            nh,nn = length(hhs), length(nts)

            # make index arrays
            ntidx = [findfirst(x->x==reg[y][hhs[i]], nts) for i=1:nh]
            if religion; relidx = [findfirst(x->x==rel[y][hhs[i]], relList) for i=1:nh] end

            # sum sample households and members by regions
            thbd = zeros(Float64, nn)   # total households by region
            tpbd = zeros(Float64, nn)   # total members of households by region
            for i=1:nh; thbd[ntidx[i]] += 1 end
            for i=1:nh; tpbd[ntidx[i]] += siz[y][hhs[i]] end
            for i=1:nn; sam[y][nts[i]] = (tpbd[i], thbd[i]) end

            # sum sample households and members by regions and by religions
            if religion
                thbdr = zeros(Float64, nn, nr)  # total households by district, by religion
                tpbdr = zeros(Float64, nn, nr)  # total members of households by district, by religion
                for i=1:nh
                    thbdr[ntidx[i],relidx[i]] += 1
                    tpbdr[ntidx[i],relidx[i]] += siz[y][hhs[i]]
                end
            end

            # calculate average monthly expenditure per capita by region
            totexp = zeros(Float64, nn)     # total expenditures by region
            if popWgh
                for i=1:nh; totexp[ntidx[i]] += inc[y][hhs[i]]*siz[y][hhs[i]]*pwgh[y][hhs[i]] end
                for i=1:nn; if nts[i] in ntcdlist[y]; ave[y][nts[i]] = totexp[i]/popList[y][n][nts[i]] end end
            else
                for i=1:nh; totexp[ntidx[i]] += inc[y][hhs[i]]*siz[y][hhs[i]] end
                for i=1:nn; ave[y][nts[i]] = totexp[i]/tpbd[i] end
            end
            # convert 'AVEpC' to annual or daily
            if period=="annual"; mmtoyy = 365/30; for i=1:nn; ave[y][nts[i]] = ave[y][ntd[i]] * mmtoyy end
            elseif period=="daily"; for i=1:nn; if nts[i] in ntcdlist[y]; ave[y][nts[i]] = ave[y][nts[i]]/30 end end
            end

            # categorize emission data
            if mode == "ie"; ec = ieHHs[y][n]
            elseif mode == "de"; ec = deHHs[y][n]
            elseif mode == "cf"; ec = cfHHs[y][n]
            end
            en = zeros(Float64, nn, nc)
            if popWgh; for i=1:nh; en[ntidx[i],:] += ec[i,:]*pwgh[y][hhs[i]] end
            else for i=1:nh; en[ntidx[i],:] += ec[i,:] end
            end
            # normalizing
            if popWgh; for i=1:nn; if nts[i] in ntcdlist[y]; for j=1:nc; en[i,j] /= popList[y][n][nts[i]] end end end
            else for i=1:nc; en[:,i] ./= tpbd end
            end
            # calculate differences
            avg = mean(en, dims=1)
            ecn = zeros(size(en))
            for i=1:size(en,2); ecn[:,i] = (en[:,i].-avg[i])/avg[i] end
            # save the results
            er[n] = en
            erd[n] = ecn
        end

        if mode == "ie"; ieReg[y] = er; ieRegDiff[y] = erd
        elseif mode == "de"; deReg[y] = er; deRegDiff[y] = erd
        elseif mode == "cf"; cfReg[y] = er; cfRegDiff[y] = er
        else println("wrong emission categorizing mode: ", mode)
        end
    end

    # for y in year
    #     if cf
    #         ieReg[y] = Dict{String, Array{Float64, 2}}()
    #         ieRegDiff[y] = Dict{String, Array{Float64, 2}}()
    #     end
    #     if de
    #         ndc = length(deCatList[y])
    #         deReg[y] = Dict{String, Array{Float64, 2}}()
    #         deRegDiff[y] = Dict{String, Array{Float64, 2}}()
    #     end
    #
    #     sam[y] = Dict{String, Tuple{Int,Int}}()
    #     ave[y] = Dict{String, Float64}()
    #
    #     for n in natList
    #         hhs = hhsList[y][n]
    #         nts = nutsList[y][n]
    #         nh = length(hhs)
    #         nn = length(nts)
    #
    #         # make index arrays
    #         ntidx = [findfirst(x->x==reg[y][hhs[i]], nts) for i=1:nh]
    #         if religion; relidx = [findfirst(x->x==rel[y][hhs[i]], relList) for i=1:nh] end
    #
    #         # sum sample households and members by regions
    #         thbd = zeros(Float64, nn)   # total households by region
    #         tpbd = zeros(Float64, nn)   # total members of households by region
    #         for i=1:nh; thbd[ntidx[i]] += 1 end
    #         for i=1:nh; tpbd[ntidx[i]] += siz[y][hhs[i]] end
    #         for i=1:nn; sam[y][nts[i]] = (tpbd[i], thbd[i]) end
    #
    #         # sum sample households and members by regions and by religions
    #         if religion
    #             thbdr = zeros(Float64, nn, nr)  # total households by district, by religion
    #             tpbdr = zeros(Float64, nn, nr)  # total members of households by district, by religion
    #             for i=1:nh
    #                 thbdr[ntidx[i],relidx[i]] += 1
    #                 tpbdr[ntidx[i],relidx[i]] += siz[y][hhs[i]]
    #             end
    #         end
    #
    #         # calculate average monthly expenditure per capita by region
    #         totexp = zeros(Float64, nn)     # total expenditures by region
    #         if popWgh
    #             for i=1:nh; totexp[ntidx[i]] += inc[y][hhs[i]]*siz[y][hhs[i]]*pwgh[y][hhs[i]] end
    #             for i=1:nn; if nts[i] in ntcdlist[y]; ave[y][nts[i]] = totexp[i]/popList[y][n][nts[i]] end end
    #         else
    #             for i=1:nh; totexp[ntidx[i]] += inc[y][hhs[i]]*siz[y][hhs[i]] end
    #             for i=1:nn; ave[y][nts[i]] = totexp[i]/tpbd[i] end
    #         end
    #         # convert 'AVEpC' to annual or daily
    #         if period=="annual"; mmtoyy = 365/30; for i=1:nn; ave[y][nts[i]] = ave[y][ntd[i]] * mmtoyy end
    #         elseif period=="daily"; for i=1:nn; if nts[i] in ntcdlist[y]; ave[y][nts[i]] = ave[y][nts[i]]/30 end end
    #         end
    #
    #         # categorize emission data
    #         if cf
    #             ec = ieHHs[y][n]
    #             en = zeros(Float64, nn, nc)
    #             if popWgh; for i=1:nh; en[ntidx[i],:] += ec[i,:]*pwgh[y][hhs[i]] end
    #             else for i=1:nh; en[ntidx[i],:] += ec[i,:] end
    #             end
    #             # normalizing
    #             if popWgh; for i=1:nn; if nts[i] in ntcdlist[y]; for j=1:nc; en[i,j] /= popList[y][n][nts[i]] end end end
    #             else for i=1:nc; en[:,i] ./= tpbd end
    #             end
    #             # calculate differences
    #             avg = mean(en, dims=1)
    #             ecn = zeros(size(en))
    #             for i=1:size(en,2); ecn[:,i] = (en[:,i].-avg[i])/avg[i] end
    #             # save the results
    #             ieReg[y][n] = en
    #             ieRegDiff[y][n] = ecn
    #         end
    #         if de
    #             dec = deHHs[y][n]
    #             den = zeros(Float64, nn, ndc)
    #             if popWgh; for i=1:nh; den[ntidx[i],:] += dec[i,:]*pwgh[y][hhs[i]] end
    #             else for i=1:nh; den[ntidx[i],:] += dec[i,:] end
    #             end
    #             # normalizing
    #             if popWgh; for i=1:nn; if nts[i] in ntcdlist[y]; for j=1:ndc; den[i,j] /= popList[y][n][nts[i]] end end end
    #             else for i=1:ndc; den[:,i] ./= tpbd end
    #             end
    #             # calculate differences
    #             davg = mean(den, dims=1)
    #             decn = zeros(size(den))
    #             for i=1:size(den,2); decn[:,i] = (den[:,i].-davg[i])/davg[i] end
    #             # save the results
    #             deReg[y][n] = den
    #             deRegDiff[y][n] = decn
    #         end
    #     end
    # end
end

# function categorizeDirectEmission(years; output="", hhsinfo=false, nutsLv=1)
#     global wgh, deSec, hhid, deCat, deCatList, siz, inc, natList
#     global directCE, deHHs
#     if isa(years, Number); years = [years] end
#
#     # categorize direct emission data
#     for y in years
#         nc = length(deCatList[y])
#         ns = length(deSec[y])
#         catidx = [if haskey(deCat[y], s); findfirst(x->x==deCat[y][s], deCatList[y]) end for s in deSec[y]]     # make an index dict
#
#         deHHs[y] = Dict{String, Array{Float64, 2}}()
#         for n in natList
#             nh = length(hhsList[y][n])
#             de = directCE[y][n]
#             dec = zeros(Float64, nh, nc)
#             # categorizing
#             for i=1:nh; for j=1:ns; if catidx[j]!=nothing; dec[i,catidx[j]] += de[j,i] end end end
#             # summing
#             for i=1:nc-1; dec[:, nc] += dec[:,i] end
#             # save the results
#             deHHs[y][n] = dec
#         end
#     end
#
#     # print the results
#     if length(output)>0
#         f = open(output, "w")
#         print(f,"Year,Nation,HHID"); for c in deCatList[years[1]]; print(f, ",", c) end
#         if hhsinfo; print(f, ",HH_size,MPCE,PopWgh") end; println(f)
#         for y in years
#             for n in natList
#                 for i = 1:length(hhsList[y][n])
#                     hh = hhsList[y][n][i]
#                     print(f, y,',',n,',',hh)
#                     for j = 1:length(deCatList[y]); print(f, ",", deHHs[y][n][i,j]) end
#                     if hhsinfo; print(f, ",",siz[y][hh],",",inc[y][hh],",",wgh[y][hh]) end
#                     println(f)
#                 end
#             end
#         end
#         close(f)
#     end
# end

function printRegionalEmission(years, outputFile; mode=[],totm=false,expm=false,popm=false,relm=false,wghm=false,denm=false,povm=false,ntweigh=false)
    # mode: "ie" indirect CE, "de" direct CE, "cf" integrated CF
    # expm: print average expenditure per capita
    # popm: print population related figures
    # hhsm: print households related figures
    # relm: print relgion related figures
    # denm: print population density
    # povm: print poverty rates

    global natList, hhsList, catList, nutsList, relList, relName, popList, ave, ieReg, deReg, cfReg
    global wgh, wghNuts

    if ntweigh; pwgh = wghNuts else pwgh = wgh end
    f = open(outputFile, "w")

    nc = length(catList)
    print(f,"Year,Nation,NUTS")
    for m in mode; for c in catList; print(f, ","*uppercase(m)*"_", c) end end
    if totm; for m in mode; print(f, ",Overall_"*uppercase(m)) end end
    if expm; print(f, ",Income") end
    if popm; print(f, ",Population") end
    if wghm; print(f, ",PopWeight") end
    # if povm; print(f, ",PovertyRatio") end
    println(f)
    for y in years, n in natList
        nts = nutsList[y][n]
        if "cf" in mode; idxs = [x for x = 1:length(nts) if !isnan(cfReg[y][n][x,end]) && cfReg[y][n][x,end]!=0]
        elseif "ie" in mode; idxs = [x for x = 1:length(nts) if !isnan(ieReg[y][n][x,end]) && ieReg[y][n][x,end]!=0]
        elseif "de" in mode; idxs = [x for x = 1:length(nts) if !isnan(deReg[y][n][x,end]) && deReg[y][n][x,end]!=0]
        end
        for i in idxs
            print(f, y,',',n,',',nts[i])
            for m in mode
                if m == "cf"; for j = 1:nc; print(f, ",", cfReg[y][n][i,j]) end
                elseif m == "ie"; for j = 1:nc; print(f, ",", ieReg[y][n][i,j]) end
                elseif m == "de"; for j = 1:nc; print(f, ",", deReg[y][n][i,j]) end
                end
            end
            if haskey(popList[y][n], nts[i])
                if totm
                    for m in mode
                        if m == "cf"; print(f, ",", cfReg[y][n][i,end]*popList[y][n][nts[i]])
                        elseif m == "ie"; print(f, ",", ieReg[y][n][i,end]*popList[y][n][nts[i]])
                        elseif m == "de"; print(f, ",", deReg[y][n][i,end]*popList[y][n][nts[i]])
                        end
                    end
                end
                if expm; if nts[i] in ntcdlist[y]; print(f, ",", ave[y][nts[i]]) else print(f, ",NA") end end
                if popm; print(f, ",", popList[y][n][nts[i]]) end
                if wghm; print(f, ",", sum([pwgh[y][h]*siz[y][h] for h in filter(x->reg[y][x]==nts[i],hhsList[y][n])])) end
                # if povm; print(f, ",", disPov[regList[i]]) end
            end
            println(f)
        end
    end
    close(f)

    # for y in years, n in natList
    #     ndc = length(deCatList[y])
    #     nts = nutsList[y][n]
    #     if cf; en = ieReg[y][n] end
    #     if de; den = deReg[y][n] end
    #     if cf; idxs = [x for x = 1:length(nts) if !isnan(en[x,end]) && en[x,end]!=0]
    #     elseif de; idxs = [x for x = 1:length(nts) if !isnan(den[x,end]) && den[x,end]!=0]
    #     end
    #     for i in idxs
    #         print(f, y,',',n,',',nts[i])
    #         if cf; for j = 1:nc; print(f, ",", en[i,j]) end end
    #         if de; for j = 1:ndc; print(f, ",", den[i,j]) end end
    #         if haskey(popList[y][n], nts[i])
    #             if totm
    #                 if cf; print(f, ",", en[i,end]*popList[y][n][nts[i]]) end
    #                 if de; print(f, ",", den[i,end]*popList[y][n][nts[i]]) end
    #             end
    #             if expm; if nts[i] in ntcdlist[y]; print(f, ",", ave[y][nts[i]]) else print(f, ",NA") end end
    #             if popm; print(f, ",", popList[y][n][nts[i]]) end
    #             if wghm; print(f, ",", sum([pwgh[y][h]*siz[y][h] for h in filter(x->reg[y][x]==nts[i],hhsList[y][n])])) end
    #             # if povm; print(f, ",", disPov[regList[i]]) end
    #         end
    #         println(f)
    #     end
    # end
    # close(f)
end

function exportRegionalEmission(years=[], tag="", outputFile=""; mode="cf", nspan=128, minmax=[], descend=false,empty=false,logarithm=false,nutsmode = "gis")
    # mode: "ie" indirect CE, "de" direct CE, "cf" integrated CF
    # nutsmode = "gis": NUTS codes follow GIS-map's NUTS (ex. DE1, DE2, DE3, ..., EL1, EL2, ...)
    # nutsmode = "hbs": NUTS codes follow HBS's NUTS (ex. DE0, DE3, DE4, ..., EL0, ...)

    global catList, nuts, nutsLv, nutsList, natList, popList, sam, ave, reg
    global gisNutsList, hbscd, gispopcdlist, giscdlist, hbspopcdlist, hbscdlist
    global ieReg, gisRegionalIe, gisRegionalIeRank, gisRegionalIePerCap, gisRegionalIeRankPerCap
    global deReg, gisRegionalDe, gisRegionalDeRank, gisRegionalDePerCap, gisRegionalDeRankPerCap
    global cfReg, gisRegionalCF, gisRegionalCFrank, gisRegionalCFperCap, gisRegionalCFrankPerCap

    nc, cdlen = length(catList), nutsLv+2
    labels, labelspc = Dict{Int, Array{String,2}}(), Dict{Int, Array{String,2}}()

    for y in years
        nts, ntslist = Dict{String, Array{String, 1}}(), Array{String, 1}()
        if nutsmode == "gis"; ntkeys = giscdlist[y]
        elseif nutsmode == "hbs"; ntkeys = hbscdlist[y]
        else println("Error: wrong NUTS mode, ", nutsmode)
        end
        for n in natList
            nts[n] = filter(x->x in ntkeys, nutsList[y][n])
            append!(ntslist, nts[n])
        end
        nn = length(ntslist)

        # making exporting table
        tb = zeros(Float64, nn, nc)     # regional CF
        tbpc = zeros(Float64, nn, nc)   # regional CF per capita
        spo = zeros(Float64, nn)        # number of sample population by region
        tpo = zeros(Float64, nn)        # total number of population by region
        aec = zeros(Float64, nn)        # average expenditure per capita by region

        for n in natList
            if mode=="ie"; ec = ieReg[y][n]; elseif mode=="de"; ec = deReg[y][n]; elseif mode=="cf"; ec = cfReg[y][n] end
            for nt in nts[n]
                hnt = hbscd[y][nt]
                gidx = findfirst(x->x==nt, ntslist)
                if nutsmode == "gis"
                    pnts = filter(x->length(x)==cdlen+1, gispopcdlist[y][nt])
                    ntidx = findfirst(x->x==hbscd[y][nt], nutsList[y][n])
                    for pnt in pnts
                        if haskey(pop[y], pnt)
                            tb[gidx,:] += ec[ntidx,:] * pop[y][pnt]
                            tpo[gidx] += pop[y][pnt]
                            aec[gidx] += ave[y][hnt] * pop[y][pnt]
                        else
                            subpnts = filter(x->length(x)==cdlen+1, gispopcdlist[y][nt])
                            for spnt in subpnts
                                if haskey(pop[y], spnt)
                                    tb[gidx,:] += ec[ntidx,:] * pop[y][spnt]
                                    tpo[gidx] += pop[y][spnt]
                                    aec[gidx] += ave[y][hnt] * pop[y][spnt]
                                end
                            end
                        end
                    end
                elseif nutsmode == "hbs"
                    ntidx = findfirst(x->x==nt, nutsList[y][n])
                    tb[gidx,:] += ec[ntidx,:] * popList[y][n][nt]
                    tpo[gidx] += popList[y][n][nt]
                    aec[gidx] += ave[y][hnt] * popList[y][n][nt]
                end
            end
        end
        # normalizing
        for i=1:nn
            aec[i] /= tpo[i]
            for j=1:nc; tbpc[i,j] = tb[i,j] / tpo[i] end
        end

        modeTag = uppercase(mode)
        filename = replace(replace(replace(outputFile,"YEAR_"=>string(y)*"_"), ".txt"=>"_"*modeTag*".txt"), ".csv"=>"_"*modeTag*".csv")
        rank, labels[y] = exportRegionalTables(replace(filename,"_OvPcTag"=>"_overall"), tag, ntslist, nspan, minmax[1], tb, logarithm, descend, tab_mode=mode)
        rankpc, labelspc[y] = exportRegionalTables(replace(filename,"_OvPcTag"=>"_percap"), tag, ntslist, nspan, minmax[2], tbpc, logarithm, descend, tab_mode=mode)

        gisTotPop[y] = tpo
        if nutsmode == "hbs"; gisSamPop[y] = spo end
        gisAvgExp[y] = aec
        gisNutsList[y] = ntslist
        if mode == "ie"
            gisRegionalIe[y] = tb
            gisRegionalIePerCap[y] = tbpc
            gisRegionalIeRank[y] = rank
            gisRegionalIeRankPerCap[y] = rankpc
        elseif mode == "de"
            gisRegionalDe[y] = tb
            gisRegionalDePerCap[y] = tbpc
            gisRegionalDeRank[y] = rank
            gisRegionalDeRankPerCap[y] = rankpc
        elseif mode == "cf"
            gisRegionalCF[y] = tb
            gisRegionalCFperCap[y] = tbpc
            gisRegionalCFrank[y] = rank
            gisRegionalCFrankPerCap[y] = rankpc
        else println("wrong emission categorizing mode: ", mode)
        end
    end

    return labels, labelspc
end

function exportRegionalTables(outputFile, tag, ntslist, nspan, minmax, tb, logarithm, descend; tab_mode="cf")
    # This function is for [exportRegionalEmission]

    global yrList, catList
    nc, nn = length(catList), length(ntslist)

    # find min. and max.: overall CF
    if length(minmax)==1; maxde = [minmax[1][2] for i=1:nc]; minde = [minmax[1][1] for i=1:nc]
    elseif length(minmax)==nc; maxde = [minmax[i][2] for i=1:nc]; minde = [minmax[i][1] for i=1:nc]
    elseif logarithm; maxde = [log10(maximum(tb[:,i])) for i=1:nc]; minde = [log10(minimum(tb[:,i])) for i=1:nc]
    else maxde = [maximum(tb[:,i]) for i=1:nc]; minde = [minimum(tb[:,i]) for i=1:nc]
    end
    replace!(minde, Inf=>0, -Inf=>0)
    # grouping by ratios; ascending order: overall CF
    span = zeros(Float64, nspan+1, nc)
    over = [maxde[i] < maximum(tb[:,i]) for i=1:nc]
    for j=1:nc
        if over[j]; span[:,j] = [[(maxde[j]-minde[j])*(i-1)/(nspan-1)+minde[j] for i=1:nspan]; maximum(tb[:,j])]
        else span[:,j] = [(maxde[j]-minde[j])*(i-1)/nspan+minde[j] for i=1:nspan+1]
        end
    end
    if logarithm; for i=1:size(span,1), j=1:nc; span[i,j] = 10^span[i,j] end end
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
        for j=1:nc, i=1:nn; rank[i,j] = nspan - rank[i,j] + 1 end
    end
    # prepare labels
    labels = Array{String, 2}(undef,nspan,nc)
    for j=1:nc
        lbstr = [string(round(span[i,j],digits=0)) for i=1:nspan+1]
        if descend; labels[:,j] = [lbstr[i+1]*"-"*lbstr[i] for i=1:nspan]
        else labels[:,j] = [lbstr[i]*"-"*lbstr[i+1] for i=1:nspan]
        end
        if over[j]; if descend; labels[1,j] = "over "*lbstr[2] else labels[nspan,j] = "over "*lbstr[nspan] end end
    end
    # exporting table: overall CF
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
    # exporting group table: overall CF
    f = open(replace(outputFile, ".csv"=>"_gr.csv"), "w")
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

    return rank, labels
end

function exportEmissionDiffRate(years=[], tag="", outputFile="", maxr=0.5, minr=-0.5, nspan=128; descend=false, empty=false)

    global catList, nutsList, gisNutsList, gisRegionalCF, gisRegionalCFperCap
    global gisRegionalCFdiff, gisRegionalCFdiffRank, gisRegionalCFdiffPerCap, gisRegionalCFdiffRankPerCap

    nc = length(catList)
    spanval = Dict{Int, Array{Float64, 2}}()
    spanvalpc = Dict{Int, Array{Float64, 2}}()

    for y in years
        ntslist = gisNutsList[y]
        gre = gisRegionalCF[y]
        grepc = gisRegionalCFperCap[y]

        filename = replace(outputFile,"YEAR_"=>string(y)*"_")
        gred, rank, spanval[y] = exportEmissionDiffTable(replace(filename,"_OvPcTag"=>"_overall"), tag, ntslist, gre, maxr, minr, nspan, descend, empty)
        gredpc, rankpc, spanvalpc[y] = exportEmissionDiffTable(replace(filename,"_OvPcTag"=>"_percap"), tag, ntslist, grepc, maxr, minr, nspan, descend, empty)

        gisRegionalCFdiff[y] = gred
        gisRegionalCFdiffRank[y] = rank
        gisRegionalCFdiffPerCap[y] = gredpc
        gisRegionalCFdiffRankPerCap[y] = rankpc
    end

    return spanval, spanvalpc
end

function exportEmissionDiffTable(outputFile, tag, ntslist, gre, maxr, minr, nspan, descend, empty)
    # this function is for [exportEmissionDiffRate]

    global catList
    nc = length(catList)

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
    if descend; for j=1:size(gred,2), i=1:size(gred,1); rank[i,j] = nspan - rank[i,j] + 1 end end

    # exporting difference table
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

    return gred, rank, spanval
end

function findMajorCity(year, ntslist, nutsmode; modnuts=false)

    global majorCity, nuts, pop, poplb, gispopcdlist, hbspopcdlist

    majorCity[year] = Dict{String, String}()
    for nt in ntslist
        mjnt = ""; mjpop = 0
        if nutsmode == "gis"; pnts = filter(x->length(x)==5, gispopcdlist[year][nt])
        elseif nutsmode == "hbs"; pnts = filter(x->length(x)==5, hbspopcdlist[year][nt])
        end
        for pnt in pnts; if haskey(pop[year], pnt) && pop[year][pnt] > mjpop; mjpop = pop[year][pnt]; mjnt = pnt end end
        if mjnt==""
            if nutsmode == "gis"; pnts = filter(x->length(x)==4, gispopcdlist[year][nt])
            elseif nutsmode == "hbs"; pnts = filter(x->length(x)==4, hbspopcdlist[year][nt])
            end
            for pnt in pnts; if haskey(pop[year], pnt) && pop[year][pnt] > mjpop; mjpop = pop[year][pnt]; mjnt = pnt end end
        end
        majorCity[year][nt] = mjnt
        if modnuts && mjnt!="" && haskey(poplb[year], mjnt); nuts[year][nt] *= " (including "* poplb[year][mjnt] *")" end
    end
end

function exportWebsiteFiles(year, path; nutsmode = "gis", rank=false, empty=false, major=false)

    global natName, nuts, pop, popList, poplb, catList, gisNutsList, gisTotPop, gisAvgExp
    global pophbscd, hbscd, gispopcdlist, hbspopcdlist, majorCity, gisCoord
    global gisRegionalCF, gisRegionalCFrank, gisRegionalCFdiffPerCap, gisRegionalCFdiffRankPerCap
    if isa(year, Number); year = [year] end

    for y in year
        tbntlist = gisNutsList[y]
        if y == 2010; if nutsmode == "gis"; exceptNt = ["FR0","DE0","EL0"] elseif nutsmode == "hbs"; exceptNt = ["FR0"] end
        elseif y == 2015; exceptNt = []
        end
        ntslist = filter(x->!(x in exceptNt), tbntlist)

        gre, grer, gredpc, gredrpc = gisRegionalCF[y], gisRegionalCFrank[y], gisRegionalCFdiffPerCap[y], gisRegionalCFdiffRankPerCap[y]

        # find major city of NUTS1
        if major; findMajorCity(y, ntslist, nutsmode, modnuts = true) end

        # print center file
        mkpath(path*string(y))
        f = open(path*string(y)*"/centers.csv", "w")
        println(f, "\"NO\",\"GID2CODE\",\"PNAME\",\"DNAME\",\"x\",\"y\"")
        cnt = 1
        for nt in ntslist
            println(f,"\"",cnt,"\",\"",nt,"\",\"",natName[nt[1:2]],"\",\"",nuts[y][nt],"\",\"",gisCoord[y][nt][1],"\",\"",gisCoord[y][nt][2],"\"")
            cnt += 1
        end
        close(f)

        # print english file
        f = open(path*string(y)*"/english.txt", "w")
        println(f, "KEY_CODE\tEN_NAME")
        for nt in ntslist
            if nt[3] == '0' && nt[1:2] != "DE";
                if '(' in nuts[y][nt]; mjcity = '(' * split(nuts[y][nt], '(')[2] else mjcity = "" end
                println(f, nt, "\t", natName[nt[1:2]] * mjcity)
            else println(f, nt, "\t", nuts[y][nt],", ",natName[nt[1:2]])
            end
        end
        close(f)

        # print english_name file
        f = open(path*string(y)*"/english_match.txt", "w")
        println(f, "KEY_CODE\tSTATE_CODE\tSTATE\tDISTRICT")
        for nt in ntslist; println(f, nt, "\t", nt[1:2], "\t", natName[nt[1:2]], "\t", nuts[y][nt]) end
        close(f)

        # print ALLP file
        f = open(path*string(y)*"/ALLP.txt", "w")
        println(f, "ALL\tALLP")
        catidx = findfirst(x->giscatlab[x]=="All", catList)
        for nt in ntslist; println(f, nt, "\t", grer[findfirst(x->x==nt, tbntlist), catidx]) end
        close(f)

        # print CF files
        for j=1:length(catList)
            mkpath(path*string(y)*"/CFAV/")
            f = open(path*string(y)*"/CFAV/"*"CFAV_"*giscatlab[catList[j]]*".txt","w")
            print(f, "KEY_CODE\tSTATE\tDISTRICT\tSTATE_NAME\tDISTRICT_NAME\tGENHH_ALLP\tGENHH_APPPC")
            if catList[j]=="Total" || catList[j]=="All"; println(f, "\tANEXPPC\tPOP")
            else println(f)
            end
            for i=1:length(ntslist)
                nt = ntslist[i]
                tbidx = findfirst(x->x==nt, tbntlist)
                print(f, nt,"\t",nt[1:2],"\t",nt,"\t",natName[nt[1:2]],"\t",nuts[y][nt],"\t")
                printfmt(f, "{:f}", gre[tbidx,j]); print(f, "\t",gre[tbidx,j]/gisTotPop[y][tbidx])
                if catList[j]=="Total" || catList[j]=="All"; println(f,"\t",gisAvgExp[y][tbidx],"\t",convert(Int, gisTotPop[y][tbidx]))
                else println(f)
                end
            end
            if empty
                # for not-covered NUTS data
            end
            close(f)

            mkpath(path*string(y)*"/CFAC/")
            f = open(path*string(y)*"/CFAC/"*"CFAC_"*giscatlab[catList[j]]*".txt","w")
            println(f, "KEY_CODE\tSTATE\tDISTRICT\tSTATE_NAME\tDISTRICT_NAME\tGENHH_ALLPPC")
            for i=1:length(ntslist)
                nt = ntslist[i]
                tbidx = findfirst(x->x==nt, tbntlist)
                print(f, nt,"\t",nt[1:2],"\t",nt,"\t",natName[nt[1:2]],"\t",nuts[y][nt],"\t")
                if rank; println(f, gredrpc[tbidx,j]) else println(f, gredpc[tbidx,j]) end
            end
            if empty
                # for not-covered NUTS data
            end
            close(f)
        end
    end
end

function buildWebsiteFolder(years, centerpath, outputpath; nutsmode = "gis", rank = false)

    global natList, catList, nuts, natName, natA3, gisNutsList, giscatlab
    global gisRegionalIe, gisRegionalIeRank, gisRegionalCFdiffPerCap, gisRegionalCFdiffRankPerCap
    global gisTotPop, gisAvgExp

    cenfile = "centers.csv"
    engfile = "english_match.txt"
    allfile = "ALLP.txt"

    for y in years
        tbntlist = gisNutsList[y]
        if nutsmode == "gis"; ntslist = filter(x->!(x in ["FR0","DE0","EL0"]), tbntlist)
        elseif nutsmode == "hbs"; ntslist = filter(x->!(x in ["FR0"]), tbntlist)
        end

        gre = gisRegionalIe[y]
        grer = gisRegionalIeRank[y]
        gredpc = gisRegionalCFdiffPerCap[y]
        gredrpc = gisRegionalCFdiffRankPerCap[y]
        centers = Dict{String, Array{Array{String, 1}, 1}}()    # {nation, center data}

        # find major city of NUTS1
        if major; findMajorCity(ntslist, nutsmode, modnuts = true) end

        # read center data
        f = open(centerpath * "centers_" * string(y) * ".csv")
        titleLine = readline(f)
        for l in eachline(f)
            s = string.(split(l, ','))
            n = replace(s[2], "\""=>"")[1:2]
            if !haskey(centers, n); centers[n] = Array{Array{String, 1}, 1}() end
            push!(centers[n], s[2:end])
        end
        close(f)

        # build web-folders
        for n in natList
            nts = sort(filter(x->x[1:2]==n, ntslist))
            fp = outputpath*string(y)*"/data/"*natA3[n]*'/'
            mkpath(fp)

            # print INI file
            f = open(outputpath*string(y)*"/data/config."*natA3[n]*".ini", "w")
            println(f, "name = ", natName[n])
            println(f, "default_area = ", nts[1])
            println(f, "zoom = ", 7)
            println(f, "year = ", y)
            println(f, "kmz_url = http://data.spatialfootprint.com.s3-website-ap-northeast-1.amazonaws.com/", lowercase(natA3[n]))
            close(f)

            # print center files
            f = open(fp*cenfile, "w")
            println(f, titleLine)
            i = 1
            for ct in centers[n]
                print(f, "\"",string(i),"\"")
                for c in ct; print(f, ",",c) end
                println(f)
                i += 1
            end
            close(f)

            # print english_name file
            f = open(fp*engfile, "w")
            println(f, "KEY_CODE\tSTATE_CODE\tSTATE\tDISTRICT")
            for nt in nts; println(f, nt, "\t", n, "\t", natName[n], "\t", nuts[y][nt]) end
            close(f)

            # print ALLP file
            f = open(fp*allfile, "w")
            println(f, "ALL\tALLP")
            catidx = findfirst(x->giscatlab[x]=="All", catList)
            for nt in nts
                ntidx = findfirst(x->x==nt, tbntlist)
                println(f, nt, "\t", grer[ntidx,catidx])
            end
            close(f)

            # print CF files
            for j=1:length(catList)
                mkpath(fp*"CFAV/")
                f = open(fp*"CFAV/CFAV_"*giscatlab[catList[j]]*".txt","w")
                print(f, "KEY_CODE\tSTATE\tDISTRICT\tSTATE_NAME\tDISTRICT_NAME\tGENHH_ALLP\tGENHH_APPPC")
                if giscatlab[catList[j]]=="All"; println(f, "\tANEXPPC\tPOP")
                else println(f)
                end
                for nt in nts
                    ntidx = findfirst(x->x==nt, tbntlist)
                    print(f, nt,"\t",n,"\t",nt,"\t",natName[n],"\t",nuts[y][nt],"\t")
                    printfmt(f, "{:f}", gre[ntidx,j]); print(f, "\t",gre[ntidx,j]/gisTotPop[y][ntidx])
                    if giscatlab[catList[j]]=="All"; println(f,"\t",gisAvgExp[y][ntidx],"\t",convert(Int, gisTotPop[y][ntidx]))
                    else println(f)
                    end
                end
                close(f)

                mkpath(fp*"CFAC/")
                f = open(fp*"CFAC/CFAC_"*giscatlab[catList[j]]*".txt","w")
                println(f, "KEY_CODE\tSTATE\tDISTRICT\tSTATE_NAME\tDISTRICT_NAME\tGENHH_ALLPPC")
                for nt in nts
                    ntidx = findfirst(x->x==nt, tbntlist)
                    print(f, nt,"\t",n,"\t",nt,"\t",natName[n],"\t",nuts[y][nt],"\t")
                    if rank; println(f, gredrpc[ntidx,j]) else println(f, gredpc[ntidx,j]) end
                end
                close(f)
            end
        end
    end
end

function analyzeCategoryComposition(year, output="")
    global sec, hhid, cat, catlist
    global indirectCE, ieHHs

    nhc = 5 # number of high composition sectors

    nh = length(hhid)
    ns = length(sec)
    nc = length(catlist)

    e = indirectCE[year]         # {India sectors, hhid}}
    ec = ieHHs[year]     # {hhid, category}

    te = [sum(e[i,:]) for i=1:ns]
    tec = [sum(ec[:,i]) for i=1:nc]
    # make index dictionaries
    indCat = [findfirst(x->x==cat[year][s], catlist) for s in sec]

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
    global ieHHs, emissionsRel

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
    ec = ieHHs[year]
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
    global ieHHs, emissionsInc

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
    ec = ieHHs[year]
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
    global ieHHs, emissionsRng

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
    ec = ieHHs[year]
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
    global ieHHs, emissionsLev

    if !absintv && length(intv) == 0; intv = [0.25,0.25,0.25,0.25]
    elseif !absintv && sum(intv) != 1; sumi = sum(intv); intv /= sumi
    end

    ec = ieHHs[year]
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
    ec = ieHHs[year]
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
    global ieHHs, emissionsIncRel

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
    ec = ieHHs[year]
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
    global ieHHs, emissionCostDis

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
    eh = ieHHs[year]
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
    global ieHHs, emissionCostDis

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
    eh = ieHHs[year]
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
    global ieHHs, emissionCostDis

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
    eh = ieHHs[year]
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

function initiate()
    hhsList = Dict{Int, Dict{String, Array{String, 1}}}()   # Household ID: {year, {nation, {hhid}}}
    sec = Array{String, 1}()                # Consumption products' or services' sectors
    secName = Dict{String, String}()        # sector name dictionary: {sector code, name}

    cat = Dict{Int, Dict{String, String}}()         # category dictionary: {sector code, category}
    nat = Dict{Int, Dict{String, String}}()         # hhid's nation: {year, {hhid, nation code}}
    reg = Dict{Int, Dict{String, String}}()         # hhid's NUTS: {year, {hhid, NUTS code}}
    typ = Dict{Int, Dict{String, String}}()         # hhid's sector type, urban or rural: {year, {hhid, "urban" or "rural"}}
    siz = Dict{Int, Dict{String, Int}}()            # hhid's family size: {year, {hhid, number of members}}
    eqs = Dict{Int, Dict{String, Float64}}()        # hhid's family equivalent size (OECD scale): {year, {hhid, number of members}}
    meqs = Dict{Int, Dict{String, Float64}}()       # hhid's family equivalent size (modified OECD scale): {year, {hhid, number of members}}
    inc = Dict{Int, Dict{String, Float64}}()        # hhid's income: {year, {hhid, total income}}
    exp = Dict{Int, Dict{String, Float64}}()        # hhid's domestic expenditure: {year, {hhid, total domestic expenditure}}
    pds = Dict{Int, Dict{String, Int}}()            # hhid region's population density: {year, {hhid, district's population density}}
    rel = Dict{Int, Dict{String, Int}}()            # hhid's religion: {year, {hhid, religion code}}
    wgh = Dict{Int, Dict{String, Float64}}()        # hhid's weight: {year, {hhid, weight}}
    wghNuts = Dict{Int, Dict{String, Float64}}()    # hhid's NUTS weight: {year, {hhid, weight}}

    nutsLv = 0      # NUTS level
    nuts = Dict{Int, Dict{String, String}}()       # NUTS: {year, {code, label}}
    nutsList = Dict{Int, Dict{String, Array{String, 1}}}()      # NUTS code list: {year, {nation_code, {NUTS_code}}}
    pop = Dict{Int, Dict{String, Float64}}()        # Population: {year, {NUTS_code, population}}
    popList = Dict{Int, Dict{String, Dict{String, Float64}}}()  # Population list: {year, {nation_code, {NUTS_code, population}}}
    poplb = Dict{Int, Dict{String, String}}()       # populaton NUTS label: {year, {NUTS_code, NUTS_label}}

    directCE = Dict{Int, Dict{String, Array{Float64, 2}}}()     # direct carbon emission: {year, {nation, {table}}}
    indirectCE = Dict{Int, Dict{String, Array{Float64, 2}}}()   # indirect carbon emission: {year, {nation, {table}}}
    integratedCF = Dict{Int, Dict{String, Array{Float64, 2}}}() # carbon footprint: {year, {nation, {table}}}

    yrList = Array{Int, 1}()        # year list
    catList = Array{String, 1}()    # category list
    deCatList = Dict{Int, Array{String, 1}}()  # DE category list: {year, {DE category}}
    natList = Array{String, 1}()    # nation list
    natName = Dict{String,String}() # nation's code and full-name: {nation code, full-name}
    natA3 = Dict{String,String}()   # nation's A3: {nation code, A3}

    ieHHs = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied indirect emission by household: {year, {nation, {hhid, category}}}
    ieReg = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied indirect emission by region: {year, {nation, {region, category}}}
    ieRegDiff = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied indirect emission differences by region: {year, {nation, {region, category}}}
    deHHs = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied direct emission by household: {year, {nation, {hhid, category}}}
    deReg = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied direct emission by region: {year, {nation, {region, category}}}
    deRegDiff = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied direct emission differences by region: {year, {nation, {region, category}}}
    cfHHs = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied carbon footprint by household: {year, {nation, {hhid, category}}}
    cfReg = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied carbon footprint by region: {year, {nation, {region, category}}}
    cfRegDiff = Dict{Int, Dict{String, Array{Float64, 2}}}() # categozied carbon footprint differences by region: {year, {nation, {region, category}}}

    gisNutsList = Dict{Int, Array{String, 1}}()           # NUTS list: {year, {region(hbscd)}}
    gisRegionalIe = Dict{Int, Array{Float64, 2}}()        # categozied indirect emission by district: {year, {region(hbscd), category}}
    gisRegionalIeRank = Dict{Int, Array{Int, 2}}()        # categozied indirect emission rank by district: {year, {region(hbscd), category}}
    gisRegionalIePerCap = Dict{Int, Array{Float64, 2}}()  # categozied indirect emission per capita by district: {year, {region(hbscd), category}}
    gisRegionalIeRankPerCap = Dict{Int, Array{Int, 2}}()  # categozied indirect emission per capita rank by district: {year, {region(hbscd), category}}
    gisRegionalDe = Dict{Int, Array{Float64, 2}}()        # categozied direct emission by district: {year, {region(hbscd), DE category}}
    gisRegionalDePerCap = Dict{Int, Array{Float64, 2}}()  # categozied direct emission per capita by district: {year, {region(hbscd), DE category}}
    gisRegionalDeRank = Dict{Int, Array{Int, 2}}()        # categozied direct emission rank by district: {year, {region(hbscd), DE category}}
    gisRegionalDeRankPerCap = Dict{Int, Array{Int, 2}}()  # categozied direct emission per capita rank by district: {year, {region(hbscd), DE category}}

    gisRegionalCF = Dict{Int, Array{Float64, 2}}()        # categozied carbon footprint by district: {year, {region(hbscd), category}}
    gisRegionalCFrank = Dict{Int, Array{Int, 2}}()        # categozied carbon footprint rank by district: {year, {region(hbscd), category}}
    gisRegionalCFperCap = Dict{Int, Array{Float64, 2}}()  # categozied carbon footprint per capita by district: {year, {region(hbscd), category}}
    gisRegionalCFrankPerCap = Dict{Int, Array{Int, 2}}()  # categozied carbon footprint per capita rank by district: {year, {region(hbscd), category}}
    gisRegionalCFdiff = Dict{Int, Array{Float64, 2}}()    # differences of categozied carbon footprint by district: (emission-mean)/mean, {year, {district(GID), category}}
    gisRegionalCFdiffRank = Dict{Int, Array{Int, 2}}()    # difference ranks of categozied carbon footprint by district: (emission-mean)/mean, {year, {district(GID), category}}
    gisRegionalCFdiffPerCap = Dict{Int, Array{Float64, 2}}()  # differences of categozied carbon footprint per capita by district: (emission-mean)/mean, {year, {district(GID), category}}
    gisRegionalCFdiffRankPerCap = Dict{Int, Array{Int, 2}}()  # difference ranks of categozied carbon footprint per capita by district: (emission-mean)/mean, {year, {district(GID), category}}

    gisTotPop = Dict{Int, Array{Float64, 1}}()      # GIS version, total population by NUTS
    gisSamPop = Dict{Int, Array{Float64, 1}}()      # GIS version, total sample members by NUTS
    gisAvgExp = Dict{Int, Array{Float64, 1}}()      # GIS version, average expenditure by NUTS

    regList = Dict{Int, Array{String, 1}}() # district list

    sam = Dict{Int, Dict{String, Tuple{Int,Int}}}() # sample population and households by districct: {district code, (population, number of households)}
    ave = Dict{Int, Dict{String, Float64}}()        # average annual expenditure per capita, USD/yr: {district code, mean Avg.Exp./cap/yr}

    popcd = Dict{Int, Dict{String, String}}()       # concordance NUTS code: {year, {NUTS code, replaced population NUTS code}}
    pophbscd = Dict{Int, Dict{String, String}}()    # concordance NUTS code: {year, {population NUTS code, HBS NUTS code}}
    hbscd = Dict{Int, Dict{String, String}}()       # concordance NUTS code: {year, {NUTS code, HBS NUTS code}}
    popgiscd = Dict{Int, Dict{String, String}}()    # concordance NUTS code: {year, {population NUTS code, GIS NUTS code}}
    popcdlist = Dict{Int, Array{String, 1}}()       # Population NUTS code list: {year, Population NUTS code}
    ntcdlist = Dict{Int, Array{String, 1}}()        # NUTS code list: {year, NUTS code}
    gispopcdlist = Dict{Int, Dict{String, Array{String, 1}}}()   # Population NUTS code list: {year, {GIS_NUTS coded, {Population NUTS code}}}
    giscdlist = Dict{Int, Array{String, 1}}()       # GIS NUTS code list: {year, GIS NUTS code}
    hbspopcdlist = Dict{Int, Dict{String, Array{String, 1}}}()   # Population NUTS code list: {year, {HBS_NUTS coded, {Population NUTS code}}}
    hbscdlist = Dict{Int, Array{String, 1}}()       # HBS NUTS code list: {year, HBS NUTS code}
    majorCity = Dict{Int, Dict{String, String}}()   # major city of NUTS: {year, {NUTS_upper, major sub-NUTS}}
    gisCoord = Dict{Int, Dict{String, Tuple{Float64, Float64}}}()   # GIS coordinates: {year, {GIS_NUTS, {X, Y}}}
    giscatlab = Dict{String, String}()              # Web-exporting category matching table: {Category label in program, in web-site files}

end

end
