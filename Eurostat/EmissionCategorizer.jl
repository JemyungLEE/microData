module EmissionCategorizer

# Developed date: 3. Aug. 2020
# Last modified date: 19. May. 2022
# Subject: Categorize EU households' carbon footprints
# Description: Read household-level CFs and them by consumption category, district, expenditure-level, and etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using XLSX
using Statistics
using Formatting: printfmt

hhsList = Dict{Int, Dict{String, Array{String, 1}}}()   # Household ID: {year, {nation, {hhid}}}
sec = Dict{Int, Array{String, 1}}()                     # Consumption products' or services' sectors: {year, {sector}}
secName = Dict{Int, Dict{String, String}}()             # sector name dictionary: {year, {sector code, name}}

# hhid -> nation(2digit)_hhid: some HHIDs are duplicated across multiple countries
cat = Dict{Int, Dict{String, String}}()         # category dictionary: {year, {sector code, category}}
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

nutsLv = 0                                              # NUTS level
nuts = Dict{Int, Dict{String, String}}()                # NUTS: {year, {code, label}}
nutsList = Dict{Int, Dict{String, Array{String, 1}}}()  # NUTS code list: {year, {nation_code, {NUTS_code}}}
nuts_intg = Dict{Int, Dict{String, String}}()           # integrated NUTS codes: {target_year, {target_NUTS, concording_NUTS}}
nuts_intg_list = Array{String, 1}()                     # integrated NUTS list


pop = Dict{Int, Dict{String, Float64}}()                # Population: {year, {NUTS_code, population}}
pops_ds = Dict{Int, Dict{String, Dict{Int, Float64}}}() # Population by population density: {year, {NUTS_code, {density, population}}}
pops_ds_hbs = Dict{Int, Dict{String, Dict{Int, Float64}}}() # Population by population density: {year, {HBS NUTS_code, {density, population}}}
popList = Dict{Int, Dict{String, Dict{String, Float64}}}()  # Population list: {year, {nation_code, {NUTS_code, population}}}
poplb = Dict{Int, Dict{String, String}}()                   # populaton NUTS label: {year, {NUTS_code, NUTS_label}}

directCE = Dict{Int, Dict{String, Array{Float64, 2}}}()     # direct carbon emission: {year, {nation, {table}}}
indirectCE = Dict{Int, Dict{String, Array{Float64, 2}}}()   # indirect carbon emission: {year, {nation, {table}}}
integratedCF = Dict{Int, Dict{String, Array{Float64, 2}}}() # carbon footprint: {year, {nation, {table}}}

yrList = Array{Int, 1}()                    # year list
catList = Array{String, 1}()                # category list
deCatList = Dict{Int, Array{String, 1}}()   # DE category list: {year, {DE category}}
natList = Dict{Int, Array{String, 1}}()     # nation list: {year, {nation}}
regList = Dict{Int, Array{String, 1}}()     # district list
natName = Dict{String,String}()             # nation's code and full-name: {nation code, full-name}
natA3 = Dict{String,String}()               # nation's A3: {nation code, A3}

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
ntsmp = Dict{Int,Dict{String,Dict{String,Array{Float64,1}}}}()  # NUTS sample size:{year,{nation,{NUTS,sample number{total,dense,mid,sparse,none}}}}
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

    global hhsList, natList, catList, natName, siz, wgh, wghNuts, ieHHs, deHHs, cfHHs
    if nuts_mode; w = wghNuts else w = wgh end

    nn = length(natList[year])
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

    ci = findfirst(x -> x == "Total", catList)
    for i=1:nn
        n = natList[year][i]
        for j = 1:length(hhsList[year][n])
            h = hhsList[year][n][j]
            ie = ieHHs[year][n][j,ci]
            de = deHHs[year][n][j,ci]
            cf = cfHHs[year][n][j,ci]

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
        print(f, natList[year][i],"\t",length(hhsList[year][natList[year][i]]),"\t",natsam[i],"\t", pop[year][natList[year][i]],"\t",natwgh[i])
        print(f, "\t",natcf[i],"\t",natcfpc[i],"\t",natcfph[i],"\t",natcfpeqs[i],"\t",natcfpmeqs[i])
        print(f, "\t", natie[i], "\t", natiepc[i], "\t", natde[i], "\t", natdepc[i])
        println(f)
    end
    close(f)
end

function readCategoryData(inputFile, year, ntlv=0; subCategory="", except=[])

    global sec, deSec, cat, gid, nam, pop, poplb, gidData, misDist
    global natList, natName, natA3, nuts, nutslList, pop, popList, catList
    global popcd, pophbscd, hbscd, gisCoord, popgiscd, popcdlist, ntcdlist
    global giscdlist, gispopcdlist, giscatlab, hbspopcdlist, hbscdlist
    global nutsLv = ntlv
    nuts_yr = Dict(2010=>"2010", 2015=>"2013")
    pop_yr = Dict(2010=>"2021", 2015=>"2021")

    xf = XLSX.readxlsx(inputFile)
    if isa(year, Number); year = [year] end
    cdlen = ntlv+2

    sh = xf["Nation"]
    for y in year
        natList[y] = Array{String, 1}()
        for r in XLSX.eachrow(sh)
            if XLSX.row_number(r)>1 && !ismissing(r[1]) && !(string(r[1]) in natList[y])
                push!(natList[y], string(r[1]))
                natName[natList[y][end]] = string(r[2])
                natA3[natList[y][end]] = string(r[3])
            end
        end
    end


    for y in year
        nuts[y], nutsList[y] = Dict{String, String}(), Dict{String, Array{String, 1}}()
        pop[y], poplb[y], popList[y] = Dict{String, Float64}(), Dict{String, String}(), Dict{String, Dict{String, Float64}}()
        for n in natList[y]; nutsList[y][n], popList[y][n] = Array{String, 1}(), Dict{String, Float64}() end
        pophbscd[y], popcd[y], hbscd[y] = Dict{String, String}(), Dict{String, String}(), Dict{String, String}()
        popgiscd[y], popcdlist[y] = Dict{String, String}(), Array{String, 1}()
        ntcdlist[y], giscdlist[y], hbscdlist[y] = Array{String, 1}(), Array{String, 1}(), Array{String, 1}()
        gispopcdlist[y], hbspopcdlist[y] = Dict{String, Array{String, 1}}(), Dict{String, Array{String, 1}}()
        gisCoord[y] = Dict{String, Tuple{Float64, Float64}}()

        sec[y], secName[y], cat[y] = Array{String, 1}(), Dict{String, String}(), Dict{String, String}()
        tb = xf["Sector_" * string(y)][:]
        ci = findfirst(x -> x == (length(subCategory) == 0 ? "Category" : subCategory * "_category"), tb[1,:])
        for i in filter(x -> !ismissing(tb[x,1]), collect(2:size(tb,1)))
            s = string(tb[i,1])
            push!(sec[y], s)
            secName[y][s] = string(tb[i,2])
            if !ismissing(tb[i,ci]) && !(string(tb[i,ci]) in except); cat[y][s] = string(tb[i,ci]) end
        end

        tb = xf["NUTS" * nuts_yr[y]][:]
        for i in filter(x -> !ismissing(tb[x,1]), collect(2:size(tb,1)))
            lv, ntcd, n = parse(Int, string(tb[i,3])), string(tb[i,1]), string(tb[i,4])
            if length(ntcd)==lv+2; nuts[y][ntcd] = ('/' in tb[i,2] ? strip(split(tb[i,2], '/')[1]) : strip(tb[i,2]))
            else println("NUTS level error: ", ntcd, "\t", lv)
            end
            if n in natList[y] && lv == ntlv && ntcd[end] != 'Z'; push!(nutsList[y][n], ntcd) end
        end

        # # add aggregated regions for region-code absent nations
        # # agg_reg = ["FR"]
        # agg_reg = natList[y]
        # for n in filter(x -> x in natList[y], agg_reg)
        #     rg = n
        #     for i=1:ntlv; rg *= "0" end
        #     if !(rg in nutsList[y][n])
        #         push!(nutsList[y][n], rg)
        #         nuts[y][rg] = nuts[y][n]
        #     end
        # end

        tb = xf["Conc" * nuts_yr[y]][:]
        hidx = findfirst(x -> x == "HBS_code", tb[1,:])
        pidx = findfirst(x -> x == "Population_code", tb[1,:])
        for i in filter(x -> !ismissing(tb[x,1]), collect(2:size(tb,1)))
            ntcd = string(tb[i,1])
            popcd[y][ntcd] = string(tb[i,pidx])
            hbscd[y][ntcd] = string(tb[i,hidx])
            if !(hbscd[y][ntcd] in ntcdlist[y]); push!(ntcdlist[y], hbscd[y][ntcd]) end
        end

        tb = xf["PopCd" * pop_yr[y]][:]
        hbs_i, gis_i = [findfirst(x -> x == lb * string(y), tb[1,:]) for lb in ["NUTS_", "GIS_"]]
        for i in filter(x -> !ismissing(tb[x,1]), collect(2:size(tb,1)))
            nt_map = string(tb[i,1])
            if !ismissing(tb[i,hbs_i]); pophbscd[y][nt_map] = string(tb[i,hbs_i]) end
            if !ismissing(tb[i,gis_i]); popgiscd[y][nt_map] = string(tb[i,gis_i]) end

            # nt_map, nt_hbs, nt_gis = string(tb[i,1]), string(tb[i,hbs_i]), string(tb[i,gis_i])
            # pophbscd[y][nt_map], popgiscd[y][nt_map] = nt_hbs, nt_gis

            # if !(nt_map in popcdlist[y]); push!(popcdlist[y], nt_map) end
            # if !haskey(hbspopcdlist[y], nt_hbs); hbspopcdlist[y][nt_hbs] = Array{String, 1}() end
            # push!(hbspopcdlist[y][nt_hbs], nt_map)
            # if !haskey(gispopcdlist[y], nt_gis); gispopcdlist[y][nt_gis] = Array{String, 1}() end
            # push!(gispopcdlist[y][nt_gis], nt_map)
        end
        popcdlist[y] = sort(collect(keys(pophbscd[y])))
        for hc in unique(collect(values(pophbscd[y]))); hbspopcdlist[y][hc] = sort(filter(x -> pophbscd[y][x] == hc, popcdlist[y])) end
        for gc in unique(collect(values(popgiscd[y]))); gispopcdlist[y][gc] = sort(filter(x -> popgiscd[y][x] == gc, popcdlist[y])) end

        # hbscdlist[y] = sort(filter(x->length(x)==cdlen && haskey(nutsList[y], x[1:2]) && x in nutsList[y][x[1:2]], collect(keys(hbspopcdlist[y]))))
        hbscdlist[y] = sort(filter(x->length(x)==cdlen, collect(keys(hbspopcdlist[y]))))
        giscdlist[y] = sort(filter(x->length(x)==cdlen, collect(keys(gispopcdlist[y]))))
    end

    tb = xf["Pop_mod_NUTS3"][:]
    yi = [findfirst(x -> x == y, tb[1,:]) for y in year]
    for ri = 2:size(tb,1)
        ntcd, ntlb = strip.(split(replace(tb[ri, 1], ",Total,Total,Number"=>""), " - "))    # code, label
        if '/' in ntlb; ntlb = strip(split(ntlb, '/')[1]) end
        n = ntcd[1:2]
        for i=1:length(year)
            y = year[i]
            if n in natList[y]
                s = strip(replace(string(tb[ri, yi[i]]), ['b','d','e','p',',']=>""))
                if tryparse(Float64, s) != nothing && !haskey(pop[y], ntcd); pop[y][ntcd] = parse(Float64, s) end
                if !haskey(poplb[y], ntcd); poplb[y][ntcd] = haskey(nuts[y], ntcd) ? nuts[y][ntcd] : ntlb end
            end
        end
    end

    for y in year
        sh = xf["GIS_coor" * nuts_yr[y]]
        for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1; gisCoord[y][r[1]] = (r[2], r[3]) end end
    end

    sh = xf["GIS_cat"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1; giscatlab[string(r[1])] = string(r[2]) end end
    close(xf)

    for y in year
        if length(catList) == 0; catList = sort(unique(values(cat[y])))
        elseif sort(filter(x -> x != "Total", catList)) != sort(unique(values(cat[y]))); println("Categories are different in ", year, "\t", sort(catList), "\t", sort(unique(values(cat[y]))))
        end
    end
    if length(subCategory) == 0 && !("Total" in catList); push!(catList, "Total")
    elseif length(subCategory) > 0 && !(subCategory in catList); push!(catList, subCategory)
    end
end

function readPopulation(years, inputFile; nuts_lv = 1, adjust = true, reverse_order = false)
    # read nuts_lv region's population if there is a corresponding NUTS_2010 or 2013 matching NUTS_2021 of population data
    # if there is not, go further level (NUTS 1 -> NUTS 2 -> NUTS 3)
    # if reverse_order = true: go upper level (NUTS 3 -> NUTS 2 -> NUTS 1)

    global natList, popList, regList, nutsList, pop, pophbscd

    xf = XLSX.readxlsx(inputFile)
    if isa(years, Number); years = [years] end
    nt_len = nuts_lv + 2

    tb = xf["Pop_mod_NUTS3"][:]
    for y in years
        popList[y] = Dict{String, Dict{String, Float64}}()
        for n in natList[y]; popList[y][n] =  Dict(nutsList[y][n] .=> 0.0) end
    end

    tb = tb[2:end,:]
    nr = size(tb,1)
    pnts = [strip(split(tb[ri, 1], " - ", limit = 2)[1]) for ri = 1:nr]

    if !reverse_order
        chk = [false for ri = 1:nr]
        for i = 1:length(years)
            y = years[i]
            for cd_len = nt_len:5, ri in filter(x -> !chk[x] && length(pnts[x]) == cd_len, 1:nr)
                n, nt = pnts[ri][1:2], pnts[ri]
                if haskey(pophbscd[y], nt) && pophbscd[y][nt] in regList[y] && pophbscd[y][nt] in nutsList[y][n] && haskey(pop[y], nt) && pop[y][nt] > 0
                    popList[y][n][pophbscd[y][nt]] += pop[y][nt]
                    chk[filter(x -> startswith(pnts[x], nt), 1:nr)] .= true
                end
            end
        end
    else
        for y in years, cd_len = 5:-1:nt_len
            for ri in findall(x -> length(x) == cd_len && haskey(pophbscd[y], x) && pophbscd[y][x] in regList[y] && popList[y][x[1:2]][pophbscd[y][x]] == 0, pnts)
                nt = pnts[ri]
                if haskey(pop[y], nt); popList[y][nt[1:2]][pophbscd[y][nt]] += pop[y][nt] end
            end
        end
    end

    if adjust
        for y in years, n in natList[y]
            tot_pop = sum(collect(values(popList[y][n])))
            gap = pop[y][n] - tot_pop
            if abs(gap) > 0; for nt in collect(keys(popList[y][n])); popList[y][n][nt] += gap * popList[y][n][nt] / tot_pop end end
        end
    end
    close(xf)
end

function readPopGridded(year, inputFile; nuts_lv = [0], adjust = false, tag = ["_dense", "_inter", "_spars", "_total"])
    # 1:Densely populated (at least 500), 2:Intermediate (between 100 and 499)
    # 3:Sparsely populated (less than 100), 4:Total, (Unit: inhabitants/km2)

    global pop, pops_ds, pops_ds_hbs, hhsList, hbscd, popList
    ntag = length(tag)
    if isa(nuts_lv, Number); nuts_lv = [nuts_lv] end
    nuts_yr = Dict(2010=>"2010", 2015=>"2013")
    gp_tag = Dict(0 => "Pop_GP", 1 => "Pop_GP_LV1")

    xf = XLSX.readxlsx(inputFile)
    if isa(year, Number); year = [year] end

    for lv in nuts_lv, y in year
        tb = xf[replace(gp_tag[lv], "GP" => "GP" * nuts_yr[y])][:]
        if !haskey(pops_ds, y); pops_ds[y] = Dict{String, Dict{Int, Float64}}() end
        if !haskey(pops_ds_hbs, y); pops_ds_hbs[y] = Dict{String, Dict{Int, Float64}}() end

        idx = [findfirst(x -> x == string(y) * t, tb[1,:]) for t in tag]
        for ri = 2:size(tb,1)
            nt = string(tb[ri, 1])
            if haskey(hhsList[y], nt[1:2])
                pops_ds[y][nt] = Dict{Int, Float64}()
                for i = 1:ntag; pops_ds[y][nt][i] = tb[ri,idx[i]] end
            end
        end
        for nt in collect(keys(pops_ds[y]))
            nt_hbs = hbscd[y][nt]
            if !haskey(pops_ds_hbs[y], nt_hbs); pops_ds_hbs[y][nt_hbs] = Dict(1:ntag .=> 0.0) end
            for i = 1:ntag; pops_ds_hbs[y][nt_hbs][i] += pops_ds[y][nt][i] end
        end

        if adjust
            for nt in sort(filter(x -> length(x) == lv+2 ,collect(keys(pops_ds_hbs[y]))))
                n = nt[1:2]
                tot_gp = sum([pops_ds_hbs[y][nt][i] for i = 1:3])
                if lv == 0 && haskey(pop[y], nt)
                    gap = pop[y][nt] - tot_gp
                    for i = 1:ntag; pops_ds_hbs[y][nt][i] += gap * pops_ds_hbs[y][nt][i] / tot_gp end
                elseif haskey(popList[y], n) && haskey(popList[y][n], nt) && popList[y][n][nt] > 0
                    gap = popList[y][n][nt] - tot_gp
                    for i = 1:ntag; pops_ds_hbs[y][nt][i] += gap * pops_ds_hbs[y][nt][i] / tot_gp end
                else println("gridded population adjustment error: ", y, " year, ", nt)
                end
            end
        end
    end
    close(xf)
end

function setCategory(list::Array{String,1})
    global catList
    if sort(catList)==sort(list); catList = list
    else println("Category items are different.\n",sort(catList),"\n",sort(list))
    end
end

function readHouseholdData(inputFile; period="annual", sampleCheck=false, alter=false, remove_nt=false, remove_z=false)
    # period: "annual"(default), "monthly", or "daily"

    global yrList, hhsList, natList, nutsList, regList, relList, disSta, nutsLv, gisCoord
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
        nat[year][hh], reg[year][hh] = s[2], hbscd[year][s[4][1:(nutsLv+2)]]
        siz[year][hh] = parse(Int,s[5])
        eqs[year][hh], meqs[year][hh] = parse(Float64,s[12]), parse(Float64,s[13])
        wgh[year][hh], pds[year][hh] = parse(Float64,s[6]), parse(Int,s[11])
        inc[year][hh], exp[year][hh] = parse(Float64,s[7]), parse(Float64,s[9])
    end
    close(f)

    for y in yrList
        filter!(x->x in unique(collect(values(nat[y]))), natList[y])
        regList[y] = sort(unique(collect(values(reg[y]))))

        if alter
            gis_nts = sort(collect(keys(gisCoord[y])))
            nts = filter(x -> !(x in gis_nts) && x[end] != 'Z', regList[y])
            for nt in nts
                ntz = nt[1:end-1] * "Z"
                for h in filter(x -> reg[y][x] == nt, hhsList[y][nt[1:2]]); reg[y][h] = ntz end
            end
            regList[y] = sort(unique(collect(values(reg[y]))))
        end
        if remove_z
            for nt in filter(x -> x[end] == 'Z', regList[y])
                hhz_ls = filter(x -> reg[y][x] == nt, hhsList[y][nt[1:2]])
                filter!(x -> !(x in hhz_ls), hhsList[y][nt[1:2]])
                filter!(x -> !(x.first in hhz_ls), nat[y])
                filter!(x -> !(x.first in hhz_ls), reg[y])
                filter!(x -> !(x.first in hhz_ls), typ[y])
                filter!(x -> !(x.first in hhz_ls), siz[y])
                filter!(x -> !(x.first in hhz_ls), eqs[y])
                filter!(x -> !(x.first in hhz_ls), meqs[y])
                filter!(x -> !(x.first in hhz_ls), pds[y])
                filter!(x -> !(x.first in hhz_ls), rel[y])
                filter!(x -> !(x.first in hhz_ls), wgh[y])
            end
            filter!(x -> x[end] != 'Z', regList[y])
        end
        if remove_nt; for n in natList[y]; filter!(x->x in regList[y], nutsList[y][n]) end end

        # convert household's income and expenditure data period
        if period=="daily"; cvr = 1/365 elseif period=="monthly" cvr = 30/365 end
        if period in ["daily", "monthly"]
            for n in collect(keys(hhsList[y]))
                for h in hhsList[y][n]; inc[y][h] = inc[y][h] * cvr; exp[y][h] = exp[y][h] * cvr end
            end
        end

        # check sample numbers by district's population density
        if sampleCheck
            samples = zeros(Int, length(regList[y]), 4)    # {total, high_dens, middle_dens, low_dens}
            for n in collect(keys(hhsList[y]))
                for h in hhsList[y][n]
                    idx = findfirst(x -> x == hbscd[y][reg[y][h]], regList[y])
                    samples[idx,1] += siz[y][h]
                    samples[idx,pds[h]+1] += siz[y][h]
                end
            end
            f = open(Base.source_dir() * "/data/extracted/SampleNumber_" * string(y) * ".txt","w")
            println(f,"District\tAll\tHigh_density\tMiddle_density\tLow_density")
            for i=1:length(regList[y]); print(f,regList[y][i]); for sn in samples[i,:]; println(f,"\t",sn) end end
            close(f)
        end
    end
end

function readEmissionData(year, nations, inputFiles; mode = "ie")
    # mode: [ie] indirect emission, [de] direct emission

    global sec, hhsList, indirectCE, directCE
    if lowercase(mode) == "ie"; indirectCE[year] = Dict{String, Array{Float64, 2}}()
    elseif lowercase(mode) == "de"; directCE[year] = Dict{String, Array{Float64, 2}}()
    end

    ns = length(sec[year])
    nn = length(nations)

    if length(inputFiles) == nn
        for i = 1:nn
            n = nations[i]
            nh = length(hhsList[year][n])
            e = zeros(Float64, ns, nh)

            f = open(inputFiles[i])
            hhs_ls = string.(split(readline(f), '\t'))[2:end]

            hidxs = [findfirst(x -> x== h[4:end] , hhs_ls) for h in hhsList[year][n]]
            if any(hidxs .== nothing); println("HHs index mismatch error: ", count(x -> x == nothing, hidxs)) end

            for l in eachline(f)
                l = string.(split(l, '\t'))
                e[findfirst(x->x==string(l[1]), sec[year]),:] = map(x->parse(Float64,x), l[2:end])[hidxs]
            end
            if lowercase(mode) == "ie"; indirectCE[year][n] = e
            elseif lowercase(mode) == "de"; directCE[year][n] = e
            end
            close(f)
        end
    else println("Sizes of nation list (", nn,") and emission files (", length(inputFiles),") do not match")
    end
end

function readExpenditureData(year, nations, inputFiles)

    # global sec, hhsList, indirectCE
    # indirectCE[year] = Dict{String, Array{Float64, 2}}()
    # ns,nn = length(sec[year]), length(nations)
    #
    # if length(inputFiles) == nn
    #     for i = 1:nn
    #         n = nations[i]
    #         nh = length(hhsList[year][n])
    #         f = open(inputFiles[i])
    #         readline(f)
    #         e = zeros(Float64, ns, nh)
    #         for l in eachline(f)
    #             l = split(l, '\t')
    #             e[findfirst(x->x==string(l[1]), sec[year]),:] = map(x->parse(Float64,x), l[2:end])
    #         end
    #         indirectCE[year][n] = e
    #         close(f)
    #     end
    # else println("Sizes of nation list (", nn,") and emission files (", length(inputFiles),") do not match")
    # end
end

function integrateCarbonFootprint(year; mode="cf")

    global natList, hhsList, sec, directCE, indirectCE
    global integratedCF[year] = Dict{String, Array{Float64, 2}}()

    ie, de = indirectCE[year], directCE[year]

    for n in filter(x -> haskey(hhsList[year], x), natList[year])
        nh = length(hhsList[year][n])
        if mode == "cf"; integratedCF[year][n] = ie[n] + de[n]
        elseif mode == "ie"; integratedCF[year][n] = ie[n]
        elseif mode == "de"; integratedCF[year][n] = de[n]
        end
    end
end

function importIntegratedNUTS(nutsDict, nutsList)
    global nuts_intg = nutsDict
    global nuts_intg_list = nutsList
end

function calculateTemporalGap(target_year, base_year, outputFile, nations=[]; mode="cf", nspan=128, minmax=[], descend=false,logarithm=false, tag="NUTS")

    ty, by, int_yr = target_year, base_year, 10000 * base_year + target_year
    if length(nations) == 0; nats = natList[by]
    elseif isa(nations, String); nats = [nations]
    elseif isa(nations, Array{String, 1}); nats = nations
    end

    global nuts_intg, nuts_intg_list, catList, popList, gisNutsList, hbscd, ave, gisTotPop, gisAvgExp
    global gisRegionalIe, gisRegionalIeRank, gisRegionalIePerCap, gisRegionalIeRankPerCap
    global gisRegionalDe, gisRegionalDeRank, gisRegionalDePerCap, gisRegionalDeRankPerCap
    global gisRegionalCF, gisRegionalCFrank, gisRegionalCFperCap, gisRegionalCFrankPerCap

    # global nuts, nutsLv, nutsList, natList, sam, reg
    # global gispopcdlist, giscdlist, hbspopcdlist, hbscdlist

    gisNutsList[int_yr] = nuts_intg_list
    nil, nn, nc = nuts_intg_list, length(nuts_intg_list), length(catList)
    labels, labelspc = Dict{Int, Array{String,2}}(), Dict{Int, Array{String,2}}()

    if mode == "ie"; gre, grer, grepc, grerpc = gisRegionalIe, gisRegionalIeRank, gisRegionalIePerCap, gisRegionalIeRankPerCap
    elseif mode == "de"; gre, grer, grepc, grerpc = gisRegionalDe, gisRegionalDeRank, gisRegionalDePerCap, gisRegionalDeRankPerCap
    elseif mode == "cf"; gre, grer, grepc, grerpc = gisRegionalCF, gisRegionalCFrank, gisRegionalCFperCap, gisRegionalCFrankPerCap
    end

    gre[int_yr], grepc[int_yr] = zeros(Float64, nn, nc), zeros(Float64, nn, nc)
    gisPop = Dict(yr => zeros(Float64, nn) for yr in [ty, by])
    gisAve = Dict(yr => zeros(Float64, nn) for yr in [ty, by])
    gisEpc = Dict(yr => zeros(Float64, nn, nc) for yr in [ty, by])

    for i = filter(x -> gisNutsList[ty][x][1:2] in nats, 1:length(gisNutsList[ty]))
        nt = gisNutsList[ty][i]
        n = nt[1:2]
        nt_i = nuts_intg[ty][nt]
        nidx = findfirst(x -> x == nt_i, nil)
        gre[int_yr][nidx,:] += gre[ty][i,:]
        gisPop[ty][nidx] += popList[ty][n][nt]
        gisAve[ty][nidx] += ave[ty][nt] * popList[ty][n][nt]
        gisEpc[ty][nidx,:] += gre[ty][i,:]
    end
    gisAve[ty] ./= gisPop[ty]
    gisEpc[ty] ./= gisPop[ty]

    for i = filter(x -> gisNutsList[by][x][1:2] in nats && gisNutsList[by][x] in nil, 1:length(gisNutsList[by]))
        nt = gisNutsList[by][i]
        n = nt[1:2]
        nidx = findfirst(x -> x == nt, nil)
        gre[int_yr][nidx,:] -= gre[by][i,:]
        gisPop[by][nidx] += popList[by][n][nt]
        gisAve[by][nidx] += ave[by][nt] * popList[by][n][nt]
        gisEpc[by][nidx,:] += gre[by][i,:]
    end
    gisAve[by] ./= gisPop[by]
    gisEpc[by] ./= gisPop[by]

    grepc[int_yr] = gisEpc[ty] - gisEpc[by]
    gisTotPop[int_yr] = gisPop[ty] - gisPop[by]
    gisAvgExp[int_yr] = gisAve[ty] - gisAve[by]

    modeTag = uppercase(mode)
    filename = replace(replace(replace(outputFile,"YEAR_"=>string(int_yr)*"_"), ".txt"=>"_"*modeTag*".txt"), ".csv"=>"_"*modeTag*".csv")
    grer[int_yr], labels[int_yr] = exportRegionalTables(replace(filename,"_OvPcTag"=>"_overall"), tag, nil, nspan, minmax[1], gre[int_yr], logarithm, descend, tab_mode=mode)
    grerpc[int_yr], labelspc[int_yr] = exportRegionalTables(replace(filename,"_OvPcTag"=>"_percap"), tag, nil, nspan, minmax[2], grepc[int_yr], logarithm, descend, tab_mode=mode)

    return labels, labelspc
end

function calculateNutsPopulationWeight(;year = 0, pop_dens = false, adjust = true)

    global yrList, natList, regList, hhsList, nutsList, popList, pops_ds_hbs, wghNuts, hbscd
    global siz, reg, pds, ntpop, ntsmp, ntwgh
    if year > 0; yrs = [year] else yrs = yrList end

    # count region samples
    typeidx = Dict(1 => 2, 2 => 3, 3 => 4, 9 => 5)
    for y in yrs
        ntsmp[y] = Dict{String, Dict{String, Array{Float64,1}}}()
        for n in filter(x -> haskey(hhsList[y], x), natList[y])
            ntsmp[y][n] = Dict{String, Array{Float64,1}}()
            smp = ntsmp[y][n]
            for nt in nutsList[y][n]
                smp[nt] = zeros(Float64, 5)
                for h in filter(x -> reg[y][x] == nt, hhsList[y][n]); smp[nt][[1, typeidx[pds[y][h]]]] .+= siz[y][h] end
            end
        end

        if adjust
            for nt in filter(x -> x[end] == 'Z', regList[y])
                n = nt[1:2]
                smp = ntsmp[y][n]
                nts = filter(x -> x[end] != 'Z', nutsList[y][n])

                smp_tot = [sum([smp[x][i] for x in nts]) for i = 1:5]
                smp_tmp = zeros(Float64, 5)
                for h in filter(x -> reg[y][x] == nt, hhsList[y][n])
                    smp_tmp[typeidx[pds[y][h]]] += siz[y][h]
                    smp_tmp[1] += siz[y][h]
                end
                for x in nts, i = 1:5; smp[x][i] += (smp_tot[i] > 0 ? smp_tmp[i] * smp[x][i] / smp_tot[i] : 0.0) end
            end
        end
    end

    # calculate weights
    typeidx_r = Dict(2 => 1, 3 => 2, 4 => 3, 5 => 9)
    for y in yrs
        ntwgh[y] = Dict{String, Dict{String, Array{Int,1}}}()
        for n in filter(x -> haskey(hhsList[y], x), natList[y])
            ntwgh[y][n] = Dict{String, Array{Int,1}}()
            wgh, smp, pls = ntwgh[y][n], copy(ntsmp[y][n]), popList[y][n]
            nts = nutsList[y][n]

            for nt in nts
                wgh[nt] = zeros(Float64, 5)
                wgh[nt][1] = pls[nt] / smp[nt][1]
            end
            if adjust
                for nt in filter(x -> x[1:2] == n && x[end] == 'Z', regList[y])
                    wgh[nt] = zeros(Float64, 5)
                    nts = filter(x -> x[end] != 'Z', nutsList[y][n])
                    smp_tot = sum([smp[x][1] for x in nts])
                    wgh[nt][1] = sum([wgh[x][1] * smp[x][1] for x in nts]) / smp_tot
                end
            end

            if pop_dens
                for nt in nts
                    pbd = copy(pops_ds_hbs[y][nt])
                    chk = [false for i=1:3]
                    for i in filter(x -> smp[nt][x] > 0 && pbd[typeidx_r[x]] == 0, 2:4)
                        smp[nt][5] += smp[nt][i]
                        smp[nt][i] = 0
                        chk[typeidx_r[i]] = true
                    end

                    for i in filter(x -> smp[nt][x] == 0 && pbd[typeidx_r[x]] > 0, 2:4)
                        pidx = filter(x->x != i, 2:4)
                        pbd_tot = sum([pbd[typeidx_r[i_p]] for i_p in pidx])
                        for i_p in pidx
                            pbd[typeidx_r[i_p]] += (pbd_tot > 0 ? pbd[typeidx_r[i]] * pbd[typeidx_r[i_p]] / pbd_tot : 0.0)
                        end
                        pbd[typeidx_r[i]] = 0
                    end

                    if smp[nt][5] > 0 && sum(smp[nt][2:4]) > 0
                        smp_adj = [smp[nt][i] + (smp[nt][i] / smp[nt][1] * smp[nt][5]) for i = 2:4]
                        wgh[nt][2:4] = [(smp_adj[i] > 0 ?  pbd[i] / smp_adj[i] : 0) for i = 1:3]
                        wgh[nt][5] = sum(wgh[nt][2:4] .* smp[nt][2:4]) / smp[nt][1]
                    elseif smp[nt][5] > 0 && sum(smp[nt][2:4]) == 0; wgh[nt][5] = wgh[nt][1]
                    else wgh[nt][2:4] = [(smp[nt][i] > 0 ? pbd[typeidx_r[i]] / smp[nt][i] : 0) for i = 2:4]
                    end
                    for i in filter(x -> chk[typeidx_r[x]], 2:4); wgh[nt][i] = wgh[nt][5] end
                end

                if adjust
                    for nt in filter(x -> x[1:2] == n && x[end] == 'Z', regList[y])
                        nts = filter(x -> x[end] != 'Z', nutsList[y][n])
                        smp_tot = [sum([smp[x][i] for x in nts]) for i = 1:5]
                        wgh[nt] = [sum([wgh[x][i] * smp[x][i] for x in nts]) / smp_tot[i] for i = 1:5]
                    end
                end
            end
        end
    end

    # allocate household weight
    for y in yrs
        wghNuts[y] = Dict{String, Float64}()
        for n in filter(x -> haskey(hhsList[y], x), natList[y])
            for hh in hhsList[y][n]
                wghNuts[y][hh] = ntwgh[y][n][reg[y][hh]][pop_dens ? typeidx[pds[y][hh]] : 1]
            end
        end
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
    for y in year; sl[y] = sec[y] end

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
        for n in natList[y]
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
        for y in year, n in natList[y], i = 1:length(hhsList[y][n])
            hh = hhsList[y][n][i]
            print(f, y,',',n,',',hh)
            for j = 1:length(catList); print(f, ",", ht[y][n][i,j]) end
            if hhsinfo; print(f, ",",siz[y][hh],",",inc[y][hh],",",wgh[y][hh]) end
            println(f)
        end
        close(f)
    end
end

function categorizeRegionalEmission(year=[]; mode = "cf", nutsLv=1, period="annual", adjust=true, popWgh=false, ntweigh=false, religion=false)
    # mode: "ie" indirect CE, "de" direct CE, "cf" integrated CF
    # period: "annual", "monthly", or "daily"
    # religion: [true] categorize districts' features by religions

    global hhsList, natList, cat, reg, siz, inc, exp, sam, ave, rel, pop, wgh, wghNuts, catList, nutsList, relList, pophbscd
    global ieHHs, ieReg, ieRegDiff, deHHs, deReg, deRegDiff, cfHH, cfReg, cfRegDiff
    global deCat, deCatList

    if isa(year, Number); year = [year] end
    if ntweigh; pwgh = wghNuts else pwgh = wgh end
    nc = length(catList)
    nr = length(relList)

    for y in year
        er, erd = Dict{String, Array{Float64, 2}}(), Dict{String, Array{Float64, 2}}()
        sam[y], ave[y] = Dict{String, Tuple{Int,Int}}(), Dict{String, Float64}()

        for n in natList[y]
            hhs, nts = hhsList[y][n], nutsList[y][n]
            nh,nn = length(hhs), length(nts)

            # make index arrays
            ntidx = [findfirst(x -> x == reg[y][hhs[i]], nts) for i=1:nh]
            nh_idx = filter(i -> ntidx[i] != nothing, 1:nh)
            zidx = filter(x -> reg[y][hhs[x]][end] == 'Z' , 1:nh)
            if religion; relidx = [findfirst(x->x==rel[y][hhs[i]], relList) for i=1:nh] end

            # sum sample households and members by regions
            thbd = zeros(Float64, nn)   # total households by region
            tpbd = zeros(Float64, nn)   # total members of households by region

            for i in nh_idx
                thbd[ntidx[i]] += 1
                tpbd[ntidx[i]] += siz[y][hhs[i]]
            end
            for i=1:nn; sam[y][nts[i]] = (tpbd[i], thbd[i]) end

            # sum sample households and members by regions and by religions
            if religion
                thbdr = zeros(Float64, nn, nr)  # total households by district, by religion
                tpbdr = zeros(Float64, nn, nr)  # total members of households by district, by religion
                for i in nh_idx
                    thbdr[ntidx[i],relidx[i]] += 1
                    tpbdr[ntidx[i],relidx[i]] += siz[y][hhs[i]]
                end
            end

            # calculate average monthly expenditure per capita by region
            totexp = zeros(Float64, nn)     # total expenditures by region
            if popWgh
                for i in nh_idx; totexp[ntidx[i]] += exp[y][hhs[i]] * siz[y][hhs[i]] * pwgh[y][hhs[i]] end
                if adjust && length(zidx) > 0
                    zexp = sum([exp[y][hhs[i]] * siz[y][hhs[i]] * pwgh[y][hhs[i]] for i in zidx])
                    for i=1:nn; totexp[i] += zexp * popList[y][n][nts[i]] / pop[y][n] end
                end
                for i=1:nn; ave[y][nts[i]] = totexp[i] / popList[y][n][nts[i]] end
            else
                for i in nh_idx; totexp[ntidx[i]] += exp[y][hhs[i]] * siz[y][hhs[i]] end
                for i=1:nn; ave[y][nts[i]] = totexp[i] / tpbd[i] end
            end
            # convert 'AVEpC' to annual or daily
            if period=="monthly"; yytomm = 30/365; for i=1:nn; ave[y][nts[i]] = ave[y][ntd[i]] * yytomm end
            elseif period=="daily"; for i=1:nn; ave[y][nts[i]] = ave[y][nts[i]] / 365 end
            end

            # categorize emission data
            if mode == "ie"; ec = ieHHs[y][n]
            elseif mode == "de"; ec = deHHs[y][n]
            elseif mode == "cf"; ec = cfHHs[y][n]
            end
            en = zeros(Float64, nn, nc)
            if popWgh
                for i in nh_idx; en[ntidx[i],:] += ec[i,:]*pwgh[y][hhs[i]] end
                if adjust && length(zidx) > 0
                    ez = [sum([ec[i,j]*pwgh[y][hhs[i]] for i in zidx]) for j = 1:nc]
                    for i=1:nn, j=1:nc; en[i,j] += ez[j] * popList[y][n][nts[i]] / pop[y][n] end
                end
            else for i in nh_idx; en[ntidx[i],:] += ec[i,:] end
            end
            # normalizing
            if popWgh; for i=1:nn, j=1:nc; en[i,j] /= popList[y][n][nts[i]] end
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
end

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
    for y in years, n in natList[y]
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
                if expm; print(f, ",", ave[y][nts[i]]) end
                if popm; print(f, ",", popList[y][n][nts[i]]) end
                if wghm; print(f, ",", sum([pwgh[y][h]*siz[y][h] for h in filter(x->hbscd[y][reg[y][x]]==nts[i],hhsList[y][n])])) end
                # if povm; print(f, ",", disPov[regList[i]]) end
            end
            println(f)
        end
    end
    close(f)
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
        for n in natList[y]
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

        for n in natList[y]
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

        if length(outputFile) > 0
            modeTag = uppercase(mode)
            filename = replace(replace(replace(outputFile,"YEAR_"=>string(y)*"_"), ".txt"=>"_"*modeTag*".txt"), ".csv"=>"_"*modeTag*".csv")
            rank, labels[y] = exportRegionalTables(replace(filename,"_OvPcTag"=>"_overall"), tag, ntslist, nspan, minmax[1], tb, logarithm, descend, tab_mode=mode)
            rankpc, labelspc[y] = exportRegionalTables(replace(filename,"_OvPcTag"=>"_percap"), tag, ntslist, nspan, minmax[2], tbpc, logarithm, descend, tab_mode=mode)
        else rank, rankpc = zeros(Int, 0, 0), zeros(Int, 0, 0)
        end

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

function exportWebsiteFiles(year, path; nutsmode = "hbs", rank=false, empty=false, major=false)

    global natName, nuts, pop, popList, poplb, catList, gisNutsList, gisTotPop, gisAvgExp
    global pophbscd, hbscd, gispopcdlist, hbspopcdlist, majorCity, gisCoord
    global gisRegionalCF, gisRegionalCFrank, gisRegionalCFdiffPerCap, gisRegionalCFdiffRankPerCap
    if isa(year, Number); year = [year] end

    for y in year
        tbntlist = gisNutsList[y]
        if y == 2010; if nutsmode == "gis"; exceptNt = ["FR0","DE0","EL0"] elseif nutsmode == "hbs"; exceptNt = ["FR0"] end
        elseif y == 2015; exceptNt = []
        else exceptNt = []
        end
        ntslist = filter(x->!(x in exceptNt), tbntlist)

        gre, grer, gredpc, gredrpc = gisRegionalCF[y], gisRegionalCFrank[y], gisRegionalCFdiffPerCap[y], gisRegionalCFdiffRankPerCap[y]

        if y > 9999; by = trunc(Int, y/10000) else by = y end

        # find major city of NUTS1
        if major; findMajorCity(by, ntslist, nutsmode, modnuts = true) end

        # print center file
        mkpath(path*string(y))
        f = open(path*string(y)*"/centers.csv", "w")
        println(f, "\"KEY_CODE\",\"EN_NAME\",\"JA_NAME\",\"COUNTRY\",\"x\",\"y\"")

        for nt in ntslist
            println(f,"\"",nt,"\",\"",nuts[by][nt],"\",\"",nuts[by][nt],"\",\"",natName[nt[1:2]],"\",\"",gisCoord[by][nt][1],"\",\"",gisCoord[by][nt][2],"\"")
        end
        close(f)

        # print english file
        f = open(path*string(y)*"/english.txt", "w")
        println(f, "KEY_CODE\tEN_NAME")
        for nt in ntslist
            if nt[3] == '0' && nt[1:2] != "DE";
                if '(' in nuts[by][nt]; mjcity = '(' * split(nuts[by][nt], '(')[2] else mjcity = "" end
                println(f, nt, "\t", natName[nt[1:2]] * mjcity)
            else println(f, nt, "\t", nuts[by][nt],", ",natName[nt[1:2]])
            end
        end
        close(f)

        # print english_name file
        f = open(path*string(y)*"/english_match.txt", "w")
        println(f, "KEY_CODE\tSTATE_CODE\tSTATE\tDISTRICT")
        for nt in ntslist; println(f, nt, "\t", nt[1:2], "\t", natName[nt[1:2]], "\t", nuts[by][nt]) end
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
                print(f, nt,"\t",nt[1:2],"\t",nt,"\t",natName[nt[1:2]],"\t",nuts[by][nt],"\t")
                printfmt(f, "{:f}", gre[tbidx,j]); print(f, "\t",gre[tbidx,j]/gisTotPop[y][tbidx])
                if catList[j]=="Total" || catList[j]=="All"
                    println(f,"\t",gisAvgExp[y][tbidx],"\t",convert(Int, round(gisTotPop[y][tbidx], digits=0)))
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
                print(f, nt,"\t",nt[1:2],"\t",nt,"\t",natName[nt[1:2]],"\t",nuts[by][nt],"\t")
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
        for n in natList[y]
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

function initVars(; year = [], nation = [], clear_all = false)

    global hhsList, nat, reg, siz, eqs, meqs, typ, inc, exp, rel, wgh, pds
    global directCE, indirectCE, integratedCF

    if isa(year, Number); year = [year] end
    if isa(nation, String); nation = [nation] end

    if length(year) == 0
        if clear_all
            nat = Dict{Int, Dict{String, String}}()
            reg = Dict{Int, Dict{String, String}}()
            typ = Dict{Int, Dict{String, String}}()
            siz = Dict{Int, Dict{String, Int}}()
            eqs = Dict{Int, Dict{String, Float64}}()
            meqs = Dict{Int, Dict{String, Float64}}()
            inc = Dict{Int, Dict{String, Float64}}()
            exp = Dict{Int, Dict{String, Float64}}()
            pds = Dict{Int, Dict{String, Int}}()
            rel = Dict{Int, Dict{String, Int}}()
            wgh = Dict{Int, Dict{String, Float64}}()
        end
        if length(nation) == 0
            directCE = Dict{Int, Dict{String, Array{Float64, 2}}}()
            indirectCE = Dict{Int, Dict{String, Array{Float64, 2}}}()
            integratedCF = Dict{Int, Dict{String, Array{Float64, 2}}}()
        else
            for y in collect(keys(hhsList)), n in nation
                if haskey(directCE, y); directCE[y][n] = Array{Float64, 2}(undef, 0, 0) end
                if haskey(indirectCE, y); indirectCE[y][n] = Array{Float64, 2}(undef, 0, 0) end
                if haskey(integratedCF, y); integratedCF[y][n] = Array{Float64, 2}(undef, 0, 0) end
            end
        end
    else
        for y in year
            if clear_all
                nat[y] = Dict{String, String}()
                reg[y] = Dict{String, String}()
                typ[y] = Dict{String, String}()
                siz[y] = Dict{String, Int}()
                eqs[y] = Dict{String, Float64}()
                meqs[y] = Dict{String, Float64}()
                inc[y] = Dict{String, Float64}()
                exp[y] = Dict{String, Float64}()
                pds[y] = Dict{String, Int}()
                rel[y] = Dict{String, Int}()
                wgh[y] = Dict{String, Float64}()
            end
            if length(nation) == 0
                directCE[y] = Dict{String, Array{Float64, 2}}()
                indirectCE[y] = Dict{String, Array{Float64, 2}}()
                integratedCF[y] = Dict{String, Array{Float64, 2}}()
            else
                for y in year, n in nation
                    if haskey(directCE, y); directCE[y][n] = Array{Float64, 2}(undef, 0, 0) end
                    if haskey(indirectCE, y); indirectCE[y][n] = Array{Float64, 2}(undef, 0, 0) end
                    if haskey(integratedCF, y); integratedCF[y][n] = Array{Float64, 2}(undef, 0, 0) end
                end
            end
        end
    end
end

function sortHHsByCF(years, nations, outputFileTag; mode = "cf")

    global yrList, hhsList, catList, natList, nat, siz
    global ieHHs, deHHs, cfHHs, wghNuts

    mkpath(rsplit(outputFileTag, '/', limit = 2)[1])

    if length(years) == 0; yrs = yrList
    elseif isa(years, Int); yrs = [years]
    elseif isa(years, Array{Int, 1}); yrs = years
    end

    if mode == "ie"; ce = ieHHs
    elseif mode == "de"; ce = deHHs
    elseif mode == "cf"; ce = cfHHs
    end

    ci = findfirst(x -> x == "Total", catList)

    for y in yrs
        if length(nations) == 0; nats = natList[y]
        elseif isa(nations, String); nats = [nations]
        elseif isa(nations, Array{String, 1}); nats = nations
        end
        for n in nats
            outputFile = replace(replace(outputFileTag, "YEAR_" => string(y) * "_"), "NATION_" => n * "_")
            f = open(outputFile, "w")
            println(f, "Year\tNation\tHHID\tSize\tWeight\tPosition\tCF_per_capita\tWeight_accumulated\tCF_accumulated")
            hhl, hhe = hhsList[y][n], ce[y][n]
            nh = length(hhl)

            cfpc = [hhe[i,ci] / siz[y][hhl[i]] for i=1:nh]
            si = sort(collect(1:nh), by = x -> cfpc[x])
            wgh_sum = sum([wghNuts[y][h] * siz[y][h] for h in hhl])
            wgh_acc, cf_acc = 0, 0
            for i in si
                hh = hhl[i]
                wgh_acc += siz[y][hh] * wghNuts[y][hh]
                cf_acc += cfpc[i] * siz[y][hhl[i]] * wghNuts[y][hh]
                println(f, y, "\t", n, "\t", hh, "\t", siz[y][hh], "\t", wghNuts[y][hh], "\t", wgh_acc / wgh_sum, "\t", cfpc[i], "\t", wgh_acc, "\t", cf_acc)
            end
            close(f)
        end
    end
end

function sortHHsByStatus(years, nations, outputFileTag = ""; mode = "cf", except=[], sort_mode="cfpc")
    # sort_mode: "income", "income_pc", "cfpc"

    global yrList, hhsList, catList, natList, nat, siz, inc
    global ieHHs, deHHs, cfHHs, wghNuts

    mkpath(rsplit(outputFileTag, '/', limit = 2)[1])

    if length(years) == 0; yrs = yrList
    elseif isa(years, Int); yrs = [years]
    elseif isa(years, Array{Int, 1}); yrs = years
    end

    if mode == "ie"; ce = ieHHs
    elseif mode == "de"; ce = deHHs
    elseif mode == "cf"; ce = cfHHs
    end

    hh_pos = Dict{Int, Dict{String, Dict{String, Float64}}}()
    cat_list = filter(x -> x != "Total", catList)
    ci = findfirst(x -> x == "Total", catList)
    cis = [findfirst(x -> x == ct, catList) for ct in cat_list]

    for y in yrs
        hh_pos[y] = Dict{String, Dict{String, Float64}}()
        if length(nations) == 0; nats = natList[y]
        elseif isa(nations, String); nats = [nations]
        elseif isa(nations, Array{String, 1}); nats = nations
        end
        filter!(x -> !(x in except), nats)

        for n in nats
            hh_pos[y][n] = Dict{String, Float64}()
            hhl, hhe = hhsList[y][n], ce[y][n]
            nh = length(hhl)

            if sort_mode == "income"; si = sort(collect(1:nh), by = x -> inc[y][hhl[x]])
            elseif sort_mode == "income_pc"; si = sort(collect(1:nh), by = x -> inc[y][hhl[x]]/siz[y][hhl[x]])
            elseif sort_mode == "cfpc"; si = sort(collect(1:nh), by = x -> hhe[x,ci]/siz[y][hhl[x]])
            else println("wrong sorting_mode: ", sort_mode)
            end
            wgh_sum = sum([wghNuts[y][h] * siz[y][h] for h in hhl])
            wgh_acc, cf_acc = 0, 0

            if length(outputFileTag) > 0
                outputFile = replace(replace(outputFileTag, "YEAR_" => string(y) * "_"), "NATION_" => n * "_")
                f = open(outputFile, "w")
                print(f, "Year\tNation\tHHID\tSize\tWeight\tPosition\tIncome\tCF_Total")
                for ct in cat_list; print(f, "\tCF_", ct) end
                println(f,"\tWeight_accumulated\tCF_accumulated")
            end

            for i in si
                hh = hhl[i]
                wgh_acc += siz[y][hh] * wghNuts[y][hh]
                cf_acc += hhe[i,ci] * wghNuts[y][hh]
                hh_pos[y][n][hh] = wgh_acc / wgh_sum

                if length(outputFileTag) > 0
                    print(f, y, "\t", n, "\t", hh, "\t", siz[y][hh], "\t", wghNuts[y][hh], "\t", hh_pos[y][n][hh], "\t", inc[y][hh], "\t", hhe[i,ci])
                    for ct_i = 1:length(cat_list); print(f, "\t", hhe[i, cis[ct_i]]) end
                    println(f, "\t", wgh_acc, "\t", cf_acc)
                end
            end
            if length(outputFileTag) > 0; close(f) end
        end
    end

    return hh_pos
end

function profilingEmissionByIncome(years=[], nations=[], outputFile=""; n_top = 10, mode = "cf", period="annual", sort_mode="cfpc", adjust=true, popWgh=false, ntweigh=false, except=[])
    # mode: "ie" indirect CE, "de" direct CE, "cf" integrated CF
    # period: "annual", "monthly", or "daily"

    global yrList, hhsList, catList, natList, nat, siz, inc, sec, secName
    global indirectCE, directCE, integratedCF, wghNuts

    if length(years) == 0; yrs = yrList
    elseif isa(years, Int); yrs = [years]
    elseif isa(years, Array{Int, 1}); yrs = years
    end

    if mode == "ie"; e = indirectCE
    elseif mode == "de"; e = directCE
    elseif mode == "cf"; e = integratedCF
    end

    cat_list = filter(x -> x != "Total", catList)
    nc = length(cat_list)
    ci = findfirst(x -> x == "Total", catList)
    cis = [findfirst(x -> x == ct, catList) for ct in cat_list]

    intervals = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    n_int = length(intervals) - 1

    pos = Dict{Int, Dict{String, Array{Float64, 1}}}()
    profile = Dict{Int, Dict{String, Array{Array{Array{Tuple{String, Float64}, 1}, 1}, 1}}}()
    for y in yrs
        if length(nations) == 0; nats = natList[y]
        elseif isa(nations, String); nats = [nations]
        elseif isa(nations, Array{String, 1}); nats = nations
        end
        filter!(x -> !(x in except), nats)
        pos[y] = Dict{String, Array{Float64, 1}}()
        profile[y] = Dict{String, Array{Array{Array{Tuple{String, Float64}, 1}, 1}, 1}}()
        ns = length(sec[y])
        catidx = [if haskey(cat[y], s); findfirst(x->x==cat[y][s], catList) end for s in sec[y]]

        for n in nats
            hhl = hhsList[y][n]
            nh = length(hhl)
            pos[y][n] = zeros(Float64, nh)
            profile[y][n] = Array{Array{Array{Tuple{String, Float64}, 1}, 1}, 1}()

            if sort_mode == "income"; si = sort(collect(1:nh), by = x -> inc[y][hhl[x]])
            elseif sort_mode == "income_pc"; si = sort(collect(1:nh), by = x -> inc[y][hhl[x]]/siz[y][hhl[x]])
            elseif sort_mode == "cfpc"; si = sort(collect(1:nh), by = x -> ce[y][n][x,ci]/siz[y][hhl[x]])
            else println("wrong sorting_mode: ", sort_mode)
            end
            wgh_sum = sum([wghNuts[y][h] * siz[y][h] for h in hhl])
            wgh_acc = 0
            for i in si
                hh = hhl[i]
                wgh_acc += siz[y][hh] * wghNuts[y][hh]
                pos[y][n][i] = wgh_acc / wgh_sum
            end

            e_wgh = zeros(Float64, ns, nh)
            for i = 1:nh; e_wgh[:,i] = e[y][n][:,i] .* wghNuts[y][hhl[i]] end

            he_lv = [zeros(Float64, ns) for i=1:n_int]
            for i = 1:n_int
                push!(profile[y][n], Array{Array{Tuple{String, Float64}, 1}, 1}())
                if i < n_int; idxs = findall(x -> intervals[i] <= pos[y][n][x] < intervals[i+1], 1:nh)
                elseif i == n_int; idxs = findall(x -> intervals[i] <= pos[y][n][x], 1:nh)
                end

                he_lv[i] = vec(sum(e_wgh[:,idxs], dims = 2))

                pop_int = sum([siz[y][hhl[hi]] * wghNuts[y][hhl[hi]] for hi in idxs])

                for j = 1:nc
                    he_lv_cat = [(catidx[k] == j ? he_lv[i][k] : 0) for k = 1:ns]
                    he_lv_sum = sum(he_lv_cat)
                    ssi = sort(collect(1:ns), by = x -> he_lv_cat[x], rev = true)
                    push!(profile[y][n][i], [(sec[y][ssi[k]], he_lv_cat[ssi[k]]/he_lv_sum) for k = 1:n_top])
                end
            end
        end
    end

    f = open(outputFile, "w")
    println(f, "Year\tNation\tHH_level\tCategory\tProfile")
    for y in yrs
        if length(nations) == 0; nats = natList[y]
        elseif isa(nations, String); nats = [nations]
        elseif isa(nations, Array{String, 1}); nats = nations
        end
        for n in nats, i = 1:n_int, j = 1:nc
            print(f, y, "\t", n, "\t", i, "\t", cat_list[j])
            for k = 1:n_top
                print(f, "\t", secName[y][profile[y][n][i][j][k][1]]," (", round(profile[y][n][i][j][k][2]*100, digits=1),")")
            end
            println(f)
        end
    end
    close(f)
end

function profilingExpenditureByIncome(year, nations=[], expFile = "", outputFile=""; n_top = 10, period="annual", sort_mode="cfpc", adjust=true, popWgh=false, ntweigh=false, except=[])
    # period: "annual", "monthly", or "daily"

    global yrList, hhsList, catList, natList, nat, siz, inc, sec, secName
    global indirectCE, directCE, integratedCF, wghNuts

    y = year

    if length(nations) == 0; nats = natList[y]
    elseif isa(nations, String); nats = [nations]
    elseif isa(nations, Array{String, 1}); nats = nations
    end
    ns = length(sec[y])

    e = Dict{Int, Dict{String, Array{Float64, 2}}}() # expenditure: {year, {nation, {sector, hhid}}}
    e[y] = Dict{String, Array{Float64, 2}}()
    cnt = Dict(nats .=> 0)
    for n in nats; e[y][n] = zeros(Float64, ns, length(hhsList[y][n])) end
    f = open(expFile)
    readline(f)
    nat_prev = ""
    yr = string(year)
    for l in eachline(f)
        s = string.(split(l, ','))
        if s[1] == yr
            n = s[2]
            cnt[n] += 1
            if hhsList[y][n][cnt[n]] == n * "_" * s[3]; e[y][n][:,cnt[n]] = map(x -> parse(Float64, x), s[4:end])
            else println("HHID mismatch: ", hhsList[y][n][cnt[n]], ", ", s[2]*"_"*s[3])
            end
        end
    end
    close(f)



    cat_list = filter(x -> x != "Total", catList)
    nc = length(cat_list)
    ci = findfirst(x -> x == "Total", catList)
    cis = [findfirst(x -> x == ct, catList) for ct in cat_list]

    intervals = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    n_int = length(intervals) - 1

    pos = Dict{Int, Dict{String, Array{Float64, 1}}}()
    profile = Dict{Int, Dict{String, Array{Array{Array{Tuple{String, Float64}, 1}, 1}, 1}}}()
    filter!(x -> !(x in except), nats)
    pos[y] = Dict{String, Array{Float64, 1}}()
    profile[y] = Dict{String, Array{Array{Array{Tuple{String, Float64}, 1}, 1}, 1}}()
    ns = length(sec[y])
    catidx = [if haskey(cat[y], s); findfirst(x->x==cat[y][s], catList) end for s in sec[y]]

    for n in nats
        hhl = hhsList[y][n]
        nh = length(hhl)
        pos[y][n] = zeros(Float64, nh)
        profile[y][n] = Array{Array{Array{Tuple{String, Float64}, 1}, 1}, 1}()

        if sort_mode == "income"; si = sort(collect(1:nh), by = x -> inc[y][hhl[x]])
        elseif sort_mode == "income_pc"; si = sort(collect(1:nh), by = x -> inc[y][hhl[x]]/siz[y][hhl[x]])
        elseif sort_mode == "cfpc"; si = sort(collect(1:nh), by = x -> ce[y][n][x,ci]/siz[y][hhl[x]])
        else println("wrong sorting_mode: ", sort_mode)
        end
        wgh_sum = sum([wghNuts[y][h] * siz[y][h] for h in hhl])
        wgh_acc = 0
        for i in si
            hh = hhl[i]
            wgh_acc += siz[y][hh] * wghNuts[y][hh]
            pos[y][n][i] = wgh_acc / wgh_sum
        end

        e_wgh = zeros(Float64, ns, nh)
        for i = 1:nh; e_wgh[:,i] = e[y][n][:,i] .* wghNuts[y][hhl[i]] end

        he_lv = [zeros(Float64, ns) for i=1:n_int]
        for i = 1:n_int
            push!(profile[y][n], Array{Array{Tuple{String, Float64}, 1}, 1}())
            if i < n_int; idxs = findall(x -> intervals[i] <= pos[y][n][x] < intervals[i+1], 1:nh)
            elseif i == n_int; idxs = findall(x -> intervals[i] <= pos[y][n][x], 1:nh)
            end

            he_lv[i] = vec(sum(e_wgh[:,idxs], dims = 2))

            pop_int = sum([siz[y][hhl[hi]] * wghNuts[y][hhl[hi]] for hi in idxs])

            for j = 1:nc
                he_lv_cat = [(catidx[k] == j ? he_lv[i][k] : 0) for k = 1:ns]
                he_lv_sum = sum(he_lv_cat)
                ssi = sort(collect(1:ns), by = x -> he_lv_cat[x], rev = true)
                push!(profile[y][n][i], [(sec[y][ssi[k]], he_lv_cat[ssi[k]]/he_lv_sum) for k = 1:n_top])
            end
        end
    end

    f = open(outputFile, "w")
    println(f, "Year\tNation\tHH_level\tCategory\tProfile")
    for n in nats, i = 1:n_int, j = 1:nc
        print(f, y, "\t", n, "\t", i, "\t", cat_list[j])
        for k = 1:n_top
            print(f, "\t", secName[y][profile[y][n][i][j][k][1]]," (", round(profile[y][n][i][j][k][2]*100, digits=1),")")
        end
        println(f)
    end
    close(f)
end

end
