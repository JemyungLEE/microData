module EmissionCategorizer

# Developed date: 17. May. 2021
# Last modified date: 25. May. 2021
# Subject: Categorize households' carbon footprints
# Description: Read household-level indirect and direct carbon emissions,  integrate them to be CF,
#              and categorize the CFs by consumption category, district, expenditure-level, and etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

include("MicroDataReader.jl")

using XLSX
using Statistics
using Formatting: printfmt
using .MicroDataReader
mdr = MicroDataReader

yr_list = Array{Int, 1}()       # year list: {YYYY}
nat_list = Array{String, 1}()   # nation list: {A3}
cat_list = Array{String, 1}()   # category list
rel_list = Array{String, 1}()   # religion list

hh_list = Dict{Int, Dict{String, Array{String, 1}}}()           # Household ID: {year, {nation, {hhid}}}
sc_list = Dict{Int, Dict{String, Array{String, 1}}}()           # commodity code list: {year, {nation A3, {code}}}
households = Dict{Int, Dict{String, Dict{String, mdr.household}}}() # household dict: {year, {nation A3, {hhid, household}}}
sectors = Dict{Int, Dict{String, Dict{String, mdr.commodity}}}()    # expenditure sector: {year, {nation A3, {code, commodity (of mdr)}}}
sc_cat = Dict{Int, Dict{String, Dict{String, String}}}()        # CES/HBS sector-category link dict: {year, {nation, {sector_code, category}}}

regions = Dict{Int, Dict{String, Dict{String, String}}}()       # region code-name: {year, {nation A3, {code, region}}}
prov_list = Dict{Int, Dict{String, Array{String, 1}}}()         # province code list: {year, {nation A3, {code}}}
dist_list = Dict{Int, Dict{String, Array{String, 1}}}()         # district code list: {year, {nation A3, {code}}}
dist_prov = Dict{Int, Dict{String, Dict{String, String}}}()     # district's province: {year, {nation A3, {district code, province code}}}

pops = Dict{Int, Dict{String, Dict{String, Float64}}}()         # population: {year, {nation, {region_code, population}}}
pop_wgh = Dict{Int, Dict{String, Dict{String, Float64}}}()      # population weight: {year, {nation, {region_code, weight}}}
pops_ur = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}()      # urban/rural population: {year, {nation, {region_code, (urban, rural)}}
pop_ur_wgh = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}()   # urban/rural population weight: {year, {nation, {region_code, (urban, rural)}}

hh_curr = Dict{Int, Dict{String, Array{String, 1}}}()            # currency unit for household values (income or expenditure): {year, {nation, {currency}}}
hh_period = Dict{Int, Dict{String, Array{String, 1}}}()          # period for household values (income or expenditure): {year, {nation, {period}}}
exp_curr = Dict{Int, Dict{String, Array{String, 1}}}()           # currency unit for expenditure values: {year, {nation, {currency}}}
exp_period = Dict{Int, Dict{String, Array{String, 1}}}()         # period for expenditure values: {year, {nation, {period}}}

directCE = Dict{Int, Dict{String, Array{Float64, 2}}}()         # direct carbon emission: {year, {nation, {CES/HBS sector, household}}}
indirectCE = Dict{Int, Dict{String, Array{Float64, 2}}}()       # indirect carbon emission: {year, {nation, {CES/HBS sector, household}}}
integratedCF = Dict{Int, Dict{String, Array{Float64, 2}}}()     # carbon footprint = DE + IE: {year, {nation, {CES/HBS sector, household}}}

ieHHs = Dict{Int, Dict{String, Array{Float64, 2}}}()        # categozied indirect emission by household: {year, {nation, {hhid, category}}}
deHHs = Dict{Int, Dict{String, Array{Float64, 2}}}()        # categozied direct emission by household: {year, {nation, {hhid, category}}}
cfHHs = Dict{Int, Dict{String, Array{Float64, 2}}}()        # categozied carbon footprint by household: {year, {nation, {hhid, category}}}
ieReg = Dict{Int, Dict{String, Array{Float64, 2}}}()        # categozied indirect emission by region: {year, {nation, {region, category}}}
deReg = Dict{Int, Dict{String, Array{Float64, 2}}}()        # categozied direct emission by region: {year, {nation, {region, category}}}
cfReg = Dict{Int, Dict{String, Array{Float64, 2}}}()        # categozied carbon footprint by region: {year, {nation, {region, category}}}
ieRegDev = Dict{Int, Dict{String, Array{Float64, 2}}}()     # categozied indirect emission deviation from mean by region: {year, {nation, {region, category}}}
deRegDev = Dict{Int, Dict{String, Array{Float64, 2}}}()     # categozied direct emission deviation from mean by region: {year, {nation, {region, category}}}
cfRegDev = Dict{Int, Dict{String, Array{Float64, 2}}}()     # categozied carbon footprint deviation from mean by region: {year, {nation, {region, category}}}

reg_sample = Dict{Int, Dict{String, Dict{String, Tuple{Int,Int}}}}()    # sample population and households by districct: {year, {nation, {district code, (sample population, number of households)}}}
reg_avgExp = Dict{Int, Dict{String, Dict{String, Float64}}}()           # average annual expenditure per capita, USD/yr: {year, {nation, {district code, mean Avg.Exp./cap/yr}}}
reg_sample_ur = Dict{Int, Dict{String, Dict{String, Array{Tuple{Int,Int}, 1}}}}()   # sample population and households by districct: {year, {nation, {district code, {(sample population, number of households): [urban, rural]}}}}
reg_avgExp_ur = Dict{Int, Dict{String, Dict{String, Array{Float64, 1}}}}()  # average annual expenditure per capita, USD/yr: {year, {nation, {district code, {mean Avg.Exp./cap/yr: [urban, rural]}}}}

# GIS data
majorCity = Dict{Int, Dict{String, String}}()                   # major city in the region: {year, {Upper_region_code, major_city_code}}
gisCoord = Dict{Int, Dict{String, Tuple{Float64, Float64}}}()   # GIS coordinates: {year, {region_code, {X_longitude, Y_latitude}}}
gisCatLab = Dict{String, String}()              # Web-exporting category matching table: {Category label in program, in web-site files}

gisRegList = Dict{Int, Array{String, 1}}()      # GIS region list: {year, {region code}}
gisPop = Dict{Int, Array{Float64, 1}}()         # GIS version, total population by region
gisSample = Dict{Int, Array{Float64, 1}}()      # GIS version, total sample members by region
gisAvgExp = Dict{Int, Array{Float64, 1}}()      # GIS version, average expenditure by region

gisRegIe = Dict{Int, Array{Float64, 2}}()       # categozied indirect emission by region: {year, {category, region(hbscd)}}
gisRegDe = Dict{Int, Array{Float64, 2}}()       # categozied direct emission by region: {year, {DE category, region(hbscd)}}
gisRegCf = Dict{Int, Array{Float64, 2}}()       # categozied carbon footprint by region: {year, {category, region(hbscd)}}
gisRegIeRank = Dict{Int, Array{Int, 2}}()       # categozied indirect emission rank by region: {year, {category, region(hbscd)}}
gisRegDeRank = Dict{Int, Array{Int, 2}}()       # categozied direct emission rank by region: {year, {DE category, region(hbscd)}}
gisRegCfRank = Dict{Int, Array{Int, 2}}()       # categozied carbon footprint rank by region: {year, {category, region(hbscd)}}
gisRegIePerCap = Dict{Int, Array{Float64, 2}}() # categozied indirect emission per capita by region: {year, {category, region(hbscd)}}
gisRegDePerCap = Dict{Int, Array{Float64, 2}}() # categozied direct emission per capita by region: {year, {DE category, region(hbscd)}}
gisRegCfPerCap = Dict{Int, Array{Float64, 2}}() # categozied carbon footprint per capita by region: {year, {category, region(hbscd)}}
gisRegIeRankPerCap = Dict{Int, Array{Int, 2}}() # categozied indirect emission per capita rank by region: {year, {category, region(hbscd)}}
gisRegDeRankPerCap = Dict{Int, Array{Int, 2}}() # categozied direct emission per capita rank by region: {year, {DE category, region(hbscd)}}
gisRegCfRankPerCap = Dict{Int, Array{Int, 2}}() # categozied carbon footprint per capita rank by region: {year, {category, region(hbscd)}}
gisRegCfDiff = Dict{Int, Array{Float64, 2}}()   # differences of categozied carbon footprint by region: (emission-mean)/mean, {year, {category, district(GID)}}
gisRegCfDiffRank = Dict{Int, Array{Int, 2}}()   # difference ranks of categozied carbon footprint by region: (emission-mean)/mean, {year, {category, district(GID)}}
gisRegCfDiffPerCap = Dict{Int, Array{Float64, 2}}() # differences of categozied carbon footprint per capita by region: (emission-mean)/mean, {year, {category, district(GID)}}
gisRegCfDiffRankPerCap = Dict{Int, Array{Int, 2}}() # difference ranks of categozied carbon footprint per capita by region: (emission-mean)/mean, {year, {category, district(GID)}}

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

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

function importMicroData(mdata::Module)

    global yr_list, nat_list
    global hh_list = mdata.hh_list
    global sc_list = mdata.sc_list
    global households = mdata.households
    global sectors = mdata.sectors

    global regions = mdata.regions
    global prov_list = mdata.prov_list
    global dist_list = mdata.dist_list
    global dist_prov = mdata.dist_prov

    global pops = mdata.pops
    global pop_wgh = mdata.pop_wgh
    global pops_ur = mdata.pops_ur
    global pop_ur_wgh = mdata.pop_ur_wgh

    global hh_curr = mdata.hh_curr
    global hh_period = mdata.hh_period
    global exp_curr = mdata.exp_curr
    global exp_period = mdata.exp_period

    yrs = sort(collect(keys(hh_list)))
    for y in yrs
        if !(y in yr_list); push!(yr_list, y) end
        nats = sort(collect(keys(hh_list[y])))
        for n in nats; if !(n in nat_list); push!(nat_list, n) end end
    end
    sort!(yr_list)
    sort!(nat_list)
end

function importEmissionData(edata::Module)

    global yr_list, nat_list
    global directCE = edata.directCE
    global indirectCE = edata.indirectCE

    for ce in [directCE, indirectCE]
        yrs = sort(collect(keys(ce)))
        neg_yrs = filter(x->!(x in yrs), yr_list)
        abs_yrs = filter(x->!(x in yr_list), yrs)
        if length(neg_yrs) > 0; println("Emission data does not contain $neg_yrs year data") end
        if length(abs_yrs) > 0; println("Micro-data does not contain $abs_yrs year data") end
        for y in yrs
            nts = sort(collect(keys(ce[y])))
            neg_nat = filter(x->!(x in nts), nat_list)
            abs_nat = filter(x->!(x in nat_list), nts)
            if length(neg_nat) > 0; println("Emission data does not contain $neg_nat nation data in $y year") end
            if length(abs_nat) > 0; println("Micro-data does not contain $abs_nat nation data in $y year") end
        end
    end
end

function readEmissionData(year, nation, inputFile; mode = "ie")

    global hh_list, sc_list, directCE, indirectCE
    hl, sl = hh_list[year][nation], sc_list[year][nation]
    ns, nh = length(sl), length(hl)

    em = zeros(Float64, ns, nh)

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    hhs = string.(strip.(split(readline(f), f_sep)[2:end]))

    if hhs == hl; hi = 1:nh
    elseif sort(hhs) == sort(hl); hi = [findfirst(x->x==hh, hl) for hh in hhs]
    else println("Error: Emission matrix's household list and micro-data's household list are not same.")
    end

    scs = Array{String, 1}()
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        push!(scs, s[1])
        si = findfirst(x->x==s[1], sl)
        em[si,hi] = [parse(Float64, x) for x in s[2:end]]
    end
    close(f)

    if mode == "ie"
        if !haskey(indirectCE, year); indirectCE[year] = Dict{String, Array{Float64, 2}}() end
        indirectCE[year][nation] = em
    elseif mode == "de"
        if !haskey(directCE, year); directCE[year] = Dict{String, Array{Float64, 2}}() end
        directCE[year][nation] = em
    else println("Wrong emission print mode: $mode")
    end
end

function integrateCarbonFootprint(; years=[], nations=[], mode="cf")

    global yr_list, nat_list, hh_list, sc_list, directCE, indirectCE, integratedCF

    if length(years) > 0; yrs = years else yrs = yr_list end
    if length(nations) > 0; nts = nations else nts = nat_list end

    for y in yrs; if !haskey(integratedCF, y); integratedCF[y] = Dict{String, Array{Float64, 2}}() end end
    for y in yrs, n in nts
        integratedCF[y][n] = zeros(Float64, length(sc_list[y][n]), length(hh_list[y][n]))
        if lowercase(mode) in ["cf", "ie"]; integratedCF[y][n] += indirectCE[y][n] end
        if lowercase(mode) in ["cf", "de"]; integratedCF[y][n] += directCE[y][n] end
    end
end

function setCategory(year, nation; subgroup = "", except::Array{String, 1}=[])  # Note: sub-grouping parts should be added

    global sectors, sc_list, cat_list, sc_cat
    ss, sl = sectors[year][nation], sc_list[year][nation]
    if !haskey(sc_cat, year); sc_cat[year] = Dict{String, Dict{String, String}}() end
    if !haskey(sc_cat[year], nation); sc_cat[year][nation] = Dict{String, String}() end
    cat_list = Array{String, 1}()
    sc_ct = sc_cat[year][nation]

    if length(subgroup) == 0
        for sc in sl
            c = ss[sc].category
            if !(c in cat_list); push!(cat_list, c) end
            sc_ct[sc] = c
        end
        cat_list = sort(filter(x->!(x in except), cat_list))
    end
end

function categorizeHouseholdEmission(years=[], nations=[]; mode="cf", output="", hhsinfo=false)

    global yr_list, nat_list, hh_list, sc_list, sc_cat, cat_list, households
    global directCE, indirectCE, integratedCF, ieHHs, deHHs, cfHHs

    agg_label = "Total"
    if !(agg_label in cat_list); push!(cat_list, agg_label) end
    nc = length(cat_list)
    if isa(years, Number); years = [years] end
    if isa(nations, String); nations = [nations] end

    if lowercase(mode) == "ie"; em = indirectCE; ce = ieHHs
    elseif lowercase(mode) == "de"; em = directCE; ce = deHHs
    elseif lowercase(mode) == "cf"; em = integratedCF; ce = cfHHs
    else println("wrong emission categorizing mode: ", mode)
    end

    for y in yr_list, n in nat_list
        hl, sl, scct = hh_list[y][n], sc_list[y][n], sc_cat[y][n]
        nh, ns = length(hl), length(sl)
        eh = em[y][n]
        ec = zeros(Float64, nh, nc)
        for i = 1:nc-1
            si = [findfirst(x->x == s,  sl) for s in filter(x -> scct[x] == cat_list[i], sl)]
            ec[:,i] = sum(eh[si,:], dims=1)
            ec[:,nc] += ec[:,i]
        end
        if !haskey(ce, y); ce[y] = Dict{String, Array{Float64, 2}}() end
        ce[y][n] = ec
    end

    # print the results
    if length(output)>0
        f_sep = getValueSeparator(output)
        f = open(output, "w")
        print(f, "Year", f_sep, "Nation", f_sep, "HHID"); for c in cat_list; print(f, f_sep, c) end
        if hhsinfo; print(f, f_sep, "Size", f_sep, "Income", f_sep, "Expend", f_sep, "PopWgh") end; println(f)
        for y in years, n in nat_list
            hl, hhs = hh_list[y][n], households[y][n]
            for i = 1:length(hl)
                hh = hl[i]
                print(f, y, f_sep, n, f_sep, hh)
                for j = 1:nc; print(f, f_sep, ce[y][n][i,j]) end
                if hhsinfo; print(f, f_sep, hhs[hh].size, f_sep, hhs[hh].totinc, f_sep, hhs[hh].totexp, f_sep, hhs[hh].popwgh) end
                println(f)
            end
        end
        close(f)
    end
end

function categorizeRegionalEmission(years=[], nations=[]; mode = "cf", period="year", popwgh=false, region = "district", ur = false, religion=false)
    # mode: "ie" indirect CE, "de" direct CE, "cf" integrated CF
    # period: "day", "week", "month", or "year"
    # region: "province" or "district"
    # religion: [true] categorize districts' features by religions

    global yr_list, nat_list, hh_list, sc_list, sc_cat, cat_list, rel_list, prov_list, dist_list
    global households, pops, pops_ur, hh_period, reg_sample, reg_avgExp
    global ieHHs, deHHs, cfHHs, ieReg, deReg, cfReg, ieRegDev, deRegDev, cfRegDev
    pr_unts = Dict("day" => 1, "week" => 7,"month" => 30, "year" => 365)

    em_tabs = Dict{Int, Dict{String, Array{Any, 1}}}()  # eReg, eRegDev, eReg_ur, eRegDev_ur, eReg_rel, eRegDev_rel, eReg_rel_ur, eRegDev_rel_ur

    if religion
        rsam_rel = Dict{Int, Dict{String, Dict{String, Array{Tuple{Int,Int}, 1}}}}()
        ravg_rel = Dict{Int, Dict{String, Dict{String, Array{Float64, 1}}}}()
        if ur
            rsam_rel_ur = Dict{Int, Dict{String, Dict{String, Array{ Array{Tuple{Int,Int}, 1}, 1}}}}()
            ravg_rel_ur = Dict{Int, Dict{String, Dict{String, Array{ Array{Float64, 1}, 1}}}}()
        end
    end

    if ur; reg_type = ["urban", "rural"]; n_rtyp = length(reg_type) end   # region type: "urban"/"rural" or "dense"/"intermediate"/"sparce"
    nc = length(cat_list)
    if isa(years, Number); years = [years] end
    if isa(nations, String); nations = [nations] end

    if lowercase(mode) == "ie"; eh = ieHHs; ereg = ieReg; eregdev = ieRegDev
    elseif lowercase(mode) == "de"; eh = deHHs; ereg = deReg; eregdev = deRegDev
    elseif lowercase(mode) == "cf"; eh = cfHHs; ereg = cfReg; eregdev = cfRegDev
    else println("wrong emission categorizing mode: ", mode)
    end

    if region == "district"; reg_list = dist_list
    elseif region == "province"; reg_list = prov_list
    else println("Wrong region mode: $region")
    end

    for y in years, n in nations
        if !haskey(em_tabs, y); em_tabs[y] = Dict{String, Array{Any, 1}}() end
        if !haskey(reg_sample, y); reg_sample[y] = Dict{String, Dict{String, Tuple{Int,Int}}}() end
        if !haskey(reg_avgExp, y); reg_avgExp[y] = Dict{String, Dict{String, Float64}}() end
        if ur && !haskey(reg_sample_ur, y); reg_sample_ur[y] = Dict{String, Dict{String, Array{Tuple{Int,Int}, 1}}}() end
        if ur && !haskey(reg_avgExp_ur, y); reg_avgExp_ur[y] = Dict{String, Dict{String, Array{Float64, 1}}}() end
        if religion; rsam_rel[y] = Dict{String, Dict{String, Array{Tuple{Int,Int}, 1}}}() end
        if religion; ravg_rel[y] = Dict{String, Dict{String, Array{Float64, 1}}}() end
        if ur && religion; rsam_rel_ur[y] = Dict{String, Dict{String, Array{ Array{Tuple{Int,Int}, 1}, 1}}}() end
        if ur && religion; ravg_rel_ur[y] = Dict{String, Dict{String, Array{ Array{Float64, 1}, 1}}}() end

        hhs, cl, hl, rl= households[y][n], cat_list, hh_list[y][n], reg_list[y][n]
        nc, nh, nr = length(cl), length(hl), length(rl)

        rsam = reg_sample[y][n] = Dict{String, Tuple{Int,Int}}()
        ravg = reg_avgExp[y][n] = Dict{String, Float64}()
        if ur; rsam_ur = reg_sample_ur[y][n] = Dict{String, Array{Tuple{Int,Int}, 1}}() end
        if ur; ravg_ur = reg_avgExp_ur[y][n] = Dict{String, Array{Float64, 1}}() end
        if religion; rsam_rel[y][n] = Dict{String, Array{Tuple{Int,Int}, 1}}() end
        if religion; ravg_rel[y][n] = Dict{String, Array{Float64, 1}}() end
        if ur && religion; rsam_rel_ur[y][n] = Dict{String, Array{ Array{Tuple{Int,Int}, 1}, 1}}() end
        if ur && religion; ravg_rel_ur[y][n] = Dict{String, Array{ Array{Float64, 1}, 1}}() end

        # make region index arrays
        if region == "district"; regidx = [filter(i->hhs[hl[i]].district == d, 1:nh) for d in rl]
        elseif region == "province"; regidx = [filter(i->hhs[hl[i]].province == p, 1:nh) for p in rl]
        end
        if ur; regidx_ur = [[filter(i->hhs[hl[i]].regtype == r_typ, ri) for r_typ in reg_type] for ri in regidx] end
        if religion; relidx = [[filter(i->hhs[hl[i]].rel == rlg, ri) for rlg in rel_list] for ri in regidx] end
        if ur && religion; relidx_ur = [[[filter(i->hhs[hl[i]].regtype == r_typ, idx) for r_typ in reg_type] for idx in r_idx] for r_idx in relidx] end

        # sum sample households and members by regions
        thbr = [length(idxs) for idxs in regidx]                            # total households by region
        tpbr = [sum([hhs[hl[i]].size for i in idxs]) for idxs in regidx]    # total sample population by region
        for i = 1:nr; rsam[rl[i]] = (tpbr[i], thbr[i]) end
        if ur
            thbr_ur = [[length(idxs) for idxs in idxs_ur] for idxs_ur in regidx_ur]
            tpbr_ur = [[sum([hhs[hl[i]].size for i in idxs]) for idxs in idxs_ur] for idxs_ur in regidx_ur]
            for i = 1:nr; rsam_ur[rl[i]] = collect(zip(thbr_ur[i], tpbr_ur[i])) end
        end

        # sum sample households and members by regions and religions
        if religion
            n_rel = length(rel_list)
            thbrr = [[length(idxs) for idxs in idxs_rel] for idxs_rel in relidx]    # total households by religion and region
            tpbrr = [[sum([hhs[hl[i]].size for i in idxs]) for idxs in idxs_rel] for idxs_rel in relidx]    # total members of households by religion and region
            for i = 1:nr; rsam_rel[y][n][rl[i]] = collect(zip(thbrr[i], tpbrr[i])) end
            if ur
                thbrr_ur = [[[length(idxs) for idxs in idxs_ur] for idxs_ur in idxs_rel] for idxs_rel in relidx_ur]
                tpbrr_ur = [[[sum([hhs[hl[i]].size for i in idxs]) for idxs in idxs_ur] for idxs_ur in idxs_rel] for idxs_rel in relidx_ur]
                for i = 1:nr; rsam_rel_ur[y][n][rl[i]] = [collect(zip(thbrr_ur[i][j], tpbrr_ur[i][j])) for j = 1:n_rel] end
            end
        end

        # calculate average expenditure per capita by region
        if popwgh
            totexp = [sum([hhs[hl[i]].totexp * hhs[hl[i]].popwgh for i in idxs]) for idxs in regidx]
            pw = [sum([hhs[hl[i]].popwgh for i in idxs]) for idxs in regidx]
            for i=1:nr; ravg[rl[i]] = totexp[i]/pw[i] end
        else
            totexp = [sum([hhs[hl[i]].totexp for i in idxs]) for idxs in regidx]
            for i=1:nr; ravg[rl[i]] = totexp[i]/tpbr[i] end
        end
        if ur && popwgh
            totexp_ur = [[sum([hhs[hl[i]].totexp * hhs[hl[i]].popwgh for i in idxs_ur]) for idxs_ur in r_idxs] for r_idxs in regidx_ur]
            pw_ur = [[sum([hhs[hl[i]].popwgh for i in idxs_ur]) for idxs_ur in r_idxs] for r_idxs in regidx_ur]
            for i=1:nr; ravg_ur[rl[i]] = [totexp_ur[i][j]/pw_ur[i][j] for j = 1:n_rtyp] end
        elseif ur && !popwgh
            totexp_ur = [[sum([hhs[hl[i]].totexp for i in idxs]) for idxs in idxs_ur] for idxs_ur in regidx_ur]
            for i=1:nr; ravg_ur[rl[i]] = [totexp_ur[i][j]/tpbr_ur[i][j] for j = 1:n_rtyp] end
        end
        if religion && popwgh
            totexp_rel = [[sum([hhs[hl[i]].totexp * hhs[hl[i]].popwgh for i in rel_idxs]) for rel_idxs in r_idxs] for r_idxs in relidx]
            pw_rel = [[sum([hhs[hl[i]].popwgh for i in rel_idxs]) for rel_idxs in r_idxs] for r_idxs in relidx]
            for i=1:nr; ravg_rel[y][n][rl[i]] = [totexp_rel[i][j]/pw_rel[i][j] for j = 1:n_rel] end
        elseif religion && !popwgh
            totexp_rel = [[sum([hhs[hl[i]].totexp for i in rel_idxs]) for rel_idxs in r_idxs] for r_idxs in relidx]
            for i=1:nr; ravg_rel[y][n][rl[i]] = [totexp_rel[i][j]/tpbrr[i][j] for j = 1:n_rel] end
        end
        if ur && religion && popwgh
            totexp_rel_ur = [[[sum([hhs[hl[i]].totexp * hhs[hl[i]].popwgh for i in ur_idx]) for ur_idx in rel_idxs] for rel_idxs in r_idxs] for r_idxs in relidx_ur]
            pw_rel_ur = [[[sum([hhs[hl[i]].popwgh for i in ur_idx]) for ur_idx in rel_idxs] for rel_idxs in r_idxs] for r_idxs in relidx_ur]
            for i=1:nr; ravg_rel_ur[y][n][rl[i]] = [[totexp_rel_ur[i][j][k]/pw_rel_ur[i][j][k] for k = 1:n_rtyp] for j = 1:n_rel] end
        elseif ur && religion && !popwgh
            totexp_rel_ur = [[[sum([hhs[hl[i]].totexp for i in ur_idx]) for ur_idx in rel_idxs] for rel_idxs in r_idxs] for r_idxs in relidx_ur]
            for i=1:nr; ravg_rel_ur[y][n][rl[i]] = [[totexp_rel_ur[i][j][k]/tpbrr_ur[i][j][k] for k = 1:n_rtyp] for j = 1:n_rel] end
        end

        # convert periods of average expenditure values
        if period != hh_period[y][n][1]
            pr_ex_r = pr_unts[period] / pr_unts[hh_period[y][n][1]]
            for r in rl
                ravg[r] *= pr_ex_r
                if ur; for i = 1:n_rtyp; ravg_ur[r][i] *= pr_ex_r end end
                if religion; for i = 1:n_rel; ravg_rel[r][i] *= pr_ex_r end end
                if ur && religion; for i = 1:n_rel, j = 1:n_rtyp; ravg_rel_ur[r][i][j] *= pr_ex_r end end
            end
        end

        # categorize emission data
        if mode == "ie"; ec = ieHHs[y][n]   # [nh, nc]
        elseif mode == "de"; ec = deHHs[y][n]
        elseif mode == "cf"; ec = cfHHs[y][n]
        end
        erc = zeros(Float64, nr, nc)
        if popwgh; for i = 1:nr; erc[i,:] = sum([ec[hi,:] * hhs[hl[hi]].popwgh for hi in regidx[i]]) end
        else for i = 1:nr; erc[i,:] = sum(ec[regidx[i],:], dims=1) end
        end
        if ur
            erc_ur = [zeros(Float64, nr, nc) for i=1:n_rtyp]
            if popwgh; for i = 1:nr, j = 1:n_rtyp; erc_ur[j][i,:] = sum([ec[hi[j],:] * hhs[hl[hi[j]]].popwgh for hi in regidx_ur[i]]) end
            elseif !popwgh; for i = 1:nr, j = 1:n_rtyp; erc_ur[j][i,:] = sum(ec[regidx_ur[i][j],:], dims=1) end
            end
        end
        if religion
            erc_rel = [zeros(Float64, nr, nc) for i=1:n_rel]
            if popwgh; for i = 1:nr, j = 1:n_rel; erc_rel[j][i,:] = sum([ec[hi[j],:] * hhs[hl[hi[j]]].popwgh for hi in relidx[i]]) end
            elseif !popwgh; for i = 1:nr, j = 1:n_rel; erc_rel[j][i,:] = sum(ec[relidx[i][j],:], dims=1) end
            end
        end
        if ur && religion
            erc_rel_ur = [[zeros(Float64, nr, nc) for i=1:n_rtyp] for j=1:n_rel]
            if popwgh; for i=1:nr, j=1:n_rel, k=1:n_rtyp; erc_rel_ur[j][k][i,:] = sum([ec[hi[j][k],:] * hhs[hl[hi[j][k]]].popwgh for hi in relidx_ur[i]]) end
            elseif !popwgh; for i=1:nr, j=1:n_rel, k=1:n_rtyp; erc_rel_ur[j][k][i,:] = sum(ec[relidx_ur[i][j][k],:], dims=1) end
            end
        end
        # normalizing
        if popwgh; for i=1:nr, j=1:nc; erc[i,j] /= pw[i] end
        else for i=1:nc; erc[:,i] ./= tpbr end
        end
        if ur && popwgh; for i=1:nr, j=1:nc, k=1:n_rtyp; erc_ur[k][i,j] /= pw_ur[i][k] end
        elseif ur && !popwgh; for i=1:nc, j=1:n_rtyp; erc_ur[j][:,i] ./= tpbr_ur[:][j] end
        end
        if religion && popwgh; for i=1:nr, j=1:nc, k=1:n_rel; erc_rel[k][i,j] /= pw_rel[i][k] end
        elseif religion && !popwgh; for i=1:nc, j=1:n_rel; erc_rel[j][:,i] ./= tpbrr[:][j] end
        end
        if ur && religion && popwgh; for i=1:nr, j=1:nc, k=1:n_rel, l=n_rtyp; erc_rel_ur[k][l][i,j] /= pw_rel_ur[i][k][l] end
        elseif ur && religion && !popwgh; for i=1:nc, j=1:n_rel, k=n_rtyp; erc_rel_ur[j][k][:,i] ./= tpbrr_ur[:][j][k] end
        end
        # calculate differences
        avg = mean(erc, dims=1)
        edrc = zeros(Float64, size(erc))
        for i=1:nc; edrc[:,i] = (erc[:,i].-avg[i])/avg[i] end
        if ur
            avg_ur = [mean(erc_ur[i], dims=1) for i = 1:n_rtyp]
            edrc_ur = [zeros(Float64, size(erc_ur[i])) for i = 1:n_rtyp]
            for i=1:n_rtyp, j=1:nc; edrc_ur[i][:,j] = (erc_ur[i][:,j].-avg_ur[i][j])/avg_ur[i][j] end
        end
        if religion
            avg_rel = [mean(erc_rel[i], dims=1) for i = 1:n_rel]
            edrc_rel = [zeros(Float64, size(erc_rel[i])) for i = 1:n_rel]
            for i=1:n_rel, j=1:nc; edrc_rel[i][:,j] = (erc_rel[i][:,j].-avg_rel[i][j])/avg_rel[i][j] end
        end
        if ur && religion
            avg_rel_ur = [[mean(erc_rel_ur[i][j], dims=1) for j = 1:n_rtyp] for i = 1:n_rel]
            edrc_rel_ur = [[zeros(Float64, size(erc_rel_ur[i][j])) for j = 1:n_rtyp] for i = 1:n_rel]
            for i=1:n_rel, j=1:n_rtyp, k=1:nc; edrc_rel_ur[i][j][:,k] = (erc_rel_ur[i][j][:,k].-avg_rel_ur[i][j][k])/avg_rel_ur[i][j][k] end
        end
        # save the results
        if !haskey(ereg, y); ereg[y] = Dict{String, Array{Float64, 2}}() end
        if !haskey(eregdev, y); eregdev[y] = Dict{String, Array{Float64, 2}}() end
        ereg[y][n] = erc
        eregdev[y][n] = edrc

        # for returning tables
        em_tabs[y][n] = Array{Any, 1}(undef, 8)
        em_tabs[y][n][1:2] = [erc, edrc]
        if ur; em_tabs[y][n][3:4] = [erc_ur, edrc_ur] end
        if religion; em_tabs[y][n][5:6] = [erc_rel, edrc_rel] end
        if ur && religion; em_tabs[y][n][7:8] = [erc_rel_ur, edrc_rel_ur] end
    end
end

function printRegionalEmission(years=[], nations=[], outputFile=""; region = "district", mode=["cf","de","ie"], popwgh=false, ur=false, religion=false)

    global yr_list, nat_list, sc_list, sc_cat, cat_list, regions, rel_list, prov_list, dist_list, dist_prov
    global pops, pops_ur, pop_wgh, pop_ur_wgh, reg_sample, reg_avgExp, reg_avgExp_ur
    global ieReg, deReg, cfReg, ieRegDev, deRegDev, cfRegDev
    if isa(years, Number); years = [years] end
    if isa(nations, String); nations = [nations] end

    ce_chk = [x in mode for x in ["cf","de","ie"]]
    items = ["Pr_code", "Province", "Ds_code", "District"]
    items = [items; "Pop"]; if ur; items = [items; ["Pop_urban", "Pop_rural"]] end
    items = [items; "Exp"]; if ur; items = [items; ["Exp_urban", "Exp_rural"]] end
    if popwgh; items = [items; "PopWgh"]; if ur; items = [items; ["PopWgh_urban", "PopWgh_rural"]] end end
    for m in mode; items = [items; [uppercase(m)*"_"*c for c in cat_list]] end

    if region == "district"; reg_list = dist_list
    elseif region == "province"; reg_list = prov_list
    else println("Wrong region mode: $region")
    end

    nc = length(cat_list)

    f_sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")

    print(f, "Year",f_sep,"Nation")
    for it in items; print(f, f_sep, it) end
    println(f)

    nc = length(cat_list)
    for y in years, n in nations
        regs, rl, pr, ds, dp = regions[y][n], reg_list[y][n], prov_list[y][n], dist_list[y][n], dist_prov[y][n]

        for i = 1:length(rl)
            r = rl[i]
            print(f, y, f_sep, n)
            if r in pr; print(f, f_sep, r, f_sep, regs[r], f_sep, f_sep)
            elseif r in ds; print(f, f_sep, dp[r], f_sep, regs[dp[r]], f_sep, r, f_sep, regs[r])
            end
            print(f, f_sep, pops[y][n][r])
            if ur; print(f, f_sep, pops_ur[y][n][r][1], f_sep, pops_ur[y][n][r][2]) end
            print(f, f_sep, reg_avgExp[y][n][r])
            if ur; print(f, f_sep, reg_avgExp_ur[y][n][r][1], f_sep, reg_avgExp_ur[y][n][r][2]) end
            if popwgh;
                print(f, f_sep, pop_wgh[y][n][r])
                if ur; print(f, f_sep, pop_ur_wgh[y][n][r][1], f_sep, pop_ur_wgh[y][n][r][2]) end
            end
            if ce_chk[1]; for j = 1:nc; print(f, f_sep, cfReg[y][n][i,j]) end end
            if ce_chk[2]; for j = 1:nc; print(f, f_sep, deReg[y][n][i,j]) end end
            if ce_chk[3]; for j = 1:nc; print(f, f_sep, ieReg[y][n][i,j]) end end
            println(f)
        end
    end

    close(f)
end

function makeNationalSummary(year, outputFile)

    global hh_list, natList, natName, siz, wgh
    global indirectCE, directCE, integratedCF

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
        for j = 1:length(hh_list[year][n])
            h = hh_list[year][n][j]
            ie = sum(indirectCE[year][n][:,j])
            de = sum(directCE[year][n][:,j])
            cf = sum(integratedCF[year][n][:,j])
            natsam[i] += siz[year][h]
            nateqs[i] += eqs[year][h]
            natmeqs[i] += meqs[year][h]
            natwgh[i] += wgh[year][h] * siz[year][h]
            natcf[i] += wgh[year][h] * cf
            natcfpc[i] += cf
            natcfph[i] += cf
            natcfpeqs[i] += cf
            natcfpmeqs[i] += cf
            natie[i] += wgh[year][h] * ie
            natiepc[i] += ie
            natde[i] += wgh[year][h] * de
            natdepc[i] += de
        end
        natcfpc[i] /= natsam[i]
        natcfph[i] /= length(hh_list[year][n])
        natcfpeqs[i] /= nateqs[i]
        natcfpmeqs[i] /= natmeqs[i]
        natiepc[i] /= natsam[i]
        natdepc[i] /= natsam[i]
    end

    f = open(outputFile, "w")
    println(f, "Nation\tHHs\tMMs\tWeights\tCF_overall\tCF_percapita\tCF_perhh\tCF_pereqs\tCF_permeqs\tIE_overall\tIE_percapita\tDE_overall\tDE_percapita")
    for i=1:nn
        print(f, natList[i],"\t",length(hh_list[year][natList[i]]),"\t",natsam[i],"\t",natwgh[i])
        print(f, "\t",natcf[i],"\t",natcfpc[i],"\t",natcfph[i],"\t",natcfpeqs[i],"\t",natcfpmeqs[i])
        print(f, "\t", natie[i], "\t", natiepc[i], "\t", natde[i], "\t", natdepc[i])
        println(f)
    end
    close(f)
end

function exportRegionalEmission(years=[], tag="", outputFile=""; mode="cf", nspan=128, minmax=[], descend=false,empty=false,logarithm=false,nutsmode = "gis")
    # mode: "ie" indirect CE, "de" direct CE, "cf" integrated CF
    # nutsmode = "gis": NUTS codes follow GIS-map's NUTS (ex. DE1, DE2, DE3, ..., EL1, EL2, ...)
    # nutsmode = "hbs": NUTS codes follow HBS's NUTS (ex. DE0, DE3, DE4, ..., EL0, ...)

    global cat_list, nuts, nutsLv, nutsList, natList, popList, reg_sample, reg_avgExp, reg
    global gisRegList, hbscd, gispopcdlist, giscdlist, hbspopcdlist, hbscdlist
    global ieReg, gisRegIe, gisRegIeRank, gisRegIePerCap, gisRegIeRankPerCap
    global deReg, gisRegDe, gisRegDeRank, gisRegDePerCap, gisRegDeRankPerCap
    global cfReg, gisRegCf, gisRegCfRank, gisRegCfPerCap, gisRegCfRankPerCap

    nc = length(cat_list)
    labels = Dict{Int, Array{String,2}}()
    labelspc = Dict{Int, Array{String,2}}()

    cdlen = nutsLv+2

    for y in years
        nts = Dict{String, Array{String, 1}}()
        ntslist = Array{String, 1}()
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
                            aec[gidx] += reg_avgExp[y][hnt] * pop[y][pnt]
                        else
                            subpnts = filter(x->length(x)==cdlen+1, gispopcdlist[y][nt])
                            for spnt in subpnts
                                if haskey(pop[y], spnt)
                                    tb[gidx,:] += ec[ntidx,:] * pop[y][spnt]
                                    tpo[gidx] += pop[y][spnt]
                                    aec[gidx] += reg_avgExp[y][hnt] * pop[y][spnt]
                                end
                            end
                        end
                    end
                elseif nutsmode == "hbs"
                    ntidx = findfirst(x->x==nt, nutsList[y][n])
                    tb[gidx,:] += ec[ntidx,:] * popList[y][n][nt]
                    tpo[gidx] += popList[y][n][nt]
                    aec[gidx] += reg_avgExp[y][hnt] * popList[y][n][nt]
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

        gisPop[y] = tpo
        if nutsmode == "hbs"; gisSample[y] = spo end
        gisAvgExp[y] = aec
        gisRegList[y] = ntslist
        if mode == "ie"
            gisRegIe[y] = tb
            gisRegIePerCap[y] = tbpc
            gisRegIeRank[y] = rank
            gisRegIeRankPerCap[y] = rankpc
        elseif mode == "de"
            gisRegDe[y] = tb
            gisRegDePerCap[y] = tbpc
            gisRegDeRank[y] = rank
            gisRegDeRankPerCap[y] = rankpc
        elseif mode == "cf"
            gisRegCf[y] = tb
            gisRegCfPerCap[y] = tbpc
            gisRegCfRank[y] = rank
            gisRegCfRankPerCap[y] = rankpc
        else println("wrong emission categorizing mode: ", mode)
        end
    end

    return labels, labelspc
end

function exportRegionalTables(outputFile, tag, ntslist, nspan, minmax, tb, logarithm, descend; tab_mode="cf")
    # This function is for [exportRegionalEmission]

    global yr_list, cat_list
    nc = length(cat_list)
    nn = length(ntslist)

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
    print(f, tag); for c in cat_list; print(f,",",c) end; println(f)
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
    print(f, tag); for c in cat_list; print(f,",",c) end; println(f)
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

    global cat_list, nutsList, gisRegList, gisRegCf, gisRegCfPerCap
    global gisRegCfDiff, gisRegCfDiffRank, gisRegCfDiffPerCap, gisRegCfDiffRankPerCap

    nc = length(cat_list)
    spanval = Dict{Int, Array{Float64, 2}}()
    spanvalpc = Dict{Int, Array{Float64, 2}}()

    for y in years
        ntslist = gisRegList[y]
        gre = gisRegCf[y]
        grepc = gisRegCfPerCap[y]

        filename = replace(outputFile,"YEAR_"=>string(y)*"_")
        gred, rank, spanval[y] = exportEmissionDiffTable(replace(filename,"_OvPcTag"=>"_overall"), tag, ntslist, gre, maxr, minr, nspan, descend, empty)
        gredpc, rankpc, spanvalpc[y] = exportEmissionDiffTable(replace(filename,"_OvPcTag"=>"_percap"), tag, ntslist, grepc, maxr, minr, nspan, descend, empty)

        gisRegCfDiff[y] = gred
        gisRegCfDiffRank[y] = rank
        gisRegCfDiffPerCap[y] = gredpc
        gisRegCfDiffRankPerCap[y] = rankpc
    end

    return spanval, spanvalpc
end

function exportEmissionDiffTable(outputFile, tag, ntslist, gre, maxr, minr, nspan, descend, empty)
    # this function is for [exportEmissionDiffRate]

    global cat_list
    nc = length(cat_list)

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
    print(f, tag); for c in cat_list; print(f,",",c) end; println(f)
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
    print(f, tag); for c in cat_list; print(f,",",c) end; println(f)
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
    end
    for nt in ntslist
        mjcity = majorCity[year][nt]
        if modnuts && mjcity!="" && haskey(poplb[year], mjcity); nuts[year][nt] *= " (including "* poplb[year][mjcity] *")" end
    end
end

function exportWebsiteFiles(years, path; nutsmode = "gis", rank=false, empty=false, major=false)

    global natName, nuts, pop, popList, poplb, cat_list, gisRegList, gisPop, gisAvgExp
    global pophbscd, hbscd, gispopcdlist, hbspopcdlist, majorCity, gisCoord
    global gisRegCf, gisRegCfRank, gisRegCfDiffPerCap, gisRegCfDiffRankPerCap

    for y in years
        tbntlist = gisRegList[y]
        if nutsmode == "gis"; ntslist = filter(x->!(x in ["FR0","DE0","EL0"]), tbntlist)
        elseif nutsmode == "hbs"; ntslist = filter(x->!(x in ["FR0"]), tbntlist)
        end

        gre = gisRegCf[y]
        grer = gisRegCfRank[y]
        gredpc = gisRegCfDiffPerCap[y]
        gredrpc = gisRegCfDiffRankPerCap[y]

        # find major city of NUTS1
        if major; findMajorCity(y, ntslist, nutsmode, modnuts = true) end

        # print center file
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
        catidx = findfirst(x->gisCatLab[x]=="All", cat_list)
        for nt in ntslist; println(f, nt, "\t", grer[findfirst(x->x==nt, tbntlist), catidx]) end
        close(f)

        # print CF files
        for j=1:length(cat_list)
            mkpath(path*string(y)*"/CFAV/")
            f = open(path*string(y)*"/CFAV/"*"CFAV_"*gisCatLab[cat_list[j]]*".txt","w")
            print(f, "KEY_CODE\tSTATE\tDISTRICT\tSTATE_NAME\tDISTRICT_NAME\tGENHH_ALLP\tGENHH_APPPC")
            if cat_list[j]=="Total" || cat_list[j]=="All"; println(f, "\tANEXPPC\tPOP")
            else println(f)
            end
            for i=1:length(ntslist)
                nt = ntslist[i]
                tbidx = findfirst(x->x==nt, tbntlist)
                print(f, nt,"\t",nt[1:2],"\t",nt,"\t",natName[nt[1:2]],"\t",nuts[y][nt],"\t")
                printfmt(f, "{:f}", gre[tbidx,j]); print(f, "\t",gre[tbidx,j]/gisPop[y][tbidx])
                if cat_list[j]=="Total" || cat_list[j]=="All"; println(f,"\t",gisAvgExp[y][tbidx],"\t",convert(Int, gisPop[y][tbidx]))
                else println(f)
                end
            end
            if empty
                # for not-covered NUTS data
            end
            close(f)

            mkpath(path*string(y)*"/CFAC/")
            f = open(path*string(y)*"/CFAC/"*"CFAC_"*gisCatLab[cat_list[j]]*".txt","w")
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

    global natList, cat_list, nuts, natName, natA3, gisRegList, gisCatLab
    global gisRegIe, gisRegIeRank, gisRegCfDiffPerCap, gisRegCfDiffRankPerCap
    global gisPop, gisAvgExp

    cenfile = "centers.csv"
    engfile = "english_match.txt"
    allfile = "ALLP.txt"

    for y in years
        tbntlist = gisRegList[y]
        if nutsmode == "gis"; ntslist = filter(x->!(x in ["FR0","DE0","EL0"]), tbntlist)
        elseif nutsmode == "hbs"; ntslist = filter(x->!(x in ["FR0"]), tbntlist)
        end

        gre = gisRegIe[y]
        grer = gisRegIeRank[y]
        gredpc = gisRegCfDiffPerCap[y]
        gredrpc = gisRegCfDiffRankPerCap[y]
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
            catidx = findfirst(x->gisCatLab[x]=="All", cat_list)
            for nt in nts
                ntidx = findfirst(x->x==nt, tbntlist)
                println(f, nt, "\t", grer[ntidx,catidx])
            end
            close(f)

            # print CF files
            for j=1:length(cat_list)
                mkpath(fp*"CFAV/")
                f = open(fp*"CFAV/CFAV_"*gisCatLab[cat_list[j]]*".txt","w")
                print(f, "KEY_CODE\tSTATE\tDISTRICT\tSTATE_NAME\tDISTRICT_NAME\tGENHH_ALLP\tGENHH_APPPC")
                if gisCatLab[cat_list[j]]=="All"; println(f, "\tANEXPPC\tPOP")
                else println(f)
                end
                for nt in nts
                    ntidx = findfirst(x->x==nt, tbntlist)
                    print(f, nt,"\t",n,"\t",nt,"\t",natName[n],"\t",nuts[y][nt],"\t")
                    printfmt(f, "{:f}", gre[ntidx,j]); print(f, "\t",gre[ntidx,j]/gisPop[y][ntidx])
                    if gisCatLab[cat_list[j]]=="All"; println(f,"\t",gisAvgExp[y][ntidx],"\t",convert(Int, gisPop[y][ntidx]))
                    else println(f)
                    end
                end
                close(f)

                mkpath(fp*"CFAC/")
                f = open(fp*"CFAC/CFAC_"*gisCatLab[cat_list[j]]*".txt","w")
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

    global hhid, dis, siz, inc, pop, reg_sample, disPov
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
    else for i=1:nd; povr[i] /= reg_sample[regList[i]][1] end
    end

    for i=1:nd; disPov[regList[i]] = povr[i] end
end

function categorizeDistrictByEmissionLevel(year, normMode = 0, intv=[])
                                            # intv: proportions between invervals of highest to lowest
                                            # normmode: [1]per capita, [2]per household, [3]basic information
    global hhid, reg_sample, cat_list, regList
    global emissionsDis
    ed = emissionsDis[year]

    if length(intv) == 0; intv = [0.25,0.25,0.25,0.25]
    elseif sum(intv) != 1; sumi = sum(intv); intv /= sumi
    end

    nh = length(hhid)
    nc = length(cat_list)
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
        tp[indDis[i]] += reg_sample[regList[i]][1]
        th[indDis[i]] += reg_sample[regList[i]][2]
    end
    # normalizing
    if normMode == 1; for i=1:nc; edl[:,i] ./= tp end
    elseif normMode == 2 ;for i=1:nc; edl[:,i] ./= th end
    # basic information
    elseif normMode == 3; ed[:,1], ed[:,2] = tp[:], th[:]
    end

    emissionsDisLev[year] = edl

    return edl, cat_list, regList
end

function categorizeHouseholdByReligion(year, normMode=0; sqrRt=false, popWgh=false, wghmode="district")
                                            # normmode: [1]per capita, [2]per houehold, [3]basic information
    global wghSta, wghDis; if wghmode=="state"; wgh=wghSta elseif wghmode=="district"; wgh=wghDis end
    global hhid, cat, dis, siz, rel, cat_list, relList
    global ieHHs, emissionsRel

    nh = length(hhid)
    nc = length(cat_list)
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

    return er, cat_list, relList, tpbr, thbr, twpbr
end

function categorizeHouseholdByIncome(year,intv=[],normMode=0; sqrRt=false,absIntv=false,perCap=false,desOrd=false,popWgh=false,wghmode="district")
                                            # intv: proportions between invervals of highest to lowest
                                            # absIntv: if "true", then intv[] is a list of income values, descending order
                                            # normmode: [1]per capita, [2]per houehold
                                            # perCap: [true] per capita mode, [false] per household mode
                                            # desOrd: [true] descening order of 'intv[]', [false] ascending order
    global wghSta, wghDis; if wghmode=="state"; wgh=wghSta elseif wghmode=="district"; wgh=wghDis end
    global sec, hhid, cat, dis, siz, inc, cat_list, incList
    global ieHHs, emissionsInc

    if !absIntv && length(intv) == 0; intv = [0.25,0.5,0.75,1.00]
    elseif sort!(intv)[end] != 1; intv /= intv[end]
    end

    nh = length(hhid)
    nc = length(cat_list)
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

    return ei, cat_list, incList, tpbi, thbi, twpbi, indInc
end

function categorizeHouseholdByExpRange(year,rng=[],normMode=0; perCap=false,popWgh=false,absRng=false,absSpn=false,over=0.1,less=0.1,wghmode="district")
                                            # absRng: [true] apply absolute range, [false] apply population ratio range
                                            # rng: standard values of ranges for grouping
                                            # over/less: range from 'stdv', ratios of samples, househods, or population
                                            # normMode: [1]per capita, [2]per houehold
                                            # perCap: [true] per capita mode, [false] per household mode
    global wghSta, wghDis; if wghmode=="state"; wgh=wghSta elseif wghmode=="district"; wgh=wghDis end
    global sec, hhid, cat, dis, siz, inc, cat_list, incList
    global ieHHs, emissionsRng

    sort!(rng)

    nh = length(hhid)
    nc = length(cat_list)
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

    return er, cat_list, thber, tpber, twpber, rsidx, expOrder
end

function categorizeHouseholdByEmissionLevel(year, intv=[], normMode = 0; squareRoot = false, absintv=false)
                                                    # intv: proportions between invervals of highest to lowest
                                                    # normmode: [1]per capita, [2]per houehold
    global hhid, sec, cat_list, levList, cat, dis, siz
    global ieHHs, emissionsLev

    if !absintv && length(intv) == 0; intv = [0.25,0.25,0.25,0.25]
    elseif !absintv && sum(intv) != 1; sumi = sum(intv); intv /= sumi
    end

    ec = ieHHs[year]
    nh = length(hhid)
    nc = length(cat_list)
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
    global sec, hhid, cat, dis, siz, inc, rel, cat_list, incList, relList
    global ieHHs, emissionsIncRel

    if !absIntv && length(intv) == 0; intv = [0.25,0.5,0.75,1.00]
    elseif sort!(intv)[end] != 1; intv /= intv[end]
    end

    nh = length(hhid)
    nc = length(cat_list)
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

    return eir, cat_list, incList, tpbir, thbir, twpbir
end

function estimateEmissionCostByDistrict(year,expIntv=[],normMode=0; perCap=false,popWgh=false,desOrd=false,name=false,output="",exportFile=[],wghmode="district")
    # Estimate poverty alleviation leading CF (CO2 emissions) increase by district
    # considering expenditure-level, and distric consumption patterns
    # intv: proportions between invervals of highest to lowest
    # normmode: [1]per capita, [2]per houehold
    # perCap: [true] per capita mode, [false] per household mode

    global wghSta, wghDis; if wghmode=="state"; wgh=wghSta elseif wghmode=="district"; wgh=wghDis end
    global sec, hhid, siz, cat, dis, typ, inc, rel, pop, reg_sample
    global staList, regList, cat_list, incList, relList, typList, disSta
    global ieHHs, emissionCostDis

    nh = length(hhid)
    ns = length(staList)
    nd = length(regList)
    nt = length(typList)+1
    nc = length(cat_list)
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
            for i=1:nd; ecpc[:,i,:] ./= reg_sample[regList[i]][1] end
        else
            for i=2:ne; for j=1:nh; if idxExp[j]<i
                cost = avg[idxDis[j],idxTyp[j],i,:] - eh[j,:]
                if cost>0; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= reg_sample[regList[i]][2] end
        end
    else
        if normMode == 1; for i=1:ne-1; for j=1:nh;
            if idxExp[j]>i
                cost = avg[idxDis[j],idxTyp[j],i,:]*siz[hhid[j]] - eh[j,:]
                if cost[end]>0 ; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= reg_sample[regList[i]][1] end
        else for i=1:ne-1; for j=1:nh;
            if idxExp[j]>i
                cost = avg[idxDis[j],idxTyp[j],i,:] - eh[j,:]
                if cost>0; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= reg_sample[regList[i]][2] end
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
            print(f,"District"); for j=1:nc; print(f,",",cat_list[j]) end; println(f)
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
            print(f,"District"); for j=1:nc; print(f,",",cat_list[j]) end; println(f)
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
        global reg_avgExp, gid, gidData

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
            spo[idx] += reg_sample[regList[i]][1]
            shh[idx] += reg_sample[regList[i]][2]
            tpo[idx] += pop[regList[i]][1]
            thh[idx] += pop[regList[i]][2]
            aec[idx] += reg_avgExp[regList[i]]*reg_sample[regList[i]][1]
        end
        for i=1:ng; aec[i] /= spo[i] end
        if normMode==1; for i=1:ng; for j=1:nc; tcpc[i,j] = tc[i,j]/tpo[i] end end
        elseif normMode==2; for i=1:ng; for j=1:nc; tcpc[i,j] = tc[i,j]/thh[i] end end
        end

        # exporting total emission cost
        f = open(replace(exportFile[1],".csv"=>"_total.csv"), "w")
        print(f, exportFile[2]); for c in cat_list; print(f,",",c) end; println(f)
        for i = 1:size(tc, 1)
            print(f, gidList[i])
            for j = 1:size(tc, 2); print(f, ",", tc[i,j]) end
            println(f)
        end
        close(f)
        # exporting emission cost per capita
        f = open(replace(exportFile[1],".csv"=>"_perCap.csv"), "w")
        print(f, exportFile[2]); for c in cat_list; print(f,",",c) end; println(f)
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
    global sec, hhid, siz, cat, dis, typ, inc, rel, pop, reg_sample
    global staList, regList, cat_list, incList, relList, typList, disSta
    global ieHHs, emissionCostDis

    nh = length(hhid)
    ns = length(staList)
    nd = length(regList)
    nt = length(typList)+1
    nc = length(cat_list)
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
            for i=1:nd; ecpc[:,i,:] ./= reg_sample[regList[i]][1] end
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
            for i=1:nd; ecpc[:,i,:] ./= reg_sample[regList[i]][2] end
        end
    else    # incomplete code
        if normMode == 1; for i=1:nr; for j=1:nh;
            if idxExp[j]>i
                cost = avg[idxDis[j],idxTyp[j],i,:]*siz[hhid[j]] - eh[j,:]
                if cost[end]>0 ; ecpc[i,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= reg_sample[regList[i]][1] end
        else for i=1:nr; for j=1:nh;
            if idxExp[j]>i
                cost = avg[idxDis[j],idxTyp[j],i,:] - eh[j,:]
                if cost>0; ecpc[i,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= reg_sample[regList[i]][2] end
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
            print(f,"District"); for j=1:nc; print(f,",",cat_list[j]) end
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
            print(f,"District"); for j=1:nc; print(f,",",cat_list[j]) end
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
        global reg_avgExp, gid, gidData

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
            spo[idx] += reg_sample[regList[i]][1]
            shh[idx] += reg_sample[regList[i]][2]
            tpo[idx] += pop[regList[i]][1]
            thh[idx] += pop[regList[i]][2]
            aec[idx] += reg_avgExp[regList[i]]*reg_sample[regList[i]][1]
        end
        for i=1:ng; aec[i] /= spo[i] end
        if normMode==1; for i=1:ng; for j=1:nc; tcpc[i,j] = tc[i,j]/tpo[i] end end
        elseif normMode==2; for i=1:ng; for j=1:nc; tcpc[i,j] = tc[i,j]/thh[i] end end
        end

        # exporting total emission cost
        f = open(replace(exportFile[1],".csv"=>"_total.csv"), "w")
        print(f, exportFile[2]); for c in cat_list; print(f,",",c) end; println(f)
        for i = 1:size(tc, 1)
            print(f, gidList[i])
            for j = 1:size(tc, 2); print(f, ",", tc[i,j]) end
            println(f)
        end
        close(f)
        # exporting emission cost per capita
        f = open(replace(exportFile[1],".csv"=>"_perCap.csv"), "w")
        print(f, exportFile[2]); for c in cat_list; print(f,",",c) end; println(f)
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
    global sec, hhid, siz, cat, dis, typ, inc, rel, pop, reg_sample
    global staList, regList, cat_list, incList, relList, typList, disSta
    global ieHHs, emissionCostDis

    nh = length(hhid)
    ns = length(staList)
    nd = length(regList)
    nt = length(typList)+1
    nc = length(cat_list)
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
            for i=1:nd; ecpc[:,i,:] ./= reg_sample[regList[i]][1] end
        else
            for i=2:ne; for j=1:nh; if idxExp[j]<i
                cost = avg[idxDis[j],idxTyp[j],i,idxRel[j],:] - eh[j,:]
                if cost>0; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= reg_sample[regList[i]][2] end
        end
    else
        if normMode == 1; for i=1:ne-1; for j=1:nh;
            if idxExp[j]>i
                cost = avg[idxDis[j],idxTyp[j],i,idxRel[j],:]*siz[hhid[j]] - eh[j,:]
                if cost[end]>0 ; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= reg_sample[regList[i]][1] end
        else for i=1:ne-1; for j=1:nh;
            if idxExp[j]>i
                cost = avg[idxDis[j],idxTyp[j],i,idxRel[j],:] - eh[j,:]
                if cost>0; ecpc[i-1,idxDis[j],:] += cost end
            end end end
            for i=1:nd; ecpc[:,i,:] ./= reg_sample[regList[i]][2] end
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
            print(f,"District"); for j=1:nc; print(f,",",cat_list[j]) end; println(f)
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
            print(f,"District"); for j=1:nc; print(f,",",cat_list[j]) end; println(f)
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
        global reg_avgExp, gid, gidData

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
            spo[idx] += reg_sample[regList[i]][1]
            shh[idx] += reg_sample[regList[i]][2]
            tpo[idx] += pop[regList[i]][1]
            thh[idx] += pop[regList[i]][2]
            aec[idx] += reg_avgExp[regList[i]]*reg_sample[regList[i]][1]
        end
        for i=1:ng; aec[i] /= spo[i] end
        if normMode==1; for i=1:ng; for j=1:nc; tcpc[i,j] = tc[i,j]/tpo[i] end end
        elseif normMode==2; for i=1:ng; for j=1:nc; tcpc[i,j] = tc[i,j]/thh[i] end end
        end

        # exporting total emission cost
        f = open(replace(exportFile[1],".csv"=>"_total.csv"), "w")
        print(f, exportFile[2]); for c in cat_list; print(f,",",c) end; println(f)
        for i = 1:size(tc, 1)
            print(f, gidList[i])
            for j = 1:size(tc, 2); print(f, ",", tc[i,j]) end
            println(f)
        end
        close(f)
        # exporting emission cost per capita
        f = open(replace(exportFile[1],".csv"=>"_perCap.csv"), "w")
        print(f, exportFile[2]); for c in cat_list; print(f,",",c) end; println(f)
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

    global cat_list, relList, relName, emissionsRel
    er = emissionsRel[year]

    f = open(outputFile, "w")

    print(f,"Religion")
    for c in cat_list; print(f, ",", c) end
    if length(tpbr)>0; print(f,",Pop.") end
    if length(thbr)>0; print(f,",HH.") end
    if length(twpbr)>0; print(f,",WghPop.") end
    println(f)
    for i = 1:length(relList)
        print(f, relName[i])
        for j = 1:length(cat_list); print(f, ",", er[i,j]) end
        if length(tpbr)>0; print(f,",",tpbr[i]) end
        if length(thbr)>0; print(f,",",thbr[i]) end
        if length(twpbr)>0; print(f,",",twpbr[i]) end
        println(f)
    end

    close(f)
end

function printEmissionByIncome(year, outputFile, intv=[], tpbi=[], thbi=[], twpbi=[]; absIntv=false, desOrd=false, relative=false)

    global cat_list, incList, emissionsInc
    ei = emissionsInc[year]
    ni = length(intv); if absIntv; ni += 1 end
    nc = length(cat_list)

    f = open(outputFile, "w")
    print(f,"Exp_Lv,Value")
    for c in cat_list; print(f, ",", c) end
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
        for c in cat_list; print(f, ",", c) end
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

    global cat_list, emissionsRng
    er = emissionsRng[year]
    nr = size(rsidx,1)

    f = open(outputFile, "w")

    print(f,"Expenditure,Range")
    for c in cat_list; print(f, ",", c) end
    if length(tpber)>0; print(f,",Pop.") end
    if length(thber)>0; print(f,",HH.") end
    if length(twpber)>0; print(f,",WghPop.") end
    println(f)
    for i = 1:nr
        print(f,round(inc[hhid[order[rsidx[i,1]]]],digits=2))
        print(f,",",round(inc[hhid[order[rsidx[i,2]]]],digits=2),"-",round(inc[hhid[order[rsidx[i,3]]]],digits=2))
        for j = 1:length(cat_list); print(f, ",", er[i,j]) end
        if length(tpber)>0; print(f,",",tpber[i]) end
        if length(thber)>0; print(f,",",thber[i]) end
        if length(twpber)>0; print(f,",",twpber[i]) end
        println(f)
    end

    close(f)
end

function printEmissionByIncomeByReligion(year, outputFile, intv=[], tpbir=[], thbir=[], twpbir=[]; absIntv=false, desOrd=false)

    global cat_list, incList, relList, emissionsIncRel
    eir = emissionsIncRel[year]
    nc = length(cat_list)
    nr = length(relList)
    ni = length(intv); if absIntv; ni += 1 end
    relName = ["Hinduism","Islam","Christianity","Sikhism","Jainism","Buddhism","Zoroastrianism", "Others", "None"]

    f = open(outputFile, "w")

    for i = 1:nr
        println(f, relName[i])
        print(f,"Exp_Lv")
        for c in cat_list; print(f, ",", c) end
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

    global cat_list, regList, emissionsDisLev
    edl = emissionsDisLev[year]

    f = open(outputFile, "w")

    print(f,"CF_Lv")
    for c in cat_list; print(f, ",", c) end
    println(f)
    for i = 1:length(intv)
        print(f, "\t<", trunc(Int, sum(intv[1:i])*100),"%")
        for j = 1:length(cat_list); print(f, ",", edl[i,j]) end
        println(f)
    end

    close(f)
end

function printEmissionByHhsEmLev(year, outputFile, intv=[])

    global cat_list, regList, emissionsLev
    el = emissionsLev[year]

    f = open(outputFile, "w")

    print(f,"CF_Lv")
    for c in cat_list; print(f, ",", c) end
    println(f)
    for i = 1:length(intv)
        print(f, "\t<", trunc(Int, sum(intv[1:i])*100),"%")
        for j = 1:length(cat_list); print(f, ",", el[i,j]) end
        println(f)
    end

    close(f)
end

end
