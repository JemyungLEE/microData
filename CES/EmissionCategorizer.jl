module EmissionCategorizer

# Developed date: 17. May. 2021
# Last modified date: 29. May. 2021
# Subject: Categorize households' carbon footprints
# Description: Read household-level indirect and direct carbon emissions,  integrate them to be CF,
#              and categorize the CFs by consumption category, district, expenditure-level, and etc.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

include("MicroDataReader.jl")

using Statistics
using Formatting: printfmt
using .MicroDataReader
mdr = MicroDataReader

yr_list = Array{Int, 1}()       # year list: {YYYY}
nat_list = Array{String, 1}()   # nation list: {A3}
cat_list = Array{String, 1}()   # category list
rel_list = Array{String, 1}()   # religion list
pr_unts = Dict("day" => 1, "week" => 7,"month" => 30, "year" => 365)

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
majorCity = Dict{Int, Dict{String, Dict{String, String}}}()     # major city in the region: {year, {nation, {Upper_region_code, major_city_code}}}
gisCoord = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}() # GIS coordinates: {year, {nation, {region_code, {X_longitude, Y_latitude}}}}
gisCatLab = Dict{String, String}()              # Web-exporting category matching table: {Category label in program, in web-site files}

gisRegList = Dict{Int, Dict{String, Array{String, 1}}}()    # GIS region list: {year, {nation, {region code}}}
gisRegLabel= Dict{Int, Dict{String, Dict{String,String}}}() # GIS region label: {year, {nation, {region_code, region_label}}}
gisRegProv = Dict{Int, Dict{String, Dict{String, Tuple{String, String}}}}() # GIS upper region's code and label: {year, {nation, {region_code, {upper region_code, label}}}}
gisRegConc = Dict{Int, Dict{String, Array{Float64, 2}}}()   # GIS-CES/HBS region concordance weight: {year, {nation, {gis_code, ces/hbs_code}}}

gisPop = Dict{Int, Dict{String, Array{Float64, 1}}}()       # GIS version, population by region: {year, {nation, {population}}}
gisSample = Dict{Int, Dict{String, Array{Float64, 1}}}()    # GIS version, total samples by region: {year, {nation, {samples}}}
gisAvgExp = Dict{Int, Dict{String, Array{Float64, 1}}}()    # GIS version, average expenditure by region: {year, {nation, {average expenditure}}}

ieRegGIS = Dict{Int, Dict{String, Array{Float64, 2}}}()     # categozied indirect emission by region: {year, {nation, {region, category}}}
deRegGIS = Dict{Int, Dict{String, Array{Float64, 2}}}()     # categozied direct emission by region: {year, {nation, {region, category}}}
cfRegGIS = Dict{Int, Dict{String, Array{Float64, 2}}}()     # categozied carbon footprint by region: {year, {nation, {region, category}}}
ieRegRankGIS = Dict{Int, Dict{String, Array{Int, 2}}}()     # categozied indirect emission rank by region: {year, {nation, {region, category}}}
deRegRankGIS = Dict{Int, Dict{String, Array{Int, 2}}}()     # categozied direct emission rank by region: {year, {nation, {region, category}}}
cfRegRankGIS = Dict{Int, Dict{String, Array{Int, 2}}}()     # categozied carbon footprint rank by region: {year, {nation, {region, category}}}
ieRegPcGIS = Dict{Int, Dict{String, Array{Float64, 2}}}()   # categozied indirect emission per capita by region: {year, {nation, {region, category}}}
deRegPcGIS = Dict{Int, Dict{String, Array{Float64, 2}}}()   # categozied direct emission per capita by region: {year, {nation, {region, category}}}
cfRegPcGIS = Dict{Int, Dict{String, Array{Float64, 2}}}()   # categozied carbon footprint per capita by region: {year, {nation, {region, category}}}
ieRegPcRankGIS = Dict{Int, Dict{String, Array{Int, 2}}}()   # categozied indirect emission per capita rank by region: {year, {nation, {region, category}}}
deRegPcRankGIS = Dict{Int, Dict{String, Array{Int, 2}}}()   # categozied direct emission per capita rank by region: {year, {nation, {region, category}}}
cfRegPcRankGIS = Dict{Int, Dict{String, Array{Int, 2}}}()   # categozied carbon footprint per capita rank by region: {year, {nation, {region, category}}}

ieRegDevGIS = Dict{Int, Dict{String, Array{Float64, 2}}}()  # categozied indirect emission deviation from mean by region: {year, {nation, {region, category}}}
ieRegDevRankGIS = Dict{Int, Dict{String, Array{Int, 2}}}()  # categozied indirect emission deviation from mean rank by region: {year, {nation, {region, category}}}
ieRegDevPcGIS = Dict{Int, Dict{String, Array{Float64, 2}}}()# categozied indirect emission per capita deviation from mean by region: {year, {nation, {region, category}}}
ieRegDevPcRankGIS = Dict{Int, Dict{String, Array{Int, 2}}}()# categozied indirect emission per capita deviation from mean rank by region: {year, {nation, {region, category}}}
deRegDevGIS = Dict{Int, Dict{String, Array{Float64, 2}}}()  # categozied direct emission deviation from mean by region: {year, {nation, {region, category}}}
deRegDevRankGIS = Dict{Int, Dict{String, Array{Int, 2}}}()  # categozied direct emission deviation from mean rank by region: {year, {nation, {region, category}}}
deRegDevPcGIS = Dict{Int, Dict{String, Array{Float64, 2}}}()# categozied direct emission per capita deviation from mean by region: {year, {nation, {region, category}}}
deRegDevPcRankGIS = Dict{Int, Dict{String, Array{Int, 2}}}()# categozied direct emission per capita deviation from mean rank by region: {year, {nation, {region, category}}}
cfRegDevGIS = Dict{Int, Dict{String, Array{Float64, 2}}}()  # categozied carbon footprint deviation from mean by region: {year, {nation, {region, category}}}
cfRegDevRankGIS = Dict{Int, Dict{String, Array{Int, 2}}}()  # categozied carbon footprint deviation from mean rank by region: {year, {nation, {region, category}}}
cfRegDevPcGIS = Dict{Int, Dict{String, Array{Float64, 2}}}()# categozied carbon footprint per capita deviation from mean by region: {year, {nation, {region, category}}}
cfRegDevPcRankGIS = Dict{Int, Dict{String, Array{Int, 2}}}()# categozied carbon footprint per capita deviation from mean rank by region: {year, {nation, {region, category}}}

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

function setCategory(year, nation; subgroup = "", except=[])  # Note: sub-grouping parts should be added

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

    global yr_list, nat_list, hh_list, sc_list, sc_cat, cat_list, rel_list, prov_list, dist_list, pr_unts
    global households, pops, pops_ur, hh_period, reg_sample, reg_avgExp
    global ieHHs, deHHs, cfHHs, ieReg, deReg, cfReg, ieRegDev, deRegDev, cfRegDev

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

function makeNationalSummary(years=[], nations=[], outputFile=""; region = "district")

    global yr_list, nat_list, hh_list, households, pops, pops_ur, pop_wgh, pop_ur_wgh
    global ieHHs, deHHs, cfHHs
    if length(years) == 0; years = yr_list elseif isa(years, Number); years = [years] end
    if length(nations) == 0; nations = nat_list elseif isa(nations, String); nations = [nations] end

    ny, nn = length(years), length(nations)
    natsam = zeros(Int, ny, nn)
    natwgh = zeros(Float64, ny, nn)
    natie = zeros(Float64, ny, nn)      # Overall IE
    natiepc = zeros(Float64, ny, nn)    # IE per capita
    natde = zeros(Float64, ny, nn)      # Overall DE
    natdepc = zeros(Float64, ny, nn)    # DE per capita
    natcf = zeros(Float64, ny, nn)      # Overall CF
    natcfpc = zeros(Float64, ny, nn)    # CF per capita

    for i = 1:ny, j = 1:nn
        y, n, hhs = years[i], nations[j], households[y][n]
        for k = 1:length(hh_list[y][n])
            h = hh_list[y][n][k]
            if region == "district"; pw = pop_wgh[y][h][hhs[h].district]
            elseif region == "province"; pw = pop_wgh[y][h][hhs[h].province]
            end
            ie = sum(ieHHs[y][n][k,:])
            de = sum(deHHs[y][n][k,:])
            cf = sum(cfHHs[y][n][k,:])
            natsam[i,j] += hhs[h].size
            natwgh[i,j] += pw * hhs[h].size
            natcf[i,j] += pw * cf
            natcfpc[i,j] += cf
            natie[i,j] += pw * ie
            natiepc[i,j] += ie
            natde[i,j] += pw * de
            natdepc[i,j] += de
        end
        natcfpc[i,j] /= natsam[i,j]
        natiepc[i,j] /= natsam[i,j]
        natdepc[i,j] /= natsam[i,j]
    end

    f_sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")
    items = ["Year","Nation","HHs","MMs","Weights","CF_overall","CF_percapita","IE_overall","IE_percapita","DE_overall","DE_percapita"]
    for it in items; print(f, it, f_sep) end; println(f)
    for i = 1:ny, j = 1:nn
        y, n = years[i], nations[j]
        print(f, y, f_sep, n, f_sep, length(hh_list[y][n]), f_sep, natsam[i,j], f_sep, natwgh[i,j])
        print(f, f_sep, natcf[i,j], f_sep, natcfpc[i,j], f_sep, natie[i,j], f_sep, natiepc[i,j], f_sep, natde[i,j], f_sep, natdepc[i,j])
        println(f)
    end
    close(f)
end

function readGISinfo(years=[], nations=[], regionFile="", gisCatFile="")

    global yr_list, nat_list, gisCoord, gisRegList, gisRegLabel, gisRegProv, gisCatLab
    if length(years) == 0; years = yr_list elseif isa(years, Number); years = [years] end
    if length(nations) == 0; nations = nat_list elseif isa(nations, String); nations = [nations] end

    essential = ["Region_ID", "Region_label", "Province_ID", "Province_label", "Nation_A3", "X_coord", "Y_coord"]

    str = Array{Array{String, 1}, 1}()
    f_sep = getValueSeparator(regionFile)
    f = open(regionFile)
    title = string.(strip.(split(readline(f), f_sep)))
    i = [findfirst(x->x==t, title) for t in essential]
    for l in eachline(f); push!(str, string.(strip.(split(l, f_sep)))[i]) end
    close(f)

    for y in years, n in nations
        coord = Dict{String, Tuple{Float64, Float64}}()
        reg_list, reg_lab, reg_prov = Array{String, 1}(), Dict{String, String}(), Dict{String, Tuple{String, String}}()

        for s in str
            if s[5] == n
                if !(s[1] in reg_list); push!(reg_list, s[1]) end
                if !haskey(reg_lab, s[1]); reg_lab[s[1]] = s[2] end
                coord[s[1]] = (parse(Float64, s[6]), parse(Float64, s[7]))
                reg_prov[s[1]] = (s[3], s[4])
            end
        end

        if !haskey(gisCoord, y); gisCoord[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
        if !haskey(gisRegList, y); gisRegList[y] = Dict{String, Array{String, 1}}() end
        if !haskey(gisRegLabel, y); gisRegLabel[y] = Dict{String, Dict{String, String}}() end
        if !haskey(gisRegProv, y); gisRegProv[y] = Dict{String, Dict{String, Tuple{String, String}}}() end
        gisCoord[y][n], gisRegList[y][n], gisRegLabel[y][n], gisRegProv[y][n] = coord, sort(reg_list), reg_lab, reg_prov
    end

    f_sep = getValueSeparator(gisCatFile)
    f = open(gisCatFile)
    title = string.(strip.(split(readline(f), f_sep)))
    i = [findfirst(x->x==t, title) for t in ["CES/HBS_category", "GIS_label"]]
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        gisCatLab[s[i[1]]] = s[i[2]]
    end
    close(f)
end

function buildGISconc(years=[], nations=[], gisConcFile=""; region = "district", remove=false)

    global yr_list, nat_list, dist_list, prov_list, gisRegList, gisRegConc
    if length(years) == 0; years = yr_list elseif isa(years, Number); years = [years] end
    if length(nations) == 0; nations = nat_list elseif isa(nations, String); nations = [nations] end
    if region == "district"; reg_list = dist_list elseif region == "province"; reg_list = prov_list end

    links = Array{Tuple{String, String, Float64}, 1}()
    gis_id = Array{String, 1}()
    f_sep = getValueSeparator(gisConcFile)
    f = open(gisConcFile)
    title = string.(strip.(split(readline(f), f_sep)))
    i = [findfirst(x->x==t, title) for t in ["GIS_region_ID", "CES/HBS_region_ID", "Weight"]]
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        push!(links, (s[i[1]], s[i[2]], parse(Float64, s[i[3]])))
        push!(gis_id, s[i[1]])
    end
    close(f)

    for y in years, n in nations
        rl, grl = reg_list[y][n], gisRegList[y][n]
        if remove; grl = gisRegList[y][n] = filter(x->x in gis_id, grl) end
        nr, ngr = length(rl), length(grl)
        conc = zeros(Float64, ngr, nr)
        for l in links; conc[findfirst(x->x==l[1], grl), findfirst(x->x==l[2], rl)] += l[3] end
        for i = 1:nr; conc[:,i] /= sum(conc[:,i]) end
        if !haskey(gisRegConc, y); gisRegConc[y] = Dict{String, Array{Float64, 2}}() end
        gisRegConc[y][n] = conc
    end

end

function exportRegionalEmission(years=[],nations=[],tag="",outputFile=""; region="district",mode=["cf"],nspan=128,minmax=[],descend=false,empty=false,logarithm=false)

    global yr_list, nat_list, cat_list, prov_list, dist_list, pops, ieReg, deReg, cfReg, ieRegDev, deRegDev, cfRegDev
    global reg_sample, reg_avgExp, gisRegList, gisRegConc, gisPop, gisSample, gisAvgExp
    global ieRegGIS, deRegGIS, cfRegGIS, ieRegRankGIS, deRegRankGIS, cfRegRankGIS
    global ieRegPcGIS, deRegPcGIS, cfRegPcGIS, ieRegPcRankGIS, deRegPcRankGIS, cfRegPcRankGIS
    global cfRegDevGIS, cfRegDevRankGIS, cfRegDevPcGIS, cfRegDevPcRankGIS

    nc = length(cat_list)
    if length(years) == 0; years = yr_list elseif isa(years, Number); years = [years] end
    if length(nations) == 0; nations = nat_list elseif isa(nations, String); nations = [nations] end
    if region == "district"; reg_list = dist_list elseif region == "province"; reg_list = prov_list end
    if isa(mode, String); mode = [mode] end
    labels, labelspc = Dict{Int, Dict{String, Array{String,2}}}(), Dict{Int, Dict{String, Array{String,2}}}()

    for y in years, n in nations
        rl, grl, con = reg_list[y][n], gisRegList[y][n], gisRegConc[y][n]
        nr, ngr = length(rl), length(grl)

        # r_pop = [pops[y][n][r] for r in rl]
        # r_sam = [reg_sample[y][n][r][1] for r in rl]
        # r_ave = [reg_avgExp[y][n][r] * pops[y][n][r] for r in rl]

        g_sam = con * [reg_sample[y][n][r][1] for r in rl]
        g_pop = con * [pops[y][n][r] for r in rl]
        g_ave = con * [reg_avgExp[y][n][r] * pops[y][n][r] for r in rl]

        if !haskey(gisPop, y); gisPop[y] = Dict{String, Array{Float64, 1}}() end
        if !haskey(gisSample, y); gisSample[y] = Dict{String, Array{Float64, 1}}() end
        if !haskey(gisAvgExp, y); gisAvgExp[y] = Dict{String, Array{Float64, 1}}() end
        gisPop[y][n], gisSample[y][n], gisAvgExp[y][n] = g_pop, g_sam, g_ave

        for m in mode
            if m == "ie"; ec = ieReg[y][n]; elseif m == "de"; ec = deReg[y][n]; elseif m == "cf"; ec = cfReg[y][n] end

            g_ce = con * ec     # {ngr, nc}
            g_cepc = g_ce ./ g_pop
            g_ave = g_ave ./ g_pop

            if m == "ie"
                if !haskey(ieRegGIS, y); ieRegGIS[y] = Dict{String, Array{Float64, 2}}() end
                if !haskey(ieRegPcGIS, y); ieRegPcGIS[y] = Dict{String, Array{Float64, 2}}() end
                ieRegGIS[y][n], ieRegPcGIS[y][n] = g_ce, g_cepc
            elseif m == "de"
                if !haskey(deRegGIS, y); deRegGIS[y] = Dict{String, Array{Float64, 2}}() end
                if !haskey(deRegPcGIS, y); deRegPcGIS[y] = Dict{String, Array{Float64, 2}}() end
                deRegGIS[y][n], deRegPcGIS[y][n] = g_ce, g_cepc
            elseif m == "cf"
                if !haskey(cfRegGIS, y); cfRegGIS[y] = Dict{String, Array{Float64, 2}}() end
                if !haskey(cfRegPcGIS, y); cfRegPcGIS[y] = Dict{String, Array{Float64, 2}}() end
                cfRegGIS[y][n], cfRegPcGIS[y][n] = g_ce, g_cepc
            end

            if !haskey(labels, y); labels[y] = Dict{String, Array{String,2}}() end
            if !haskey(labelspc, y); labelspc[y] = Dict{String, Array{String,2}}() end
            fns = rsplit(outputFile, '.', limit=2)
            filename = replace(fns[1],"YEAR_"=>string(y)*"_") * "_" * uppercase(m) * "." * fns[2]
            rank, labels[y][n] = exportRegionalTables(replace(filename,"_OvPcTag"=>"_overall"), tag, rl, nspan, minmax[1], g_ce, logarithm, descend)
            rankpc, labelspc[y][n] = exportRegionalTables(replace(filename,"_OvPcTag"=>"_percap"), tag, rl, nspan, minmax[2], g_cepc, logarithm, descend)

            if m == "ie"
                if !haskey(ieRegRankGIS, y); ieRegRankGIS[y] = Dict{String, Array{Int, 2}}() end
                if !haskey(ieRegPcRankGIS, y); ieRegPcRankGIS[y] = Dict{String, Array{Int, 2}}() end
                ieRegRankGIS[y][n], ieRegPcRankGIS[y][n] = rank, rankpc
            elseif m == "de"
                if !haskey(deRegRankGIS, y); deRegRankGIS[y] = Dict{String, Array{Int, 2}}() end
                if !haskey(deRegPcRankGIS, y); deRegPcRankGIS[y] = Dict{String, Array{Int, 2}}() end
                deRegRankGIS[y][n], deRegPcRankGIS[y][n] = rank, rankpc
            elseif m == "cf"
                if !haskey(cfRegRankGIS, y); cfRegRankGIS[y] = Dict{String, Array{Int, 2}}() end
                if !haskey(cfRegPcRankGIS, y); cfRegPcRankGIS[y] = Dict{String, Array{Int, 2}}() end
                cfRegRankGIS[y][n], cfRegPcRankGIS[y][n] = rank, rankpc
            end
        end
    end

    return labels, labelspc
end

function exportRegionalTables(outputFile="", tag="", reg_list=[], nspan=[], minmax=[], ce=[], logarithm=false, descend=false)
    # This function is for [exportRegionalEmission]

    global cat_list
    nc, ngr = length(cat_list), length(reg_list)

    # find min. and max.: overall CF
    if length(minmax)==1; maxce = [minmax[1][2] for i=1:nc]; mince = [minmax[1][1] for i=1:nc]
    elseif length(minmax)==nc; maxce = [minmax[i][2] for i=1:nc]; mince = [minmax[i][1] for i=1:nc]
    elseif logarithm; maxce = [log10(maximum(ce[:,i])) for i=1:nc]; mince = [log10(minimum(ce[:,i])) for i=1:nc]
    else maxce = [maximum(ce[:,i]) for i=1:nc]; mince = [minimum(ce[:,i]) for i=1:nc]
    end
    replace!(mince, Inf=>0, -Inf=>0)
    # grouping by ratios; ascending order: overall CF
    span = zeros(Float64, nspan+1, nc)
    over = [maxce[i] < maximum(ce[:,i]) for i=1:nc]
    for j = 1:nc
        if over[j]; span[:,j] = [[(maxce[j]-mince[j])*(i-1)/(nspan-1)+mince[j] for i=1:nspan]; maximum(ce[:,j])]
        else span[:,j] = [(maxce[j]-mince[j])*(i-1)/nspan+mince[j] for i=1:nspan+1]
        end
    end
    if logarithm; for i=1:size(span,1), j=1:nc; span[i,j] = 10^span[i,j] end end
    # grouping by rank; ascending order
    rank = zeros(Int, ngr, nc)
    for i = 1:ngr, j = 1:nc
        if ce[i,j]>=span[end-1,j]; rank[i,j] = nspan
        elseif ce[i,j] <= span[1,j]; rank[i,j] = 1
        else rank[i,j] = findfirst(x->x>=ce[i,j],span[:,j]) - 1
        end
    end
    # for descending order, if "descend == true"
    if descend
        for i = 1:nc; span[:,i] = reverse(span[:,i]) end
        for i = 1:ngr, j = 1:nc; rank[i,j] = nspan - rank[i,j] + 1 end
    end
    # prepare labels
    labels = Array{String, 2}(undef, nspan, nc)
    for j = 1:nc
        lbstr = [string(round(span[i,j],digits=0)) for i=1:nspan+1]
        if descend; labels[:,j] = [lbstr[i+1]*"-"*lbstr[i] for i=1:nspan]
        else labels[:,j] = [lbstr[i]*"-"*lbstr[i+1] for i=1:nspan]
        end
        if over[j]; if descend; labels[1,j] = "over "*lbstr[2] else labels[nspan,j] = "over "*lbstr[nspan] end end
    end
    # exporting table: overall CF
    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f_sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")
    print(f, tag); for c in cat_list; print(f, f_sep, c) end; println(f)
    for i = 1:ngr
        print(f, reg_list[i])
        for j = 1:nc; print(f, f_sep, ce[i,j]) end
        println(f)
    end
    if empty
        # for empty GIS region
    end
    close(f)
    # exporting group table: overall CF
    fns = rsplit(outputFile, '.', limit=2)
    f = open(fns[1]*"_gr."*fns[2], "w")
    print(f, tag); for c in cat_list; print(f, f_sep, c) end; println(f)
    for i = 1:ngr
        print(f, reg_list[i])
        for j = 1:nc; print(f, f_sep, rank[i,j]) end
        println(f)
    end
    if empty
        # for empty GIS region
    end
    close(f)

    return rank, labels
end

function exportEmissionDevRate(years=[], nations=[], tag="", outputFile=""; mode=["cf"], maxr=0.5, minr=-0.5, nspan=128, descend=false, empty=false)

    global yr_list, nat_list, cat_list, gisRegList
    global ieRegGIS, deRegGIS, cfRegGIS, ieRegPcGIS, deRegPcGIS, cfRegPcGIS
    global ieRegDevGIS, deRegDevGIS, cfRegDevGIS, ieRegDevPcGIS, deRegDevPcGIS, cfRegDevPcGIS
    global ieRegDevRankGIS, deRegDevRankGIS, cfRegDevRankGIS, ieRegDevPcRankGIS, deRegDevPcRankGIS, cfRegDevPcRankGIS

    nc = length(cat_list)
    if length(years) == 0; years = yr_list elseif isa(years, Number); years = [years] end
    if length(nations) == 0; nations = nat_list elseif isa(nations, String); nations = [nations] end

    spanval = Dict{Int, Dict{String, Array{Float64, 2}}}()
    spanvalpc = Dict{Int, Dict{String, Array{Float64, 2}}}()

    for y in years, n in nations, m in mode
        rl = gisRegList[y][n]
        if m == "ie"; gre, grepc = ieRegGIS[y][n], ieRegPcGIS[y][n]
        elseif m == "de"; gre, grepc = deRegGIS[y][n], deRegPcGIS[y][n]
        elseif m == "cf"; gre, grepc = cfRegGIS[y][n], cfRegPcGIS[y][n]
        end
        if !haskey(spanval, y); spanval[y] = Dict{String, Array{Float64, 2}}() end
        if !haskey(spanvalpc, y); spanvalpc[y] = Dict{String, Array{Float64, 2}}() end

        println(m)

        fns = rsplit(outputFile, '.', limit=2)
        filename = replace(fns[1],"YEAR_"=>string(y)*"_") * "_" * uppercase(m) * "." * fns[2]
        gred, rank, spanval[y][n] = exportEmissionDevTable(replace(filename,"_OvPcTag"=>"_overall"), tag, rl, gre, maxr, minr, nspan, descend, empty)
        gredpc, rankpc, spanvalpc[y][n] = exportEmissionDevTable(replace(filename,"_OvPcTag"=>"_percap"), tag, rl, grepc, maxr, minr, nspan, descend, empty)

        if m == "ie"
            if !haskey(ieRegDevGIS, y); ieRegDevGIS[y] = Dict{String, Array{Float64, 2}}() end
            if !haskey(ieRegDevPcGIS, y); ieRegDevPcGIS[y] = Dict{String, Array{Float64, 2}}() end
            if !haskey(ieRegDevRankGIS, y); ieRegDevRankGIS[y] = Dict{String, Array{Int, 2}}() end
            if !haskey(ieRegDevPcRankGIS, y); ieRegDevPcRankGIS[y] = Dict{String, Array{Int, 2}}() end
            ieRegDevGIS[y][n], ieRegDevPcGIS[y][n] = gred, gredpc
            ieRegDevRankGIS[y][n], ieRegDevPcRankGIS[y][n] = rank, rankpc
        elseif m == "de"
            if !haskey(deRegDevGIS, y); deRegDevGIS[y] = Dict{String, Array{Float64, 2}}() end
            if !haskey(deRegDevPcGIS, y); deRegDevPcGIS[y] = Dict{String, Array{Float64, 2}}() end
            if !haskey(deRegDevRankGIS, y); deRegDevRankGIS[y] = Dict{String, Array{Int, 2}}() end
            if !haskey(deRegDevPcRankGIS, y); deRegDevPcRankGIS[y] = Dict{String, Array{Int, 2}}() end
            deRegDevGIS[y][n], deRegDevPcGIS[y][n] = gred, gredpc
            deRegDevRankGIS[y][n], deRegDevPcRankGIS[y][n] = rank, rankpc
        elseif m == "cf"
            if !haskey(cfRegDevGIS, y); cfRegDevGIS[y] = Dict{String, Array{Float64, 2}}() end
            if !haskey(cfRegDevPcGIS, y); cfRegDevPcGIS[y] = Dict{String, Array{Float64, 2}}() end
            if !haskey(cfRegDevRankGIS, y); cfRegDevRankGIS[y] = Dict{String, Array{Int, 2}}() end
            if !haskey(cfRegDevPcRankGIS, y); cfRegDevPcRankGIS[y] = Dict{String, Array{Int, 2}}() end
            cfRegDevGIS[y][n], cfRegDevPcGIS[y][n] = gred, gredpc
            cfRegDevRankGIS[y][n], cfRegDevPcRankGIS[y][n] = rank, rankpc
        end
    end

    return spanval, spanvalpc
end

function exportEmissionDevTable(outputFile, tag, reg_list, gre, maxr, minr, nspan, descend, empty)
    # this function is for [exportEmissionDiffRate]

    global cat_list
    nc, ngr = length(cat_list), length(reg_list)

    # calculate difference rates
    avg = mean(gre, dims=1)
    gred = zeros(ngr, nc)
    for i=1:nc; if avg[i]>0; gred[:,i] = (gre[:,i].-avg[i])/avg[i] end end

    # grouping by ratios; ascending order
    span = [(maxr-minr)*(i-1)/(nspan-2)+minr for i=1:nspan-1]
    spanval = zeros(Float64, nspan, nc)
    for i=1:nc
        spanval[1:end-1,i] = span[:].*avg[i].+avg[i]
        spanval[end,i] = spanval[end-1,i]
    end

    rank = zeros(Int, ngr, nc)
    for i = 1:ngr, j = 1:nc
        if gred[i,j]>=maxr; rank[i,j] = nspan
        else rank[i,j] = findfirst(x->x>gred[i,j],span)
        end
    end
    # for descending order, if "descend == true".
    if descend; for i = 1:ngr, j = 1:nc; rank[i,j] = nspan - rank[i,j] + 1 end end

    # exporting difference table
    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f_sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")
    print(f, tag); for c in cat_list; print(f, f_sep, c) end; println(f)
    for i = 1:ngr
        print(f, reg_list[i])
        for j = 1:nc; print(f, f_sep, gred[i,j]) end
        println(f)
    end
    if empty
        # for not-covered NUTS data
    end
    close(f)

    # exporting difference group table
    fns = rsplit(outputFile, '.', limit=2)
    f = open(fns[1]*"_gr."*fns[2], "w")
    print(f, tag); for c in cat_list; print(f, f_sep,c) end; println(f)
    for i = 1:ngr
        print(f, reg_list[i])
        for j = 1:nc; print(f, f_sep, rank[i,j]) end
        println(f)
    end
    if empty
        # for not-covered NUTS data
    end
    close(f)

    return gred, rank, spanval
end

function exportWebsiteFiles(years=[], nations=[], path=""; mode=["cf"], rank=false, empty=false)

    global yr_list, nat_list, cat_list, gisCoord, gisRegList, gisRegLabel, gisRegProv, gisPop, gisSample, gisAvgExp
    global ieRegGIS, deRegGIS, cfRegGIS, ieRegRankGIS, deRegRankGIS, cfRegRankGIS
    global ieRegDevPcGIS, deRegDevPcGIS, cfRegDevPcGIS, ieRegDevPcRankGIS, deRegDevPcRankGIS, cfRegDevPcRankGIS

    if length(years) == 0; years = yr_list elseif isa(years, Number); years = [years] end
    if length(nations) == 0; nations = nat_list elseif isa(nations, String); nations = [nations] end
    nc = length(cat_list)

    for y in years, n in nations, m in mode
        gr, lab, crd, pr, gp, ave = gisRegList[y][n], gisRegLabel[y][n], gisCoord[y][n], gisRegProv[y][n], gisPop[y][n], gisAvgExp[y][n]
        nr = length(gr)

        if m == "ie"; gre, grer, gredpc, gredrpc = ieRegGIS[y][n], ieRegRankGIS[y][n], ieRegDevPcGIS[y][n], ieRegDevPcRankGIS[y][n]
        elseif m == "de"; gre, grer, gredpc, gredrpc = deRegGIS[y][n], deRegRankGIS[y][n], deRegDevPcGIS[y][n], deRegDevPcRankGIS[y][n]
        elseif m == "cf"; gre, grer, gredpc, gredrpc = cfRegGIS[y][n], cfRegRankGIS[y][n], cfRegDevPcGIS[y][n], cfRegDevPcRankGIS[y][n]
        end

        f_path = path * n *"/" * string(y) * "/" * uppercase(m) * "/"
        mkpath(f_path)

        # print center file
        f = open(f_path * "centers.csv", "w")
        println(f, "\"NO\",\"GID2CODE\",\"PNAME\",\"DNAME\",\"x\",\"y\"")
        cnt = 1
        for r in gr
            println(f, "\"", cnt, "\",\"", r, "\",\"", pr[r][2], "\",\"", lab[r], "\",\"", crd[r][1], "\",\"", crd[r][2], "\"")
            cnt += 1
        end
        close(f)

        # print english file
        f = open(f_path * "english.txt", "w")
        println(f, "KEY_CODE\tEN_NAME")
        for r in gr; println(f, r, "\t", lab[r], ", ", pr[r][2]) end
        close(f)

        # print english_name file
        f = open(f_path * "english_match.txt", "w")
        println(f, "KEY_CODE\tSTATE_CODE\tSTATE\tDISTRICT")
        for r in gr; println(f, r, "\t", pr[r][1], "\t", pr[r][2], "\t", lab[r]) end
        close(f)

        # print ALLP file
        f = open(f_path * "ALLP.txt", "w")
        println(f, "ALL\tALLP")
        ci = findfirst(x->gisCatLab[x]=="All", cat_list)
        for i = 1:nr
            r = gr[i]
            println(f, r, "\t", grer[i, ci])
        end
        close(f)

        # print CF files
        for j = 1:nc
            mkpath(f_path * "CFAV/")
            f = open(f_path * "CFAV/" * "CFAV_" * gisCatLab[cat_list[j]] * ".txt", "w")
            print(f, "KEY_CODE\tSTATE\tDISTRICT\tSTATE_NAME\tDISTRICT_NAME\tGENHH_ALLP\tGENHH_APPPC")
            if cat_list[j]=="Total" || cat_list[j]=="All"; println(f, "\tANEXPPC\tPOP")
            else println(f)
            end
            for i = 1:nr
                r = gr[i]
                tbidx = i
                print(f, r, "\t", pr[r][1], "\t", r, "\t", pr[r][2], "\t", lab[r],"\t")
                printfmt(f, "{:f}", gre[i,j]); print(f, "\t", gre[i,j]/gp[i])
                if cat_list[j]=="Total" || cat_list[j]=="All"; println(f, "\t",ave[i], "\t", convert(Int, gp[i]))
                else println(f)
                end
            end
            if empty
                # for not-covered NUTS data
            end
            close(f)

            mkpath(f_path * "CFAC/")
            f = open(f_path * "CFAC/" * "CFAC_" * gisCatLab[cat_list[j]] * ".txt", "w")
            println(f, "KEY_CODE\tSTATE\tDISTRICT\tSTATE_NAME\tDISTRICT_NAME\tGENHH_ALLPPC")
            for i = 1:nr
                r = gr[i]
                print(f, r, "\t", pr[r][1], "\t", r, "\t", pr[r][2], "\t", lab[r], "\t")
                if rank; println(f, gredrpc[i,j]) else println(f, gredpc[i,j]) end
            end
            if empty
                # for not-covered NUTS data
            end
            close(f)
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
