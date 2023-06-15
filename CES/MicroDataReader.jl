# SPDX-FileCopyrightText: © 2021 Jemyung Lee <jemyung81@gmail.com>
# SPDX-License-Identifier: GPL-3.0

module MicroDataReader

# Developed date: 17. Mar. 2021
# Last modified date: 15. Jun. 2023
# Subject: Household consumption expenditure survey microdata reader
# Description: read consumption survey microdata and store household, member, and expenditure data
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using Statistics

mutable struct expenditure
    code::String        # product or service item code
    value::Float64      # total consumption monetary value
    quantity::Float64   # total consumption physical quantity
    valUnit::String     # total value's currency unit; ex. 'USD' or 'EUR'
    qntUnit::String     # total quantity's unit; ex. 'kg', 'litre', or 'number'
    period::Int16       # total value's consuming period; generally 7, 30, or 365 days

    expenditure(c, v=0, q=0, vu="", qu="", p=0) = new(c, v, q, vu, qu, p)
end

mutable struct member
    hhid::String        # household id
    age::Int16          # age
    gen::Int8           # gender: [1]male, [2]female, [9]not available
    nat::String         # nationality (A3)
    mar::Int8           # marital status: [1]single, [2]married, [9]not available
    edu::Int8           # education level: [1]non literate, [2]literate, [3]below primary, [4]primary, [5]below secondary, [6]secondary,
                        #                  [7]higher secondary, [8]tertiary, [8]graduate, [9]post gradutate/above, [10]not specified, [99]not available
    rel::Int8           # relationship with head: (of reference person and/or of the spouse) [1]reference person, [2]spouse or partner,
                        #                         [3]child, [4]parent, [5]relative, [6]no family relationship, [9]not available,
    occ::String         # occupation
    inc::Float64        # income amount
    incUnit::String     # income unit: ex. 'USD/month'

    member(id, ag=0, ge=0, na="", ma=0, ed=0, re=0, oc="", in=0, iu="") = new(id, ag, ge, na, ma, ed, re, oc, in, iu)
end

mutable struct household
    hhid::String        # household identification no.
    date::String        # survey date: YYYYMM, ex) 2015. Jan = '201501'
    province::String    # province/state code
    district::String    # district/city code
    regtype::String     # region type: urban/rural or dense/intermediate/sparce
    size::Int16         # household size
    age::Int16          # head's age
    rel::String         # head's religion
    occ::String         # head's occupation
    edu::String         # head's education level
    totexp::Float64	    # hh's total expenditure
    totinc::Float64     # hh's total income
    totexppc::Float64	# hh's total expenditure per capita
    totincpc::Float64   # hh's total income per capita
    unit::String        # currency unit of expenditure and income: ex. 'PPP_USD/month'
    popwgh::Float64     # population weight: (urban/rural) province/state/distric/city population represented by one member

    aggexp::Float64     # aggregate total expenditure (for validation)

    members::Array{member,1}        # household member(s)
    expends::Array{expenditure,1}   # consumed product or sevice items

    group::String       # (survey) group

    household(hi,da="",pr="",di="",rt="",sz=0,ag=0,rl="",oc="",ed="",te=0,ti=0,tep=0,tip=0,ut="",pw=0,ae=0,mm=[],ex=[],sv="") = new(hi,da,pr,di,rt,sz,ag,rl,oc,ed,te,ti,tep,tip,ut,pw,ae,mm,ex,sv)
end

mutable struct commodity
    code::String        # commodity code in the survey data
    sector::String      # commodity's sector or label
    category::String    # main category, ex) Food, Electricity, Gas, Other energy, Public transport, Private transport, Medical care, Education, Consumable goods, Durable goods, Other services
    subCategory::String # sub-category, ex) Food related: Grain, Vegetable, Fruit, Dairy, Beef, Pork, Poultry, Other meat, Fish, Alcohol, Other beverage, Confectionery, Restaurant, Other food, etc
                        #                   Energy-related: Electricity, Gas, Wood, Dung cake, Kerosene, Coal, Petrol, Diesel, Biogas, Other fuel, etc.
                        #                   Transport-related: Road (private), Road (public), Rail, Air, Water, Other, etc.
    unit_cur::String        # commodity's currency unit, ex) USD, EUR, etc.
    unit_qnt::String        # commodity's physical unit, ex) liter, kg, m^3, etc.
    coicop::String      # commodity corresponding COICOP code
    group::String       # corresponding (survey) group

    commodity(cod, sec="", cat="", subcat="", unt_cr = "", unt_qt = "", coi="", gr="") = new(cod, sec, cat, subcat, unt_cr, unt_qt, coi, gr)
end

global households = Dict{Int, Dict{String, Dict{String, household}}}()  # household dict: {year, {nation A3, {hhid, household}}}
global sectors = Dict{Int, Dict{String, Dict{String, commodity}}}()     # expenditure sector: {year, {nation A3, {code, commodity}}}
global hh_list = Dict{Int, Dict{String, Array{String, 1}}}()            # hhid list: {year, {nation A3, {hhid}}}
global sc_list = Dict{Int, Dict{String, Array{String, 1}}}()            # commodity code list: {year, {nation A3, {code}}}
global gr_list = Dict{Int, Dict{String, Array{String, 1}}}()            # group (survey type) list: {year, {nation A3, {group}}}

global regions = Dict{Int, Dict{String, Dict{String, String}}}()        # region code-name: {year, {nation A3, {code, region}}}
global prov_list = Dict{Int, Dict{String, Array{String, 1}}}()          # province code list: {year, {nation A3, {code}}}
global dist_list = Dict{Int, Dict{String, Array{String, 1}}}()          # district code list: {year, {nation A3, {code}}}
global dist_prov = Dict{Int, Dict{String, Dict{String, String}}}()      # district's province: {year, {nation A3, {district code, province code}}}
global region_modified = Dict{Int, Dict{String, Dict{String,String}}}() # revising regoin code list: {year, {nation A3, {original code, modified_code}}}

global pops = Dict{Int, Dict{String, Dict{String, Float64}}}()          # population: {year, {nation, {region_code, population}}}
global pops_ur = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}()   # urban/rural population: {year, {nation, {region_code, (urban, rural)}}
global pop_wgh = Dict{Int, Dict{String, Dict{String, Float64}}}()       # population weight: {year, {nation, {region_code, weight}}}
global pop_ur_wgh = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}()    # urban/rural population weight: {year, {nation, {region_code, (urban, rural)}}
global pop_gr_wgh = Dict{Int, Dict{String, Dict{String, Array{Float64, 1}}}}()    # population weight by group(or survey type): {year, {nation, {region_code, {group}}}}

global hh_curr = Dict{Int, Dict{String, Array{String, 1}}}()            # currency unit for household values (income or expenditure): {year, {nation, {currency}}}
global hh_period = Dict{Int, Dict{String, Array{String, 1}}}()          # period for household values (income or expenditure): {year, {nation, {period}}}
global exp_curr = Dict{Int, Dict{String, Array{String, 1}}}()           # currency unit for expenditure values: {year, {nation, {currency}}}
global exp_period = Dict{Int, Dict{String, Array{String, 1}}}()         # period for expenditure values: {year, {nation, {period}}}
global exchange_rate = Dict{String, Dict{String, Float64}}()            # exchange rate: {currency("original/target", ex."EUR/USD", {year(string), converting_rate}}

global expMatrix = Dict{Int, Dict{String, Array{Float64, 2}}}()         # expenditure matrix: {year, {nation, {hhid, commodity}}}
global qntMatrix = Dict{Int, Dict{String, Array{Float64, 2}}}()         # qunatity matrix: {year, {nation, {hhid, commodity}}}

global cpi_list = Dict{Int, Dict{String, Array{String, 1}}}()           # Consumption price indexes: {year, {nation, {COICOP_category}}}
global cpis = Dict{Int, Dict{String, Dict{String, Float64}}}()          # Consumption price indexes: {year, {nation, {COICOP_category, CPI}}}
global scl_rate = Dict{Int, Dict{String, Dict{String, Float64}}}()      # CPI scaling rate: {year, {nation, {HBS code, rate}}}

global pr_unts = Dict(1=>"day", 7=>"week", 30=>"month", 365=>"year")    # period units
global pr_scl = Dict("year"=>365.0, "month"=>30.0, "week"=>7.0, "day"=>1.0, "annual"=>365, "monthly"=>30.0, "weekly"=>7.0, "daily"=>1.0)    # period scales

function appendCommoditySectorData(year, nation, cmm_data)
    global sectors, sc_list

    push!(sc_list[year][nation], cmm_data[1])
    sectors[year][nation][cmm_data[1]] = commodity(cmm_data[1], cmm_data[2], cmm_data[3], "", cmm_data[4], cmm_data[5], cmm_data[6], cmm_data[7])
end

function appendExpenditureData(year, nation, id, exp_value)
    global households
    hh = households[year][nation][id]

    push!(hh.expends, expenditure(exp_value[1], exp_value[2], exp_value[3], exp_value[4], exp_value[5], exp_value[6]))
end

function appendMemberData(year, nation, id, mm_value)
    global households
    hh = households[year][nation][id]

    push!(hh.members, member(id, mm_value[1],mm_value[2],mm_value[3],mm_value[4],mm_value[5],mm_value[6],mm_value[7],mm_value[8],mm_value[9]))
end

function reviseHouseholdData(year, nation, id, hh_value)
    global households
    hh = households[year][nation][id]

    if length(id)>0 && length(hh_value)==18
        da = hh_value[1]; if length(da)>0; hh.date = da end
        pr = hh_value[2]; if length(pr)>0; hh.province = pr end
        di = hh_value[3]; if length(di)>0; hh.district = di end
        rt = hh_value[4]; if length(rt)>0; hh.regtype = rt end
        sz = hh_value[5]; if sz>0; hh.size = sz end
        ag = hh_value[6]; if ag>0; hh.age = ag end
        rl = hh_value[7]; if length(rl)>0; hh.rel = rl end
        oc = hh_value[8]; if length(oc)>0; hh.occ = oc end
        ed = hh_value[9]; if length(ed)>0; hh.edu = ed end
        te = hh_value[10]; if te>0; hh.totexp = te end
        ti = hh_value[11]; if ti>0; hh.totinc = ti end
        tep = hh_value[12]; if tep>0; hh.totexppc = tep end
        tip = hh_value[13]; if tip>0; hh.totincpc = tip end
        ut = hh_value[14]; if length(ut)>0; hh.unit = ut end
        pw = hh_value[15]; if pw>0; hh.popwgh = pw end
        mm = hh_value[16]; if length(mm)>0; hh.members = mm end
        ex = hh_value[17]; if length(ex)>0; hh.expends = ex end
        st = hh_value[18]; if length(st)>0; hh.group = st end
    else println("HHID is absent or incomplete HH values.")
    end
end

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

function readIndexFile(year, nation, indexFile; err_display=false)

    titles = Array{String, 1}()
    idxs = Array{Array{String, 1}, 1}()         # essential index values
    idxs_opt = Array{Array{String, 1}, 1}()     # optional index values
    idxstart = idxend = metaend = titlend = optchk = yrchk = natchk = false

    val_sep = getValueSeparator(indexFile)
    f = open(indexFile)
    for l in eachline(f)
        if length(strip(l, ['=', ' ', '\t'])) == 0; if !idxstart; idxstart = true else idxend = true end
        elseif length(strip(l, ['-', ' ', '\t'])) == 0 && idxstart
            if !metaend; metaend = true
            elseif !titlend; titlend = true
            elseif !optchk; optchk = true; titlend = false
            end
        elseif !titlend
            l = strip(l)
            if l[1] != '[' && l[end] != ']'; titles = strip.(split(l, val_sep)) end
        elseif idxstart && !idxend && titlend
            s = strip.(split(l, val_sep))
            tag = lowercase(filter(x->x!=':', s[1]))
            if !metaend
                if tag == "year" && parse(Int, s[2]) == year; yrchk = true
                elseif tag == "nation" && lowercase(s[2]) == lowercase(nation); natchk = true
                elseif tag == "a3" && uppercase(s[2]) == nation; natchk = true
                end
            elseif !optchk
                push!(idxs, [tag; s[2:end]])
                if length(s) < 5; push!(idxs[end], "") end
            else
                push!(idxs_opt, [tag; s[2:end]])
                if length(s) < 5; push!(idxs_opt[end], "") end
            end
        # else println("No applying condition: ", l)
        end
    end
    if err_display
        if !yrchk; println("Year does not match: ", year, ", ", indexFile) end
        if !natchk; println("Nation does not match: ", nation, ", ", indexFile) end
        if !idxstart || !idxend || !metaend || !titlend; println("Wrong index file format: ", indexFile) end
    end
    close(f)

    return yrchk, natchk, idxs, idxs_opt
end

function readPopulation(year, nation, populationFile)

    global pops, pops_ur

    if !haskey(pops, year); pops[year] = Dict{String, Dict{String, Float64}}() end
    if !haskey(pops[year], nation); pops[year][nation] = Dict{String, Float64}() end
    if !haskey(pops_ur, year); pops_ur[year] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
    if !haskey(pops_ur[year], nation); pops_ur[year][nation] = Dict{String, Tuple{Float64, Float64}}() end
    pop = pops[year][nation]
    pop_ur = pops_ur[year][nation]

    yrchk, natchk, idxs = readIndexFile(year, nation, populationFile)
    for s in idxs
        pop[s[1]] = parse(Float64, s[3])
        if length(s)==5 && length(s[4])>0 && length(s[5])>0
            pop_ur[s[1]] = tuple(parse(Float64, s[4]), parse(Float64, s[5]))
        end
    end
end

function readRegion(year, nation, regionFile; region_revised_file = "")

    global regions, prov_list, dist_list, dist_prov, region_modified

    if !haskey(regions, year); regions[year] = Dict{String, Dict{String, String}}() end
    if !haskey(regions[year], nation); regions[year][nation] = Dict{String, String}() end
    if !haskey(prov_list, year); prov_list[year] = Dict{String, Array{String, 1}}() end
    if !haskey(dist_list, year); dist_list[year] = Dict{String, Array{String, 1}}() end
    if !haskey(dist_prov, year); dist_prov[year] = Dict{String, Dict{String, String}}() end
    if !haskey(dist_prov[year], nation); dist_prov[year][nation] = Dict{String, String}() end

    rg = regions[year][nation]
    dp = dist_prov[year][nation]
    pl = Array{String, 1}()
    dl = Array{String, 1}()

    yrchk, natchk, idxs = readIndexFile(year, nation, regionFile)
    for s in idxs
        rg[s[1]] = s[6]
        dp[s[1]] = s[3]
        if !haskey(rg, s[3]); rg[s[3]] = s[4] end
        if !(s[1] in dl); push!(dl, s[1]) end
        if !(s[3] in pl); push!(pl, s[3]) end
    end

    prov_list[year][nation] = sort(pl)
    dist_list[year][nation] = sort(dl)

    if length(region_revised_file) > 0
        if !haskey(region_modified, year); region_modified[year] = Dict{String, Dict{String, String}}() end
        if !haskey(region_modified[year], nation); region_modified[year][nation] = Dict{String, String}() end

        yrchk, natchk, idxs = readIndexFile(year, nation, region_revised_file)
        for s in idxs
            region_modified[year][nation][s[2]] = s[3]
            if lowercase(s[1]) in ["city", "district"]; filter!(x -> x != s[2], dist_list[year][nation])
            elseif lowercase(s[1]) in ["state", "province"]; filter!(x -> x != s[2], prov_list[year][nation])
            end
        end
    end
end

function readMicroData(year, nation, microdataPath, hhIdxFile, mmIdxFile, cmmIdxFile, expIdxFile; hhid_sec = "hhid",
                        skip_title = true, periodFiltering = false, ignoreException = true, region_modify = false,
                        visible = false, date_type = "standard")

    if visible; print(" (") end
    for idxfile in filter(x->length(x)>0, [hhIdxFile, mmIdxFile, cmmIdxFile, expIdxFile])
        # read microdata index file
        yrchk, natchk, idxs, ipdx_o = readIndexFile(year, nation, idxfile)

        # read microdata contents
        if idxfile == hhIdxFile
            if visible; print("hhs") end
            readHouseholdData(year, nation, [idxs;ipdx_o], microdataPath, hhid_sec = hhid_sec, skip_title = skip_title, region_modify = region_modify, date_type = date_type)
            # idxs: {sector, position, type, file, tag}
        elseif idxfile == mmIdxFile
            if visible; print(", mms") end
            readMemberData(year, nation, idxs, microdataPath, hhid_sec = hhid_sec, skip_title = skip_title)
            # idxs: {sector, position, type, file, tag}
        elseif idxfile == cmmIdxFile
            if visible; print(", cmm") end
            readCommoditySectors(year, nation, idxs)
            # idxs: {code, sector, category, currency_unit, quantity_unit, coicop}
        elseif idxfile == expIdxFile
            if visible; print(", exp") end
            readExpenditureData(year, nation, idxs, microdataPath, skip_title = skip_title, periodFiltering = periodFiltering, ignoreException = ignoreException)
            # idxs: {}
        end
    end
    if visible; print(" )") end
end

function readHouseholdData(year, nation, indices, microdataPath; hhid_sec = "hhid", skip_title = true, region_modify = false, date_type = "standard")
    # date_type: "standard"=YYYYMMDD, "YYMMDD", "DDMMYY" => "date" items are stored as a standard type

    global households, hh_list, hh_curr, hh_period, region_modified

    sectors = ["survey_date", "province/state", "district/city", "region_type", "hh_size", "head_age", "head_religion", "head_occupation", "head_education", "expenditure", "income", "exp_percap", "inc_percap", "currency_unit", "pop_weight", "agg_exp"]
    nsec = length(sectors)
    dat_sec = ["survey_date"]
    int_sec = ["hh_size", "head_age"]
    flo_sec = ["pop_weight"]
    acc_sec = ["expenditure", "income", "exp_percap", "inc_percap"]
    acc_scale = 1.0

    dt_chk = (date_type == "standard")
    dt_idxs = Dict("YYMMDD" => [1:2, 3:4, 5:6], "DDMMYY" => [5:6, 3:4, 1:2])
    if !dt_chk; yyidx, mmidx, ddidx = dt_idxs[date_type] end

    if !haskey(households, year)
        households[year] = Dict{String, Dict{String, household}}()
        hh_list[year] = Dict{String, Array{String, 1}}()
    end
    if !haskey(households[year], nation)
        households[year][nation] = Dict{String, household}()
        hh_list[year][nation] = Array{String, 1}()
    end
    if !haskey(hh_curr, year); hh_curr[year] = Dict{String, Array{String, 1}}() end
    if !haskey(hh_period, year); hh_period[year] = Dict{String, Array{String, 1}}() end
    if !haskey(hh_curr[year], nation); hh_curr[year][nation] = Array{String, 1}() end
    if !haskey(hh_period[year], nation); hh_period[year][nation] = Array{String, 1}() end

    # analyze index data
    mdFiles = Dict{String, Dict{String, Tuple{Int, String, String}}}()  # {microdata_file, {data_sector, {position, type}}}
    for idx in indices      # indices: {sector, position, type, file, tag}
        if length(idx[2])>0 && length(idx[4])>0
            mf = idx[4]
            if !haskey(mdFiles, mf); mdFiles[mf] = Dict{String, Tuple{Int, String, String}}() end
            mdFiles[mf][idx[1]] = (parse(Int, idx[2]), idx[3], idx[5])      # (position, type, tag)
        elseif length(idx[3])>0 && lowercase(idx[1]) == "currency_unit"
            currency, period = string.(strip.(split(idx[3], '/', limit = 2)))
            if '/' in period
                period, scale = string.(strip.(split(period, '/')))
                acc_scale = parse(Float64, scale)
            elseif '*' in period
                period, scale = string.(strip.(split(period, '*')))
                acc_scale = 1.0 / parse(Float64, scale)
            end
            if !(currency in hh_curr[year][nation]); push!(hh_curr[year][nation], currency) end
            if !(period in hh_period[year][nation]); push!(hh_period[year][nation], period) end
        end
    end
    mfs = sort(collect(keys(mdFiles)))
    for mf in mfs; if !haskey(mdFiles[mf], hhid_sec); println(mf, "does not contain HHID sector.") end end
    hh_unit = hh_curr[year][nation][1] * "/" * hh_period[year][nation][1]

    # read household data: all the microdata files should contain HHID values
    if region_modify; reg_mod = region_modified[year][nation] end
    for mf in mfs
        mfd = mdFiles[mf]
        hhid_pos = mfd[hhid_sec][1]
        hhid_tag = mfd[hhid_sec][3]
        mdf_sep = getValueSeparator(mf)
        f = open(microdataPath * mf)
        if skip_title; readline(f) end  # skip title line
        for l in eachline(f)
            s = strip.(split(l, mdf_sep))
            hhid = hhid_tag * s[hhid_pos]
            if !haskey(households[year][nation], hhid)
                households[year][nation][hhid] = household(hhid)
                push!(hh_list[year][nation], hhid)
            end
            hh_vals = ["","","","",0,0,"","","",0,0,0,0,"",0,[],[],""]

            for i = 1:nsec
                c = sectors[i]
                if haskey(mfd, c)
                    val = s[mfd[c][1]]
                    if length(val) > 0
                        if c in int_sec; hh_vals[i] = parse(Int, val)
                        elseif c in acc_sec; hh_vals[i] = parse(Float64, val) * acc_scale
                        elseif c in flo_sec; hh_vals[i] = parse(Float64, val)
                        elseif !dt_chk && c in dat_sec
                            if length(val) == 5 && date_type == "DDMMYY"; val = "0" * val end
                            yy, mm, dd = val[yyidx], val[mmidx], val[ddidx]
                            if length(yy) == 2; yy = "20" * yy end
                            hh_vals[i] = yy * mm * dd
                        else hh_vals[i] = val
                        end
                    end
                end
            end
            hh_vals[14] = hh_unit

            if region_modify
                if haskey(reg_mod, hh_vals[2]); hh_vals[2] = reg_mod[hh_vals[2]] end
                if haskey(reg_mod, hh_vals[3]); hh_vals[3] = reg_mod[hh_vals[3]] end
            end
            reviseHouseholdData(year, nation, hhid, hh_vals)
        end
        close(f)
    end
end

function filterData(year, nation; group = false, region = "district", quantity = false, filter_region = true, filter_hhs = true, filter_exp = true)
    # filter_region: filter regions correspond to households' region
    # filter_hhs: filter households correspond to region list
    # filter_exp: filter expenditure matrix rows correspond to houseohld list

    global households, hh_list, gr_list, prov_list, dist_list, expMatrix, qntMatrix
    hhs, hl, em = households[year][nation], hh_list[year][nation], expMatrix[year][nation]
    if quantity; qm = qntMatrix[year][nation] end

    if filter_region
        prl = filter(p -> findfirst(x -> hhs[x].province == p, hl) != nothing, prov_list[year][nation])
        dsl = filter(d -> findfirst(x -> hhs[x].district == d, hl) != nothing, dist_list[year][nation])
        if group
            for g in gr_list[year][nation]
                filter!(p -> findfirst(x -> hhs[x].province == p && hhs[x].group == g, hl) != nothing, prl)
                filter!(d -> findfirst(x -> hhs[x].district == d && hhs[x].group == g, hl) != nothing, dsl)
            end
        end
        prov_list[year][nation], dist_list[year][nation] = prl, dsl
    end

    if filter_hhs
        if region == "province"; hidx = findall(h -> hhs[h].province in prl, hl) end
        if region == "district"; hidx = findall(h -> hhs[h].district in dsl, hl) end
        hh_list[year][nation]  = hl[hidx]

        if filter_exp
            expMatrix[year][nation] = em[hidx, :]
            if quantity; qntMatrix[year][nation] = qm[hidx, :] end
        end
    end
end

function findLostRegion(year, nation)

    global households, hh_list, prov_list, dist_list

    hhs, hhl, prl, dsl = households[year][nation], hh_list[year][nation], prov_list[year][nation], dist_list[year][nation]

    lost_prs, lost_dss = Array{String, 1}(), Array{String, 1}()
    for h in hhl
        pr, ds = hhs[h].province, hhs[h].district
        if !(pr in prl) && !(pr in lost_prs); push!(lost_prs, pr) end
        if !(ds in dsl) && !(ds in lost_dss); push!(lost_dss, ds) end
    end

    if length(lost_prs) > 0
        println("\n[Lost province(s)]")
        for lpr in lost_prs; println(lpr) end
    end
    if length(lost_dss) > 0
        println("\n[Lost district(s)]")
        for lds in lost_dss; println(lds) end
    end
end

function readMemberData(year, nation, indices, microdataPath; hhid_sec = "hhid", skip_title = true)

    global households

    sectors = ["age", "gender", "nationality", "head_relation", "marital_status", "education_level", "occupation", "income", "income_unit"]
    nsec = length(sectors)
    int_sec = ["age", "gender", "marital_status", "education_level", "head_relation"]
    acc_sec = ["income"]
    acc_scale = 1.0

    # analyze index data
    mdFiles = Dict{String, Dict{String, Tuple{Int, String, String}}}()  # {microdata_file, {data_sector, {position, type tag}}}
    for idx in indices      # indices: {sector, position, type, file, tag}
        if length(idx[2])>0 && length(idx[4])>0
            mf = idx[4]
            if !haskey(mdFiles, mf); mdFiles[mf] = Dict{String, Tuple{Int, String, String}}() end
            mdFiles[mf][idx[1]] = (parse(Int, idx[2]), idx[3], idx[5])      # (position, type, tag)
        elseif length(idx[3])>0 && lowercase(idx[1]) == "currency_unit"
            period = string.(strip.(split(idx[3], '/', limit = 2)))[end]
            if '/' in period; acc_scale = parse(Float64, string.(strip.(split(period, '/')))[end])
            elseif '*' in period; acc_scale = 1.0 / parse(Float64, string.(strip.(split(period, '/')))[end])
            end
        end
    end

    mfs = sort(collect(keys(mdFiles)))
    for mf in mfs; if !haskey(mdFiles[mf], hhid_sec); println(mf, "does not contain HHID sector.") end end

    # read household member data: all the microdata files should contain HHID values
    for mf in mfs
        mfd = mdFiles[mf]
        hhid_pos = mfd[hhid_sec][1]
        hhid_tag = mfd[hhid_sec][3]
        mdf_sep = getValueSeparator(mf)
        f = open(microdataPath * mf)
        if skip_title; readline(f) end  # skip title line
        for l in eachline(f)
            s = strip.(split(l, mdf_sep))
            hhid = hhid_tag * s[hhid_pos]
            mm_vals = [0, 0, "", 0, 0, 0, "", 0, ""]
            for i = 1:nsec
                c = sectors[i]
                if c in mfd
                    val = s[mfd[c][1]]
                    if lenght(val) > 0
                        if c in int_sec; mm_vals[i] = parse(Int, val)
                        elseif c in acc_sec; mm_vals[i] = parse(Float64, val) * acc_scale
                        else mm_vals[i] = val
                        end
                    end
                end
            end
            appendMemberData(year, nation, hhid, mm_vals)
        end
        close(f)
    end
end

function readCommoditySectors(year, nation, indices)
    # indices: [Code, Sector, Main_category, Unit_currency, Unit_quantity, COICOP, Survey]

    global sectors, sc_list

    if !haskey(sectors, year)
        sectors[year] = Dict{String, Dict{String, commodity}}()
        sc_list[year] = Dict{String, Array{String, 1}}()
    end
    if !haskey(sectors[year], nation)
        sectors[year][nation] = Dict{String, commodity}()
        sc_list[year][nation] = Array{String, 1}()
    end

    sec = sectors[year][nation]
    scl = sc_list[year][nation]
    for s in indices
        if length(s[1]) == 0; println("\nEmpty codes: ", s[1], ", ", s[3], ", ", sec[s[1]])
        elseif !(s[1] in scl); appendCommoditySectorData(year, nation, [s; ["" for i=1:(7 - length(s))]])
        else println("\nDuplicated codes: ", s[1], ", ", s[3], ", ", sec[s[1]])
        end
    end
end

function readExpenditureData(year, nation, indices, microdataPath; skip_title = true, periodFiltering = false, ignoreException = true)

    global households, hh_list, sectors, sc_list, exp_curr, exp_period, pr_unts
    if !haskey(exp_curr, year) exp_curr[year] = Dict{String, Array{String, 1}}() end
    if !haskey(exp_period, year) exp_period[year] = Dict{String, Array{String, 1}}() end
    if !haskey(exp_curr[year], nation) exp_curr[year][nation] = Array{String, 1}() end
    if !haskey(exp_period[year], nation) exp_period[year][nation] = Array{String, 1}() end

    cm = sectors[year][nation]
    sl = sc_list[year][nation]
    units = ["day", "week", "month", "year"]
    op_code = Dict{Tuple{String, String, String} , Tuple{String, String}}() # optional index list: {(microdata_file, category, code_position), (start_code, end_code)}

    # analyze index data
    # mdFiles: {microdata_file, {category, {data_tag, hhid_position, code_position, period(days), value_position, value_unit, quantity_position, quantity_unit, value_scale, quantity_scale}}}
    mdFiles = Dict{String, Dict{String, Tuple{String, Int, Int, Int, Int, String, Int, String, Float64, Float64}}}()
    mdFiles_op = Dict{String, Dict{String, Tuple{String, Int, Int, Int, Int, String, Int, String, Float64, Float64}}}()  # optional
    for idx in indices      # index: {category, hhid_position, code_position, period(days), value_position, value_unit, file, data_tag, quantity_position, quantity_unit}
        val_scale, qnt_scale = 1.0, 1.0
        if all(length.(idx[[1, 2, 3, 4, 7]]).>0) && (all(length.(idx[[5, 6]]).>0) || all(length.(idx[[9, 10]]).>0))
            if length(idx) in [7,8]
                li = findlast(x->tryparse(Float64, x) != nothing, idx)
                if li == 5; idx = [idx; (length(idx) == 7 ? ["", "0", ""] : ["0", ""])]
                elseif li == 7; idx = [idx[1:4]; ["0", ""] ;idx[5:end]]
                end
            end
            if length(idx) == 10
                idxs = [[idx[1], parse(Int, idx[2]), idx[3], parse(Int, idx[4])]; idx[5:end]]
                for i in [5, 9]; if length(idxs[i]) > 0; idxs[i] = parse(Int, idxs[i]); else idxs[i] = 0 end end
                mf = idxs[7]
                for i in [6, 10]
                    if '/' in idx[i]
                        idx[i], scale = string.(strip.(split(idx[i] , '/')))
                        i == 6 ? (val_scale = parse(Float64, scale)) : (qnt_scale = parse(Float64, scale))
                    elseif '*' in idx[i]
                        idx[i], scale = string.(strip.(split(idx[i], '*')))
                        i == 6 ? (val_scale /= parse(Float64, scale)) : (qnt_scale /= parse(Float64, scale))
                    end
                end
                if '_' in idx[3]    # optional indexing
                    idxs[3], st_cd, ed_cd = string.(strip.(split(idx[3], '_')))
                    op_code[(mf, idx[1], idxs[3])] = (st_cd, ed_cd)
                    mdf = mdFiles_op
                else mdf = mdFiles
                end
                idxs[3] = parse(Int, idxs[3])
                if !haskey(mdf, mf); mdf[mf] = Dict{String, Tuple{String, Int, Int, Int, Int, String, Int, String, Float64, Float64}}() end
                mdf[mf][idxs[1]] = tuple([idxs[8]; idxs[2:6]; idxs[9:end]; [val_scale]; [qnt_scale]]...)
            else println("Expenditure index content length error: ", year, ", ", nation, "\t", idx)
            end
            if !(idxs[6] in exp_curr[year][nation]); push!(exp_curr[year][nation], idxs[6]) end
            if !(pr_unts[idxs[4]] in exp_period[year][nation]); push!(exp_period[year][nation], pr_unts[idxs[4]]) end
        else println("Expenditure index lacks essential contents: ", year, ", ", nation, "\t", idx)
        end
    end
    mfs = sort(collect(keys(mdFiles)))

    # read expenditure data: all the microdata files should contain HHID and Code values
    exp_str_idx, mfd_str_idx = [1, 4, 5], [3, 6, 8]     # for string values
    n_str_idx = length(mfd_str_idx)

    exp_num_idx, mfd_num_idx = [2, 3], [5, 7]           # for numeric values
    n_num_idx = length(mfd_num_idx)

    for mf in mfs
        mfd = mdFiles[mf]
        sec = collect(keys(mfd))
        mdf_sep = getValueSeparator(mf)
        pre_mfd = []
        op_chk = haskey(mdFiles_op, mf)
        for sc in sec
            dup_chk = false
            op_chk = op_chk && haskey(mdFiles_op[mf], sc)
            if op_chk
                mfd_op = mdFiles_op[mf][sc]
                st_cd, ed_cd = op_code[(mf, sc, string(mfd_op[3]))]
                op_val_scale, op_qnt_scale = mfd_op[9], mfd_op[10]
            end
            for pre_sc in pre_mfd; if mfd[sc] == mfd[pre_sc]; dup_chk = true end end
            if !dup_chk
                push!(pre_mfd, sc)
                hhid_tag, hhid_pos, code_pos, val_scale, qnt_scale = mfd[sc][1], mfd[sc][2], mfd[sc][3], mfd[sc][9], mfd[sc][10]
                # {data_tag, hhid_position, code_position, period(days), value_position, value_unit, quantity_position, quantity_unit, value_scale, quantity_scale}
                f = open(microdataPath * mf)
                if skip_title; readline(f) end  # skip title line
                for l in eachline(f)
                    s = strip.(split(l, mdf_sep))
                    if !ignoreException || s[code_pos] in sl
                        hhid = hhid_tag * s[hhid_pos]
                        exp_vals = ["", 0, 0, "", "", 0]    #{code, value, quantity, value_unit, quantity_unit, period(days)}
                        # for string values
                        for i = 1:n_str_idx
                            val = mfd[sc][mfd_str_idx[i]]
                            if isa(val, Number); val = s[val] end
                            if length(val) > 0; exp_vals[exp_str_idx[i]] = val end
                        end
                        c_unit, q_unit = cm[exp_vals[1]].unit_cur, cm[exp_vals[1]].unit_qnt
                        if length(c_unit)>0; exp_vals[4] = c_unit end
                        if length(q_unit)>0; exp_vals[5] = q_unit end
                        # for numeric values
                        for i = 1:n_num_idx
                            val = mfd[sc][mfd_num_idx[i]] > 0 && length(s[mfd[sc][mfd_num_idx[i]]]) > 0 ? parse(Float64, s[mfd[sc][mfd_num_idx[i]]]) : 0
                            if val > 0; exp_vals[exp_num_idx[i]] = val * (mfd_num_idx[i] == 5 ? val_scale : qnt_scale) end
                        end
                        # for period-value
                        if mfd[sc][4] > 0; exp_vals[6] = mfd[sc][4]
                        else println("Error: period data does not exist. ", mf, ", ", sc, ", ", exp_vals[1])
                        end
                        if op_chk && st_cd <= exp_vals[1] <= ed_cd && mfd_op[4] > 0
                            exp_vals[4], exp_vals[5], exp_vals[6]  = mfd_op[6], mfd_op[8], mfd_op[4]
                            exp_vals[2] *= op_val_scale / val_scale
                            exp_vals[3] *= op_qnt_scale / qnt_scale
                        end
                        if periodFiltering
                            pr_str = pr_unts[exp_vals[6]]
                            pr_unit = split(exp_vals[5], '/')[end]
                            if pr_unit in units
                                if pr_unit != pr_str; exp_vals[3] = 0 end
                                exp_vals[5] = replace(replace(exp_vals[5], pr_str=>""), "/"=>"")
                            end
                        end
                        # append data
                        if !(exp_vals[2] == exp_vals[3] == 0); appendExpenditureData(year, nation, hhid, exp_vals) end
                    end
                end
                close(f)
            end
        end
    end
end

function buildExpenditureMatrix(year, nation; transpose = false, period = 365, quantity = false)
    # build an expenditure matrix as period (init: 365) days consumption monetary values
    # [row]: household, [column]: commodity

    global households, hh_list, sc_list, expMatrix, qntMatrix, exp_period, pr_unts
    if !haskey(expMatrix, year); expMatrix[year] = Dict{String, Array{Float64, 2}}() end
    if quantity && !haskey(qntMatrix, year); qntMatrix[year] = Dict{String, Array{Float64, 2}}() end

    if haskey(pr_unts, period); exp_period[year][nation] = Array{String, 1}([pr_unts[period]])
    else println("Expenditure matrix's period $period is not listed in ", pr_unts)
    end

    hhs = households[year][nation]
    row = hh_list[year][nation]
    col = sc_list[year][nation]
    nr = length(row)
    nc = length(col)

    mat = zeros(Float64, nr, nc)

    if quantity; qmat = zeros(Float64, nr, nc) end
    rowErr = zeros(Int, nr)
    colErr = zeros(Int, nc)

    # make expenditure matrix
    for ri = 1:nr
        h = hhs[row[ri]]
        total = 0
        for he in h.expends
            if he.code in col
                ci = findfirst(x -> x==he.code, col)
                if he.value > 0
                    val = he.value
                    if he.period != period; val *= period / he.period end
                    mat[ri,ci] += val
                    total += val
                end
                if quantity && he.quantity > 0
                    qnt = he.quantity
                    if he.period != period; qnt *= period / he.period end
                    qmat[ri,ci] += qnt
                end
                if he.value <= 0 && he.quantity <= 0; rowErr[ri] += 1; colErr[ci] += 1 end
            end
        end
        h.aggexp = total
    end
    if transpose
        mat = transpose(mat)
        qmat = transpose(qmat)
        row, nr, rowErr, col, nc, colErr = col, nc, colErr, row, nr, rowErr
    end
    expMatrix[year][nation] = mat
    if quantity; qntMatrix[year][nation] = qmat end

    return mat, row, col, rowErr, colErr
end

function exchangeExpCurrency(year, exchangeYear, nation, org_curr, exRateFile; target_curr = "USD", hhs_exp = true, hhs_info = true, exp_mat = false)

    global households, hh_list, exchange_rate, exp_curr, hh_curr

    if org_curr == target_curr; return end

    hhs = households[year][nation]
    hhl = hh_list[year][nation]
    er = exchange_rate[target_curr * "/" * org_curr] = Dict{String, Float64}()
    if exp_mat; em = expMatrix[year][nation] end
    if all(exp_curr[year][nation] .== org_curr); exp_chk = true; exp_curr[year][nation] = Array{String, 1}([target_curr])
    else exp_chk = false; println("Expenditure's original currency ", exp_curr[year][nation], " mismatch with $org_curr")
    end
    if hhs_info && all(hh_curr[year][nation] .== org_curr); hh_chk = true; hh_curr[year][nation] = Array{String, 1}([target_curr])
    elseif hhs_info; hh_chk = false; println("HH info's original currency ", hh_curr[year][nation], " mismatch with $org_curr")
    end

    # read exchange rate
    f_sep = getValueSeparator(exRateFile)
    f = open(exRateFile)
    readline(f)
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        trg, org = split(s[2], '/')
        trg_scl, org_scl = tryparse(Float64, filter(isdigit, trg)), tryparse(Float64, filter(isdigit, org))
        if trg_scl != nothing; trg = filter(!isdigit, trg) else trg_scl = 1.0 end
        if org_scl != nothing; org = filter(!isdigit, org) else org_scl = 1.0 end
        if (trg, org) == (target_curr, org_curr); er[s[1]] = parse(Float64, s[3]) * trg_scl / org_scl
        elseif (trg, org) == (org_curr, target_curr); er[s[1]] = 1.0/parse(Float64, s[3]) / trg_scl * org_scl
        end
    end
    if length(er) == 0; println("No exchange data for $org_curr to $target_curr") end
    close(f)

    # exchange the expenditure currency
    yr = string(exchangeYear)
    if !haskey(er, yr)
        mms = filter(x->(length(x)==6 && x[1:4]==yr), collect(keys(er)))
        length(mms)>0 ? er[yr] = mean([er[mm] for mm in mms]) : println("Error: no $yr correspoding exchange rate.")
    end
    for hh in hhl
        if length(hhs[hh].date)==0; r = er[yr]
        else
            h = hhs[hh]
            lh = length(h.date); mmidx = 5:lh; yyidx = 1:4
            if haskey(er, h.date[mmidx]); r = er[h.date[mmidx]]
            elseif haskey(er, h.date[yyidx]); r = er[h.date[yyidx]]
            else println("Exchange rate error: no exchange rate data for ", h.date)
            end
        end
        if hhs_exp && exp_chk
            hhs[hh].aggexp *= r
            for he in hhs[hh].expends
                he.value *= r
                he.valUnit = target_curr
            end
        end
        if exp_mat && exp_chk; em[findfirst(x -> x == hh, hhl), :] *= r end
        if hhs_info && hh_chk
            hhs[hh].totexp *= r
            hhs[hh].totinc *= r
            hhs[hh].totexppc *= r
            hhs[hh].totincpc *= r
            if org_curr == strip(string(split(hhs[hh].unit, '/')[1])); hhs[hh].unit = replace(hhs[hh].unit, org_curr => target_curr) end
        end
    end
end

function convertAvgExpToPPP(year, nation, pppConvRate; inverse=false)
    # PPP rates: can be a file path that contains PPP rates, a constant value of
    #            a nation's currency to USD (normally) currency exchange rate (USD/A3), or a set of values of Dict[MMYY] or Dict[YY]

    global households, hh_list
    hhs = households[year][nation]
    hhl = hh_list[year][nation]

    # read converting rate from the recieved file if 'pppFile' is 'String'
    if typeof(pppConvRate) <:AbstractString
        ppp = Dict{String, Float64}()
        f_sep = getValueSeparator(pppConvRate)
        f = open(pppConvRate)
        readline(f)
        for l in eachline(f); s = split(l, f_sep); ppp[s[1]] = parse(Float64, s[2]) end
        close(f)
        if inverse; for x in collect(keys(ppp)); ppp[x] = 1/ppp[x] end end
        pppConvRate = ppp
    end

    if typeof(pppConvRate) <: Number
        for hh in hhl
            hhs[hh].totexp /= pppConvRate
            hhs[hh].totinc /= pppConvRate
            hhs[hh].totexppc /= pppConvRate
            hhs[hh].totincpc /= pppConvRate
            for mm in hhs[hh].members; mm.inc /= pppConvRate end
        end
    elseif typeof(pppConvRate) <: AbstractDict
        yr = string(year)
        if !(yr in collect(keys(pppConvRate)))
            rates = [pppConvRate[mm] for mm in filter(x->(length(x)==6 && x[1:4]==yr), collect(keys(pppConvRate)))]
            pppConvRate[yr] = sum(rates) / length(rates)
        end
        for hh in hhl
            if length(hhs[hh].date)==0; ppp = pppConvRate[yr]
            else
                lh = length(h.date); mmidx = 5:lh; yyidx = 1:4
                if haskey(pppConvRate, h.date[mmidx]); ppp = pppConvRate[h.date[mmidx]]
                elseif haskey(pppConvRate, h.date[yyidx]); ppp = pppConvRate[h.date[yyidx]]
                else println("PPP converting rate error: no exchange rate data for ", h.date)
                end
            end
            hhs[hh].totexp /= ppp
            hhs[hh].totinc /= ppp
            hhs[hh].totexppc /= ppp
            hhs[hh].totincpc /= ppp
            for mm in hhs[hh].members; mm.inc /= ppp end
        end
    end
end

function calculatePopWeight(year, nation, outputFile=""; ur_wgh = false, district=true, province=false, hhs_wgh = false, gr_wgh = false)
    # filter: [true] exclude households not correspond to the region list

    global regions, prov_list, dist_list, pops, pops_ur, pop_wgh, pop_ur_wgh, pop_gr_wgh
    global households, hh_list, gr_list

    prl, dsl = prov_list[year][nation], dist_list[year][nation]
    rl = Array{String, 1}()
    if province; append!(rl, prl) end
    if district; append!(rl, dsl) end

    hl = hh_list[year][nation]
    hhs = households[year][nation]
    pop = pops[year][nation]
    smp = Dict{String, Int}()           # Province sample size, {regin code, sample number}
    wgh = Dict{String, Float64}()       # Province population weight, {region code, weight}

    if ur_wgh
        pop_ur = pops_ur[year][nation]
        smp_ur = Dict{String, Tuple{Int,Int}}()             # Urban/rural province sample size, {regin code, (urban, rural)}
        wgh_ur = Dict{String, Tuple{Float64, Float64}}()    # Urban/rural province population weight, {region code, (urban, rural)}
    end
    if gr_wgh
        gl = gr_list[year][nation]
        ng = length(gl)
        smp_gr = Dict{String, Array{Int,1}}()               # regional sample size by group, {regin code, {group}}
        wgh_gr = Dict{String, Array{Float64, 1}}()          # population weight by group, {region code, {group}}
    end

    # count sample number
    for r in rl
        smp[r] = 0
        if ur_wgh
            smp_ur[r] = Tuple{Int, Int}()
            wgh_ur[r] = Tuple{Float64, Float64}()
        end
        if gr_wgh
            smp_gr[r] = zeros(Int, ng)
            wgh_gr[r] = zeros(Float64, ng)
        end
    end
    if province; for h in hl; smp[hhs[h].province] += hhs[h].size end end
    if district; for h in hl; smp[hhs[h].district] += hhs[h].size end end

    if ur_wgh
        for h in hl
            if hhs[h].regtype == "urban"
                if province; smp_ur[hhs[h].province][1] += hhs[h].size end
                if district; smp_ur[hhs[h].district][1] += hhs[h].size end
            elseif hhs[h].regtype == "rural"
                if province; smp_ur[hhs[h].province][2] += hhs[h].size end
                if district; smp_ur[hhs[h].district][2] += hhs[h].size end
            end
        end
    end
    if gr_wgh
        for i = 1:ng
            g = gl[i]
            for h in filter(x -> hhs[x].group == g, hl)
                if province; smp_gr[hhs[h].province][i] += hhs[h].size end
                if district; smp_gr[hhs[h].district][i] += hhs[h].size end
            end
        end
    end

    # calculate weights
    for r in rl; wgh[r] = pop[r]/smp[r] end
    if ur_wgh
        for r in rl, i = 1:2
            if pop_ur[r][i]>0 && smp_ur[r][i]>0; wgh_ur[r][i] = pop_ur[r][i] / smp_ur[r][i] end
        end
    end
    if gr_wgh
        for r in rl, i = 1:ng
            if pop[r]>0 && smp_gr[r][i]>0; wgh_gr[r][i] = pop[r] / smp_gr[r][i] end
        end
    end

    if !haskey(pop_wgh, year); pop_wgh[year] = Dict{String, Dict{String, Float64}}() end
    pop_wgh[year][nation] = wgh
    if ur_wgh
        if !haskey(pop_ur_wgh, year); pop_ur_wgh[year] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
        pop_ur_wgh[year][nation] = wgh_ur
    end
    if gr_wgh
        if !haskey(pop_gr_wgh, year); pop_gr_wgh[year] = Dict{String, Dict{String, Array{Float64, 1}}}() end
        pop_gr_wgh[year][nation] = wgh_gr
    end

    # assign household weight
    if hhs_wgh
        for h in hl
            hh = hhs[h]
            hh_reg = (district ? hh.district : (province ? hh.province : ""))

            if gr_wgh
                igr = findfirst(x -> x == hh.group, gl)
            elseif ur_wgh
                if hh.regtype in ["urban", "sparce"]; iur = 1
                elseif hh.regtype in ["rural", "intermediate"]; iur = 2
                elseif hh.regtype in ["dense"]; iur=3
                end
            end

            hh.popwgh = (gr_wgh ? wgh_gr[hh_reg][igr] : (ur_wgh ? wgh_ur[hh_reg][iur] : wgh[hh_reg]))
        end
    end

    # print population weights
    if length(outputFile)>0
        mkpath(rsplit(outputFile, '/', limit = 2)[1])
        f = open(outputFile, "w")
        print(f, "Code\tRegion\tPopulation_weight")
        if ur_wgh; print(f, "\tWeight_urban\tWeight_rural") end
        if gr_wgh; for g in gl; print(f, "\tWeight_", g) end end
        println(f)
        for r in rl
            print(f, r, "\t", regions[year][nation][r], "\t", pop_wgh[year][nation][r])
            if ur_wgh; for i = 1:2; print(f, "\t", pop_ur_wgh[year][nation][r][i]) end end
            if gr_wgh; for i = 1:ng; print(f, "\t", pop_gr_wgh[year][nation][r][i]) end end
            println(f)
        end
        close(f)
    end
end

function scalingExpByCPI(year, nation, targetYear, cpiCodeFile, statFile, sec_linkFile, reg_linkFile; top_sector = "", period="year", region="district", revHH=true, revMat=false)
    # scaling hh expenditure to the target year based on Consumer Price Index
    # period: "year" or "month"
    # region: "nation", "province" or "district"

    global hh_list, sc_list, prov_list, dist_list, dist_prov, households, pops, sectors, expMatrix
    global cpi_list, cpis, scl_rate

    for y in [year, targetYear]
        if !haskey(cpi_list, y); cpi_list[y] = Dict{String, Array{String, 1}}() end
        if !haskey(cpis, y); cpis[y] = Dict{String, Dict{String, Float64}}() end
        if !haskey(scl_rate, y); scl_rate[y] = Dict{String, Dict{String, Float64}}() end
    end

    y, n = year, nation
    hl, sl, hhs, pop = hh_list[y][n], sc_list[y][n], households[y][n], pops[y][n]
    ds_pr = dist_prov[y][n]
    ces_sectors = sectors[y][n]     # CES/HBS micro-data sectors: {code, commodity}
    if region == "district"; rl = dist_list[y][n]
    elseif region == "province"; rl = prov_list[y][n]
    end

    cpi_sec = Array{String, 1}()                    # CPI sector code list
    cpi_reg = Array{String, 1}()                    # CPI region code list
    cpi_per = Array{String, 1}()                    # CPI data period list: XXXX = year, XXXXMM = year(XXXX), month(MM)
    ces_cpi_link = Dict{String, String}()           # CES/HBS - CPI sector code matching: {micro-data, CPI}
    ces_cpi_region = Dict{String, String}()         # CES/HBS - CPI region code matching: {micro-data, CPI}
    cpi_sectors = Dict{String, String}()            # CPI sectors: {code, sector}
    cpi_vals = Array{Float64, 3}(undef,0,0,0)       # CPI values: {region, sector, period(year/month)}

    scl_rate_yr = Array{Float64, 2}(undef,0,0)      # annual average scaling rate = target_year_CPI/current_year_CPI: {region, sector}
    scl_rate_mth = Array{Float64, 3}(undef,0,0,0)   # monthly average scaling rate = target_year_CPI/current_year_CPI: {region, sector, month}

    # read CPI sectors
    val_sep = getValueSeparator(cpiCodeFile)
    f = open(cpiCodeFile)
    readline(f)
    for l in eachline(f)
        s = string.(strip.(split(l, val_sep)))
        cpi_sectors[s[1]] = s[2]
    end
    cpi_sec = sort(collect(keys(cpi_sectors)))
    close(f)
    cpi_list[y][n] = cpi_sec

    # read CPI-CES/HBS sector linkages
    if length(sec_linkFile) > 0
        val_sep = getValueSeparator(sec_linkFile)
        f = open(sec_linkFile)
        t_label = ["cpi_code", "cpi_sector", lowercase(n)*"_code", lowercase(n)*"_sector"]
        title = string.(strip.(split(lowercase(readline(f)), val_sep)))
        si = [findfirst(x -> x == tlb, title) for tlb in t_label]
        ces_cod_chk = Array{String, 1}()
        for l in eachline(f)
            s = string.(strip.(split(l, val_sep)))
            c_ces, c_cpi  = s[si[3]], s[si[1]]
            ces_cpi_link[c_ces] = c_cpi
            if !haskey(cpi_sectors, c_cpi); println("CPI sectors do not contain code $c_cpi") end
            if !haskey(ces_sectors, c_ces); println("CES/HBS sectors do not contain code $c_ces")
            elseif strip(ces_sectors[c_ces].sector) != strip(s[si[4]])
                println("CES/HBS sector does not match with $c_ces\t", s[si[4]], "\t", ces_sectors[c_ces].sector)
            end
            if !(c_ces in ces_cod_chk); push!(ces_cod_chk, c_ces) else println("Duplicated CES/HBS sector code: $c_ces") end
        end
        if !issubset(sort(ces_cod_chk), sort(sl)); println("CSI code matching list's CES/HBS codes doen't match with micro-data's") end
        close(f)
    else ces_sectors = Dict(sl .=> sl)
    end

    # read CPI values
    val_sep = getValueSeparator(statFile)
    f = open(statFile)
    cpi_per = string.(strip.(split(readline(f), val_sep)[5:end]))
    for l in eachline(f)
        reg_cod = string.(strip.(split(l, val_sep, limit=2)[1]))
        if !(reg_cod in cpi_reg); push!(cpi_reg, reg_cod) end
    end
    seek(f, 0)
    nr, ns, np = length(cpi_reg), length(cpi_sec), length(cpi_per)
    cpi_vals = zeros(Float64, nr, ns, np)
    readline(f)
    for l in eachline(f)
        s = string.(strip.(split(l, val_sep)))
        regidx = findfirst(x->x==s[1], cpi_reg)
        codidx = findfirst(x->x==s[3], cpi_sec)
        cpi_vals[regidx, codidx, :] = [parse(Float64, x) for x in s[5:end]]
    end
    close(f)

    # read CPI_CES/HBS region code linkanges
    if length(reg_linkFile) > 0
        val_sep = getValueSeparator(reg_linkFile)
        f = open(reg_linkFile)
        t_label = ["cpi_city_code", "hbs_city_code"]
        title = string.(strip.(split(lowercase(readline(f)), val_sep)))
        ti = [findfirst(x -> x == tlb, title) for tlb in t_label]
        ces_cod_chk = Array{String, 1}()
        for l in eachline(f)
            s = string.(strip.(split(l, val_sep)))
            r_cpi, r_ces  = s[ti[1]], s[ti[2]]
            ces_cpi_region[r_ces] = r_cpi
            if !(r_cpi in cpi_reg); println("CPI regions do not contain code ", r_cpi) end
            if !(r_ces in rl); println("CES/HBS regions do not contain code ", r_ces) end
            if !(r_ces in ces_cod_chk); push!(ces_cod_chk, r_ces) else println("Duplicated CES/HBS region code: $r_ces") end
        end
        if !issubset(sort(ces_cod_chk), sort(rl)); println("CPI region matching list's CES/HBS codes doen't match with micro-data's") end
        close(f)
    elseif nr == 1 && cpi_reg[1] == n; ces_cpi_region = Dict(rl .=> cpi_reg[1])
    elseif issubset(rl, cpi_reg); ces_cpi_region = Dict(rl .=> rl)
    else println("\nIncorrect CPI-CES region code linkanges: ", cpi_reg)
    end

    # assign CPIs
    yrs = unique(sort([p[1:4] for p in cpi_per]))
    mths = filter(x->length(x)==6, cpi_per)
    ny, nm = length(yrs), length(mths)
    cpi_vals_yr, cpi_vals_mth = zeros(Float64, nr, ns, ny), zeros(Float64, nr, ns, nm)
    for i = 1:np
        p = cpi_per[i]
        if length(p)==4
            yi = findfirst(x->x==p, yrs)
            for j=1:nr, k=1:ns; cpi_vals_yr[j, k, yi] = cpi_vals[j, k, i] end
        elseif length(p)==6
            mi = findfirst(x->x==p, mths)
            for j=1:nr, k=1:ns; cpi_vals_mth[j, k, mi] = cpi_vals[j, k, i] end
        else println("Period is not year nor month: ", p)
        end
    end
    for i = 1:nr, j = 1:ns, k = 1:ny
        if cpi_vals_yr[i, j, k] == 0
            cpi_vals_yr[i, j, k] = mean(filter(x->x>0, cpi_vals_mth[i, j, findall(x->x[1:4]==yrs[k], mths)]))
        end
    end

    for yr in [year, targetYear]
        yi = findfirst(x -> x == yr, yrs)
        if nr == 1 && cpi_reg[1] == n
            cpis[yr][n] = Dict(cpi_sec[j] => cpi_vals_yr[1, j, yi] for j = 1:ns)
        else
            cpi_nat, tot_pop = zeros(Float64, ns), 0
            for i = 1:nr
                p = pops[yr][n][cpi_ces_region[cpi_reg[i]]]
                for j = 1:ns; cpi_nat[j] += cpi_vals_yr[i, j, yi] * p end
                tot_pop += p
            end
            cpi_nat ./= tot_pop
            cpis[yr][n] = Dict(cpi_sec[j] => cpi_nat[j] for j = 1:ns)
        end
    end

    # calculate scaling rates
    cy, ty = string(year), string(targetYear)
    target_mths, current_mths = filter(x->x[1:4]==ty, mths), filter(x->x[1:4]==cy, mths)
    nm = length(target_mths)
    scl_rate_yr, scl_rate_mth = zeros(Float64, nr, ns), zeros(Float64, nr, ns, nm)
    cyi, tyi = findfirst(x->x==cy, yrs), findfirst(x->x==ty, yrs)
    cmi, tmi = [findfirst(x->x==cm, mths) for cm in current_mths], [findfirst(x->x==tm, mths) for tm in target_mths]
    for i = 1:nr, j = 1:ns
        scl_rate_yr[i, j] = cpi_vals_yr[i, j, tyi] / cpi_vals_yr[i, j, cyi]
        for k = 1:nm; scl_rate_mth[i, j, k] = cpi_vals_mth[i, j, tmi[k]] / cpi_vals_mth[i, j, cmi[k]] end
    end
    for i = 1:nr, j = 1:ns
        if isnan(scl_rate_yr[i, j]) || isinf(scl_rate_yr[i, j])
            println("scale rate is not available: ",cy, ", ", ty, ", ", cpi_reg[i], ", ", cpi_sec[j])
        end
        for k = 1:nm; if isnan(scl_rate_mth[i, j, k]) || isinf(scl_rate_mth[i, j, k])
            println("scale rate is not available: ", current_mths[k], ", ", target_mths[k], ", ", cpi_reg[i], ", ", cpi_sec[j])
        end end
    end

    # determine regional scaling rates
    pr_code, ds_code = prov_list[y][n], dist_list[y][n]
    npr, nds = length(pr_code), length(ds_code)
    scl_yr_nt, scl_mth_nt = zeros(Float64, ns), zeros(Float64, ns, nm)
    scl_yr_pr, scl_mth_pr = zeros(Float64, npr, ns), zeros(Float64, npr, ns, nm)
    scl_yr_ds, scl_mth_ds = zeros(Float64, nds, ns), zeros(Float64, nds, ns, nm)
    cpi_dist = filter(x-> haskey(ces_cpi_region, x), ds_code)
    cpi_prov = filter(x-> haskey(ces_cpi_region, x), pr_code)

    for i = 1:npr
        prc = pr_code[i]
        if prc in cpi_prov
            pidx = findfirst(x -> x == ces_cpi_region[prc], cpi_reg)
            for j = 1:ns; scl_yr_pr[i, j] = scl_rate_yr[pidx, j] end
            for j = 1:ns, k = 1:nm; scl_mth_pr[i, j, k] = scl_rate_mth[pidx, j, k] end
        else
            pr_dist = filter(x -> ds_pr[x] == prc, cpi_dist)
            pr_pop = sum([pop[ds] for ds in pr_dist])
            didx = [findfirst(x -> x == ds, cpi_reg) for ds in pr_dist]
            for j = 1:ns, dsi in didx
                ds = cpi_reg[dsi]
                scl_yr_pr[i, j] += scl_rate_yr[dsi, j] * pop[ds]
                for k = 1:nm; scl_mth_pr[i, j, k] += scl_rate_mth[dsi, j, k] * pop[ds] end
            end
            for j = 1:ns
                if scl_yr_pr[i, j] > 0; scl_yr_pr[i, j] /= pr_pop end
                for k = 1:nm; if scl_mth_pr[i, j, k] > 0; scl_mth_pr[i, j, k] /= pr_pop end end
            end
        end
    end
    ntidx = findfirst(x -> x == nation, cpi_reg)
    if ntidx != nothing
        for i = 1:ns; scl_yr_nt[i] = scl_rate_yr[ntidx, i] end
        for i = 1:ns, j = 1:nm; scl_mth_nt[i, j] = scl_rate_mth[ntidx, i, j] end
    else
        pr_pop = [haskey(pop,prc) && pop[prc]>0 ? pop[prc] : sum([pop[ds] for ds in filter(x->ds_pr[x]==prc, cpi_dist)]) for prc in pr_code]
        for j = 1:ns
            ntpop = 0
            for i = 1:npr
                if scl_yr_pr[i, j] > 0
                    ntpop += pr_pop[i]
                    scl_yr_nt[j] += scl_yr_pr[i, j] * pr_pop[i]
                end
            end
            scl_yr_nt[j] /= ntpop
            for k = 1:nm
                ntpop = 0
                for i = 1:npr
                    if scl_mth_pr[i, j, k] > 0
                        ntpop += pr_pop[i]
                        scl_mth_nt[j, k] += scl_mth_pr[i, j, k] * pr_pop[i]
                    end
                end
                scl_mth_nt[j, k] /= ntpop
            end
        end
    end
    for i = 1:npr, j = 1:ns
        if scl_yr_pr[i, j] == 0; scl_yr_pr[i, j] = scl_yr_nt[j] end
        for k = 1:nm; if scl_mth_pr[i, j, k] == 0; scl_mth_pr[i, j, k] = scl_mth_nt[j, k] end end
    end

    for i = 1:nds
        dsc = ds_code[i]
        if dsc in cpi_dist
            didx = findfirst(x -> x == ces_cpi_region[dsc], cpi_reg)
            for j = 1:ns; scl_yr_ds[i, j] = scl_rate_yr[didx, j] end
            for j = 1:ns, k = 1:nm; scl_mth_ds[i, j, k] = scl_rate_mth[didx, j, k] end
        else
            pidx = findfirst(x -> x == ds_pr[dsc], pr_code)
            if pidx == nothing; println("Province code matching error: ", dsc, ", ", ds_pr[dsc]) end
            for j = 1:ns; scl_yr_ds[i, j] = scl_yr_pr[pidx, j] end
            for j = 1:ns, k = 1:nm; scl_mth_ds[i, j, k] = scl_mth_pr[pidx, j, k] end
        end
    end

    # scaling expenditure matrix
    if revMat; nsl = length(sl); em = expMatrix[year][nation] end
    sec_idx = Dict(c => findfirst(x->x==ces_cpi_link[c], cpi_sec) for c in filter(x->haskey(ces_cpi_link, x), sl))
    if top_sector in cpi_sec
        tsi = findfirst(x -> x == top_sector, cpi_sec)
        for s in filter(x -> !haskey(ces_cpi_link, x), sl); sec_idx[s] = tsi end
    end

    if period == "year"; if region == "province"; scl = scl_yr_pr; elseif region =="district"; scl = scl_yr_ds end
    elseif period == "month"; if region == "province"; scl = scl_mth_pr; elseif region =="district"; scl = scl_mth_ds end
    end
    for i = 1:length(hl)
        hh = hhs[hl[i]]
        if region == "province"; reg_idx = findfirst(x->x== hh.province, pr_code)
        elseif region =="district"; reg_idx = findfirst(x->x== hh.district, ds_code)
        end
        if period == "month" && length(hh.date)>=6; midx = findfirst(x->x==hh.date[1:6], current_mths) else midx = 1 end
        if revHH
            hh.aggexp = 0.0
            for he in hh.expends
                he.value *= scl[reg_idx, sec_idx[he.code], midx]
                hh.aggexp += he.value
            end
        end
        if revMat; for j=1:nsl; em[i, j] *= scl[reg_idx, sec_idx[sl[j]], midx] end end
    end
end

function scalingExpByGDP(year, nation, statFile; wgh_mode="district", ur_dist=false)
    # scaling the expenditures to mitigate the differences between national expenditure accounts (GDP) and micro-data

    global hh_list, households, pop_wgh, pop_ur_wgh
    hl = hh_list[year][nation]
    hhs = households[year][nation]
    pw = pop_wgh[year][nation]
    pwur = pop_ur_wgh[year][nation]

    # read national expenditure accounts
    nat_exp = 0
    val_sep = getValueSeparator(statFile)
    f = open(statFile)
    readline(f)
    for l in eachline(f)
        s = strip.(split(l, val_sep))
        if parse(Int, s[1]) == year && s[2] == nation; nat_exp = parse(Float64, s[3]) end
    end
    close(f)

    # calculate nation-level total expenditure
    md_exp = 0
    if !ur_dist
        for h in hl
            hh = hhs[h]
            if wgh_mode=="district"; md_exp += hh.aggexp * pw[hh.district]
            elseif wgh_mode=="province"; md_exp += hh.aggexp * pw[hh.province]
            end
        end
    elseif ur_dist
        for h in hl
            hh = hhs[h]
            if hh.regtype == "urban"; uridx = 1 else hh.regtype == "rural"; uridx = 2 end
            if wgh_mode=="district"; md_exp += hh.aggexp * pwur[hh.district][uridx]
            elseif wgh_mode=="province"; md_exp += hh.aggexp * pwur[hh.province][uridx]
            end
        end
    end

    # scaling expenditure
    scl_rate = nat_exp / md_exp
    for h in hl
        hhs[h].aggexp *= scl_rate
        for e in hhs[h].expends; e.value *= scl_rate end
    end
end

function printRegionData(year, nation, outputFile; region = "district", ur = false)

    global regions, prov_list, dist_list, dist_prov, pops, pops_ur, pop_wgh, pop_ur_wgh
    rg = regions[year][nation]
    ps, pw = pops[year][nation], pop_wgh[year][nation]
    if ur; psur, pwur = pops_ur[year][nation], pop_ur_wgh[year][nation] end

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    if region == "province"
        rl = prov_list[year][nation]
        print(f, "Code\tProvince_ID\tProvince_name\tPopulation\tWeight")
        # print(f, "Code\tCode_State/Province\tState/Province\tPopulation\tWeight")
    elseif region == "district"
        rl, dp = dist_list[year][nation], dist_prov[year][nation]
        print(f, "Code\tProvince_ID\tProvince_name\tCity_ID\tCity_name\tPopulation\tWeight")
        # print(f, "Code\tCode_State/Province\tState/Province\tCode_District/City\tDistrict/City\tPopulation\tWeight")
    end
    if ur; print(f, "\tPop_urban\tPop_rural\tWgh_urban\tWgh_rural") end
    println(f)
    count = 0
    for r in rl
        if region == "province"; print(f, r, "\t", r, "\t", rg[r], "\t", ps[r], "\t", pw[r])
        elseif region == "district"; print(f, r, "\t", dp[r], "\t", rg[dp[r]], "\t", r, "\t" , rg[r], "\t", ps[r], "\t", pw[r])
        end
        if ur; print(f, "\t", psur[r][1], "\t", psur[r][2], "\t", pwur[r][1], "\t", pwur[r][2]) end
        println(f)
        count += 1
    end
    close(f)
    println("$count regions' data is printed.")
end

function printCommoditySectors(year, nation, outputFile)

    global sc_list, sectors
    sl = sc_list[year][nation]
    sc = sectors[year][nation]

    count = 0
    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    println(f, "Code\tSector\tMain_category\tUnit\tUnit_quantity")
    for s in sl
        println(f, s, "\t", sc[s].sector, "\t", sc[s].category, "\t", sc[s].unit_cur, "\t", sc[s].unit_qnt)
        count += 1
    end
    close(f)
    println("$count commodities' data is printed.")
end

function printHouseholdData(year, nation, outputFile; hh_wgh=false, tot_inc = false, ur_dist = false, religion = false, surv_date = false, surv_type = false)

    global households, hh_list, regions, pop_wgh, pop_ur_wgh, exp_curr, exp_period

    hhs = households[year][nation]
    # rg = regions[year][nation]
    # pw = pop_wgh[year][nation]
    # if ur_dist; pwur = pop_ur_wgh[year][nation] end

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    count = 0

    print(f, "HHID\tProvince_ID\tCity_ID\tHH_size\tTotal_exp\tTot_exp_unit")
    # print(f, "HHID\tCode_province/state\tCode_district/city\tHH_size\tTotal_exp\tTot_exp_unit")
    if hh_wgh; print(f, "\tPop_wgh_percap") end
    if tot_inc; print(f, "\tTotal_inc") end
    if ur_dist; print(f, "\tRegion_type") end
    if religion; print(f, "\tReligion") end
    if surv_date; print(f, "\tStart_date\tEnd_date") end
    if surv_type; print(f, "\tSurvey_type") end
    println(f)

    for h in hh_list[year][nation]
        hh = hhs[h]
        print(f, hh.hhid, "\t", hh.province, "\t", hh.district, "\t", hh.size, "\t", hh.totexp, "\t", hh.unit)
        if hh_wgh; print(f, "\t", hh.popwgh) end
        if tot_inc; print(f, "\t", hh.totinc) end
        if ur_dist; print(f, "\t", hh.regtype) end
        if religion; print(f, "\t", hh.rel) end
        if surv_date; print(f, "\t", hh.date) end
        if surv_type; print(f, "\t", hh.group) end
        println(f)
        count += 1
    end
    close(f)
    println("$count households' data is printed.")
end

function printExpenditureData(year, nation, outputFile; quantity = false)

    global households, hh_list
    hhs = households[year][nation]
    hl = hh_list[year][nation]

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    count = 0; err_cnt = 0

    print(f, "HHID\tCommodity_code\tExpended_value\tCurrency_unit")
    if quantity; print(f, "\tConsumed_quantity\tQuantity_unit") end
    println(f, "\tConsumption_period")
    for h in hl
        hh = hhs[h]
        for e in hh.expends
            print(f, hh.hhid, "\t", e.code, "\t", e.value, "\t", e.valUnit)
            if quantity; print(f, "\t", e.quantity, "\t", e.qntUnit) end
            println(f, "\t", e.period)
            count += 1
            if quantity && e.value == e.quantity == 0; err_cnt += 1 end
        end
    end
    close(f)
    if quantity; print("$err_cnt errors, ") end
    println("$count items' data is printed.")
end

function printExpenditureMatrix(year, nation, outputFile = ""; quantity = false, rowErr = [], colErr = [])

    global households, hh_list, sc_list, expMatrix, qntMatrix

    hhs = households[year][nation]
    mat = Array{Array{Float64, 2}, 1}()
    push!(mat, expMatrix[year][nation])
    if quantity; push!(mat, qntMatrix[year][nation]) end
    row, col = hh_list[year][nation], sc_list[year][nation]
    nr, nc = length(row), length(col)
    nre, nce = length(rowErr), length(colErr)

    f_tag = ["", "_qnt"]
    f_sep = getValueSeparator(outputFile)
    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    for i = 1:length(mat)
        m = mat[i]
        f = open(replace(outputFile, ".txt" => f_tag[i] * ".txt"), "w")
        print(f, "HHID/Sector")
        for c in col; print(f, f_sep, c) end
        # if nre > 0; print(f, f_sep, "Row_error") end
        println(f)
        for ri = 1:nr
            print(f, row[ri])
            for ci = 1:nc; print(f, f_sep, m[ri,ci]) end
            # if nre > 0; print(f, f_sep, rowErr[ri]) end
            println(f)
        end
        # if nce > 0
        #     print(f, "Column error")
        #     for ci = 1:nc; print(f, f_sep, colErr[ci]) end
        #     print(f, f_sep, sum(colErr))
        # end
        println(f)
        close(f)
    end
    println("$nr by $nc expenditure matrix is printed.")
end

function initVars()
    global households = Dict{Int, Dict{String, Dict{String, household}}}()
    global sectors = Dict{Int, Dict{String, Dict{String, commodity}}}()
    global hh_list = Dict{Int, Dict{String, Array{String, 1}}}()
    global sc_list = Dict{Int, Dict{String, Array{String, 1}}}()
    global regions = Dict{Int, Dict{String, Dict{String, String}}}()
    global prov_list = Dict{Int, Dict{String, Array{String, 1}}}()
    global dist_list = Dict{Int, Dict{String, Array{String, 1}}}()
    global dist_prov = Dict{Int, Dict{String, Dict{String, String}}}()
    global pops = Dict{Int, Dict{String, Dict{String, Float64}}}()
    global pops_ur = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}()
    global pop_wgh = Dict{Int, Dict{String, Dict{String, Float64}}}()
    global pop_ur_wgh = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}()
    global exchange_rate = Dict{String, Dict{String, Float64}}()
    global expMatrix = Dict{Int, Dict{String, Array{Float64, 2}}}()
    global qntMatrix = Dict{Int, Dict{String, Array{Float64, 2}}}()
end

function exportCommodityUnit(year, nation)

    global sc_list, sectors
    sc = sectors[year][nation]
    qnt_unit = Array{String, 1}()

    for s in sc_list[year][nation]; push!(qnt_unit, sc[s].unit_qnt) end

    return qnt_unit
end

function readExtractedRegionData(year, nation, inputFile; key_district = false, merged_key = false, legacy_mode = true, ignore = false)
    # legacy_mode: [true] apply the previous item label for the previously made data (should be removed after all revisions)

    global regions, prov_list, dist_list, dist_prov, pops, pops_ur, pop_wgh, pop_ur_wgh
    essential = ["Code", "Province_ID", "Province_name", "City_ID", "City_name", "Population"]
    essential_lag = ["Code", "Code_State/Province", "State/Province", "Code_District/City", "District/City", "Population"]
    optional = ["Weight"]
    ur_title = ["Pop_urban", "Pop_rural", "Wgh_urban", "Wgh_rural"]

    essential, essential_lag, optional, ur_title = lowercase.(essential), lowercase.(essential_lag), lowercase.(optional), lowercase.(ur_title)

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = lowercase.(string.(strip.(split(readline(f), f_sep))))
    if issubset(essential, title); i = [findfirst(x->x==item, title) for item in essential]
    elseif legacy_mode && issubset(essential_lag, title); i = [findfirst(x->x==item, title) for item in essential_lag]
    else println(inputFile, " household file does not contain all essential data.")
    end
    io = [findfirst(x->x==item, title) for item in optional]
    iu = [findfirst(x->x==item, title) for item in ur_title]
    op_chk = all(io.!=nothing)
    ur_p_chk, ur_w_chk = all(iu[[1,2]].!=nothing), all(iu[[3,4]].!=nothing)
    if !haskey(regions, year)
        regions[year] = Dict{String, Dict{String, String}}()
        pops[year] = Dict{String, Dict{String, Float64}}()
        pop_wgh[year] = Dict{String, Dict{String, Float64}}()
        prov_list[year] = Dict{String, Array{String, 1}}()
        dist_list[year] = Dict{String, Array{String, 1}}()
        dist_prov[year] = Dict{String, Dict{String, String}}()
        if ur_p_chk; pops_ur[year] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
        if ur_w_chk; pop_ur_wgh[year] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
    end
    regions[year][nation] = Dict{String, String}()
    pops[year][nation] = Dict{String, Float64}()
    pop_wgh[year][nation] = Dict{String, Float64}()
    prov_list[year][nation] = Array{String, 1}()
    dist_list[year][nation] = Array{String, 1}()
    dist_prov[year][nation] = Dict{String, String}()
    if ur_p_chk; pops_ur[year][nation] = Dict{String, Tuple{Float64, Float64}}() end
    if ur_w_chk; pop_ur_wgh[year][nation] = Dict{String, Tuple{Float64, Float64}}() end

    count = 0
    r_cds = Array{String, 1}()
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        if ignore || all(length.(s) .> 0)
            pr_code = s[i[2]]
            ds_code = (merged_key ? pr_code * "_" * s[i[4]] : s[i[4]])
            r_cd = key_district ? ds_code : s[i[1]]
            push!(prov_list[year][nation], pr_code)
            push!(dist_list[year][nation], ds_code)
            dist_prov[year][nation][ds_code] = pr_code
            rg = regions[year][nation]
            if !haskey(rg, pr_code) rg[pr_code] = s[i[3]] end
            if !haskey(rg, ds_code) rg[ds_code] = s[i[5]] end
            pops[year][nation][r_cd] = parse(Float64, s[i[6]])
            if op_chk && s[io[1]] != ""; pop_wgh[year][nation][r_cd] = parse(Float64, s[io[1]]) end
            if ur_p_chk && !(s[iu[1]]==s[iu[2]]==""); pops_ur[year][nation][r_cd] = (parse(Float64, s[iu[1]]), parse(Float64, s[iu[2]])) end
            if ur_w_chk && !(s[iu[3]]==s[iu[4]]==""); pop_ur_wgh[year][nation][r_cd] = (parse(Float64, s[iu[3]]), parse(Float64, s[iu[4]])) end
            push!(r_cds, r_cd)
            count += 1
        end
    end
    close(f)
    print("$count regions")
    if length(r_cds) != length(unique(r_cds))
        println("\nDuplicated region key codes: ", length(r_cds) - length(unique(r_cds)), " duplicated")
    end
end

function readExtractedHouseholdData(year, nation, inputFile; period = "year", merged_key = false, skip_empty = true, legacy_mode = true)
    # skip_empty: [true] exclude household that does not have district code
    # legacy_mode: [true] apply the previous item label for the previously made data (should be removed after all revisions)
    # split_mode: [true] split household data for multiple-group hhs

    global households, hh_list, regions, hh_curr, hh_period, pr_scl
    essential = ["HHID", "Province_ID", "City_ID", "HH_size", "Total_exp", "Currency"]
    essential_lag = ["HHID", "Code_province/state", "Code_district/city", "HH_size", "Total_exp", "Tot_exp_unit"]
    optional = ["Pop_wgh_percap", "Total_inc", "Region_type", "Religion", "Start_date", "End_date", "Survey"]

    essential, essential_lag, optional = lowercase.(essential), lowercase.(essential_lag), lowercase.(optional)

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = lowercase.(string.(strip.(split(readline(f), f_sep))))
    if issubset(essential, title); i = [findfirst(x->x==t, title) for t in [essential;optional]]
    elseif legacy_mode && issubset(essential_lag, title); i = [findfirst(x->x==t, title) for t in [essential_lag;optional]]
    else println(inputFile, " household file does not contain all essential data.")
    end
    if !haskey(households, year)
        households[year] = Dict{String, Dict{String, household}}()
        hh_list[year] = Dict{String, Array{String, 1}}()
    end
    if !haskey(hh_curr, year); hh_curr[year] = Dict{String, Array{String, 1}}() end
    if !haskey(hh_period, year); hh_period[year] = Dict{String, Array{String, 1}}() end

    hhs = households[year][nation] = Dict{String, household}()
    hl = hh_list[year][nation] = Array{String, 1}()
    hhc = hh_curr[year][nation] = Array{String, 1}()
    hhp = hh_period[year][nation] = Array{String, 1}()

    chk_pw, chk_ic, chk_rt, chk_sd, chk_sv = [isa(i[oi], Int) for oi in [7, 8, 9, 11, 13]]
    # chk_pw, chk_ic, chk_rt, chk_sd, chk_sv = isa(i[7], Int), isa(i[8], Int), isa(i[9], Int), isa(i[11], Int), isa(i[13], Int)
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        if skip_empty && length(s[i[3]]) == 0; continue end
        hhid = s[i[1]]
        push!(hl, hhid)
        hhs[hhid] = household(hhid)
        hh_vals = ["","","","",0,0,"","","",0,0,0,0,"",0,[],[],""]
        ds_code = (merged_key ? s[i[2]] * "_" * s[i[3]] : s[i[3]])
        hh_vals[[2,3,5,10,14]] = [s[i[2]], ds_code, parse(Int, s[i[4]]), parse(Float64, s[i[5]]), s[i[6]]]
        hh_vals[4] = chk_rt ? s[i[9]] : ""
        hh_vals[1] = chk_sd ? s[i[11]] : string(year)
        hh_vals[15] = chk_pw && s[i[7]] != "" ? parse(Float64, s[i[7]]) : 0.0
        hh_vals[11] = chk_ic && s[i[8]] != "" ? parse(Float64, s[i[8]]) : 0.0
        hh_vals[18] = chk_sv ? s[i[13]] : ""

        currency = s[i[6]]
        crr_scl = tryparse(Float64, filter(isdigit, currency))
        if crr_scl != nothing
            for vi = 10:13; hh_vals[vi] *= crr_scl end
            currency = filter(!isdigit, currency)
        end
        pru = string(rsplit(currency, '/', limit=2)[end])
        if haskey(pr_scl, pru)
            crr_scl = pr_scl[period] / pr_scl[pru]
            for vi = 10:13; hh_vals[vi] *= crr_scl end
            currency = string(rsplit(currency, '/', limit=2)[1])
        end
        reviseHouseholdData(year, nation, hhid, hh_vals)

        if !(currency in hhc); push!(hhc, currency) end
    end
    if !(period in hhp); push!(hhp, period) end
    close(f)
    print(length(hl), " households")
end

function readExtractedMemberData(year, nation, inputFile)

    global households, hh_list
    essential = ["HHID", "Age", "Gender"]

    essential = lowercase.(essential)

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = lowercase.(string.(strip.(split(readline(f), f_sep))))
    if issubset(essential, title); i = [findfirst(x->x==et, title) for et in essential]
    else println(inputFile, " member file does not contain all essential data.")
    end
    hhs = households[year][nation]

    count = 0
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        hhid = s[i[1]]
        push!(hhs[hhid].members, member(hhid, parse(Int16, s[i[2]]), parse(Int8, s[i[3]])))
        count += 1
    end
    close(f)
    print("$count members")
end

function readExtractedExpenditureData(year, nation, inputFile; quantity = false, deflet_unit = true)
    # deflet_unit: remove the scale in the expenditure currency unit (ex, 1000USD -> USD)

    global households, hh_list, exp_curr, exp_period, pr_unts
    essential = ["HHID", "Commodity_code", "Expended_value", "Currency_unit", "Consumption_period"]
    if quantity; essential = [essential[1:4]; ["Consumed_quantity", "Quantity_unit"] ; essential[5:end]] end
    if !haskey(exp_curr, year) exp_curr[year] = Dict{String, Array{String, 1}}() end
    if !haskey(exp_period, year) exp_period[year] = Dict{String, Array{String, 1}}() end
    if !haskey(exp_curr[year], nation) exp_curr[year][nation] = Array{String, 1}() end
    if !haskey(exp_period[year], nation) exp_period[year][nation] = Array{String, 1}() end

    if deflet_unit; exp_curr[year][nation] = unique([filter(!isdigit, ec) for ec in exp_curr[year][nation]]) end

    hhs = households[year][nation]
    e_curr, e_per = exp_curr[year][nation], exp_period[year][nation]

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = string.(strip.(split(readline(f), f_sep)))
    if issubset(essential, title); i = [findfirst(x->x==et, title) for et in essential]
    else println(inputFile, " expenditure file does not contain all essential data.")
    end

    count = 0
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        hhid = s[i[1]]
        e_val = length(s[i[3]]) > 0 ? parse(Float64, s[i[3]]) : 0.0
        c_unt = string(split(s[i[4]], '/', limit = 2)[1])
        c_scl = tryparse(Float64, filter(isdigit, c_unt))
        if c_scl != nothing
            e_val *= c_scl
            c_unt = filter(!isdigit, c_unt)
        end
        exp_vals = [s[i[2]], e_val, 0, c_unt, "", parse(Int16, s[i[end]])]
        if quantity
            q_val = length(s[i[5]]) > 0 ? parse(Float64, s[i[5]]) : 0.0
            q_unt = string(split(s[i[6]], '/', limit = 2)[1])
            q_scl = tryparse(Float64, filter(isdigit, c_unt))
            if q_scl != nothing
                q_val *= q_scl
                q_unt = filter(!isdigit, q_unt)
            end
            exp_vals[[3, 5]] = [q_val, q_unt]
        end

        appendExpenditureData(year, nation, hhid, exp_vals)
        if !(exp_vals[4] in e_curr); push!(e_curr, exp_vals[4]) end
        if !(pr_unts[exp_vals[6]] in e_per); push!(e_per, pr_unts[exp_vals[6]]) end
        count += 1
    end
    close(f)
    print("$count expenditures")
end

function readExtractedSectorData(year, nation, itemfile)

    global sc_list, sectors, exp_curr, gr_list
    if !haskey(sc_list, year); sc_list[year] = Dict{String, Array{String, 1}}() end
    if !haskey(sectors, year); sectors[year] = Dict{String, Dict{String, commodity}}() end
    if !haskey(exp_curr, year); exp_curr[year] = Dict{String, Array{String, 1}}() end
    if !haskey(gr_list, year); gr_list[year] = Dict{String, Array{String, 1}}() end
    sl = sc_list[year][nation] = Array{String, 1}()
    sc = sectors[year][nation] = Dict{String, commodity}()
    ec = exp_curr[year][nation] = Array{String, 1}()
    gl = gr_list[year][nation] = Array{String, 1}()

    essential = ["Code", "Sector", "Main_category", "Unit"]
    optional = ["Unit_quantity", "Survey"]

    essential, optional = lowercase.(essential), lowercase.(optional)

    f_sep = getValueSeparator(itemfile)
    count = 0
    f = open(itemfile)
    title = lowercase.(string.(strip.(split(readline(f), f_sep))))
    if issubset(essential, title); i = [findfirst(x->x==t, title) for t in [essential;optional]]
    else println(itemfile, " commodity sector file does not contain all essential data.")
    end
    chk_uq, chk_gr = isa(i[5], Int), isa(i[6], Int)
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        push!(sl, s[i[1]])
        u_qnt = chk_uq ? s[i[5]] : ""
        grp = chk_gr ? s[i[6]] : ""
        sc[s[i[1]]] = commodity(s[i[1]], s[i[2]], s[i[3]], "", s[i[4]], u_qnt, "", grp)
        if length(s[i[4]]) > 0 && !(s[i[4]] in ec); push!(ec, s[i[4]]) end
        if chk_gr && length(grp) > 0 && !(grp in gl); push!(gl, grp) end
        count += 1
    end
    close(f)
    print("$count sectors")
end

function readExtractedExpenditureMatrix(year, nation, inputFile; quantity = false)

    global hh_list, sc_list, expMatrix, qntMatrix, sectors, exp_curr

    if !haskey(expMatrix, year); expMatrix[year] = Dict{String, Array{Float64, 2}}() end
    if !haskey(exp_curr, year); exp_curr[year] = Dict{String, Array{String, 1}}() end
    if !haskey(exp_curr[year], nation); exp_curr[year][nation] = Array{String, 1}() end
    sl, hl = sc_list[year][nation], hh_list[year][nation]
    ns, nh = length(sl), length(hl)
    em = expMatrix[year][nation] = zeros(Float64, nh, ns)
    d_set = [[inputFile, em, "expenditure"]]

    if quantity
        if !haskey(qntMatrix, year); qntMatrix[year] = Dict{String, Array{Float64, 2}}() end
        qm = qntMatrix[year][nation] = zeros(Float64, nh, ns)
        nat_curr = rsplit(rsplit(inputFile, '.', limit = 2)[1], '_', limit = 2)[end]
        append!(d_set, [replace(nat_curr => "qnt", inputFile), qm, "quantity"])
    end

    for ds in d_set
        ifile, mat, tag = ds

        f_sep = getValueSeparator(ifile)
        f = open(ifile)
        codes = string.(strip.(split(readline(f), f_sep)[2:end]))
        if issubset(sl, codes); i = [findfirst(x->x==sc, codes) for sc in sl]
        else println(inputFile, " " * tag * " matrix file does not contain all essential data.")
        end

        hhs = Array{String, 1}()
        for l in eachline(f)
            s = string.(strip.(split(l, f_sep)))
            hi = findfirst(x->x==s[1], hl)
            vals = parse.(Float64, s[2:1+ns])
            if hi !=nothing
                push!(hhs, s[1])
                mat[hi,:] = vals
            end
        end
        close(f)
        if sort(hhs) != sort(hl)
            println(inputFile, " " * tag * " matrix file does not contain all household data.")
        end
    end

    exp_c = exp_curr[year][nation]
    if length(exp_c) == 0
        exp_c = [string(rsplit(rsplit(inputFile, '.', limit = 2)[1], '_', limit = 2)[2])]
    elseif length(exp_c) == 1
        crr_scl = tryparse(Float64, filter(isdigit, exp_c[1]))
        if crr_scl != nothing; em .*= crr_scl end
    elseif length(exp_c) > 1
        sc = sectors[year][nation]
        crr_scl = [tryparse(Float64, filter(isdigit, sc[s].unit_cur)) for s in sl]
        if any(crr_scl .!= nothing)
            crr_scl = [scl == nothing ? 1.0 : scl for scl in crr_scl]
            em .*= crr_scl'
        end
    end
    exp_curr[year][nation] = [filter(!isdigit, c) for c in exp_c]
end

function reshapeCommoditySectors(year, nation; except = ["None", "Taxes"], hhs_reshape = false)

    global sectors, sc_list, expMatrix, qntMatrix, households, hh_list
    y, n = year, nation
    sec, sl, hl, hhs = sectors[y][n], sc_list[y][n], hh_list[y][n], households[y][n]
    ns = length(sl)
    except = lowercase.(except)

    rs_idx = filter(x -> !(lowercase(sec[sl[x]].category) in except), 1:ns)

    sc_list[y][n] = sl[rs_idx]
    filter!(x -> x.first in sl, sec)
    if haskey(expMatrix, y) && haskey(expMatrix[y], n); expMatrix[y][n] = expMatrix[y][n][:,rs_idx] end
    if haskey(qntMatrix, y) && haskey(qntMatrix[y], n); qntMatrix[y][n] = qntMatrix[y][n][:,rs_idx] end

    if hhs_reshape
        for h in filter(x -> length(hhs[x].expends) > 0, hl)
            filter!(x -> x.code in sl, hhs[h].expends)
            hhs[h].aggexp = sum([he.value for he in hhs[h].expends])
        end
    end
end

end
