module MicroDataReader

# Developed date: 17. Mar. 2021
# Last modified date: 2. Jun. 2021
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

    household(hi,da="",pr="",di="",rt="",sz=0,ag=0,rl="",oc="",ed="",te=0,ti=0,tep=0,tip=0,ut="",pw=0,ae=0,mm=[],ex=[]) = new(hi,da,pr,di,rt,sz,ag,rl,oc,ed,te,ti,tep,tip,ut,pw,ae,mm,ex)
end

mutable struct commodity
    code::String        # commodity code in the survey data
    sector::String      # commodity's sector or label
    category::String    # main category, ex) Food, Electricity, Gas, Other energy, Public transport, Private transport, Medical care, Education, Consumable goods, Durable goods, Other services
    subCategory::String # sub-category, ex) Food related: Grain, Vegetable, Fruit, Dairy, Beef, Pork, Poultry, Other meat, Fish, Alcohol, Other beverage, Confectionery, Restaurant, Other food, etc
                        #                   Energy-related: Electricity, Gas, Wood, Dung cake, Kerosene, Coal, Petrol, Diesel, Biogas, Other fuel, etc.
                        #                   Transport-related: Road (private), Road (public), Rail, Air, Water, Other, etc.
    unit::String        # commodity's physical unit, ex) liter, kg, m^3, etc.
    coicop::String      # commodity corresponding COICOP code

    commodity(cod, sec="", cat="", subcat="", unt = "", coi="") = new(cod, sec, cat, subcat, unt, coi)
end

global households = Dict{Int, Dict{String, Dict{String, household}}}()  # household dict: {year, {nation A3, {hhid, household}}}
global sectors = Dict{Int, Dict{String, Dict{String, commodity}}}()     # expenditure sector: {year, {nation A3, {code, commodity}}}
global hh_list = Dict{Int, Dict{String, Array{String, 1}}}()            # hhid list: {year, {nation A3, {hhid}}}
global sc_list = Dict{Int, Dict{String, Array{String, 1}}}()            # commodity code list: {year, {nation A3, {code}}}

global regions = Dict{Int, Dict{String, Dict{String, String}}}()        # region code-name: {year, {nation A3, {code, region}}}
global prov_list = Dict{Int, Dict{String, Array{String, 1}}}()          # province code list: {year, {nation A3, {code}}}
global dist_list = Dict{Int, Dict{String, Array{String, 1}}}()          # district code list: {year, {nation A3, {code}}}
global dist_prov = Dict{Int, Dict{String, Dict{String, String}}}()      # district's province: {year, {nation A3, {district code, province code}}}

global pops = Dict{Int, Dict{String, Dict{String, Float64}}}()          # population: {year, {nation, {region_code, population}}}
global pops_ur = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}()   # urban/rural population: {year, {nation, {region_code, (urban, rural)}}
global pop_wgh = Dict{Int, Dict{String, Dict{String, Float64}}}()       # population weight: {year, {nation, {region_code, weight}}}
global pop_ur_wgh = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}()    # urban/rural population weight: {year, {nation, {region_code, (urban, rural)}}

global hh_curr = Dict{Int, Dict{String, Array{String, 1}}}()            # currency unit for household values (income or expenditure): {year, {nation, {currency}}}
global hh_period = Dict{Int, Dict{String, Array{String, 1}}}()          # period for household values (income or expenditure): {year, {nation, {period}}}
global exp_curr = Dict{Int, Dict{String, Array{String, 1}}}()           # currency unit for expenditure values: {year, {nation, {currency}}}
global exp_period = Dict{Int, Dict{String, Array{String, 1}}}()         # period for expenditure values: {year, {nation, {period}}}
global exchange_rate = Dict{String, Dict{String, Float64}}()            # exchange rate: {currency("original/target", ex."EUR/USD", {year(string), converting_rate}}

global expMatrix = Dict{Int, Dict{String, Array{Float64, 2}}}()         # expenditure matrix: {year, {nation, {hhid, commodity}}}
global qntMatrix = Dict{Int, Dict{String, Array{Float64, 2}}}()         # qunatity matrix: {year, {nation, {hhid, commodity}}}

global pr_unts = Dict(1=>"day", 7=>"week", 30=>"month", 365=>"year")    # period units

function appendCommoditySectorData(year, nation, cmm_data)
    global sectors, sc_list

    push!(sc_list[year][nation], cmm_data[1])
    sectors[year][nation][cmm_data[1]] = commodity(cmm_data[1], cmm_data[2], cmm_data[3], "", cmm_data[4], cmm_data[5])
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

    if length(id)>0 && length(hh_value)==17
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

function readRegion(year, nation, regionFile)

    global regions, prov_list, dist_list, dist_prov

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
end

function readMicroData(year, nation, microdataPath, hhIdxFile, mmIdxFile, cmmIdxFile, expIdxFile; hhid_sec = "hhid", periodFiltering = false, ignoreException = true)

    for idxfile in filter(x->length(x)>0, [hhIdxFile, mmIdxFile, cmmIdxFile, expIdxFile])
        # read microdata index file
        yrchk, natchk, idxs, ipdx_o = readIndexFile(year, nation, idxfile)

        # read microdata contents
        if idxfile == hhIdxFile; readHouseholdData(year, nation, [idxs;ipdx_o], microdataPath, hhid_sec = hhid_sec)  # idxs: {sector, position, type, file, tag}
        elseif idxfile == mmIdxFile; readMemberData(year, nation, idxs, microdataPath, hhid_sec = hhid_sec)          # idxs: {sector, position, type, file, tag}
        elseif idxfile == cmmIdxFile; readCommoditySectors(year, nation, idxs)   # idxs: {code, coicop, sector, entity, category, sub_category}
        elseif idxfile == expIdxFile; readExpenditureData(year, nation, idxs, microdataPath, periodFiltering = periodFiltering, ignoreException = ignoreException)    # idxs: {}
        end
    end
end

function readHouseholdData(year, nation, indices, microdataPath; hhid_sec = "hhid")

    global households, hh_list, hh_curr, hh_period

    sectors = ["survey_date", "province/state", "district/city", "region_type", "hh_size", "head_age", "head_religion", "head_occupation", "head_education", "expenditure", "income", "exp_percap", "inc_percap", "currency_unit", "pop_weight", "agg_exp"]
    nsec = length(sectors)
    int_sec = ["hh_size", "head_age"]
    flo_sec = ["expenditure", "income", "exp_percap", "inc_percap", "pop_weight"]

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
            currency, period = string.(strip.(split(idx[3], '/')))
            if !(currency in hh_curr[year][nation]); push!(hh_curr[year][nation], currency) end
            if !(period in hh_period[year][nation]); push!(hh_period[year][nation], period) end
        end
    end
    mfs = sort(collect(keys(mdFiles)))
    for mf in mfs; if !haskey(mdFiles[mf], hhid_sec); println(mf, "does not contain HHID sector.") end end
    hh_unit = hh_curr[year][nation][1] * "/" * hh_period[year][nation][1]

    # read household data: all the microdata files should contain HHID values
    for mf in mfs
        mfd = mdFiles[mf]
        hhid_pos = mfd[hhid_sec][1]
        hhid_tag = mfd[hhid_sec][3]
        mdf_sep = getValueSeparator(mf)
        f = open(microdataPath * mf)
        readline(f)     # read title line
        for l in eachline(f)
            s = strip.(split(l, mdf_sep))
            hhid = hhid_tag * s[hhid_pos]
            if !haskey(households[year][nation], hhid)
                households[year][nation][hhid] = household(hhid)
                push!(hh_list[year][nation], hhid)
            end
            hh_vals = ["","","","",0,0,"","","",0,0,0,0,"",0,[],[]]
            for i = 1:nsec
                c = sectors[i]
                if haskey(mfd, c)
                    val = s[mfd[c][1]]
                    if length(val) > 0
                        if c in int_sec; hh_vals[i] = parse(Int, val)
                        elseif c in flo_sec; hh_vals[i] = parse(Float64, val)
                        else hh_vals[i] = val
                        end
                    end
                end
            end
            hh_vals[14] = hh_unit
            reviseHouseholdData(year, nation, hhid, hh_vals)
        end
        close(f)
    end
end

function readMemberData(year, nation, indices, microdataPath; hhid_sec = "hhid")

    global households

    sectors = ["age", "gender", "nationality", "head_relation", "marital_status", "education_level", "occupation", "income", "income_unit"]
    nsec = length(sectors)
    int_sec = ["age", "gender", "marital_status", "education_level", "head_relation"]
    flo_sec = ["income"]

    # analyze index data
    mdFiles = Dict{String, Dict{String, Tuple{Int, String, String}}}()  # {microdata_file, {data_sector, {position, type tag}}}
    for idx in indices      # indices: {sector, position, type, file, tag}
        mf = idx[4]
        if !haskey(mdFiles, mf); mdFiles[mf] = Dict{String, Tuple{Int, String, String}}() end
        mdFiles[mf][idx[1]] = (parse(Int, idx[2]), idx[3], idx[5])      # (position, type, tag)
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
        readline(f)     # read title line
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
                        elseif c in flo_sec; mm_vals[i] = parse(Float64, val)
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
        if !(s[1] in scl); appendCommoditySectorData(year, nation, [s; ["" for i=1:(5 - length(s))]])
        else println("Duplicated codes: ", s[1], "\t", s[3], ", ", sec[s[1]])
        end
    end
end

function readExpenditureData(year, nation, indices, microdataPath; periodFiltering = false, ignoreException = true)

    global households, hh_list, sectors, sc_list, exp_curr, exp_period, pr_unts
    if !haskey(exp_curr, year) exp_curr[year] = Dict{String, Array{String, 1}}() end
    if !haskey(exp_period, year) exp_period[year] = Dict{String, Array{String, 1}}() end
    if !haskey(exp_curr[year], nation) exp_curr[year][nation] = Array{String, 1}() end
    if !haskey(exp_period[year], nation) exp_period[year][nation] = Array{String, 1}() end

    cm = sectors[year][nation]
    sl = sc_list[year][nation]
    units = ["day", "week", "month", "year"]

    # analyze index data
    # mdFiles: {microdata_file, {category, {data_tag, hhid_position, code_position, period(days), value_position, value_unit, quantity_position, quantity_unit}}}
    mdFiles = Dict{String, Dict{String, Tuple{String, Int, Int, Int, Int, String, Int, String}}}()
    for idx in indices      # index: {category, hhid_position, code_position, period(days), value_position, value_unit, file, data_tag, quantity_position, quantity_unit}
        if all(length.(idx[[1, 2, 3, 4, 7]]).>0) && all(length.(idx[[5, 9]]).>0)
            if length(idx) == 7 || length(idx) == 8
                li = findlast(x->tryparse(Float64, x) != nothing, idx)
                if li == 5; idx = [idx; ["","0", ""]]
                elseif li == 7; idx = [idx[1:4]; ["0", ""] ;idx[5:end]]
                end
            end
            if length(idx) == 10
                idxs = [[idx[1], parse(Int, idx[2]), parse(Int, idx[3]), parse(Int, idx[4])]; idx[5:end]]
                for i in [5, 9]; if length(idxs[i]) > 0; idxs[i] = parse(Int, idxs[i]); else idxs[i] = 0 end end
                mf = idxs[7]
                if !haskey(mdFiles, mf); mdFiles[mf] = Dict{String, Tuple{String, Int, Int, Int, Int, String, Int, String}}() end
                mdFiles[mf][idxs[1]] = tuple([idxs[8]; idxs[2:6]; idxs[9:end]]...)
            else println("Expenditure index content length error: ", year, ", ", nation, "\t", idx)
            end
            if !(idxs[6] in exp_curr[year][nation]); push!(exp_curr[year][nation], idxs[6]) end
            if !(pr_unts[idxs[4]] in exp_period[year][nation]); push!(exp_period[year][nation], pr_unts[idxs[4]]) end
        else println("Expenditure index lacks essential contents: ", year, ", ", nation, "\t", idx)
        end
    end
    mfs = sort(collect(keys(mdFiles)))

    # read expenditure data: all the microdata files should contain HHID and Code values
    for mf in mfs
        mfd = mdFiles[mf]
        sec = collect(keys(mfd))
        mdf_sep = getValueSeparator(mf)
        pre_mfd = []
        for sc in sec
            dup_chk = false
            for pre_sc in pre_mfd; if mfd[sc] == mfd[pre_sc]; dup_chk = true end end
            if !dup_chk
                push!(pre_mfd, sc)
                hhid_tag, hhid_pos, code_pos = mfd[sc][1], mfd[sc][2], mfd[sc][3]
                # {data_tag, hhid_position, code_position, period(days), value_position, value_unit, quantity_position, quantity_unit}
                f = open(microdataPath * mf)
                readline(f)     # read title line
                for l in eachline(f)
                    s = strip.(split(l, mdf_sep))
                    if !ignoreException || s[code_pos] in sl
                        hhid = hhid_tag * s[hhid_pos]
                        exp_vals = ["", 0, 0, "", "", 0]
                        # for string values
                        exp_idx = [1, 4, 5]
                        mfd_idx = [3, 6, 8]
                        for i = 1:length(mfd_idx)
                            val = mfd[sc][mfd_idx[i]]
                            if isa(val, Number); val = s[val] end
                            if length(val) > 0; exp_vals[exp_idx[i]] = val end
                        end
                        unt = cm[exp_vals[1]].unit
                        if length(unt)>0; exp_vals[5] = unt end
                        # for numeric values
                        exp_idx = [2, 3]
                        mfd_idx = [5, 7]
                        for i = 1:length(mfd_idx)
                            if mfd[sc][mfd_idx[i]]>0; val = parse(Float64, s[mfd[sc][mfd_idx[i]]]) else val = 0 end
                            if val > 0; exp_vals[exp_idx[i]] = val end
                        end
                        # for period-value
                        if mfd[sc][4] > 0; exp_vals[6] = mfd[sc][4] end
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

function exchangeExpCurrency(year, exchangeYear, nation, org_curr, exRateFile; target_curr = "USD", hhs_exp = true, hhs_info = false)

    global households, hh_list, exchange_rate, exp_curr, hh_curr
    hhs = households[year][nation]
    hhl = hh_list[year][nation]
    er = exchange_rate[target_curr * "/" * org_curr] = Dict{String, Float64}()
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
        if (trg, org) == (target_curr, org_curr); er[s[1]] = parse(Float64, s[3])
        elseif (trg, org) == (org_curr, target_curr); er[s[1]] = 1.0/parse(Float64, s[3])
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
            lh = length(h.date); mmidx = 5:lh; yyidx = 1:4
            if haskey(er, h.date[mmidx]); r = er[h.date[mmidx]]
            elseif haskey(er, h.date[yyidx]); r = er[h.date[yyidx]]
            else println("Exchange rate error: no exchange rate data for ", h.date)
            end
        end
        if hhs_exp && exp_chk
            for he in hhs[hh].expends
                he.value *= r
                he.valUnit = target_curr
            end
        end
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

function calculatePopWeight(year, nation, outputFile=""; ur_wgh = false, district=true, province=false, hhs_wgh = false)

    global regions, prov_list, dist_list, pops, pops_ur, pop_wgh, pop_ur_wgh
    global households, hh_list

    rl = Array{String, 1}()
    if province; append!(rl, prov_list[year][nation]) end
    if district; append!(rl, dist_list[year][nation]) end

    hl = hh_list[year][nation]
    hhs = households[year][nation]
    pop = pops[year][nation]
    pop_ur = pops_ur[year][nation]
    smp = Dict{String, Int}()           # Province sample size, {regin code, sample number}
    wgh = Dict{String, Float64}()       # Province population weight, {region code, weight}
    smp_ur = Dict{String, Tuple{Int,Int}}()             # Urban/rural province sample size, {regin code, (urban, rural)}
    wgh_ur = Dict{String, Tuple{Float64, Float64}}()    # Urban/rural province population weight, {region code, (urban, rural)}

    # count sample number
    for r in rl; smp[r] = 0 end
    for h in hl
        if province; smp[hhs[h].province] += hhs[h].size end
        if district; smp[hhs[h].district] += hhs[h].size end
    end
    if ur_wgh
        for h in hl
            if hhs[h].regtype == "urban"; smp_ur[hhs[h].province][1] += hhs[h].size
            elseif hhs[h].regtype == "rural"; smp_ur[hhs[h].province][2] += hhs[h].size
            end
        end
    end

    # calculate weights
    for r in rl; wgh[r] = pop[r]/smp[r] end
    if ur_wgh
        for r in rl, i=1:2
            if pop_ur[r][i]>0 && smp_ur[r][i]>0; wgh_ur[r][i] = pop_ur[r][i] / smp_ur[r][i] end
        end
    end
    if !haskey(pop_wgh, year); pop_wgh[year] = Dict{String, Dict{String, Float64}}() end
    pop_wgh[year][nation] = wgh
    if ur_wgh && !haskey(pop_ur_wgh, year)
        pop_ur_wgh[year] = Dict{String, Dict{String, Tuple{Float64, Float64}}}()
        pop_ur_wgh[year][nation] = wgh_ur
    end

    # assign household weight
    if hhs_wgh
        for h in hl
            if ur_wgh
                if hhs[h].regtype == "urban"; iur = 1 elseif hhs[h].regtype == "rural"; iur = 2
                elseif hhs[h].regtype=="sparce"; iur=1 elseif hhs[h].regtype=="intermediate"; iur=2 elseif hhs[h].regtype=="dense"; iur=3
                end
            end
            if district; ur_wgh ? hhs[h].popwgh = wgh_ur[hhs[h].district][iur] : hhs[h].popwgh = wgh[hhs[h].district]
            elseif province; ur_wgh ? hhs[h].popwgh = wgh_ur[hhs[h].province][iur] : hhs[h].popwgh = wgh[hhs[h].province]
            end
        end
    end

    # print population weights
    if length(outputFile)>0
        f = open(outputFile, "w")
        println(f, "Code\tRegion\tPopulation_weight");
        for r in rl; println(f, r, "\t", regions[year][nation][r], "\t", pop_wgh[year][nation][r]) end
        close(f)
    end
end

function scalingExpByCPI(year, nation, cpiCodeFile, statFile, linkFile, targetYear; period="year", region="district", revHH=true, revMat=false)
    # scaling hh expenditure to the target year based on Consumer Price Index
    # period: "year" or "month"
    # region: "province" or "district"

    global hh_list, sc_list, prov_list, dist_list, dist_prov, households, pops, pops_ur, expMatrix
    hl = hh_list[year][nation]
    sl = sc_list[year][nation]
    hhs = households[year][nation]
    pop = pops[year][nation]
    # pop_ur = pops_ur[year][nation]
    ds_pr = dist_prov[year][nation]
    ces_sectors = sectors[year][nation]     # CES/HBS micro-data sectors: {code, commodity}

    cpi_sec = Array{String, 1}()            # CPI sector code list
    cpi_reg = Array{String, 1}()            # CPI region code list
    cpi_per = Array{String, 1}()            # CPI data period list: XXXX = year, XXXXMM = year(XXXX), month(MM)
    ces_cpi_link = Dict{String, String}()   # CES/HBS - CPI code matching: {micro-data sector code, CPI sector code}
    cpi_sectors = Dict{String, String}()    # CPI sectors: {code, sector}
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

    # read CPI-CES/HBS linkages
    val_sep = getValueSeparator(linkFile)
    f = open(linkFile)
    readline(f)
    ces_cod_chk = Array{String, 1}()
    for l in eachline(f)
        s = string.(strip.(split(l, val_sep)))
        ces_cpi_link[s[3]] = s[1]
        if !haskey(cpi_sectors, s[1]); println("CPI sectors do not contain code ", s[1]) end
        if !haskey(ces_sectors, s[3]); println("CES/HBS sectors do not contain code ", s[3])
        elseif strip(ces_sectors[s[3]].sector) != strip(s[4])
            println("CES/HBS sector does not match with ", s[3], "\t", s[4], "\t", ces_sectors[s[3]].sector)
        end
        if !(s[3] in ces_cod_chk); push!(ces_cod_chk, s[3]) else println(s[3], "Duplicated CES/HBS sector code: ", s[3]) end
    end
    if !issubset(sort(ces_cod_chk), sort(sl)); println("CSI code matching list's CES/HBS codes doen't match with micro-data's") end
    close(f)

    # read CPIs
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
    pr_code, ds_code = prov_list[year][nation], dist_list[year][nation]
    npr, nds = length(pr_code), length(ds_code)
    scl_yr_nt, scl_mth_nt = zeros(Float64, ns), zeros(Float64, ns, nm)
    scl_yr_pr, scl_mth_pr = zeros(Float64, npr, ns), zeros(Float64, npr, ns, nm)
    scl_yr_ds, scl_mth_ds = zeros(Float64, nds, ns), zeros(Float64, nds, ns, nm)
    cpi_dist = filter(x->x in ds_code, cpi_reg)

    for i = 1:npr
        prc = pr_code[i]
        pidx = findfirst(x->x==prc, cpi_reg)
        if pidx != nothing
            for j = 1:ns; scl_yr_pr[i, j] = scl_rate_yr[pidx, j] end
            for j = 1:ns, k = 1:nm; scl_mth_pr[i, j, k] = scl_rate_mth[pidx, j, k] end
        else
            pr_dist = filter(x->ds_pr[x]==prc, cpi_dist)
            pr_pop = sum([pop[ds] for ds in pr_dist])
            didx = [findfirst(x->x==ds, cpi_reg) for ds in pr_dist]
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
    ntidx = findfirst(x->x==nation, cpi_reg)
    if ntidx != nothing
        for i = 1:ns; scl_yr_nt[i] = scl_rate_yr[ntidx, i] end
        for i = 1:ns, j = 1:nm; scl_mth_nt[i, j] = scl_rate_mth[ntidx, i, j] end
    else
        pr_pop = [haskey(pop, prc)&&pop[prc]>0 ? pop[prc] : sum([pop[ds] for ds in filter(x->ds_pr[x]==prc, cpi_dist)]) for prc in pr_code]
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
        didx = findfirst(x->x==dsc, cpi_reg)
        if didx != nothing
            for j = 1:ns; scl_yr_ds[i, j] = scl_rate_yr[didx, j] end
            for j = 1:ns, k = 1:nm; scl_mth_ds[i, j, k] = scl_rate_mth[didx, j, k] end
        else
            pidx = findfirst(x->x==ds_pr[dsc], pr_code)
            if pidx == nothing; println("Province code matching error: ", dsc, ", ", ds_pr[dsc]) end
            for j = 1:ns; scl_yr_ds[i, j] = scl_yr_pr[pidx, j] end
            for j = 1:ns, k = 1:nm; scl_mth_ds[i, j, k] = scl_mth_pr[pidx, j, k] end
        end
    end

    # scaling expenditure matrix
    if revMat; nsl = length(sl); em = expMatrix[year][nation] end
    sl_ex = filter(x->haskey(ces_cpi_link, x), sl)
    sec_idx = Dict(c => findfirst(x->x==ces_cpi_link[c], cpi_sec) for c in sl_ex)
    if period == "year"; if region == "province"; scl = scl_yr_pr; elseif region =="district"; scl = scl_yr_ds end
    elseif period == "month"; if region == "province"; scl = scl_mth_pr; elseif region =="district"; scl = scl_mth_ds end
    end
    for i = 1:length(hl)
        hh = hhs[hl[i]]
        if region == "province"; reg_idx = findfirst(x->x== hh.province, pr_code)
        elseif region =="district"; reg_idx = findfirst(x->x== hh.district, ds_code)
        end
        if period == "month" && length(hh.date)>=6; midx = findfirst(x->x==hh.date[1:6], current_mths) else midx = 1 end
        if revHH; for he in hh.expends; if he.code in sl_ex; he.value *= scl[reg_idx, sec_idx[he.code], midx] end end end
        if revMat; for j=1:nsl; if sl[j] in sl_ex; em[i, j] *= scl[reg_idx, sec_idx[sl[j]], midx] end end end
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
    for h in hl; for e in hhs[h].expends; e.value *= scl_rate end end
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
        print(f, "Code\tCode_State/Province\tState/Province\tPopulation\tWeight")
    elseif region == "district"
        rl, dp = dist_list[year][nation], dist_prov[year][nation]
        print(f, "Code\tCode_State/Province\tState/Province\tCode_District/City\tDistrict/City\tPopulation\tWeight")
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
    println(f, "Code\tSector\tMain_category\tUnit")
    for s in sl; println(f, s, "\t", sc[s].sector, "\t", sc[s].category, "\t", sc[s].unit); count += 1 end
    close(f)
    println("$count commodities' data is printed.")
end

function printHouseholdData(year, nation, outputFile; prov_wgh=false, dist_wgh=false, ur_dist = false, surv_date = false)

    global households, hh_list, regions, pop_wgh, pop_ur_wgh, exp_curr, exp_period

    hhs = households[year][nation]
    rg = regions[year][nation]
    pw = pop_wgh[year][nation]
    if ur_dist; pwur = pop_ur_wgh[year][nation] end
    exp_unit = exp_curr[year][nation][1] * "/" * exp_period[year][nation][1]

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    count = 0

    print(f, "HHID\tCode_province/state\tProvince/State\tCode_district/city\tDistrict/City")
    if ur_dist; print(f, "\tUrban/Rural") end
    if surv_date; print(f,"\tSurvey_date") end
    print(f, "\tHH size\tTotal_exp\tTot_exp_unit\tAgg_exp\tAgg_exp_unit")
    if prov_wgh; print(f, "\tProv_PopWgh") end
    if dist_wgh; print(f, "\tDist_PopWgh") end
    println(f)
    for h in hh_list[year][nation]
        hh = hhs[h]
        if !haskey(rg, hh.province) && haskey(rg, hh.district); println(hh.province, "\t", hh.district) end
        print(f, hh.hhid, "\t", hh.province, "\t", rg[hh.province], "\t", hh.district, "\t", rg[hh.district])
        if ur_dist; print(f, "\t", hh.regtype); if hh.regtype == "urban"; uridx=1 elseif hh.regtype == "rural"; uridx=2 end end
        if surv_date; print(f, "\t", hh.date) end
        print(f, "\t", hh.size, "\t", hh.totexp, "\t", hh.unit, "\t", hh.aggexp, "\t", exp_unit)
        if prov_wgh; if ur_dist; print(f, "\t", pwur[hh.province][uridx]) else print(f, "\t", pw[hh.province]) end end
        if dist_wgh; if ur_dist; print(f, "\t", pwur[hh.district][uridx]) else print(f, "\t", pw[hh.district]) end end
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
    mat = [expMatrix[year][nation]]
    if quantity; push!(mat, qntMatrix[year][nation]) end
    row, col = hh_list[year][nation], sc_list[year][nation]
    nr, nc = length(row), length(col)
    nre, nce = length(rowErr), length(colErr)

    f_tag = ["_exp", "_qnt"]
    f_sep = getValueSeparator(outputFile)
    f_name = rsplit(outputFile, '.', limit=2)
    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    for i = 1:length(mat)
        m = mat[i]
        f = open(f_name[1] * f_tag[i] * "." * f_name[2], "w")
        for c in col; print(f, f_sep, c) end
        if nre > 0; print(f, f_sep, "Row_error") end
        println(f)
        for ri = 1:nr
            print(f, row[ri])
            for ci = 1:nc; print(f, f_sep, m[ri,ci]) end
            if nre > 0; print(f, f_sep, rowErr[ri]) end
            println(f)
        end
        if nce > 0
            print(f, "Column error")
            for ci = 1:nc; print(f, f_sep, colErr[ci]) end
            print(f, f_sep, sum(colErr))
        end
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

function readPrintedRegionData(year, nation, inputFile)

    global regions, prov_list, dist_list, dist_prov, pops, pops_ur, pop_wgh, pop_ur_wgh
    essential = ["Code", "Code_State/Province", "State/Province", "Code_District/City", "District/City"]
    optional = ["Population", "Weight"]
    ur_title = ["Pop_urban", "Pop_rural", "Wgh_urban", "Wgh_rural"]

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = string.(strip.(split(readline(f), f_sep)))
    if issubset(essential, title)
        i = [findfirst(x->x==item, title) for item in essential]
        io = [findfirst(x->x==item, title) for item in optional]
        iu = [findfirst(x->x==item, title) for item in ur_title]
        op_chk = all(io.!=nothing)
        ur_chk = all(iu.!=nothing)
    else println(inputFile, " household file does not contain all essential data.")
    end
    if !haskey(regions, year)
        regions[year] = Dict{String, Dict{String, String}}()
        pops[year] = Dict{String, Dict{String, Float64}}()
        pop_wgh[year] = Dict{String, Dict{String, Float64}}()
        prov_list[year] = Dict{String, Array{String, 1}}()
        dist_list[year] = Dict{String, Array{String, 1}}()
        dist_prov[year] = Dict{String, Dict{String, String}}()
        if ur_chk
            pops_ur[year] = Dict{String, Dict{String, Tuple{Float64, Float64}}}()
            pop_ur_wgh[year] = Dict{String, Dict{String, Tuple{Float64, Float64}}}()
        end
    end
    regions[year][nation] = Dict{String, String}()
    pops[year][nation] = Dict{String, Float64}()
    pop_wgh[year][nation] = Dict{String, Float64}()
    prov_list[year][nation] = Array{String, 1}()
    dist_list[year][nation] = Array{String, 1}()
    dist_prov[year][nation] = Dict{String, String}()
    if ur_chk
        pops_ur[year][nation] = Dict{String, Tuple{Float64, Float64}}()
        pop_ur_wgh[year][nation] = Dict{String, Tuple{Float64, Float64}}()
    end

    count = 0
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        r_cd = s[i[1]]
        push!(prov_list[year][nation], s[i[2]])
        push!(dist_list[year][nation], s[i[4]])
        dist_prov[year][nation][s[i[4]]] = s[i[2]]
        rg = regions[year][nation]
        if !haskey(rg, s[i[2]]) rg[s[i[2]]] = s[i[3]] end
        if !haskey(rg, s[i[4]]) rg[s[i[4]]] = s[i[5]] end
        if op_chk
            pops[year][nation][r_cd] = parse(Float64, s[io[1]])
            pop_wgh[year][nation][r_cd] = parse(Float64, s[io[2]])
        end
        if ur_chk
            pops_ur[year][nation][r_cd] = (parse(Float64, s[iu[1]]), parse(Float64, s[iu[2]]))
            pop_ur_wgh[year][nation][r_cd] = (parse(Float64, s[iu[3]]), parse(Float64, s[iu[4]]))
        end
        count += 1
    end
    close(f)
    print(" read $count regions")
end

function readPrintedHouseholdData(year, nation, inputFile)

    global households, hh_list, regions, hh_curr, hh_period, exp_curr, exp_period
    essential = ["HHID", "Code_province/state", "Province/State", "Code_district/city", "District/City", "HH size"]
    optional = ["Total_exp", "Agg_exp", "Urban/Rural", "Survey_date", "Tot_exp_unit", "Prov_PopWgh", "Dist_PopWgh", "Agg_exp_unit"]

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = string.(strip.(split(readline(f), f_sep)))
    if issubset(essential, title); i = [findfirst(x->x==t, title) for t in [essential;optional]]
    else println(inputFile, " household file does not contain all essential data.")
    end
    if !haskey(households, year)
        households[year] = Dict{String, Dict{String, household}}()
        hh_list[year] = Dict{String, Array{String, 1}}()
    end
    if !haskey(hh_curr, year); hh_curr[year] = Dict{String, Array{String, 1}}() end
    if !haskey(hh_period, year); hh_period[year] = Dict{String, Array{String, 1}}() end
    if !haskey(exp_curr, year); exp_curr[year] = Dict{String, Array{String, 1}}() end
    if !haskey(exp_period, year); exp_period[year] = Dict{String, Array{String, 1}}() end


    rg = regions[year][nation]
    hhs = households[year][nation] = Dict{String, household}()
    hl = hh_list[year][nation] = Array{String, 1}()
    hhc = hh_curr[year][nation] = Array{String, 1}()
    hhp = hh_period[year][nation] = Array{String, 1}()
    exp_c = exp_curr[year][nation] = Array{String, 1}()
    exp_p = exp_period[year][nation] = Array{String, 1}()

    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        hhid = s[i[1]]
        push!(hl, hhid)
        hhs[hhid] = household(hhid)
        hh_vals = ["","","","",0,0,"","","",0,0,0,0,"",0,[],[]]
        hh_vals[[2,3,5,10]] = [s[i[2]], s[i[4]], parse(Int, s[i[6]]), parse(Float64, s[i[7]])]
        for oi = 9:11; if isa(i[oi], Int); hh_vals[[4, 1, 14][oi-8]] = s[i[oi]] end end
        if isa(i[13], Int); hh_vals[15] = parse(Float64, s[i[13]]) elseif isa(i[12], Int); hh_vals[15] = parse(Float64, s[i[12]]) end
        reviseHouseholdData(year, nation, hhid, hh_vals)
        currency, period = string.(strip.(split(s[i[11]], '/')))
        if !(currency in hhc); push!(hhc, currency) end
        if !(period in hhp); push!(hhp, period) end
        currency, period = string.(strip.(split(s[i[14]], '/')))
        if !(currency in exp_c); push!(exp_c, currency) end
        if !(period in exp_p); push!(exp_p, period) end
    end
    close(f)
    print(" read ", length(hl), " households")
end

function readPrintedMemberData(year, nation, inputFile)

    global households, hh_list
    essential = ["HHID", "Age", "Gender"]

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = string.(strip.(split(readline(f), f_sep)))
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
    print(" read $count members")
end

function readPrintedExpenditureData(year, nation, inputFile; quantity = false)

    global households, hh_list, exp_curr, exp_period, pr_unts
    essential = ["HHID", "Commodity_code", "Expended_value", "Currency_unit", "Consumption_period"]
    if quantity; essential = [essential[1:4]; ["Consumed_quantity", "Quantity_unit"] ; essential[5:end]] end
    if !haskey(exp_curr, year) exp_curr[year] = Dict{String, Array{String, 1}}() end
    if !haskey(exp_period, year) exp_period[year] = Dict{String, Array{String, 1}}() end
    if !haskey(exp_curr[year], nation) exp_curr[year][nation] = Array{String, 1}() end
    if !haskey(exp_period[year], nation) exp_period[year][nation] = Array{String, 1}() end
    e_curr, e_per = exp_curr[year][nation], exp_period[year][nation]

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = string.(strip.(split(readline(f), f_sep)))
    if issubset(essential, title); i = [findfirst(x->x==et, title) for et in essential]
    else println(inputFile, " expenditure file does not contain all essential data.")
    end
    hhs = households[year][nation]

    count = 0
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        hhid = s[i[1]]
        exp_vals = [s[i[2]], parse(Float64, s[i[3]]), 0, s[i[4]], "", parse(Int16, s[i[end]])]
        if quantity; exp_vals[[3, 5]] = [parse(Float64, s[i[5]]), s[i[6]]] end

        appendExpenditureData(year, nation, hhid, exp_vals)
        if !(exp_vals[4] in e_curr); push!(e_curr, exp_vals[4]) end
        if !(pr_unts[exp_vals[6]] in e_per); push!(e_per, pr_unts[exp_vals[6]]) end
        count += 1
    end
    close(f)
    print(" read $count expenditures")
end

function readPrintedSectorData(year, nation, itemfile)

    global sc_list, sectors
    if !haskey(sc_list, year); sc_list[year] = Dict{String, Array{String, 1}}() end
    if !haskey(sectors, year); sectors[year] = Dict{String, Dict{String, commodity}}() end
    sl = sc_list[year][nation] = Array{String, 1}()
    sc = sectors[year][nation] = Dict{String, commodity}()

    essential = ["Code", "Sector", "Main_category", "Unit"]
    f_sep = getValueSeparator(itemfile)
    count = 0
    f = open(itemfile)
    title = string.(strip.(split(readline(f), f_sep)))
    if issubset(essential, title); i = [findfirst(x->x==et, title) for et in essential]
    else println(itemfile, " commodity sector file does not contain all essential data.")
    end
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        push!(sl, s[i[1]])
        sc[s[i[1]]] = commodity(s[i[1]], s[i[2]], s[i[3]], s[i[4]], "")
        count += 1
    end
    close(f)
    print(" read $count sectors")
end

function readPrintedExpenditureMatrix(year, nation, inputFile)

    global hh_list, sc_list, expMatrix

    if !haskey(expMatrix, year); expMatrix[year] = Dict{String, Array{Float64, 2}}() end
    sl, hl = sc_list[year][nation], hh_list[year][nation]
    ns, nh = length(sl), length(hl)
    em = expMatrix[year][nation] = zeros(Float64, nh, ns)

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    codes = string.(strip.(split(readline(f), f_sep)[2:end]))
    if issubset(sl, codes); i = [findfirst(x->x==sc, codes) for sc in sl]
    else println(inputFile, " expenditure matrix file does not contain all essential data.")
    end

    hhs = Array{String, 1}()
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        hi = findfirst(x->x==s[1], hl)
        vals = [parse(Float64, v) for v in s[2:1+ns]]
        if hi !=nothing
            push!(hhs, s[1])
            em[hi,:] = vals[i]
        end
    end
    close(f)
    if sort(hhs) != sort(hl); println(inputFile, " expenditure matrix file does not contain all household data.") end
end

end
