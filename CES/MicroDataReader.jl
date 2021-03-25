hhmodule MicroDataReader

# Developed date: 17. Mar. 2021
# Last modified date: 24. Mar. 2021
# Subject: Household consumption expenditure survey microdata reader
# Description: read consumption survey microdata and store household, member, and expenditure data
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

mutable struct expenditure
    code::String        # product or service item code
    value::Float64      # total consumption value
    quantity::Float64   # total consumption quantity
    valUnit::String     # total value's currency unit; ex. 'USD' or 'EUR'
    qntUnit::String     # total quantity's unit; ex. 'kg', 'litre', or 'number'
    period::Int16       # total value's consuming period; generally 7, 30, or 365 days

    expenditure(c, v=0, q=0, vu="", qu="", p=0) = new(c, v, q, vu, qu, p)
end

mutable struct member
    hhid::String    # household id
    age::Int16      # age
    gen::Int8       # gender: [1]male, [2]female, [9]not available
    nat::String     # nationality (A3)
    mar::Int8       # marital status: [1]single, [2]married, [9]not available
    edu::Int8       # education level: [1]non literate, [2]literate, [3]below primary, [4]primary, [5]below secondary, [6]secondary,
                    #                  [7]higher secondary, [8]tertiary, [8]graduate, [9]post gradutate/above, [10]not specified, [99]not available
    rel::Int8       # relationship with head: (of reference person and/or of the spouse) [1]reference person, [2]spouse or partner,
                    #                         [3]child, [4]parent, [5]relative, [6]no family relationship, [9]not available,
    occ::String     # occupation
    inc::Float64    # income amount
    incUnit::String # income unit: ex. 'USD/month'

    member(id, ag=0, ge=0, na="", ma=0, ed=0, re=0, oc="", in=0, iu="") = new(id, ag, ge, na, ma, ed, re, oc, in, iu)
end

mutable struct household
    hhid::String        # household identification no.
    date::String        # survey date
    state::String       # state code
    province::String    # province code
    district::String    # district code
    city::String        # city code
    regtype::String     # region type: urban/rural or sparce/intermediate/dense
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
    popWgh::Float64     # population weight: (urban/rural) province/state/distric/city population represented by one member

    members::Array{member,1}        # household member(s)
    expends::Array{expenditure,1}   # consumed product or sevice items

    household(hi,da="",sa="",pr="",di="",ct="",rt="",sz=0,ag=0,rl="",oc="",ed="",te=0,ti=0,tep=0,tip=0,ut="",pw=0,mm=[],ex=[]) = new(hi,da,sa,pr,di,ct,rt,sz,ag,rl,oc,ed,te,ti,tep,tip,ut,pw,mm,ex)
end

mutable struct commodity
    code::String        # commodity code in the survey data
    sector::String      # commodity's sector or label
    entity::String      # propduct or service
    category::String    # main category, ex) Food, Electricity, Gas, Other energy, Public transport, Private transport, Medical care, Education, Consumable goods, Durable goods, Other services
    subCategory::String # sub-category, ex) Food related: Grain, Vegetable, Fruit, Dairy, Beef, Pork, Poultry, Other meat, Fish, Alcohol, Other beverage, Confectionery, Restaurant, Other food, etc
                        #                   Energy-related: Electricity, Gas, Wood, Dung cake, Kerosene, Coal, Petrol, Diesel, Biogas, Other fuel, etc.
                        #                   Transport-related: Road (private), Road (public), Rail, Air, Water, Other, etc.
    coicop::String      # commodity corresponding COICOP code

    commodity(cod, sec="", ent="", cat="", subcat="", coi="") = new(cod, sec, ent, cat, subcat, coi)
end

global households = Dict{Int, Dict{String, Dict{String, household}}}()  # household dict: {year, {nation A3, {hhid, household}}}
global sectors = Dict{Int, Dict{String, Dict{String, commodity}}}()        # expenditure sector: {year, {nation A3, {code, commodity}}}
global hh_list = Dict{Int, Dict{String, Array{String, 1}}}()            # hhid list: {year, {nation A3, {hhid}}}
global sc_list = Dict{Int, Dict{String, Array{String, 1}}}()            # commodity code list: {year, {nation A3, {code}}}

function appendCommoditySectorData(year, nation, cmm_data)
    global sectors, sc_list

    push!(sc_list[year][nation], cmm_data[1])
    sectors[year][nation][cmm_data[1]] = commodity(cmm_data[1], cmm_data[2], cmm_data[3], cmm_data[4], cmm_data[5], cmm_data[6])
end

function appendExpenditureData(year, nation, id, exp_value)
    global households
    hh = households[year][nation][id]

    push!(hh.expends, expenditure(exp_value[1], exp_value[2], exp_value[3], exp_value[4], exp_value[5], exp_value[6]))
end

function appendMemberData(year, nation, id, mm_value)
    global households
    hh = households[year][nation][id]

    push!(hh.members, member(id, mm_value[1],mm_value[2],mm_value[3],mm_value[4],mm_value[5],mm_value[6],mm_value[7],mm_value[8],mm_value[9])
end

function reviseHouseholdData(year, nation, id, hh_value)
    global households
    hh = households[year][nation][id]

    if length(id)>0 && length(hh_value)==19
        da = hh_value[1]; if length(da)>0; hh.date = da end
        sa = hh_value[2]; if length(sa)>0; hh.state = sa end
        pr = hh_value[3]; if length(pr)>0; hh.province = pr end
        di = hh_value[4]; if length(di)>0; hh.district = di end
        ct = hh_value[5]; if length(ct)>0; hh.city = ct end
        rt = hh_value[6]; if length(rt)>0; hh.regtype = rt end
        sz = hh_value[7]; if sz>0; hh.size = sz end
        ag = hh_value[8]; if ag>0; hh.age = ag end
        rl = hh_value[9]; if length(rl)>0; hh.rel = rl end
        oc = hh_value[10]; if length(oc)>0; hh.occ = oc end
        ed = hh_value[11]; if length(ed)>0; hh.edu = ed end
        te = hh_value[12]; if te>0; hh.totexp = te end
        ti = hh_value[13]; if ti>0; hh.totinc = ti end
        tep = hh_value[14]; if tep>0; hh.totexppc = tep end
        tip = hh_value[15]; if tip>0; hh.totincpc = tip end
        ut = hh_value[16]; if length(ut)>0; hh.unit = ut end
        pw = hh_value[17]; if pw>0; hh.popWgh = pw end
        mm = hh_value[18]; if length(mm)>0; hh.members = mm end
        ex = hh_value[19]; if length(ex)>0; hh.expends = ex end
    else println("HHID is absent or incomplete HH values.")
    end
end

function readMicroData(year, nation, microdataPath, hhIdxFile, mmIdxFile, cmmIdxFile, expIdxFile)

    for idxfile in [hhIdxFile, mmIdxFile, cmmIdxFile, expIdxFile]
        # read microdata index file
        idxs = Array{Array{String, 1}, 1}()
        idxstart = idxend = metaend = yrchk = natchk = false
        fext = idxfile[findlast(isequal('.'), idxfile)+1:end]
        if fext == "csv"; idxf_sep = ',' elseif fext == "tsv" || fext == "txt"; idxf_sep = '\t' end
        f = open(idxfile)
        for l in eachline(f)
            if startswith(l) == "========================"; if !idxstart; idxstart = true else idxend = true end
            elseif startswith(l) == "------------------------" && idxstart; metaend = true
            elseif idxstart && !idxend
                s = strip.(split(l, idxf_sep))
                tag = lowercase(filter(x->x!=':', s[1]))
                if !metaend
                    if tag == 'year' && parse(Int, s[2]) == year; yrchk = true
                    elseif tag == 'nation' && lowercase(s[2]) == lowercase(nation); natchk = true
                    elseif tag == 'a3' && uppercase(s[2]) == nation; natchk = true
                    end
                else
                    push!(idxs, [tag; s[2:end]])
                    if length(s) < 5; push!(idxs[end], "") end
                end

            end
        end
        if !idxstart || !idxend || !metaend; println("Wrong index file format: ", indexFile) end
        if !yrchk; println("Year does not match: ", year, ", ", indexFile) end
        if !natchk; println("Nation does not match: ", nation, ", ", indexFile) end
        close(f)

        # read microdata contents
        if idxfile == hhIdxFile; readHouseholdData(year, nation, idxs)          # idxs: {sector, position, type, file, tag}
        elseif idxfile == mmIdxFile; readMemberData(year, nation, idxs)         # idxs: {sector, position, type, file, tag}
        elseif idxfile == cmmIdxFile; readCommoditySectors(year, nation, idxs)  # idxs: {code, coicop, sector, entity, category, sub_category}
        elseif idxfile == expIdxFile; readExpenditureData(year, nation, idxs)   # idxs: {}
        end
    end
end

function readHouseholdData(year, nation, indices; hhid_sec = "hhid")

    global households, hh_list

    sectors = ["survey_date", "state", "province", "district", "city", "region_type", "hh_size", "head_age", "head_religion", "head_occupation", "head_education", "expenditure", "income", "exp_percap", "inc_percap", "currency_unit", "pop_weight"]
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

    # analyze index data
    mdFiles = Dict{String, Dict{String, Tuple{Int, String, String}}}()  # {microdata_file, {data_sector, {position, type}}}
    for idx in indices      # indices: {sector, position, type, file, tag}
        mf = idx[4]
        if !haskey(mdFiles, mf); mdFiles[mf] = Dict{String, Tuple{Int, String, String}}() end
        mdFiles[mf][idx[1]] = (parse(Int, idx[2]), idx[3], idx[5])      # (position, type, tag)
    end
    mfs = sort(collect(keys(mdFiles)))
    for mf in mfs; if !haskey(mdFiles[mf], hhid_sec); println(mf, "does not contain HHID sector.") end end

    # read household data: all the microdata files should contain HHID values
    for mf in mfs
        mfd = mdFiles[mf]
        hhid_pos = mfd[hhid_sec][1]
        hhid_tag = mfd[hhid_sec][3]
        fext = mf[findlast(isequal('.'), mf)+1:end]
        if fext == "csv"; mdf_sep = ',' elseif fext == "tsv" || fext == "txt"; mdf_sep = '\t' end
        f = open(mf)
        for l in eachline(f)
            s = strip.(split(l, mdf_sep))
            hhid = hhid_tag * s[hhid_pos]
            if !haskey(households[year][nation], hhid)
                households[year][nation][hhid] = household[hhid]
                push!(hh_list[year][nation], hhid)
            end
            hh_vals = ["","","","","","",0,0,"","","",0,0,0,0,"",0,[],[]]
            for i = 1:nsec
                c = sectors[i]
                if c in mfd
                    val = s[mfd[c][1]]
                    if lenght(val) > 0
                        if c in int_sec; hh_vals[i] = parse(Int, val)
                        elseif c in flo_sec; hh_vals[i] = parse(Float64, val)
                        else hh_vals[i] = val
                        end
                    end
                end
            end
            reviseHouseholdData(year, nation, hhid, hh_vals)
        end
        close(f)
    end
end

function readMemberData(year, nation, indices; hhid_sec = "hhid")

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
        fext = mf[findlast(isequal('.'), mf)+1:end]
        if fext == "csv"; mdf_sep = ',' elseif fext == "tsv" || fext == "txt"; mdf_sep = '\t' end
        f = open(mf)
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
        if !(s[1] in scl); appendCommoditySectorData(year, nation, [s; ["" for i=1:(6 - length(s))]])
        else println("Duplicated codes: ", s[1], "\t", s[3], ", ", sec[s[1]])
        end
    end
end

function readExpenditureData(year, nation, indices; build_table=true)

    global households, hh_list, sectors, sc_list

    # analyze index data
    # mdFiles: {microdata_file, {category, {data_tag, hhid_position, code_position, period(days), value_position, value_unit, quantity_position, quantity_unit}}}
    mdFiles = Dict{String, Dict{String, Tuple{String, Int, Int, Int, Int, String, Int, String}}}()
    for idx in indices      # index: {category, hhid_position, code_position, period(days), value_position, value_unit, file, data_tag, quantity_position, quantity_unit}
        if length(idx) == 8
            li = findlast(x->tryparse(Float64, x) != nothing, idx)
            if li == 5; idx = [idx; ["0", ""]]
            elseif li == 7; idx = [idx[1:4]; "0", "", ;idx[5:end]]
            end
        end
        valchk = true; for i in [1, 2, 3, 4, 7]; if length(idx[i])==0; chk = false end end
        if length(idx) == 10 && valchk
            idxs = [[idx[1], parse(Int, idx[2]), parse(Int, idx[3]), parse(Int, idx[4])]; idx[5:end]]
            for i in [5, 9]; if length(idxs[i]) > 0; idxs[i] = parse(Int, idxs[i]); else idxs[i] = 0 end end
            mf = idxs[7]
            if !haskey(mdFiles, mf); mdFiles[mf] = Dict{String, Tuple{Int, Int, Int, Int, String, Int, String}}() end
            mdFiles[mf][idxs[1]] = tuple([idxs[8]; idxs[2:6]; idxs[9:end]]...)
        else println("Expenditure index content error: ", year, ", ", nation)
        end
    end
    mfs = sort(collect(keys(mdFiles)))

    # read expenditure data: all the microdata files should contain HHID and Code values
    for mf in mfs
        mfd = mdFiles[mf]
        sec = collect(keys(mfd))
        fext = mf[findlast(isequal('.'), mf)+1:end]
        if fext == "csv"; mdf_sep = ',' elseif fext == "tsv" || fext == "txt"; mdf_sep = '\t' end

        pre_mfd = []
        for sc in sec
            dup_chk = false
            for pre_sc in pre_mfd; if mfd[sc] == mfd[pre_sc]; dup_chk = true end end
            if !dup_chk
                push!(pre_mfd, sc)
                hhid_tag, hhid_pos, code_pos = mfd[sc][1], mfd[sc][2], mfd[sc][3]
                # {data_tag, hhid_position, code_position, period(days), value_position, value_unit, quantity_position, quantity_unit}
                f = open(mf)
                for l in eachline(f)
                    s = strip.(split(l, mdf_sep))
                    hhid = hhid_tag * s[hhid_pos]
                    exp_vals = ["", 0, 0, "", "", 0]
                    for i = 1:3
                        exp_idx = [1, 4, 5]
                        mfd_idx = [3, 6, 8]
                        val = s[mfd[sc][mfd_idx[i]]]
                        if length(val) > 0; exp_vals[exp_idx[i]] = val end
                    end
                    for i = 1:3
                        exp_idx = [2, 3, 6]
                        mfd_idx = [5, 7, 4]
                        val = s[mfd[sc][mfd_idx[i]]]
                        if val > 0; exp_vals[exp_idx[i]] = val end
                    end
                    appendExpenditureData(year, nation, hhid, exp_vals)
                end
                close(f)
            end
        end
    end
end

function buildExpenditureMatrix(outputFile = "")
    # build expenditure matrix with 365 days consumption monetary values

    row = sort(collect(keys(households)))
    col = sort(collect(keys(categories)))

    mat = zeros(Float64, length(row), length(col))
    rowErr = zeros(Int32, length(row))
    colErr = zeros(Int32, length(col))

    # make expenditure matrix, row: households, col: expenditure categories
    for r = 1:length(row)
        h = households[row[r]]
        total = 0
        for i in h.items
            if haskey(categories, i.code)
                c = findfirst(x -> x==i.code, col)
                if i.value > 0
                    val = i.value
                    if i.period != 365; val *= 365 / i.period end
                    mat[r,c] += val
                    total += val
                else rowErr[r] += 1; colErr[c] += 1
                end
            end
        end
        h.totExp = total
        h.totExpMrp = h.mpceMrp * h.size * 12
    end

    # print expenditure matrix as a file
    if length(outputFile) > 0
        f = open(outputFile, "w")

        for c in col; print(f, "\t", c) end
        print(f, "\tRow_error")
        println(f)
        for r = 1:length(row)
            print(f, row[r])
            for c = 1:length(col); print(f, "\t", mat[r,c]) end
            print(f, "\t", rowErr[r])
            println(f)
        end
        print(f, "Column error")
        for c = 1:length(col); print(f, "\t", colErr[c]) end
        print(f, "\t", sum(colErr))
        println(f)
        close(f)
    end

    return mat, row, col, rowErr, colErr
end

function exchangeExpCurrency(exchangeRate; inverse=false)
                                # exchangeRate: can be a file path that contains excahnge rates,
                                #               a constant value of Rupees to USD currency exchange rate
                                #               , or a set of values of Dict[MMYY] or Dict[YY]
    global households

    # read exchange rate from the recieved file if 'exchangeRate' is 'String'
    if typeof(exchangeRate) <:AbstractString
        er = Dict{String, Float64}()
        f = open(exchangeRate)
        readline(f)
        for l in eachline(f); s = split(l, '\t'); er[s[1]] = parse(Float64,s[2]) end
        close(f)
        if inverse; for x in collect(keys(er)); er[x] = 1/er[x] end end
        exchangeRate = er
    end

    # if 'exchangeRate' is a constant rate
    if typeof(exchangeRate) <: Number
        for h in collect(values(households))
            for i in h.items; i.value *= exchangeRate; i.homeVal *= exchangeRate end
        end
    # if 'exchangeRate' is a set of rates
    elseif typeof(exchangeRate) <: AbstractDict
        stdRate = exchangeRate[minimum(filter(x->length(x)==2, collect(keys(exchangeRate))))]
        for h in collect(values(households))
            if length(h.date)==0; er = stdRate
            else
                lh = length(h.date); idxmm = lh-3:lh; idxyy = lh-1:lh
                if haskey(exchangeRate, h.date[idxmm]); er = exchangeRate[h.date[idxmm]]
                elseif haskey(exchangeRate, h.date[idxyy]); er = exchangeRate[h.date[idxyy]]
                else println("Exchange rate error: no exchange rate data for ", h.date)
                end
            end
            for i in h.items; i.value *= er; i.homeVal *= er end
        end
    end
end

function convertMpceToPPP(pppRate; inverse=false)
                                    # PPP rates: can be a file path that contains PPP rates,
                                    #            a constant value of India Rupees/USD PPP converting rate
                                    #            , or a set of values of Dict[MMYY] or Dict[YY]
    global households

    # read converting rate from the recieved file if 'pppFile' is 'String'
    if typeof(pppRate) <:AbstractString
        ppp = Dict{String, Float64}()
        f = open(pppRate)
        readline(f)
        for l in eachline(f); s = split(l, '\t'); ppp[s[1]] = parse(Float64,s[2]) end
        close(f)
        pppRate = ppp
    end

    if typeof(pppRate) <: Number; for h in collect(values(households)); h.mpceMrp /= pppRate end
    elseif typeof(pppRate) <: AbstractDict
        stdRate = pppRate[minimum(filter(x->length(x)==2, collect(keys(pppRate))))]
        for h in collect(values(households))
            if length(h.date)==0; h.mpceMrp /= stdRate
            else
                lh = length(h.date); idxmm = lh-3:lh; idxyy = lh-1:lh
                if haskey(ppp, h.date[idxmm]); h.mpceMrp /= pppRate[h.date[idxmm]]
                elseif haskey(ppp, h.date[idxyy]); h.mpceMrp /= pppRate[h.date[idxyy]]
                else println("Exchange rate error: no exchange rate data for ", h.date)
                end
            end
        end
    end
end

function readCategory(inputFile)
    f= open(inputFile)

    readline(f)
    for l in eachline(f)
        s = split(l, '\t')
        categories[s[1]] = s[2]
    end

    close(f)
end

function calculateStatePopulationWeight(populationFile)

    global households

    stalist = Array{String, 1}()                    # State code list
    stapop = Dict{String, Tuple{Int, Int, Int}}()   # State population, {State code, population{total, rural, urban}}
    stasmp = Dict{String, Array{Int, 1}}()          # State sample size, {State code, sample number{total, rural, urban}}
    stawgh = Dict{String, Array{Float64, 1}}()      # State population weight, {State code, population{total, rural, urban}}

    # read population data
    f = open(populationFile)
    readline(f)
    for l in eachline(f)
        s = split(l, ",")
        stapop[s[1]] = (parse(Int, s[4]), parse(Int, s[6]), parse(Int, s[8]))
        stasmp[s[1]] = zeros(Int, 3)
    end
    close(f)
    stalist = sort(collect(keys(stapop)))

    # count sample number
    for h in collect(values(households))
        if h.sector == "1"; stidx=1         # rural
        elseif h.sector == "2"; stidx=2     # urban
        else println("HH sector error: not \"urban\" nor \"rural\"")
        end
        stasmp[h.state][1] += h.size
        stasmp[h.state][stidx+1] += h.size
    end

    # calculate weights
    for st in stalist
        stawgh[st] = zeros(Float64, 3)
        for i=1:3; stawgh[st][i] = stapop[st][i]/stasmp[st][i] end
    end

    for h in collect(values(households))
        if h.sector=="1"; h.popStaWgh = stawgh[h.state][2]
        elseif h.sector=="2"; h.popStaWgh = stawgh[h.state][3]
        else println("Household ",h," sector is wrong")
        end
    end
end

function calculateDistrictPopulationWeight(populationFile, concordanceFile)

    global households

    disList = Array{String, 1}()                # Survey district code list
    dispop = Dict{String, Array{Int, 1}}()      # District population, {Survey district code, population{total, rural, urban}}
    dissmp = Dict{String, Array{Int, 1}}()      # District sample size, {Survey district code, sample number{total, rural, urban}}
    diswgh = Dict{String, Array{Float64, 1}}()  # District population weight, {Survey district code, population{total, rural, urban}}

    disConc = Dict{String, String}()     # Population-Survey concordance, {Statistics district code, Survey district code}

    # read concordance data
    xf = XLSX.readxlsx(concordanceFile)
    sh = xf["IND_pop"]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r)>1 && !ismissing(r[3]); disConc[string(r[2])]=string(r[3]) end
    end
    disList = sort(collect(values(disConc)))

    # read population data
    f = open(populationFile)
    readline(f)
    for l in eachline(f)
        s = split(l, ",")
        discode = string(s[2])
        if haskey(disConc, discode)
            if string(s[4])=="Total"; dispop[disConc[discode]] = [parse(Int, s[6]), 0, 0]
            elseif string(s[4])=="Rural"; dispop[disConc[discode]][2] = parse(Int, s[6])
            elseif string(s[4])=="Urban"; dispop[disConc[discode]][3] = parse(Int, s[6])
            else println(discode, " does not have \"Total\", \"Urban\" nor \"Rural\" data")
            end
        end
    end
    close(f)
    for dc in disList; dissmp[dc] = zeros(Int, 3) end

    # count sample number
    for h in collect(values(households))
        if h.sector == "1"; stidx=1         # rural
        elseif h.sector == "2"; stidx=2     # urban
        else println("HH sector error: not \"urban\" nor \"rural\"")
        end
        dissmp[h.district][1] += h.size
        dissmp[h.district][stidx+1] += h.size
    end

    # calculate weights
    for dc in disList
        diswgh[dc] = zeros(Float64, 3)
        diswgh[dc][1] = dispop[dc][1]/dissmp[dc][1]
        if dissmp[dc][2]==0 && dissmp[dc][3]==0; println(dc," does not have samples in both rural and urban areas.")
        elseif dissmp[dc][2]==0; diswgh[dc][2]=0; diswgh[dc][3]=dispop[dc][1]/dissmp[dc][3]
        elseif dissmp[dc][3]==0; diswgh[dc][3]=0; diswgh[dc][2]=dispop[dc][1]/dissmp[dc][2]
        end
    end

    for h in collect(values(households))
        if h.sector=="1"; h.popDisWgh = diswgh[h.district][2]
        elseif h.sector=="2"; h.popDisWgh = diswgh[h.district][3]
        else println("Household ",h," sector is wrong")
        end
    end
end

function convertHouseholdData(outputFile = "")
    id = String[]; dat = String[]; fsu = String[]; sta = String[]; dis = String[]
    sec = String[]; siz = Int16[]; mrp = Float64[]; rel = Int8[]
    tot = Float64[]; totMrp = Float64[]

    avg = Float16[]; mal = Int32[]; fem = Int32[]; chi = Int32[]; mid = Int32[]; old = Int32[]

    for hhid in sort(collect(keys(households)))
        h = households[hhid]

        agesum = 0
        total = 0
        male = 0
        female = 0
        child = 0
        grownup = 0
        aged = 0

        for m in h.members
            total += 1
            agesum += m.age
            if m.sex == 1; male += 1
            elseif m.sex == 2; female +=1
            end
            if 0<= m.age < 15; child += 1
            elseif 15<= m.age <65; grownup += 1
            elseif m.age >= 65; aged += 1
            end
        end

        push!(id, h.id)
        push!(dat, h.date)
        push!(fsu, h.fsu)
        push!(sta, h.state)
        push!(dis, h.district)
        push!(sec, h.sector)
        push!(siz, h.size)
    #    push!(urp, h.mpceUrp)
        push!(mrp, h.mpceMrp)
        push!(rel, h.religion)
        push!(tot, h.totExp)
        push!(totMrp, h.totExpMrp)
        push!(avg, agesum / total)
        push!(mal, male)
        push!(fem, female)
        push!(chi, child)
        push!(mid, grownup)
        push!(old, aged)
    end

    df = DataFrame(HHID=id, Survey_date=dat, FSU=fsu, State=sta, District=dis, Sector=sec, Size=siz,
                    MPCE_Mrp=mrp, Religion=rel, Total_Exp=tot, Total_Exp_Mrp=totMrp,
                    Avg_age=avg, Male=mal, Female=fem, Children=chi, Grownups=mid, Aged_persons=old)

    if length(outputFile) > 0
        f = open(outputFile, "w")
        println(f, df)
        close(f)
    end

    return df
end


function printHouseholdData(outputFile; addCds=false, addPov=false, addWgh=false)
    f = open(outputFile, "w")
    sd = ["rural", "urban"]
    count = 0

    print(f, "HHID\tSurvey Date\tFSU\tState\tDistrict\tSector\tHH size\tMPCE_MRP\tReligion\tTot_exp(yr)")
    if addCds; print(f, "\tState_code\tDistrict_code\tStratum\tSubstratum_No\tFOD_Sub_Region") end
    if addPov; print(f, "\tPovLine") end
    if addWgh; print(f, "\tStaPopWeight\tDisPopWeight") end
    println(f)
    for hhid in sort(collect(keys(households)))
        h = households[hhid]
        sector = sd[parse(Int8, h.sector)]
        print(f, h.id,"\t",h.date,"\t",h.fsu,"\t",h.state,"\t",h.district,"\t",sector,"\t",h.size,"\t",h.mpceMrp,"\t",h.religion,"\t",h.totExp)
        if addCds; print(f, "\t", h.regCds[1],"\t",h.regCds[2],"\t",h.regCds[3],"\t",h.regCds[4],"\t",h.regCds[5]) end
        if addPov; if h.pov; print(f, "\tunder") else print(f, "\tover") end end
        if addWgh; print(f, "\t",h.popStaWgh,"\t",h.popDisWgh) end
        println(f)
        count += 1
    end
    close(f)

    println("$count households' data is printed.")
end

function printMemberData(outputFile)
    f = open(outputFile, "w")
    count = 0

    sd = ["male", "female"]
    md = ["never married", "currently married", "widowed", "divorced/seperated"]
    ed = ["non literate", "literate without formal schooling: through EGC/NFEC/AEC",
            "literate without formal schooling: through TLC", "literate without formal schooling: Others",
            "literate with formal schooling:below primary", "literate with formal schooling:primary",
            "literate with formal schooling:middle", "literate with formal schooling:secondary", "",
            "literate with formal schooling:higher secondary",
            "literate with formal schooling:diploma/certificate course",
            "literate with formal schooling:graduate",
            "literate with formal schooling:post graduate and above"]

    println(f, "HHID\tAge\tSex\tMarital status\tEducation")
    for hhid in sort(collect(keys(households)))
        h = households[hhid]
        for m in h.members
            if m.sex>0; sex = sd[m.sex]; else sex = ""; end
            if m.mar>0; mar = md[m.mar]; else mar = ""; end
            if m.edu>0; edu = ed[m.edu]; else edu = ""; end

            println(f, h.id,"\t",m.age,"\t",sex,"\t",mar,"\t",edu)
            count += 1
        end
    end
    close(f)

    println("$count members' data is printed.")
end

function printMicroData(outputFile)
    f = open(outputFile, "w")
    count = 0

    sd = ["only purchase", "only home-grown stock", "both purchase and home-grown stock",
          "only free collection", "only exchange of goods and services", "only gifts/charities", "", "", "others"]

    println(f, "HHID\tItem code\tTotal value\tTotal quantity\tConsumption period\tConsumption source\tHome produce value\tHome produce quantity")
    for hhid in sort(collect(keys(households)))
        h = households[hhid]
        for i in h.items
            if i.source>0; source = sd[i.source]; else source = ""; end

            println(f, h.id,"\t",i.code,"\t",i.value,"\t",i.quantity,"\t",i.period,"\t",source,"\t",i.homeVal,"\t",i.homeQt)
            count += 1
        end
    end
    close(f)

    println("$count items' data is printed.")
end

function initVars()
    global households = Dict{String, household}()
    global categories = Dict{String, String}()
end

end
