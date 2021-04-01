module MicroDataReader

# Developed date: 17. Mar. 2021
# Last modified date: 30. Mar. 2021
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
    coicop::String      # commodity corresponding COICOP code

    commodity(cod, sec="", cat="", subcat="", coi="") = new(cod, sec, cat, subcat, coi)
end

global households = Dict{Int, Dict{String, Dict{String, household}}}()  # household dict: {year, {nation A3, {hhid, household}}}
global sectors = Dict{Int, Dict{String, Dict{String, commodity}}}()     # expenditure sector: {year, {nation A3, {code, commodity}}}
global hh_list = Dict{Int, Dict{String, Array{String, 1}}}()            # hhid list: {year, {nation A3, {hhid}}}
global sc_list = Dict{Int, Dict{String, Array{String, 1}}}()            # commodity code list: {year, {nation A3, {code}}}

global regions = Dict{Int, Dict{String, Dict{String, String}}}()        # expenditure sector: {year, {nation A3, {code, region}}}
global prov_list = Dict{Int, Dict{String, Array{String, 1}}}()          # province code list: {year, {nation A3, {code}}}
global dist_list = Dict{Int, Dict{String, Array{String, 1}}}()          # district code list: {year, {nation A3, {code}}}

global pops = Dict{Int, Dict{String, Dict{String, Float64}}}()          # population: {year, {nation, {region_code, population}}}
global pops_ur = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}()   # urban/rural population: {year, {nation, {region_code, (urban, rural)}}
global pop_wgh = Dict{Int, Dict{String, Dict{String, Float64}}}()       # population weight: {year, {nation, {region_code, weight}}}
global pop_ur_wgh = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}()    # urban/rural population weight: {year, {nation, {region_code, (urban, rural)}}

global expMatrix = Dict{Int, Dict{String, Array{Float64, 2}}}()         # expenditure matrix: {year, {nation, {hhid, commodity}}}

function appendCommoditySectorData(year, nation, cmm_data)
    global sectors, sc_list

    push!(sc_list[year][nation], cmm_data[1])
    sectors[year][nation][cmm_data[1]] = commodity(cmm_data[1], cmm_data[2], cmm_data[3], "", cmm_data[4])
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


function readIndexFile(year, nation, indexFile; err_display=false)

    titles = Array{String, 1}()
    idxs = Array{Array{String, 1}, 1}()         # essential index values
    idxs_opt = Array{Array{String, 1}, 1}()     # optional index values
    idxstart = idxend = metaend = titlend = optchk = yrchk = natchk = false
    fext = indexFile[findlast(isequal('.'), indexFile)+1:end]
    if fext == "csv"; val_sep = ',' elseif fext == "tsv" || fext == "txt"; val_sep = '\t' end
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

    global regions, prov_list, dist_list

    if !haskey(regions, year); regions[year] = Dict{String, Dict{String, String}}() end
    if !haskey(regions[year], nation); regions[year][nation] = Dict{String, String}() end
    if !haskey(prov_list, year); prov_list[year] = Dict{String, Array{String, 1}}() end
    if !haskey(dist_list, year); dist_list[year] = Dict{String, Array{String, 1}}() end

    rg = regions[year][nation]
    pl = Array{String, 1}()
    dl = Array{String, 1}()

    yrchk, natchk, idxs = readIndexFile(year, nation, regionFile)
    for s in idxs
        rg[s[1]] = s[6]
        if !haskey(rg, s[3]); rg[s[3]] = s[4] end
        if !(s[1] in dl); push!(dl, s[1]) end
        if !(s[3] in pl); push!(pl, s[3]) end
    end

    prov_list[year][nation] = sort(pl)
    dist_list[year][nation] = sort(dl)
end

function readMicroData(year, nation, microdataPath, hhIdxFile, mmIdxFile, cmmIdxFile, expIdxFile)

    for idxfile in filter(x->length(x)>0, [hhIdxFile, mmIdxFile, cmmIdxFile, expIdxFile])
        # read microdata index file
        yrchk, natchk, idxs, ipdx_o = readIndexFile(year, nation, idxfile)

        # read microdata contents
        if idxfile == hhIdxFile; readHouseholdData(year, nation, [idxs;ipdx_o], microdataPath)  # idxs: {sector, position, type, file, tag}
        elseif idxfile == mmIdxFile; readMemberData(year, nation, idxs, microdataPath)          # idxs: {sector, position, type, file, tag}
        elseif idxfile == cmmIdxFile; readCommoditySectors(year, nation, idxs)   # idxs: {code, coicop, sector, entity, category, sub_category}
        elseif idxfile == expIdxFile; readExpenditureData(year, nation, idxs, microdataPath)    # idxs: {}
        end
    end
end

function readHouseholdData(year, nation, indices, microdataPath; hhid_sec = "hhid")

    global households, hh_list

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

    # analyze index data
    mdFiles = Dict{String, Dict{String, Tuple{Int, String, String}}}()  # {microdata_file, {data_sector, {position, type}}}
    for idx in indices      # indices: {sector, position, type, file, tag}
        if length(idx[2])>0 && length(idx[4])>0
            mf = idx[4]
            if !haskey(mdFiles, mf); mdFiles[mf] = Dict{String, Tuple{Int, String, String}}() end
            mdFiles[mf][idx[1]] = (parse(Int, idx[2]), idx[3], idx[5])      # (position, type, tag)
        end
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
        fext = mf[findlast(isequal('.'), mf)+1:end]
        if fext == "csv"; mdf_sep = ',' elseif fext == "tsv" || fext == "txt"; mdf_sep = '\t' end
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
        if !(s[1] in scl); appendCommoditySectorData(year, nation, [s; ["" for i=1:(4 - length(s))]])
        else println("Duplicated codes: ", s[1], "\t", s[3], ", ", sec[s[1]])
        end
    end
end

function readExpenditureData(year, nation, indices, microdataPath)

    global households, hh_list, sectors, sc_list

    # analyze index data
    # mdFiles: {microdata_file, {category, {data_tag, hhid_position, code_position, period(days), value_position, value_unit, quantity_position, quantity_unit}}}
    mdFiles = Dict{String, Dict{String, Tuple{String, Int, Int, Int, Int, String, Int, String}}}()
    for idx in indices      # index: {category, hhid_position, code_position, period(days), value_position, value_unit, file, data_tag, quantity_position, quantity_unit}
        valchk = true; for i in [1, 2, 3, 4, 7]; if length(idx[i])==0; chk = false end end
        if valchk
            if length(idx) == 7 || length(idx) == 8
                li = findlast(x->tryparse(Float64, x) != nothing, idx)
                if li == 5; idx = [idx; ["0", ""]]
                elseif li == 7; idx = [idx[1:4]; ["0", ""] ;idx[5:end]]
                end
            elseif length(idx) == 10
                idxs = [[idx[1], parse(Int, idx[2]), parse(Int, idx[3]), parse(Int, idx[4])]; idx[5:end]]
                for i in [5, 9]; if length(idxs[i]) > 0; idxs[i] = parse(Int, idxs[i]); else idxs[i] = 0 end end
                mf = idxs[7]
                if !haskey(mdFiles, mf); mdFiles[mf] = Dict{String, Tuple{Int, Int, Int, Int, String, Int, String}}() end
                mdFiles[mf][idxs[1]] = tuple([idxs[8]; idxs[2:6]; idxs[9:end]]...)
            else println("Expenditure index content length error: ", year, ", ", nation, "\t", idx)
            end
        else println("Expenditure index lacks essential contents: ", year, ", ", nation, "\t", idx)
        end
    end
    mfs = sort(collect(keys(mdFiles)))

    # for m in mfs
    #     for s in collect(keys(mdFiles[m]))
    #         println(m,"\t",s,"\t",mdFiles[m][s])
    #     end
    # end

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
                f = open(microdataPath * mf)
                readline(f)     # read title line
                for l in eachline(f)
                    s = strip.(split(l, mdf_sep))
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
                    # for numeric values
                    exp_idx = [2, 3]
                    mfd_idx = [5, 7]
                    for i = 1:length(mfd_idx)
                        val = parse(Float64, s[mfd[sc][mfd_idx[i]]])
                        if val > 0; exp_vals[exp_idx[i]] = val end
                    end
                    # for period-value
                    if mfd[sc][4] > 0; exp_vals[6] = mfd[sc][4] end

                    appendExpenditureData(year, nation, hhid, exp_vals)
                end
                close(f)
            end
        end
    end
end

function buildExpenditureMatrix(year, nation, outputFile = ""; transpose = false, period = 365, print_err = false)
    # build an expenditure matrix as period (init: 365) days consumption monetary values
    # [row]: household, [column]: commodity

    global households, hh_list, sc_list, expMatrix
    if !haskey(expMatrix, year); expMatrix[year] = Dict{String, Array{Float64, 2}}() end

    hhs = households[year][nation]
    row = hh_list[year][nation]
    col = sc_list[year][nation]
    nr = length(row)
    nc = length(col)

    mat = zeros(Float64, nr, nc)
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
                else rowErr[ri] += 1; colErr[ci] += 1
                end
            end
        end
        h.aggexp = total
    end
    if transpose
        mat = transpose(mat)
        row, nr, rowErr, col, nc, colErr = col, nc, colErr, row, nr, rowErr
    end
    expMatrix[year][nation] = mat

    # print expenditure matrix
    if length(outputFile) > 0
        fext = outputFile[findlast(isequal('.'), outputFile)+1:end]
        if fext == "csv"; f_sep = ',' elseif fext == "tsv" || fext == "txt"; f_sep = '\t' end
        f = open(outputFile, "w")
        for c in col; print(f, f_sep, c) end
        if print_err; print(f, f_sep, "Row_error") end
        println(f)
        for ri = 1:nr
            print(f, row[ri])
            for ci = 1:nc; print(f, f_sep, mat[ri,ci]) end
            if print_err; print(f, f_sep, rowErr[ri]) end
            println(f)
        end
        if print_err
            print(f, "Column error")
            for ci = 1:nc; print(f, f_sep, colErr[ci]) end
            print(f, f_sep, sum(colErr))
        end
        println(f)
        close(f)
    end

    return mat, row, col, rowErr, colErr
end

function exchangeExpCurrency(year, nation, exchangeRate; inverse=false)
    # exchangeRate: can be a file path that contains excahnge rates, a constant value of
    #               a nation's currency to USD (normally) currency exchange rate (USD/A3), or a set of values of Dict[MMYY] or Dict[YY]

    global households, hh_list
    hhs = households[year][nation]
    hhl = hh_list[year][nation]

    # read exchange rate from the recieved file if 'exchangeRate' is 'String'
    if typeof(exchangeRate) <:AbstractString
        fext = exchangeRate[findlast(isequal('.'), exchangeRate)+1:end]
        if fext == "csv"; f_sep = ',' elseif fext == "tsv" || fext == "txt"; f_sep = '\t' end
        er = Dict{String, Float64}()
        f = open(exchangeRate)
        readline(f)
        for l in eachline(f); s = split(l, f_sep); er[s[1]] = parse(Float64, s[2]) end
        close(f)
        if inverse; for x in collect(keys(er)); er[x] = 1/er[x] end end
        exchangeRate = er
    end

    # exchange the expenditure currency
    if typeof(exchangeRate) <: Number
        for hh in hhl; for he in hhs[hh].expends; he.value *= exchangeRate end end
    elseif typeof(exchangeRate) <: AbstractDict
        yr = string(year)
        if !(yr in collect(keys(exchangeRate)))
            rates = [exchangeRate[mm] for mm in filter(x->(length(x)==6 && x[1:4]==yr), collect(keys(exchangeRate)))]
            exchangeRate[yr] = sum(rates) / length(rates)
        end
        for hh in hhl
            if length(hhs[hh].date)==0; er = exchangeRate[yr]
            else
                lh = length(h.date); mmidx = 5:lh; yyidx = 1:4
                if haskey(exchangeRate, h.date[mmidx]); er = exchangeRate[h.date[mmidx]]
                elseif haskey(exchangeRate, h.date[yyidx]); er = exchangeRate[h.date[yyidx]]
                else println("Exchange rate error: no exchange rate data for ", h.date)
                end
            end
            for he in hhs[hh].expends; he.value *= er end
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
        fext = exchangeRate[findlast(isequal('.'), exchangeRate)+1:end]
        if fext == "csv"; f_sep = ',' elseif fext == "tsv" || fext == "txt"; f_sep = '\t' end
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

function calculatePopWeight(year, nation, outputFile=""; ur_wgh = false, district=true, province=false)

    global regions, prov_list, dist_list, pops, pops_ur, pop_wgh
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

    # print population weights
    if length(outputFile)>0
        f = open(outputFile, "w")
        println(f, "Code\tRegion\tPopulation_weight");
        for r in rl; println(f, r, "\t", regions[year][nation][r], "\t", pop_wgh[year][nation][r]) end
        close(f)
    end
end

function printHouseholdData(year, nation, outputFile; prov_wgh=false, dist_wgh=false, ur_dist = false)

    global households, hh_list, regions, pop_wgh, pop_ur_wgh

    hhs = households[year][nation]
    rg = regions[year][nation]
    pw = pop_wgh[year][nation]
    if ur_dist; pwur = pop_ur_wgh[year][nation] end

    f = open(outputFile, "w")
    count = 0

    print(f, "HHID\tProvince/State\tDistrict/City")
    if ur_dist; print(f, "\tUrban/Rural") end
    print(f, "\tHH size\tTotal_exp\tAgg_exp")
    if prov_wgh; print(f, "\tProv_PopWgh") end
    if dist_wgh; print(f, "\tDist_PopWgh") end
    println(f)
    for h in hh_list[year][nation]
        hh = hhs[h]


        if !haskey(rg, hh.province) && haskey(rg, hh.district); println(hh.province, "\t", hh.district) end


        print(f, hh.hhid , "\t", rg[hh.province], "\t", rg[hh.district])
        if ur_dist; print(f, hh.regtype); if hh.regtype == "urban"; uridx=1 elseif hh.regtype == "rural"; uridx=2 end end
        print(f, "\t", hh.size, "\t", hh.totexp, "\t", hh.aggexp)
        if prov_wgh; if ur_dist; print(f, "\t", pwur[hh.province][uridx]) else print(f, "\t", pw[hh.province]) end end
        if dist_wgh; if ur_dist; print(f, "\t", pwur[hh.district][uridx]) else print(f, "\t", pw[hh.district]) end end
        println(f)
        count += 1
    end
    close(f)

    println("$count households' data is printed.")
end

function printExpenditureData(year, nation, outputFile)

    global households, hh_list
    hhs = households[year][nation]
    hl = hh_list[year][nation]

    f = open(outputFile, "w")
    count = 0

    println(f, "HHID\tCommodity_code\tExpended_value\tConsumption_periody")
    for h in hl
        hh = hhs[h]
        for e in hh.expends
            println(f, hh.hhid,"\t",e.code,"\t",e.value,"\t",e.period)
            count += 1
        end
    end
    close(f)

    println("$count items' data is printed.")
end

function initVars()
    global households = Dict{Int, Dict{String, Dict{String, household}}}()
    global sectors = Dict{Int, Dict{String, Dict{String, commodity}}}()
end

end
