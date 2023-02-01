module GroupEstimator

# Developed date: 17. Oct. 2022
# Last modified date: 31. Jan. 2023
# Subject: Group statistics estimator
# Description: Calculate statistics of each household group
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

include("MicroDataReader.jl")

using Statistics
using .MicroDataReader
mdr = MicroDataReader

mutable struct group
    label::String    # group label
    type::String    # group classification type: "s" (String) or "v" (Value)
    sort::String    # sort of group (available only if type == "s")
    range::Tuple{Float64, Float64}  # range (=<, <) of the group (available only if type == "v")
    unit::String    # unit of the range (available only if type == "v")

    group(l, t, s="", r=(0.0,0.0)) = new(l, t, s, r)
end

yr_list = Array{Int, 1}()       # year list: {YYYY}
nat_list = Array{String, 1}()   # nation list: {A3}
cat_list = Array{String, 1}()   # category list
reg_list = Dict{Int, Dict{String, Array{String, 1}}}()  # region list: {year, {nation A3, {region (gis id)}}}
ces_gis_link = Dict{Int, Dict{String, Dict{String,String}}}() # GIS ID corresponds CES/HBS region: {year, {nation, {ces/hbs region, gis id}}}


gr_list = Dict{String, Array{String, 1}}()              # group label list: {nation A3, {label}}
groups = Dict{String, Dict{String, group}}()            # groups: {nation A3, {label, group}}
qt_list = Dict{Int, Dict{String, Array{String, 1}}}()   # questionary code list: {year, {nation A3, {code}}}
responses = Dict{Int, Dict{String, Array{Bool, 2}}}()# questionary responses: {year, {nation, {household, question}}}

hh_list = Dict{Int, Dict{String, Array{String, 1}}}()   # Household ID: {year, {nation A3, {hhid}}}
households = Dict{Int, Dict{String, Dict{String, mdr.household}}}() # household dict: {year, {nation A3, {hhid, household}}}
hhs_cf = Dict{Int, Dict{String, Array{Float64, 2}}}()   # Household categorized CF: {year, {nation A3, {household, category}}}

region_pov = Dict{Int, Dict{String, Array{Float64, 1}}}()   # region's poverty rate: {year, {nation A3, {region}}}
region_cf = Dict{Int, Dict{String, Array{Float64, 1}}}() # region CF: {year, {nation A3, {region}}}
region_cfpc = Dict{Int, Dict{String, Array{Float64, 1}}}() # region pre capita CF: {year, {nation A3, {region}}}
region_rsp = Dict{Int, Dict{String, Array{Float64, 2}}}() # region's response rates: {year, {nation A3, {region, questionary}}}
region_inc = Dict{Int, Dict{String, Array{Float64, 1}}}() # region's average income: {year, {nation A3, {region}}}
region_exp = Dict{Int, Dict{String, Array{Float64, 1}}}() # region's average expenditure: {year, {nation A3, {region}}}
region_cf_gr = Dict{Int, Dict{String, Array{Float64, 3}}}() # region's CF by group: {year, {nation A3, {region, group, (all/yes/no)}}}
region_cfpc_gr = Dict{Int, Dict{String, Array{Float64, 3}}}() # region's per capita CF by group: {year, {nation A3, {region, group, (all/yes/no)}}}
region_stt_gr = Dict{Int, Dict{String, Dict{String, Array{Float64, 3}}}}()  # region's response rates by group: {year, {nation A3, {questionary code, {region, group, (all/yes/no)}}}}
region_inc_gr = Dict{Int, Dict{String, Array{Float64, 3}}}() # region's average income by group: {year, {nation A3, {region, group, (all/yes/no)}}}
region_exp_gr = Dict{Int, Dict{String, Array{Float64, 3}}}() # region's average expenditure by group: {year, {nation A3, {region, group, (all/yes/no)}}}

period_unit = Dict("year" => 365.0, "month" => 30.0, "week" => 7.0, "day" => 1.0, "annual" => 365.0, "monthly" => 30.0, "weekly" => 7.0, "daily" => 1.0)

function importData(mdata::Module, em_data::Module; mode = "cf")
    # mode: "cf" (carbon footprint), "ie" (embedded emission), or "de" (direct emission)

    global yr_list, nat_list
    global hh_list, households = mdata.hh_list, deepcopy(mdata.households)
    global hhs_cf, region_cf, region_cfpc
    global cat_list, reg_list, ces_gis_link = em_data.cat_list, em_data.gisRegList, em_data.gisRegDist

    yrs = sort(collect(keys(hh_list)))
    for y in yrs
        if !(y in yr_list); push!(yr_list, y) end
        nats = sort(collect(keys(hh_list[y])))
        for n in nats; if !(n in nat_list); push!(nat_list, n) end end
    end
    sort!(yr_list)
    sort!(nat_list)

    if mode == "cf"
        hhs_cf = em_data.cfHHs
        region_cf = em_data.cfRegGIS
        region_cfpc = em_data.cfRegPcGIS
    elseif mode == "ie"
        hhs_cf = em_data.ieHHs
        region_cf = em_data.ieRegGIS
        region_cfpc = em_data.ieRegPcGIS
    elseif mode == "de"
        hhs_cf = em_data.deHHs
        region_cf = em_data.deRegGIS
        region_cfpc = em_data.deRegPcGIS
    else println("Wrong emission mode: ", mode)
    end
end

function setGroups(nation, threshold = [0, 1.9, 3.0, 5.0], label = ["pov", "low", "mid", "high"]; type = "v", currency = "USD", period = "day", percap = true)
    # type: group classification type: "s" (String) or "v" (Value)

    global gr_list, groups
    grl, grs = Array{String, 1}(), Dict{String, group}()

    if type == "s"
        for (tr, lb) in zip(threshold, label)
            push!(grl, lb)
            grs[lb] = group(lb, type, tr)
        end
    elseif type == "v"
        for i = 1:length(label)-1
            push!(grl, label[i])
            grs[label[i]] = group(label[i], type, "", (threshold[i], threshold[i+1]))
        end
        push!(grl, label[end])
        grs[label[end]] = group(label[end], type, "", (threshold[end], 0))
        unit_label = currency * "/" * period * (percap ? "/percap" : "")
        for gr in grs; gr.unit = unit_label end
    else println("Incorrect group threshold type: ", type)
    end

    gr_list[nation] = grl
    groups[nation] = grs
end

function estimateRegionState(year, nation; group_mode = true, poverty_label = "pov", weight_mode = true, region_mode = "district", income_mode = true)
    # region_mode: "district" or "province"

    global hh_list, reg_list, gr_list, qt_list, period_unit, cat_list
    global groups, households, responses, hhs_cf
    global region_pov, region_cf, region_inc, region_exp, region_cf_gr, region_inc_gr, region_exp_gr

    y, n = year, nation
    rl, hhl, hhs, cg_lnk = reg_list[y][n], hh_list[y][n], households[y][n], ces_gis_link[y][n]
    pr, grl, gtl = period_unit, gr_list[n], qt_list[y][n]
    hcf = hhs_cf[y][n]
    nr, nh, ng, nq = length(rl), length(hhl), length(grl), length(qtl)
    w_t, w_rg, w_gr = 0.0, zeros(Float64, nr), zeros(Float64, nr, ng, 3)
    cf, inc, exp, pvr = zeros(Float64, nr), zeros(Float64, nr), zeros(Float64, nr), zeros(Float64, nr)
    cf_gr, inc_gr, exp_gr = zeros(Float64, nr, ng, nq, 3), zeros(Float64, nr, ng, nq, 3), zeros(Float64, nr, ng, nq, 3)

    if !haskey(region_cf, y); region_cf[y] = Dict{String, Array{Float64, 1}}() end
    if !haskey(region_inc, y); region_inc[y] = Dict{String, Array{Float64, 1}}() end
    if !haskey(region_exp, y); region_exp[y] = Dict{String, Array{Float64, 1}}() end
    if !haskey(region_pov, y); region_pov[y] = Dict{String, Array{Float64, 1}}() end

    ci = findfirst(x -> x == "total", cat_list)
    p_l = groups[nation][poverty_label].range[2]

    for hi = 1:nh
        h = hhl[hi]
        hh = hhs[h]
        sz = hh.size
        if region_mode == "district"; r = cg_lnk[hh.district]
        elseif region_mode == "province"; r = cg_lnk[hh.province]
        else println("\nIncorrect region_mode: ", region_mode)
        end
        ri = findfirst(x -> x == r, rl)
        w = (weight_mode ? hh.popwgh : 1.0)
        w_sz = w * sz
        w_t += w_sz
        w_rg[ri] += w_sz
        cf[ri] += hcf[hi, ci] * w
        inc[ri] += hh.totinc * w
        exp[ri] += hh.totexp * w
        pvr[ri] += ((income_mode ? hh.totinc : hh.totexp) < p_l ? w_sz : 0.0)


    end
    cf ./= w_rg
    inc ./= w_rg
    exp ./= w_rg
    pvr ./= w_rg

    region_cf[y][n] = cf
    region_inc[y][n] = inc
    region_exp[y][n] = exp
    region_pov[y][n] = pvr
end

function readResponseData(year, nation, inputFile; qst_label =[], qst_res = [], hhid_label = "", district_mode = true, weight_mode = true)

    global hh_list, qt_list, reg_list, households, responses, region_rsp, ces_gis_link
    y, n = year, nation
    hhs, cg_lnk = households[y][n], ces_gis_link[y][n]

    if !haskey(region_rsp, y); region_rsp[y] = Dict{String, Array{Float64, 2}}() end
    if !haskey(qt_list, y); qt_list[y] = Dict{String, Array{String, 1}}() end
    qt_list[y][n] = qst_label[:]

    nr, nq = length(reg_list[y][n]), length(qst_label)
    r_idx = Dict([reg_list[y][n][i] => i for i = 1:nr])

    reg_rate = zeros(Float64, nr, nq)
    reg_samp = zeros(Float64, nr)

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = string.(strip.(split(readline(f), f_sep)))
    hi = findfirst(x -> x == hhid_label, title)
    qi = [findfirst(x -> x == ql, title) for ql in qst_label]
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        hh = hhs[s[hi]]
        ri = r_idx[cg_lnk[(district_mode ? hh.district : hh.province)]]
        w = (weight_mode ? hh.size * hh.weight : 1.0)
        for qi = 1:nq
            qs_er = qst_res[i]
            if s[qi[i]] in (isa(qs_er, Array{}) ? qs_er : [qs_er]); reg_rate[ri, qi] += w end
        end
        reg_samp[ri] += w
    end
    close(f)
    reg_rate ./= reg_samp

    region_rsp[y][n] = reg_rate
end

function printResponseStatistics(year, nation, outputFile)

    global reg_list, qt_list, region_rsp
    y, n = year, nation
    rgl, qtl, reg_rate = reg_list[y][n], qt_list[y][n], region_rsp[y][n]
    rn, qn = length(rgl), length(qtl)

    f_sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")
    for qt in qtl; print(f, f_sep, qt) end; println(f)
    for i = 1:rn
        print(f, rgl[i])
        for j = 1:qn; print(f, f_sep, reg_rate[i,j]) end
        println(f)
    end
    println(f)

    close(f)
end

function readQuestionaryCodes(inputFile; code_label = "code")

    global qt_list

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = string.(strip.(split(readline(f), f_sep)))
    ci = findfirst(x->x==code_label, title)
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        push!(qt_list, s[ci])
    end
    close(f)
end

# function readGroupData(inputFile)
#     # "Type" = "s" (String) or "v" (Value)
#
#     global gr_list
#     sectors = ["Group", "Code", "Type", "Range", "Unit"]
#
#     f_sep = getValueSeparator(inputFile)
#     f = open(inputFile)
#     title = string.(strip.(split(readline(f), f_sep)))
#     gi, ci, ti, ri, ui = [findfirst(x->x==label, title) for label in sectors]
#     for l in eachline(f)
#         s = string.(strip.(split(l, f_sep)))
#         gr, tp = s[gi], lowercase(s[ti])
#         if !haskey(gr_list, gr); gr_list[gr] = Array{group, 1}() end
#         push!(gr_list[gr], group(s[ci], tp, (tp == "s" ? s[ri] : ""), (tp == "v" ?  : (0.0,0.0))))
#
#     end
#     close(f)
#
# end

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

end
