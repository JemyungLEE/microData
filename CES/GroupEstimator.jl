# SPDX-FileCopyrightText: © 2022 Jemyung Lee <jemyung81@gmail.com>
# SPDX-License-Identifier: GPL-3.0

module GroupEstimator

# Developed date: 17. Oct. 2022
# Last modified date: 7. Feb. 2023
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
gr_thrs = Dict{Int, Dict{String, Array{Float64, 1}}}()  # group thresholds: {year, {nation, {threshold}}}
gr_unit = Dict{Int, Dict{String, String}}()             # unit of group thresholds: year, {nation, unit}}
groups = Dict{String, Dict{String, group}}()            # groups: {nation A3, {label, group}}
qt_list = Dict{Int, Dict{String, Array{String, 1}}}()   # questionary code list: {year, {nation A3, {code}}}
responses = Dict{Int, Dict{String, Array{Bool, 2}}}()   # questionary responses: {year, {nation, {household, question}}}

hh_list = Dict{Int, Dict{String, Array{String, 1}}}()   # Household ID: {year, {nation A3, {hhid}}}
households = Dict{Int, Dict{String, Dict{String, mdr.household}}}() # household dict: {year, {nation A3, {hhid, household}}}
hhs_cf = Dict{Int, Dict{String, Array{Float64, 2}}}()   # Household categorized CF: {year, {nation A3, {household, category}}}

region_pov = Dict{Int, Dict{String, Array{Float64, 1}}}()   # region's poverty rate: {year, {nation A3, {region}}}
region_rsp = Dict{Int, Dict{String, Array{Float64, 2}}}()   # region's response rates: {year, {nation A3, {region, questionary}}}
region_cf = Dict{Int, Dict{String, Array{Float64, 1}}}()    # region per capita CF: {year, {nation A3, {region}}}
region_inc = Dict{Int, Dict{String, Array{Float64, 1}}}()   # region's average per capita income: {year, {nation A3, {region}}}
region_exp = Dict{Int, Dict{String, Array{Float64, 1}}}()   # region's average per capita expenditure: {year, {nation A3, {region}}}
region_rsp_gr = Dict{Int, Dict{String, Array{Float64, 3}}}()    # region's response rates by group: {year, {nation A3, {region, group, questionary}}}
region_cf_gr = Dict{Int, Dict{String, Array{Float64, 4}}}()     # region's per capita CF by group: {year, {nation A3, {region, group, questionary, (all/yes/no)}}}
region_inc_gr = Dict{Int, Dict{String, Array{Float64, 4}}}()    # region's average per capita income by group: {year, {nation A3, {region, group, questionary, (all/yes/no)}}}
region_exp_gr = Dict{Int, Dict{String, Array{Float64, 4}}}()    # region's average per capita expenditure by group: {year, {nation A3, {region, group, questionary, (all/yes/no)}}}

nation_pov = Dict{Int, Dict{String, Float64}}() # nation's poverty rate: {year, {nation A3, rate}}
nation_rsp = Dict{Int, Dict{String, Array{Float64, 1}}}()   # nation's response rates: {year, {nation A3, {questionary}}}
nation_cf = Dict{Int, Dict{String, Float64}}()  # nation per capita CF: {year, {nation A3, cf}}
nation_inc = Dict{Int, Dict{String, Float64}}() # nation's average per capita income: {year, {nation A3, income}}
nation_exp = Dict{Int, Dict{String, Float64}}() # nation's average per capita expenditure: {year, {nation A3, expenditure}}
nation_rsp_gr = Dict{Int, Dict{String, Array{Float64, 2}}}()    # nation's response rates by group: {year, {nation A3, {group, questionary}}}
nation_cf_gr = Dict{Int, Dict{String, Array{Float64, 3}}}()     # nation's per capita CF by group: {year, {nation A3, {group, questionary, (all/yes/no)}}}
nation_inc_gr = Dict{Int, Dict{String, Array{Float64, 3}}}()    # nation's average per capita income by group: {year, {nation A3, {group, questionary, (all/yes/no)}}}
nation_exp_gr = Dict{Int, Dict{String, Array{Float64, 3}}}()    # nation's average per capita expenditure by group: {year, {nation A3, {group, questionary, (all/yes/no)}}}


period_unit = Dict("year" => 365.0, "month" => 30.0, "week" => 7.0, "day" => 1.0, "annual" => 365.0, "monthly" => 30.0, "weekly" => 7.0, "daily" => 1.0)

function importData(mdata::Module, em_data::Module; mode = "cf")
    # mode: "cf" (carbon footprint), "ie" (embedded emission), or "de" (direct emission)

    global yr_list, nat_list, hhs_cf
    global hh_list, households = mdata.hh_list, deepcopy(mdata.households)
    global cat_list, reg_list, ces_gis_link = em_data.cat_list, em_data.gisRegList, em_data.gisRegDist

    yrs = sort(collect(keys(hh_list)))
    for y in yrs
        if !(y in yr_list); push!(yr_list, y) end
        nats = sort(collect(keys(hh_list[y])))
        for n in nats; if !(n in nat_list); push!(nat_list, n) end end
    end
    sort!(yr_list)
    sort!(nat_list)

    if mode == "cf"; hhs_cf = em_data.cfHHs
    elseif mode == "ie"; hhs_cf = em_data.ieHHs
    elseif mode == "de"; hhs_cf = em_data.deHHs
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

function setGroupThresholds(year, nation, threshold = [1.9, 3.0, 5.0], label = ["pov", "low", "mid", "high"]; unit = "USD/day/cap")

    global gr_list[nation] = label
    global gr_thrs[year] = Dict(nation => threshold)
    global gr_unit[year] = Dict(nation => unit)
end

function convertGroupThresholds(year, nation; conv_unit = "USD/year/cap")

    global gr_thrs, gr_unit, period_unit

    grt, gru = gr_thrs[year][nation], gr_unit[year][nation]
    or_crr, or_prd, or_scl = string.(split(gru, '/'))
    co_crr, co_prd, co_scl = string.(split(conv_unit, '/'))

    if or_crr != co_crr
        # not prepared
        println("\nCurrency change mode is not prepared.")
    end
    if or_prd != co_prd; grt .*= period_unit[co_prd] / period_unit[or_prd] end
    if or_scl != co_scl; println("Scale of convert units are different: ", co_scl, " (convert), ", or_scl, "(current)") end

    gr_thrs[year][nation], gr_unit[year][nation] = grt, conv_unit
end

function estimateRegionState(year, nation; group_mode = true, weight_mode = true, region_mode = "district", income_mode = true, category = "total")
    # region_mode: "district" or "province"

    global hh_list, reg_list, gr_list, qt_list, period_unit, cat_list, gr_thrs
    global households, responses, hhs_cf
    global region_pov, region_cf, region_inc, region_exp
    global region_rsp, region_rsp_gr, region_cf_gr, region_inc_gr, region_exp_gr

    y, n = year, nation
    rl, hhl, hhs, cg_lnk = reg_list[y][n], hh_list[y][n], households[y][n], ces_gis_link[y][n]
    pr, grl, qtl, grt, hhs_rsp = period_unit, gr_list[n], qt_list[y][n], gr_thrs[y][n], responses[y][n]
    hcf = hhs_cf[y][n]
    nr, nh, ng, nq = length(rl), length(hhl), length(grl), length(qtl)
    w_t, w_rg, w_grp = 0.0, zeros(Float64, nr), zeros(Float64, nr, ng, nq, 3)
    cf, inc, exp, pvr = zeros(Float64, nr), zeros(Float64, nr), zeros(Float64, nr), zeros(Float64, nr)
    rsp_rg, rsp_grp = zeros(Float64, nr, nq), zeros(Float64, nr, ng, nq)
    cf_grp, inc_grp, exp_grp = zeros(Float64, nr, ng, nq, 3), zeros(Float64, nr, ng, nq, 3), zeros(Float64, nr, ng, nq, 3)

    if !haskey(region_cf, y); region_cf[y] = Dict{String, Array{Float64, 1}}() end
    if !haskey(region_inc, y); region_inc[y] = Dict{String, Array{Float64, 1}}() end
    if !haskey(region_exp, y); region_exp[y] = Dict{String, Array{Float64, 1}}() end
    if !haskey(region_pov, y); region_pov[y] = Dict{String, Array{Float64, 1}}() end
    if !haskey(region_rsp, y); region_rsp[y] = Dict{String, Array{Float64, 2}}() end
    if !haskey(region_rsp_gr, y); region_rsp_gr[y] = Dict{String, Array{Float64, 3}}() end
    if !haskey(region_cf_gr, y); region_cf_gr[y] = Dict{String, Array{Float64, 4}}() end
    if !haskey(region_inc_gr, y); region_inc_gr[y] = Dict{String, Array{Float64, 4}}() end
    if !haskey(region_exp_gr, y); region_exp_gr[y] = Dict{String, Array{Float64, 4}}() end

    ci = findfirst(x -> lowercase(x) == lowercase(category), cat_list)
    p_l = grt[1]

    for hi = 1:nh
        h = hhl[hi]
        hh = hhs[h]
        sz = hh.size
        h_inc = (income_mode ? hh.totinc : hh.totexp)
        gi = (h_inc >= grt[end] ? ng : findfirst(x -> x > h_inc, grt))
        if region_mode == "district"; r = cg_lnk[hh.district]
        elseif region_mode == "province"; r = cg_lnk[hh.province]
        else println("\nIncorrect region_mode: ", region_mode)
        end
        ri = findfirst(x -> x == r, rl)
        w = (weight_mode ? hh.popwgh : 1.0)
        w_sz = w * sz
        w_cf = w * hcf[hi, ci]
        w_inc = w * hh.totinc
        w_exp = w * hh.totexp

        w_t += w_sz
        w_rg[ri] += w_sz
        cf[ri] += w_cf
        inc[ri] += w_inc
        exp[ri] += w_exp
        pvr[ri] += (h_inc < p_l ? w_sz : 0.0)

        for qi = 1:nq
            h_rsp = hhs_rsp[hi, qi]
            rsp_idx = [1, (h_rsp ? 2 : 3)]
            w_grp[ri, gi, qi, rsp_idx] .+= w_sz
            rsp_rg[ri, qi] += (h_rsp ? w_sz : 0.0)
            rsp_grp[ri, gi, qi] += (h_rsp ? w_sz : 0.0)
            cf_grp[ri, gi, qi, rsp_idx] .+= w_cf
            inc_grp[ri, gi, qi, rsp_idx] .+= w_inc
            exp_grp[ri, gi, qi, rsp_idx] .+= w_exp
        end
    end
    cf ./= w_rg
    inc ./= w_rg
    exp ./= w_rg
    pvr ./= w_rg
    rsp_rg ./= w_rg
    rsp_grp ./= w_grp[:, :, :, 1]
    cf_grp ./= w_grp
    inc_grp ./= w_grp
    exp_grp ./= w_grp

    region_cf[y][n], region_inc[y][n], region_exp[y][n], region_pov[y][n] = cf, inc, exp, pvr
    region_rsp[y][n], region_rsp_gr[y][n] = rsp_rg, rsp_grp
    region_cf_gr[y][n], region_inc_gr[y][n], region_exp_gr[y][n] = cf_grp, inc_grp, exp_grp
end

function estimateNationState(year, nation; group_mode = true, weight_mode = true, income_mode = true, category = "total")
    # region_mode: "district" or "province"

    global hh_list, reg_list, gr_list, qt_list, period_unit, cat_list, gr_thrs
    global households, responses, hhs_cf
    global nation_pov, nation_cf, nation_inc, nation_exp
    global nation_rsp, nation_rsp_gr, nation_cf_gr, nation_inc_gr, nation_exp_gr

    y, n = year, nation
    hhl, hhs = hh_list[y][n], households[y][n]
    pr, grl, qtl, grt, hhs_rsp = period_unit, gr_list[n], qt_list[y][n], gr_thrs[y][n], responses[y][n]
    hcf = hhs_cf[y][n]
    nh, ng, nq = length(hhl), length(grl), length(qtl)
    w_t, w_grp = 0.0, zeros(Float64, ng, nq, 3)
    cf, inc, exp, pvr = 0.0, 0.0, 0.0, 0.0
    rsp, rsp_grp = zeros(Float64, nq), zeros(Float64, ng, nq)
    cf_grp, inc_grp, exp_grp = zeros(Float64, ng, nq, 3), zeros(Float64, ng, nq, 3), zeros(Float64, ng, nq, 3)

    if !haskey(nation_cf, y); nation_cf[y] = Dict{String, Float64}() end
    if !haskey(nation_inc, y); nation_inc[y] = Dict{String, Float64}() end
    if !haskey(nation_exp, y); nation_exp[y] = Dict{String, Float64}() end
    if !haskey(nation_pov, y); nation_pov[y] = Dict{String, Float64}() end
    if !haskey(nation_rsp, y); nation_rsp[y] = Dict{String, Array{Float64, 1}}() end
    if !haskey(nation_rsp_gr, y); nation_rsp_gr[y] = Dict{String, Array{Float64, 2}}() end
    if !haskey(nation_cf_gr, y); nation_cf_gr[y] = Dict{String, Array{Float64, 3}}() end
    if !haskey(nation_inc_gr, y); nation_inc_gr[y] = Dict{String, Array{Float64, 3}}() end
    if !haskey(nation_exp_gr, y); nation_exp_gr[y] = Dict{String, Array{Float64, 3}}() end

    ci = findfirst(x -> lowercase(x) == lowercase(category), cat_list)
    p_l = grt[1]

    for hi = 1:nh
        h = hhl[hi]
        hh = hhs[h]
        sz = hh.size
        h_inc = (income_mode ? hh.totinc : hh.totexp)
        gi = (h_inc >= grt[end] ? ng : findfirst(x -> x > h_inc, grt))
        w = (weight_mode ? hh.popwgh : 1.0)
        w_sz = w * sz
        w_cf = w * hcf[hi, ci]
        w_inc = w * hh.totinc
        w_exp = w * hh.totexp

        w_t += w_sz
        cf += w_cf
        inc += w_inc
        exp += w_exp
        pvr += (h_inc < p_l ? w_sz : 0.0)

        for qi = 1:nq
            h_rsp = hhs_rsp[hi, qi]
            rsp_idx = [1, (h_rsp ? 2 : 3)]
            w_grp[gi, qi, rsp_idx] .+= w_sz
            rsp[qi] += (h_rsp ? w_sz : 0.0)
            rsp_grp[gi, qi] += (h_rsp ? w_sz : 0.0)
            cf_grp[gi, qi, rsp_idx] .+= w_cf
            inc_grp[gi, qi, rsp_idx] .+= w_inc
            exp_grp[gi, qi, rsp_idx] .+= w_exp
        end
    end
    cf /= w_t
    inc /= w_t
    exp /= w_t
    pvr /= w_t
    rsp /= w_t
    rsp_grp ./= w_grp[:, :, 1]
    cf_grp ./= w_grp
    inc_grp ./= w_grp
    exp_grp ./= w_grp

    nation_cf[y][n], nation_inc[y][n], nation_exp[y][n], nation_pov[y][n] = cf, inc, exp, pvr
    nation_rsp[y][n], nation_rsp_gr[y][n] = rsp, rsp_grp
    nation_cf_gr[y][n], nation_inc_gr[y][n], nation_exp_gr[y][n] = cf_grp, inc_grp, exp_grp
end

function readResponseData(year, nation, inputFile; qst_label =[], qst_res = [], hhid_label = "", weight_mode = true)

    global hh_list, qt_list, households, responses

    y, n = year, nation
    hhl, hhs = hh_list[y][n], households[y][n]
    nh, nq = length(hhl), length(qst_label)

    if !haskey(qt_list, y); qt_list[y] = Dict{String, Array{String, 1}}() end
    if !haskey(responses, y); responses[y] = Dict{String, Array{Bool, 2}}() end
    qt_list[y][n] = qst_label[:]

    qrs = [[] for i = 1:nq]
    for i = 1:nq
        qr = qst_res[i]
        qrs[i] = (isa(qr, Array{}) ? qr : [qr])
    end

    hhs_rsp = zeros(Bool, nh, nq)

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = string.(strip.(split(readline(f), f_sep)))
    hti = findfirst(x -> x == hhid_label, title)
    qti = [findfirst(x -> x == ql, title) for ql in qst_label]
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        h = s[hti]
        hi = findfirst(x -> x == h, hhl)
        for qi = 1:nq; if s[qti[qi]] in qrs[qi]; hhs_rsp[hi, qi] = true end end
    end
    close(f)

    responses[y][n] = hhs_rsp
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

function printRegionStatus(year, nation, outputFile)

    global reg_list, region_pov, region_cf, region_inc, region_exp

    y, n = year, nation
    rl, r_pov, r_cf, r_inc, r_exp = reg_list[y][n], region_pov[y][n], region_cf[y][n], region_inc[y][n], region_exp[y][n]
    nr = length(rl)

    f_sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")
    print(f, "Region")
    print(f, f_sep, "CF per capita", f_sep, "Income per capita", f_sep, "Expenditure per capita", f_sep, "Poverty rate")
    println(f)
    for i = 1:nr
        print(f, rl[i])
        print(f, f_sep, r_cf[i], f_sep, r_inc[i], f_sep, r_exp[i], f_sep, r_pov[i])
        println(f)
    end

    close(f)
end

function printRegionalGroupStatus(year, nation, outputFile)

    global gr_list, gr_thrs, reg_list, region_rsp_gr, region_cf_gr, region_inc_gr, region_exp_gr
    y, n = year, nation
    gl, gtl, rl = gr_list[n], gr_thrs[y][n], reg_list[y][n]
    cf_grp, inc_grp, exp_grp = region_cf_gr[y][n], region_inc_gr[y][n], region_exp_gr[y][n]
    ng, nr = length(gl), length(rl)

    gr_label = ["" for i = 1:ng]
    for i = 1:ng-1; gr_label[i] = gl[i] * " (less " * string(gtl[i]) * ")" end
    gr_label[end] = gl[end] * " (over or same" * string(gtl[end]) * ")"

    st_label = ["_cf", "_inc", "_exp"]

    f_sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")
    print(f, "Region")
    for i = 1:ng; print(f, f_sep, gr_label[i]) end
    println(f)
    for stt in zip(st_label, [cf_grp, inc_grp, exp_grp])
        st_lab, stt_grp = stt
        for i = 1:nr
            print(f, rl[i] * st_lab)
            for j = 1:ng; print(f, f_sep, stt_grp[i, j, 1, 1]) end
            println(f)
        end
    end

    close(f)
end

function printRegionalGroupStatusByResponse(year, nation, outputFile)

    global gr_list, gr_thrs, reg_list, region_rsp, region_rsp_gr, region_cf_gr, region_inc_gr, region_exp_gr
    y, n = year, nation
    gl, gtl, ql, rl = gr_list[n], gr_thrs[y][n], qt_list[y][n], reg_list[y][n]
    rsp_rg, rsp_grp = region_rsp[y][n], region_rsp_gr[y][n]
    cf_grp, inc_grp, exp_grp = region_cf_gr[y][n], region_inc_gr[y][n], region_exp_gr[y][n]
    ng, nr, nq = length(gl), length(rl), length(ql)

    gr_label = ["" for i = 1:ng]
    for i=1:ng-1; gr_label[i] = "_" * gl[i] * " (less " * string(gtl[i]) * ")" end
    gr_label[end] = "_" * gl[end] * " (over or same" * string(gtl[end]) * ")"

    st_label = ["_cf", "_inc", "_exp"]
    rs_label = ["_all", "_yes", "_no"]

    f_sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")
    print(f, "Region")
    for i = 1:nq; print(f, f_sep, ql[i]) end
    println(f)
    for i = 1:nr
        print(f, rl[i] * "_res")
        for j = 1:nq; print(f, f_sep, rsp_rg[i, j]) end
        println(f)
    end
    for i = 1:nr, j = 1:ng
        print(f, rl[i] * gr_label[j] * "_res")
        for k = 1:nq; print(f, f_sep, rsp_grp[i, j, k]) end
        println(f)
    end
    for stt in zip(st_label, [cf_grp, inc_grp, exp_grp])
        st_lab, stt_grp = stt
        for i = 1:nr, j = 1:ng, k = 1:3
            print(f, rl[i] * gr_label[j] * st_lab * rs_label[k])
            for l = 1:nq; print(f, f_sep, stt_grp[i, j, l, k]) end
            println(f)
        end
    end
    close(f)
end

function printNationalGroupStatus(year, nation, outputFile)

    global gr_list, gr_thrs, reg_list
    global nation_pov, nation_cf, nation_inc, nation_exp
    global nation_rsp_gr, nation_cf_gr, nation_inc_gr, nation_exp_gr
    y, n = year, nation
    gl, gtl = gr_list[n], gr_thrs[y][n]
    n_pov, n_cf, n_inc, n_exp = nation_pov[y][n], nation_cf[y][n], nation_inc[y][n], nation_exp[y][n]
    cf_grp, inc_grp, exp_grp = nation_cf_gr[y][n], nation_inc_gr[y][n], nation_exp_gr[y][n]
    ng = length(gl)

    gr_label = ["" for i = 1:ng]
    for i = 1:ng-1; gr_label[i] = gl[i] * " (less " * string(gtl[i]) * ")" end
    gr_label[end] = gl[end] * " (over or same" * string(gtl[end]) * ")"

    st_label = ["_cf", "_inc", "_exp"]

    f_sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")
    print(f, "Status")
    print(f, f_sep, nation)
    for i = 1:ng; print(f, f_sep, gr_label[i]) end
    println(f)
    print(f, nation * "_pov", f_sep, n_pov)
    for i = 1:ng; print(f, f_sep) end
    println(f)
    for stt in zip(st_label, [n_cf, n_inc, n_exp], [cf_grp, inc_grp, exp_grp])
        st_lab, stt_nat, stt_grp = stt
        print(f, nation * st_lab)
        print(f, f_sep, stt_nat)
        for i = 1:ng; print(f, f_sep, stt_grp[i, 1, 1]) end
        println(f)
    end
    close(f)
end

function printNationalGroupStatusByResponse(year, nation, outputFile)

    global gr_list, gr_thrs, nation_rsp, nation_rsp_gr, nation_cf_gr, nation_inc_gr, nation_exp_gr
    y, n = year, nation
    gl, gtl, ql = gr_list[n], gr_thrs[y][n], qt_list[y][n]
    rsp, rsp_grp = nation_rsp[y][n], nation_rsp_gr[y][n]
    cf_grp, inc_grp, exp_grp = nation_cf_gr[y][n], nation_inc_gr[y][n], nation_exp_gr[y][n]
    ng, nq = length(gl), length(ql)

    gr_label = ["" for i = 1:ng]
    for i=1:ng-1; gr_label[i] = "_" * gl[i] * " (less " * string(gtl[i]) * ")" end
    gr_label[end] = "_" * gl[end] * " (over or same" * string(gtl[end]) * ")"

    st_label = ["_cf", "_inc", "_exp"]
    rs_label = ["_all", "_yes", "_no"]

    f_sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")
    print(f, "Status")
    for i = 1:nq; print(f, f_sep, ql[i]) end
    println(f)
    print(f, nation * "_res")
    for i = 1:nq; print(f, f_sep, rsp[i]) end
    println(f)
    for i = 1:ng
        print(f, nation * gr_label[i] * "_res")
        for j = 1:nq; print(f, f_sep, rsp_grp[i, j]) end
        println(f)
    end
    for stt in zip(st_label, [cf_grp, inc_grp, exp_grp])
        st_lab, stt_grp = stt
        for i = 1:ng, j = 1:3
            print(f, nation * gr_label[i] * st_lab * rs_label[j])
            for k = 1:nq; print(f, f_sep, stt_grp[i, k, j]) end
            println(f)
        end
    end
    close(f)
end

function printRegionalGroupCompared(year, nation, outputFile)

    global gr_list, gr_thrs, reg_list, region_cf_gr
    y, n = year, nation
    gl, gtl, ql, rl = gr_list[n], gr_thrs[y][n], qt_list[y][n], reg_list[y][n]
    cf_grp = region_cf_gr[y][n]
    ng, nr, nq = length(gl), length(rl), length(ql)

    gr_label = ["" for i = 1:ng]
    for i=1:ng-1; gr_label[i] = "_" * gl[i] * " (less " * string(gtl[i]) * ")" end
    gr_label[end] = "_" * gl[end] * " (over or same" * string(gtl[end]) * ")"

    f_sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")
    print(f, "Region")
    for i = 1:nq; print(f, f_sep, ql[i]) end
    println(f)
    for i = 1:nr, j = 1:ng
        print(f, rl[i] * gr_label[j] * "_cf")
        for k = 1:nq; print(f, f_sep, cf_grp[i, j, k, 2] < cf_grp[i, j, k, 3]) end
        println(f)
    end
    close(f)
end

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

end
