# Developed date: 17. Oct. 2022
# Last modified date: 27. Oct. 2022
# Subject: Group statistics estimator
# Description: Calculate statistics of each household group
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

module GroupEstimator

include("MicroDataReader.jl")

using Statistics
using .MicroDataReader
mdr = MicroDataReader

mutable struct group
    code::String    # group code
    type::String    # group classification type: "s" (String) or "v" (Value)
    sort::String    # sort of group (available only if type == "s")
    range::Tuple{Float64, Floate64}     # range of the group (available only if type == "v")
    unit::String    # unit of the range (available only if type == "v")

    group(c, t, s="", r=(0.0,0.0)) = new(c, t, s, r)
end

yr_list = Array{Int, 1}()       # year list: {YYYY}
nat_list = Array{String, 1}()   # nation list: {A3}
cat_list = Array{String, 1}()   # category list
reg_list = Dict{Int, Dict{String, Array{String, 1}}}()    # region list: {year, {nation A3, {region}}}

gr_list = Dict{String, Dict{String, Array{group, 1}}}()   # group code list: {nation A3, {group, {code}}}
qt_list = Dict{Int, Dict{String, Array{String, 1}}}()  # questionary code list: {year, {nation A3, {code}}}

hh_list = Dict{Int, Dict{String, Array{String, 1}}}() # Household ID: {year, {nation A3, {hhid}}}
households = Dict{Int, Dict{String, Dict{String, mdr.household}}}() # household dict: {year, {nation A3, {hhid, household}}}

region_pov = Dict{Int, Dict{String, Dict{String, Float64}}}()   # region's poverty rate: {year, {nation A3, {region, rate}}}
region_stt = Dict{Int, Dict{String, Array{Float64, 2}}}()   # region's response rates: {year, {nation A3, {region, questionary code}}}

function importMicroData(mdata::Module)

    global yr_list, nat_list
    global hh_list, reg_list, households = mdata.hh_list, mdata.dist_list, mdata.households

    yrs = sort(collect(keys(hh_list)))
    for y in yrs
        if !(y in yr_list); push!(yr_list, y) end
        nats = sort(collect(keys(hh_list[y])))
        for n in nats; if !(n in nat_list); push!(nat_list, n) end end
    end
    sort!(yr_list)
    sort!(nat_list)
end

function readResponseData(Int year, String nation, String inputFile; qst_label =[], qst_res = [], hhid_label = "")

    global hh_list, qt_list, reg_list, households, region_stt
    y, n = year, nation
    hhs = households[y][n]


    if !haskkey(region_stt, y); region_stt[y] = Dict{String, Dict{String, Dict{String, Float64}}}() end
    if !haskkey(region_stt[y], n); region_stt[y][n] = Dict{String, Dict{String, Float64}}() end

    if !haskkey(qt_list, y); qt_list[y] = Dict{String, Array{String, 1}}() end
    qt_list[y][n] = qst_label[:]

    nr, nq = length(reg_list[y][n]), length(qst_label)
    r_idx = Dict([reg_list[y][n][i] => i for i = 1:nr])

    reg_rate = zeros(Float64, nr, nq)
    reg_samp = zeros(Float64, nr)

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = string.(strip.(split(readline(f), f_sep)))
    hi = findfirst(x->x==hhid_label, title)
    qi = [findfirst(x->x==ql, title) for ql in qst_label]
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        hh = hhs[s[hi]]
        ri = r_idx[hh.district]
        for i = 1:nq; if s[qi[i]] == qst_res[i]; reg_rate[ri, i] += 1.0 end end
        reg_samp[ri] += hh.size
    end
    reg_rate ./= reg_samp
    close(f)

    region_stt[y][n] = reg_rate
end

function printResponseStatistics(Int year, String nation, String outputFile)

    global reg_list, qt_list, region_stt
    y, n = year, nation
    rgl, qtl, reg_rate = reg_list[y][n], qt_list[y][n], region_stt[y][n]
    rn, qn = lenth(rgl), length(qtl)

    f_sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")
    for qt in qtl; print(f, f_sep, qt) end
    for i = 1:nr
        print(f, rgl[i])
        for j = 1:qn; print(f, f_sep, reg_rate[i,j]) end
        println(f)
    end
    println(f)

    close(f)
end

function readQuestionaryCodes(String inputFile; code_label = "code")

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

function readGroupData(String inputFile)
    # "Type" = "s" (String) or "v" (Value)

    global gr_list
    sectors = ["Group", "Code", "Type", "Range", "Unit"]

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = string.(strip.(split(readline(f), f_sep)))
    gi, ci, ti, ri, ui = [findfirst(x->x==label, title) for label in sectors]
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        gr, tp = s[gi], lowercase(s[ti])
        if !haskey(gr_list, gr); gr_list[gr] = Array{group, 1}() end
        push!(gr_list[gr], group(s[ci], tp, (tp == "s" ? s[ri] : ""), (tp == "v" ?  : (0.0,0.0)))

    end
    close(f)

end

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end
