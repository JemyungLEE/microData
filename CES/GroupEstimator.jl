# Developed date: 17. Oct. 2022
# Last modified date: 19. Oct. 2022
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

gr_list = Dict{String, Array{group, 1}}()  # group code list: {group, {code}}
qt_list = Array{String, 1}()    # questionary code list

hh_list = Dict{Int, Dict{String, Array{Tuple{String,String}, 1}}}()           # Household ID: {year, {nation, {hhid}}}
households = Dict{Int, Dict{String, Dict{String, mdr.household}}}() # household dict: {year, {nation A3, {hhid, household}}}

function importMicroData(mdata::Module)

    global yr_list, nat_list
    global hh_list = mdata.hh_list
    global households = mdata.households

    yrs = sort(collect(keys(hh_list)))
    for y in yrs
        if !(y in yr_list); push!(yr_list, y) end
        nats = sort(collect(keys(hh_list[y])))
        for n in nats; if !(n in nat_list); push!(nat_list, n) end end
    end
    sort!(yr_list)
    sort!(nat_list)
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
    # "Type" = "String" or "Value"

    global gr_list

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)
    title = string.(strip.(split(readline(f), f_sep)))
    gi, ci, ti, ri, ui = [findfirst(x->x==label, title) for label in ["Group", "Code", "Type", "Range", "Unit"]]
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
