# SPDX-FileCopyrightText: © 2021 Jemyung Lee <jemyung81@gmail.com>
# SPDX-License-Identifier: GPL-3.0

module EmissionEstimator

# Developed date: 26. Apr. 2021
# Last modified date: 12. Sep. 2023
# Subject: Calculate household carbon emissions
# Description: Calculate direct and indirect carbon emissions by analyzing
#              Customer Expenditure Survey (CES) or Household Budget Survey (HBS) micro-data.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX
using LinearAlgebra
using SparseArrays
using Statistics

mutable struct tables
    year::Int
    t::Array{Float64, 2}    # T matrix, MRIO tables
    v::Array{Float64, 2}    # V matrix, value added
    y::Array{Float64, 2}    # Y matrix, final demand
    q::Array{Float64, 2}    # Q matrix, indicators (satellite accounts)

    tables(year, nt, nv, ny, nq) =  new(year, zeros(nt, nt), zeros(nv, nt), zeros(nt, ny), zeros(nq, nt))
end

mutable struct idx      # index structure
    nation::String      # A3 abbreviation
    entity::String
    sector::String

    idx(nat::String, ent::String, sec::String) = new(nat, ent, sec)
end

mutable struct ind      # indicator structure
    code::String
    name::String
    item::String

    ind(cd::String, na::String, it::String) = new(cd, na, it)
end

global natA3 = Dict{String, String}()      # Nation name's A3 abbreviation, {Nation, A3}
global natList = Array{String, 1}()        # Nation A3 list
global natName = Dict{String, String}()    # Nation names: {A3, Nation}

global hh_list = Dict{Int, Dict{String, Array{String, 1}}}()   # Household ID list: {year, {nation A3, {hhid}}}
global sc_list = Dict{Int, Dict{String, Array{String, 1}}}()   # commodity code list: {year, {nation A3, {code}}}
global sc_unit = Dict{Int, Dict{String, Array{String, 1}}}()   # commodity units: {year, {nation A3, {unit: kg/week, liter/month, etc.}}}
global hhExp = Dict{Int, Dict{String, Array{Float64, 2}}}()    # Household expenditure matrix: {year, {nation, {commodity, hhid}}}
global hhCmm = Dict{Int, Dict{String, Array{Float64, 2}}}()    # Household physical consumption matrix: {year, {nation, {commodity, hhid}}}

# indirect carbon emission variables
global ti = Array{idx, 1}()                # index T
global vi = Array{idx, 1}()                # index V
global yi = Array{idx, 1}()                # index Y
global qi = Array{ind, 1}()                # index Q
global lti = []                            # inversed Leontief matrix
global eoraExp = Dict{Int, Dict{String, Array{Float64, 2}}}()   # transformed households expenditure: {year, {nation, {Eora sectors, households}}}
global mTables = Dict{Int, tables}()       # {Year, tables}
global concMat = Dict{Int, Dict{String, Array{Float64, 2}}}()   # Concordance matrix: {year, {nation, {Eora sectors, CES/HBS sectors}}}
global concMatWgh = Dict{Int, Dict{String, Array{Float64, 2}}}()# Weighted concordance matrix: {year, {nation, {Eora sectors, CES/HBS sectors}}}
global indirectCE = Dict{Int, Dict{String, Array{Float64, 2}}}()# indirect carbon emission: {year, {nation, {CES/HBS sector, households}}}

# direct carbon emission variables
global concMatDe = Dict{Int, Dict{String, Array{Float64, 2}}}()        # Concordance matrix for direct emission: {year, {nation, {DE sectors, {CES/HBS sectors}}}
global de_sectors = Dict{Int, Dict{String, Array{String, 1}}}() # direct emission sectors: {year, {nation, {DE sector}}}
global de_units = Dict{Int, Dict{String, Array{String, 1}}}()   # units of converting rates: {year, {nation, {DE category, DE unit}}}
global directCE = Dict{Int, Dict{String, Array{Float64, 2}}}()  # direct carbon emission: {year, {nation, {}}}

# direct carbon emission variables
global de_list = Array{String, 1}()                             # DE IEA sector label list
global de_pr_link = Dict{Int, Dict{String, Array{String, 1}}}() # DE sector - IEA price sector link: {year, {price item, {De item}}}
global de_emits = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA emission: {year, {nation, {emission}}}
global de_enbal = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA energy balance: {year, {nation, {energy}}}
global de_massc = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA mass to volume converting rate: {year, {nation, {barrels/tonne}}}
global de_enerc = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA mass to energy converting rate: {year, {nation, {KJ/Kg}}}

global de_price = Dict{Int, Dict{String, Array{Float64, 1}}}()  # fuel prices: {year, {nation, {prices}}}
global de_price_unit = Dict{Int, Dict{String, Array{String, 1}}}()  # fuel price unit: {year, {nation, {units}}}
global de_price_item = Dict{Int, Array{String, 1}}()            # IEA price items: {year, {items}}
global de_intens = Dict{Int, Dict{String, Array{Float64, 1}}}() # IEA price items' CO2 intensity: {year, {nation, {tCO2/unit(USD)}}}
global de_energy = Dict{Int, Dict{String, Array{Float64, 1}}}() # IEA price items' CO2 energy balance: {year, {nation, {TJ}}}

global sc_conv_qnt = Dict{Int, Dict{String, Array{Float64, 1}}}()   # CES/HBS sectors' converting rate from physical untit to Kg
global de_intens_qnt = Dict{Int, Dict{String, Array{Float64, 1}}}() # IEA sectors' CO2 intensity: {year, {nation, {tCO2/unit(Kg, litre, etc.)}}}


function readIndex(indexFilePath)

    global natA3, ti, vi, yi, qi = Dict{String, String}(), Array{idx, 1}(), Array{idx, 1}(), Array{idx, 1}(), Array{ind, 1}()

    f = open(indexFilePath*"a3.csv"); readline(f)
    for l in eachline(f); l = string.(split(replace(l,"\""=>""), ',')); natA3[l[1]] = l[2] end
    close(f)
    f = open(indexFilePath*"index_t.csv"); readline(f)
    for l in eachline(f); l=string.(strip.(split(replace(l,"\""=>""), ',', limit = 5))); push!(ti, idx(l[3],l[4],l[5])) end
    close(f)
    f = open(indexFilePath*"index_v.csv"); readline(f)
    for l in eachline(f); l=string.(strip.(split(replace(l,"\""=>""), ',', limit = 5))); push!(vi, idx(l[3],l[4],l[5])) end
    close(f)
    f = open(indexFilePath*"index_y.csv"); readline(f)
    for l in eachline(f); l=string.(strip.(split(replace(l,"\""=>""), ',', limit = 5))); push!(yi, idx(l[3],l[4],l[5])) end
    close(f)
    f = open(indexFilePath*"index_q.csv"); readline(f)
    for l in eachline(f); l=string.(strip.(split(replace(l,"\""=>""), ',', limit = 5))); push!(qi, ind(l[2],l[3],l[4])) end
    close(f)
end

function readIOTables(year, tfile, vfile, yfile, qfile)

    global mTables

    tb = tables(year, length(ti), length(vi), length(yi), length(qi))

    f = open(tfile)
    i = 1; for l in eachline(f); tb.t[i,:] = [parse(Float64, x) for x in split(l, ',')]; i += 1 end
    close(f)
    f = open(vfile)
    i = 1; for l in eachline(f); tb.v[i,:] = [parse(Float64, x) for x in split(l, ',')]; i += 1 end
    close(f)
    f = open(yfile)
    i = 1; for l in eachline(f); tb.y[i,:] = [parse(Float64, x) for x in split(l, ',')]; i += 1 end
    close(f)
    f = open(qfile)
    i = 1; nt = length(ti); for l in eachline(f); tb.q[i,:] = [parse(Float64, x) for x in split(l, ',')][1:nt]; i += 1 end
    close(f)

    mTables[year] = tb
end

function rearrangeMRIOtables(year; qmode = "PRIMAP")
    # Note: Q-table's indexing should be re-arranged corresponding to the changed Eora tables

    global natList

    if qmode == "I_CHG_CO2"
        ql = deleteat!(collect(10:64), [12,13,43,44,45,46,47])  # I-GHG-CO2 emissions (Gg)
    elseif qmode == "PRIMAP"
        ql = filter(x -> startswith(qi[x].item, "PRIMAP|KYOTOGHGAR4|TOTAL") && endswith(qi[x].item, "GgCO2eq"), 1:length(qi))
    end

    global ti = ti[filter(x -> uppercase(ti[x].nation) != "ROW", 1:length(ti))]
    global vi = vi[filter(x -> uppercase(vi[x].nation) != "ROW", 1:length(vi))]
    global yi = yi[filter(x -> uppercase(yi[x].nation) != "ROW", 1:length(yi))]
    global qi = qi[ql]

    for t in ti; if !(t.nation in natList); push!(natList, t.nation) end end

    tb = mTables[year]

    nt = length(ti)
    tb.t = tb.t[1:nt, 1:nt]
    tb.v = tb.v[1:length(vi), 1:nt]
    tb.y = tb.y[1:nt, 1:length(yi)]
    tb.q = tb.q[ql, 1:nt]
end

# function rearrangeIndex(; qmode = "")
#
#     if qmode == "I_CHG_CO2"; ql = deleteat!(collect(10:64), [12,13,43,44,45,46,47]) # I-GHG-CO2 emissions (Gg)
#     elseif qmode == "PRIMAP"; ql = [2570]                                           # PRIMAP|KYOTOGHGAR4|TOTAL|GgCO2eq
#     else println("Define Q_table mode.")
#     end
#
#     global ti = ti[1:end-1]
#     global vi = vi[1:end-6]
#     global yi = yi[1:end-6]
#     global qi = qi[ql]
#
#     global natList
#     for t in ti; if !(t.nation in natList); push!(natList, t.nation) end end
# end
#
# function rearrangeTables(year; qmode = "")
#     global ti, vi, yi, qi, mTables
#     if qmode == "I_CHG_CO2"; ql = deleteat!(collect(10:64), [12,13,43,44,45,46,47]) # I-GHG-CO2 emissions (Gg)
#     elseif qmode == "PRIMAP"; ql = [2570]                                           # PRIMAP|KYOTOGHGAR4|TOTAL|GgCO2eq
#     else println("Define Q_table mode.")
#     end
#
#     tb = mTables[year]
#
#     nt = length(ti)
#     tb.t = tb.t[1:nt, 1:nt]
#     tb.v = tb.v[1:length(vi), 1:nt]
#     tb.y = tb.y[1:nt, 1:length(yi)]
#     tb.q = tb.q[ql, 1:nt]
# end

function getDomesticData(year, nation, hhid_list, sector_list, expMat, qntMat; cmmUnit = [])

    global hh_list, sc_list, sc_unit, hhExp, hhCmm

    if !haskey(hh_list, year); hh_list[year] = Dict{String, Array{String, 1}}() end
    if !haskey(sc_list, year); sc_list[year] = Dict{String, Array{String, 1}}() end
    if !haskey(sc_unit, year); sc_unit[year] = Dict{String, Array{String, 1}}() end
    if !haskey(hhExp, year); hhExp[year] = Dict{String, Array{Float64, 2}}() end
    if !haskey(hhCmm, year); hhCmm[year] = Dict{String, Array{Float64, 2}}() end

    hl = hh_list[year][nation] = hhid_list
    sl = sc_list[year][nation] = sector_list
    if size(expMat,1) == length(sl) && size(expMat,2) == length(hl); global hhExp[year][nation] = expMat
    elseif size(expMat,2) == length(sl) && size(expMat,1) == length(hl); global hhExp[year][nation] = transpose(expMat)
    else println(", Matrices sizes don't match: expMat ", size(expMat), ", hhid ", size(hl), ", sec", size(sl))
    end
    if length(qntMat) > 0
        if size(qntMat,1) == length(sl) && size(qntMat,2) == length(hl); global hhCmm[year][nation] = qntMat
        elseif size(qntMat,2) == length(sl) && size(qntMat,1) == length(hl); global hhCmm[year][nation] = transpose(qntMat)
        else println(", Matrices sizes don't match: qntMat ", size(qntMat), ", hhid ", size(hl), ", sec", size(sl))
        end
    end

    if length(cmmUnit) > 0; sc_unit[year][nation] = cmmUnit end
end

function readEmissionIntensity(year, nation, sectorFile, intensityFile; quantity = false, emit_unit = "tCO2", curr_unit = "USD")

    global de_sectors, de_intens, de_units, de_energy, de_intens_qnt
    if !haskey(de_sectors, year); de_sectors[year] = Dict{String, Array{String, 1}}() end
    if !haskey(de_intens, year); de_intens[year] = Dict{String, Array{Float64, 1}}() end
    if !haskey(de_intens_qnt, year); de_intens_qnt[year] = Dict{String, Array{Float64, 1}}() end
    if !haskey(de_energy, year); de_energy[year] = Dict{String, Array{Float64, 1}}() end
    if !haskey(de_units, year); de_units[year] = Dict{String, Array{String, 1}}() end
    if !haskey(de_sectors[year], nation); de_sectors[year][nation] = Array{String, 1}() end
    if !haskey(de_intens[year], nation); de_intens[year][nation] = Array{Float64, 1}() end
    if !haskey(de_intens_qnt[year], nation); de_intens_qnt[year][nation] = Array{Float64, 1}() end
    if !haskey(de_energy[year], nation); de_energy[year][nation] = Array{Float64, 1}() end
    if !haskey(de_units[year], nation); de_units[year][nation] = Array{String, 1}() end

    dsl = de_sectors[year][nation]
    # qnt_units = ["liter", "kg", "m^3"]
    qnt_units = ["kg"]

    f_sep = getValueSeparator(sectorFile)
    f = open(sectorFile)
    ci = findfirst(x->x=="Sector", string.(strip.(split(readline(f), f_sep))))
    for l in eachline(f)
        c = string.(strip.(split(l, f_sep)))[ci]
        if !(c in dsl); push!(dsl, c) end
    end
    close(f)

    nds = length(dsl)
    dit, den, dun = zeros(Float64, nds), zeros(Float64, nds), Array{String, 1}(undef, nds)
    f_sep = getValueSeparator(intensityFile)
    essential = ["Year", "Nation", "DE_code", "DE_sector", "Factor", "F_unit", "Emission", "E_unit"]
    f = open(intensityFile)
    title = string.(strip.(split(readline(f), f_sep)))
    i = [findfirst(x->x==et, title) for et in essential]
    ener_unit = ""
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        if s[i[1]] == string(year) && s[i[2]] == nation
            unit = string.(strip.(split(s[i[6]], '/')))
            if emit_unit == unit[1] && ((!quantity && curr_unit == unit[2]) || (quantity && lowercase(unit[2]) in qnt_units))
                di = findfirst(x->x==s[i[4]], dsl)
                dit[di], dun[di] = parse(Float64, s[i[5]]), s[i[6]]
                if ener_unit == ""; ener_unit = s[i[8]] end
                if ener_unit == s[i[8]]; den[di] = parse(Float64, s[i[7]]) end
            end
        end
    end
    de_energy[year][nation] = den
    if !quantity; de_intens[year][nation], de_units[year][nation] = dit, dun
    else de_intens_qnt[year][nation], de_units[year][nation] = dit, dun
    end
    close(f)
end

function setNationDict(nat_dict)
    global natA3, natName
    if isa(nat_dict, Dict); natName = nat_dict
    elseif length(natA3) > 0; natName = Dict(value => key for (key, value) in natA3)
    else println("nat_dict is not dictionary form.")
    end
end

function readDirectEmissionData(year, nation, filepath; output_path = "", output_tag = "", integrate = false, cpi_scaling = false, cpi_base = 0, cpi_vals = [])

    global natName, de_list, de_pr_link, de_emits, de_enbal, de_massc, de_enerc, de_price, de_price_unit

    if isa(nation, String); nat_code = [nation] end
    nat_name = [natName[x] for x in nat_code]
    nn = length(nat_code)
    cpi_scaling = cpi_scaling && (year != cpi_base)

    emit_year = (year < 2018 ? (year >= 1990 ? year : 1990) : 2017) # NOTE: This part should be removed when the IEA 2018 datasets are prepared

    sector_file = filepath * "Emission_sectors.txt"
    emiss_road_file = filepath * "emission/road/" * "Emission_road_ktCO2_" * string(emit_year) * ".xlsx"
    emiss_res_file = filepath * "emission/residential/" * "Emission_residential_ktCO2_" * string(emit_year) * ".xlsx"
    enbal_road_file = filepath * "energy_balance/road/" * "EnergyBalance_road_TJ_" * string(emit_year) * ".xlsx"
    enbal_res_file = filepath  * "energy_balance/residential/"* "EnergyBalance_residential_TJ_" * string(emit_year) * ".xlsx"
    mass_conv_file = filepath * "barrels_per_tonne/" * "Barrels_per_tonne_" * string(emit_year) * ".xlsx"
    ener_conv_file = filepath * "kj_per_kg/" * "KJ_per_Kg_" * string(emit_year) * ".xlsx"
    price_trans_file = filepath * "energy_price/"* "EnergyPrice_transport.xlsx"
    price_other_file = filepath * "energy_price/" * "EnergyPrice_others.xlsx"

    emiss_file = filepath * "Emission_ktCO2_" * string(year) * ".txt"
    enbal_file = filepath * "Energy_balance_TJ_" * string(year) * ".txt"

    de_emits[year], de_enbal[year] = Dict{String, Array{Float64, 1}}(), Dict{String, Array{Float64, 1}}()
    de_massc[year], de_enerc[year] = Dict{String, Array{Float64, 1}}(), Dict{String, Array{Float64, 1}}()
    de_price[year], de_price_unit[year] = Dict{String, Array{Float64, 1}}(), Dict{String, Array{String, 1}}()
    de_price_item[year], de_pr_link[year] = Array{String, 1}(), Dict{String, Array{String, 1}}()

    f = open(sector_file)
    readline(f)
    for l in eachline(f)
        s = string.(strip.(split(l, '\t')))
        push!(de_list, s[2])
        if length(s[3]) > 0
            if !haskey(de_pr_link[year], s[3]); de_pr_link[year][s[3]] = Array{String, 1}() end
            push!(de_pr_link[year][s[3]], s[2])
        end
    end
    close(f)
    nds = length(de_list)

    for n in nat_code
        de_emits[year][n], de_enbal[year][n] = zeros(Float64, nds), zeros(Float64, nds)
        de_massc[year][n], de_enerc[year][n] = zeros(Float64, nds), zeros(Float64, nds)
    end

    de_files = [(emiss_res_file, de_emits[year]), (enbal_res_file, de_enbal[year]), (mass_conv_file, de_massc[year]), (ener_conv_file, de_enerc[year])]
    if integrate; append!(de_files, [(emiss_road_file, de_emits[year]), (enbal_road_file, de_enbal[year])]) end
    for (filename, de_data) in de_files

        sht_label = string(year)
        if year >= 2018; sht_label ="2017"  # NOTE: This part should be removed when the IEA 2018 datasets are prepared
        elseif year < 1990; sht_label ="1990"
        end

        de_idx = []
        xf = XLSX.readxlsx(filename)
        tb = xf[sht_label][:]
        de_idx = [findfirst(x->x==d, [replace(v, "Memo: "=>"") for v in tb[3,:]]) for d in de_list]
        for i = 1:nn, j = 1:nds
            ncd = nat_code[i]
            n_idx = findfirst(x->x==nat_name[i], tb[5:end,1]) + 4
            if isa(de_idx[j], Number) && isa(tb[n_idx, de_idx[j]], Number); de_data[ncd][j] += tb[n_idx, de_idx[j]] end
        end
        close(xf)
    end

    price_item = ["Gasoline", "Diesel"]
    npi = length(price_item)
    for n in nat_code; de_price[year][n], de_price_unit[year][n] = zeros(Float64, npi), ["" for i=1:npi] end
    xf = XLSX.readxlsx(price_trans_file)
    tb = xf["Transport"][:]
    yr_idx = findfirst(x->x==string(emit_year), tb[1,3:end]) + 2
    gas_lab = ["Regular motor gasoline", "Mid-grade motor gasoline", "High-grade motor gasoline"]
    die_lab = "Automotive diesel"
    unit_lab = "Total price (USD/litre)"
    for i = 1:nn
        n, n_name = nat_code[i], nat_name[i]
        n_idx = findfirst(x->!ismissing(x) && x==n_name, tb[:,1])
        for j = 1:3
            pri = tb[n_idx+1+j*5, yr_idx]
            if isa(pri, Number) && pri > 0 && strip(tb[n_idx+(j-1)*5, 2]) == gas_lab[j] && strip(tb[n_idx+1+(j-1)*5, 3]) == unit_lab
                de_price[year][n][1] = pri
                de_price_unit[year][n][1] = "USD/litre"
                break
            end
        end
        pri = tb[n_idx+16, yr_idx]
        if isa(pri, Number) && pri > 0 && strip(tb[n_idx+15, 2]) == die_lab && strip(tb[n_idx+16, 3]) == unit_lab
            de_price[year][n][2] = pri
            de_price_unit[year][n][2] = "USD/litre"
        end
    end
    close(xf)

    xf = XLSX.readxlsx(price_other_file)
    tb = xf["OtherProducts"][:]
    unit_lab = "Total price (USD/unit)"
    sorts = Array{String, 1}()
    pri_ot = Dict{Int, Dict{String, Dict{String, Tuple{Float64, String}}}}()   # {year, {nation, {fuel_sort, {price, unit}}}}
    pri_ot[year] = Dict{String, Dict{String, Tuple{Float64, String}}}()
    yr_idx = Dict(year => findfirst(x->x==string(emit_year), tb[1,3:end]) + 2)
    yrs = [year]
    if cpi_scaling
        push!(yrs, cpi_base)
        yr_idx[cpi_base] = findfirst(x->x==string(cpi_base), tb[1,3:end]) + 2
        pri_ot[cpi_base] = Dict{String, Dict{String, Tuple{Float64, String}}}()
    end
    for rt in ["Residential","Transport"], i in filter(x->!ismissing(tb[x,1]) && rsplit(tb[x,1], '.', limit=2)[2] == rt, collect(3:size(tb)[1]))
        n, s, t = split(tb[i,1], '.')
        ft, ut = strip.(rsplit(s, ('(',')'), limit=3))
        for y in yrs
            pri = tb[i+1, yr_idx[y]]
            if isa(pri, Number) && tb[i+1, 2] == unit_lab && !(haskey(pri_ot[y], n) && haskey(pri_ot[y][n], ft))
                if !(ft in sorts) && y == year; push!(sorts, ft) end
                if !haskey(pri_ot[y], n); pri_ot[y][n] = Dict{String, Tuple{Float64, String}}() end
                if occursin(" ", ut) && tryparse(Float64, split(ut, " ")[1]) != nothing
                    scl, ut = split(ut, " ")
                    pri /= parse(Float64, scl)
                end
                pri_ot[y][n][ft] = (pri, "USD/" * (ut == "litres" ? "litre" : ut))
            end
        end
    end
    close(xf)

    cpi_codes = Dict("Coal" => "CP0454", "Charcoal" => "CP0454", "Kerosene" => "CP0453")
    sort!(sorts); ns = length(sorts)

    all_avg, grp_avg = Dict(yrs .=> [zeros(Float64, ns) for i=1:length(yrs)]), Dict(yrs .=> [zeros(Float64, ns) for i=1:length(yrs)])
    all_cnt, grp_cnt = Dict(yrs .=> [zeros(Int, ns) for i=1:length(yrs)]), Dict(yrs .=> [zeros(Int, ns) for i=1:length(yrs)])
    grp_lst = Dict(yrs .=> [Dict{String, Array{String, 1}}() for i=1:length(yrs)])  # {year, {item, {nations}}}
    all_unit = ["" for i = 1:ns]
    for y in yrs, n in collect(keys(pri_ot[y])), i = 1:ns
        s = sorts[i]
        if haskey(pri_ot[y][n], s)
            if !haskey(grp_lst[y], s); grp_lst[y][s] = Array{String, 1}() end
            all_avg[y][i] += pri_ot[y][n][s][1]; all_cnt[y][i] += 1
            if n in nat_name; grp_avg[y][i] += pri_ot[y][n][s][1]; grp_cnt[y][i] += 1; push!(grp_lst[y][s], n) end
            if all_unit[i] == ""; all_unit[i] = pri_ot[y][n][s][2]
            elseif all_unit[i] != pri_ot[y][n][s][2]; println("Price units are different: ", s, "\t", all_unit[i], "\t", pri_ot[y][n][s][2])
            end
        end
    end
    for y in yrs; all_avg[y] ./= all_cnt[y] end
    for y in yrs; grp_avg[y] ./= grp_cnt[y] end

    n_ps = length(price_item)
    append!(price_item, sorts)
    for n in nat_code; append!(de_price[year][n], zeros(Float64, ns)); append!(de_price_unit[year][n], ["" for i=1:ns]) end
    for i = 1:nn, j = 1:ns
        n, n_name, s = nat_code[i], nat_name[i], sorts[j]
        price, unit = 0, ""
        if haskey(pri_ot[year], n_name) && haskey(pri_ot[year][n_name], s); price, unit = pri_ot[year][n_name][s]
        elseif grp_avg[year][j] > 0; price, unit = grp_avg[year][j], all_unit[j]
        elseif cpi_scaling && grp_avg[cpi_base][j] > 0
            g_avg = Array{Float64, 1}()
            for ng in grp_lst[cpi_base][s]
                nc = nat_code[findfirst(x -> x == ng, nat_name)]
                cpi_rate = cpi_vals[year][nc][cpi_codes[s]] / cpi_vals[cpi_base][nc][cpi_codes[s]]
                if !isa(cpi_rate, Number) && cpi_rate > 0
                    cpi_rate = cpi_vals[year][nc][cpi_codes[s][1:end-1]] / cpi_vals[cpi_base][nc][cpi_codes[s][1:end-1]]
                end
                push!(g_avg, pri_ot[cpi_base][ng][s][1] * cpi_rate)
            end
            grp_avg[year][j] = mean(g_avg)
            price, unit = grp_avg[year][j], all_unit[j]
        elseif all_avg[year][j] > 0; price, unit = all_avg[year][j], all_unit[j]
        else println(sorts[j], " does not have any price data in ", year)
        end
        de_price[year][n][j+n_ps], de_price_unit[year][n][j+n_ps] = price, unit
    end
    de_price_item[year] = price_item

    if length(output_path)>0
        mkpath(output_path)
        for (filename, de_data) in [(emiss_file, de_emits[year]), (enbal_file, de_enbal[year]), (mass_conv_file, de_massc[year]), (ener_conv_file, de_enerc[year])]
            fn = rsplit(replace(filename, filepath=>output_path), '.', limit=2)
            f = open(fn[1] * "_" * output_tag * ".txt", "w")
            for s in de_list; print(f, "\t", s) end; println(f)
            for n in nat_code; print(f, n); for i = 1:nds; print(f, "\t", de_data[n][i]) end; println(f) end
            close(f)
        end
        f = open(output_path * "Energy_price_" * string(year) * ".txt", "w")
        for pit in price_item; print(f, "\t", pit) end; println(f)
        for n in nat_code; print(f, n); for i = 1:length(price_item); print(f, "\t", de_price[year][n][i], " ", de_price_unit[year][n][i]) end; println(f) end
        close(f)
        f = open(output_path * "Energy_price_others_" * string(year) * ".txt", "w")
        for s in sorts; print(f, "\t", s) end; println(f)
        for n in sort(collect(keys(pri_ot[year])))
            print(f, n); for s in sorts; print(f, "\t", haskey(pri_ot[year][n], s) ? string(pri_ot[year][n][s][1])*" "*pri_ot[year][n][s][2] : "") end; println(f)
        end
        close(f)
    end
end

function exchangeEmCurrency(year, exchFile; target = "EUR", origin = "USD", output = "")

    global de_price, de_price_unit, de_price_item
    ni = length(de_price_item[year])

    er = Dict{Int, Float64}()
    f = open(exchFile)
    s = split(readline(f), '\t')
    er_idx = findfirst(x->x== target * "_" * origin, s)
    if er_idx != nothing; for l in eachline(f); s = split(l, '\t'); er[parse(Int,s[1])] = parse(Float64,s[er_idx]) end
    else
        er_idx = findfirst(x->x== origin * "_" * target, s)
        if er_idx != nothing; for l in eachline(f); s = split(l, '\t'); er[parse(Int,s[1])] = 1.0 / parse(Float64,s[er_idx]) end
        else println("No matching exchange rate from ", origin, " to ", target)
        end
    end
    close(f)

    if haskey(er, year); er = er[year] else println(year," year exchange rate is not on the list.") end
    nats = sort(collect(keys(de_price[year])))
    for n in nats, i = 1:ni
        curr, unit = split(de_price_unit[year][n][i], '/')
        if curr == origin
            de_price[year][n][i] *= er
            de_price_unit[year][n][i] = target * "/" * unit
        end
    end

    if length(output)>0
        f = open(output, "w")
        for pit in de_price_item[year]; print(f, "\t", pit) end; println(f)
        for n in nats; print(f, n); for i = 1:ni; print(f, "\t", de_price[year][n][i], " ", de_price_unit[year][n][i]) end; println(f) end
        close(f)
    end
end

function calculateEmissionRates(year; output = "", currency = "USD", quantity = false)

    global de_code, de_list, de_pr_link, de_emits, de_enbal, de_massc, de_enerc, de_price, de_price_unit, de_price_item
    global de_intens, de_energy, de_intens_qnt

    emi, bal, mas, ene = de_emits[year], de_enbal[year], de_massc[year], de_enerc[year]
    pri_itm, pri, pri_unt = de_price_item[year], de_price[year], de_price_unit[year]
    de_intens[year] = Dict{String, Array{Float64, 1}}()
    de_energy[year] = Dict{String, Array{Float64, 1}}()
    nats = sort(collect(keys(emi)))
    nn, nd, np = length(nats), length(de_list), length(pri_itm)
    intens = zeros(Float64, nn, np)         # tCO2 / USD
    tot_emits = zeros(Float64, nn, np)      # TJ
    emit_rate = zeros(Float64, nn, nd)      # ktCO2 / TJ
    if quantity
        de_intens_qnt[year] = Dict{String, Array{Float64, 1}}()
        intens_qnt = zeros(Float64, nn, nd) # tCO2 / kg
        ng_di  = findfirst(x -> lowercase(x) == "natural gas", de_list)
        lng_di = findfirst(x -> lowercase(x) == "natural gas liquids", de_list)
    end

    for i = 1:nn, j = 1:nd
        n = nats[i]
        if emi[n][j] > 0 && bal[n][j] > 0; emit_rate[i,j] =  emi[n][j] / bal[n][j] end
        if quantity; intens_qnt[i,j] = emit_rate[i,j] * (j != ng_di ? ene[n][j] : ene[n][lng_di]) / 10^6 end
    end

    dp_idx = [haskey(de_pr_link[year], pit) ? [findfirst(x->x==l, de_list) for l in de_pr_link[year][pit]] : [] for pit in pri_itm]
    for i = 1:nn, j in filter(x->haskey(de_pr_link[year], pri_itm[x]), collect(1:np))
        n = nats[i]
        curr, unit = split(pri_unt[n][j], '/')
        unit = lowercase(unit)
        if curr == currency
            dp_link = de_pr_link[year][pri_itm[j]]
            if unit in ["tonne", "t", "kg"]
                if length(dp_link) == 1
                    intens[i, j] = emit_rate[i, dp_idx[j][1]] * ene[n][dp_idx[j][1]] / pri[n][j] / 10^3
                    tot_emits[i, j] = bal[n][dp_idx[j][1]]
                elseif length(dp_link) > 1
                    idxs = filter(x->(emi[n][x]>0 && bal[n][x]>0 && ene[n][x]>0), dp_idx[j])
                    emit_r = [emit_rate[i, k] * ene[n][k] for k in idxs]
                    emit_t = [bal[n][k] for k in idxs]
                    if length(idxs)>0
                        tot_emits[i, j] = sum(emit_t)
                        intens[i, j] = sum(emit_r .* emit_t) / tot_emits[i, j] / pri[n][j] / 10^3
                    end
                end
                if unit in ["kg"]; intens[i, j] /= 10^3 end
            elseif unit in ["litre", "liter", "litres", "liters", "l"]
                if length(dp_link) == 1
                    intens[i, j] = emit_rate[i, dp_idx[j][1]] * ene[n][dp_idx[j][1]] / mas[n][dp_idx[j][1]] / pri[n][j] / 10^3  / 158.99
                    tot_emits[i, j] = bal[n][dp_idx[j][1]]
                elseif length(dp_link) > 1
                    idxs = filter(x->(emi[n][x]>0 && bal[n][x]>0 && ene[n][x]>0 && mas[n][x]>0), dp_idx[j])
                    emit_r = [emit_rate[i, k] * ene[n][k] / mas[n][k] for k in idxs]
                    emit_t = [bal[n][k] for k in idxs]
                    if length(idxs)>0
                        tot_emits[i, j] = sum(emit_t)
                        intens[i, j] = sum(emit_r .* emit_t) / tot_emits[i, j] / pri[n][j] / 10^3  / 158.99
                    end
                end
            elseif unit in ["mwh"]
                if length(dp_link) == 1
                    intens[i, j] = emit_rate[i, dp_idx[j][1]] / pri[n][j] * 10^3 / 277.78
                    tot_emits[i, j] = bal[n][dp_idx[j][1]]
                elseif length(dp_link) > 1
                    idxs = filter(x->(emi[n][x]>0 && bal[n][x]>0), dp_idx[j])
                    emit_r = [emit_rate[i, k] for k in idxs]
                    emit_t = [bal[n][k] for k in idxs]
                    if length(idxs)>0
                        tot_emits[i, j] = sum(emit_t)
                        intens[i, j] = sum(emit_r .* emit_t) / tot_emits[i, j] / pri[n][j] * 10^3 / 277.78
                    end
                end
            else println("Unit ", unit, " is not in the list.", )
            end
        else println("Currency ", curr, " is not: ", currency)
        end
    end
    for i = 1:nn
        idx_avail = filter(x -> intens[i, x] > 0, collect(1:np))
        if length(idx_avail) > 0
            int_avg = sum(intens[i, idx_avail] .* tot_emits[i, idx_avail]) ./ sum(tot_emits[i, idx_avail])
            for j in filter(x -> intens[i, x]==0, collect(1:np)); intens[i, j] = int_avg end
        end
        de_intens[year][nats[i]] = intens[i, :]
        de_energy[year][nats[i]] = tot_emits[i, :]
        if quantity; de_intens_qnt[year][nats[i]] = intens_qnt[i, :] end
    end

    if length(output) > 0
        mkpath(rsplit(output, '/', limit = 2)[1])
        f = open(output, "w")
        pri_idx = filter(x->haskey(de_pr_link[year], pri_itm[x]), collect(1:np))
        for i in pri_idx; print(f, "\t", pri_itm[i]) end; println(f)
        for i = 1:nn; print(f, nats[i]); for j in pri_idx; print(f, "\t", intens[i, j]) end; println(f) end
        close(f)
        if quantity
            f = open(replace(output, currency * ".txt" => "Kg.txt"), "w")
            for i in nd; print(f, "\t", de_list[i]) end; println(f)
            for i = 1:nn; print(f, nats[i]); for j = 1:nd; print(f, "\t", intens_qnt[i, j]) end; println(f) end
            close(f)
        end
    end
end

function printEmissionConvRates(year, outputFile; emit_unit = "tCO2", curr_unit = "USD", qnt_unit = "kg", quantity = false)

    global de_list, de_pr_link, de_price_item, de_intens, de_energy, de_intens_qnt, de_enbal
    nats = sort(collect(keys(de_intens[year])))
    sep = getValueSeparator(outputFile)

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    de_idxs = filter(x->haskey(de_pr_link[year], de_price_item[year][x]), collect(1:length(de_price_item[year])))
    ecr_sets = [[outputFile, de_price_item[year], de_intens[year], de_energy[year], de_idxs, curr_unit]]

    if quantity
        fname, fext = rsplit(outputFile, '.', limit = 2)
        push!(ecr_sets, [fname * "_qnt." * fext, de_list, de_intens_qnt[year], de_enbal[year], collect(1:length(de_list)), qnt_unit])
    end

    for (filename, dlst, dint, den, didx, unt) in ecr_sets
        strs = Array{String, 1}()
        if isfile(filename)
            f = open(filename)
            yr_str = string(year)
            readline(f)
            for l in eachline(f)
                yr, na = string.(strip.(split(l, sep, limit = 3)))[1:2]
                if yr != yr_str || !(na in nats); push!(strs, l) end
            end
        end

        for n in nats, j in didx
            l = string(year) * "\t" * n * "\t" * string(j) * "\t" * dlst[j] * "\t"
            l *= string(dint[n][j]) * "\t" * emit_unit * "/" * unt * "\t" * string(den[n][j]) * "\tTJ"
            push!(strs, l)
        end

        tags = [string.(strip.(split(l, sep, limit = 3)))[1:2] for l in strs]
        yl, nl = map(x->x[1], tags), map(x->x[2], tags)
        y_ord = sortperm(yl)
        strs, yl, nl = strs[y_ord], yl[y_ord], nl[y_ord]
        n_ord = sortperm(nl)
        strs, yl, nl = strs[n_ord], yl[n_ord], nl[n_ord]

        f = open(filename, "w")
        println(f, "Year\tNation\tDE_code\tDE_sector\tFactor\tF_unit\tEmission\tE_unit")
        for l in strs; println(f, l) end
        close(f)
    end
end

function buildDeConcMat(year, nation, concFile; norm = false, output = "", energy_wgh = false)

    global natList, sc_list, de_sectors, de_energy, concMatDe

    dsl, scl, den = de_sectors[year][nation], sc_list[year][nation], de_energy[year][nation]
    nd, nc = length(dsl), length(scl)
    nats = collect(keys(de_sectors[year]))
    if !haskey(concMatDe, year); concMatDe[year] = Dict{String, Array{Float64, 2}}() end
    cmat = concMatDe[year][nation] = zeros(Float64, nd, nc)
    if energy_wgh; ene_avg = [mean(filter(x->x>0, [de_energy[year][n][j] for n in nats])) for j=1:nd] end

    essential = ["Nation", "DE_sector", "CES/HBS_code", "Weight"]
    f_sep = getValueSeparator(concFile)
    f = open(concFile)
    title = string.(strip.(split(readline(f), f_sep)))
    i = [findfirst(x->x==et, title) for et in essential]
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        if s[i[1]] == nation
            decode, cescode, wgh = s[i[2]], s[i[3]], parse(Float64, s[i[4]])
            di, ci = findfirst(x->x==decode, dsl), findfirst(x->x==cescode, scl)
            if di != nothing && ci != nothing
                if energy_wgh; cmat[di, ci] += wgh * (den[di]>0 ? den[di] : (ene_avg[di]>0 ? ene_avg[di] : 1.0))
                else cmat[di, ci] += wgh
                end
            else println("Incorrect DE ($decode) or CES ($cescode) code")
            end
        else println("Nation does not match: ", naiton, ", and ", s[i[1]])
        end
    end
    close(f)

    if norm
        for i = 1:nc
            ces_sum = sum(cmat[:,i])
            if ces_sum > 0; cmat[:, i] /= ces_sum end
        end
    end

    if length(output)>0
        mkpath(rsplit(output, '/', limit = 2)[1])
        f = open(output, "w")
        for sc in scl; print(f, "\t", sc) end; println(f)
        for i = 1:nd
            print(f, dsl[i])
            for j = 1:nc; print(f, "\t", cmat[i, j]) end
            println(f)
        end
        close(f)
    end

    return cmat
end

function readDeConcMat(year, nation, concMatFile; norm = false, output = "", energy_wgh = false, float_mode = false)

    global natList, sc_list, de_sectors, de_energy, concMatDe

    dsl, scl, den = de_sectors[year][nation], sc_list[year][nation], de_energy[year][nation]
    nd, nc = length(dsl), length(scl)
    nats = collect(keys(de_sectors[year]))
    if !haskey(concMatDe, year); concMatDe[year] = Dict{String, Array{Float64, 2}}() end
    cmat = concMatDe[year][nation] = zeros(Float64, nd, nc)
    if energy_wgh; ene_avg = [mean(filter(x->x>0, [de_energy[year][n][j] for n in nats])) for j=1:nd] end

    f_sep = getValueSeparator(concMatFile)
    f = open(concMatFile)
    sec = string.(strip.(split(readline(f), f_sep)[2:end]))
    si = [findfirst(x -> x == c, sec) for c in scl]

    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        if s[1] in dsl
            di = findfirst(x -> x == s[1], dsl)
            if energy_wgh
                cmat[di, :] .+= (float_mode ? parse.(Float64, s[2:end][si]) : parse.(Int, s[2:end][si])) .* (den[di]>0 ? den[di] : (ene_avg[di]>0 ? ene_avg[di] : 1.0))
            else
                cmat[di, :] .+= (float_mode ? parse.(Float64, s[2:end][si]) : parse.(Int, s[2:end][si]))
            end
        end
    end
    close(f)

    if norm
        for i = 1:nc
            ces_sum = sum(cmat[:,i])
            if ces_sum > 0; cmat[:, i] /= ces_sum end
        end
    end

    if length(output)>0
        mkpath(rsplit(output, '/', limit = 2)[1])
        f = open(output, "w")
        for sc in scl; print(f, "\t", sc) end; println(f)
        for i = 1:nd
            print(f, dsl[i])
            for j = 1:nc; print(f, "\t", cmat[i, j]) end
            println(f)
        end
        close(f)
    end

    return cmat
end

function calculateQuantityConvRate(year, nation, de_conc_file; qnt_unit = "kg")

    # NOTE: 1 Scm (standard cubic meter) = 7.350x10^-4 tonne LNG (from IEA, http://wds.iea.org/wds/pdf/gas_documentation.pdf)

    global sc_list, sc_unit, sc_conv_qnt, de_list, de_enbal, de_massc, de_enerc
    bal, mas, ene = de_enbal[year], de_massc[year], de_enerc[year]


    scl, scu = sc_list[year][nation], sc_unit[year][nation]
    if !haskey(sc_conv_qnt, year); sc_conv_qnt[year] = Dict{String, Array{Float64, 1}}() end
    ces_de_links = Dict{String, Array{String, 1}}()
    ns = length(scl)
    scc = zeros(Float64, ns)

    f_sep = getValueSeparator(de_conc_file)
    f = open(de_conc_file)
    title = string.(strip.(split(readline(f), f_sep)))
    i = [findfirst(x->x==et, title) for et in ["Nation", "DE_sector", "CES/HBS_code"]]
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        if s[i[1]] == nation
            decode, cescode = s[i[2]], s[i[3]]
            if !haskey(ces_de_links, cescode); ces_de_links[cescode] = Array{String, 1}() end
            if !(decode in ces_de_links[cescode]); push!(ces_de_links[cescode], decode) end
        else println("Nation does not match: ", naiton, ", and ", s[i[1]])
        end
    end
    close(f)

    ng_di  = findfirst(x -> lowercase(x) == "natural gas", de_list)
    lng_di = findfirst(x -> lowercase(x) == "natural gas liquids", de_list)

    for si in filter(x -> haskey(ces_de_links, scl[x]) && length(scu[x]) > 0, 1:ns)
        s = scl[si]
        qnt = length(scu[si]) > 0 ? lowercase.(string.(split(scu[si], '/')))[1] : ""

        if qnt in ["kg"]; cr = 1.0
        elseif qnt in ["tonne", "t"]; cr = 1000.0
        elseif qnt in ["litre", "liter", "litres", "liters", "l"]
            dec = ces_de_links[s]
            if length(dec) == 1
                di = findfirst(x -> x == dec[1], de_list)
                if di == ng_di && mas[nation][di] == 0; di = lng_di end
                cr = 1000.0 / mas[nation][di] / 158.99
            else
                cr, tb = 0.0, 0.0
                for di in [findfirst(x -> x == ds, de_list) for ds in dec]
                    if di == ng_di && mas[nation][di] == 0; di = lng_di end
                    if mas[nation][di] > 0
                        cr += 1000.0 / mas[nation][di] / 158.99 * bal[nation][di]
                        tb += bal[nation][di]
                    end
                end
                cr /= tb
            end
        elseif qnt in ["m^3", "m3", "cbm"]
            dec = ces_de_links[s]
            cr = 0.735  # only correct for natural gas
            if length(filter(x -> !(lowercase(x) in ["natural gas", "natural gas liquids"]), dec)) > 0
                println("CBM unit conversion does not support: ", dec)
            end
        else println("No matching qunatity unit: ", qnt, ", of sector: ", s)
        end
        scc[si] = cr
    end

    sc_conv_qnt[year][nation] = scc
end

function calculateDirectEmission(year, nation; quantity = false, sparseMat = false, enhance = false, full = false)

    global sc_list, de_sectors, hh_list, hhExp, directCE, concMatDe, de_intens, de_intens_qnt, sc_conv_qnt
    hl, sl, dsl = hh_list[year][nation], sc_list[year][nation], de_sectors[year][nation]
    if !haskey(concMatDe, year); concMatDe[year] = Dict{String, Array{Float64, 2}}() end
    cmn = concMatDe[year][nation]
    dit = (quantity ? transpose(de_intens_qnt[year][nation]) : transpose(de_intens[year][nation]))
    he = (quantity ? hhCmm[year][nation] : hhExp[year][nation])
    nh, ns, nds = length(hl), length(sl), length(dsl)
    if !haskey(directCE, year); directCE[year] = Dict{String, Array{Float64, 2}}() end

    de = zeros(Float64, ns, nh)
    if full || enhance; dit_cmn = dit * cmn end

    if quantity
        scc = sc_conv_qnt[year][nation]

        if full
            for i = 1:ns
                hes = zeros(Float64, ns, nh)
                hes[i,:] = he[i,:]
                de[i,:] = dit_cmn * scc[i] * hes
            end
        else
            for i = 1:ns, j = 1:nh, k = 1:nds; de[i,j] += dit[k] * cmn[k,i] * scc[i] * he[i,j] end
        end
    else
        if sparseMat
            dits = dropzeros(sparse(dit))
            for i = 1:ns
                hes, cmns = zeros(Float64, ns, nh), zeros(Float64, nds, ns)
                hes[i,:] = he[i,:]
                cmns[:,i] = cmn[:,i]
                hes, cmns = dropzeros(sparse(hes)), dropzeros(sparse(cmns))
                de[i,:] = dits * cmns * hes
            end
        elseif enhance; for i = 1:ns; de[i,:] = dit_cmn[i] * he[i,:] end
        elseif full
            for i = 1:ns
                hes = zeros(Float64, ns, nh)
                hes[i,:] = he[i,:]
                de[i,:] = dit_cmn * hes
            end
        else for i = 1:ns, j = 1:nh, k = 1:nds; de[i,j] += dit[k] * cmn[k,i] * he[i,j] end
        end
    end

    directCE[year][nation] = de
end

function buildWeightedConcMat(year, eoraYear, natA3; con_mat=[], con_mat_file="", normalize=false, output="", sum_ouput="")
    # concordance matrix (Eora, Nation)

    global concMat, concMatWgh, mTables, natList, sc_list, ti, yi
    sl = sc_list[year][natA3]
    tb = mTables[eoraYear]
    ns, nn, nti = length(sl), length(natList), length(ti)
    cMat = zeros(Float64, 0, 0)

    # get final demand of nation 'natA3'
    ye = tb.y[:,findfirst(x->x.nation==natA3 && x.sector=="Household final consumption P.3h", yi)]

    # check if a nation's account have 'Commodities' entities
    chk = Dict{String, Bool}()      # {A3, (true: has 'Commodities', false: only 'Industries')}
    for t in ti
        if !haskey(chk, t.nation); chk[t.nation] = false end
        if !chk[t.nation] && t.entity == "Commodities"; chk[t.nation] = true end
    end

    # count number of 'Industries' entities by nation: 0 means having only 'Industries' or only 'Commodities'.
    cnt = zeros(Int, nn)
    for t in ti; if chk[t.nation] && t.entity == "Industries"; cnt[findfirst(x->x==t.nation, natList)] += 1 end end

    if length(con_mat) > 0
        # assemble concordance matrices
        cMat = zeros(Float64, 0, ns)
        for i = 1:nn
            if cnt[i]>0; cMat = vcat(cMat, zeros(Float64, cnt[i], ns)) end
            cMat = vcat(cMat, con_mat[natList[i]])
        end
    elseif length(con_mat_file) > 0
        cMat, sumMat = zeros(Float64, nti, ns), zeros(Float64, 0, ns)

        f_sep = getValueSeparator(con_mat_file)
        f = open(con_mat_file)
        codes = string.(strip.(split(readline(f), f_sep)[3:end]))
        if issubset(sl, codes); i = [findfirst(x->x==sc, codes) for sc in sl]
        else println(inputFile, " expenditure matrix file does not contain all essential data.")
        end
        cnt_l, country = 0, ""
        for l in eachline(f)
            cnt_l += 1
            s = string.(strip.(split(l, f_sep)))

            if country != s[1]
                country = s[1]
                ni = findfirst(x -> x==s[1], natList)
                if cnt[ni]>0; cnt_l += cnt[ni] end
            end

            if (ti[cnt_l].nation, lowercase(ti[cnt_l].sector)) == (s[1], lowercase(s[2]))
                conc = parse.(Int, s[3:end])
                cMat[cnt_l, :] = conc
            else println("Eora concordance matrix index error: at row ", cnt_l, ", ", ti[cnt_l].nation, ", ",  ti[cnt_l].sector, ", ",  s[1], ", ",  s[2])
            end
        end
        close(f)

        enat, eora_nats, conc_sum = "", Array{String, 1}(), zeros(Float64, ns)
        for i = 1:nti
            n = ti[i].nation
            if enat != n
                enat = n
                push!(eora_nats, n)
                if i > 1; sumMat = vcat(sumMat, conc_sum') end
                conc_sum = zeros(Float64, ns)
            end
            conc_sum .+= cMat[i,:]
        end
        sumMat = vcat(sumMat, conc_sum')

        if length(sum_ouput) > 0
            mkpath(rsplit(sum_ouput, '/', limit = 2)[1])
            f = open(sum_ouput, "w")
            for s in sl; print(f, "\t", s) end; println(f)
            for i = 1:length(eora_nats)
                print(f, eora_nats[i])
                for j = 1:ns; print(f, "\t", sumMat[i,j]) end
                println(f)
            end
            close(f)
        end

        if normalize
            for i = 1:nti
                n = ti[i].nation
                ni = findfirst(x -> x == n, eora_nats)
                cMat[i,:] ./= sumMat[ni,:]
            end
        end
    else println("Concordance matrix building error: no given conc_mat")
    end
    if !haskey(concMat, year); concMat[year] = Dict{String, Array{Float64, 2}}() end
    concMat[year][natA3] = copy(cMat)

    # reflect Eora final demand accounts' ratios
    for j = 1:ns
        cMat[:, j] .*= ye
        tsum = sum(cMat[:,j])
        cMat[:, j] /= tsum
    end
    if !haskey(concMatWgh, year); concMatWgh[year] = Dict{String, Array{Float64, 2}}() end
    concMatWgh[year][natA3] = copy(cMat)

    # print concordance matrix
    if length(output)>0
        mkpath(rsplit(output, '/', limit = 2)[1])
        f = open(output, "w")
        print(f,"Nation,Entity,Sector");for i=1:ns; print(f,",",sl[i]) end; println(f)
        for i=1:size(concMatWgh[year][natA3],1)
            print(f,ti[i].nation,",",ti[i].entity,",",ti[i].sector)
            for j=1:size(concMatWgh[year][natA3],2); print(f,",",concMatWgh[year][natA3][i,j]) end; println(f)
        end
        close(f)
    end
end

function calculateLeontief(year)

    global mTables, ti, vi, yi, qi, lti
    nt = length(ti)
    tb = mTables[year]

    # x = sum(tb.t, dims = 1) +  sum(tb.v, dims = 1)  # calculate X
    x = transpose(sum(tb.t, dims = 2) +  sum(tb.y, dims = 2))  # calculate X

    f = sum(tb.q, dims = 1) ./ x                    # calculate EA
    lt = Matrix{Float64}(I, nt, nt)                 # calculate Leontief matrix
    for i = 1:nt; for j = 1:nt; lt[i,j] -= tb.t[i,j] / x[j] end end
    lti = inv(lt)
    for i = 1:nt; lti[i,:] *= f[i] end
end

function calculateIndirectEmission(cesYear, eoraYear, nation; sparseMat = false, enhance = false, full = false, elapChk = 0)

    global indirectCE, mTables, concMatWgh, lti
    global hh_list, sc_list, hhExp
    global ti, vi, yi, qi

    tb = mTables[eoraYear]
    hl, sl, em = hh_list[cesYear][nation], sc_list[cesYear][nation], hhExp[cesYear][nation]
    nt, nh, ns = length(ti), length(hl), length(sl)
    cmat = concMatWgh[cesYear][nation]
    if !haskey(indirectCE, cesYear); indirectCE[cesYear] = Dict{String, Array{Float64, 2}}() end

    # calculate emission, by CES/HBS micro-data sectors, by Eora T matrix sectors
    e = zeros(Float64, ns, nh)

    st = time()     # check start time
    if enhance || full; lti_conc = lti * cmat end
    for i = 1:ns
        if enhance; e[i,:] = sum(lti_conc[:,i] * transpose(em[i,:]), dims=1)
        else
            hce = zeros(Float64, ns, nh)
            hce[i,:] = em[i,:]

            if sparseMat
                hceS = dropzeros(sparse(hce))
                hce = []
                concMatS = zeros(Float64, nt, ns)
                concMatS[:,i] = cmat[:,i]
                concMatS = dropzeros(sparse(concMatS))
                conc_hce = Array(concMatS * hceS)
                concMatS = hceS = []
                ebe = lti * conc_hce
            elseif full; ebe = lti_conc * hce
            else ebe = lti * cmat * hce       # household emission by Eora sectors
            end
            e[i,:] = sum(ebe, dims=1)       # calculate total emission (=sum of Eora emissions) of each nation sector
        end

        if elapChk > 0   # check elapsed and remained time
            elap = floor(Int, time() - st)
            (eMin, eSec) = fldmod(elap, 60)
            (eHr, eMin) = fldmod(eMin, 60)
            (rMin, rSec) = fldmod(floor(Int, (elap / i) * (ns - i)), 60)
            (rHr, rMin) = fldmod(rMin, 60)

            if i%elapChk == 0
                print("\n\t$i/$ns iterations, $eHr:$eMin:$eSec elapsed, $rHr:$rMin:$rSec remained")
            end
        end
    end

    indirectCE[cesYear][nation] = e
end

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

function printEmissions(year, nation, outputFile; mode = "ie")

    global hh_list, sc_list, directCE, indirectCE
    hl, sl= hh_list[year][nation], sc_list[year][nation]
    ns, nh = length(sl), length(hl)

    if mode == "ie"; em = indirectCE[year][nation]
    elseif mode == "de"; em = directCE[year][nation]
    else println("Wrong emission print mode: $mode")
    end

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    for h in hl; print(f, "\t", h) end
    println(f)
    for i = 1:ns
        print(f, sl[i])
        for j = 1:nh; print(f, "\t", em[i,j]) end
        println(f)
    end

    close(f)
end

end
