module EmissionEstimator

# Developed date: 29. Jul. 2020
# Last modified date: 7. Aug. 2021
# Subject: Calculate EU households carbon emissions
# Description: Calculate emissions by analyzing Eurostat Household Budget Survey (HBS) micro-data.
#              Transform HH consumptions matrix to nation by nation matrix of Eora form.
#              Intermidiate results can be utilized for the global carbon footprint mapping.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX
using LinearAlgebra
using SparseArrays
using Statistics

mutable struct tables
    year::Int16
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

global abb = Dict{String, String}()    # Nation name's A3 abbreviation, {Nation, A3}
global euA2 = Dict{String, String}()   # EU nation name's A2 abbreviation, {Nation, A2}
global euA3 = Dict{String, String}()   # EU nation name's A3 abbreviation, {Nation, A3}
global natList = Array{String, 1}()    # Nation A3 list
global ti = Array{idx, 1}()     # index T
global vi = Array{idx, 1}()     # index V
global yi = Array{idx, 1}()     # index Y
global qi = Array{ind, 1}()     # index Q

global sec = Dict{Int, Array{String, 1}}()          # Household expenditure sectors: {year, {sector}}
global hhid = Dict{Int, Array{String, 1}}()         # Household ID: {year, {hhid}}
global hhExp = Dict{Int, Array{Float64, 2}}()       # Households enpenditure: {year, {Nation sectors, households}}
global concMat = Dict{Int, Array{Float64, 2}}()     # Assembled concordance matrix {Eora sectors, Nation sectors}
global concMatWgh = Dict{Int, Array{Float64, 2}}()  # Weighted concordance matrix {Eora sectors, Nation sectors}
global concMatDe = Dict{Int, Array{Float64, 2}}()   # concordance matrix sets for direct emission

global lti = []                            # Inversed Leontief matrix
global eoraExp = Array{Float64, 2}         # Transformed households expenditure, {Eora sectors, households}
global mTables = Dict{Int16, tables}()     # {Year, tables}
global emissions = Dict{Int16, Array{Float64, 2}}()    # carbon footprint

# direct carbon emission variables
global directCE = Dict{Int, Array{Float64, 2}}()        # direct carbon emission
global deSecList = Dict{Int, Array{String, 1}}()        # direct emission sectors: {year, {DE sector}}
global deHbsList = Dict{Int, Array{String, 1}}()        # direct emission HBS sectors: {year, {HBS code}}
global deHbsSec = Dict{Int, Dict{String, String}}()     # HBS sector - DE sector link: {year, {HBS code, DE sector}}
global dePerEUR = Dict{Int, Dict{String, Dict{String, Float64}}}() # converting rate from EUR to CO2: {year, {Nation, {DE category, tCO2/EUR}}}
global deIdx = Dict{Int, Dict{String, Array{Int, 1}}}() # Direct emission sector matched EU expenditure sector index: {year, {DE sector, {HBS expenditure index}}}

global de_list = Array{String, 1}()                             # DE IEA sector label list
global de_pr_link = Dict{Int, Dict{String, Array{String, 1}}}() # DE sector - IEA price sector link: {year, {price item, {De item}}}
global de_emits = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA emission: {year, {nation, {emission}}}
global de_enbal = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA energy balance: {year, {nation, {energy}}}
global de_massc = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA mass to volume converting rate: {year, {nation, {barrels/tonne}}}
global de_enerc = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA mass to energy converting rate: {year, {nation, {KJ/Kg}}}

global de_price = Dict{Int, Dict{String, Array{Float64, 1}}}()  # fuel prices: {year, {nation, {prices}}}
global de_price_unit = Dict{Int, Dict{String, Array{String, 1}}}()  # fuel price unit: {year, {nation, {units}}}
global de_price_item = Dict{Int, Array{String, 1}}()            # IEA price items: {year, {items}}
global de_intens = Dict{Int, Dict{String, Array{Float64, 1}}}() # IEA price items' CO2 intensity: {year, {nation, {tCO2/EUR}}}
global de_energy = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA price items' CO2 energy balance: {year, {nation, {TJ}}}
global de_sectors = Dict{Int, Array{String, 1}}()           # direct emission sectors: {year, {DE sector}}

function readIndexXlsx(inputFile; revised = false)

    global abb, ti, vi, yi, qi
    if !revised
        xf = XLSX.readxlsx(inputFile)
        sh = xf["A3"]
        for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; abb[r[1]] = r[2] end end
        sh = xf["index_t"]
        for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; push!(ti, idx(r[3], r[4], r[5])) end end
        sh = xf["index_v"]
        for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; push!(vi, idx(r[3], r[4], r[5])) end end
        sh = xf["index_y"]
        for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; push!(yi, idx(r[3], r[4], r[5])) end end
        sh = xf["index_q"]
        for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; push!(qi, ind(r[2], r[3], r[4])) end end
        close(xf)
    else
        f = open(inputFile*"a3.csv"); readline(f)
        for l in eachline(f); l = string.(split(replace(l,"\""=>""), ',')); abb[l[1]] = l[2] end
        close(f)
        f = open(inputFile*"index_t.csv"); readline(f)
        for l in eachline(f); l=string.(split(replace(l,"\""=>""), ',')); push!(ti, idx(l[3],l[4],l[5])) end
        close(f)
        f = open(inputFile*"index_v.csv"); readline(f)
        for l in eachline(f); l=string.(split(replace(l,"\""=>""), ',')); push!(vi, idx(l[3],l[4],l[5])) end
        close(f)
        f = open(inputFile*"index_y.csv"); readline(f)
        for l in eachline(f); l=string.(split(replace(l,"\""=>""), ',')); push!(yi, idx(l[3],l[4],l[5])) end
        close(f)
        f = open(inputFile*"index_q.csv"); readline(f)
        for l in eachline(f); l=string.(split(replace(l,"\""=>""), ',')); push!(qi, ind(l[2],l[3],l[4])) end
        close(f)
    end
end

function readIOTables(year, tfile, vfile, yfile, qfile)

    global mTables

    tb = tables(year, length(ti), length(vi), length(yi), length(qi))

    f = open(tfile)
    i = 1
    for l in eachline(f); tb.t[i,:] = [parse(Float64, x) for x in split(l, ',')]; i += 1 end
    close(f)
    f = open(vfile)
    i = 1
    for l in eachline(f); tb.v[i,:] = [parse(Float64, x) for x in split(l, ',')]; i += 1 end
    close(f)
    f = open(yfile)
    i = 1
    for l in eachline(f); tb.y[i,:] = [parse(Float64, x) for x in split(l, ',')]; i += 1 end
    close(f)
    f = open(qfile)
    i = 1
    nt = length(ti)
    for l in eachline(f); tb.q[i,:] = [parse(Float64, x) for x in split(l, ',')][1:nt]; i += 1 end
    close(f)

    mTables[year] = tb
end

function rearrangeIndex(; qmode = "")

    if qmode == "" || qmode == "I_CHG_CO2"; ql = deleteat!(collect(10:64), [12,13,43,44,45,46,47])  # I-GHG-CO2 emissions (Gg)
    elseif qmode == "PRIMAP"; ql = [2570]                                                           # PRIMAP|KYOTOGHGAR4|TOTAL|GgCO2eq
    end

    global ti = ti[1:end-1]
    global vi = vi[1:end-6]
    global yi = yi[1:end-6]
    global qi = qi[ql]

    global natList
    for t in ti; if !(t.nation in natList); push!(natList, t.nation) end end
end

function rearrangeTables(year; qmode = "")
    global ti, vi, yi, qi, mTables
    if qmode == "" || qmode == "I_CHG_CO2"; ql = deleteat!(collect(10:64), [12,13,43,44,45,46,47])  # I-GHG-CO2 emissions (Gg)
    elseif qmode == "PRIMAP"; ql = [2570]                                                           # PRIMAP|KYOTOGHGAR4|TOTAL|GgCO2eq
    end

    tb = mTables[year]

    nt = length(ti)
    tb.t = tb.t[1:nt, 1:nt]
    tb.v = tb.v[1:length(vi), 1:nt]
    tb.y = tb.y[1:nt, 1:length(yi)]
    tb.q = tb.q[ql, 1:nt]
end

function getSectorData(year, sector, subst_sector = [])
    if length(subst_sector)>0 && haskey(subst_sector, year); global sec[year] = vcat(sector[year], subst_sector[year])
    else; global sec[year] = sector[year]
    end
end

function getDomesticData(year, nation, expendMat, householdID)
    global hhid[year] = householdID[year][nation]
    expMat = expendMat[year][nation]

    if size(expMat,1) == length(sec[year]) && size(expMat,2) == length(hhid[year]); global hhExp[year] = expMat
    elseif size(expMat,2) == length(sec[year]) && size(expMat,1) == length(hhid[year]); global hhExp[year] = transpose(expMat)
    else println("Matrices sizes don't match: expMat,", size(expMat), "\thhid,", size(hhid[year]), "\tsec,", size(sec[year]))
    end
end

# function readEmissionRates(year, categoryFile, convertingFile)
#
#     global sec, deSecList, dePerEUR, deIdx, euA3, euA2
#     nats = []
#
#     # read DE-Expenditure matching sectors
#     deSecList[year] = Array{String, 1}()
#     deHbsList[year] = Array{String, 1}()
#     deHbsSec[year] = Dict{String, String}()
#     deIdx[year] = Dict{String, Array{Int, 1}}()
#     xf = XLSX.readxlsx(categoryFile)
#     sh = xf["Nation"]
#     for r in XLSX.eachrow(sh)
#         if XLSX.row_number(r) > 1
#             push!(nats, r[2])
#             euA2[r[2]] = r[1]
#             euA3[r[2]] = r[3]
#         end
#     end
#     sh = xf["DE_sector"]
#     for r in XLSX.eachrow(sh)
#         if XLSX.row_number(r) > 1 && r[1] == year
#             if !(r[2] in deSecList[year])
#                 push!(deSecList[year], r[2])
#                 deIdx[year][r[2]] = Array{Int, 1}()
#             end
#             if r[3] in sec[year]
#                 if !(r[3] in deHbsList[year]); push!(deHbsList[year], r[3]) end
#                 deHbsSec[year][r[3]] = r[2]
#                 push!(deIdx[year][r[2]], findfirst(x->x==r[3], sec[year]))
#             end
#         end
#     end
#     sort!(deHbsList[year])
#     close(xf)
#
#     # read DE/EUR exchange rates
#     dePerEUR[year] = Dict{String, Dict{String, Float64}}()
#     for n in nats; dePerEUR[year][euA2[n]] = Dict{String, Float64}() end
#     f = open(convertingFile); readline(f)
#     for l in eachline(f)
#         s = string.(split(l, '\t'))
#         if parse(Int, s[1]) == year && s[3] in deSecList[year]; dePerEUR[year][euA2[s[2]]][s[3]] = parse(Float64, s[4]) end
#     end
#     close(f)
# end

# function calculateDirectEmission(year, nation)
#
#     global sec, hhid, hhExp, directCE, deSecList, deHbsList, deHbsSec, dePerEUR
#     ns = length(deHbsList[year])
#     nh = length(hhid[year])
#     de = zeros(Float64, ns, nh)
#
#     for i = 1:ns
#         s = deHbsList[year][i]
#         des = deHbsSec[year][s]     # HBS code corresponding DE sector
#         idx = findfirst(x->x==s, sec[year])
#         der = dePerEUR[year][nation][des]
#         for j = 1:nh; de[i,j] = der * hhExp[year][idx,j] end
#     end
#     directCE[year] = de
# end

function readEmissionData(year, nat_dict, filepath; output_path = "", output_tag = "", integrate = false)

    global de_list, de_pr_link, de_emits, de_enbal, de_massc, de_enerc, de_price, de_price_unit
    nat_code = sort(collect(keys(nat_dict)))
    nat_name = [x == "SK" ? "Slovak Republic" : nat_dict[x] for x in nat_code]
    nn = length(nat_code)

    sector_file = filepath * "Emission_sectors.txt"
    emiss_file = filepath * "Emission_ktCO2_" * string(year) * ".xlsx"
    emiss_road_file = filepath * "Emission_road_ktCO2_" * string(year) * ".xlsx"
    emiss_res_file = filepath * "Emission_residential_ktCO2_" * string(year) * ".xlsx"
    enbal_file = filepath * "EnergyBlanaces_TJ_" * string(year) * ".xlsx"
    enbal_road_file = filepath * "EnergyBlanaces_road_TJ_" * string(year) * ".xlsx"
    enbal_res_file = filepath * "EnergyBlanaces_residential_TJ_" * string(year) * ".xlsx"
    mass_conv_file = filepath * "Barrels_per_Tonne_EU_" * string(year) * ".xlsx"
    ener_conv_file = filepath * "KJ_per_Kg_EU_" * string(year) * ".xlsx"
    price_trans_file = filepath * "EnergyPrices_Transport.xlsx"
    price_other_file = filepath * "EnergyPrices_OtherProducts.xlsx"

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
        de_idx = []
        xf = XLSX.readxlsx(filename)
        tb = xf[string(year)][:]
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
    yr_idx = findfirst(x->x==string(year), tb[1,3:end]) + 2
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
    yr_idx = findfirst(x->x==string(year), tb[1,3:end]) + 2
    sorts = Array{String, 1}()
    pri_ot = Dict{String, Dict{String, Tuple{Float64, String}}}()   # {nation, {fuel_sort, {price, unit}}}
    for i in filter(x->!ismissing(tb[x,1]) && rsplit(tb[x,1], '.', limit=2)[2] in ["Residential","Transport"], collect(3:size(tb)[1]))
        n, s, t = split(tb[i,1], '.')
        ft, ut = strip.(rsplit(s, ('(',')'), limit=3))
        pri = tb[i+1, yr_idx]
        if isa(pri, Number) && tb[i+1, 2] == unit_lab
            if !(ft in sorts); push!(sorts, ft) end
            if !haskey(pri_ot, n); pri_ot[n] = Dict{String, Tuple{Float64, String}}() end
            if occursin(" ", ut) && tryparse(Float64, split(ut, " ")[1]) != nothing
                scl, ut = split(ut, " ")
                pri /= parse(Float64, scl)
            end
            pri_ot[n][ft] = (pri, "USD/" * ut)
        end
    end
    close(xf)

    sort!(sorts); ns = length(sorts)
    all_avg, grp_avg, all_cnt, grp_cnt = zeros(Float64, ns), zeros(Float64, ns), zeros(Int, ns), zeros(Int, ns)
    all_unit = ["" for i = 1:ns]
    for n in collect(keys(pri_ot)), i = 1:ns
        s = sorts[i]
        if haskey(pri_ot[n], s)
            all_avg[i] += pri_ot[n][s][1]; all_cnt[i] += 1
            if n in nat_name; grp_avg[i] += pri_ot[n][s][1]; grp_cnt[i] += 1 end
            if all_unit[i] == ""; all_unit[i] = pri_ot[n][s][2]
            elseif all_unit[i] != pri_ot[n][s][2]; println("Price units are different: ", s, "\t", all_unit[i], "\t", pri_ot[n][s][2])
            end
        end
    end
    all_avg ./= all_cnt
    grp_avg ./= grp_cnt

    n_ps = length(price_item)
    append!(price_item, sorts)
    for n in nat_code; append!(de_price[year][n], zeros(Float64, ns)); append!(de_price_unit[year][n], ["" for i=1:ns]) end
    for i = 1:nn, j = 1:ns
        n, n_name, s = nat_code[i], nat_name[i], sorts[j]
        if haskey(pri_ot, n_name) && haskey(pri_ot[n_name], s); de_price[year][n][j+n_ps], de_price_unit[year][n][j+n_ps] = pri_ot[n_name][s]
        elseif grp_avg[j] > 0; de_price[year][n][j+n_ps], de_price_unit[year][n][j+n_ps] = grp_avg[j], all_unit[j]
        elseif all_avg[j] > 0; de_price[year][n][j+n_ps], de_price_unit[year][n][j+n_ps] = all_avg[j], all_unit[j]
        else println(sorts[j], " does not have any price data in ", year)
        end
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
        f = open(output_path * "Price_" * string(year) * "_EU.txt", "w")
        for pit in price_item; print(f, "\t", pit) end; println(f)
        for n in nat_code; print(f, n); for i = 1:length(price_item); print(f, "\t", de_price[year][n][i], " ", de_price_unit[year][n][i]) end; println(f) end
        close(f)
        f = open(output_path * "Price_others_" * string(year) * ".txt", "w")
        for s in sorts; print(f, "\t", s) end; println(f)
        for n in sort(collect(keys(pri_ot)))
            print(f, n); for s in sorts; print(f, "\t", haskey(pri_ot[n], s) ? string(pri_ot[n][s][1])*" "*pri_ot[n][s][2] : "") end; println(f)
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

function calculateEmissionRates(year; output = "", currency = "EUR")

    global de_code, de_list, de_pr_link, de_emits, de_enbal, de_massc, de_enerc, de_price, de_price_unit, de_price_item, de_intens, de_energy

    emi, bal, mas, ene = de_emits[year], de_enbal[year], de_massc[year], de_enerc[year]
    pri_itm, pri, pri_unt = de_price_item[year], de_price[year], de_price_unit[year]
    de_intens[year] = Dict{String, Array{Float64, 1}}()
    de_energy[year] = Dict{String, Array{Float64, 1}}()
    nats = sort(collect(keys(emi)))
    nn, nd, np = length(nats), length(de_list), length(pri_itm)
    intens = zeros(Float64, nn, np)         # tCO2 / EUR
    tot_emits = zeros(Float64, nn, np)      # TJ
    emit_rate = zeros(Float64, nn, nd)      # ktCO2 / TJ

    for i = 1:nn, j = 1:nd
        n = nats[i]
        if emi[n][j] > 0 && bal[n][j] > 0; emit_rate[i,j] =  emi[n][j] / bal[n][j] end
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
    for j = 1:np
        int_avg = mean(filter(x -> x>0, intens[:, j]))
        for i in filter(x -> intens[x, j]==0, collect(1:nn)); intens[i, j] = int_avg end
    end
    for i = 1:nn
        de_intens[year][nats[i]] = intens[i, :]
        de_energy[year][nats[i]] = tot_emits[i, :]
    end

    if length(output) > 0
        mkpath(rsplit(output, '/', limit = 2)[1])
        f = open(output, "w")
        pri_idx = filter(x->haskey(de_pr_link[year], pri_itm[x]), collect(1:np))
        for i in pri_idx; print(f, "\t", pri_itm[i]) end; println(f)
        for i = 1:nn; print(f, nats[i]); for j in pri_idx; print(f, "\t", intens[i, j]) end; println(f) end
        close(f)
    end
end

function printEmissionConvRates(year, outputFile; emit_unit = "tCO2", curr_unit = "EUR")

    global de_pr_link, de_price_item, de_intens, de_energy
    nats = sort(collect(keys(de_intens[year])))

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    println(f, "Year\tNation\tDE_code\tDE_sector\tFactor\tUnit\tEmission\tUnit")
    for n in nats, j in filter(x->haskey(de_pr_link[year], de_price_item[year][x]), collect(1:length(de_price_item[year])))
        print(f, year, "\t", n, "\t", j, "\t", de_price_item[year][j], "\t")
        println(f, de_intens[year][n][j], "\t", emit_unit, "/", curr_unit, "\t", de_energy[year][n][j], "\t", "TJ")
    end
    close(f)
end

function readEmissionIntensity(year, nations, sectorFile, intensityFile; emit_unit = "tCO2", curr_unit = "EUR")

    global de_sectors, de_intens, de_price_unit, de_energy
    if !haskey(de_sectors, year); de_sectors[year] = Array{String, 1}() end
    if !haskey(de_intens, year); de_intens[year] = Dict{String, Array{Float64, 1}}() end
    if !haskey(de_price_unit, year); de_price_unit[year] = Dict{String, Array{String, 1}}() end
    if !haskey(de_energy, year); de_energy[year] = Dict{String, Array{Float64, 1}}() end
    dsl = de_sectors[year]

    f_sep = getValueSeparator(sectorFile)
    f = open(sectorFile)
    si = findfirst(x->x=="Sector", string.(strip.(split(readline(f), f_sep))))
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))[si]
        if !(s in dsl); push!(dsl, s) end
    end
    close(f)

    nd = length(dsl)
    for n in nations; de_intens[year][n], de_price_unit[year][n], de_energy[year][n] = zeros(Float64, nd), ["" for i=1:nd], zeros(Float64, nd) end

    f_sep = getValueSeparator(intensityFile)
    essential = ["Year", "Nation", "DE_code", "DE_sector", "Factor", "Unit", "Emission", "Unit"]
    f = open(intensityFile)
    title = string.(strip.(split(readline(f), f_sep)))
    i = [findfirst(x->x==et, title) for et in essential]
    ener_unit = ""
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        if s[i[1]] == string(year) && s[i[2]] in nations
            unit = string.(strip.(split(s[i[6]], '/')))
            if [emit_unit, curr_unit] == unit
                di = findfirst(x->x==s[i[4]], dsl)
                de_intens[year][s[i[2]]][di], de_price_unit[year][s[i[2]]][di] = parse(Float64, s[i[5]]), s[i[6]]
                if ener_unit == ""; ener_unit = s[i[8]] end
                if ener_unit == s[i[8]]; de_energy[year][s[i[2]]][di] = parse(Float64, s[i[7]]) end
            end
        end
    end
    close(f)
end

function buildDeConcMat(year, nation, deCodeFile, concFile; norm = false, energy_wgh = false, output = "", group = "EU")

    global sec, de_sectors, de_energy, concMatDe

    nd, nc = length(de_sectors[year]), length(sec[year])
    concMatDe[year] = zeros(Float64, nd, nc)
    nats = collect(keys(de_energy[year]))
    ene_avg = [mean(filter(x->x>0, [de_energy[year][n][j] for n in nats])) for j=1:nd]

    essential = ["Nation", "DE_sector", "CES/HBS_code", "Weight"]
    f_sep = getValueSeparator(concFile)
    f = open(concFile)
    title = string.(strip.(split(readline(f), f_sep)))
    i = [findfirst(x->x==et, title) for et in essential]

    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        if s[i[1]] in [nation, group]
            desec, cescode, wgh = s[i[2]], s[i[3]], parse(Float64, s[i[4]])
            di, ci = findfirst(x->x==desec, de_sectors[year]), findfirst(x->x==cescode, sec[year])
            if energy_wgh; concMatDe[year][di, ci] += wgh * (de_energy[year][nation][di]>0 ? de_energy[year][nation][di] : ene_avg[di])
            else concMatDe[year][di, ci] += wgh
            end
        end
    end
    close(f)

    if norm
        for i = 1:nc
            ces_sum = sum(concMatDe[year][:,i])
            if ces_sum > 0; concMatDe[year][:, i] /= ces_sum end
        end
    end

    if length(output)>0
        mkpath(rsplit(output, '/', limit = 2)[1])
        f = open(output, "w")
        for sc in sec[year]; print(f, "\t", sc) end; println(f)
        for i = 1:nd
            print(f, de_sectors[year][i])
            for j = 1:nc; print(f, "\t", concMatDe[year][i, j]) end
            println(f)
        end
        close(f)
    end

    return concMatDe[year]
end

function calculateDirectEmission(year, nation; sparseMat = false, enhance = false, full = false)

    global sec, de_sectors, hhid, hhExp, directCE, concMatDe, de_intens, de_price_unit
    hl, sl, dsl, cmn, he = hhid[year], sec[year], de_sectors[year], concMatDe[year], hhExp[year]
    dit = transpose(de_intens[year][nation])
    nh, ns, nds = length(hl), length(sl), length(dsl)

    de = zeros(Float64, ns, nh)
    if full || enhance; dit_cmn = dit * cmn end
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
    directCE[year] = de
end

function assembleConcMat(year, conMat)

    global concMat, mTables
    global natList, sec, ti
    tb = mTables[year]
    ns = length(sec[year])

    # check whether a nation's account have 'Commodities' entities
    chk = Dict{String, Bool}()      # {A3, (true: has 'Commodities', false: only 'Industries')}
    for t in ti
        if !haskey(chk, t.nation); chk[t.nation] = false end
        if !chk[t.nation] && t.entity == "Commodities"; chk[t.nation] = true end
    end

    # count number of 'Industries' entities by nation: 0 means having only 'Industries' or only 'Commodities'.
    cnt = zeros(Int, length(natList))
    for t in filter(x->chk[x.nation] && x.entity == "Industries", ti)
        cnt[findfirst(x->x==t.nation, natList)] += 1
    end

    # assemble concordance matrices
    cMat = zeros(Float64, 0, ns)
    for i = 1:length(natList)
        if cnt[i]>0; cMat = vcat(cMat, zeros(Float64, cnt[i], ns)) end
        cMat = vcat(cMat, conMat[natList[i]])
    end

    concMat[year] = cMat

    return concMat[year], ti, sec[year]
end

function buildWeightedConcMat(year, nat; output="") # feasical year, nation A3, concordance matrix (Eora, Nation)

    global concMat, concMatWgh, mTables
    global natList, sec, ti, yi
    tb = mTables[year]
    ns = length(sec[year])
    cMat = concMat[year][:,:]

    # get final demand of nation 'nat'
    ye = tb.y[:,findfirst(x->x.nation==nat && x.sector=="Household final consumption P.3h", yi)]

    # reflect ratios of Eora final demand accounts
    for j = 1:ns
        cMat[:, j] .*= ye
        tsum = sum(cMat[:,j])
        cMat[:, j] /= tsum
    end

    concMatWgh[year] = cMat

    # print concordance matrix
    mkpath(rsplit(output, '/', limit = 2)[1])
    if length(output)>0
        f = open(output, "w")
        print(f,"Nation,Entity,Sector");for i=1:ns; print(f,",",sec[year][i]) end; println(f)
        for i=1:size(concMatWgh[year],1)
            print(f,ti[i].nation,",",ti[i].entity,",",ti[i].sector)
            for j=1:size(concMatWgh[year],2); print(f,",",concMatWgh[year][i,j]) end; println(f)
        end
        close(f)
    end

    return concMatWgh[year], ti, sec[year]
end

# function buildWeightedConcMat(year, nat, conMat; output="") # feasical year, nation A3, concordance matrix (Eora, Nation)
#
#     global concMatWgh, mTables
#     global natList, sec, ti, yi
#     tb = mTables[year]
#     ns = length(sec[year])
#
#     # get final demand of nation 'nat'
#     if length(nat) > 0; ye = tb.y[:,findfirst(x->x.nation==nat && x.sector=="Household final consumption P.3h", yi)] end
#
#     # check whether a nation's account have 'Commodities' entities
#     chk = Dict{String, Bool}()      # {A3, (true: has 'Commodities', false: only 'Industries')}
#     for t in ti
#         if !haskey(chk, t.nation); chk[t.nation] = false end
#         if !chk[t.nation] && t.entity == "Commodities"; chk[t.nation] = true end
#     end
#
#     # count number of 'Industries' entities by nation: 0 means having only 'Industries' or only 'Commodities'.
#     cnt = zeros(Int, length(natList))
#     for t in ti
#         if chk[t.nation] && t.entity == "Industries"
#             cnt[findfirst(x->x==t.nation, natList)] += 1
#         end
#     end
#
#     # assemble concordance matrices
#     cMat = zeros(Float64, 0, ns)
#     for i = 1:length(natList)
#         if cnt[i]>0; cMat = vcat(cMat, zeros(Float64, cnt[i], ns)) end
#         cMat = vcat(cMat, conMat[natList[i]])
#     end
#
#     # reflect ratios of Eora final demand accounts
#     if length(nat) > 0
#         for j = 1:ns
#             cMat[:, j] .*= ye
#             tsum = sum(cMat[:,j])
#             cMat[:, j] /= tsum
#         end
#     end
#
#     concMatWgh[year] = cMat
#
#     # print concordance matrix
#     mkpath(rsplit(output, '/', limit = 2)[1])
#     if length(output)>0
#         f = open(output, "w")
#         print(f,"Nation,Entity,Sector");for i=1:ns; print(f,",",sec[year][i]) end; println(f)
#         for i=1:size(concMatWgh[year],1)
#             print(f,ti[i].nation,",",ti[i].entity,",",ti[i].sector)
#             for j=1:size(concMatWgh[year],2); print(f,",",concMatWgh[year][i,j]) end; println(f)
#         end
#         close(f)
#     end
#
#     return concMatWgh[year], ti, sec[year]
# end

function calculateLeontief(year)

    global mTables, ti, vi, yi, qi, lti
    nt = length(ti)
    tb = mTables[year]

    x = sum(tb.t, dims = 1) +  sum(tb.v, dims = 1)  # calculate X
    f = sum(tb.q, dims = 1) ./ x                    # calculate EA
    lt = Matrix{Float64}(I, nt, nt)                 # calculate Leontief matrix
    for i = 1:nt; for j = 1:nt; lt[i,j] -= tb.t[i,j] / x[j] end end
    lti = inv(lt)
    for i = 1:nt; lti[i,:] *= f[i] end
end


function calculateIndirectEmission(year, sparseMat = false, elapChk = 0)

    global emissions, mTables, concMatWgh, lti
    global sec, hhid, hhExp
    global ti, vi, yi, qi

    tb = mTables[year]
    nt = length(ti)
    ns = length(sec[year])
    nh = length(hhid[year])

    # calculate emission, by Eurostat micro-data sectors, by Eora T matrix sectors
    e = zeros(Float64, ns, nh)

    if sparseMat
        concMatS = SparseArrays.sortSparseMatrixCSC!(sparse(concMatWgh[year]), sortindices=:doubletranspose)
        ltiS = SparseArrays.sortSparseMatrixCSC!(sparse(lti), sortindices=:doubletranspose)
        concMatWgh[year] = []
        lti = []
    end

    st = time()     # check start time
    for i = 1:ns
        hce = zeros(Float64, ns, nh)
        hce[i,:] = hhExp[year][i,:]

        if sparseMat
            hceS = SparseArrays.sortSparseMatrixCSC!(sparse(hce), sortindices=:doubletranspose)
            hce = []
            ebe = ltiS * concMatS * hceS
        else ebe = lti * concMatWgh[year] * hce       # household emission by Eora sectors
        end
        e[i,:] = sum(ebe, dims=1)       # calculate total emission (=sum of Eora emissions) of each nation sector

        if elapChk > 0   # check elapsed and remained time
            elap = floor(Int, time() - st)
            (eMin, eSec) = fldmod(elap, 60)
            (eHr, eMin) = fldmod(eMin, 60)
            (rMin, rSec) = fldmod(floor(Int, (elap / i) * (ns - i)), 60)
            (rHr, rMin) = fldmod(rMin, 60)

            if i%elapChk == 0
                println(i,"/",ns," iterations, ",eHr,":",eMin,":",eSec," elapsed, ",rHr,":",rMin,":",rSec," remained")
            end
        end
    end

    emissions[year] = e

    return e, sec[year], hhid[year]
end

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

function printIndirectEmissions(year, outputFile)

    global sec, hhid, emissions

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    e = emissions[year]

    ns = length(sec[year])
    nh = length(hhid[year])

    for h in hhid[year]; print(f, "\t", h) end
    println(f)
    for i = 1:ns
        print(f, sec[year][i])
        for j = 1:nh; print(f, "\t", e[i,j]) end
        println(f)
    end

    close(f)
end

# function printDirectEmissions(year, outputFile)
#
#     global de_sectors, hhid, directCE
#
#     mkpath(rsplit(outputFile, '/', limit = 2)[1])
#     f = open(outputFile, "w")
#     de = directCE[year]
#     ns, nh = length(de_sectors[year]), length(hhid[year])
#
#     for h in hhid[year]; print(f, "\t", h) end
#     println(f)
#     for i = 1:ns
#         print(f, de_sectors[year][i])
#         for j = 1:nh; print(f, "\t", de[i,j]) end
#         println(f)
#     end
#     close(f)
# end

function printEmissions(year, outputFile; mode = "ie")

    global hhid, sec, directCE, indirectCE
    hl, sl = hhid[year], sec[year]
    ns, nh = length(sl), length(hl)

    if lowercase(mode) == "ie"; em = indirectCE[year]
    elseif lowercase(mode) == "de"; em = directCE[year]
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

function initiate()
    global abb = Dict{String, String}()    # Nation name's A3 abbreviation, {Nation, A3}
    global euA2 = Dict{String, String}()   # EU nation name's A2 abbreviation, {Nation, A2}
    global euA3 = Dict{String, String}()   # EU nation name's A3 abbreviation, {Nation, A3}
    global natList = Array{String, 1}()    # Nation A3 list
    global ti = Array{idx, 1}()     # index T
    global vi = Array{idx, 1}()     # index V
    global yi = Array{idx, 1}()     # index Y
    global qi = Array{ind, 1}()     # index Q

    global sec = Dict{Int, Array{String, 1}}()          # Household expenditure sectors: {year, {sector}}
    global hhid = Dict{Int, Array{String, 1}}()         # Household ID: {year, {hhid}}
    global hhExp = Dict{Int, Array{Float64, 2}}()       # Households enpenditure: {year, {Nation sectors, households}}
    global concMat = Dict{Int, Array{Float64, 2}}()     # Assembled concordance matrix {Eora sectors, Nation sectors}
    global concMatWgh = Dict{Int, Array{Float64, 2}}()  # Weighted concordance matrix {Eora sectors, Nation sectors}
    global concMatDe = Dict{Int, Array{Float64, 2}}()   # concordance matrix sets for direct emission

    global lti = []                            # Inversed Leontief matrix
    global eoraExp = Array{Float64, 2}         # Transformed households expenditure, {Eora sectors, households}
    global mTables = Dict{Int16, tables}()     # {Year, tables}
    global emissions = Dict{Int16, Array{Float64, 2}}()    # carbon footprint

    # direct carbon emission variables
    global directCE = Dict{Int, Array{Float64, 2}}()        # direct carbon emission
    global deSecList = Dict{Int, Array{String, 1}}()        # direct emission sectors: {year, {DE sector}}
    global deHbsList = Dict{Int, Array{String, 1}}()        # direct emission HBS sectors: {year, {HBS code}}
    global deHbsSec = Dict{Int, Dict{String, String}}()     # HBS sector - DE sector link: {year, {HBS code, DE sector}}
    global dePerEUR = Dict{Int, Dict{String, Dict{String, Float64}}}() # converting rate from EUR to CO2: {year, {Nation, {DE category, tCO2/EUR}}}
    global deIdx = Dict{Int, Dict{String, Array{Int, 1}}}() # Direct emission sector matched EU expenditure sector index: {year, {DE sector, {HBS expenditure index}}}

    global de_list = Array{String, 1}()                             # DE IEA sector label list
    global de_pr_link = Dict{Int, Dict{String, Array{String, 1}}}() # DE sector - IEA price sector link: {year, {price item, {De item}}}
    global de_emits = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA emission: {year, {nation, {emission}}}
    global de_enbal = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA energy balance: {year, {nation, {energy}}}
    global de_massc = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA mass to volume converting rate: {year, {nation, {barrels/tonne}}}
    global de_enerc = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA mass to energy converting rate: {year, {nation, {KJ/Kg}}}

    global de_price = Dict{Int, Dict{String, Array{Float64, 1}}}()  # fuel prices: {year, {nation, {prices}}}
    global de_price_unit = Dict{Int, Dict{String, Array{String, 1}}}()  # fuel price unit: {year, {nation, {units}}}
    global de_price_item = Dict{Int, Array{String, 1}}()            # IEA price items: {year, {items}}
    global de_intens = Dict{Int, Dict{String, Array{Float64, 1}}}() # IEA price items' CO2 intensity: {year, {nation, {tCO2/EUR}}}
    global de_energy = Dict{Int, Dict{String, Array{Float64, 1}}}()  # IEA price items' CO2 energy balance: {year, {nation, {TJ}}}
    global de_sectors = Dict{Int, Array{String, 1}}()           # direct emission sectors: {year, {DE sector}}
end

end
