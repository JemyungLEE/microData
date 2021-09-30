module EmissionDecomposer

# Developed date: 27. Jul. 2021
# Last modified date: 28. Sep. 2021
# Subject: Decompose EU households' carbon footprints
# Description: Process for Input-Output Structural Decomposition Analysis
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

include("MicroDataReader.jl")
include("EmissionEstimator.jl")
include("EmissionCategorizer.jl")

using XLSX
using Statistics
using LinearAlgebra
using .MicroDataReader
using .EmissionEstimator
using .EmissionCategorizer

mdr = MicroDataReader
ee = EmissionEstimator
ec = EmissionCategorizer

mutable struct factors
    f::Array{Float64, 1}        # Emission factors: {Eora t-index}
    l::Array{Float64, 2}        # Leontief matrix: {Eora t-index, Eora t-index}
    p::Array{Float64, 1}        # Populations: {region}
    cepc::Array{Float64, 1}     # Consumption expenditures per capita: {region}
    cspf::Array{Float64, 2}     # Regional household expenditure profile: {Eora t-index, region}
    de::Array{Float64, 1}       # direct emission: {region}

    cspfbc::Array{Float64, 3}   # Regional household expenditure profile by expenditure category: {Eora t-index, category, region}
    cpbc::Array{Float64, 2}     # Consumption profile by expenditure category: {category, region}

    cepcbc::Array{Array{Float64, 3}, 1} # Consumption expenditures per capita by expenditure category: {category, {category, category, region}}

    function factors(nr=0, nt=0; nc=0, factor_by_cat = false, f=zeros(0), l=zeros(0,0), p=zeros(0), exp_pc=zeros(0), exp_prof=zeros(0,0), de=zeros(0), int_share=zeros(0,0,0), dom_share=zeros(0,0), exp_pc_cat=Array{Array{Float64,3},1}())
        if nr > 0 && nt > 0 && length(f) == length(l) == length(p) == 0
            if nc == length(exp_pc) == length(exp_prof) == length(de) == 0
                new(zeros(nt), zeros(nt, nt), zeros(nr), zeros(nr), zeros(nt, nr), zeros(nr), zeros(0,0,0), zeros(0,0), exp_pc_cat)
            elseif !factor_by_cat && nc > 0 && length(exp_pc) == length(de) == length(int_share) == length(dom_share) == 0
                new(zeros(nt), zeros(nt, nt), zeros(nr), zeros(nr), zeros(0, 0), zeros(nr), zeros(nt, nc, nr), zeros(nc, nr), exp_pc_cat)
            elseif factor_by_cat && nc > 0 && length(de) == length(int_share) == length(exp_pc_cat) == 0
                new(zeros(nt), zeros(nt, nt), zeros(nr), zeros(0), zeros(0, 0), zeros(nr), zeros(nt, nc, nr), zeros(0, 0), [zeros(nc, nc, nr) for i=1:nc])
            end
        else new(f, l, p, exp_pc, exp_prof, de, int_share, dom_share, exp_pc_cat)
        end
    end
end

# mutable struct factors
#     f::Array{Float64, 1}        # Emission factors: {Eora t-index}
#     l::Array{Float64, 2}        # Leontief matrix: {Eora t-index, Eora t-index}
#     p::Array{Float64, 1}        # Populations: {region}
#     cepc::Array{Float64, 1}     # Consumption expenditures per capita: {region}
#     cspf::Array{Float64, 2}     # Regional household expenditure profile: {Eora t-index, region}
#     de::Array{Float64, 1}       # direct emission: {region}
#
#     function factors(nr=0, nt=0; f=Array{Float64,1}(), l=Array{Float64,2}(undef,0,0), p=Array{Float64,1}(), exp_pc=Array{Float64,1}(), exp_prof=Array{Float64,2}(undef,0,0), de=Array{Float64,1}())
#         if nr>0 && nt>0 && length(f)==0 && length(l)==0 && length(p)==0 && length(exp_pc)==0 && length(exp_prof)==0 && length(de)==0
#             new(zeros(nt), zeros(nt, nt), zeros(nr), zeros(nr), zeros(nt, nr), zeros(nr))
#         else new(f, l, p, exp_pc, exp_prof, de)
#         end
#     end
# end

mutable struct delta
    nth::Int                    # denotes this delta is for the n_th factor
    n_f::Int                    # number of factors
    w::Float64                  # weight of delta_factor
    subs::Array{Int, 1}         # list of subscripts for the components of other factors except n_th
    d::Array{Float64, 1}        # delta emission

    function delta(nth, n_factor; sub_list = [], weight = 0.0, delta_value = [])
        if length(sub_list) == 0; new(nth, n_factor, weight, zeros(Int, n_factor-1), delta_value)
        elseif length(sub_list) == n_factor-1; new(nth, n_factor, weight, sub_list, delta_value)
        else println("subscript list's length does not match number of factors - 1.")
        end
    end
end

global yr_list = Array{Int, 1}()            # year list: {YYYY}
global nat_list = Array{String, 1}()        # nation list: {A2}
global nat_name = Dict{String, String}()    # nation names: {Nation code, Name}
global cat_list = Array{String, 1}()        # category list
global pr_unts = Dict("day" => 1, "week" => 7, "month" => 30, "year" => 365)

global sc_list = Dict{Int, Array{String, 1}}()                              # HBS sectors: {year, {sector}}
global hh_list = Dict{Int, Dict{String, Array{String, 1}}}()                # Household ID: {year, {nation, {hhid}}}
global households = Dict{Int, Dict{String, Dict{String, mdr.household}}}()  # household dict: {year, {nation, {hhid, household}}}
global exp_table = Dict{Int, Dict{String, Array{Float64, 2}}}()             # household expenditure table: {year, {nation, {hhid, category}}}
global sc_cat = Dict{Int, Dict{String, String}}()                           # category dictionary: {year, {sector code, category}}
# global exp_table_conv = Dict{Int, Dict{String, Array{Float64, 2}}}()        # Base-year price converted household expenditure table: {year, {nation, {hhid, category}}}

global nuts = Dict{Int, Dict{String, String}}()                     # NUTS: {year, {code, label}}
global nutsByNat = Dict{Int, Dict{String, Array{String, 1}}}()      # NUTS code list: {year, {nation_code, {NUTS_code}}}
global nuts_list = Dict{Int, Array{String, 1}}()                    # NUTS code list: {year, {NUTS_code}}
global pops = Dict{Int, Dict{String, Float64}}()                    # Population: {year, {NUTS_code, population}}
global pop_list = Dict{Int, Dict{String, Dict{String, Float64}}}()  # Population list: {year, {nation_code, {NUTS_code, population}}}
global pop_label = Dict{Int, Dict{String, String}}()                # populaton NUTS label: {year, {NUTS_code, NUTS_label}}
global pop_linked_cd = Dict{Int, Dict{String, String}}()            # concordance NUTS code: {year, {NUTS code, replaced population NUTS code}}

global pop_density = Dict{Int, Dict{String, Float64}}()             # Population density: {year, {NUTS_code, density}}
global pops_ds = Dict{Int, Dict{String, Dict{Int, Float64}}}()      # Population by population density: {year, {NUTS_code, {density, population}}}

global cpi_list = Dict{Int, Dict{String, Array{String, 1}}}()       # Consumption price indexes: {year, {nation, {COICOP_category}}}
global cpis = Dict{Int, Dict{String, Dict{String, Float64}}}()      # Consumption price indexes: {year, {nation, {COICOP_category, CPI}}}
global scl_rate = Dict{Int, Dict{String, Dict{String, Float64}}}()  # CPI scaling rate: {year, {nation, {HBS code, rate}}}
global conc_mat = Dict{Int, Array{Float64, 2}}()                    # Assembled concordance matrix {Eora sectors, Nation sectors}
global conc_mat_wgh = Dict{Int, Dict{String, Array{Float64, 2}}}()  # Weighted concordance matrix: {year, {nation, {Eora sectors, Nation sectors}}}
global mrio_idxs = Array{ee.idx, 1}()                               # index T
global mrio_tabs = Dict{Int, ee.tables}()                           # MRIO tables: {Year, MRIO tables (t, v ,y , q)}
global mrio_tabs_conv = Dict{Int, Dict{String, ee.tables}}()        # Base-year price converted MRIO tables: {Year, {natoin, MRIO tables (t, v ,y , q)}}
# global conc_mat_conv = Dict{Int, Dict{String, Array{Float64, 2}}}() # Base-year price converted concordance matrix: {year, {nation, {Eora sectors, Nation sectors}}}

global nt_wgh = Dict{Int, Dict{String, Float64}}()                  # hhid's NUTS weight: {year, {hhid, weight}}
global in_emiss = Dict{Int, Dict{String, Array{Float64, 2}}}()      # indirect carbon emission: {year, {nation, {table}}}
global di_emiss = Dict{Int, Dict{String, Array{Float64, 2}}}()      # direct carbon emission: {year, {nation, {table}}}

global l_factor = Dict{Int, Array{Float64, 2}}()                    # Leontief matrix: {year, {Eora t-index, Eora t-index}}
global sda_factors = Dict{Int, Dict{String, factors}}()             # SDA factors: {year, {nation, factors}}
global dltByNat = Dict{String, Dict{Int, Any}}()                    # delta by factor, by nation: {nation, {factor, {target_year - base_year}}}
global deltas = Dict{Tuple{Int,Int}, Dict{String, Dict{String, Array{Float64, 1}}}}()  # Deltas of elements: {(target_year, base_year), {nation, {region, {factor}}}}
global ieByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()        # indirect CF by nation, NUTS: {year, {nation, {nuts}}}
global deByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()        # direct CF by nation, NUTS: {year, {nation, {nuts}}}
global popByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()       # population by nation, NUTS: {year, {nation, {nuts}}}
global expPcByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()     # total expenditures by nation, NUTS: {year, {nation, {nuts}}}

global ci_ie = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}()   # confidence intervals of indirect emission: {year, {nation, {NUTS, {lower, upper}}}}
global ci_de = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}()   # confidence intervals of indirect emission: {year, {nation, {NUTS, {lower, upper}}}}

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

function getPopGridded(year, inputFile; tag = ["_dense", "_inter", "_spars", "_total"])
    # 1:Densely populated (at least 500), 2:Intermediate (between 100 and 499)
    # 3:Sparsely populated (less than 100), 4:Total, (Unit: inhabitants/km2)

    global nat_list, pops, pops_ds
    ntag = length(tag)

    xf = XLSX.readxlsx(inputFile)
    if isa(year, Number); year = [year] end

    tb = xf["Pop_GP"][:]
    for y in year
        pops_ds[y] = Dict{String, Dict{Int, Float64}}()
        idx = [findfirst(x -> x == string(y) * t, tb[1,:]) for t in tag]
        for ri = 2:size(tb,1)
            nt = string(tb[ri, 1])
            pops_ds[y][nt] = Dict{Int, Float64}()
            for i = 1:ntag; pops_ds[y][nt][i] = tb[ri,idx[i]] end
        end
    end
    close(xf)
end

function readPopDensity(year, densityFile)

    global pop_density, pops
    pop_density[year] = Dict{String, Float64}()

    sep = getValueSeparator(densityFile)
    f = open(densityFile)
    idx = findfirst(x -> string(year) == x, strip.(string.(split(readline(f), sep))))

    for l in eachline(f)
        s = strip.(string.(split(l, sep)))
        nt = rsplit(s[1], ',', limit = 2)[end]
        if haskey(pops[year], nt) && tryparse(Float64, s[idx]) != nothing; pop_density[year][nt] = parse(Float64, s[idx]) end
    end
    close(f)
end

function filterPopByDensity(year; nuts_lv = 3)
    # 1:Densely populated (at least 500 inhabitants/km2),  2:Intermediate (between 100 and 499 inhabitants/km2)
    # 3:Sparsely populated (less than 100 inhabitants/km2), 9:Not specified

    global nat_list, pop_density, pops, pops_ds
    pops_ds[year] = Dict{String, Dict{Int, Float64}}()

    nlv = nuts_lv + 2
    nt_list = sort(collect(keys(pops[year])))

    for nt in nuts_list[year]
        nt_str = nt
        while nt_str[end] == '0'; nt_str = nt_str[1:end-1] end

        nl = length(nt_str)
        if nl < nlv
            if !haskey(pops_ds[year], nt_str); pops_ds[year][nt_str] = Dict(1 => 0, 2 => 0, 3 => 0) end
            nts = filter(x -> length(x) == nlv && x[1:nl] == nt_str, nt_list)

            for n in filter(x -> haskey(pop_density[year], x), nts)
                ds = pop_density[year][n]
                if ds >= 500; pops_ds[year][nt_str][1] += pops[year][n]
                elseif ds >= 100; pops_ds[year][nt_str][2] += pops[year][n]
                elseif ds > 0; pops_ds[year][nt_str][3] += pops[year][n]
                end
            end
        else println(nt, " is excluded in filtering by population density: NUTS_level ", nuts_lv)
        end
    end
end

function printPopByDens(year, outputFile)

    global pops, pops_ds

    nts = sort(collect(keys(pops_ds[year])))
    sep = getValueSeparator(outputFile)
    f = open(outputFile, "w")
    println(f, "NUTS", sep, "Population", sep, "Pop_dense", sep, "Pop_inter", sep, "Pop_sparse")

    println()

    for nt in nts;
        nt_p = nt
        while nt_p[end] == '0'; nt_p = nt_p[1:end-1] end
        println(f, nt, sep, pops[year][nt_p], sep, pops_ds[year][nt][1], sep, pops_ds[year][nt][2], sep, pops_ds[year][nt][3]) end
    close(f)
end

function filterNations()

    global hh_list, nat_list

    n_list = Array{String, 1}()
    for y in sort(collect(keys(hh_list)))
        if length(n_list) == 0; n_list = sort(collect(keys(hh_list[y])))
        else filter!(x -> x in sort(collect(keys(hh_list[y]))), n_list)
        end
    end
    nat_list = n_list
    return n_list
end

function filterNonPopDens(year, nations = []; pop_dens = [1,2,3])

    global nat_list, hh_list, households, nutsByNat
    if isa(year, Number); year = [year] end
    if isa(pop_dens, Number); pop_dens = [pop_dens] end
    if length(nations) == 0; nats = nat_list else nats = nations end

    nats_flt = nats[:]

    for y in year, n in nats
        hhs, nts = hh_list[y][n], nutsByNat[y][n]
        nh = length(hhs)
        nts_flt = nts[:]
        for r in nts
            idxs = filter(x -> households[y][n][hhs[x]].nuts1 == r, 1:nh)
            for pd in pop_dens
                if length(filter(x -> households[y][n][hhs[x]].popdens == pd, idxs)) == 0; filter!(x -> !(x in [r]), nts_flt) end
            end
        end
        if length(nts_flt) == 0
            filter!(x -> !(x in [n]), nats_flt)
            delete!(nutsByNat[y], n)
        else nutsByNat[y][n] = nts_flt
        end
    end
    if length(nations) == 0; nat_list = nats_flt end

    return nats_flt
end

function detectNations(file_path, target_year, base_year; factor_file_tag = "_factors.txt")

    nats = Dict{Int, Array{String, 1}}()
    for y in [target_year, base_year]
        nats[y] = Array{String, 1}()
        for f in readdir(file_path * string(y) * "/")
            if startswith(f, string(y)) && endswith(f, factor_file_tag); push!(nats[y], f[6:7]) end
        end
    end

    global nat_list = filter(x -> x in nats[base_year], nats[target_year])
end

function importData(; hh_data::Module, mrio_data::Module, cat_data::Module, nations = [])

    global yr_list, nat_name = hh_data.year_list, hh_data.nationNames
    global hh_list, households, exp_table, scl_rate, cpis = hh_data.hhsList, hh_data.mdata, hh_data.expTable, hh_data.sclRate, hh_data.cpis
    global mrio_idxs, mrio_tabs, sc_list, conc_mat = mrio_data.ti, mrio_data.mTables, mrio_data.sec, mrio_data.concMat
    global nt_wgh, in_emiss, di_emiss = cat_data.wghNuts, cat_data.indirectCE, cat_data.directCE
    global cat_list, nuts, sc_cat = cat_data.catList, cat_data.nuts, cat_data.cat
    global pops, pop_list, pop_label, pop_linked_cd = cat_data.pop, cat_data.popList, cat_data.poplb, cat_data.popcd
    global nat_list = length(nations) > 0 ? nations : hh_data.nations

end

function storeNUTS(year; cat_data::Module)

    global nuts_list, nutsByNat

    nuts_list[year] = Array{String, 1}()
    nutsByNat[year] = Dict{String, Array{String, 1}}()
    for n in nat_list
        nutsByNat[year][n] = filter(x->x in cat_data.giscdlist[year], cat_data.nutsList[year][n])
        append!(nuts_list[year], nutsByNat[year][n])
    end
end

function storeNutsWeight(; year = 0)

    global yr_list, nat_list, households, hh_list, nt_wgh
    if year > 0; yrs = [year] else yrs = yr_list end

    for y in yrs, hh in collect(keys(nt_wgh[y]))
        n, h = hh[1:2], hh[4:end]
        if h in hh_list[y][n]; households[y][n][h].weight_nt = nt_wgh[y][hh]
        else println(h, " household is not in the ", y, " year's list")
        end
    end
end

function storeConcMat(year, nation, concMat)

    global conc_mat_wgh

    if !haskey(conc_mat_wgh, year); conc_mat_wgh[year] = Dict{String, Array{Float64, 2}}() end
    conc_mat_wgh[year][nation] = concMat
end

function readMrioTable(year, mrioPath, file_tag)

    f = open(mrioPath * string(year) * "/" * string(year) * file_tag)
    tb = Array{Array{Float64, 1}, 1}()
    for l in eachline(f); push!(tb, [parse(Float64, x) for x in split(l, ',')]) end
    close(f)
    nr, nc = length(tb), length(tb[1])
    mrio_tb = zeros(Float64, nr, nc)
    for i = 1:nr; mrio_tb[i,:] = tb[i][:] end

    return mrio_tb
end

function setMrioTables(year, mrioPath; t="_eora_t.csv", tax="_eora_t_tax.csv", sub="_eora_t_sub.csv", v="_eora_v.csv", y="_eora_y.csv")

    t_bp = readMrioTable(year, mrioPath, t)
    t_tax = readMrioTable(year, mrioPath, tax)
    t_sub = readMrioTable(year, mrioPath, sub)
    v_bp = readMrioTable(year, mrioPath, v)
    y_bp = readMrioTable(year, mrioPath, y)

    return t_bp, t_tax, t_sub, v_bp, y_bp
end

function convertTable(year, nation, base_year, mrioPath; total_cp = "CP00", t_bp, t_tax, t_sub, v_bp, y_bp)
    # double deflation method

    global sc_list, scl_rate, conc_mat, mrio_tabs, mrio_tabs_conv, cpis, mrio_idxs
    sclr, mrio, tidx = scl_rate[year][nation], mrio_tabs[year], mrio_idxs

    if !haskey(mrio_tabs_conv, year); mrio_tabs_conv[year] = Dict{String, ee.tables}() end
    avg_scl = cpis[year][nation][total_cp] / cpis[year][nation][total_cp]
    cvr_conc = [sclr[c] for c in sc_list[year]]
    cmat = conc_mat[year]

    row_idx, col_idx = Array{Int, 1}(), Array{Int, 1}()
    nti = length(tidx)

    ind_nat, com_nat = Array{String, 1}(), Array{String, 1}()
    for i = 1:nti
        if tidx[i].entity == "Industries" && !(tidx[i].nation in ind_nat); push!(ind_nat, tidx[i].nation)
        elseif tidx[i].entity == "Commodities" && !(tidx[i].nation in com_nat); push!(com_nat, tidx[i].nation)
        end
    end
    col_idx = filter(i -> (tidx[i].nation in ind_nat && tidx[i].nation in com_nat) && tidx[i].entity == "Commodities", 1:nti)
    row_idx = filter(i -> i in col_idx || !(tidx[i].nation in ind_nat) || !(tidx[i].nation in com_nat), 1:nti)

    cvr_mrio = [r > 0 ? r : avg_scl for r in (sum(cmat .* cvr_conc', dims = 2) ./ sum(cmat, dims = 2))]

    mrio_conv = ee.tables(year, size(mrio.t, 1), size(mrio.v, 1), size(mrio.y, 2), size(mrio.q, 1))
    mrio_conv.t = t_bp[:,:]
    mrio_conv.t[:,col_idx] .*= cvr_mrio[col_idx]'
    mrio_conv.t[row_idx,:] .*= cvr_mrio[row_idx]
    mrio_conv.y = y_bp[:,:]
    mrio_conv.y[row_idx,:] .*= cvr_mrio[row_idx]
    mrio_conv.v  = v_bp[:,:]

    t_all = t_bp + t_tax + t_sub
    t_all[:,col_idx] .*= cvr_mrio[col_idx]'
    t_all[row_idx,:] .*= cvr_mrio[row_idx]

    x_out = vec(sum(mrio_conv.t, dims=2) + sum(mrio_conv.y, dims=2))
    x_in = vec(sum(t_all, dims=1))
    d_sum = x_out - x_in
    v_sum = vec(sum(mrio_conv.v, dims=1))
    r_sum = d_sum ./ v_sum
    mrio_conv.v .*= r_sum'
    nv = size(mrio_conv.v, 1)
    for i in filter(x -> d_sum[x] == v_sum[x] == 0, 1:nti); mrio_conv.v[:,i] = zeros(Float64, nv) end
    for i in filter(x -> abs(d_sum[x]) >0 && v_sum[x] == 0, 1:nti); mrio_conv.v[:,i] = [d_sum[i] / nv for j = 1:nv] end

    nti, nvi, nyi = size(mrio.t, 1), size(mrio.v, 1), size(mrio.y, 2)
    mrio_conv.t = mrio_conv.t[1:nti, 1:nti]
    mrio_conv.y = mrio_conv.y[1:nti, 1:nyi]
    mrio_conv.v = mrio_conv.v[1:nvi, 1:nti]
    mrio_conv.q = mrio.q[:,:]
    mrio_tabs_conv[year][nation] = mrio_conv
end

# function convertTable(year, nation, base_year; total_cp = "CP00")
#     # double deflation method
#
#     global sc_list, scl_rate, conc_mat, mrio_tabs, mrio_tabs_conv, cpis, mrio_idxs
#     sclr, mrio, tidx = scl_rate[year][nation], mrio_tabs[year], mrio_idxs
#
#     if !haskey(mrio_tabs_conv, year); mrio_tabs_conv[year] = Dict{String, ee.tables}() end
#     avg_scl = cpis[year][nation][total_cp] / cpis[year][nation][total_cp]
#     cvr_conc = [sclr[c] for c in sc_list[year]]
#     cmat = conc_mat[year]
#
#     row_idx, col_idx = Array{Int, 1}(), Array{Int, 1}()
#     nti = length(tidx)
#
#     ind_nat, com_nat = Array{String, 1}(), Array{String, 1}()
#     for i = 1:nti
#         if tidx[i].entity == "Industries" && !(tidx[i].nation in ind_nat); push!(ind_nat, tidx[i].nation)
#         elseif tidx[i].entity == "Commodities" && !(tidx[i].nation in com_nat); push!(com_nat, tidx[i].nation)
#         end
#     end
#     col_idx = filter(i -> (tidx[i].nation in ind_nat && tidx[i].nation in com_nat) && tidx[i].entity == "Commodities", 1:nti)
#     row_idx = filter(i -> i in col_idx || !(tidx[i].nation in ind_nat) || !(tidx[i].nation in com_nat), 1:nti)
#
#     cvr_mrio = [r > 0 ? r : avg_scl for r in (sum(cmat .* cvr_conc', dims = 2) ./ sum(cmat, dims = 2))]
#
#     mrio_conv = ee.tables(year, size(mrio.t, 1), size(mrio.v, 1), size(mrio.y, 2), size(mrio.q, 1))
#     mrio_conv.t = mrio.t[:,:]
#     mrio_conv.t[:,col_idx] .*= cvr_mrio[col_idx]'
#     mrio_conv.t[row_idx,:] .*= cvr_mrio[row_idx]
#     mrio_conv.y = mrio.y .* cvr_mrio
#
#     row_sum = vec(sum(mrio_conv.t, dims=2) + sum(mrio_conv.y, dims=2))
#     col_sum = vec(sum(mrio_conv.t, dims=1))
#     d_sum = row_sum - col_sum
#     v_sum = vec(sum(mrio.v, dims=1))
#     r_sum = d_sum ./ v_sum
#     mrio_conv.v = mrio.v .* r_sum'
#     nv = size(mrio.v, 1)
#     for i in filter(x -> d_sum[x] == v_sum[x] == 0, 1:nti); mrio_conv.v[:,i] = zeros(Float64, nv) end
#     for i in filter(x -> abs(d_sum[x]) >0 && v_sum[x] == 0, 1:nti); mrio_conv.v[:,i] = [d_sum[i] / nv for j = 1:nv] end
#     mrio_conv.q = mrio.q
#     mrio_tabs_conv[year][nation] = mrio_conv
# end

function calculateLeontief(mrio_table)

    tb = mrio_table
    nt = size(tb.t, 1)

    x = sum(tb.t, dims = 1) +  sum(tb.v, dims = 1)      # calculate X
    lt = Matrix{Float64}(I, nt, nt)                     # calculate Leontief matrix
    for i = 1:nt; for j = 1:nt; lt[i,j] -= tb.t[i,j] / x[j] end end
    lti = inv(lt)

    return lti   # Leontied matrix
end

function calculateIntensity(mrio_table)

    tb = mrio_table
    nt = size(tb.t, 1)

    x = sum(tb.t, dims = 1) +  sum(tb.v, dims = 1)                  # calculate X
    f = [x[i] > 0 ? sum(tb.q, dims = 1)[i] / x[i] : 0.0 for i=1:nt] # calculate EA

    return f   # Leontied matrix, emission factor (intensity)
end

function decomposeFactors(year, baseYear, nation = "", mrioPath = ""; visible = false, pop_nuts3 = true, pop_dens = 0, mode="penta")
    # 1:Densely populated (at least 500), 2:Intermediate (between 100 and 499)
    # 3:Sparsely populated (less than 100)

    # 5-factors: f * L * p * tot_ce_pc * [con * hbs_profile] + DE
    # 6-factors: f * L * p * tot_ce_pc * [con * hbs_profile_by_category * hbs_proportion_by_category] + DE

    global mrio_tabs, mrio_tabs_conv, conc_mat_wgh, sda_factors, di_emiss, l_factor, cat_list, sc_cat, sc_list
    global nat_list, nutsByNat, hh_list, pops, pop_list, pop_linked_cd, pops_ds, pop_list_ds
    pt_sda, hx_sda, cat_sda = "penta", "hexa", "categorized"

    if isa(year, Number); year = [year] end
    if length(nation) == 0; nats = nat_list else nats = [nation] end
    if mode in [hx_sda, cat_sda]; nc = length(cat_list) end

    for y in year
        if y != baseYear; t_bp, t_tax, t_sub, v_bp, y_bp = setMrioTables(y, mrioPath) end
        for n in nats
            if visible; print("\t", n) end
            etab, cmat = exp_table[y][n], conc_mat_wgh[y][n]
            hhs, de, nts = hh_list[y][n], di_emiss[y][n], nutsByNat[y][n]
            nh = length(hhs)

            ft = factors()
            if y == baseYear
                if !haskey(l_factor, y); l_factor[y] = calculateLeontief(mrio_tabs[y]) end
                mrio, ft.l = mrio_tabs[y], l_factor[y]
            else
                convertTable(y, n, baseYear, mrioPath, t_bp = t_bp, t_tax = t_tax, t_sub = t_sub, v_bp = v_bp, y_bp = y_bp)
                mrio = mrio_tabs_conv[y][n]
                ft.l = calculateLeontief(mrio)
            end
            ft.f = calculateIntensity(mrio)
            nr, nt = length(nts), size(mrio.t, 1)
            if mode in [hx_sda, cat_sda]
                scl, sct = sc_list[y], sc_cat[y]
                nc, ns = length(cat_list), length(scl)
                ct_idx = [findall(x -> haskey(sct,x) && sct[x] == c, scl) for c in cat_list]
            end
            ft_p, ft_de = zeros(nr), zeros(nr)
            if mode == pt_sda; ft_cepc, ft_cspf = zeros(nr), zeros(nt, nr)
            elseif mode == hx_sda; ft_cepc, ft_cspf, ft_cpbc = zeros(nr), zeros(nt, nc, nr), zeros(nc, nr)
            elseif mode == cat_sda; ft_cepc, ft_cspf = [zeros(nc, nc, nr) for i=1:nc], zeros(nt, nc, nr)
            end

            for r in nts
                if pop_nuts3 && !(pop_dens in [1,2,3]); p_reg = pop_list[y][n][r]
                else
                    r_p = pop_linked_cd[y][r]
                    while r_p[end] == '0'; r_p = r_p[1:end-1] end
                    p_reg = (pop_dens in [1,2,3] ? pops_ds[y][r_p][pop_dens] : pops[y][r_p])
                end
                ri = findfirst(x -> x == r, nts)
                idxs = filter(x -> households[y][n][hhs[x]].nuts1 == r, 1:nh)
                if pop_dens in [1,2,3]; filter!(x -> households[y][n][hhs[x]].popdens == pop_dens, idxs) end

                wg_reg = [households[y][n][h].weight_nt for h in hh_list[y][n][idxs]]
                wg_sum = sum(wg_reg .* [households[y][n][h].size for h in hh_list[y][n][idxs]])
                etb_wg = wg_reg .* etab[idxs, :]

                if mode == pt_sda
                    et_sum = sum(etb_wg, dims=1)
                    ce_tot = sum(et_sum)
                    ce_pf = et_sum ./ ce_tot
                    ft_p[ri] = p_reg
                    ft_cepc[ri] = ce_tot / wg_sum
                    ft_cspf[:,ri] = cmat * ce_pf'
                elseif mode in [hx_sda, cat_sda]
                    et_sum = vec(sum(etb_wg, dims=1))
                    ct_pf = [sum(et_sum[ci]) for ci in ct_idx]
                    ce_tot = sum(ct_pf)
                    ce_pf = zeros(Float64, ns, nc)
                    for i = 1:nc; ce_pf[ct_idx[i], i] = (ct_pf[i] > 0 ? (et_sum[ct_idx[i]] ./ ct_pf[i]) : [0 for x = 1:length(ct_idx[i])]) end
                    ft_cspf[:,:,ri] = cmat * ce_pf
                    ft_p[ri] = p_reg

                    if mode == hx_sda
                        ft_cepc[ri] = ce_tot / wg_sum
                        ft_cpbc[:,ri] = ct_pf ./ ce_tot
                    elseif mode == cat_sda
                        for i = 1:nc, j = 1:nc; ft_cepc[i][j,j,ri] = (i == j ? ct_pf[i] / wg_sum : 1.0) end
                    end
                end
                ft_de[ri] = (sum(de[:, idxs], dims=1) * wg_reg)[1] / wg_sum * p_reg
            end
            if mode == pt_sda; ft.p, ft.cepc, ft.cspf, ft.de = ft_p, ft_cepc, ft_cspf, ft_de
            elseif mode == hx_sda; ft.p, ft.cepc, ft.cspfbc, ft.cpbc, ft.de = ft_p, ft_cepc, ft_cspf, ft_cpbc, ft_de
            elseif mode == cat_sda; ft.p, ft.cepcbc, ft.cspfbc, ft.de = ft_p, ft_cepc, ft_cspf, ft_de
            else println("SDA mode error: ", mode)
            end

            if !haskey(sda_factors, y); sda_factors[y] = Dict{String, factors}() end
            sda_factors[y][n] = ft
            mrio, etab, cmat = [], [], []
        end
        if y != baseYear; t_bp, t_tax, t_sub, v_bp, y_bp = [], [], [], [], [] end
    end
end

function prepareDeltaFactors(target_year, base_year; nation = "", mode = "penta")

    global nat_list, sda_factors, dltByNat, cat_list
    if length(nation) == 0; nats = nat_list else nats = [nation] end
    if mode == "categorized"; nc = length(cat_list) end

    for n in nats
        t_ft, b_ft = sda_factors[target_year][n], sda_factors[base_year][n]

        dltByNat[n] = Dict{Int, Any}()
        dltByNat[n][1] = t_ft.f - b_ft.f
        dltByNat[n][2] = t_ft.l - b_ft.l
        dltByNat[n][3] = t_ft.p - b_ft.p

        if mode == "penta"
            dltByNat[n][4] = t_ft.cepc - b_ft.cepc
            dltByNat[n][5] = t_ft.cspf - b_ft.cspf
        elseif mode == "hexa"
            dltByNat[n][4] = t_ft.cepc - b_ft.cepc
            dltByNat[n][5] = t_ft.cspfbc - b_ft.cspfbc
            dltByNat[n][6] = t_ft.cpbc - b_ft.cpbc
        elseif mode == "categorized"
            for i = 1:nc; dltByNat[n][i+3] = t_ft.cepcbc[i] - b_ft.cepcbc[i] end
            dltByNat[n][nc+4] = t_ft.cspfbc - b_ft.cspfbc
        else println("SDA mode error: ", mode)
        end
    end
end

function calculateDeltaFactors(target_year, base_year, nation, delta_factor, sub_list)

    global sda_factors, dltByNat, cat_list

    yrs = [target_year, base_year]
    fts = sda_factors[target_year][nation], sda_factors[base_year][nation]
    subs = [sub_list[1:delta_factor-1]; 1; sub_list[delta_factor:end]]

    if mode == "penta"
        var = [fts[subs[1]].f, fts[subs[2]].l, fts[subs[3]].p, fts[subs[4]].cepc, fts[subs[5]].cspf]
        var[delta_factor] = dltByNat[nation][delta_factor]
        ie = vec(sum(var[1] .* var[2] * ((var[3] .* var[4])' .* var[5]), dims=1))
    elseif mode == "hexa"
        var = [fts[subs[1]].f, fts[subs[2]].l, fts[subs[3]].p, fts[subs[4]].cepc, fts[subs[5]].cspfbc, fts[subs[6]].cpbc]
        var[delta_factor] = dltByNat[nation][delta_factor]
        nt, nr = size(var[2], 1), size(var[3], 1)
        cspf = zeros(Float64, nt, nr)
        for i = 1:nr; cspf[:,i] = var[5][:,:,i] * var[6][:,i] end
        ie = vec(sum(var[1] .* var[2] * ((var[3] .* var[4])' .* cspf), dims=1))
    elseif mode == "categorized"
        var = [fts[subs[1]].f, fts[subs[2]].l, fts[subs[3]].p]
        for i = 1:nc; push!(var, fts[subs[i+3]].cepcbc[i]) end
        push!(var, fts[subs[nc+4]].cspfbc)
        var[delta_factor] = dltByNat[nation][delta_factor]
        nt, nr, nc = size(var[2], 1), size(var[3], 1), cat_list
        ie = zeros(Float64, nr, nc)
        for i = 1:nr
            ft_ce = zeros(Float64, nt, nr)
            ft_cepcbc = Matrix(1.0I, nc, nc)
            for j = 1:nc; ft_cepcbc *= var[i+3][:,:,i] end
            ft_cepfbc = var[nc+4][:,:,i] * ft_cepcbc
            ie[i,:] = vec(sum(var[1] .* var[2] * (var[3] .* ft_cepfbc), dims=1))
        end
    else println("SDA mode error: ", mode)
    end

    return ie
end

function calculateEmission(year, nation; mode = "penta")
    global sda_factors, nutsByNat, cat_list
    ft = sda_factors[year][nation]
    if mode in ["hexa", "categorized"]; nt, nr, nc = size(ft.cspfbc, 1), length(nutsByNat[year][nation]), length(cat_list) end

    if mode == "penta"
        ie = vec(sum(ft.f .* ft.l * ((ft.p .* ft.cepc)' .* ft.cspf), dims=1))
    elseif mode == "hexa"
        ft_cspf = zeros(Float64, nt, nr)
        for i = 1:nr; ft_cspf[:,i] = ft.cspfbc[:,:,i] * ft.cpbc[:,i] end
        ie = vec(sum(ft.f .* ft.l * ((ft.p .* ft.cepc)' .* ft_cspf), dims=1))
    elseif mode == "categorized"
        ie = zeros(Float64, nr, nc)
        for i = 1:nr
            ft_ce = zeros(Float64, nt, nr)
            ft_cepcbc = Matrix(1.0I, nc, nc)
            for j = 1:nc; ft_cepcbc *= ft.cepcbc[j][:,:,i] end
            ft_cepfbc = ft.cspfbc[:,:,i] * ft_cepcbc
            ie[i,:] = vec(sum(ft.f .* ft.l * (ft.p[i] .* ft_cepfbc), dims=1))
        end
    else println("SDA mode error: ", mode)
    end

    return ie
end

function generateAllCombination(subs_list, n_factor; elements = [0,1])

    subs = Array{Array{Int, 1}, 1}()

    if length(subs_list) == 0; subs = [[e] for e in elements]
    else for sl in subs_list, e in elements; push!(subs, [sl; e]) end
    end

    if length(subs[1]) == n_factor - 1; return subs
    else generateAllCombination(subs, n_factor, elements = elements)
    end
end

function structuralAnalysis(target_year, base_year, nation, n_factor)

    global deltas, nutsByNat, sda_factors, ieByNat, deByNat, popByNat, totExpByNat
    if !haskey(deltas, (target_year, base_year)); deltas[(target_year, base_year)] = Dict{String, Dict{String, Array{Float64, 1}}}() end
    if !haskey(ieByNat, target_year); ieByNat[target_year] = Dict{String, Array{Float64, 1}}() end
    if !haskey(deByNat, target_year); deByNat[target_year] = Dict{String, Array{Float64, 1}}() end
    if !haskey(ieByNat, base_year); ieByNat[base_year] = Dict{String, Array{Float64, 1}}() end
    if !haskey(deByNat, base_year); deByNat[base_year] = Dict{String, Array{Float64, 1}}() end
    if !haskey(popByNat, target_year); popByNat[target_year] = Dict{String, Array{Float64, 1}}() end
    if !haskey(popByNat, base_year); popByNat[base_year] = Dict{String, Array{Float64, 1}}() end
    if !haskey(expPcByNat, target_year); expPcByNat[target_year] = Dict{String, Array{Float64, 1}}() end
    if !haskey(expPcByNat, base_year); expPcByNat[base_year] = Dict{String, Array{Float64, 1}}() end

    deltas[(target_year, base_year)][nation] = Dict{String, Array{Float64, 1}}()

    for nt in nutsByNat[target_year][nation]; deltas[(target_year, base_year)][nation][nt] = zeros(Float64, n_factor) end
    # if nutsByNat[target_year][nation] != nutsByNat[base_year][nation]
    #     println("NUTS lists are not consistent: ", nation, "\t", target_year, " ", nutsByNat[target_year][nation], "\t", base_year, " ", nutsByNat[base_year][nation])
    # end

    nf = n_factor
    nk = nf - 1
    nn = length(nutsByNat[target_year][nation])
    dlt_repo = Array{Array{delta, 1}, 1}()

    wghs = Dict(0:nk .=> [factorial(nk - k) * factorial(k) for k = 0:nk])
    subs_list = generateAllCombination(Array{Int, 1}(), n_factor, elements = [1,2])
    wgh_subs = Array{Tuple{Float64, Array{Int, 1}}, 1}()

    for sl in subs_list; push!(wgh_subs , (wghs[count(x->x==2, sl)], sl)) end

    for i = 1:nf
        tot_wgh, dlts = 0.0, zeros(Float64, nn)
        dlt_list = Array{delta, 1}()
        for (wgh, sl) in wgh_subs
            dlt_vec = vec(calculateDeltaFactors(target_year, base_year, nation, i, sl))
            dlt = delta(i, n_factor, sub_list = sl, weight = wgh, delta_value = dlt_vec)
            dlts .+= wgh .* dlt_vec
            push!(dlt_list, dlt)
            tot_wgh += wgh
        end
        dlts ./= tot_wgh
        for j = 1:nn; deltas[(target_year, base_year)][nation][nutsByNat[target_year][nation][j]][i] = dlts[j]  end
        push!(dlt_repo, dlt_list)
    end

    ieByNat[target_year][nation], ieByNat[base_year][nation] = vec(calculateEmission(target_year,nation)), vec(calculateEmission(base_year,nation))
    deByNat[target_year][nation], deByNat[base_year][nation] = sda_factors[target_year][nation].de, sda_factors[base_year][nation].de
    popByNat[target_year][nation], popByNat[base_year][nation] = sda_factors[target_year][nation].p, sda_factors[base_year][nation].p
    expPcByNat[target_year][nation], expPcByNat[base_year][nation] = sda_factors[target_year][nation].cepc, sda_factors[base_year][nation].cepc

    return dlt_repo
end

function estimateConfidenceIntervals(year, nation = []; iter = 10000, ci_rate = 0.95, resample_size = 0, replacement = false, pop_dens = 0, visible = false)
    # bootstrap method
    # ci_per: confidence interval percentage
    # replacement: [0] sampling with replacement
    # if resample_size: [0] resample_size = sample_size

    global nat_list, nutsByNat, hh_list, pops, pop_list, pop_linked_cd, pops_ds
    global ci_ie, ci_de, in_emiss, di_emiss, ieByNat, deByNat

    if resample_size == 0; replacement = true end
    if isa(year, Number); year = [year] end
    if length(nation) == 0; nats = nat_list else nats = nation end

    for y in year, n in nats
        if visible; println(" ", y," ", n) end

        if !haskey(ci_ie, y); ci_ie[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
        if !haskey(ci_de, y); ci_de[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
        if !haskey(ci_ie[y], n); ci_ie[y][n] = Dict{String, Tuple{Float64, Float64}}() end
        if !haskey(ci_de[y], n); ci_de[y][n] = Dict{String, Tuple{Float64, Float64}}() end
        if !haskey(ieByNat, y); ieByNat[y] = Dict{String, Array{Float64, 1}}() end
        if !haskey(deByNat, y); deByNat[y] = Dict{String, Array{Float64, 1}}() end

        nts, hhs = nutsByNat[y][n], hh_list[y][n]
        ie, de = vec(sum(in_emiss[y][n], dims=1)), vec(sum(di_emiss[y][n], dims=1))
        nh, nnt = length(hhs), length(nts)
        ieByNat[y][n], deByNat[y][n] = zeros(Float64, nnt), zeros(Float64, nnt)

        for ri = 1:nnt
            r = nts[ri]
            r_p = pop_linked_cd[y][r]
            while r_p[end] == '0'; r_p = r_p[1:end-1] end
            p_reg = pop_dens in [1,2,3] ? pops_ds[y][r_p][pop_dens] : pops[y][r_p]

            idxs = filter(x -> households[y][n][hhs[x]].nuts1 == r, 1:nh)
            if pop_dens in [1,2,3]; filter!(x -> households[y][n][hhs[x]].popdens == pop_dens, idxs) end
            wg_reg = [households[y][n][h].weight_nt for h in hh_list[y][n][idxs]]

            nsam = (resample_size == 0 ? length(idxs) : resample_size)
            ie_vals, de_vals = zeros(Float64, iter), zeros(Float64, iter)

            for i = 1:iter
                if replacement; re_idx = [trunc(Int, nsam * rand())+1 for x = 1:nsam]
                else re_idx = sortperm([rand() for x = 1:nsam])
                end
                wg_sum = sum(wg_reg[re_idx] .* [households[y][n][h].size for h in hh_list[y][n][idxs[re_idx]]])

                ie_vals[i] = sum(ie[idxs[re_idx]] .* wg_reg[re_idx]) / wg_sum * p_reg
                de_vals[i] = sum(de[idxs[re_idx]] .* wg_reg[re_idx]) / wg_sum * p_reg
            end
            sort!(ie_vals)
            sort!(de_vals)

            l_idx, u_idx = trunc(Int, (1 - ci_rate) / 2 * iter) + 1, trunc(Int, ((1 - ci_rate) / 2 + ci_rate) * iter) + 1
            ci_ie[y][n][r] = (ie_vals[l_idx], ie_vals[u_idx])
            ci_de[y][n][r] = (de_vals[l_idx], de_vals[u_idx])

            wg_sum = sum(wg_reg .* [households[y][n][h].size for h in hh_list[y][n][idxs]])

            ieByNat[y][n][ri] = sum(ie[idxs] .* wg_reg) / wg_sum * p_reg
            deByNat[y][n][ri] = sum(de[idxs] .* wg_reg) / wg_sum * p_reg
        end
    end
end

function printConfidenceIntervals(year, outputFile, nation = []; pop_dens = 0, ci_rate = 0.95)

    global nat_list, nutsByNat, hh_list, pops, pop_list, pop_linked_cd, pops_ds, ci_ie, ci_de, ieByNat, deByNat, in_emiss, di_emiss
    if isa(year, Number); year = [year] end
    if length(nation) == 0; nats = nat_list else nats = nation end

    dens_label = Dict(0 => "all", 1 => "densely", 2 => "inter", 3 => "sparsely")
    low_lab, upp_lab = string((1 - ci_rate) / 2), string((1 - ci_rate) / 2 + ci_rate)

    f = open(outputFile, "w")
    print(f, "Year\tNation\tNUTS\tDensity\tSamples\t")
    print(f, "Overall_IE\tIE_CI_", low_lab, "\tIE_CI_", upp_lab, "\tOverall_DE\tDE_CI_", low_lab, "\tDE_CI_", upp_lab)
    print(f, "\tIE_per_capita\tDE_per_capita\tPopulation\tTotal_weight")
    println(f)

    for y in year, n in nats
        nts, hhs, nh = nutsByNat[y][n], hh_list[y][n], length(hh_list[y][n])
        ie, de = vec(sum(in_emiss[y][n], dims=1)), vec(sum(di_emiss[y][n], dims=1))

        for ri = 1:length(nts)
            r = nts[ri]
            r_p = pop_linked_cd[y][r]
            while r_p[end] == '0'; r_p = r_p[1:end-1] end
            p_reg = pop_dens in [1,2,3] ? pops_ds[y][r_p][pop_dens] : pops[y][r_p]

            idxs = filter(x -> households[y][n][hhs[x]].nuts1 == r, 1:nh)
            if pop_dens in [1,2,3]; filter!(x -> households[y][n][hhs[x]].popdens == pop_dens, idxs) end

            wg_reg = [households[y][n][h].weight_nt for h in hh_list[y][n][idxs]]
            wg_sum = sum(wg_reg .* [households[y][n][h].size for h in hh_list[y][n][idxs]])

            print(f, y, "\t", n, "\t", r, "\t", dens_label[pop_dens], "\t", length(idxs))
            print(f, "\t", ieByNat[y][n][ri], "\t", ci_ie[y][n][r][1], "\t", ci_ie[y][n][r][2])
            print(f, "\t", deByNat[y][n][ri], "\t", ci_de[y][n][r][1], "\t", ci_de[y][n][r][2])
            print(f, "\t", sum(ie[idxs] .* wg_reg) / wg_sum, "\t", sum(de[idxs] .* wg_reg) / wg_sum, "\t", p_reg, "\t", wg_sum)
            println(f)
        end
    end
    close(f)
end

function printNUTS(year, outputPath)

    global sda_factors, nutsByNat

    f_path = outputPath * string(year) * "/"
    mkpath(f_path)
    nats = sort(collect(keys(sda_factors[year])))
    f = open(f_path * string(year) * "_nuts.txt", "w")
    println(f, "Nation\tNUTS")
    for n in nats; print(f, n); for nt in nutsByNat[year][n]; print(f, "\t", nt) end; println(f) end
    close(f)
end

function printLmatrix(year, outputPath; nation="", base_year=0)

    global sda_factors, l_factor

    f_path = outputPath * string(year) * "/"
    mkpath(f_path)
    if length(nation) == 0; f_name, l = string(year) * "_L_factor.txt", l_factor[year]
    elseif base_year == 0; f_name, l = string(year) * "_" * nation * "_L_factor.txt", sda_factors[year][nation].l
    else f_name, l = string(year) * "_" * nation * "_L_factor_" * string(base_year) * "_price.txt", sda_factors[year][nation].l
    end
    f = open(f_path * f_name, "w")
    for i = 1:size(l,1); print(f, l[i,1]); for j = 2:size(l,2); print(f, "\t", l[i,j]) end; println(f) end
    close(f)
end

function printFactors(outputPath; year=0, nation="")

    global sda_factors, l_factor, nutsByNat
    if year == 0; yrs = sort(collect(keys(sda_factors))) else yrs = [year] end

    for y in yrs
        if length(nation) == 0; nats = sort(collect(keys(sda_factors[y]))) else nats = [nation] end
        f_path = outputPath * string(y) * "/"
        mkpath(f_path)

        for n in nats
            ft = sda_factors[y][n]
            f = open(f_path * string(y) * "_" * n * "_factors.txt", "w")
            println(f, "[f]")
            nt = size(ft.f, 1)
            print(f, ft.f[1]); for i = 2:nt; print(f, "\t", ft.f[i]) end; println(f)
            println(f, "[p]")
            nr = size(ft.p, 1)
            print(f, ft.p[1]); for i = 2:nr; print(f, "\t", ft.p[i]) end; println(f)
            println(f, "[cepc]")
            if nr == size(ft.cepc, 1); print(f, ft.cepc[1]); for i = 2:nr; print(f, "\t", ft.cepc[i]) end; println(f)
            else println(n, ": P matrix and CEPC vector sizes do not match.")
            end
            println(f, "[cspf]")
            if nt == size(ft.cspf,1) && nr == size(ft.cspf,2)
                # for i = 1:nt; print(f, ft.cspf[i,1]); for j = 2:nr; print(f, "\t", ft.cspf[i,j]) end; println(f) end
                for i = 1:nr; print(f, ft.cspf[1,i]); for j = 2:nt; print(f, "\t", ft.cspf[j,i]) end; println(f) end
            else println(n, ": CSPF matrix and F & P vectors sizes do not match.")
            end
            println(f, "[de]")
            if nr == size(ft.de, 1); print(f, ft.de[1]); for i = 2:nr; print(f, "\t", ft.de[i]) end; println(f)
            else println(n, ": P matrix and DE vector sizes do not match.")
            end
            close(f)
        end
    end
end

function readNUTS(year, inputPath)

    global nuts_list, nutsByNat
    nuts_list[year], nutsByNat[year] = Array{String, 1}(), Dict{String, Array{String, 1}}()

    f_path = inputPath * string(year) * "/"
    f = open(f_path * string(year) * "_nuts.txt")
    readline(f)
    for l in eachline(f)
        s = string.(strip.(split(l, "\t")))
        if s[1] in nat_list
            if !haskey(nutsByNat[year], s[1]); nutsByNat[year][s[1]] = Array{String, 1}() end
            append!(nutsByNat[year][s[1]], filter(x->!(x in nutsByNat[year][s[1]]), s[2:end]))
            append!(nuts_list[year], filter(x->!(x in nuts_list[year]), nutsByNat[year][s[1]]))
        end
    end
    close(f)
end

function readLmatrix(year, inputFile)

    f = open(inputFile)
    l = parse.(Float64, string.(split(strip(readline(f)), "\t")))
    nt = length(l)
    l_mat = zeros(Float64, nt, nt)
    l_mat[1,:] = l
    for i = 2:nt; l_mat[i,:] = parse.(Float64, string.(split(strip(readline(f)), "\t"))) end
    close(f)

    return l_mat
end

function readFactors(year, base_year, inputPath; nation = "", visible = false)

    global sda_factors, nat_list, l_factor
    sda_factors[year] = Dict{String, factors}()
    if length(nation) == 0; nats = nat_list else nats = [nation] end

    f_path = inputPath * string(year) * "/"

    if year == base_year; l_factor[year] = readLmatrix(year, f_path * string(year) * "_L_factor.txt") end
    for n in nat_list
        if visible; print(", ", n) end
        nt, nr, f_mat, p_mat, cepc_mat, cspf_mat, de_mat = 0, 0, [], [], [], [], []
        f = open(f_path * string(year) * "_" * n * "_factors.txt")
        if string(strip(readline(f))) == "[f]"
            f_mat = parse.(Float64, string.(split(strip(readline(f)), "\t")))
            nt = length(f_mat)
        else println("F-factor data error.")
        end
        if string(strip(readline(f))) == "[p]"
            p_mat = parse.(Float64, string.(split(strip(readline(f)), "\t")))
            nr = length(p_mat)
        else println("P-factor data error.")
        end
        if string(strip(readline(f))) == "[cepc]"
            cepc_mat = parse.(Float64, string.(split(strip(readline(f)), "\t")))
            if length(cepc_mat) != nr; println("CEPC size ", length(cepc_mat)," does not match ", nr) end
        else println("CEPC-factor data error.")
        end
        if string(strip(readline(f))) == "[cspf]"
            cspf_mat = zeros(Float64, nt, nr)
            # for i = 1:nt; cspf_mat[i,:] = parse.(Float64, string.(split(strip(readline(f)), "\t"))) end
            for i = 1:nr; cspf_mat[:,i] = parse.(Float64, string.(split(strip(readline(f)), "\t"))) end
            if size(cspf_mat, 1) != nt || size(cspf_mat, 2) != nr; println("CSPF size ", size(cspf_mat)," does not match ", nt, " x ", nr) end
        else println("CSPF-factor data error.")
        end
        if string(strip(readline(f))) == "[de]"
            de_mat = parse.(Float64, string.(split(strip(readline(f)), "\t")))
            if length(de_mat) != nr; println("DE size ", length(de_mat)," does not match ", nr) end
        else println("DE-factor data error.")
        end
        close(f)
        if year == base_year; l_mat = l_factor[year]
        else l_mat = readLmatrix(year, f_path * string(year) * "_" * n * "_L_factor_" * string(base_year) * "_price.txt")
        end
        if size(l_mat, 1) != nt; println("L-matrix and F-factor size do not match: ", size(l_mat), "\t", size(f_mat)) end

        sda_factors[year][n] = factors(f = f_mat, l = l_mat, p = p_mat, exp_pc = cepc_mat, exp_prof = cspf_mat, de = de_mat)
    end
end

function printDelta(outputFile; cf_print = true, st_print = true, mode = "penta")

    global nutsByNat, deltas, ieByNat, deByNat, popByNat, expPcByNat, cat_list

    vs = getValueSeparator(outputFile)
    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    print(f, "Target_year", vs, "Base_year", vs, "Nation", vs, "NUTS")
    print(f, vs, "f", vs, "L", vs, "p")
    if mode == "penta"; print(f, vs, "exp_per_cap", vs, "exp_profile")
    elseif mode == "hexa"; print(f, vs, "exp_per_cap", vs, "exp_profile", vs, "exp_cat")
    elseif mode == "categorized"
        for c in cat_list; print(f, vs, c, "_exp_per_cap") end
        print(f, vs, "exp_profile")
    else println("SDA mode error: ", mode)
    end
    print(f, vs, "de", vs, "total_delta")
    if cf_print; print(f, vs, "Measured_delta", vs, "Target_year_IE", vs, "Base_year_IE", vs, "Target_year_DE", vs, "Base_year_DE") end
    if st_print; print(f, vs, "Population_target", vs, "Population_base", vs, "Total_exp_target", vs, "Total_exp_base") end
    println(f)
    for yrs in sort(collect(keys(deltas))), n in sort(collect(keys(deltas[yrs])))
        if cf_print
            t_ie, t_de = ieByNat[yrs[1]][n], deByNat[yrs[1]][n]
            b_ie, b_de = ieByNat[yrs[2]][n], deByNat[yrs[2]][n]
            t_p, t_epc = popByNat[yrs[1]][n], expPcByNat[yrs[1]][n]
            b_p, b_epc = popByNat[yrs[2]][n], expPcByNat[yrs[2]][n]
        end
        nn = length(nutsByNat[yrs[1]][n])
        for i = 1:nn
            nt = nutsByNat[yrs[1]][n][i]
            print(f, yrs[1], vs, yrs[2], vs, n, vs, nt)
            for d in deltas[yrs][n][nt]; print(f, vs, d) end
            print(f, vs, t_de[i] - b_de[i])
            print(f, vs, sum(deltas[yrs][n][nt]))
            if cf_print; print(f, vs, t_ie[i] - b_ie[i], vs, t_ie[i], vs, b_ie[i], vs, t_de[i], vs, b_de[i]) end
            if st_print; print(f, vs, t_p[i], vs, b_p[i], vs, t_epc[i], vs, b_epc[i]) end
            println(f)
        end

    end
    close(f)
end

function clearFactors(; year = 0, nation = "")

    global sda_factors, households, exp_table, mrio_tabs_conv, conc_mat_wgh, di_emiss

    if year == 0; yrs = sort(collect(keys(sda_factors))) else yrs = [year] end
    for y in yrs
        if length(nation) == 0; nats = sort(collect(keys(sda_factors[y]))) else nats = [nation] end
        for n in nats
            households[y][n] = Dict{String, mdr.household}()
            exp_table[y][n] = Array{Float64, 2}(undef, 0, 0)
            conc_mat_wgh[y][n] = Array{Float64, 2}(undef, 0, 0)
            sda_factors[y][n] = factors()
            if haskey(di_emiss, y); di_emiss[y][n] = Array{Float64, 2}(undef, 0, 0)  end
            if haskey(mrio_tabs_conv, y); mrio_tabs_conv[y][n] = ee.tables() end
        end
    end
end

function testSDA(year, output)

    global sda_factors, nat_list, nutsBuNat
    cf = Dict{String, Array{Float64, 2}}()

    for n in nat_list, nt in nutsBuNat[year][n]
        ft = sda_factors[year][n]
        cf[n] = sum(ft.f .* ft.l * ((ft.p .* ft.cepc)' .* ft.cspf), dims=1)' + ft.de
    end

    for n in nat_list
        println(n ,"\t", nutsBuNat[year][n], "\t", cf[n], "\t", cf[n]./sda_factors[year][n].p)
        println(n, "\t", sda_factors[year][n].de, "\t", sda_factors[year][n].de./sda_factors[year][n].p)
        println(n, "\t", sda_factors[year][n].p)
        println(n, "\t", sda_factors[year][n].cepc)
    end

end

end
