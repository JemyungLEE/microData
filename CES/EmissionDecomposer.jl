# SPDX-FileCopyrightText: © 2023 Jemyung Lee <jemyung81@gmail.com>
# SPDX-License-Identifier: GPL-3.0

module EmissionDecomposer

# Developed date: 25. Apr. 2023
# Last modified date: 4. Jul. 2023
# Subject: Decompose household carbon footprints
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

    cspfbc::Array{Array{Float64, 2}, 1} # Regional household expenditure profile by expenditure category: {region, {Eora t-index, category}}
    cpbc::Array{Float64, 2}             # Consumption profile by expenditure category: {category, region}

    cepcbc::Array{Array{Array{Float64, 2}, 1}, 1}   # Consumption expenditures per capita by expenditure category: {category, {region, {category, category}}}

    function factors(nr=0, nt=0; nc=0, factor_by_cat = false, f=zeros(0), l=zeros(0,0), p=zeros(0), exp_pc=zeros(0), exp_prof=zeros(0,0), de=zeros(0), int_share=Array{Array{Float64, 2}, 1}(), dom_share=zeros(0,0), exp_pc_cat=Array{Array{Array{Float64, 3}, 1}, 1}())
        if nr > 0 && nt > 0 && length(f) == length(l) == length(p) == 0
            if nc == length(exp_pc) == length(exp_prof) == length(de) == 0
                new(zeros(nt), zeros(nt, nt), zeros(nr), zeros(nr), zeros(nt, nr), zeros(nr), int_share, zeros(0,0), exp_pc_cat)
            elseif !factor_by_cat && nc > 0 && length(exp_pc) == length(de) == length(int_share) == length(dom_share) == 0
                new(zeros(nt), zeros(nt, nt), zeros(nr), zeros(nr), zeros(0, 0), zeros(nr), [zeros(nt, nc) for i=1:nr], zeros(nc, nr), exp_pc_cat)
            elseif factor_by_cat && nc > 0 && length(de) == length(int_share) == length(exp_pc_cat) == 0
                new(zeros(nt), zeros(nt, nt), zeros(nr), zeros(0), zeros(0, 0), zeros(nr), [zeros(nt, nc) for i=1:nr], zeros(0, 0), [[zeros(nc, nc) for j=1:nr] for i=1:nc])
            end
        else new(f, l, p, exp_pc, exp_prof, de, int_share, dom_share, exp_pc_cat)
        end
    end
end

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

yr_list = Array{Int, 1}()            # year list: {YYYY}
nat_list = Array{String, 1}()        # nation list: {A3}
cat_list = Array{String, 1}()        # category list
pr_unts = Dict("day" => 1, "week" => 7, "month" => 30, "year" => 365, "annual" => 365, "monthly" => 30, "weekly" => 7, "daily" => 1)

sc_list = Dict{Int, Dict{String, Array{String, 1}}}()               # commodity code list: {year, {nation A3, {code}}}
hh_list = Dict{Int, Dict{String, Array{String, 1}}}()               # Household ID: {year, {nation, {hhid}}}
gr_list = Dict{Int, Dict{String, Array{String, 1}}}()               # group (survey type) list: {year, {nation, {group}}}
households = Dict{Int, Dict{String, Dict{String, mdr.household}}}() # household dict: {year, {nation, {hhid, household}}}
sectors = Dict{Int, Dict{String, Dict{String, mdr.commodity}}}()    # expenditure sector: {year, {nation A3, {code, commodity (of mdr)}}}
exp_table = Dict{Int, Dict{String, Array{Float64, 2}}}()            # household expenditure table: {year, {nation, {hhid, category}}}
sc_cat = Dict{Int, Dict{String, Dict{String, String}}}()            # category dictionary: {year, {nation, {sector code, category}}}

reg_list = Dict{Int, Dict{String, Array{String, 1}}}()      # GIS region list: {year, {nation, {gis_id}}}
regions = Dict{Int, Dict{String, Dict{String, String}}}()   # region code-name: {year, {nation A3, {code, region}}}
reg_linked = Dict{Int, Dict{String, Dict{String, String}}}()# CES/HBS region corresponding GIS region code: {year, {nation, {CES/HBS region, GIS code}}}

hh_ie = Dict{Int, Dict{String, Array{Float64, 1}}}()        # indirect carbon emission: {year, {nation, {hhid}}}
hh_de = Dict{Int, Dict{String, Array{Float64, 1}}}()        # direct carbon emission: {year, {nation, {hhid}}}
hh_cf = Dict{Int, Dict{String, Array{Float64, 1}}}()        # categorized carbon footprint by household: {year, {nation, {hhid (AA_HHID)}}}
pops = Dict{Int, Dict{String, Float64}}()                   # Population: {year, {region, population}}

cpis = Dict{Int, Dict{String, Dict{String, Float64}}}()     # Consumption price indexes: {year, {nation, {COICOP_category, CPI}}}
scl_rate = Dict{Int, Dict{String, Dict{String, Float64}}}() # CPI scaling rate: {year, {nation, {HBS code, rate}}}
conc_mat_wgh = Dict{Int, Dict{String, Array{Float64, 2}}}() # Weighted concordance matrix: {year, {nation, {Eora sectors, Nation sectors}}}
conc_mat = Dict{Int, Dict{String, Array{Float64, 2}}}()     # Concordance matrix: {year, {nation, {Eora sectors, Nation sectors}}}
mrio_idxs = Array{ee.idx, 1}()                              # index T
mrio_tabs = Dict{Int, ee.tables}()                          # MRIO tables: {Year, MRIO tables (t, v ,y , q)}
mrio_tabs_conv = Dict{Int, Dict{String, ee.tables}}()       # Base-year price converted MRIO tables: {Year, {natoin, MRIO tables (t, v ,y , q)}}

l_factor = Dict{Int, Array{Float64, 2}}()                   # Leontief matrix: {year, {Eora t-index, Eora t-index}}
f_factor = Dict{Int, Array{Float64, 1}}()                   # Intensity in Leontief matrix: {year, {Eora t-index}}
sda_factors = Dict{Int, Dict{String, factors}}()            # SDA factors: {year, {nation, factors}}
dltByNat = Dict{String, Dict{Int, Any}}()                   # delta by factor, by nation: {nation, {factor, {target_year - base_year}}}
deltas = Dict{Tuple{Int,Int}, Dict{String, Dict{String, Array{Float64, 1}}}}()  # Deltas of elements: {(target_year, base_year), {nation, {region, {factor}}}}
ieByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()      # indirect CF by nation, region: {year, {nation, {region}}}
deByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()      # direct CF by nation, region: {year, {nation, {region}}}
cfByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()      # total CF by nation, region: {year, {nation, {region}}}
cfByCat = Dict{Int, Dict{String, Array{Float64, 2}}}()      # categozied carbon footprint by region: {year, {nation, {region, category}}}
popByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()     # population by nation, region: {year, {nation, {region}}}
expPcByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()   # total expenditures by nation, region: {year, {nation, {region}}}

ci_ie = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}() # confidence intervals of indirect emission: {year, {nation, {region, {lower, upper}}}}
ci_de = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}() # confidence intervals of direct emission: {year, {nation, {region, {lower, upper}}}}
ci_cf = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}() # confidence intervals of CF: {year, {nation, {region, {lower, upper}}}}
ci_cfpc = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}() # confidence intervals of per capita CF: {year, {nation, {region, {lower, upper}}}}
ci_sda = Dict{Tuple{Int,Int}, Dict{String, Dict{String, Array{Tuple{Float64, Float64}, 1}}}}()   # confidence intervals of SDA factors: {(target_year, base_year), {nation, {region, {factor (lower, upper)}}}}

samples_gr = Dict{String, Dict{Int, []}}()

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

function importData(; hh_data::Module, mrio_data::Module, cat_data::Module, conc_data::Module, cat_filter = true, region = "district", exp_mode = "gis")
    # exp_mode: "gis" use GIS map's region classification, "ces" or "hbs" use CES/HBS region classification

    global yr_list = cat_data.yr_list
    global hh_list, households, exp_table, sectors = hh_data.hh_list, hh_data.households, hh_data.expMatrix, hh_data.sectors
    global gr_list = hh_data.gr_list
    global scl_rate, cpis = hh_data.scl_rate, hh_data.cpis
    global mrio_idxs, mrio_tabs, sc_list = mrio_data.ti, mrio_data.mTables, mrio_data.sc_list
    global conc_mat, conc_mat_wgh = mrio_data.concMat, mrio_data.concMatWgh
    global cat_list, sc_cat = cat_data.cat_list, cat_data.sc_cat

    if lowercase(exp_mode) == "gis"
        global reg_list, regions, reg_linked, cfByCat = cat_data.gisRegList, cat_data.gisRegID, cat_data.gisRegDist, cat_data.cfRegGIS
        for y in yr_list
            global pops[y] = Dict{String, Array{Float64, 1}}()
            for n in collect(keys(cat_data.gisPop[y])); pops[y][n] = Dict(reg_list[y][n] .=> cat_data.gisPop[y][n]) end
        end
    elseif lowercase(exp_mode) in ["ces", "hbs", "ces/hbs", "hbs/ces"]
        global regions, pops, cfByCat = cat_data.regions, cat_data.pops, cat_data.cfReg
        if region == "district"; global reg_list = cat_data.dist_list
        elseif region == "province"; global reg_list = cat_data.prov_list
        end
    end

    if cat_filter; filter!(x -> !(lowercase(x) in ["total", "all"]), cat_list) end

    hh_ie = Dict{Int, Dict{String, Array{Float64, 1}}}()        # indirect carbon emission: {year, {nation, {hhid}}}
    hh_de = Dict{Int, Dict{String, Array{Float64, 1}}}()        # direct carbon emission: {year, {nation, {hhid}}}
    hh_cf = Dict{Int, Dict{String, Array{Float64, 1}}}()        # categorized carbon footprint by household: {year, {nation, {hhid (AA_HHID)}}}


    for y in yr_list
        hh_ie[y], hh_de[y], hh_cf[y] = Dict{String, Array{Float64,1}}(), Dict{String, Array{Float64,1}}(), Dict{String, Array{Float64,1}}()
        for n in collect(keys(cat_data.cfHHs[y]))
            hh_ie[y][n] = vec(sum(cat_data.indirectCE[y][n], dims=1))
            hh_de[y][n] = vec(sum(cat_data.directCE[y][n], dims=1))
            hh_cf[y][n] = cat_data.cfHHs[y][n][:, end]
        end
    end
end

function storeConcMat(year, nation, concMat; conc_mat_nw = [])

    global conc_mat_wgh, conc_mat

    if !haskey(conc_mat_wgh, year); conc_mat_wgh[year] = Dict{String, Array{Float64, 2}}() end
    conc_mat_wgh[year][nation] = concMat
    if length(conc_mat_nw) > 0
        if !haskey(conc_mat, year); conc_mat[year] = Dict{String, Array{Float64, 2}}() end
        conc_mat[year] = conc_mat_nw
    end

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
    avg_scl = cpis[base_year][nation][total_cp] / cpis[year][nation][total_cp]
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

function calculateLeontief(mrio_table)

    tb = mrio_table
    nt = size(tb.t, 1)

    # x = sum(tb.t, dims = 1) +  sum(tb.v, dims = 1)      # calculate X

    x = vec(sum(tb.t, dims = 2) + sum(tb.y, dims = 2))       # calculate X
    lt = Matrix{Float64}(I, nt, nt)                     # calculate Leontief matrix
    for i = 1:nt; for j = 1:nt; lt[i,j] -= tb.t[i,j] / x[j] end end
    lti = inv(lt)

    return lti, x   # Leontief inverse matrix
end

function calculateIntensity(mrio_table; x = [])

    tb = mrio_table
    nt = size(tb.t, 1)

    # x = sum(tb.t, dims = 1) +  sum(tb.v, dims = 1)                  # calculate X

    if length(x) == 0; x = vec(sum(tb.t, dims = 2) + sum(tb.y, dims = 2)) end    # calculate X
    f = [x[i] > 0 ? sum(tb.q, dims = 1)[i] / x[i] : 0.0 for i=1:nt] # calculate EA

    return f    # emission factor (intensity)
end

function decomposeFactors(year, baseYear, nation = "", mrioPath = ""; mode="penta", region = "district", exp_mode = "gis", group = false)
    # group: [true] apply population weightfor each household (survey) group

    # 5-factors: f * L * p * tot_ce_pc * [con * hbs_profile] + DE
    # 6-factors: f * L * p * tot_ce_pc * [con * hbs_profile_by_category * hbs_proportion_by_category] + DE

    global mrio_tabs, mrio_tabs_conv, conc_mat_wgh, sda_factors, hh_de, cat_list, sc_cat, sc_list, gr_list
    global l_factor, f_factor
    global nat_list, reg_list, hh_list, pops, reg_linked

    pt_sda, hx_sda, cat_sda = "penta", "hexa", "categorized"

    if isa(year, Number); year = [year] end
    nats = (length(nation) == 0 ? nat_list : (isa(nation, String) ? [nation] : nation))
    nc = length(cat_list)
    rc_mode = (lowercase(exp_mode) == "gis")

    for y in year
        if y != baseYear; t_bp, t_tax, t_sub, v_bp, y_bp = setMrioTables(y, mrioPath) end
        if !haskey(ieByNat, y); ieByNat[y] = Dict{String, Array{Float64, 1}}() end
        if !haskey(deByNat, y); deByNat[y] = Dict{String, Array{Float64, 1}}() end
        if !haskey(popByNat, y); popByNat[y] = Dict{String, Array{Float64, 1}}() end
        if !haskey(expPcByNat, y); expPcByNat[y] = Dict{String, Array{Float64, 1}}() end
        for n in nats
            etab, cmat = exp_table[y][n], conc_mat_wgh[y][n]
            hhs, hl, de, rl, sl = households[y][n], hh_list[y][n], hh_de[y][n], reg_list[y][n], sc_list[y][n]
            nh, ns = length(hl), length(sl)
            if group; gl = gr_list[y][n]; ng = length(gl) end
            if rc_mode; r_lnk = reg_linked[y][n] end

            ft = factors()
            if y == baseYear
                if !haskey(l_factor, y)
                    l_factor[y], ft_x = calculateLeontief(mrio_tabs[y])
                    f_factor[y] = calculateIntensity(mrio, ft_x)
                end
                mrio, ft.l, ft.f = mrio_tabs[y], l_factor[y], f_factor[y]
            else
                convertTable(y, n, baseYear, mrioPath, t_bp = t_bp, t_tax = t_tax, t_sub = t_sub, v_bp = v_bp, y_bp = y_bp)
                mrio = mrio_tabs_conv[y][n]
                ft.l, ft_x = calculateLeontief(mrio)
                ft.f = calculateIntensity(mrio, ft_x)
            end
            nr, nt = length(rl), size(mrio.t, 1)
            if mode in [hx_sda, cat_sda]
                sct = sc_cat[y][n]
                ct_idx = [findall(x -> haskey(sct,x) && sct[x] == c, sl) for c in cat_list]
            end
            ft_p, ft_de = zeros(nr), zeros(nr)
            if mode == pt_sda; ft_cepc, ft_cspf = zeros(nr), zeros(nt, nr)
            elseif mode == hx_sda; ft_cepc, ft_cspfbc, ft_cpbc = zeros(nr), [zeros(nt, nc) for i=1:nr], zeros(nc, nr)
            elseif mode == cat_sda; ft_cepc, ft_cepcbc, ft_cspfbc = zeros(nr), [[zeros(nc, nc) for j=1:nc] for i=1:nr], [zeros(nt, nc) for i=1:nr]
            end

            hhs_wgs = [hhs[h].popwgh for h in hl]
            hhs_siz = [hhs[h].size for h in hl]
            hhs_wsz = hhs_wgs .* hhs_siz

            for ri = 1:nr
                r = rl[ri]
                p_reg = pops[y][n][r]

                if region == "district"
                    idxs = (rc_mode ? filter(x -> r_lnk[hhs[hl[x]].district] == r, 1:nh) : filter(x -> hhs[hl[x]].district == r, 1:nh))
                elseif region == "province"
                    idxs = (rc_mode ? filter(x -> r_lnk[hhs[hl[x]].province] == r, 1:nh) : filter(x -> hhs[hl[x]].province == r, 1:nh))
                end

                if !group
                    wg_reg = hhs_wgs[idxs]
                    wg_sum = sum(hhs_wsz[idxs])
                    etb_wg = wg_reg .* etab[idxs, :]
                    ft_de[ri] = (de[idxs] * wg_reg)[1] / wg_sum * p_reg
                elseif group
                    idxs_gr = [filter(x -> hhs[hl[x]].group == g, idxs) for g in gl]
                    wg_reg_gr = [hhs_wgs[idx_ls] for idx_ls in idxs_gr]
                    etb_wg_gr = [wg_reg_gr[gi] .* etab[idxs_gr[gi], :] for gi = 1:ng]
                    wg_sum_gr = [sum(hhs_wsz[idx_ls]) for idx_ls in idxs_gr]
                    ft_de[ri] = sum([de[idxs_gr[gi]] * wg_reg_gr[gi])[1] / wg_sum_gr[gi] for gi = 1:ng]) * p_reg
                end
                ft_p[ri] = p_reg

                if mode == pt_sda
                    if !group
                        et_sum = sum(etb_wg, dims=1)
                        ce_tot = sum(et_sum)
                        ce_pf = et_sum ./ ce_tot
                        ft_cepc[ri] = ce_tot / wg_sum
                    elseif group
                        et_sum_gr = [sum(etb_wg_gr[gi], dims=1) for gi = 1:ng]
                        ce_tot_gr = [sum(et_sum_gr[gi]) for gi = 1:ng]
                        ce_pf = sum(et_sum_gr ./ wg_sum_gr)
                        ce_tot = sum(ce_pf)
                        ce_pf ./= ce_tot
                        ft_cepc[ri] = sum(ce_tot_gr ./ wg_sum_gr)
                    end
                    ft_cspf[:,ri] = cmat * ce_pf'
                elseif mode in [hx_sda, cat_sda]
                    if !group
                        et_sum = vec(sum(etb_wg, dims=1))
                        ct_pf = [sum(et_sum[ci]) for ci in ct_idx]
                        ce_tot = sum(ct_pf)
                        ce_pf = zeros(Float64, ns, nc)
                        for i = 1:nc; if ct_pf[i] > 0; ce_pf[ct_idx[i], i] = et_sum[ct_idx[i]] ./ ct_pf[i] end end
                        ft_cepc[ri] = ce_tot / wg_sum
                        if mode == hx_sda; ft_cpbc[:,ri] = ct_pf ./ ce_tot
                        elseif mode == cat_sda; for i = 1:nc, j = 1:nc; ft_cepcbc[ri][i][j,j] = (i == j ? ct_pf[i] / wg_sum : 1.0) end
                        end
                    elseif group
                        et_sum_gr = [vec(sum(etb_wg_gr[gi], dims=1)) for gi = 1:ng]
                        ct_pf_gr = [[sum(et_sum_gr[gi][ci]) for ci in ct_idx] for gi = 1:ng]
                        ce_tot_gr = [sum(ct_pf_gr[gi]) for gi = 1:ng]
                        ce_pf = zeros(Float64, ns, nc)
                        for i = 1:nc
                            for gi = 1:ng
                                if ct_pf_gr[gi][i] > 0; ce_pf[ct_idx[i], i] += et_sum_gr[gi][ct_idx[i]] ./ wg_sum_gr[gi] end
                            end
                            ce_pf_tot = sum(ce_pf[ct_idx[i], i])
                            ce_pf[ct_idx[i], i] ./= ce_pf_tot
                        end
                        ft_cepc[ri] = sum(ce_tot_gr ./ wg_sum_gr)

                        if mode == hx_sda
                            ft_cpbc[:,ri] = sum(ct_pf_gr ./ wg_sum_gr)
                            ct_pf_tot = sum(ft_cpbc[:,ri])
                            ft_cpbc[:,ri] ./= ct_pf_tot
                        elseif mode == cat_sda
                            for i = 1:nc, j = 1:nc
                                ft_cepcbc[ri][i][j,j] = (i == j ? sum([ct_pf_gr[gi][i] / wg_sum_gr[gi] for gi = 1:ng]) : 1.0)
                            end
                        end
                    end
                    ft_cspfbc[ri] = cmat * ce_pf
                end
            end
            if mode == pt_sda; ft.p, ft.cepc, ft.cspf, ft.de = ft_p, ft_cepc, ft_cspf, ft_de
            elseif mode == hx_sda; ft.p, ft.cepc, ft.cspfbc, ft.cpbc, ft.de = ft_p, ft_cepc, ft_cspfbc, ft_cpbc, ft_de
            elseif mode == cat_sda; ft.p, ft.cepc, ft.cepcbc, ft.cspfbc, ft.de = ft_p, ft_cepc, ft_cepcbc, ft_cspfbc, ft_de
            else println("SDA mode error: ", mode)
            end

            if !haskey(sda_factors, y); sda_factors[y] = Dict{String, factors}() end
            sda_factors[y][n] = ft
            ieByNat[y][n] = vec(calculateEmission(y, n, mode = mode))
            deByNat[y][n], popByNat[y][n], expPcByNat[y][n] = ft.de, ft.p, ft.cepc

            mrio, etab, cmat = Float64[;;], Float64[;;], Float64[;;]
        end
        if y != baseYear; t_bp, t_tax, t_sub, v_bp, y_bp = Float64[;;], Float64[;;], Float64[;;], Float64[;;], Float64[;;] end
    end
end

function decomposeFactorsByGroup(year, baseYear, nation = "", mrioPath = ""; mode="penta", region = "district", bndr_mode = "percap", exp_mode = "gis",
                                cf_intv = [], inc_intv = [], hpos_cf = [], hpos_inc = [], cf_bndr = [], inc_bndr = [], group = false)
    # group: [true] apply population weightfor each household (survey) group

    # mode = penta: f * L * p * tot_ce_pc * [con * hbs_profile] + DE
    # mode = hexa:  f * L * p * tot_ce_pc * [con * hbs_profile_by_category * hbs_proportion_by_category] + DE

    # grouping = cf_intv: CF per capita intervals (stacked proportion)
    # grouping = inc_intv: income per capita intervals (stacked proportion)
    # grouping = cf_bndr: CF per capita boundaries (absolute limits)
    # grouping = inc_bndr: income per capita boundaries (absolute limits)
    # bndr_mode: "percap" per capita income or CF for boundary vales, "hhs" household income or CF

    global mrio_tabs, mrio_tabs_conv, conc_mat_wgh, sda_factors, hh_de, cat_list, sc_cat, sc_list, gr_list
    global l_factor, f_factor
    global nat_list, reg_list, hh_list, pops, hh_cf, reg_linked

    pt_sda, hx_sda, cat_sda = "penta", "hexa", "categorized"

    if isa(year, Number); year = [year] end
    nats = (length(nation) == 0 ? nat_list : (isa(nation, String) ? [nation] : nation))
    nc = length(cat_list)
    rc_mode = (lowercase(exp_mode) == "gis")

    n_cf, n_inc = length(cf_intv), length(inc_intv)
    n_cfb, n_incb = length(cf_bndr), length(inc_bndr)
    n_gr = 1 + n_cf + n_inc + n_cfb + n_incb

    if n_cf > 0
        cf_intv_lb = ["_bottom_" * string(ceil(Int, cf_intv[1] * 100)) * "%"]
        append!(cf_intv_lb, ["_middle_" * string(ceil(Int, (cf_intv[i] - cf_intv[i-1]) * 100)) * "%" for i = 2:n_cf-1])
        append!(cf_intv_lb, ["_top_" * string(ceil(Int, (cf_intv[end] - cf_intv[end-1]) * 100)) * "%"])
    end
    if n_inc > 0
        inc_intv_lb = ["_bottom_" * string(ceil(Int, inc_intv[1] * 100)) * "%"]
        append!(inc_intv_lb, ["_middle_" * string(ceil(Int, (inc_intv[i] - inc_intv[i-1]) * 100)) * "%" for i = 2:n_inc-1])
        append!(inc_intv_lb, ["_top_" * string(ceil(Int, (inc_intv[end] - inc_intv[end-1]) * 100)) * "%"])
    end
    if n_cfb > 0
        cf_bndr_lb = [[""];["_" * string(round(cf_bndr[i], digits = 2)) * "≤" for i = 2:n_cfb]]
        for i = 1:n_cfb-1; cf_bndr_lb[i] *= "_<" * string(round(cf_bndr[i+1], digits = 2)) end
    end
    if n_incb > 0
        inc_bndr_lb = [[""];["_" * string(round(inc_bndr[i]/1000, digits = 2)) * "≤" for i = 2:n_incb]]
        for i = 1:n_incb-1; inc_bndr_lb[i] *= "_<" * string(round(inc_bndr[i+1]/1000, digits = 2)) end
    end

    for y in year
        if y != baseYear; t_bp, t_tax, t_sub, v_bp, y_bp = setMrioTables(y, mrioPath) end
        if !haskey(ieByNat, y); ieByNat[y] = Dict{String, Array{Float64, 1}}() end
        if !haskey(deByNat, y); deByNat[y] = Dict{String, Array{Float64, 1}}() end
        if !haskey(popByNat, y); popByNat[y] = Dict{String, Array{Float64, 1}}() end
        if !haskey(expPcByNat, y); expPcByNat[y] = Dict{String, Array{Float64, 1}}() end
        for n in nats
            etab, cmat = exp_table[y][n], conc_mat_wgh[y][n]
            hhs, hl, de, rl, sl = households[y][n], hh_list[y][n], hh_de[y][n], reg_list[y][n][:], sc_list[y][n]
            nh, ns = length(hl), length(sl)
            if group; gl = gr_list[y][n]; ng = length(gl) end
            if rc_mode; r_lnk = reg_linked[y][n] end

            ft = factors()
            if y == baseYear
                if !haskey(l_factor, y)
                    l_factor[y], ft_x = calculateLeontief(mrio_tabs[y])
                    f_factor[y] = calculateIntensity(mrio, ft_x)
                end
                mrio, ft.l, ft.f = mrio_tabs[y], l_factor[y], f_factor[y]
            else
                convertTable(y, n, baseYear, mrioPath, t_bp = t_bp, t_tax = t_tax, t_sub = t_sub, v_bp = v_bp, y_bp = y_bp)
                mrio = mrio_tabs_conv[y][n]
                ft.l, ft_x = calculateLeontief(mrio)
                ft.f = calculateIntensity(mrio, ft_x)
            end
            nr, nt = length(rl), size(mrio.t, 1)
            n_tg = nr * n_gr
            append!(reg_list[y][n], [r * "_CF" * cf_intv_lb[gi] for r in rl, gi = 1:n_cf])
            append!(reg_list[y][n], [r * "_inc" * inc_intv_lb[gi] for r in rl, gi = 1:n_inc])
            append!(reg_list[y][n], [r * "_CF" * cf_bndr_lb[gi] for r in rl, gi = 1:n_cfb])
            append!(reg_list[y][n], [r * "_inc" * inc_bndr_lb[gi] for r in rl, gi = 1:n_incb])

            if mode in [hx_sda, cat_sda]
                sct = sc_cat[y][n]
                ct_idx = [findall(x -> haskey(sct,x) && sct[x] == c, sl) for c in cat_list]
            end
            ft_p, ft_de = zeros(n_tg), zeros(n_tg)
            if mode == pt_sda; ft_cepc, ft_cspf = zeros(n_tg), zeros(nt, n_tg)
            elseif mode == hx_sda; ft_cepc, ft_cspfbc, ft_cpbc = zeros(n_tg), [zeros(nt, nc) for i=1:n_tg], zeros(nc, n_tg)
            elseif mode == cat_sda; ft_cepc, ft_cepcbc, ft_cspfbc = zeros(n_tg), [[zeros(nc, nc) for j=1:nc] for i=1:n_tg], [zeros(nt, nc) for i=1:n_tg]
            end

            intv_dataset, bndr_dataset = [], []
            if n_cf > 0; push!(intv_dataset, (n_cf, hpos_cf[y][n], cf_intv)) end
            if n_inc > 0; push!(intv_dataset, (n_inc, hpos_inc[y][n], inc_intv)) end
            if n_cfb > 0; push!(bndr_dataset, (n_cfb, Dict(hl .=> hh_cf[y][n][:]), cf_bndr)) end
            if n_incb > 0; push!(bndr_dataset, (n_incb, Dict(h => hhs[h].totinc for h in hl), inc_bndr)) end

            hhs_wgs = [hhs[h].popwgh for h in hl]
            hhs_siz = [hhs[h].size for h in hl]
            hhs_wsz = hhs_wgs .* hhs_siz

            for ri = 1:nr
                r = rl[ri]

                if region == "district"
                    r_idxs = (rc_mode ? filter(x -> r_lnk[hhs[hl[x]].district] == r, 1:nh) : filter(x -> hhs[hl[x]].district == r, 1:nh))
                elseif region == "province"
                    r_idxs = (rc_mode ? filter(x -> r_lnk[hhs[hl[x]].province] == r, 1:nh) : filter(x -> hhs[hl[x]].province == r, 1:nh))
                end

                idx_lst = [r_idxs]
                for (n_lv, hpos, intv) in intv_dataset
                    idxs = [filter(x -> hpos[hl[x]] < intv[1], r_idxs)]
                    append!(idxs, [filter(x -> intv[i] <= hpos[hl[x]] < intv[i+1], r_idxs) for i = 1:n_lv-2])
                    push!(idxs, filter(x -> intv[end-1] <= hpos[hl[x]], r_idxs))
                    append!(idx_lst, idxs)
                end
                for (n_bd, hh_val, bndr) in bndr_dataset
                    if bndr_mode == "percap"
                        idxs = [filter(x -> hh_val[hl[x]] / hhs[hl[x]].size < bndr[2], r_idxs)]
                        append!(idxs, [filter(x -> bndr[i] <= hh_val[hl[x]] / hhs[hl[x]].size < bndr[i+1], r_idxs) for i = 2:n_bd-1])
                        push!(idxs, filter(x -> hh_val[hl[x]] / hhs[hl[x]].size >= bndr[end], r_idxs))
                    elseif bndr_mode == "hhs"
                        idxs = [filter(x -> hh_val[hl[x]] < bndr[2], r_idxs)]
                        append!(idxs, [filter(x -> bndr[i] <= hh_val[hl[x]] < bndr[i+1], r_idxs) for i = 2:n_bd-1])
                        push!(idxs, filter(x -> hh_val[hl[x]] >= bndr[end], r_idxs))
                    end
                    append!(idx_lst, idxs)
                end

                for gri = 1:n_gr
                    rgi = (gri == 1 ? ri : nr + (n_gr - 1) * (ri - 1) + gri-1)
                    idxs = idx_lst[gri]
                    if length(idxs) > 0
                        if !group
                            wg_reg = hhs_wgs[idxs]
                            wg_sum = sum(hhs_wsz[idxs])
                            p_reg = (gri == 1 ? pops[y][n][r] : wg_sum)
                            etb_wg = wg_reg .* etab[idxs, :]
                            ft_de[rgi] = (de[idxs] * wg_reg)[1] / wg_sum * p_reg
                        elseif group
                            idxs_gr = [filter(x -> hhs[hl[x]].group == g, idxs) for g in gl]
                            wg_reg_gr = [hhs_wgs[idx_ls] for idx_ls in idxs_gr]
                            etb_wg_gr = [wg_reg_gr[gi] .* etab[idxs_gr[gi], :] for gi = 1:ng]
                            wg_sum_gr = [sum(hhs_wsz[idx_ls]) for idx_ls in idxs_gr]
                            p_reg = (gri == 1 ? pops[y][n][r] : (all(wg_sum_gr .> 0) ? mean(wg_sum_gr) : 0.0))
                            ft_de[rgi] = sum([de[idxs_gr[gi]] * wg_reg_gr[gi])[1] / wg_sum_gr[gi] for gi = 1:ng]) * p_reg
                        end
                        ft_p[rgi] = p_reg

                        if mode == pt_sda
                            if !group
                                et_sum = sum(etb_wg, dims=1)
                                ce_tot = sum(et_sum)
                                ce_pf = et_sum ./ ce_tot
                                ft_cepc[rgi] = ce_tot / wg_sum
                            elseif group
                                et_sum_gr = [sum(etb_wg_gr[gi], dims=1) for gi = 1:ng]
                                ce_tot_gr = [sum(et_sum_gr[gi]) for gi = 1:ng]
                                ce_pf = sum(et_sum_gr ./ wg_sum_gr)
                                ce_tot = sum(ce_pf)
                                ce_pf ./= ce_tot
                                ft_cepc[rgi] = sum(ce_tot_gr ./ wg_sum_gr)
                            end
                            ft_cspf[:,rgi] = cmat * ce_pf'
                        elseif mode in [hx_sda, cat_sda]
                            if !group
                                et_sum = vec(sum(etb_wg, dims=1))
                                ct_pf = [sum(et_sum[ci]) for ci in ct_idx]
                                ce_tot = sum(ct_pf)
                                ce_pf = zeros(Float64, ns, nc)
                                for i = 1:nc; if ct_pf[i] > 0; ce_pf[ct_idx[i], i] = et_sum[ct_idx[i]] ./ ct_pf[i] end end
                                ft_cepc[rgi] = ce_tot / wg_sum
                                if mode == hx_sda; ft_cpbc[:,rgi] = ct_pf ./ ce_tot
                                elseif mode == cat_sda; for i = 1:nc, j = 1:nc; ft_cepcbc[rgi][i][j,j] = (i == j ? ct_pf[i] / wg_sum : 1.0) end
                                end
                            elseif group
                                et_sum_gr = [vec(sum(etb_wg_gr[gi], dims=1)) for gi = 1:ng]
                                ct_pf_gr = [[sum(et_sum_gr[gi][ci]) for ci in ct_idx] for gi = 1:ng]
                                ce_tot_gr = [sum(ct_pf_gr[gi]) for gi = 1:ng]
                                ce_pf = zeros(Float64, ns, nc)
                                for i = 1:nc
                                    for gi = 1:ng
                                        if ct_pf_gr[gi][i] > 0; ce_pf[ct_idx[i], i] += et_sum_gr[gi][ct_idx[i]] ./ wg_sum_gr[gi] end
                                    end
                                    ce_pf_tot = sum(ce_pf[ct_idx[i], i])
                                    ce_pf[ct_idx[i], i] ./= ce_pf_tot
                                end
                                ft_cepc[rgi] = sum(ce_tot_gr ./ wg_sum_gr)

                                if mode == hx_sda
                                    ft_cpbc[:,rgi] = sum(ct_pf_gr ./ wg_sum_gr)
                                    ct_pf_tot = sum(ft_cpbc[:,rgi])
                                    ft_cpbc[:,rgi] ./= ct_pf_tot
                                elseif mode == cat_sda
                                    for i = 1:nc, j = 1:nc
                                        ft_cepcbc[rgi][i][j,j] = (i == j ? sum([ct_pf_gr[gi][i] / wg_sum_gr[gi] for gi = 1:ng]) : 1.0)
                                    end
                                end
                            end
                            ft_cspfbc[rgi] = cmat * ce_pf
                        end
                    end
                end
            end
            if mode == pt_sda; ft.p, ft.cepc, ft.cspf, ft.de = ft_p, ft_cepc, ft_cspf, ft_de
            elseif mode == hx_sda; ft.p, ft.cepc, ft.cspfbc, ft.cpbc, ft.de = ft_p, ft_cepc, ft_cspfbc, ft_cpbc, ft_de
            elseif mode == cat_sda; ft.p, ft.cepc, ft.cepcbc, ft.cspfbc, ft.de = ft_p, ft_cepc, ft_cepcbc, ft_cspfbc, ft_de
            else println("SDA mode error: ", mode)
            end

            if !haskey(sda_factors, y); sda_factors[y] = Dict{String, factors}() end
            sda_factors[y][n] = ft
            ieByNat[y][n] = vec(calculateEmission(y, n, mode = mode))
            deByNat[y][n], popByNat[y][n], expPcByNat[y][n] = ft.de, ft.p, ft.cepc

            mrio, etab, cmat = Float64[;;], Float64[;;], Float64[;;]
        end
        if y != baseYear; t_bp, t_tax, t_sub, v_bp, y_bp = Float64[;;], Float64[;;], Float64[;;], Float64[;;], Float64[;;] end
    end
end

function prepareDeltaFactors(target_year, base_year; nation = "", mode = "penta", reuse = false)

    global nat_list, sda_factors, dltByNat, cat_list, reg_list
    nats = (length(nation) == 0 ? nat_list : (isa(nation, String) ? [nation] : nation))
    if mode == "categorized"; nc = length(cat_list) end

    for n in nats
        t_ft, b_ft = sda_factors[target_year][n], sda_factors[base_year][n]
        rl = reg_list[base_year][n]
        nr = length(rl)

        if !reuse || !haskey(dltByNat, n)
            dltByNat[n] = Dict{Int, Any}()
            dltByNat[n][1] = t_ft.f - b_ft.f
            dltByNat[n][2] = t_ft.l - b_ft.l
            dltByNat[n][3] = t_ft.p - b_ft.p
        end

        if mode == "penta"
            dltByNat[n][4] = t_ft.cepc - b_ft.cepc
            dltByNat[n][5] = t_ft.cspf - b_ft.cspf
        elseif mode == "hexa"
            dltByNat[n][4] = t_ft.cepc - b_ft.cepc
            dltByNat[n][5] = [t_ft.cspfbc[ri] - b_ft.cspfbc[ri] for ri = 1:nr]
            dltByNat[n][6] = t_ft.cpbc - b_ft.cpbc
        elseif mode == "categorized"
            for i = 1:nc; dltByNat[n][i+3] = [t_ft.cepcbc[ri][i] - b_ft.cepcbc[ri][i] for ri = 1:nr] end
            dltByNat[n][nc+4] = [t_ft.cspfbc[ri] - b_ft.cspfbc[ri] for ri = 1:nr]
        else println("SDA mode error: ", mode)
        end
    end
end

function calculateDeltaFactors(target_year, base_year, nation, delta_factor, sub_list; mode = "penta", fl_mat = [])

    global sda_factors, dltByNat, cat_list

    yrs = [target_year, base_year]
    fts = sda_factors[target_year][nation], sda_factors[base_year][nation]
    subs = [sub_list[1:delta_factor-1]; 1; sub_list[delta_factor:end]]
    fl = []

    if mode == "penta"
        var = [fts[subs[1]].f, fts[subs[2]].l, fts[subs[3]].p, fts[subs[4]].cepc, fts[subs[5]].cspf]
        var[delta_factor] = dltByNat[nation][delta_factor]
        if length(fl_mat) > 0; fl = fl_mat else fl = var[1] .* var[2] end
        ie = vec(sum(fl * ((var[3] .* var[4])' .* var[5]), dims=1))
    elseif mode == "hexa"
        var = [fts[subs[1]].f, fts[subs[2]].l, fts[subs[3]].p, fts[subs[4]].cepc, fts[subs[5]].cspfbc, fts[subs[6]].cpbc]
        var[delta_factor] = dltByNat[nation][delta_factor]
        nt, nr = size(var[2], 1), size(var[3], 1)
        cspf = zeros(Float64, nt, nr)
        for i = 1:nr; cspf[:,i] = var[5][i] * var[6][:,i] end
        if length(fl_mat) > 0; fl = fl_mat else fl = var[1] .* var[2] end
        ie = vec(sum(fl * ((var[3] .* var[4])' .* cspf), dims=1))
    elseif mode == "categorized"
        nc = length(cat_list)
        var = [[fts[subs[1]].f, fts[subs[2]].l, fts[subs[3]].p]; [Array{Array{Float64, 2}, 1}() for i = 1:nc]]
        nt, nr = size(var[2], 1), size(var[3], 1)
        for i = 1:nc; var[i+3] = [fts[subs[i+3]].cepcbc[j][i] for j = 1:nr] end
        push!(var, fts[subs[nc+4]].cspfbc)
        var[delta_factor] = dltByNat[nation][delta_factor]
        if length(fl_mat) > 0; fl = fl_mat else fl = var[1] .* var[2] end

        ft_ce = zeros(Float64, nt, nr)
        for i = 1:nr
            ft_cepcbc = [var[j+3][i][j,j] for j = 1:nc]
            ft_ce[:, i] = sum(ft_cepcbc' .* var[nc+4][i], dims=2)
        end
        ie = vec(sum(fl * (var[3]' .* ft_ce)), dims=1)

        # iebc = zeros(Float64, nr, nc)
        # for i = 1:nr
        #     ft_cepcbc = Matrix(1.0I, nc, nc)
        #     for j = 1:nc; ft_cepcbc *= var[j+3][i] end
        #     ft_cepfbc = var[nc+4][i] * ft_cepcbc
        #     iebc[i,:] = vec(sum(fl * (var[3][i] .* ft_cepfbc), dims=1))
        # end
        # ie = vec(sum(iebc, dims=2))
    else println("SDA mode error: ", mode)
    end

    return ie, fl
end

function calculateEmission(year, nation; mode = "penta", fl_mat = [])
    global sda_factors, reg_list, cat_list
    ft = sda_factors[year][nation]
    if mode in ["hexa", "categorized"]; nt, nr, nc = size(ft.cspfbc[1], 1), length(reg_list[year][nation]), length(cat_list) end

    if length(fl_mat) > 0; fl = fl_mat else fl = ft.f .* ft.l end
    if mode == "penta"
        ie = vec(sum(fl * ((ft.p .* ft.cepc)' .* ft.cspf), dims=1))
    elseif mode == "hexa"
        ft_cspf = zeros(Float64, nt, nr)
        for i = 1:nr; ft_cspf[:,i] = ft.cspfbc[i] * ft.cpbc[:,i] end
        ie = vec(sum(fl * ((ft.p .* ft.cepc)' .* ft_cspf), dims=1))
    elseif mode == "categorized"
        ft_ce = zeros(Float64, nt, nr)
        for i = 1:nr
            ft_cepcbc = [ft.cepcbc[i][j][j,j] for j = 1:nc]
            ft_ce[:, i] = sum(ft_cepcbc' .* ft.cspfbc[i], dims=2)
        end
        ie = vec(sum(fl * (ft.p' .* ft_ce)), dims=1)

        # ie = zeros(Float64, nr, nc)
        # for i = 1:nr
        #     ft_cepcbc = Matrix(1.0I, nc, nc)
        #     for j = 1:nc; ft_cepcbc *= ft.cepcbc[i][j] end
        #     ft_cepfbc = ft.cspfbc[i] * ft_cepcbc
        #     ie[i,:] = vec(sum(fl * (ft.p[i] .* ft_cepfbc), dims=1))
        # end
        # ie = vec(sum(ie, dims=2))
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

function structuralAnalysis(target_year, base_year, nation; mode = "penta", fl_mats = [], reuse = false)

    global deltas, reg_list, sda_factors, cat_list
    if !haskey(deltas, (target_year, base_year)); deltas[(target_year, base_year)] = Dict{String, Dict{String, Array{Float64, 1}}}() end

    rl = reg_list[target_year][nation]

    deltas[(target_year, base_year)][nation] = Dict{String, Array{Float64, 1}}()

    n_factor = Dict("penta" => 5, "hexa" => 6, "categorized" => (4+length(cat_list)))
    nf = n_factor[mode]

    for r in rl; deltas[(target_year, base_year)][nation][r] = zeros(Float64, nf) end

    nk = nf - 1
    nr = length(rl)
    dlt_repo = Array{Array{delta, 1}, 1}()

    wghs = Dict(0:nk .=> [factorial(nk - k) * factorial(k) for k = 0:nk])
    subs_list = generateAllCombination(Array{Int, 1}(), nf, elements = [1,2])
    wgh_subs = Array{Tuple{Float64, Array{Int, 1}}, 1}()

    for sl in subs_list; push!(wgh_subs , (wghs[count(x->x==2, sl)], sl)) end
    if reuse && length(fl_mats) == 0; fl_mats = Dict((i,j) => zeros(0,0) for i=1:3, j=1:3) end

    for i = 1:nf
        tot_wgh, dlts = 0.0, zeros(Float64, nr)
        # dlt_list = Array{delta, 1}()
        for (wgh, sl) in wgh_subs
            if i == 1; fsl, lsl = 3, sl[1] elseif i == 2; fsl, lsl = sl[1], 3 else fsl, lsl = sl[1], sl[2] end

            if !reuse; fl = [] else fl = fl_mats[(fsl, lsl)] end
            ie, fl = calculateDeltaFactors(target_year, base_year, nation, i, sl, mode = mode, fl_mat = fl)
            dlt_vec = vec(ie)
            dlts .+= wgh .* dlt_vec
            # push!(dlt_list, delta(i, nf, sub_list = sl, weight = wgh, delta_value = dlt_vec))
            tot_wgh += wgh
            if reuse; fl_mats[(fsl, lsl)] = fl end
        end
        dlts ./= tot_wgh
        for j = 1:nr; deltas[(target_year, base_year)][nation][rl[j]][i] = dlts[j]  end
        # push!(dlt_repo, dlt_list)
    end

    return dlt_repo, fl_mats
end

function estimateConfidenceIntervals(year, nation = []; iter = 10000, ci_rate = 0.95, resample_size = 0, replacement = true, region = "district", exp_mode = "gis", group = false)
    # bootstrap method
    # ci_per: confidence interval percentage
    # replacement: [true] sampling with replacement
    # if resample_size: [0] resample_size = sample_size
    # group: [true] divide households by (survey) groups

    global nat_list, reg_list, hh_list, households, pops, cat_list, hh_cf
    global ci_ie, ci_de, ci_cf, ci_cfpc, hh_ie, hh_de, ieByNat, deByNat, cfByNat

    if resample_size == 0; replacement = true end
    if isa(year, Number); year = [year] end
    nats = (length(nation) == 0 ? nat_list : (isa(nation, String) ? [nation] : nation))
    nc = length(cat_list)
    rc_mode = (lowercase(exp_mode) == "gis")

    for y in year
        if !haskey(ci_ie, y); ci_ie[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
        if !haskey(ci_de, y); ci_de[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
        if !haskey(ci_cf, y); ci_cf[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
        if !haskey(ci_cfpc, y); ci_cfpc[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
        if !haskey(ieByNat, y); ieByNat[y] = Dict{String, Array{Float64, 1}}() end
        if !haskey(deByNat, y); deByNat[y] = Dict{String, Array{Float64, 1}}() end
        if !haskey(cfByNat, y); cfByNat[y] = Dict{String, Array{Float64, 1}}() end
        for n in nats
            if !haskey(ci_ie[y], n); ci_ie[y][n] = Dict{String, Tuple{Float64, Float64}}() end
            if !haskey(ci_de[y], n); ci_de[y][n] = Dict{String, Tuple{Float64, Float64}}() end
            if !haskey(ci_cf[y], n); ci_cf[y][n] = Dict{String, Tuple{Float64, Float64}}() end
            if !haskey(ci_cfpc[y], n); ci_cfpc[y][n] = Dict{String, Tuple{Float64, Float64}}() end

            rl, hl, hhs = reg_list[y][n], hh_list[y][n], households[y][n]
            if rc_mode; r_lnk = reg_linked[y][n] end
            ie, de, hcf = hh_ie[y][n], hh_de[y][n], hh_cf[y][n]
            nh, nr = length(hl), length(rl)
            if group; gl = gr_list[y][n]; ng = length(gl) end
            ieByNat[y][n], deByNat[y][n], cfByNat[y][n] = zeros(Float64, nr), zeros(Float64, nr), zeros(Float64, nr)

            hhs_wgs = [hhs[h].popwgh for h in hl]
            hhs_siz = [hhs[h].size for h in hl]
            hhs_wsz = hhs_wgs .* hhs_siz
            iew = ie .* hhs_wgs
            dew = de .* hhs_wgs
            hcfw = hcf .* hhs_wgs

            for ri = 1:nr
                r = rl[ri]
                p_reg = pops[y][n][r]

                if region == "district"
                    idxs = (rc_mode ? filter(x -> r_lnk[hhs[hl[x]].district] == r, 1:nh) : filter(x -> hhs[hl[x]].district == r, 1:nh))
                elseif region == "province"
                    idxs = (rc_mode ? filter(x -> r_lnk[hhs[hl[x]].province] == r, 1:nh) : filter(x -> hhs[hl[x]].province == r, 1:nh))
                end
                ie_vals, de_vals, cf_vals, cfpc_vals = zeros(Float64, iter), zeros(Float64, iter), zeros(Float64, iter), zeros(Float64, iter)

                if !group
                    hwg_reg = hhs_wsz[idxs]
                    iew_reg = iew[idxs]
                    dew_reg = dew[idxs]
                    hcfw_reg = hcfw[idxs]

                    nsam = (resample_size == 0 ? length(idxs) : resample_size)

                    for i = 1:iter
                        if replacement; re_idx = [trunc(Int, nsam * rand())+1 for x = 1:nsam]
                        else re_idx = sortperm([rand() for x = 1:length(idxs)])[1:nsam]
                        end
                        wg_sum = sum(hwg_reg[re_idx])
                        ie_vals[i] = sum(iew_reg[re_idx]) / wg_sum * p_reg
                        de_vals[i] = sum(dew_reg[re_idx]) / wg_sum * p_reg
                        cfpc_vals[i] = sum(hcfw_reg[re_idx]) / wg_sum
                        cf_vals[i] = cfpc_vals[i] * p_reg
                    end
                elseif group
                    idxs_gr = [filter(x -> hhs[hl[x]].group == g, idxs) for g in gl]
                    hwg_reg_gr = [hhs_wsz[il] for il in idxs_gr]
                    iew_reg_gr = [iew[il] for il in idxs_gr]
                    dew_reg_gr = [dew[il] for il in idxs_gr]
                    hcfw_reg_gr = [hcfw[il] for il in idxs_gr]

                    nsam_gr = [length(il) for il in idxs_gr]
                    if resample_size > 0; nsam_gr = trunc.(Int, nsam_gr .* (resample_size / length(idxs))) end

                    for i = 1:iter
                        if replacement; re_idx_gr = [[trunc(Int, ns_gr * rand())+1 for x = 1:ns_gr] for ns_gr in nsam_gr]
                        else re_idx_gr = [sortperm([rand() for x = 1:length(idxs_gr[gi])])[1:nsam_gr[gi]] for gi = 1:ng]
                        end
                        wg_sum_gr = [sum(hwg_reg_gr[gi][re_idx_gr[gi]]) for gi = 1:ng]
                        ie_vals[i] = sum([sum(iew_reg_gr[gi][re_idx_gr[gi]]) / wg_sum_gr[gi] for gi = 1:ng]) * p_reg
                        de_vals[i] = sum([sum(dew_reg_gr[gi][re_idx_gr[gi]]) / wg_sum_gr[gi] for gi = 1:ng]) * p_reg
                        cfpc_vals[i] = sum([sum(hcfw_reg_gr[gi][re_idx_gr[gi]]) / wg_sum_gr[gi] for gi = 1:ng])
                        cf_vals[i] = cfpc_vals[i] * p_reg
                    end
                end

                sort!(ie_vals)
                sort!(de_vals)
                sort!(cf_vals)
                sort!(cfpc_vals)

                l_idx, u_idx = trunc(Int, (1 - ci_rate) / 2 * iter) + 1, trunc(Int, ((1 - ci_rate) / 2 + ci_rate) * iter) + 1
                ci_ie[y][n][r] = (ie_vals[l_idx], ie_vals[u_idx])
                ci_de[y][n][r] = (de_vals[l_idx], de_vals[u_idx])
                ci_cf[y][n][r] = (cf_vals[l_idx], cf_vals[u_idx])
                ci_cfpc[y][n][r] = (cfpc_vals[l_idx], cfpc_vals[u_idx])

                if !group
                    wg_sum = sum(hwg_reg)
                    ieByNat[y][n][ri] = sum(iew_reg) / wg_sum * p_reg
                    deByNat[y][n][ri] = sum(dew_reg) / wg_sum * p_reg
                    cfByNat[y][n][ri] = sum(hcfw_reg) / wg_sum * p_reg
                elseif group
                    wg_sum_gr = [sum(hwg_reg_gr[gi]) for gi = 1:ng]
                    ieByNat[y][n][ri] = sum(sum.(iew_reg_gr) ./ wg_sum_gr) * p_reg
                    deByNat[y][n][ri] = sum(sum.(dew_reg_gr) ./ wg_sum_gr) * p_reg
                    cfByNat[y][n][ri] = sum(sum.(hcfw_reg_gr) ./ wg_sum_gr) * p_reg
                end
            end
        end
    end
end

function estimateSdaCi(target_year, base_year, nation = [], mrioPath = ""; iter = 10000, ci_rate = 0.95, mode="penta", region = "district", exp_mode = "gis",
                        resample_size = 0, replacement = true, visible = false, reuse = true, group = false,
                        min_itr = 1000, chk_itr = 10, err_crt = 0.0001, visible_iter = 0)
    # bootstrap method
    # ci_per: confidence interval percentage
    # replacement: [true] sampling with replacement
    # if resample_size: [0] resample_size = sample_size
    # group: [true] divide households by (survey) groups

    er_limit = err_crt          # maximum acceptable error
    iter_min = min_itr          # minimum iterations (maximum = 'iter')
    er_chk_iter = chk_itr       # assess variable change rates every 'chk_itr' interation(s)

    pt_mode, hx_mode, cat_mode = "penta", "hexa", "categorized"

    global nat_list, reg_list, hh_list, pops, sda_factors, gr_list
    global ci_ie, ci_de, ci_sda, hh_ie, hh_de, ieByNat, deByNat, exp_table, conc_mat_wgh
    global mrio_tabs, l_factor, f_factor, mrio_tabs_conv, deltas, samples_gr

    ty, by = target_year, base_year
    if resample_size == 0; replacement = true end
    if length(nation) == 0; nats = nat_list
    elseif isa(nation, String); nats = [nation]
    elseif isa(nation, Array{String, 1}); nats = nation
    end
    rc_mode = (lowercase(exp_mode) == "gis")

    for y in [ty, by]
        ci_ie[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}()
        ci_de[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}()
        ieByNat[y] = Dict{String, Array{Float64, 1}}()
        deByNat[y] = Dict{String, Array{Float64, 1}}()
    end
    ci_sda[(ty,by)] = Dict{String, Dict{String, Array{Tuple{Float64, Float64}, 1}}}()

    t_bp, t_tax, t_sub, v_bp, y_bp = setMrioTables(ty, mrioPath)
    ft_by = factors()
    if !haskey(l_factor, by)
        l_factor[by], lx_by = calculateLeontief(mrio_tabs[by])
        f_factor[by] = calculateIntensity(mrio_tabs[by], lx_by)
        lx_by = []
    end
    mrio_by, ft_by.l, ft_by.f  = mrio_tabs[by], l_factor[by], f_factor[by]
    l_idx, u_idx = trunc(Int, (1 - ci_rate) / 2 * iter) + 1, trunc(Int, ((1 - ci_rate) / 2 + ci_rate) * iter) + 1

    st = time()
    for n in nats
        if visible; print(" ", by,"_",ty ,":") end

        sda_factors[ty], sda_factors[by] = Dict{String, factors}(), Dict{String, factors}()

        ci_ie[ty][n], ci_ie[by][n] = Dict{String, Tuple{Float64, Float64}}(), Dict{String, Tuple{Float64, Float64}}()
        ci_de[ty][n], ci_de[by][n] = Dict{String, Tuple{Float64, Float64}}(), Dict{String, Tuple{Float64, Float64}}()
        ci_sda[(ty,by)][n] = Dict{String, Array{Tuple{Float64, Float64}, 1}}()

        rl_ty, hl_ty, hhs_ty = reg_list[ty][n], hh_list[ty][n], households[ty][n]
        rl_by, hl_by, hhs_by = reg_list[by][n], hh_list[by][n], households[by][n]
        ie_ty, de_ty = hh_ie[ty][n], hh_de[ty][n]
        ie_by, de_by = hh_ie[by][n], hh_de[by][n]
        nh_ty, nr_ty, nh_by, nr_by = length(hl_ty), length(rl_ty), length(hl_by), length(rl_by)
        ieByNat[ty][n], deByNat[ty][n] = zeros(Float64, nr_ty), zeros(Float64, nr_ty)
        ieByNat[by][n], deByNat[by][n] = zeros(Float64, nr_by), zeros(Float64, nr_by)
        nsam = Dict{Int, Array{Int, 1}}()

        if rc_mode; r_lnk_ty, r_lnk_by = reg_linked[ty][n], reg_linked[by][n] end
        if group; gl = gr_list[by][n]; ng = length(gl) end

        etab_ty, cmat_ty = exp_table[ty][n], conc_mat_wgh[ty][n]
        etab_by, cmat_by = exp_table[by][n], conc_mat_wgh[by][n]

        ft_ty = factors()
        convertTable(ty, n, by, mrioPath, t_bp = t_bp, t_tax = t_tax, t_sub = t_sub, v_bp = v_bp, y_bp = y_bp)
        mrio_ty = mrio_tabs_conv[ty][n]
        ft_ty.l, lx_ty = calculateLeontief(mrio_ty)
        ft_ty.f = calculateIntensity(mrio_ty, lx_ty)
        lx_ty = []

        nt_ty = size(mrio_ty.t, 1)
        nt_by = size(mrio_by.t, 1)

        ft_ty_p, ft_by_p = zeros(nr_ty), zeros(nr_by)
        idxs_ty, idxs_by = [Array{Int, 1}() for ri = 1:nr_ty], [Array{Int, 1}() for ri = 1:nr_by]
        nsam[ty], nsam[by] = zeros(Int, nr_ty), zeros(nr_by)
        if group
            nsam_gr = Dict{Int, Array{Array{Int, 1}, 1}}()
            idxs_ty_gr, idxs_by_gr = [Array{Array{Int, 1}, 1}() for ri = 1:nr_ty], [Array{Array{Int, 1}, 1}() for ri = 1:nr_by]
            nsam_gr[ty], nsam_gr[by] = [zeros(Int, ng) for ri = 1:nr_ty], [zeros(Int, ng) for ri = 1:nr_by]
        end

        ie_vals_ty, de_vals_ty = [zeros(Float64, 0) for i=1:nr_ty], [zeros(Float64, 0) for i=1:nr_ty]
        ie_vals_by, de_vals_by = [zeros(Float64, 0) for i=1:nr_by], [zeros(Float64, 0) for i=1:nr_by]
        cepc_vals, cspf_vals = [zeros(Float64, 0) for i=1:nr_by], [zeros(Float64, 0) for i=1:nr_by]
        fls = []

        hhs_wgs_ty = [hhs_ty[h].popwgh for h in hl_ty]
        hhs_wgs_by = [hhs_by[h].popwgh for h in hl_by]
        hhs_wsz_ty = hhs_wgs_ty .* [hhs_ty[h].size for h in hl_ty]
        hhs_wsz_by = hhs_wgs_by .* [hhs_by[h].size for h in hl_by]
        iew_ty = ie_ty .* hhs_wgs_ty
        iew_by = ie_by .* hhs_wgs_by
        dew_ty = de_ty .* hhs_wgs_ty
        dew_by = de_by .* hhs_wgs_by
        etw_ty = hhs_wgs_ty .* etab_ty
        etw_by = hhs_wgs_by .* etab_by

        for ri = 1:nr_by
            r = rl_by[ri]
            ft_ty_p[ri] = pops[ty][n][r]
            ft_by_p[ri] = pops[by][n][r]

            if region == "district"
                if rc_mode
                    idxs_ty[ri] = filter(x -> r_lnk_ty[hhs_ty[hl_ty[x]].district] == r, 1:nh_ty)
                    idxs_by[ri] = filter(x -> r_lnk_by[hhs_by[hl_by[x]].district] == r, 1:nh_by)
                else
                    idxs_ty[ri] = filter(x -> hhs_ty[hl_ty[x]].district == r, 1:nh_ty)
                    idxs_by[ri] = filter(x -> hhs_by[hl_by[x]].district == r, 1:nh_by)
                end
            elseif region == "province"
                if rc_mode
                    idxs_ty[ri] = filter(x -> r_lnk_ty[hhs_ty[hl_ty[x]].province] == r, 1:nh_ty)
                    idxs_by[ri] = filter(x -> r_lnk_by[hhs_by[hl_by[x]].province] == r, 1:nh_by)
                else
                    idxs_ty[ri] = filter(x -> hhs_ty[hl_ty[x]].province == r, 1:nh_ty)
                    idxs_by[ri] = filter(x -> hhs_by[hl_by[x]].province == r, 1:nh_by)
                end
            end
            nsam[ty][ri], nsam[by][ri] = (resample_size == 0 ? [length(idxs_ty[ri]), length(idxs_by[ri])] : [resample_size, resample_size])

            if group
                idxs_ty_gr[ri] = [filter(x -> hhs_ty[hl_ty[x]].group == g, idxs_ty[ri]) for g in gl]
                idxs_by_gr[ri] = [filter(x -> hhs_by[hl_by[x]].group == g, idxs_by[ri]) for g in gl]
                nsam_gr[ty][ri] = [length(il) for il in idxs_ty_gr[ri]]
                nsam_gr[by][ri] = [length(il) for il in idxs_by_gr[ri]]
                if resample_size > 0
                    nsam_gr[ty][ri]  = trunc.(Int, nsam_gr[ty][ri] .* (resample_size / length(idxs_ty_gr[ri])))
                    nsam_gr[by][ri]  = trunc.(Int, nsam_gr[by][ri] .* (resample_size / length(idxs_by_gr[ri])))
                end
            end
        end
        samples_gr[n] = (!group ? nsam : nsam_gr)

        ie_prv_l_ty, ie_prv_u_ty, ie_prv_l_by, ie_prv_u_by = zeros(Float64, nr_by), zeros(Float64, nr_by), zeros(Float64, nr_by), zeros(Float64, nr_by)
        cepc_prv_l, cepc_prv_u, cspf_prv_l, cspf_prv_u  = zeros(Float64, nr_by), zeros(Float64, nr_by), zeros(Float64, nr_by), zeros(Float64, nr_by)

        er = 1.0
        er_c = er_limit

        i = 0
        while (er > er_c && i < iter)
            i += 1

            ft_ty_de, ft_by_de = zeros(nr_ty), zeros(nr_by)
            if mode == pt_mode
                ft_ty_cepc, ft_ty_cspf, ft_ty_de = zeros(nr_ty), zeros(nt_ty, nr_ty), zeros(nr_ty)
                ft_by_cepc, ft_by_cspf, ft_by_de = zeros(nr_by), zeros(nt_by, nr_by), zeros(nr_by)
            # elseif mode == hx_sda; ft_cepc, ft_cspf, ft_cpbc = zeros(nr), zeros(nt, nc, nr), zeros(nc, nr)
            # elseif mode == cat_sda; ft_cepc, ft_cspf = [zeros(nc, nc, nr) for i=1:nc], zeros(nt, nc, nr)
            end

            for ri = 1:nr_by
                r = rl_by[ri]

                if !group
                    if replacement; re_idx_ty, re_idx_by = [[trunc(Int, ns * rand())+1 for x = 1:ns] for ns in [nsam[ty][ri], nsam[by][ri]]]
                    else re_idx_ty, re_idx_by = [sortperm([rand() for x = 1:ns]) for ns in [nsam[ty][ri], nsam[by][ri]]]
                    end

                    idt, idb = idxs_ty[ri][re_idx_ty], idxs_by[ri][re_idx_by]
                    wrt, wrb = hhs_wgs_ty[idt], hhs_wgs_by[idb]
                    wst, wsb = sum(hhs_wsz_ty[idt]), sum(hhs_wsz_by[idb])

                    push!(ie_vals_ty[ri], sum(iew_ty[idt]) / wst * ft_ty_p[ri])
                    push!(ie_vals_by[ri], sum(iew_by[idb]) / wsb * ft_by_p[ri])
                    push!(de_vals_ty[ri], sum(dew_ty[idt]) / wst * ft_ty_p[ri])
                    push!(de_vals_by[ri], sum(dew_by[idb]) / wsb * ft_by_p[ri])

                    etb_wg_ty = etw_ty[idt, :]
                    etb_wg_by = etw_by[idb, :]

                    if mode == pt_mode
                        et_sum_ty = sum(etb_wg_ty, dims=1)
                        ce_tot_ty = sum(et_sum_ty)
                        ce_pf_ty = et_sum_ty ./ ce_tot_ty
                        ft_ty_cepc[ri] = ce_tot_ty / wst
                        ft_ty_cspf[:,ri] = cmat_ty * ce_pf_ty'

                        et_sum_by = sum(etb_wg_by, dims=1)
                        ce_tot_by = sum(et_sum_by)
                        ce_pf_by = et_sum_by ./ ce_tot_by
                        ft_by_cepc[ri] = ce_tot_by / wsb
                        ft_by_cspf[:,ri] = cmat_by * ce_pf_by'
                    # elseif mode in [hx_mode, cat_mode]
                    #     et_sum = vec(sum(etb_wg, dims=1))
                    #     ct_pf = [sum(et_sum[ci]) for ci in ct_idx]
                    #     ce_tot = sum(ct_pf)
                    #     ce_pf = zeros(Float64, ns, nc)
                    #     for i = 1:nc; ce_pf[ct_idx[i], i] = (ct_pf[i] > 0 ? (et_sum[ct_idx[i]] ./ ct_pf[i]) : [0 for x = 1:length(ct_idx[i])]) end
                    #     ft_cspf[:,:,ri] = cmat * ce_pf
                    #     ft_p[ri] = p_reg
                    #
                    #     if mode == hx_mode
                    #         ft_cepc[ri] = ce_tot / wg_sum
                    #         ft_cpbc[:,ri] = ct_pf ./ ce_tot
                    #     elseif mode == cat_mode
                    #         for i = 1:nc, j = 1:nc; ft_cepc[i][j,j,ri] = (i == j ? ct_pf[i] / wg_sum : 1.0) end
                    #     end
                    end
                elseif group
                    if replacement
                        re_idx_ty_gr, re_idx_by_gr = [[[trunc(Int, ns * rand())+1 for x = 1:ns] for ns in nsam_gr[y][ri]] for y in [ty, by]]
                    else
                        re_idx_ty_gr = [sortperm([rand() for x = 1:length(idxs_ty_gr[ri][gi])])[1:nsam_gr[ty][gi]] for gi = 1:ng]
                        re_idx_by_gr = [sortperm([rand() for x = 1:length(idxs_by_gr[ri][gi])])[1:nsam_gr[by][gi]] for gi = 1:ng]
                    end

                    idt_gr, idb_gr = [idxs_ty_gr[ri][ridx] for ridx in re_idx_ty_gr], [idxs_by_gr[ri][ridx] for ridx in re_idx_by_gr]
                    wrt_gr, wrb_gr = [hhs_wgs_ty[il] for il in idt_gr], [hhs_wgs_by[il] for il in idb_gr]
                    wst_gr, wsb_gr = [sum(hhs_wsz_ty[il]) for il in idt_gr], [sum(hhs_wsz_by[il]) for il in idb_gr]

                    push!(ie_vals_ty[ri], sum([sum(iew_ty[idt_gr[gi]]) / wst_gr[gi] for gi = 1:ng]) * ft_ty_p[ri])
                    push!(ie_vals_by[ri], sum([sum(iew_by[idb_gr[gi]]) / wsb_gr[gi] for gi = 1:ng]) * ft_by_p[ri])
                    push!(de_vals_ty[ri], sum([sum(dew_ty[idt_gr[gi]]) / wst_gr[gi] for gi = 1:ng]) * ft_ty_p[ri])
                    push!(de_vals_by[ri], sum([sum(dew_by[idb_gr[gi]]) / wsb_gr[gi] for gi = 1:ng]) * ft_by_p[ri])

                    etb_wg_ty_gr = [etw_ty[il, :] for il in idt_gr]
                    etb_wg_by_gr = [etw_by[il, :] for il in idb_gr]

                    if mode == pt_mode
                        et_sum_ty_gr = [sum(etb_wg_ty_gr[gi], dims=1) for gi = 1:ng]
                        ce_tot_ty_gr = [sum(et_sum_ty_gr[gi]) for gi = 1:ng]
                        ce_pf_ty = sum(et_sum_ty_gr ./ wst_gr)
                        ce_tot_ty = sum(ce_pf_ty)
                        ce_pf_ty ./= ce_tot_ty
                        ft_ty_cepc[ri] = sum(ce_tot_ty_gr ./ wst_gr)
                        ft_ty_cspf[:,ri] = cmat_ty * ce_pf_ty'

                        et_sum_by_gr = [sum(etb_wg_by_gr[gi], dims=1) for gi = 1:ng]
                        ce_tot_by_gr = [sum(et_sum_by_gr[gi]) for gi = 1:ng]
                        ce_pf_by = sum(et_sum_by_gr ./ wsb_gr)
                        ce_tot_by = sum(ce_pf_by)
                        ce_pf_by ./= ce_tot_by
                        ft_by_cepc[ri] = sum(ce_tot_by_gr ./ wsb_gr)
                        ft_by_cspf[:,ri] = cmat_by * ce_pf_by'
                    end
                end
                ft_ty_de[ri] = de_vals_ty[ri][i]
                ft_by_de[ri] = de_vals_by[ri][i]
            end

            if mode == pt_mode
                ft_ty.p, ft_ty.cepc, ft_ty.cspf, ft_ty.de = ft_ty_p, ft_ty_cepc, ft_ty_cspf, ft_ty_de
                ft_by.p, ft_by.cepc, ft_by.cspf, ft_by.de = ft_by_p, ft_by_cepc, ft_by_cspf, ft_by_de
            # elseif mode == hx_sda; ft.p, ft.cepc, ft.cspfbc, ft.cpbc, ft.de = ft_p, ft_cepc, ft_cspf, ft_cpbc, ft_de
            # elseif mode == cat_sda; ft.p, ft.cepcbc, ft.cspfbc, ft.de = ft_p, ft_cepc, ft_cspf, ft_de
            # else println("SDA mode error: ", mode)
            end

            sda_factors[ty][n], sda_factors[by][n] = ft_ty, ft_by
            prepareDeltaFactors(ty, by, nation = n, mode = mode, reuse = reuse)
            fls = structuralAnalysis(ty, by, n, mode = mode, fl_mats = fls, reuse = reuse)[2]

            for ri = 1:nr_by
                push!(cepc_vals[ri], deltas[(ty, by)][n][rl_by[ri]][4])
                push!(cspf_vals[ri], deltas[(ty, by)][n][rl_by[ri]][5])
            end

            if i >= iter_min && i % er_chk_iter == 0
                li, ui = trunc(Int, (1 - ci_rate) / 2 * i) + 1, trunc(Int, ((1 - ci_rate) / 2 + ci_rate) * i) + 1
                current_vals, previous_vals = [], []
                for ri = 1:nr_by
                    r = rl_by[ri]
                    sort!(ie_vals_ty[ri]); sort!(ie_vals_by[ri]); sort!(de_vals_ty[ri]); sort!(de_vals_by[ri])
                    sort!(cepc_vals[ri]); sort!(cspf_vals[ri])

                    ci_ie[ty][n][r] = (ie_vals_ty[ri][li], ie_vals_ty[ri][ui])
                    ci_ie[by][n][r] = (ie_vals_by[ri][li], ie_vals_by[ri][ui])
                    ci_de[ty][n][r] = (de_vals_ty[ri][li], de_vals_ty[ri][ui])
                    ci_de[by][n][r] = (de_vals_by[ri][li], de_vals_by[ri][ui])
                    ci_sda[(ty,by)][n][r] = [(cepc_vals[ri][li], cepc_vals[ri][ui]), (cspf_vals[ri][li], cspf_vals[ri][ui])]

                    append!(current_vals, [ie_vals_ty[ri][li], ie_vals_ty[ri][ui], ie_vals_by[ri][li], ie_vals_by[ri][ui], cepc_vals[ri][li], cepc_vals[ri][ui], cspf_vals[ri][li], cspf_vals[ri][ui]])
                    append!(previous_vals, [ie_prv_l_ty[ri], ie_prv_u_ty[ri], ie_prv_l_by[ri], ie_prv_u_by[ri], cepc_prv_l[ri], cepc_prv_u[ri], cspf_prv_l[ri], cspf_prv_u[ri]])
                end
                ers = abs.((current_vals - previous_vals) ./ previous_vals)
                ers[isnan.(ers)] .= 0
                ers[isinf.(ers)] .= 0
                er = maximum(ers)

                for ri = 1:nr_by
                    r = rl_by[ri]
                    ie_prv_l_ty[ri], ie_prv_u_ty[ri] = ci_ie[ty][n][r]
                    ie_prv_l_by[ri], ie_prv_u_by[ri] = ci_ie[by][n][r]
                    cepc_prv_l[ri], cepc_prv_u[ri] = ci_sda[(ty,by)][n][r][1]
                    cspf_prv_l[ri], cspf_prv_u[ri] = ci_sda[(ty,by)][n][r][2]
                end
            end
            if visible_iter > 0 && i % visible_iter == 0; print(" ", i) end
        end

        for ri = 1:nr_by
            if !group
                idt, idb = idxs_ty[ri], idxs_by[ri]
                wg_sum_ty = sum(hhs_wsz_ty[idt])
                wg_sum_by = sum(hhs_wsz_by[idb])
                ieByNat[ty][n][ri] = sum(iew_ty[idt]) / wg_sum_ty * ft_ty_p[ri]
                ieByNat[by][n][ri] = sum(iew_by[idb]) / wg_sum_by * ft_by_p[ri]
                deByNat[ty][n][ri] = sum(dew_ty[idt]) / wg_sum_ty * ft_ty_p[ri]
                deByNat[by][n][ri] = sum(dew_by[idb]) / wg_sum_by * ft_by_p[ri]
            elseif group
                idt_gr, idb_gr = idxs_ty_gr[ri], idxs_by_gr[ri]
                wg_sum_ty_gr = [sum(hhs_wsz_ty[idt]) for idt in idt_gr]
                wg_sum_by_gr = [sum(hhs_wsz_by[idb]) for idb in idb_gr]
                ieByNat[ty][n][ri] = sum([sum(iew_ty[idt_gr[gi]]) / wg_sum_ty_gr[gi] for gi = 1:ng]) * ft_ty_p[ri]
                ieByNat[by][n][ri] = sum([sum(iew_by[idb_gr[gi]]) / wg_sum_by_gr[gi] for gi = 1:ng]) * ft_by_p[ri]
                deByNat[ty][n][ri] = sum([sum(dew_ty[idt_gr[gi]]) / wg_sum_ty_gr[gi] for gi = 1:ng]) * ft_ty_p[ri]
                deByNat[by][n][ri] = sum([sum(dew_by[idb_gr[gi]]) / wg_sum_by_gr[gi] for gi = 1:ng]) * ft_by_p[ri]
            end
        end

        if visible
            elap = floor(Int, time() - st); (eMin, eSec) = fldmod(elap, 60); (eHr, eMin) = fldmod(eMin, 60)
            print(eHr,":",eMin,":",eSec," elapsed,\t", i, " iterations")
        end
    end
end

function estimateSdaCiByGroup(target_year, base_year, nation = [], mrioPath = ""; iter = 10000, ci_rate = 0.95, mode="penta",
                            resample_size = 0, replacement = false, visible = false, reuse = true, exp_mode = "gis", group = false,
                            cf_intv = [], inc_intv = [], hpos_cf = [], hpos_inc = [], cf_bndr = [], inc_bndr = [],
                            min_itr = 1000, chk_itr = 10, err_crt = 0.0001, visible_iter = 0, bndr_mode = "percap")
    # bootstrap method
    # ci_per: confidence interval percentage
    # replacement: [true] sampling with replacement
    # if resample_size: [0] resample_size = sample_size
    # group: [true] divide households by (survey) groups

    # grouping = cf_intv: CF pre capita intervals (stacked proportion)
    # grouping = inc_intv: income per capita intervals (stacked proportion)
    # grouping = cf_bndr: CF pre capita boundaries (absolute limits)
    # grouping = inc_bndr: income per capita boundaries (absolute limits)
    # bndr_mode: "percap" per capita income or CF for boundary vales, "hhs" household income or CF

    er_limit = err_crt          # maximum acceptable error
    iter_min = min_itr          # minimum iterations (maximum = 'iter')
    er_chk_iter = chk_itr       # check every 'er_chk_iter' interation

    pt_mode, hx_mode, cat_mode = "penta", "hexa", "categorized"
    pd_tag = Dict(1 => "densly", 2 => "inter", 3 => "sparsly")

    global nat_list, reg_list, hh_list, pops, sda_factors, gr_list
    global ci_ie, ci_de, ci_sda, hh_ie, hh_de, ieByNat, deByNat, exp_table, conc_mat_wgh
    global mrio_tabs, l_factor, f_factor, mrio_tabs_conv, deltas, samples_gr, hh_cf, cat_hhl

    ty, by = target_year, base_year
    if resample_size == 0; replacement = true end
    if length(nation) == 0; nats = nat_list
    elseif isa(nation, String); nats = [nation]
    elseif isa(nation, Array{String, 1}); nats = nation
    end
    rc_mode = (lowercase(exp_mode) == "gis")

    n_cf, n_inc = length(cf_intv), length(inc_intv)
    n_cfb, n_incb = length(cf_bndr), length(inc_bndr)
    n_gr = 1 + n_cf + n_inc + n_cfb + n_incb

    if n_cf > 0
        cf_intv_lb = ["_bottom_" * string(ceil(Int, cf_intv[1] * 100)) * "%"]
        append!(cf_intv_lb, ["_middle_" * string(ceil(Int, (cf_intv[i] - cf_intv[i-1]) * 100)) * "%" for i = 2:n_cf-1])
        append!(cf_intv_lb, ["_top_" * string(ceil(Int, (cf_intv[end] - cf_intv[end-1]) * 100)) * "%"])
    end
    if n_inc > 0
        inc_intv_lb = ["_bottom_" * string(ceil(Int, inc_intv[1] * 100)) * "%"]
        append!(inc_intv_lb, ["_middle_" * string(ceil(Int, (inc_intv[i] - inc_intv[i-1]) * 100)) * "%" for i = 2:n_inc-1])
        append!(inc_intv_lb, ["_top_" * string(ceil(Int, (inc_intv[end] - inc_intv[end-1]) * 100)) * "%"])
    end
    if n_cfb > 0
        cf_bndr_lb = [[""];["_" * string(round(cf_bndr[i], digits = 2)) * "≤" for i = 2:n_cfb]]
        for i = 1:n_cfb-1; cf_bndr_lb[i] *= "_<" * string(round(cf_bndr[i+1], digits = 2)) end
    end
    if n_incb > 0
        inc_bndr_lb = [[""];["_" * string(round(inc_bndr[i]/1000, digits = 2)) * "≤" for i = 2:n_incb]]
        for i = 1:n_incb-1; inc_bndr_lb[i] *= "_<" * string(round(inc_bndr[i+1]/1000, digits = 2)) end
    end

    for y in [ty, by]
        ci_ie[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}()
        ci_de[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}()
        ieByNat[y] = Dict{String, Array{Float64, 1}}()
        deByNat[y] = Dict{String, Array{Float64, 1}}()
    end
    ci_sda[(ty,by)] = Dict{String, Dict{String, Array{Tuple{Float64, Float64}, 1}}}()

    t_bp, t_tax, t_sub, v_bp, y_bp = setMrioTables(ty, mrioPath)
    l_idx, u_idx = trunc(Int, (1 - ci_rate) / 2 * iter) + 1, trunc(Int, ((1 - ci_rate) / 2 + ci_rate) * iter) + 1

    st = time()
    for n in nats
        if visible; print(" ", by,"_",ty ,":") end

        ci_sda[(ty,by)][n] = Dict{String, Array{Tuple{Float64, Float64}, 1}}()
        ft_p = Dict{Int, Array{Float64, 1}}()
        idx_ls = Dict{Int, Array{Array{Int, 1}, 1}}()
        wg_reg = Dict{Int, Array{Float64, 1}}()
        wg_hhs = Dict{Int, Array{Float64, 1}}()
        etb_wg = Dict{Int, Array{Float64, 2}}()
        nsam = Dict{Int, Array{Int, 1}}()
        ie_vals = Dict{Int, Array{Array{Float64, 1}, 1}}()
        de_vals = Dict{Int, Array{Array{Float64, 1}, 1}}()
        ie, ie_wg = Dict{Int, Array{Float64, 1}}(), Dict{Int, Array{Float64, 1}}()
        de, de_wg = Dict{Int, Array{Float64, 1}}(), Dict{Int, Array{Float64, 1}}()
        ie_prv_l = Dict{Int, Array{Float64, 1}}()
        ie_prv_u = Dict{Int, Array{Float64, 1}}()
        cmat = Dict{Int, Array{Float64, 2}}()
        nr, nt, n_nt = 0, 0, 0
        rl = Array{String, 1}()
        if rc_mode; r_lnk_ty, r_lnk_by = reg_linked[ty][n], reg_linked[by][n] end
        if group
            gl = gr_list[by][n]
            ng = length(gl)
            nsam_gr = Dict{Int, Array{Array{Int, 1}, 1}}()
        end

        for y in [ty, by]
            sda_factors[y] = Dict{String, factors}()
            ci_ie[y][n] = Dict{String, Tuple{Float64, Float64}}()
            ci_de[y][n] = Dict{String, Tuple{Float64, Float64}}()
            rl, hl, hhs = reg_list[y][n][:], hh_list[y][n], households[y][n]
            ie[y], de[y] = hh_ie[y][n], hh_de[y][n]
            nh, nr = length(hl), length(rl)

            etab, cmat[y] = exp_table[y][n], conc_mat_wgh[y][n]
            ft = factors()
            if y != base_year
                convertTable(y, n, by, mrioPath, t_bp = t_bp, t_tax = t_tax, t_sub = t_sub, v_bp = v_bp, y_bp = y_bp)
                mrio = mrio_tabs_conv[y][n]
                ft.l, lx = calculateLeontief(mrio)
                ft.f = calculateIntensity(mrio, lx)
            elseif y == base_year
                if !haskey(l_factor, y)
                    l_factor[y], lx = calculateLeontief(mrio_tabs[y])
                    f_factor[y] = calculateIntensity(mrio_tabs[y], lx)
                    lx_by = []
                end
                mrio, ft.l, ft.f = mrio_tabs[y], l_factor[y], f_factor[y]
            end
            sda_factors[y][n] = ft

            nt = size(mrio.t, 1)
            n_nt = nr * n_gr
            append!(reg_list[y][n], [r * "_CF" * cf_intv_lb[gi] for r in rl, gi = 1:n_cf])
            append!(reg_list[y][n], [r * "_inc" * inc_intv_lb[gi] for r in rl, gi = 1:n_inc])
            append!(reg_list[y][n], [r * "_CF" * cf_bndr_lb[gi] for r in rl, gi = 1:n_cfb])
            append!(reg_list[y][n], [r * "_inc" * inc_bndr_lb[gi] for r in rl, gi = 1:n_incb])

            ft_p[y] = zeros(Float64, n_nt)
            idx_ls[y] = [Array{Int, 1}() for i = 1:nr]
            nsam[y] = zeros(Int, n_nt)
            ie_vals[y], de_vals[y] = [zeros(Float64, 0) for i=1:n_nt], [zeros(Float64, 0) for i=1:n_nt]
            ieByNat[y][n], deByNat[y][n] = zeros(Float64, n_nt), zeros(Float64, n_nt)
            if group; nsam_gr[y] = Array{Array{Int, 1}, 1}()end

            intv_dataset, bndr_dataset = [], []
            if n_cf > 0; push!(intv_dataset, (n_cf, hpos_cf[y][n], cf_intv)) end
            if n_inc > 0; push!(intv_dataset, (n_inc, hpos_inc[y][n], inc_intv)) end
            if n_cfb > 0; push!(bndr_dataset, (n_cfb, Dict(cat_hhl[y][n] .=> hh_cf[y][n][:,end]), cf_bndr)) end
            if n_incb > 0; push!(bndr_dataset, (n_incb, Dict(h => hhs[h].totinc for h in hl), inc_bndr)) end

            wg_reg[y] = [hhs[h].popwgh for h in hl]
            wg_hhs[y] = wg_reg[y] .* [hhs[h].size for h in hl]
            ie_wg[y] = ie[y] .* wg_reg[y]
            de_wg[y] = de[y] .* wg_reg[y]
            etb_wg[y] = wg_reg[y] .* etab

            for ri = 1:nr
                r = rl[ri]

                if region == "district"
                    idxs = (rc_mode ? filter(x -> r_lnk[hhs[hl[x]].district] == r, 1:nh) : filter(x -> hhs[hl[x]].district == r, 1:nh))
                elseif region == "province"
                    idxs = (rc_mode ? filter(x -> r_lnk[hhs[hl[x]].province] == r, 1:nh) : filter(x -> hhs[hl[x]].province == r, 1:nh))
                end
                idx_ls[y][ri] = idxs

                for (n_lv, hpos, intv) in intv_dataset
                    push!(idx_ls[y], filter(x -> hpos[hl[x]] < intv[1], idxs))
                    for i = 1:n_lv-2; push!(idx_ls[y], filter(x -> intv[i] <= hpos[hl[x]] < intv[i+1], idxs)) end
                    push!(idx_ls[y], filter(x -> intv[end-1] <= hpos[hl[x]], idxs))
                end
                for (n_bd, hh_val, bndr) in bndr_dataset
                    if bndr_mode == "percap"
                        push!(idx_ls[y], filter(x -> hh_val[hl[x]] / hhs[hl[x]].size < bndr[2], idxs))
                        for i = 2:n_bd-1; push!(idx_ls[y], filter(x -> bndr[i] <= hh_val[hl[x]] / hhs[hl[x]].size < bndr[i+1], idxs)) end
                        push!(idx_ls[y], filter(x -> hh_val[hl[x]] / hhs[hl[x]].size >= bndr[end], idxs))
                    elseif bndr_mode == "hhs"
                        push!(idx_ls[y], filter(x -> hh_val[hl[x]] < bndr[2], idxs))
                        for i = 2:n_bd-1; push!(idx_ls[y], filter(x -> bndr[i] <= hh_val[hl[x]] < bndr[i+1], idxs)) end
                        push!(idx_ls[y], filter(x -> hh_val[hl[x]] >= bndr[end], idxs))
                    end
                end

                for gri = 1:n_gr
                    ri = (gri == 1 ? ri : nr + (n_gr - 1) * (ri - 1) + gri-1)
                    ft_p[y][ri] = (gri == 1 ? pops[y][n][r] : sum(wg_hhs[y][idx_ls[y][ri]]))
                    nsam[y][ri] = (resample_size == 0 ? length(idx_ls[y][ri]) : resample_size)
                end
            end
            ie_prv_l[y], ie_prv_u[y] = zeros(Float64, n_nt), zeros(Float64, n_nt)
        end
        samples_gr[n] = nsam

        cepc_vals, cspf_vals = [zeros(Float64, 0) for i=1:n_nt], [zeros(Float64, 0) for i=1:n_nt]
        fls = []
        cepc_prv_l, cepc_prv_u, cspf_prv_l, cspf_prv_u  = zeros(Float64, n_nt), zeros(Float64, n_nt), zeros(Float64, n_nt), zeros(Float64, n_nt)
        er = 1.0
        er_c = er_limit

        i = 0
        while (er > er_c && i < iter)
            i += 1
            for y in [ty, by]
                ft = sda_factors[y][n]
                rl = reg_list[y][n]
                nr, nt = length(rl), size(ft.l, 1)

                ft_de = zeros(nr)
                if mode == pt_mode; ft_cepc, ft_cspf, ft_de = zeros(nr), zeros(nt, nr), zeros(nr) end

                for ri = 1:nr
                    r = rl[ri]
                    if replacement; re_idx = [trunc(Int, nsam[y][ri] * rand())+1 for x = 1:nsam[y][ri]]
                    else re_idx = sortperm([rand() for x = 1:nsam[y][ri]])
                    end

                    id = idx_ls[y][ri][re_idx]
                    if length(id) > 0
                        ws = sum(wg_hhs[y][id])
                        etw = etb_wg[y][id, :]

                        push!(ie_vals[y][ri], sum(ie_wg[y][id]) / ws * ft_p[y][ri])
                        push!(de_vals[y][ri], sum(de_wg[y][id]) / ws * ft_p[y][ri])

                        if mode == pt_mode
                            et_sum = sum(etw, dims=1)
                            ce_tot = sum(et_sum)
                            ce_pf = et_sum ./ ce_tot
                            ft_cepc[ri] = ce_tot / ws
                            ft_cspf[:,ri] = cmat[y] * ce_pf'
                        end
                        ft_de[ri] = de_vals[y][ri][i]
                    end
                end

                if mode == pt_mode; ft.p, ft.cepc, ft.cspf, ft.de = ft_p[y], ft_cepc, ft_cspf, ft_de end
            end

            prepareDeltaFactors(ty, by, nation = n, mode = mode, reuse = reuse)
            fls = structuralAnalysis(ty, by, n, mode = mode, fl_mats = fls, reuse = reuse)[2]

            for ri = 1:nr
                push!(cepc_vals[ri], deltas[(ty, by)][n][rl[ri]][4])
                push!(cspf_vals[ri], deltas[(ty, by)][n][rl[ri]][5])
            end

            sam_chk = [length(idx_ls[ty][ri]) > 0 && length(idx_ls[by][ri]) > 0 for ri = 1:nr]
            if i >= iter_min && i % er_chk_iter == 0
                li, ui = trunc(Int, (1 - ci_rate) / 2 * i) + 1, trunc(Int, ((1 - ci_rate) / 2 + ci_rate) * i) + 1
                current_vals, previous_vals = [], []
                for ri = 1:nr
                    r = rl[ri]
                    if sam_chk[ri]
                        for y in [ty, by]
                            sort!(ie_vals[y][ri])
                            sort!(de_vals[y][ri])
                            ci_ie[y][n][r] = (ie_vals[y][ri][li], ie_vals[y][ri][ui])
                            ci_de[y][n][r] = (de_vals[y][ri][li], de_vals[y][ri][ui])
                        end
                        sort!(cepc_vals[ri])
                        sort!(cspf_vals[ri])
                        ci_sda[(ty,by)][n][r] = [(cepc_vals[ri][li], cepc_vals[ri][ui]), (cspf_vals[ri][li], cspf_vals[ri][ui])]

                        append!(current_vals, [ie_vals[ty][ri][li], ie_vals[ty][ri][ui], ie_vals[by][ri][li], ie_vals[by][ri][ui], cepc_vals[ri][li], cepc_vals[ri][ui], cspf_vals[ri][li], cspf_vals[ri][ui]])
                        append!(previous_vals, [ie_prv_l[ty][ri], ie_prv_u[ty][ri], ie_prv_l[by][ri], ie_prv_u[by][ri], cepc_prv_l[ri], cepc_prv_u[ri], cspf_prv_l[ri], cspf_prv_u[ri]])
                    else
                        for y in [ty, by]; ci_ie[y][n][r], ci_de[y][n][r] = (0, 0), (0, 0) end
                        ci_sda[(ty,by)][n][r] = [(0, 0), (0, 0)]
                        append!(current_vals, [0, 0, 0, 0, 0, 0, 0, 0])
                        append!(previous_vals, [0, 0, 0, 0, 0, 0, 0, 0])
                    end
                end
                ers = abs.((current_vals - previous_vals) ./ previous_vals)
                ers[isnan.(ers)] .= 0
                ers[isinf.(ers)] .= 0
                # er = (all(x -> x == 0, ers) ? 1.0 : maximum(ers))
                er = maximum(ers)

                for ri = 1:nr
                    r = rl[ri]
                    for y in [ty, by]; ie_prv_l[y][ri], ie_prv_u[y][ri] = ci_ie[y][n][r] end
                    cepc_prv_l[ri], cepc_prv_u[ri] = ci_sda[(ty,by)][n][r][1]
                    cspf_prv_l[ri], cspf_prv_u[ri] = ci_sda[(ty,by)][n][r][2]
                end
            end
            if visible_iter > 0 && i % visible_iter == 0; print(" ", i) end
        end

        for ri = 1:nr, y in [ty, by]
            idxs = idx_ls[y][ri]
            wg_sum = sum(wg_hhs[y][idxs])
            ieByNat[y][n][ri] = sum(ie_wg[y][idxs]) / wg_sum * ft_p[y][ri]
            deByNat[y][n][ri] = sum(de_wg[y][idxs]) / wg_sum * ft_p[y][ri]
        end

        if visible
            elap = floor(Int, time() - st); (eMin, eSec) = fldmod(elap, 60); (eHr, eMin) = fldmod(eMin, 60)
            println(eHr,":",eMin,":",eSec," elapsed,\t", i, " iterations")
        end
    end
end

function printConfidenceIntervals(year, outputFile, nation = []; ci_rate = 0.95, region = "district", exp_mode = "gis")

    global nat_list, reg_list, hh_list, pops, ci_ie, ci_de, ieByNat, deByNat, hh_ie, hh_de, reg_linked
    if isa(year, Number); year = [year] end
    nats = (length(nation) == 0 ? nat_list : (isa(nation, String) ? [nation] : nation))
    rc_mode = (lowercase(exp_mode) == "gis")

    low_lab, upp_lab = string((1 - ci_rate) / 2), string((1 - ci_rate) / 2 + ci_rate)

    f = open(outputFile, "w")
    print(f, "Year\tNation\tRegion\tSamples\t")
    print(f, "Overall_IE\tIE_CI_", low_lab, "\tIE_CI_", upp_lab, "\tOverall_DE\tDE_CI_", low_lab, "\tDE_CI_", upp_lab)
    print(f, "\tIE_per_capita\tDE_per_capita\tPopulation\tTotal_weight")
    println(f)

    for y in year, n in nats
        rl, hl, hhs, nh = reg_list[y][n], hh_list[y][n], households[y][n], length(hh_list[y][n])
        if rc_mode; r_lnk = reg_linked[y][n] end
        ie, de = vec(sum(hh_ie[y][n], dims=1)), vec(sum(hh_de[y][n], dims=1))

        hhs_wgs = [hhs[h].popwgh for h in hl]
        hhs_wsz = hhs_wgs .* [hhs[h].size for h in hl]
        iew = ie .* hhs_wgs
        dew = de .* hhs_wgs

        for ri = 1:length(rl)
            r = rl[ri]
            p_reg = pops[y][n][r]

            if region == "district"
                idxs = (rc_mode ? filter(x -> r_lnk[hhs[hl[x]].district] == r, 1:nh) : filter(x -> hhs[hl[x]].district == r, 1:nh))
            elseif region == "province"
                idxs = (rc_mode ? filter(x -> r_lnk[hhs[hl[x]].province] == r, 1:nh) : filter(x -> hhs[hl[x]].province == r, 1:nh))
            end

            wg_reg = hhs_wgs[idxs]
            wg_sum = sum(hhs_wsz[idxs])

            print(f, y, "\t", n, "\t", r, "\t", length(idxs))
            print(f, "\t", ieByNat[y][n][ri], "\t", ci_ie[y][n][r][1], "\t", ci_ie[y][n][r][2])
            print(f, "\t", deByNat[y][n][ri], "\t", ci_de[y][n][r][1], "\t", ci_de[y][n][r][2])
            print(f, "\t", sum(iew[idxs]) / wg_sum, "\t", sum(dew[idxs]) / wg_sum, "\t", p_reg, "\t", wg_sum)
            println(f)
        end
    end
    close(f)
end

function printSdaCI_title(target_year, base_year, outputFile; ci_rate = 0.95, mode = "penta")

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")

    ty, by = target_year, base_year
    ll, ul = round((1 - ci_rate) / 2, digits = 3), round((1 - ci_rate) / 2 + ci_rate, digits = 3)

    print(f, "Nation\tRegion\t", ty,"_Samples\t", by,"_Samples")
    print(f, "\t", ty,"_IE\t", ty,"_IE_CI_",ll, "\t", ty,"_IE_CI_",ul, "\t", ty,"_DE\t", ty,"_DE_CI_",ll, "\t", ty,"_DE_CI_",ul)
    print(f, "\t", by,"_IE\t", by,"_IE_CI_",ll, "\t", by,"_IE_CI_",ul, "\t", by,"_DE\t", by,"_DE_CI_",ll, "\t", by,"_DE_CI_",ul)
    print(f, "\t", ty, "_IE_pc\t", ty, "_DE_pc\t", ty,"_Population")
    print(f, "\t", by, "_IE_pc\t", by, "_DE_pc\t", by,"_Population")
    if mode == "penta"; print(f, "\tCEPC_CI_", ll, "\tCEPC_CI_", ul, "\tCSPF_CI_", ll, "\tCSPF_CI_", ul) end
    println(f)

    close(f)
end

function printSdaCI_values(target_year, base_year, outputFile, nation = []; ci_rate = 0.95, mode = "penta")

    global nat_list, reg_list, hh_list, pops
    global ci_ie, ci_de, ci_sda, ieByNat, deByNat, hh_ie, hh_de, sda_factors, samples_gr

    ty, by = target_year, base_year
    ll, ul = round(((1 - ci_rate) / 2), digits = 3), round(((1 - ci_rate) / 2 + ci_rate), digits = 3)

    if length(nation) == 0; nats = nat_list
    elseif isa(nation, String); nats = [nation]
    elseif isa(nation, Array{String, 1}); nats = nation
    end

    f = open(outputFile, "a")

    for n in nats
        rl_by = reg_list[by][n]
        ie_ty, de_ty = vec(sum(hh_ie[ty][n], dims=1)), vec(sum(hh_de[ty][n], dims=1))
        ie_by, de_by = vec(sum(hh_ie[by][n], dims=1)), vec(sum(hh_de[by][n], dims=1))
        ft_ty , ft_by = sda_factors[ty][n], sda_factors[by][n]
        sam_ty, sam_by = samples_gr[n][ty], samples_gr[n][by]

        for ri = 1:length(rl_by)
            r = rl_by[ri]
            p_reg_ty, p_reg_by = ft_ty.p[ri], ft_by.p[ri]

            print(f, n, "\t", r, "\t", sam_ty[ri], "\t", sam_by[ri])
            for y in [ty, by]
                print(f, "\t", ieByNat[y][n][ri], "\t", ci_ie[y][n][r][1], "\t", ci_ie[y][n][r][2])
                print(f, "\t", deByNat[y][n][ri], "\t", ci_de[y][n][r][1], "\t", ci_de[y][n][r][2])
            end
            print(f, "\t", ieByNat[ty][n][ri] / ft_ty.p[ri], "\t", deByNat[ty][n][ri] / ft_ty.p[ri], "\t", p_reg_ty)
            print(f, "\t", ieByNat[by][n][ri] / ft_by.p[ri], "\t", deByNat[by][n][ri] / ft_by.p[ri], "\t", p_reg_by)
            print(f, "\t", ci_sda[(ty,by)][n][r][1][1], "\t", ci_sda[(ty,by)][n][r][1][2])
            print(f, "\t", ci_sda[(ty,by)][n][r][2][1], "\t", ci_sda[(ty,by)][n][r][2][2])
            println(f)
        end
    end
    close(f)
end

function printDeltaTitle(outputFile; cf_print = true, st_print = true, mode = "penta")

    global cat_list

    vs = getValueSeparator(outputFile)
    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    print(f, "Target_year", vs, "Base_year", vs, "Nation", vs, "Region")
    print(f, vs, "f", vs, "L", vs, "p")
    if mode == "penta"; print(f, vs, "exp_pc", vs, "exp_profile")
    elseif mode == "hexa"; print(f, vs, "exp_pc", vs, "exp_profile", vs, "exp_cat")
    elseif mode == "categorized"
        for c in cat_list; print(f, vs, c, "_exp_pc") end
        print(f, vs, "exp_profile")
    else println("SDA mode error: ", mode)
    end
    print(f, vs, "de", vs, "total_delta")
    if cf_print; print(f, vs, "Measured_delta", vs, "Target_year_IE", vs, "Base_year_IE", vs, "Target_year_DE", vs, "Base_year_DE") end
    if st_print; print(f, vs, "Population_target", vs, "Population_base", vs, "Total_exp_target", vs, "Total_exp_base") end
    println(f)
    close(f)
end

function printDeltaValues(outputFile, nation = ""; cf_print = true, st_print = true, mode = "penta")

    global reg_list, deltas, ieByNat, deByNat, cfByNat, popByNat, expPcByNat, cat_list

    vs = getValueSeparator(outputFile)
    f = open(outputFile, "a")
    for yrs in sort(collect(keys(deltas)))
        if nation == ""; nats = sort(collect(keys(deltas[yrs])))
        elseif isa(nation, String); nats = [nation]
        elseif isa(nation, Array{String, 1}); nats = nation
        end
        for n in nats
            if cf_print
                t_ie, t_de = ieByNat[yrs[1]][n], deByNat[yrs[1]][n]
                b_ie, b_de = ieByNat[yrs[2]][n], deByNat[yrs[2]][n]
                t_p, t_epc = popByNat[yrs[1]][n], expPcByNat[yrs[1]][n]
                b_p, b_epc = popByNat[yrs[2]][n], expPcByNat[yrs[2]][n]
            end
            nn = length(reg_list[yrs[1]][n])
            for i = 1:nn
                r = reg_list[yrs[1]][n][i]
                print(f, yrs[1], vs, yrs[2], vs, n, vs, r)
                for d in deltas[yrs][n][r]; print(f, vs, d) end
                print(f, vs, t_de[i] - b_de[i])
                print(f, vs, sum(deltas[yrs][n][r]))
                if cf_print; print(f, vs, t_ie[i] - b_ie[i], vs, t_ie[i], vs, b_ie[i], vs, t_de[i], vs, b_de[i]) end
                if st_print; print(f, vs, t_p[i], vs, b_p[i], vs, t_epc[i], vs, b_epc[i]) end
                println(f)
            end
        end
    end
    close(f)
end

function clearFactors(; year = 0, nation = "")

    global sda_factors, hh_list, households, exp_table, mrio_tabs_conv, conc_mat_wgh, hh_de

    if year == 0; yrs = sort(collect(keys(sda_factors))) else yrs = [year] end
    for y in yrs
        if length(nation) == 0; nats = sort(collect(keys(sda_factors[y])))
        elseif isa(nation, String); nats = [nation]
        elseif isa(nation, Array{String, 1}); nats = nation
        end
        for n in nats
            hh_list[y][n] = Array{String, 1}()
            households[y][n] = Dict{String, mdr.household}()
            exp_table[y][n] = Array{Float64, 2}(undef, 0, 0)
            conc_mat_wgh[y][n] = Array{Float64, 2}(undef, 0, 0)
            sda_factors[y][n] = factors()
            dltByNat[n] = Dict{Int, Any}()
            if haskey(hh_ie, y); hh_ie[y][n] = Array{Float64, 2}(undef, 0, 0)  end
            if haskey(hh_de, y); hh_de[y][n] = Array{Float64, 2}(undef, 0, 0)  end
            if haskey(mrio_tabs_conv, y); mrio_tabs_conv[y][n] = ee.tables() end
        end
    end
end

end
