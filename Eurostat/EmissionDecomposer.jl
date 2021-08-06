module EmissionDecomposer

# Developed date: 27. Jul. 2021
# Last modified date: 5. Aug. 2021
# Subject: Decompose EU households' carbon footprints
# Description: Process for Input-Output Structural Decomposition Analysis
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

include("MicroDataReader.jl")
include("EmissionEstimator.jl")
include("EmissionCategorizer.jl")

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

    factors(nr, nt) =  new(zeros(nt), zeros(nt, nt), zeros(nr), zeros(nr), zeros(nt, nr), zeros(nr))
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
global exp_table_conv = Dict{Int, Dict{String, Array{Float64, 2}}}()        # Base-year price converted household expenditure table: {year, {nation, {hhid, category}}}

global nuts = Dict{Int, Dict{String, String}}()                     # NUTS: {year, {code, label}}
global nuts_list = Dict{Int, Array{String, 1}}()                    # NUTS code list: {year, {NUTS_code}}
global pops = Dict{Int, Dict{String, Float64}}()                    # Population: {year, {NUTS_code, population}}
global pop_list = Dict{Int, Dict{String, Dict{String, Float64}}}()  # Population list: {year, {nation_code, {NUTS_code, population}}}
global pop_label = Dict{Int, Dict{String, String}}()                # populaton NUTS label: {year, {NUTS_code, NUTS_label}}

global cpi_list = Dict{Int, Dict{String, Array{String, 1}}}()       # Consumption price indexes: {year, {nation, {COICOP_category}}}
global cpis = Dict{Int, Dict{String, Dict{String, Float64}}}()      # Consumption price indexes: {year, {nation, {COICOP_category, CPI}}}
global scl_rate = Dict{Int, Dict{String, Dict{String, Float64}}}()  # CPI scaling rate: {year, {nation, {HBS code, rate}}}
global conc_mat = Dict{Int, Dict{String, Array{Float64, 2}}}()      # Concordance matrix: {year, {nation, {Eora sectors, Nation sectors}}}
global mrio_tabs = Dict{Int, ee.tables}()                           # MRIO tables: {Year, MRIO tables (t, v ,y , q)}
global conc_mat_conv = Dict{Int, Dict{String, Array{Float64, 2}}}() # Base-year price converted concordance matrix: {year, {nation, {Eora sectors, Nation sectors}}}
global mrio_tabs_conv = Dict{Int, Dict{String, ee.tables}}()        # Base-year price converted MRIO tables: {Year, {natoin, MRIO tables (t, v ,y , q)}}

global nt_wgh = Dict{Int, Dict{String, Float64}}()                  # hhid's NUTS weight: {year, {hhid, weight}}
global di_emiss = Dict{Int, Dict{String, Array{Float64, 2}}}()      # direct carbon emission: {year, {nation, {table}}}

global sda_factors = Dict{Int, Dict{String, factors}}()             # SDA factors: {year, {nation, factors}}

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

function importData(; hh_data::Module, mrio_data::Module, cat_data::Module, nations = [])

    global yr_list, nat_list, nat_name = hh_data.year_list, hh_data.nations, hh_data.nationNames
    global hh_list, households, exp_table, scl_rate = hh_data.hhsList, hh_data.mdata, hh_data.expTable, hh_data.sclRate
    global mrio_tabs, sc_list = mrio_data.mTables, mrio_data.sec
    global nt_wgh, di_emiss = cat_data.wghNuts, cat_data.directCE
    global cat_list, nuts, pops, pop_list, pop_label = cat_data.catList, cat_data.nuts, cat_data.pop, cat_data.popList, cat_data.poplb
    global nuts_list

    if length(nations)>0; nat_list = nations end

    for y in yr_list
        nuts_list[y] = Array{String, 1}()
        for n in nat_list; append!(nuts_list[y], filter(x->x in cat_data.giscdlist[y], cat_data.nutsList[y][n])) end
        # exp_table_conv[y] = Dict{String, Array{Float64, 2}}()
    end
end

function storeNutsWeight()

    global yr_list, nat_list, households, hh_list, nt_wgh

    for y in yr_list, hh in collect(keys(nt_wgh[y]))
        n, h = hh[1:2], hh[4:end]
        if h in hh_list[y][n]; households[y][n][h].weight_nt = nt_wgh[y][hh]
        else println(h, " household is not in the ", y, " year's list")
        end
    end
end

function storeConcMat(year, nation, concMat)

    global conc_mat

    if !haskey(conc_mat, year); conc_mat[year] = Dict{String, Array{Float64, 2}}() end
    conc_mat[year][nation] = concMat
end

function convertTable(year, nation)

    global sc_list, scl_rate, conc_mat, mrio_tabs, mrio_tabs_conv, exp_table, exp_table_conv
    sclr, mrio, expt = scl_rate[year][nation], mrio_tabs[year], exp_table[year][nation]

    if !haskey(mrio_tabs_conv, year); mrio_tabs_conv[year] = Dict{String, ee.tables}() end
    cvr_conc = [sclr[c] for c in sc_list[year]]
    cmat = conc_mat[year][nation]

    # exp_table_conv[year][nation] =  expt .* cvr_conc'

    cvr_mrio = sum(cmat .* cvr_conc', dims = 2) ./ sum(cmat, dims = 2)
    mrio_conv = ee.tables(year, size(mrio.t, 1), size(mrio.v, 1), size(mrio.y, 2), size(mrio.q, 1))
    mrio_conv.t = mrio.t .* cvr_mrio'       # notice about the converting direction of 'cvr_mrio': column-wise
    mrio_conv.v = mrio.v .* cvr_mrio'       # notice about the converting direction of 'cvr_mrio': column-wise
    mrio_conv.y = mrio.y .* cvr_mrio        # notice about the converting direction of 'cvr_mrio': row-wise
    mrio_conv.q = mrio.q
    mrio_tabs_conv[year][nation] = mrio_conv
end

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

function prepareSDA(year, baseYear)

    # f * L * y + DE
    # f * L * [conc * hbs_region] + DE
    # f * L * [conc * tot_ce_region * hbs_profile] + DE
    # f * L * tot_ce_region * [conc * hbs_profile] + DE
    # f * L * p * tot_ce_pc * [con * hbs_profile] + DE

    global mrio_tabs, mrio_tabs_conv, conc_mat, sda_factors, di_emiss
    global nat_list, nuts_list, hh_list, pop_list
    if isa(year, Number); year = [year] end

    for y in year, n in nat_list
        etab, cmat = exp_table[y][n], conc_mat[y][n]
        if y == baseYear; mrio = mrio_tabs[y]
        else convertTable(y,n); mrio = mrio_tabs_conv[y][n]
        end
        hhs, de = hh_list[y][n], di_emiss[y][n]
        nr, nt = length(nuts_list[y]), size(mrio.t, 1)
        ft = factors(nr, nt)
        ft.l = calculateLeontief(mrio_tabs[y])
        ft.f = calculateIntensity(mrio)

        regs = filter(x->x[1:2] == n, nuts_list[y])
        for r in regs
            p_reg = pop_list[y][n][r]
            ri = findfirst(x->x==r, nuts_list[y])
            idxs = [findfirst(x->x==hh, hhs) for hh in filter(x->households[y][n][x].nuts1 == r, hh_list[y][n])]

            wg_reg = [households[y][n][h].weight_nt for h in hh_list[y][n][idxs]]

            etb_wg = wg_reg .* etab[idxs, :]
            ce_tot = sum(etb_wg)
            ce_pf = sum(etb_wg, dims=1) ./ ce_tot

            ft.p[ri], ft.cepc[ri] = p_reg, ce_tot/p_reg
            ft.cspf[:,ri] = cmat * ce_pf'
            ft.de[ri] = (sum(de[:, idxs], dims=1) * wg_reg)[1]
        end

        if !haskey(sda_factors, y); sda_factors[y] = Dict{String, factors}() end
        sda_factors[y][n] = ft
    end
end

function testSDA(year, output)

    global sda_factors, nat_list, nuts_list
    cf = Dict{String, Array{Float64, 2}}()

    for n in nat_list, nt in filter(x->x[1:2] == n, nuts_list[year])
        ft = sda_factors[year][n]
        cf[n] = sum(ft.f .* ft.l * ((ft.p .* ft.cepc)' .* ft.cspf), dims=1) + ft.de
    end

    for n in nat_list
        nts = filter(x->x[1:2] == n, nuts_list[year])
        for i=1:length(cf[n]); println(n ,"\t", nts[i], "\t", cf[n], "\t", cf[n]./sda_factors[year][n].p) end
        println(n, "\t", sda_factors[year][n].de, "\t", sda_factors[year][n].de./sda_factors[year][n].p)
        println(n, "\t", sda_factors[year][n].p)
        println(n, "\t", sda_factors[year][n].cepc)
    end

end

end
