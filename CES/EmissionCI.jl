# SPDX-FileCopyrightText: © 2022 Jemyung Lee <jemyung81@gmail.com>
# SPDX-License-Identifier: GPL-3.0

module EmissionCI

# Developed date: 10. Mar. 2022
# Last modified date: 8. Aug. 2023
# Subject: CI estimation of household CF
# Description: Estimate confidence intervals of household carbon footprint assessed from
#               Customer Expenditure Survey (CES) or Household Budget Survey (HBS) micro-data.
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

global yr_list = Array{Int, 1}()                    # year list: {YYYY}
global nat_list = Dict{Int, Array{String, 1}}()     # nation list: {year, {A2}}
global nat_name = Dict{String, String}()            # nation names: {Nation code, Name}
global cat_list = Array{String, 1}()                # category list
global pr_unts = Dict("day" => 1, "week" => 7, "month" => 30, "year" => 365)
global web_index = Array{Tuple{String, String}, 1}()

global sc_list = Dict{Int, Dict{String, Array{String, 1}}}()        # HBS commodity code list: {year, {nation A3, {code}}}
global hh_list = Dict{Int, Dict{String, Array{String, 1}}}()        # Household ID: {year, {nation, {hhid}}}
global gr_list = Dict{Int, Dict{String, Array{String, 1}}}()        # group (survey type) list: {year, {nation A3, {group}}}
global households = Dict{Int, Dict{String, Dict{String, mdr.household}}}()  # household dict: {year, {nation, {hhid, household}}}
global exp_table = Dict{Int, Dict{String, Array{Float64, 2}}}()     # household expenditure table: {year, {nation, {hhid, category}}}
global sc_cat = Dict{Int, Dict{String, Dict{String, String}}}()     # CES/HBS sector-category link dict: {year, {nation, {sector_code, category}}}
global conc_mat = Dict{Int, Array{Float64, 2}}()                    # Concordance matrix: {year, {Eora sectors, CES/HBS sectors}}

global hh_cf = Dict{Int, Dict{String, Array{Float64, 2}}}()         # categozied carbon footprint by household: {year, {nation, {hhid, category}}}
global cat_hhl = Dict{Int, Dict{String, Array{String, 1}}}()        # Household ID: {year, {nation, {hhid}}}

global pops = Dict{Int, Dict{String, Dict{String, Float64}}}()      # population: {year, {nation, {region_code, population}}}

global in_emiss = Dict{Int, Dict{String, Array{Float64, 2}}}()      # indirect carbon emission: {year, {nation, {CES/HBS sector, household}}}
global di_emiss = Dict{Int, Dict{String, Array{Float64, 2}}}()      # direct carbon emission: {year, {nation, {CES/HBS sector, household}}}

global ieByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()       # indirect CF by region: {year, {nation, {region}}}
global deByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()       # direct CF by region: {year, {nation, {region}}}
global cfByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()       # total CF by region: {year, {nation, {region}}}
global cfByReg = Dict{Int, Dict{String, Array{Float64, 2}}}()       # categozied carbon footprint per capita by region: {year, {nation, {region, category}}}

global popByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()      # population by nation, NUTS: {year, {nation, {region}}}
global expPcByNat = Dict{Int, Dict{String, Array{Float64, 1}}}()    # total expenditures by nation, NUTS: {year, {nation, {region}}}

global ci_ie = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}() # confidence intervals of indirect emission: {year, {nation, {region, {lower, upper}}}}
global ci_de = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}() # confidence intervals of direct emission: {year, {nation, {region, {lower, upper}}}}
global ci_cf = Dict{Int, Dict{String, Dict{String, Tuple{Float64, Float64}}}}() # confidence intervals of CF: {year, {nation, {region, {lower, upper}}}}
global ci_cfpc = Dict{Int, Dict{String, Dict{String, Array{Tuple{Float64, Float64}, 1}}}}() # confidence intervals of CF/capita by category: {year, {nation, {region, {(by category) {lower, upper}}}}}

global cat_reg = Dict{Int, Dict{String, Dict{String, String}}}()    # region code-name: {year, {nation A3, {code, region}}}
global cat_prov = Dict{Int, Dict{String, Array{String, 1}}}()       # province code list: {year, {nation A3, {code}}}
global cat_dist = Dict{Int, Dict{String, Array{String, 1}}}()       # district code list: {year, {nation A3, {code}}}

global gis_reg_id = Dict{Int, Dict{String, Dict{String,String}}}()  # GIS region ID: {year, {nation, {region_code, region_ID}}}
global gis_reg_list = Dict{Int, Dict{String, Array{String, 1}}}()   # GIS region list: {year, {nation, {region code}}}
global gis_reg_conc = Dict{Int, Dict{String, Array{Float64, 2}}}()  # GIS-CES/HBS region concordance weight: {year, {nation, {gis_code, ces/hbs_code}}}
global gis_reg_dist= Dict{Int, Dict{String, Dict{String,String}}}() # GIS ID corresponds CES/HBS region: {year, {nation, {CES/HBS region, GIS id}}}

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

function readCityFileSector(webIndexFile)

    global web_index

    f = open(webIndexFile)
    for l in eachline(f)
        s = string.(split(l, '\t'))
        push!(web_index, (s[1], s[2]))
    end
    close(f)
end

function importData(; hh_data::Module, mrio_data::Module, cat_data::Module, cat_filter = true)

    global yr_list, hh_list, households, exp_table, gr_list = cat_data.yr_list, hh_data.hh_list, hh_data.households, hh_data.expMatrix, hh_data.gr_list
    global mrio_idxs, mrio_tabs, sc_list, conc_mat = mrio_data.ti, mrio_data.mTables, mrio_data.sc_list, mrio_data.concMat
    global in_emiss, di_emiss, hh_cf, cfByReg = cat_data.indirectCE, cat_data.directCE, cat_data.cfHHs, cat_data.cfReg
    global cat_list, cat_reg, cat_dist, cat_prov = cat_data.cat_list, cat_data.regions, cat_data.dist_list, cat_data.prov_list
    global sc_cat, pops = cat_data.sc_cat, cat_data.pops
    global gis_reg_list, gis_reg_id, gis_reg_conc, gis_reg_dist = cat_data.gisRegList, cat_data.gisRegID, cat_data.gisRegConc, cat_data.gisRegDist

    if cat_filter; filter!(x -> !(lowercase(x) in ["total", "all"]), cat_list) end
end

function estimateConfidenceIntervals(year, nation; iter = 10000, ci_rate = 0.95, resample_size = 0, replacement = true, boundary="district", group = false)
    # bootstrap method
    # ci_per: confidence interval percentage
    # replacement: [true] sampling with replacement
    # if resample_size: [0] resample_size = sample_size
    # group: [true] divide households by (survey) groups

    global hh_list, pops, cat_list, cat_dist, cat_prov, households, gr_list
    global ci_cf, ci_cfpc, cfByNat, hh_cf
    # global ci_ie, ci_de, in_emiss, di_emiss, ieByNat, deByNat

    y, n = year, nation

    if resample_size == 0; replacement = true end
    nc = length(cat_list)

    if !haskey(ci_ie, y); ci_ie[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
    if !haskey(ci_de, y); ci_de[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
    if !haskey(ci_cf, y); ci_cf[y] = Dict{String, Dict{String, Tuple{Float64, Float64}}}() end
    if !haskey(ci_cfpc, y); ci_cfpc[y] = Dict{String, Dict{String, Array{Tuple{Float64, Float64}, 1}}}() end
    # if !haskey(ieByNat, y); ieByNat[y] = Dict{String, Array{Float64, 1}}() end
    # if !haskey(deByNat, y); deByNat[y] = Dict{String, Array{Float64, 1}}() end
    if !haskey(cfByNat, y); cfByNat[y] = Dict{String, Array{Float64, 1}}() end

    if !haskey(ci_ie[y], n); ci_ie[y][n] = Dict{String, Tuple{Float64, Float64}}() end
    if !haskey(ci_de[y], n); ci_de[y][n] = Dict{String, Tuple{Float64, Float64}}() end
    if !haskey(ci_cf[y], n); ci_cf[y][n] = Dict{String, Tuple{Float64, Float64}}() end
    if !haskey(ci_cfpc[y], n); ci_cfpc[y][n] = Dict{String, Array{Tuple{Float64, Float64}, 1}}() end

    hhl, hhs = hh_list[y][n], households[y][n]
    hcf = hh_cf[y][n]
    nh = length(hhl)
    # ie, de = vec(sum(in_emiss[y][n], dims=1)), vec(sum(di_emiss[y][n], dims=1))

    if boundary == "district"
        reg_ls, nr = cat_dist[y][n], length(cat_dist[y][n])
        ridx = [filter(i->hhs[hhl[i]].district == d, 1:nh) for d in reg_ls]
    elseif boundary == "province"
        reg_ls, nr = cat_prov[y][n], length(cat_prov[y][n])
        ridx = [filter(i->hhs[hhl[i]].province == p, 1:nh) for p in reg_ls]
    end
    cfByNat[y][n] = zeros(Float64, nr)
    # ieByNat[y][n], deByNat[y][n] = zeros(Float64, nr), zeros(Float64, nr)

    if group
        grl = gr_list[y][n]
        ng = length(grl)
        ridx_gr = [[filter(i->hhs[hhl[i]].group == g, ridx[ri]) for g in grl] for ri = 1:nr]
    end

    wgs = [hhs[h].popwgh for h in hhl]
    hwgs = [hhs[h].size for h in hhl] .* wgs
    # iews = ie .* wgs
    # dews = de .* wgs
    hcfws = hcf .* wgs

    for ri = 1:nr
        r = reg_ls[ri]
        p_reg = pops[y][n][r]

        cf_vals = zeros(Float64, iter)
        # ie_vals, de_vals = zeros(Float64, iter), zeros(Float64, iter)
        cfpc_vals = [zeros(Float64, iter) for i = 1:nc]

        if !group
            wg_reg = wgs[ridx[ri]]
            hwg_reg = hwgs[ridx[ri]]
            # iew_reg = iews[ridx[ri]]
            # dew_reg = dews[ridx[ri]]
            hcfw_reg = hcfws[ridx[ri], :]
            nsam = (resample_size == 0 ? length(ridx[ri]) : resample_size)

            for i = 1:iter
                if replacement; re_idx = [trunc(Int, nsam * rand())+1 for x = 1:nsam]
                else re_idx = sortperm([rand() for x = 1:length(ridx[ri])])[1:nsam]
                end
                wg_sum = sum(hwg_reg[re_idx])
                # ie_vals[i] = sum(iew_reg[re_idx]) / wg_sum * p_reg
                # de_vals[i] = sum(dew_reg[re_idx]) / wg_sum * p_reg
                cf_vals[i] = sum(hcfw_reg[re_idx, end]) / wg_sum * p_reg
                cfpcs = vec(sum(hcfw_reg[re_idx, 1:nc], dims = 1)) / wg_sum
                for j = 1:nc; cfpc_vals[j][i] = cfpcs[j] end
            end
        elseif group
            wg_reg_gr = [wgs[igr] for igr in ridx_gr[ri]]
            hwg_reg_gr = [hwgs[igr] for igr in ridx_gr[ri]]
            # iew_reg_gr = [iews[igr] for igr in ridx_gr[ri]]
            # dew_reg_gr = [dews[igr] for igr in ridx_gr[ri]]
            hcfw_reg_gr = [hcfws[igr, :] for igr in ridx_gr[ri]]
            nsam_gr = [length(igr) for igr in ridx_gr[ri]]
            if resample_size > 0; nsam_gr *= resample_size / sum(nsam_gr) end

            for i = 1:iter
                if replacement; re_idx_gr = [[trunc(Int, ns_gr * rand())+1 for x = 1:ns_gr] for ns_gr in nsam_gr]
                else re_idx_gr = [sortperm([rand() for x = 1:length(ridx_gr[ri][gi])])[1:nsam_gr[gi]] for gi = 1:ng]
                end
                wg_sum_gr = [sum(hwg_reg_gr[gi][re_idx_gr[gi]]) for gi = 1:ng]

                # ie_vals[i] = sum([sum(iew_reg_gr[gi][re_idx_gr[gi]]) / wg_sum_gr[gi] for gi = 1:ng]) * p_reg
                # de_vals[i] = sum([sum(dew_reg_gr[gi][re_idx_gr[gi]]) / wg_sum_gr[gi] for gi = 1:ng]) * p_reg
                cf_vals[i] = sum([sum(hcfw_reg_gr[gi][re_idx_gr[gi], end]) / wg_sum_gr[gi] for gi = 1:ng]) * p_reg
                cfpcs = vec(sum([sum(hcfw_reg_gr[gi][re_idx_gr[gi], 1:nc], dims = 1) / wg_sum_gr[gi] for gi = 1:ng]))

                for j = 1:nc; cfpc_vals[j][i] = cfpcs[j] end
            end
        end

        # sort!(ie_vals)
        # sort!(de_vals)
        sort!(cf_vals)
        for i = 1:nc; sort!(cfpc_vals[i]) end

        l_idx, u_idx = trunc(Int, (1 - ci_rate) / 2 * iter) + 1, trunc(Int, ((1 - ci_rate) / 2 + ci_rate) * iter) + 1
        # ci_ie[y][n][r] = (ie_vals[l_idx], ie_vals[u_idx])
        # ci_de[y][n][r] = (de_vals[l_idx], de_vals[u_idx])
        ci_cf[y][n][r] = (cf_vals[l_idx], cf_vals[u_idx])
        ci_cfpc[y][n][r] = [(cfpc_vals[i][l_idx], cfpc_vals[i][u_idx]) for i = 1:nc]

        if !group
            wg_sum = sum(hwg_reg)
            # ieByNat[y][n][ri] = sum(iew_reg) / wg_sum * p_reg
            # deByNat[y][n][ri] = sum(dew_reg) / wg_sum * p_reg
            cfByNat[y][n][ri] = sum(hcfw_reg[:, end]) / wg_sum * p_reg
        elseif group
            wg_sum_gr = [sum(hwg_reg_gr[gi]) for gi = 1:ng]
            cfByNat[y][n][ri] = sum([sum(hcfw_reg_gr[gi][:, end]) / wg_sum_gr[gi] for gi = 1:ng]) * p_reg
        end
    end
end


# function printConfidenceIntervals(year, outputFile, nation = []; pop_dens = 0, ci_rate = 0.95)
#
#     global nat_list, nutsByNat, hh_list, pops, pop_list, pop_linked_cd, pops_ds, ci_ie, ci_de, ieByNat, deByNat, in_emiss, di_emiss
#     if isa(year, Number); year = [year] end
#     if length(nation) == 0; nats = nat_list else nats = nation end
#
#     dens_label = Dict(0 => "all", 1 => "densely", 2 => "inter", 3 => "sparsely")
#     low_lab, upp_lab = string((1 - ci_rate) / 2), string((1 - ci_rate) / 2 + ci_rate)
#
#     f = open(outputFile, "w")
#     print(f, "Year\tNation\tNUTS\tDensity\tSamples\t")
#     print(f, "Overall_IE\tIE_CI_", low_lab, "\tIE_CI_", upp_lab, "\tOverall_DE\tDE_CI_", low_lab, "\tDE_CI_", upp_lab)
#     print(f, "\tIE_per_capita\tDE_per_capita\tPopulation\tTotal_weight")
#     println(f)
#
#     for y in year
#         if length(nation) == 0; nats = nat_list[y] else nats = nation end
#         for n in nats
#             nts, hhs, nh = nutsByNat[y][n], hh_list[y][n], length(hh_list[y][n])
#             ie, de = vec(sum(in_emiss[y][n], dims=1)), vec(sum(di_emiss[y][n], dims=1))
#
#             for ri = 1:length(nts)
#                 r = nts[ri]
#                 # r_p = pop_linked_cd[y][r]
#                 # while r_p[end] == '0'; r_p = r_p[1:end-1] end
#                 # p_reg = pop_dens in [1,2,3] ? pops_ds[y][r_p][pop_dens] : pops[y][r_p]
#                 p_reg = (pop_dens in [1,2,3] ? pops_ds[y][r][pop_dens] : pop_list[y][n][r])
#
#                 idxs = filter(x -> households[y][n][hhs[x]].nuts1 == r, 1:nh)
#                 if pop_dens in [1,2,3]; filter!(x -> households[y][n][hhs[x]].popdens == pop_dens, idxs) end
#
#                 wg_reg = [households[y][n][h].weight_nt for h in hh_list[y][n][idxs]]
#                 wg_sum = sum(wg_reg .* [households[y][n][h].size for h in hh_list[y][n][idxs]])
#
#                 print(f, y, "\t", n, "\t", r, "\t", dens_label[pop_dens], "\t", length(idxs))
#                 print(f, "\t", ieByNat[y][n][ri], "\t", ci_ie[y][n][r][1], "\t", ci_ie[y][n][r][2])
#                 print(f, "\t", deByNat[y][n][ri], "\t", ci_de[y][n][r][1], "\t", ci_de[y][n][r][2])
#                 print(f, "\t", sum(ie[idxs] .* wg_reg) / wg_sum, "\t", sum(de[idxs] .* wg_reg) / wg_sum, "\t", p_reg, "\t", wg_sum)
#                 println(f)
#             end
#         end
#     end
#     close(f)
# end

function exportWebsiteCityFiles(year, nation, path, web_cat, cfav_file, cfac_file; boundary="district")

    global cat_list, cat_dist, cat_prov, pops, cfByNat, cfByReg, ci_cf, ci_cfpc
    global gis_reg_list, gis_reg_id, gis_reg_conc, web_index
    if isa(year, Number)
        if isa(cfav_file, String); cfav_file = Dict(year => cfav_file) end
        if isa(cfac_file, String); cfac_file = Dict(year => cfac_file) end
        year = [year]
    end
    if isa(nation, String); nats = [nation] elseif isa(nation, Array{String, 1}); nats = nation end

    web_cat_conc = Dict()
    for i = 1:length(cat_list)
        rx = r"[^a-zA-Z0-9]+"
        wc_lab = lowercase(replace(web_cat[i], rx => ""))
        web_cat_conc[web_cat[i]] = cat_list[findfirst(x -> occursin(wc_lab, lowercase(replace(x, rx => ""))), cat_list)]
    end

    cfav, cfac = Dict{Int, Dict{String, Array{String, 1}}}(), Dict{Int, Dict{String, Array{String, 1}}}()

    mkpath(path)
    nc = length(cat_list)

    reg_ls = Dict(y => Dict{String, Array{String, 1}}() for y in year)
    for y in year, n in nats
        if boundary == "district"; reg_ls[y][n] = cat_dist[y][n]
        elseif boundary == "province"; reg_ls[y][n] = cat_prov[y][n]
        end
    end

    regs = Array{String, 1}()
    for y in year, n in nats; append!(regs, [gis_reg_id[y][n][r] for r in gis_reg_list[y][n]]) end
    sort!(unique!(regs))

    strs = Dict{String, Dict{Int, String}}()

    for rg in regs
        city_file = path * rg * ".txt"
        strs[rg] = Dict{Int, String}()

        if isfile(city_file)
            f = open(city_file)
            ctitle = string.(split(readline(f), "\t"))
            widxs = [widx[1] for widx in web_index]

            if ctitle == widxs
                for l in eachline(f); strs[rg][parse(Int, string(split(l, "\t", limit = 2)[1]))] = l end
            else println("\nInconsistent web index: ", ctitle, " and ", widxs)
            end
            close(f)
        end

        # f = open(city_file, "w")
        # print(f, web_index[1][1])
        # for widx in web_index[2:end]; print(f, "\t", widx[1]) end
        # println(f)
        # close(f)
    end

    for y in year
        cfav[y], cfac[y] = Dict{String, Array{String, 1}}(), Dict{String, Array{String, 1}}()

        f = open(cfav_file[y])
        ctitle = string.(split(readline(f), ","))[2:end]
        cidx = [findfirst(x -> x == c, ctitle) for c in cat_list]
        push!(cidx, findfirst(x -> x == "Total", ctitle))
        for l in eachline(f)
            s = string.(split(l, ","))
            cfav[y][s[1]] = s[2:end][cidx]
        end
        close(f)

        f = open(cfac_file[y])
        ctitle = string.(split(readline(f), ","))[2:end]
        cidx = [findfirst(x -> x == c, ctitle) for c in cat_list]
        push!(cidx, findfirst(x -> x == "Total", ctitle))
        for l in eachline(f)
            s = string.(split(l, ","))
            cfac[y][s[1]] = s[2:end][cidx]
        end
        close(f)
    end

    for y in year, n in nats, ri = filter(x -> haskey(gis_reg_dist[y][n], reg_ls[y][n][x]) && gis_reg_dist[y][n][reg_ls[y][n][x]] in gis_reg_list[y][n], 1:length(reg_ls[y][n]))
        rg = reg_ls[y][n][ri]
        rid = gis_reg_id[y][n][gis_reg_dist[y][n][rg]]
        str = ""
        str *= string(y)
        for widx in web_index
            wsec = widx[1]
            if wsec != "YEAR"
                ws_type, ws_cat = string.(split(wsec, "_", limit = 2))
                if !(ws_cat in ["ALL", "CF"]); ci = findfirst(x -> x == web_cat_conc[ws_cat], cat_list) end

                if ws_type == "CFAV"
                    if ws_cat == "CF"; str *= "\t" * string(cfByNat[y][n][ri])
                    elseif ws_cat == "ALL"; str *= "\t" * string(cfByNat[y][n][ri] / pops[y][n][rg])
                    else str *= "\t" * string(cfByReg[y][n][ri, ci])
                    end
                elseif ws_type == "CFAC"
                    if ws_cat == "CF"; str *= "\t" * string(cfav[y][rid][end])
                    elseif ws_cat == "ALL"; str *= "\t" * string(cfac[y][rid][end])
                    else str *= "\t" * string(cfac[y][rid][ci])
                    end
                elseif ws_type == "CFAL"
                    if ws_cat == "CF"; str *= "\t" * string(ci_cf[y][n][rg][1])
                    elseif ws_cat == "ALL"; str *= "\t" * string(ci_cf[y][n][rg][1] / pops[y][n][rg])
                    else str *= "\t" * string(ci_cfpc[y][n][rg][ci][1])
                    end
                elseif ws_type == "CFAU"
                    if ws_cat == "CF"; str *= "\t" * string(ci_cf[y][n][rg][2])
                    elseif ws_cat == "ALL"; str *= "\t" * string(ci_cf[y][n][rg][2] / pops[y][n][rg])
                    else str *= "\t" * string(ci_cfpc[y][n][rg][ci][2])
                    end
                else str *= "\t" * string(widx[2])
                end
            end
        end
        strs[rid][y] = str

        # f = open(path * rid * ".txt", "a")
        # print(f, y)
        # for widx in web_index
        #     wsec = widx[1]
        #     if wsec != "YEAR"
        #         ws_type, ws_cat = string.(split(wsec, "_", limit = 2))
        #         if !(ws_cat in ["ALL", "CF"]); ci = findfirst(x -> x == web_cat_conc[ws_cat], cat_list) end
        #
        #         if ws_type == "CFAV"
        #             if ws_cat == "CF"; print(f, "\t", cfByNat[y][n][ri])
        #             elseif ws_cat == "ALL"; print(f, "\t", cfByNat[y][n][ri] / pops[y][n][rg])
        #             else print(f, "\t", cfByReg[y][n][ri, ci])
        #             end
        #         elseif ws_type == "CFAC"
        #             if ws_cat == "CF"; print(f, "\t", cfav[y][rid][end])
        #             elseif ws_cat == "ALL"; print(f, "\t", cfac[y][rid][end])
        #             else print(f, "\t", cfac[y][rid][ci])
        #             end
        #         elseif ws_type == "CFAL"
        #             if ws_cat == "CF"; print(f, "\t", ci_cf[y][n][rg][1])
        #             elseif ws_cat == "ALL"; print(f, "\t", ci_cf[y][n][rg][1] / pops[y][n][rg])
        #             else print(f, "\t", ci_cfpc[y][n][rg][ci][1])
        #             end
        #         elseif ws_type == "CFAU"
        #             if ws_cat == "CF"; print(f, "\t", ci_cf[y][n][rg][2])
        #             elseif ws_cat == "ALL"; print(f, "\t", ci_cf[y][n][rg][2] / pops[y][n][rg])
        #             else print(f, "\t", ci_cfpc[y][n][rg][ci][2])
        #             end
        #         else print(f, "\t", widx[2])
        #         end
        #     end
        # end
        # println(f)
        # close(f)
    end

    for rg in regs
        f = open(path * rg * ".txt", "w")
        print(f, web_index[1][1])
        for widx in web_index[2:end]; print(f, "\t", widx[1]) end
        println(f)
        for y in sort(collect(keys(strs[rg]))); println(f, strs[rg][y]) end
        close(f)
    end
end

end
