module EmissionEstimator

# Developed date: 26. Apr. 2021
# Last modified date: 17. May. 2021
# Subject: Calculate household carbon emissions
# Description: Calculate direct and indirect carbon emissions by analyzing
#              Customer Expenditure Survey (CES) or Household Budget Survey (HBS) micro-data.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX
using LinearAlgebra
using SparseArrays

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

natA3 = Dict{String, String}()  # Nation name's A3 abbreviation, {Nation, A3}
natList = Array{String, 1}()    # Nation A3 list

hh_list = Dict{Int, Dict{String, Array{String, 1}}}()   # Household ID list: {year, {nation A3, {hhid}}}
sc_list = Dict{Int, Dict{String, Array{String, 1}}}()   # commodity code list: {year, {nation A3, {code}}}
hhExp = Dict{Int, Dict{String, Array{Float64, 2}}}()    # Household expenditure matrix: {year, {nation, {commodity, hhid}}}
hhCmm = Dict{Int, Dict{String, Array{Float64, 2}}}()    # Household physical consumption matrix: {year, {nation, {commodity, hhid}}}

# indirect carbon emission variables
ti = Array{idx, 1}()                # index T
vi = Array{idx, 1}()                # index V
yi = Array{idx, 1}()                # index Y
qi = Array{ind, 1}()                # index Q
lti = []                            # inversed Leontief matrix
eoraExp = Dict{Int, Dict{String, Array{Float64, 2}}}()  # transformed households expenditure: {year, {nation, {Eora sectors, households}}}
mTables = Dict{Int, tables}()       # {Year, tables}
concMat = Dict{Int, Array{Float64, 2}}()  # Concordance matrix: {year, {Eora sectors, CES/HBS sectors}}
indirectCE = Dict{Int, Dict{String, Array{Float64, 2}}}()   # indirect carbon emission: {year, {nation, {CES/HBS sector, households}}}

# direct carbon emission variables
concMatDe = Dict{Int, Dict{String, Array{Float64, 2}}}()    # Concordance matrix for direct emission: {year, {nation, {DE sectors, {CES/HBS sectors}}}
de_sc_list = Dict{Int, Dict{String, Array{String, 1}}}()    # direct emission sectors: {year, {nation, {DE sector}}}
deIntens = Dict{Int, Dict{String, Array{Float64, 1}}}()     # converting rate from comodity's value (ex.USD) to tCO2: {year, {nation, {DE category, DE factor (ex.tCO2/USD)}}}
deUnits = Dict{Int, Dict{String, Array{String, 1}}}()       # units of converting rates: {year, {nation, {DE category, DE unit}}}
directCE = Dict{Int, Dict{String, Array{Float64, 2}}}()     # direct carbon emission: {year, {nation, {}}}

function readIndex(indexFilePath)

    global natA3, ti, vi, yi, qi

    f = open(indexFilePath*"a3.csv"); readline(f)
    for l in eachline(f); l = string.(split(replace(l,"\""=>""), ',')); natA3[l[1]] = l[2] end
    close(f)
    f = open(indexFilePath*"index_t.csv"); readline(f)
    for l in eachline(f); l=string.(split(replace(l,"\""=>""), ',')); push!(ti, idx(l[3],l[4],l[5])) end
    close(f)
    f = open(indexFilePath*"index_v.csv"); readline(f)
    for l in eachline(f); l=string.(split(replace(l,"\""=>""), ',')); push!(vi, idx(l[3],l[4],l[5])) end
    close(f)
    f = open(indexFilePath*"index_y.csv"); readline(f)
    for l in eachline(f); l=string.(split(replace(l,"\""=>""), ',')); push!(yi, idx(l[3],l[4],l[5])) end
    close(f)
    f = open(indexFilePath*"index_q.csv"); readline(f)
    for l in eachline(f); l=string.(split(replace(l,"\""=>""), ',')); push!(qi, ind(l[2],l[3],l[4])) end
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

function rearrangeIndex(; qmode = "")

    if qmode == "I_CHG_CO2"; ql = deleteat!(collect(10:64), [12,13,43,44,45,46,47]) # I-GHG-CO2 emissions (Gg)
    elseif qmode == "PRIMAP"; ql = [2570]                                           # PRIMAP|KYOTOGHGAR4|TOTAL|GgCO2eq
    else println("Define Q_table mode.")
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
    if qmode == "I_CHG_CO2"; ql = deleteat!(collect(10:64), [12,13,43,44,45,46,47]) # I-GHG-CO2 emissions (Gg)
    elseif qmode == "PRIMAP"; ql = [2570]                                           # PRIMAP|KYOTOGHGAR4|TOTAL|GgCO2eq
    else println("Define Q_table mode.")
    end

    tb = mTables[year]

    nt = length(ti)
    tb.t = tb.t[1:nt, 1:nt]
    tb.v = tb.v[1:length(vi), 1:nt]
    tb.y = tb.y[1:nt, 1:length(yi)]
    tb.q = tb.q[ql, 1:nt]
end

function getDomesticData(year, nation, hhid_list, sector_list, expMat, qntMap)

    global hh_list, sc_list, hhExp, hhCmm

    if !haskey(hh_list, year); hh_list[year] = Dict{String, Array{String, 1}}() end
    if !haskey(sc_list, year); sc_list[year] = Dict{String, Array{String, 1}}() end
    if !haskey(hhExp, year); hhExp[year] = Dict{String, Array{Float64, 2}}() end
    if !haskey(hhCmm, year); hhCmm[year] = Dict{String, Array{Float64, 2}}() end

    hl = hh_list[year][nation] = hhid_list
    sl = sc_list[year][nation] = sector_list
    if size(expMat,1) == length(sl) && size(expMat,2) == length(hl); global hhExp[year][nation] = expMat
    elseif size(expMat,2) == length(sl) && size(expMat,1) == length(hl); global hhExp[year][nation] = transpose(expMat)
    else println("Matrices sizes don't match: expMat,", size(expMat), "\thhid,", size(hl), "\tsec,", size(sl))
    end
    if size(qntMap,1) == length(sl) && size(qntMap,2) == length(hl); global hhCmm[year][nation] = qntMap
    elseif size(qntMap,2) == length(sl) && size(qntMap,1) == length(hl); global hhCmm[year][nation] = transpose(qntMap)
    else println("Matrices sizes don't match: expMat,", size(qntMap), "\thhid,", size(hl), "\tsec,", size(sl))
    end
end

function readEmissionIntensity(year, nation, sectorFile, intensityFile; quantity = false, emit_unit = "tCO2", curr_unit = "USD")

    global de_sc_list, deIntens
    if !haskey(de_sc_list, year); de_sc_list[year] = Dict{String, Array{String, 1}}() end
    if !haskey(deIntens, year); deIntens[year] = Dict{String, Array{Float64, 1}}() end
    if !haskey(deUnits, year); deUnits[year] = Dict{String, Array{String, 1}}() end
    if !haskey(de_sc_list[year], nation); de_sc_list[year][nation] = Array{String, 1}() end
    if !haskey(deIntens[year], nation); deIntens[year][nation] = Array{Float64, 1}() end
    if !haskey(deUnits[year], nation); deUnits[year][nation] = Array{String, 1}() end
    dsl = de_sc_list[year][nation]
    qnt_units = ["liter", "kg", "m^3"]

    f_sep = getValueSeparator(sectorFile)
    f = open(sectorFile)
    ci = findfirst(x->x=="DE_code", string.(strip.(split(readline(f), f_sep))))
    for l in eachline(f)
        c = string.(strip.(split(l, f_sep)))[ci]
        if !(c in dsl); push!(dsl, c) end
    end
    close(f)

    nds = length(dsl)
    dit, dun = zeros(Float64, nds), Array{String, 1}(undef, nds)
    f_sep = getValueSeparator(intensityFile)
    essential = ["Year", "Nation", "DE_code", "DE_sector", "Factor", "Unit"]
    f = open(intensityFile)
    title = string.(strip.(split(readline(f), f_sep)))
    i = [findfirst(x->x==et, title) for et in essential]
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        if s[i[1]] == string(year) && s[i[2]] == nation
            unit = string.(strip.(split(s[i[6]], '/')))
            if emit_unit == unit[1] && ((!quantity && curr_unit == unit[2]) || (quantity && unit[2] in qnt_units))
                di = findfirst(x->x==s[i[3]], dsl)
                dit[di], dun[di] = parse(Float64, s[i[5]]), s[i[6]]
            end
        end
    end
    deIntens[year][nation] = dit
    deUnits[year][nation] = dun
    close(f)
end

function calculateDirectEmission(year, nation, cm_de; quantity = false, sparseMat = false, enhance = false, full = false)

    global sc_list, de_sc_list, hh_list, hhExp, directCE, concMatDe, deIntens, deUnits
    hl, sl, dsl = hh_list[year][nation], sc_list[year][nation], de_sc_list[year][nation]
    if !haskey(concMatDe, year); concMatDe[year] = Dict{String, Array{Float64, 2}}() end
    cmn = concMatDe[year][nation] = cm_de
    dit = transpose(deIntens[year][nation])
    if quantity; he = hhCmm[year][nation]; else he = hhExp[year][nation] end
    nh, ns, nds = length(hl), length(sl), length(dsl)
    if !haskey(directCE, year); directCE[year] = Dict{String, Array{Float64, 2}}() end

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
    directCE[year][nation] = de
end

function buildWeightedConcMat(year, eoraYear, natA3, conMat; output="")
    # concordance matrix (Eora, Nation)

    global concMat, mTables
    global natList, sc_list, ti, yi
    sl = sc_list[year][natA3]
    tb = mTables[eoraYear]
    ns, nn = length(sl), length(natList)

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

    # assemble concordance matrices
    cMat = zeros(Float64, 0, ns)
    for i = 1:nn
        if cnt[i]>0; cMat = vcat(cMat, zeros(Float64, cnt[i], ns)) end
        cMat = vcat(cMat, conMat[natList[i]])
    end

    # reflect Eora final demand accounts' ratios
    for j = 1:ns
        cMat[:, j] .*= ye
        tsum = sum(cMat[:,j])
        cMat[:, j] /= tsum
    end
    concMat[year] = cMat

    # print concordance matrix
    if length(output)>0
        mkpath(rsplit(output, '/', limit = 2)[1])
        f = open(output, "w")
        print(f,"Nation,Entity,Sector");for i=1:ns; print(f,",",sl[i]) end; println(f)
        for i=1:size(concMat[year],1)
            print(f,ti[i].nation,",",ti[i].entity,",",ti[i].sector)
            for j=1:size(concMat[year],2); print(f,",",concMat[year][i,j]) end; println(f)
        end
        close(f)
    end
end

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

function calculateIndirectEmission(cesYear, eoraYear, nation; sparseMat = false, enhance = false, full = false, elapChk = 0)

    global indirectCE, mTables, concMat, lti
    global hh_list, sc_list, hhExp
    global ti, vi, yi, qi

    tb = mTables[eoraYear]
    hl, sl, em = hh_list[cesYear][nation], sc_list[cesYear][nation], hhExp[cesYear][nation]
    nt, nh, ns = length(ti), length(hl), length(sl)
    if !haskey(indirectCE, cesYear); indirectCE[cesYear] = Dict{String, Array{Float64, 2}}() end

    # calculate emission, by CES/HBS micro-data sectors, by Eora T matrix sectors
    e = zeros(Float64, ns, nh)

    st = time()     # check start time
    if enhance || full; lti_conc = lti * concMat[cesYear] end
    for i = 1:ns
        if enhance; e[i,:] = sum(lti_conc[:,i] * transpose(em[i,:]), dims=1)
        else
            hce = zeros(Float64, ns, nh)
            hce[i,:] = em[i,:]

            if sparseMat
                hceS = dropzeros(sparse(hce))
                hce = []
                concMatS = zeros(Float64, nt, ns)
                concMatS[:,i] = concMat[cesYear][:,i]
                concMatS = dropzeros(sparse(concMatS))
                conc_hce = Array(concMatS * hceS)
                concMatS = hceS = []
                ebe = lti * conc_hce
            elseif full; ebe = lti_conc * hce
            else ebe = lti * concMat[cesYear] * hce       # household emission by Eora sectors
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
                println(i,"/",ns," iterations, ",eHr,":",eMin,":",eSec," elapsed, ",rHr,":",rMin,":",rSec," remained")
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
