module EmissionEstimator

# Developed date: 29. Jul. 2020
# Last modified date: 14. Aug. 2020
# Subject: Calculate EU households carbon emissions
# Description: Calculate emissions by analyzing Eurostat Household Budget Survey (HBS) micro-data.
#              Transform HH consumptions matrix to nation by nation matrix of Eora form.
#              Intermidiate results can be utilized for the global carbon footprint mapping.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX
using LinearAlgebra
using SparseArrays

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

abb = Dict{String, String}()    # Nation name's A3 abbreviation, {Nation, A3}
natList = Array{String, 1}()    # Nation A3 list
ti = Array{idx, 1}()     # index T
vi = Array{idx, 1}()     # index V
yi = Array{idx, 1}()     # index Y
qi = Array{ind, 1}()     # index Q

sec = Dict{Int16, Array{String, 1}}()           # Household expenditure sectors: {year, {sector}}
hhid = Dict{Int16, Array{String, 1}}()          # Household ID: {year, {hhid}}
hhExp = Dict{Int16, Array{Float64, 2}}()        # Households enpenditure: {year, {Nation sectors, households}}
concMat = Dict{Int16, Array{Float64, 2}}()      # Concordance matrix {Eora sectors, Nation sectors}

eoraExp = Array{Float64, 2}         # Transformed households expenditure, {Eora sectors, households}
mTables = Dict{Int16, tables}()     # {Year, tables}
emissions = Dict{Int16, Array{Float64, 2}}()    # carbon footprint
directEms = Dict{Int16, Array{Float64, 2}}()    # direct emission

lti = []                            # Inversed Leontief matrix

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

function getDomesticData(year, expMat, householdID, sector)
    global hhid[year] = householdID
    global sec[year] = sector

    if size(expMat,1) == length(sec[year]) && size(expMat,2) == length(hhid[year]); global hhExp[year] = expMat
    elseif size(expMat,2) == length(sec[year]) && size(expMat,1) == length(hhid[year]); global hhExp[year] = transpose(expMat)
    else println("Matrices sizes don't match: expMat,", size(expMat), "\thhid,", size(hhid[year]), "\tsec,", size(sec[year]))
    end
end

function buildWeightedConcMat(year, nat, conMat; output="") # feasical year, nation A3, concordance matrix (Eora, Nation)

    global concMat, mTables
    global natList, sec, ti, yi
    tb = mTables[year]
    ns = length(sec[year])

    # get final demand of nation 'nat'
    ye = tb.y[:,findfirst(x->x.nation==nat && x.sector=="Household final consumption P.3h", yi)]

    # check whether a nation's account have 'Commodities' entities
    chk = Dict{String, Bool}()      # {A3, (true: has 'Commodities', false: only 'Industries')}
    for t in ti
        if !haskey(chk, t.nation); chk[t.nation] = false end
        if !chk[t.nation] && t.entity == "Commodities"; chk[t.nation] = true end
    end

    # count number of 'Industries' entities by nation: 0 means having only 'Industries' or only 'Commodities'.
    cnt = zeros(Int, length(natList))
    for t in ti
        if chk[t.nation] && t.entity == "Industries"
            cnt[findfirst(x->x==t.nation, natList)] += 1
        end
    end

    # assemble concordance matrices
    cMat = zeros(Float64, 0, ns)
    for i = 1:length(natList)
        if cnt[i]>0; cMat = vcat(cMat, zeros(Float64, cnt[i], ns)) end
        cMat = vcat(cMat, conMat[natList[i]])
    end

    # reflect ratios of Eora final demand accounts
    for j = 1:length(sec[year])
        cMat[:, j] .*= ye
        tsum = sum(cMat[:,j])
        cMat[:, j] /= tsum
    end
    concMat[year] = cMat

    # print concordance matrix
    if length(output)>0
        f = open(output, "w")
        print(f,"Nation,Entity,Sector");for i=1:ns; print(f,",",sec[year][i]) end; println(f)
        for i=1:size(concMat[year],1)
            print(f,ti[i].nation,",",ti[i].entity,",",ti[i].sector)
            for j=1:size(concMat[year],2); print(f,",",concMat[year][i,j]) end; println(f)
        end
        close(f)
    end

    return concMat[year], ti, sec[year]
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

function calculateDirectEmission(year)

    global sec, hhid, hhExp, directEms
end

function calculateCarbonFootprint(year, sparseMat = false, elapChk = 0; reuseLti = false)

    global emissions, mTables, concMat, lti
    global sec, hhid, hhExp
    global ti, vi, yi, qi

    tb = mTables[year]
    nt = length(ti)
    ns = length(sec[year])
    nh = length(hhid[year])

    # if !reuseLti || length(lti) == 0
    #     x = sum(tb.t, dims = 1) +  sum(tb.v, dims = 1)  # calculate X
    #     f = sum(tb.q, dims=1) ./ x                      # calculate EA
    #     lt = Matrix{Float64}(I, nt, nt)                 # calculate Leontief matrix
    #     for i = 1:nt; for j = 1:nt; lt[i,j] -= tb.t[i,j] / x[j] end end
    #     lti = inv(lt)
    #     for i = 1:nt; lti[i,:] *= f[i] end
    # end

    # calculate emission, by Eurostat micro-data sectors, by Eora T matrix sectors
    e = zeros(Float64, ns, nh)

    if sparseMat
        concMatS = SparseArrays.sortSparseMatrixCSC!(sparse(concMat[year]), sortindices=:doubletranspose)
        ltiS = SparseArrays.sortSparseMatrixCSC!(sparse(lti), sortindices=:doubletranspose)
        concMat[year] = []
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
        else ebe = lti * concMat[year] * hce       # household emission by Eora sectors
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

function printEmissions(year, outputFile)

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

end
