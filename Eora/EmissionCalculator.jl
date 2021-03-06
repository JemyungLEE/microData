module EmissionCalculator

# Developed date: 3. Dec. 2019
# Last modified date: 18. Dec. 2019
# Subject: Calculate carbon emissions by final demands
# Description: Read Eora MRIO T, V, Y, Q tables and
#              calculate household consumption-based carbon emissions
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX
using LinearAlgebra

mutable struct tables
    year::Int16
    t::Array{Float64, 2}    # T matrix, MRIO tables
    v::Array{Float64, 2}    # V matrix, value added
    y::Array{Float64, 2}    # Y matrix, final demand
    q::Array{Float64, 2}    # Q matrix, indicators (satellite accounts)

    function tables(year, nt, nv, ny, nq)
        tt = zeros(Float64, nt, nt)
        vt = zeros(Float64, nv, nt)
        yt = zeros(Float64, nt, ny)
        qt = zeros(Float64, nq, nt)

        new(year, tt, vt, yt, qt)
    end
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
ti = Array{idx, 1}()     # index T
vi = Array{idx, 1}()     # index V
yi = Array{idx, 1}()     # index Y
qi = Array{ind, 1}()     # index Q

mTables = Dict{Int16, tables}()     # {Year, tables}
emissions = Dict{Int16, Array{Float64, 2}}()

function readIndexXlsx(inputFile)

    global abb, ti, vi, yi, qi
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
end

function readIOTables(year, tfile, vfile, yfile, qfile)

    global mTables

    nt = length(ti)
    tb = tables(year, nt, length(vi), length(yi), length(qi))

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
    for l in eachline(f); tb.q[i,:] = [parse(Float64, x) for x in split(l, ',')][1:nt]; i += 1 end
    close(f)

    mTables[year] = tb
end

function rearrangeIndex()
    global ti = ti[1:end-1]
    global vi = vi[1:end-6]
    global yi = yi[1:end-6]
    ql = deleteat!(collect(10:64), [12,13,43,44,45,46,47])
    global qi = qi[ql]
end

function rearrangeTables(year)
    global ti, vi, yi, qi, mTables

    tb = mTables[year]
    ql = deleteat!(collect(10:64), [12,13,43,44,45,46,47])

    nt = length(ti)
    tb.t = tb.t[1:nt, 1:nt]
    tb.v = tb.v[1:length(vi), 1:nt]
    tb.y = tb.t[1:nt, 1:length(yi)]
    tb.q = tb.t[ql, 1:nt]
end

function getMatchingIndex()
    # search commodity sectors' industry-sector-match indexes
    global ti

    nat = ""
    idx = zeros(Int, length(ti))
    for i = 1:length(ti)
        if ti[i].entity == "Industries"
            idx[i] = i
            nat = ti[i].nation
        elseif ti[i].entity == "Commodities"
            fi = findfirst(x->x.nation==nat && x.entity=="Industries" && x.sector==ti[i].sector, ti)
            if fi != nothing; idx[i] = fi
            else idx[i] = -1
            end
        end
    end

    #=
    f= open("/Users/leejimac/github/microData/Eora/data/test.txt", "w")
    for i = 1:length(ti); println(f, i,"\t",ti[i].nation,"\t",ti[i].entity,"\t",ti[i].sector,"\t",idx[i]) end
    close(f)
    =#

    return idx
end

function calculateEmission(year, match = []) # match[]: commodity sectors' industry sectors correponding indexes

    global emissions, mTables, ti, vi, yi, qi
    tb = mTables[year]

    nt = length(ti)
    nv = length(vi)
    ny = length(yi)
    nq = length(qi)

    # calculate X
    x = zeros(Float64, nt)
    for j = 1:nt
        for i = 1:nt; x[j] += tb.t[i,j] end
        for i = 1:nv; x[j] += tb.v[i,j] end
    end

    # calculate EA
    f = zeros(Float64, nt)
    for j = 1:nt
        for i = 1:nq
            if length(match)==0; f[j] += tb.q[i,j]
            else f[j] += tb.q[match[i],j]
            end
        end
        f[j] /= x[j]
    end

    # calculate Leontief matrix
    lt = Matrix{Float64}(I, nt, nt)
    for i = 1:nt; for j = 1:nt; lt[i,j] -= tb.t[i,j] / x[j] end end

    # calculate emissions
    lti = inv(lt)
    for i = 1:nt; for j = 1:nt; lti[i,j] *= f[i] end end
    e = lti * tb.y[1:nt, 1:ny]

    emissions[year] = e

    #=
    fl = open("/Users/leejimac/github/microData/Eora/data/test1.txt", "w")
    for i = 1:nt; println(fl, ti[i].nation,"\t",ti[i].entity,"\t",ti[i].sector,"\t",f[i]) end
    close(fl)

    nation = "IND"
    fl = open("/Users/leejimac/github/microData/Eora/data/test2.txt", "w")
    tmp = inv(lt)
    print(fl, "\t\t")
    for i = 1:nt
        if nation == "ALL" || ti[i].nation == nation
            print(fl, "\t", ti[i].nation,"_",ti[i].entity,"_",ti[i].sector)
        end
    end
    println(fl)
    for i = 1:nt
        print(fl, ti[i].nation,"\t",ti[i].entity,"\t",ti[i].sector)
        for j = 1:nt
            if nation == "ALL" || ti[j].nation == nation
                print(fl,"\t",tmp[i,j])
            end
        end
        println(fl)
    end
    close(fl)

    fl = open("/Users/leejimac/github/microData/Eora/data/test3.txt", "w")
    tmp *= tb.y[1:nt, 1:ny]
    print(fl, "nation\tentity\tsector")
    for i = 1:ny; if nation == "ALL" || ti[i].nation == nation
        print(fl, "\t", yi[i].nation,"_",yi[i].sector)
    end end
    println(fl)
    for i = 1:nt
        print(fl, ti[i].nation,"\t",ti[i].entity,"\t",ti[i].sector)
        for j = 1:ny; if nation == "ALL" || ti[j].nation == nation; print(fl,"\t",tmp[i,j]) end end
        println(fl)
    end
    close(fl)
    =#

    return e
end

function extractHouseholdEmission(year, overwrite = false)

    global ti, yi, emissions

    ent = Dict{String, Bool}()      # {nation a3, whether have commodity entities}
    for i = 1:length(ti)
        if !haskey(ent, ti[i].nation); ent[ti[i].nation] = false
        elseif !ent[ti[i].nation] && ti[i].entity == "Commodities"; ent[ti[i].nation] = true
        end
    end

    tl = findall(x->(x.entity=="Industries"&&!ent[x.nation])||(x.entity=="Commodities"&&ent[x.nation]), ti)
    yl = findall(x->x.sector=="Household final consumption P.3h", yi)

    e = emissions[year][tl,yl]
    t = ti[tl]
    y = yi[yl]

    if overwrite
        ti = t
        yi = y
        emissions[year] = e
    end

    return e, t, y
end

function getNationEmission(year, nation, overwrite = false)

    global ti, yi, emissions

    yl = findall(x->x.nation==nation, yi)
    y = yi[yl]

    e = emissions[year][:,yl]

    if overwrite
        yi = y
        emissions[year] = e
    end

    return e, ti, y
end

function getEmissionDataset(year, nation)
    # getNationEmission function should be operated in advance

    global ti, vi, yi, emissions
    nat = String[]                              # nations' a3
    sec = Dict{String, Array{String, 1}}()      # {a3, {section list}}
    ceMat = Dict{String, Array{Float64, 1}}()   # {a3, {FD table, row: sections}}

    sl = String[]
    j = findfirst(x->x.nation==nation, yi)
    for i = 1:length(ti)
        n = ti[i].nation
        if !(n in nat)
            push!(nat, n)
            sec[n] = String[]
            ceMat[n] = Float64[]
        end

        push!(sec[n], ti[i].sector)
        push!(ceMat[n], emissions[year][i,j])
    end

    return ceMat, sec, nat, abb
end

function printEmissions(year, outputFile)

    f = open(outputFile, "w")
    e = emissions[year]

    nt = length(ti)
    ny = length(yi)

    for y in yi; print(f, "\t", y.nation) end
    println(f)
    for i = 1:nt
        print(f, ti[i].nation, "_", ti[i].sector)
        for j = 1:ny; print(f, "\t", e[i,j]) end
        println(f)
    end

    close(f)
end

function printTables(year, outputFile)

    f= open(outputFile, "w")
    tb = mTables[year]

    nt = length(ti)
    nv = length(vi)
    ny = length(yi)
    nq = length(qi)

    print(f, "[T table]")
    for i = 1:nt; print(f, "\t") end
    println(f, "\t\t[Y table]")
    for i in ti; print(f, "\t", i.nation, "_", i.sector) end
    print(f, "\t")
    for i in yi; print(f, "\t", i.nation, "_", i.sector) end
    println(f)
    for i = 1:nt
        print(f, ti[i].nation, "_", ti[i].sector)
        for j = 1:nt; print(f, "\t", tb.t[i,j]) end
        print(f, "\t")
        for j = 1:ny; print(f, "\t", tb.y[i,j]) end
        println(f)
    end
    println(f)

    println(f,"[V table]")
    for i = 1:nv
        print(f, vi[i].nation, "_", vi[i].sector)
        for j = 1:nt; print(f, "\t", tb.v[i,j]) end
        println(f)
    end

    println(f,"[Q table]")
    for i = 1:nq
        print(f, qi[i].code, "_", qi[i].item)
        for j = 1:nt; print(f, "\t", tb.q[i,j]) end
        println(f)
    end

    close(f)
end

end
