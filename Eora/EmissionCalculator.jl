module EmissionCalculator

# Developed date: 3. Dec. 2019
# Last modified date: 3. Dec. 2019
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

abb = Dict{String, String}()    # Nation name A3 abbreviation
ti = Array{idx, 1}()     # index T
vi = Array{idx, 1}()     # index V
yi = Array{idx, 1}()     # index Y
qi = Array{ind, 1}()     # index Q

mTables = Dict{Int16, tables}()
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
    global ti, vi, yi, qi

    ti = ti[1:end-1]
    vi = vi[1:end-6]
    yi = yi[1:end-6]
    ql = deleteat!(collect(10:64), [12,13,43,44,45,46,47])
    qi = qi[ql]
end

function rearrangeTables(year)

    nt = length(ti)
    nv = length(vi)
    ny = length(yi)
    nq = length(qi)

    tb = mTables[year]
    ql = deleteat!(collect(10:64), [12,13,43,44,45,46,47])

    tb.t = tb.t[1:nt, 1:nt]
    tb.v = tb.v[1:nv, 1:nt]
    tb.y = tb.t[1:nt, 1:ny]
    tb.q = tb.t[ql, 1:nt]
end

function calculateEmission(year)

    global emissions
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
        for i = 1:nq; f[j] += tb.q[i,j] end
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

    return e
end

function extractHouseholdEmission(year, overwrite = false)

    global ti, yi, emissions

    tl = findall(x->x.entity=="Industry", ti)
    yl = findall(x->x.sector=="Household final consumption P.3h"&& x.nation!="ROW", yi)

    e = emissions[year][tl,yl]
    t = ti[tl]
    y = yi[yl]

    if overwrite
        ti = t
        yi = y
        emissions[year] = e
    end

    return e
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
