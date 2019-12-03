module EmissionCalculator

# Developed date: 3. Dec. 2019
# Last modified date: 3. Dec. 2019
# Subject: Calculate carbon emissions by final demands
# Description: Read Eora MRIO T, V, Y, Q tables and
#              calculate household consumption-based carbon emissions
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX

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
emissions = Array{Float64, 2}()

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

function calculateEmission()
    global emissions

    

    return emissions
end

end
