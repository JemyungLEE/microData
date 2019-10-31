module HsDataReader

# Developed date: 21. Oct. 2019
# Last modified date: 30. Oct. 2019
# Subject: Harmonized System (HS) UN comtrade data reader
# Description: read trade data of HS classification, store the data, and return it as DataFrames
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using CSV
using DataFrames

mutable struct trade
    year::Int16
    flow::String
    reporter::String
    partner::String
    code::String        # commodity HS code
    desc::String        # commodity description
    value::Int64        # trade value
    netw::Int64         # net weight (Kg)
    unit::String        # unit
    quant::Int64           # trade quantity

    orig::String        # origin
    dest::String        # destination
    hslv::String        # HS code level: H0, H1, H2, H3
    hscd::String        # HS commodity code


    function trade(yr::Int16, fl::String, rp="", pt="", cd="", dc="", vl=0, nt=0, un="", qt=0)
        if fl == "Import" || fl == "Re-Import"
            or = pt
            dt = rp
        elseif fl == "Export" || fl == "Re-Export"
            or = rp
            dt = pt
        end
        hs = split(cd, '-')
        new(yr, fl, rp, pt, cd, dc, vl, nt, un, qt, or, dt, hs[1], hs[2])
    end
end

titles = Array{String, 1}
trades = Array{trade, 1}

function readTradeData(inputFile)
    f = open(inputFile)

    global titles = []
    global trades = []
    s = ""

    titles = [string(t) for t in [strip(t) for t in split(readline(f), ',')]]

    for l in eachline(f)
        s = split(l, r",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)")
        s = [strip(c, '"') for c in s]
        s = [string(c) for c in s]      # transfer a substring to a string array

        if length(s[1]) > 0; yr = parse(Int16,s[1])     # period
        else yr = -1 end
        if length(s[7]) > 0; tv = parse(Int64,s[7])     # trade value
        else tv = -1 end
        if length(s[8]) > 0; nw = parse(Int64,s[8])     # net weight
        else nw = -1 end
        if length(s[10]) > 0; qt = parse(Int64,s[10])   # quantity
        else qt = -1 end

        push!(trades,trade(yr, s[2], s[3], s[4], s[5], s[6], tv, nw, s[9], qt))
    end

    return trades

end

function readTradeDataCSV(inputFile)
    td = CSV.read(inputFile)

    println(CSV.describe(td))

    return trades
end

function exportDataFrames()

    tl = [replace(t, r"[()]"=> "") for t in [replace(t, " "=> "_") for t in titles]]
    append!(tl, ["Origin", "Destination", "HS_class", "HS_code"])

    yr = Int16[]
    fl = String[]
    rp = String[]
    pt = String[]
    cd = String[]
    dc = String[]
    vl = Int64[]
    nw = Int64[]
    un = String[]
    qt = Int64[]
    og = String[]
    dt = String[]
    hl = String[]
    hc = String[]

    for t in trades
        push!(yr, t.year)
        push!(fl, t.flow)
        push!(rp, t.reporter)
        push!(pt, t.partner)
        push!(cd, t.code)
        push!(dc, t.desc)
        push!(vl, t.value)
        push!(nw, t.netw)
        push!(un, t.unit)
        push!(qt, t.quant)
        push!(og, t.orig)
        push!(dt, t.dest)
        push!(hl, t.hslv)
        push!(hc, t.hscd)
    end

    df = DataFrame(year=yr,flow=fl,reporter=rp,partner=pt,code=cd,desc=dc,value=vl,netw=nw,unit=un,quant=qt,orig=og,dest=dt,hslv=hl,hscd=hc)
    names!(df, Symbol.(tl))

    return df

end

end
