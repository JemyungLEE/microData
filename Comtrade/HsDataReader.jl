module HsDataReader

# Developed date: 21. Oct. 2019
# Last modified date: 24. Oct. 2019
# Subject: Harmonized System (HS) UN comtrade data reader
# Description: read trade data of HS classification
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

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

    function trade(yr::Int16, fl::String, rp="", pt="", cd="", dc="", vl=0, nt=0, un="", qt=0)
        if fl == "Import"
            or = pt
            dt = rp
        elseif fl == "Export"
            or = rp
            dt = pt
        end
        new(yr, fl, rp, pt, cd, dc, vl, nt, un, qt, or, dt)
    end
end

trades = Array{trade, 1}

function readTradeData(inputFile)
    f = open(inputFile)

    global trades = []
    s = ""

    readline(f)
    for l in eachline(f)
        s = split(l, ',')
        s = [strip(c, '"') for c in s]
        s = [string(c) for c in s]

        push!(trades,trade(parse(Int16,s[1]),s[2],s[3],s[4],s[5],s[6],parse(Int64,s[7]),parse(Int64,s[8]),s[9],parse(Int64,s[10])))
    end

    return trades
end

function readTradeDataCSV(inputFile)
    using CSV

    td = CSV.read(inputFile)

    println(CSV.describe(td))

    return trades
end

end
