module HsDataReader

# Developed date: 21. Oct. 2019
# Last modified date: 8. Nov. 2019
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
    code::String        # commodity HS code, whole (= hslv+hscd)
    desc::String        # commodity description
    value::Int64        # trade value
    netw::Int64         # net weight (Kg)
    unit::String        # unit
    qnt::Int64          # trade quantity

    orig::String        # origin
    dest::String        # destination
    hslv::String        # HS code level: H0, H1, H2, H3
    hscd::String        # HS commodity code


    function trade(yr::Int16, fl::String, rp="", pt="", cd="", dc="", vl=0, nt=0, un="", qt=0, hl="", hc="")
        if fl == "Import" || fl == "Re-Import"
            or = pt
            dt = rp
        elseif fl == "Export" || fl == "Re-Export"
            or = rp
            dt = pt
        end

        new(yr, fl, rp, pt, cd, dc, vl, nt, un, qt, or, dt, hl, hc)
    end
end

titles = Array{String, 1}()
trades = Array{trade, 1}()
nations = Array{String, 1}()
df = DataFrame()

function readTradeData(inputFile, hsLv="")
    f = open(inputFile)

    global trades　
    global nations
    global titles = [string(t) for t in [strip(t) for t in split(readline(f), ',')]]

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
        hs = split(s[5], '-')

        if !(s[3] in nations); push!(nations, s[3]) end
        if !(s[4] in nations); push!(nations, s[4]) end
        sort!(nations)

        if length(hsLv)==0 || hsLv==hs[1]
            push!(trades,trade(yr, s[2], s[3], s[4], s[5], s[6], tv, nw, s[9], qt, hs[1], hs[2]))
        end
    end

    return trades, titles, nations

end

function analyzeStatus(rpt=[], prt=[], flow="")

    stat = zeros(Int32, 4, 3)   # row: H0, H1, H2, H3; column: level 2, 4, 6
    row = ["H0", "H1", "H2", "H3"]
    col = [2, 4, 6]

    for t in trades
        r = findfirst(x -> x==t.hslv, row)
        c = findfirst(x -> x==length(t.hscd), col)

        if length(rpt) == length(prt) == 0 || t.reporter in rpt || t.partner in prt
            if length(flow)==0; stat[r,c] += 1
            elseif flow == "Import" && t.dest in prt; stat[r,c] += 1
            elseif flow == "Export" && t.orig in rpt; stat[r,c] += 1
            end
        end
    end

    println("Comtrade status:")
    for c in col; print("\t", c) end
    println()
    for r = 1:length(row)
        print(row[r])
        for c = 1:length(col); print("\t", stat[r,c]) end
        println()
    end
end

function matchingTest(nation, outputFile)

    nls = String[]
    rptDic = Dict{String, Array{Float64,1}}()    # partner; [1]import, [2]export
    prtDic = Dict{String, Array{Float64,1}}()    # reporter;[1]import, [2]export

    for t in trades
        if t.reporter == nation
            if !(t.partner in nls); push!(nls, t.partner)
            end
            if !haskey(rptDic, t.partner); rptDic[t.partner] = zeros(2)
            end
            if t.flow == "Import"; rptDic[t.partner][1] += t.value
            elseif t.flow == "Export"; rptDic[t.partner][2] += t.value
            end
        elseif t.partner == nation
            if !(t.reporter in nls); push!(nls, t.reporter)
            end
            if !haskey(prtDic, t.reporter); prtDic[t.reporter] = zeros(2)
            end
            if t.flow == "Import"; prtDic[t.reporter][1] += t.value
            elseif t.flow == "Export"; prtDic[t.reporter][2] += t.value
            end
        end
    end

    if length(outputFile) > 0
        f = open(outputFile, "w")

        println(f, "Reporter: $nation\t\t\tPartner: $nation")
        println(f, "Partner\tImport\tExport\tReporter\tImport\tExport")
        for n in nls
            print(f, n)
            if haskey(rptDic, n); print(f, "\t", rptDic[n][1], "\t", rptDic[n][2])
            else print(f, "\t\t")
            end
            print(f, "\t", n)
            if haskey(prtDic, n); println(f, "\t", prtDic[n][1], "\t", prtDic[n][2])
            else println(f, "\t\t")
            end
        end

        close(f)
    end

    return rptDic, prtDic, nls
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
        push!(qt, t.qnt)
        push!(og, t.orig)
        push!(dt, t.dest)
        push!(hl, t.hslv)
        push!(hc, t.hscd)
    end

    global　df = DataFrame(year=yr,flow=fl,reporter=rp,partner=pt,code=cd,desc=dc,value=vl,netw=nw,unit=un,qnt=qt,orig=og,dest=dt,hslv=hl,hscd=hc)
    names!(df, Symbol.(tl))

    return df
end

function printDataFrames(outputFile)
    f = open(outputFile, "w")

    println(f, df)

    close(f)
end

end
