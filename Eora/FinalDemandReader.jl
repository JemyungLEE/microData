module FinalDemandReader

# Developed date: 18. Nov. 2019
# Last modified date: 20. Nov. 2019
# Subject: Eora final demand sector reader
# Description: Read household final demand data from Eora MRIO.
#              Classify by nation and commiodity category
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX

fdNat = Array{String, 1}()          # destination: final demand (FD) reading nation
abb = Dict{String, String}()        # abbreviation: {Eora's nation name, A3 abbriviation}
nations = Dict{String, String}()    # for matching: {Eora's nation name, Comtrade's nation name}
sections = Dict{String, Array{Tuple{String, Int}, 1}}()     # {a3, {(section, index) list}}
columns = Dict{String, Array{Tuple{String, Int}, 1}}()      # {a3, {(final demand section, index) list}}
fdTable = Array{Array{Float64, 1}, 1}()         # Eora Final Demand table
fdTables = Dict{String, Array{Float64, 2}}()    # {a3, {FD tables, row: sections, col: destinations}}

function readIndexData(inputFile)

    global abb, nations, sections

    xf = XLSX.readxlsx(inputFile)

    sh = xf["A3"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r) > 1; abb[r[1]] = strip(r[2]) end end

    sh = xf["Nation"]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r) > 1 && length(r[1])>0
            if r[2] !== missing; nations[r[1]] = strip(r[2])
            else nations[r[1]] = ""
            end
        end
    end

    sh = xf["index_t"]  # read row index of Eora Y table
    nat = ""        # a3
    entity = ""     # "Industries" or "Commodities"
    sec = []        # {(section, index) list}
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r) == 1 || r[3] == "ROW"; continue    # ignore the first and the last lines
        elseif r[4] != "Industries" && r[4] != "Commodities"
            println("Entity error: ",r[1],"\t",r[2],"\t",r[3],"\t",r[4],"\t",r[5])
        end

        if nat != r[3]
            if length(nat) > 0; sections[nat] = sec end
            nat = r[3]
            entity = r[4]
            sec = Array{Tuple{String, Int}, 1}()
        elseif nat == r[3] && entity != r[4]
            if r[4] == "Commodities"
                entity = r[4]
                sec = Array{Tuple{String, Int}, 1}()
            end
        end
        push!(sec, (strip(r[5]), r[1]))
    end

    sh = xf["index_y"]      # read column index of Eora Y table
    nat = ""    # a3
    sec = []    # {(final demand section, index) list}
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r) == 1 || r[3] == "ROW"; continue   # ignore the first and the last lines
        elseif r[4] != "Final Demand"; println("Entity error: ",r[1],"\t",r[2],"\t",r[3],"\t",r[4],"\t",r[5])
        end

        if nat != r[3]
            if length(nat) > 0; columns[nat] = sec end
            nat = r[3]
            sec = Array{Tuple{String, Int}, 1}()
        end
        push!(sec, (strip(r[5]), r[1]))
    end

    return abb, nations, sections, columns
end

function readFinalDemand(inputFile, nat=[])

    global fdNat = nat
    global fdTable, fdTables

    # read Eora final demand Y table
    f = open(inputFile)
    for l in eachline(f)
        t = []
        for r in split(l, ','); push!(t, parse(Float64, r)) end
        push!(fdTable, t)
    end
    close(f)

    # extract target nations' household final demand
    idx = []
    for n in nat; push!(idx, columns[abb[n]][1][2]) end
    for n in sort(collect(keys(sections)))
        sec = sections[n]
        fd = zeros(Float64, length(sec), length(nat))
        for r = 1:length(sec); for c = 1:length(nat); fd[r,c] = fdTable[sec[r][2]][idx[c]] end end
        fdTables[n] = fd
    end

    return fdTables
end

function test(outputFile)

    f = open(outputFile, "w")

    println(f, "[ABB]")
    for n in sort(collect(keys(abb))); println(f, n, "\t", abb[n]) end
    println(f)

    println(f, "[Nations]")
    for n in sort(collect(keys(nations))); println(f, n, "\t", nations[n]) end
    println(f)

    println(f, "[Sections]")
    for n in sort(collect(keys(sections)))
        for s in sections[n]; println(f, n, "\t", s[1], "\t", s[2]) end
    end
    println(f)

    println(f, "[Columns]")
    for n in sort(collect(keys(columns)))
        for c in columns[n]; println(f, n, "\t", c[1], "\t", c[2]) end
    end
    println(f)

    println(f, "[Final demand]")
    for n in fdNat; print(f, "\t", n) end
    println(f)
    for n in sort(collect(keys(fdTables)))
        for r = 1:length(sections[n])
            print(f, n, "\t", sections[n][r][1])
            for c = 1:length(fdNat); print(f, "\t", fdTables[n][r][c]) end
            println(f)
        end
    end
    println(f)

    close(f)
end

end
