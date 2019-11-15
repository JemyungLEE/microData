module TradeMatrixBuilder

# Developed date: 7. Nov. 2019
# Last modified date: 11. Nov. 2019
# Subject: Nation by nation trade matrix builder
# Description: Build nation by nation trade matri based on UN Comtrade data.
#              Product categories are based on the Harmonized System (HS).
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

nat = Array{String, 1}()    # Nations list
hsSec = Array{String, 1}()  # HS classification's sectors; code, description
tdMat = Dict{String, Array{Float64, 2}}()   # Trade matrix: {Nation, Matrix}
shMat = Dict{String, Array{Float64, 2}}()   # Trade-share matrix: {Nation, Matrix}

function readHsCategories(inputFile, level=[]) # level = 2, 4, or 6

    global hsSec
    f = open(inputFile)

    for l in eachline(f)
        s = split(l, '\t')
        if length(level)==0 || length(s[2]) in level; push!(hsSec, s[2]) end
    end
    close(f)
end

function buildTradeMatrix(trades, nation=[], tradeFlow="")  # flow: "Import", "Export", or "Net"
    # build trade matrix of the reporter 'nation'
    global tdMat
    global nat

    # detect nations
    for t in trades
        if !(t.reporter in nat); push!(nat, t.reporter) end
        if !(t.partner in nat); push!(nat, t.partner) end
    end
    sort!(nat)

    # build matrix
    for n in nation
        mat = zeros(length(hsSec), length(nat))
        for t in trades
            if t.reporter == n && (tradeFlow == "Net" || t.flow == tradeFlow)
                r = findfirst(x -> x==t.hscd, hsSec)
                c = findfirst(x -> x==t.partner, nat)
                if tradeFlow == "Net" && t.flow == "Export"; mat[r,c] -= t.value
                else mat[r,c] += t.value
                end
            end
        end
        tdMat[n] = mat
    end

    return tdMat, nat, hsSec
end

function buildShareMatrix(nation=[])
    # calculate the share among the given nations in the 'nation' array
    # if 'nation' array is empty (default), then calculate the share among entire nations.

    global shMat

    # detect nation index
    if length(nation) == 0
        col = 1:length(nat)
    else
        col = []
        for n in nation; push!(col, findfirst(x -> x==n, nat)) end
    end

    for n in sort(collect(keys(tdMat)))
        # total trade values by each HS sector
        total = zeros(length(hsSec))
        m = tdMat[n]
        for r = 1:length(hsSec); for c in col; total[r] += m[r,c] end end

        # trade share proportions by each HS sector
        mat = zeros(length(hsSec), length(col))
        for r = 1:length(hsSec)
            for c = 1:length(col)
                if m[r,col[c]] != 0 && total[r] != 0; mat[r,c] = m[r,col[c]] / total[r] end
            end
        end

        shMat[n] = mat
    end

    if length(nation) == 0; return shMat, nat
    else return shMat, nation
    end
end

function printTradeMatrix(outputFile)

    f= open(outputFile, "w")
    print(f, "Nation\tSection")
    for n in nat; print(f, "\t", n) end
    println(f)
    for n in sort(collect(keys(tdMat)))
        m = tdMat[n]
        for r = 1:length(hsSec)
            print(f, n, "\t", hsSec[r])
            for c = 1:length(nat)
                print(f, "\t", m[r,c])
            end
            println(f)
        end
    end
    close(f)
end

function printShareMatrix(outputFile, nation=[])

    f= open(outputFile, "w")
    print(f, "Nation\tSection")
    if length(nation) == 0; for n in nat; print(f, "\t", n) end
    else for n in nation; print(f, "\t", n) end
    end
    println(f)
    for n in sort(collect(keys(shMat)))
        m = shMat[n]
        for r = 1:size(m,1)
            print(f, n, "\t", hsSec[r])
            for c = 1:size(m,2); print(f, "\t", m[r,c]) end
            println(f)
        end
    end
    close(f)
end

end
