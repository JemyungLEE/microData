module ConcMatBuilder

# Developed date: 14. Apr. 2021
# Last modified date: 14. Apr. 2021
# Subject: Build concordance matric between MRIO and HBS/CES micro-data
# Description: read sector matching information from a XLSX/TXT/CSV file and
#              build concordance matrix bewteen converting nation and Eora accounts
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX

mutable struct sector       # category data
    source::String          # 3-digit abbreviation of nation, cf) "_SD": standard 26 categcories, "_EU": EU 61 categcories
    code::String            # classification codes of each Eora classification
    categ::String           # category of Eora industry or commodity classification
    linked::Array{String,1}     # linked converting nation's sector codes
    weight::Array{Float64,1}    # linked converting nation's sector weight

    sector(src, cod, cat) = new(src, cod, cat, [], [])
end

mutable struct nation       # nation data
    name::String
    abb::String             # abbreviation of country name
    ns::Int                 # number of sectors
    hasComEn::Bool          # wether the nation have 'Commodities'-entity data
    matchCode::String       # matching nation code of this nation. cf) "_SD": standard 26 categcories, "_EU": EU 61 categcories
    sectors::Dict{String, sector}   # Eora sectors: {code, sector}

    nation(n::String, a::String, ns::Int, has::Bool, mc::String) = new(n, a, ns, has, mc, Dict{String, sector}())
end

mutable struct conTab       # concordance tables
    conMat::Array{Int,2}    # concordance matrix
    sumEora::Array{Int,1}   # sums of eora sectors
    sumNat::Array{Int,1}    # sums of converting nation's sectors

    function conTab(eorSecNum, natSecNum)
        new(zeros(Int, eorSecNum, natSecNum), zeros(Int, eorSecNum), zeros(Int, natSecNum))
    end
end

mutable struct conTabNorm   # normalized concordance tables
    conMat::Array{Float64,2}  # concordance matrix
    sumEora::Array{Float64,1} # sums of eora sectors
    sumNat::Array{Float64,1}  # sums of converting nation's sectors

    function conTabNorm(eorSecNum, natSecNum)
        new(zeros(Int, eorSecNum, natSecNum), zeros(Int, eorSecNum), zeros(Int, natSecNum))
    end
end

totals = 0  # total sectors
names = Dict{String, String}()      # Full names, abbreviation
nations = Dict{String, nation}()    # abbreviation, nation
natCodes = Array{String, 1}()       # converting nation's code list
eorCodes = Dict{String, Array{String, 1}}() # Eora's code list: [nation(A3), [code]]
convSec = Dict{String, String}()    # converting nation's sectors; code, sector
concMat = Dict{String, conTab}()    # concordance matrix sets: abbreviation, conTab
concMatNorm = Dict{String, conTabNorm}()    # normalized concordance matrix sets: abbreviation, conTab

function readConcMatFile(natFile, sectorFile, concFile, convNat; weight=false)

    global names, nations, convSec, eorCodes

    f = open(natFile)
    readline(f)     # read title
    for l in eachline(f)
        s = string.(strip.(split(l, '\t')))
        names[s[1]] = s[2]
        nations[s[2]] = nation(s[1], s[2], parse(Int, s[3]), parse(Bool, lowercase(s[4])), "")
    end
    close(f)

    f = open(sectorFile)
    readline(f)     # read title
    pre_nat = ""
    si = 0
    for l in eachline(f)
        s = strip.(split(l, '\t'))
        n = s[3]
        if haskey(nations, n)
            if pre_nat != n
                si = 0
                eorCodes[n] = Array{String, 1}()
            end
            if s[4] == "Commodities" || (!nations[n].hasComEn && s[4] == "Industries")
                si += 1
                c = string(si)
                push!(eorCodes[n], c)
                nations[n].sectors[c] = sector(n, c, s[5])
            end
            pre_nat = n
        end
    end
    close(f)

    f = open(concFile)
    readline(f)     # read title
    for l in eachline(f)
        s = string.(strip.(split(l, '\t')))
        n, c = names[s[1]], s[2]
        push!(nations[n].sectors[c].linked, s[4])
        if length(s) > 5 && tryparse(Float64, s[6]) != nothing
            push!(nations[n].sectors[c].weight, parse(Float64, s[6]))
        else push!(nations[n].sectors[c].weight, 1.0)
        end
        if !haskey(convSec, s[4]); convSec[s[4]] = s[5] end
    end
    close(f)
end

function readXlsxData(inputFile, convNat; weight=false)

    global totals, names, nations, convSec, natCodes, eorCodes

    xf = XLSX.readxlsx(inputFile)

    # read all nations' abstract information
    sh = xf["Abstract"]
    nc = length(sh["A"])
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r) == 1; continue end
        names[r[2]] = r[3]
        nations[r[3]] = nation(r[2], r[3], r[4], r[5], r[6])
        totals += r[4]
    end

    # read converting nation's sectors
    sh = xf[convNat]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r) == 1; continue end
        if r[1] == convNat; convSec[string(r[2])] = r[3] end
    end
    natCodes = sort(collect(keys(convSec)))

    # read sector data
    for n in sort(collect(keys(nations)))
        sh = xf[nations[n].matchCode]
        eorCodes[n] = []
        lc = ""
        for r in XLSX.eachrow(sh)
            if XLSX.row_number(r) == 1
                if r[1] != "No."; println(n,": Heading error.") end
            elseif r[2] == nations[n].matchCode[2:end]
                nations[n].sectors[string(r[3])] = sector(r[2], string(r[3]), r[4])
                if !(string(r[3]) in eorCodes[n]) push!(eorCodes[n], string(r[3])) end
                lc = string(r[3])
            elseif r[2] == convNat
                if convSec[string(r[3])] == r[4]
                    push!(nations[n].sectors[lc].linked, string(r[3]))
                    if weight; push!(nations[n].sectors[lc].weight, parse(Float64, string(r[5])))
                    else push!(nations[n].sectors[lc].weight, 1)
                    end
                else
                    println(n,"\t", r[1], "\t", r[2], "\t", convSec[string(r[3])], "\t", r[4], "\tsectors do not match")
                end
            else
                println(n,"\t", r[1], "\t", r[2], "\tsource error.")
            end
        end
    end

    close(xf)

    return nations
end

function exportLinkedSectors(outputFile, nation; mrio="Eora")

    global nations, convSec, eorCodes

    f = open(outputFile, "w")
    println(f, "Nation\t"*mrio*"_code\t"*mrio*"category\t"*nation*"_code\t"*nation*"_category\tWeight")
    for n in sort(collect(keys(nations)))
        for sc in eorCodes[n]
            s = nations[n].sectors[sc]
            for i = 1:length(s.linked)
                c = s.linked[i]
                w = s.weight[i]
                println(f, nations[n].name, "\t", s.code, "\t", s.categ, "\t", c, "\t", convSec[c], "\t", w)
            end
        end
    end
    close(f)
end

function buildConMat()  # build concordance matrix for all countries

    global nations, concMat, natCodes, eorCodes

    for n in collect(keys(nations))
        concMat[n] = conTab(nations[n].ns, length(natCodes))
        for idxEor = 1:length(eorCodes[n])
            s = nations[n].sectors[eorCodes[n][idxEor]]
            for i = 1:length(s.linked)
                l = s.linked[i]
                w = s.weight[i]
                idxNat = findfirst(x -> x==l, natCodes)
                concMat[n].conMat[idxEor, idxNat] += w
                concMat[n].sumEora[idxEor] += w
                concMat[n].sumNat[idxNat] += w
            end
        end
    end

    cm = Dict{String, Array{Int,2}}()
    for n in collect(keys(concMat)); cm[n] = concMat[n].conMat end
    return cm
end

function addSubstSec(substCodes, subDict, secDict)
    # rebuild concordance matrix adding substitute sectors: {[substitute code], Dict[subst.code, [sub.code]], Dict[code, sector]}

    global nations, concMat, natCodes, eorCodes, convSec

    nc = length(natCodes)   # without substitution codes
    append!(natCodes, substCodes)
    for c in substCodes; convSec[c] = secDict[c] end
    ntc = length(natCodes)  # with substitution codes

    for n in collect(keys(nations))
        ctab = conTab(nations[n].ns, ntc)
        ctab.conMat[:,1:nc] = concMat[n].conMat
        ctab.sumNat[1:nc] = concMat[n].sumNat

        for i=1:length(substCodes)
            subcds = subDict[substCodes[i]]
            for j=1:length(subcds); ctab.conMat[:,nc+i] += concMat[n].conMat[:,findfirst(x->x==subcds[j], natCodes)] end
            ctab.sumNat[nc+i] = sum(ctab.conMat[:,nc+i])
        end
        for i=1:nations[n].ns; ctab.sumEora[i] = sum(ctab.conMat[i,:]) end

        concMat[n] = ctab
    end
end

function normConMat() # normalize concordance matrix

    global concMatNorm, natCodes, eorCodes
    nnc = length(natCodes)
    for n in collect(keys(nations))
        concMatNorm[n] = conTabNorm(nations[n].ns, nnc)
        for i = 1:nnc
            if concMat[n].sumNat[i] > 1
                for j = 1:nations[n].ns; concMatNorm[n].conMat[j, i] = concMat[n].conMat[j, i] / concMat[n].sumNat[i] end
            elseif concMat[n].sumNat[i] == 1
                for j = 1:nations[n].ns; concMatNorm[n].conMat[j, i] = concMat[n].conMat[j, i] end
            else println(n,"\tsum of ",natCodes[i]," is ", concMat[n].sumNat[i], ": concordance matrix value error")
            end
        end

        for i = 1:nnc; for j = 1:nations[n].ns; concMatNorm[n].sumNat[i] += concMatNorm[n].conMat[j, i] end end
        for i = 1:nations[n].ns; for j = 1:nnc; concMatNorm[n].sumEora[i] += concMatNorm[n].conMat[i, j] end end
    end

    cmn = Dict{String, Array{Float64,2}}()
    for n in collect(keys(concMatNorm)); cmn[n] = concMatNorm[n].conMat end
    return cmn
end

function printConMat(outputFile, convNat = ""; norm = false, categ = false)

    global natCodes, eorCodes, nations

    f = open(outputFile, "w")
    tmpEor = sort(collect(keys(names)))

    #File printing
    print(f, "Eora/"*convNat*"\t")
    for s in natCodes; print(f, "\t", s) end
    println(f, "\tSum")

    for n in tmpEor
        abb = names[n]
        for i = 1:length(eorCodes[abb])
            if categ; print(f, abb, "\t", nations[abb].sectors[eorCodes[abb][i]].categ)
            else print(f, abb, "\t", nations[abb].sectors[eorCodes[abb][i]].code)
            end
            for j = 1:length(natCodes)
                if !norm; print(f, "\t", concMat[abb].conMat[i,j])
                elseif norm; print(f, "\t", concMatNorm[abb].conMat[i,j])
                end
            end
            if !norm; println(f, "\t", concMat[abb].sumEora[i])
            elseif norm; println(f, "\t", concMatNorm[abb].sumEora[i])
            end
        end
    end

    close(f)
end

function printSumNat(outputFile, convNat = ""; norm = false)

    global natCodes

    f = open(outputFile, "w")
    tmpEor = sort(collect(keys(names)))

    #File printing
    print(f, "Eora/"*convNat)
    for s in natCodes; print(f, "\t", s) end
    println(f)

    for n in sort(collect(keys(nations)))
        print(f, n)
        if !norm; for s in concMat[n].sumNat; print(f, "\t", s) end
        elseif norm; for s in concMatNorm[n].sumNat; print(f, "\t", s) end
        end
        println(f)
    end

    close(f)
end

end
