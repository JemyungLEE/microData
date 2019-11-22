module XLSXextractor

import XLSX

# Developed date: 1. Oct. 2019
# Last modified date: 22. Nov. 2019
# Subject: XLSX data extractor and Concordance matrix builder
# Description: read sector matching information from a XLSX file and build concordance matrix
#              bewteen converting nation and Eora accounts
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

mutable struct sector       # category data
    source::String          # 3-digit abbreviation of nation, cf) "_SD": standard 26 categcories, "_EU": EU 61 categcories
    code::Int               # classification codes of each Eora classification
    categ::String           # category of Eora industry or commodity classification
    linked::Array{Int16,1}  # linked converting nation's sector codes

    function sector(src, cod, cat)
        new(src, cod, cat, [])
    end
end

mutable struct nation       # nation data
    name::String
    abb::String             # abbreviation of country name
    ns::Int                 # number of sectors
    hasComEn::Bool          # wether the nation have 'Commodities'-entity data
    matchCode::String       # matching nation code of this nation. cf) "_SD": standard 26 categcories, "_EU": EU 61 categcories
    sectors::Array{sector,1}    # Eora sectors

    nation(n::String, a::String, ns::Int, has::Bool, mc::String) = new(n, a, ns, has, mc, [])
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
    conMat::Array{Float32,2}  # concordance matrix
    sumEora::Array{Float32,1} # sums of eora sectors
    sumNat::Array{Float32,1}  # sums of converting nation's sectors

    function conTabNorm(eorSecNum, natSecNum)
        new(zeros(Int, eorSecNum, natSecNum), zeros(Int, eorSecNum), zeros(Int, natSecNum))
    end
end

totals = 0  # total sectors
names = Dict{String, String}()      # Full names, abbreviation
nations = Dict{String, nation}()    # abbreviation, nation
convSec = Dict{Int16, String}()     # converting nation's sectors; code, sector
concMat = Dict{String, conTab}()    # concordance matrix sets: abbreviation, conTab
concMatNorm = Dict{String, conTabNorm}()    # normalized concordance matrix sets: abbreviation, conTab

function readXlsxData(inputFile, convNat)

    global totals, nations, convSec

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
        if r[1] == convNat; convSec[r[2]] = r[3] end
    end

    # read sector data
    for n in sort(collect(keys(nations)))
        sh = xf[nations[n].matchCode]
        for r in XLSX.eachrow(sh)
            if XLSX.row_number(r) == 1
                if r[1] != "No."; println(n,": Heading error.") end
            elseif r[2] == nations[n].matchCode[2:end]
                push!(nations[n].sectors, sector(r[2], r[3], r[4]))
            elseif r[2] == convNat
                if convSec[r[3]] == r[4]
                    push!(nations[n].sectors[end].linked, r[3])
                else
                    println(n,"\t", r[1], "\t", r[2], "\t", convSec[r[3]], "\t", r[4], "\tsectors do not match")
                end
            else
                println(n,"\t", r[1], "\t", r[2], "\tsource error.")
            end

        end
    end

    # check read data
    #=
    for n in sort(collect(keys(nations)))
        for c in nations[n].sectors
            print(n,"\t",c.source,"\t",c.code,"\t",c.categ)
            for s in c.linked; print("\t", s) end
            println()
        end
    end
    =#

    close(xf)

    return nations
end

function buildConMat()  # build concordance matrix for all countries in the XLSX file

    global concMat
    tmpSec = sort(collect(keys(convSec)))

    for n in collect(keys(nations))
        concMat[n] = conTab(nations[n].ns, length(tmpSec))
        for s in nations[n].sectors
            idxEor = s.code
            for l in s.linked
                idxNat = findfirst(x -> x==l, tmpSec)
                concMat[n].conMat[idxEor, idxNat] += 1
                concMat[n].sumEora[idxEor] += 1
                concMat[n].sumNat[idxNat] += 1
            end
        end
    end

    cm = Dict{String, Array{Int,2}}()
    for n in collect(keys(concMat)); cm[n] = concMat[n].conMat end
    return cm
end

function normConMat() # normalize concordance matrix

    global concMatNorm

    for n in collect(keys(nations))
        concMatNorm[n] = conTabNorm(nations[n].ns, length(keys(convSec)))
        for i = 1:length(keys(convSec))
            if concMat[n].sumNat[i] > 1
                for j = 1:nations[n].ns
                    concMatNorm[n].conMat[j, i] = concMat[n].conMat[j, i] / concMat[n].sumNat[i]
                end
            elseif concMat[n].sumNat[i] == 1
                for j = 1:nations[n].ns
                    concMatNorm[n].conMat[j, i] = concMat[n].conMat[j, i]
                end
            else println(n,"\t", concMat[n].sumNat[i], "\tconcordance matrix value error")
            end
        end

        for i = 1:length(keys(convSec))
            for j = 1:nations[n].ns
                concMatNorm[n].sumNat[i] += concMatNorm[n].conMat[j, i]
            end
        end

        for i = 1:nations[n].ns
            for j = 1:length(keys(convSec))
                concMatNorm[n].sumEora[i] += concMatNorm[n].conMat[i, j]
            end
        end
    end

    cmn = Dict{String, Array{Float32,2}}()
    for n in collect(keys(concMatNorm)); cmn[n] = concMatNorm[n].conMat end
    return cmn
end

function printConMat(outputFile, convNat = "", norm = false)
    f = open(outputFile, "w")
    tmpSec = sort(collect(keys(convSec)))
    tmpEor = sort(collect(keys(names)))

    #File printing
    print(f, "Eora/"*convNat*"\t")
    for s in tmpSec
        print(f, "\t", s)
    end
    println(f, "\tSum")

    for n in tmpEor
        abb = names[n]
        for i = 1:length(nations[abb].sectors)
            print(f, abb, "\t", nations[abb].sectors[i].code)
            for j = 1:length(tmpSec)
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

function printSumNat(outputFile, convNat = "", norm = false)
    f = open(outputFile, "w")
    tmpSec = sort(collect(keys(convSec)))
    tmpEor = sort(collect(keys(names)))

    #File printing
    print(f, "Eora/"*convNat)
    for s in tmpSec
        print(f, "\t", s)
    end
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
