module ConcMatBuilder

# Developed date: 29. Jul. 2021
# Last modified date: 25. Aug. 2021
# Subject: Build concordance matric between MRIO and HBS micro-data
# Description: read sector matching information from a XLSX/TXT/CSV file and
#              build concordance matrix bewteen converting nation and Eora accounts
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX

mutable struct sector       # category data
    source::String          # 3-digit abbreviation of nation, cf) "_SD": standard 26 categcories, "_EU": EU 61 categcories
    code::String            # classification codes of each Eora classification
    categ::String           # category of Eora industry or commodity classification
    linked::Array{String,1} # linked converting nation's sector codes

    sector(src, cod, cat) = new(src, cod, cat, [])
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
    conMat::Array{Float64,2}  # concordance matrix
    sumEora::Array{Float64,1} # sums of eora sectors
    sumNat::Array{Float64,1}  # sums of converting nation's sectors

    function conTabNorm(eorSecNum, natSecNum)
        new(zeros(Float64, eorSecNum, natSecNum), zeros(Float64, eorSecNum), zeros(Float64, natSecNum))
    end
end

totals = 0  # total sectors
names = Dict{String, String}()                  # Full names, abbreviation
nations = Dict{String, nation}()                # abbreviation, nation
natCodes = Dict{Int, Array{String, 1}}()        # converting nation's code list: {year, {codes}}
eorCodes = Dict{String, Array{String, 1}}()     # Eora's code list: [nation(A3), [code]]
convSec = Dict{Int, Dict{String, String}}()     # converting nation's sectors; {year, {code, sector}}
concMat = Dict{Int, Dict{String, conTab}}()     # concordance matrix sets: {year, {abbreviation, conTab}}
concMatNorm = Dict{Int, Dict{String, conTabNorm}}() # normalized concordance matrix sets: {year, {abbreviation, conTab}}

function readXlsxData(year, inputFile, convNat; nat_label = "")

    global totals, names, nations, convSec, natCodes, eorCodes
    natCodes[year] = Array{String, 1}()
    convSec[year] = Dict{String, String}()
    if length(nat_label) == 0; nat_label = convNat end

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
        if r[1] == nat_label; convSec[year][string(r[2])] = string(r[3]) end
    end
    natCodes[year] = sort(collect(keys(convSec[year])))

    # read sector data
    for n in sort(collect(keys(nations)))
        sh = xf[nations[n].matchCode]
        eorCodes[n] = []
        for r in XLSX.eachrow(sh)
            if XLSX.row_number(r) == 1
                if r[1] != "No."; println(n,": Heading error.") end
            elseif r[2] == nations[n].matchCode[2:end]
                push!(nations[n].sectors, sector(string(r[2]), string(r[3]), string(r[4])))
                if !(string(r[3]) in eorCodes[n]) push!(eorCodes[n], string(r[3])) end
            elseif r[2] == nat_label
                if convSec[year][string(r[3])] == string(r[4])
                    push!(nations[n].sectors[end].linked, string(r[3]))
                else
                    println(n,"\t", r[1], "\t", r[2], "\t", convSec[year][string(r[3])], "\t", r[4], "\tsectors do not match")
                end
            else
                println(n,"\t", r[1], "\t", r[2], "\t", nat_label, "\tsource error.")
            end
        end
    end

    close(xf)

    return nations
end

function exportLinkedSectors(year, outputFile, nation; mrio="Eora")

    global nations, convSec

    f = open(outputFile, "w")
    println(f, "Nation\t"*mrio*"_code\t"*mrio*"category\t"*nation*"_code\t"*nation*"_category")
    for n in sort(collect(keys(nations))), s in nations[n].sectors, c in s.linked
        println(f, nations[n].name,"\t",s.code,"\t",s.categ,"\t",c,"\t",convSec[year][c])
    end
    close(f)

end

function buildConMat(year)  # build concordance matrix for all countries in the XLSX file

    global nations, concMat, natCodes, eorCodes
    concMat[year] = Dict{String, conTab}()

    for n in collect(keys(nations))
        concMat[year][n] = conTab(nations[n].ns, length(natCodes[year]))
        for idxEor = 1:length(nations[n].sectors)
            s = nations[n].sectors[idxEor]
            for l in s.linked
                idxNat = findfirst(x -> x==l, natCodes[year])
                concMat[year][n].conMat[idxEor, idxNat] += 1
                concMat[year][n].sumEora[idxEor] += 1
                concMat[year][n].sumNat[idxNat] += 1
            end
        end
    end

    cm = Dict{String, Array{Int,2}}()
    for n in collect(keys(concMat[year])); cm[n] = concMat[year][n].conMat end
    return cm
end

function addSubstSec(year, substCodes, subDict, secDict; exp_table = [])
    # rebuild concordance matrix adding substitute sectors: {[substitute code], Dict[subst.code, [sub.code]], Dict[code, sector]}
    global nations, concMat, natCodes, eorCodes, convSec

    nc = length(natCodes[year])   # without substitution codes
    append!(natCodes[year], substCodes[year])
    for c in substCodes[year]; convSec[year][c] = secDict[year][c] end
    ntc = length(natCodes[year])  # with substitution codes

    if size(exp_table, 1) == nc; wgh = sum(exp_table, 2) ./ sum(exp_table)
    elseif size(exp_table, 2) == nc; wgh = sum(exp_table, 1) ./ sum(exp_table)
    else wgh = ones(nc)
    end

    for n in collect(keys(nations))
        ctab = conTab(nations[n].ns, ntc)
        ctab.conMat[:,1:nc] = concMat[year][n].conMat
        ctab.sumNat[1:nc] = concMat[year][n].sumNat

        for i=1:length(substCodes[year])
            subcds = subDict[year][substCodes[year][i]]
            for j=1:length(subcds)
                idx = findfirst(x->x==subcds[j], natCodes[year])
                ctab.conMat[:,nc+i] += concMat[year][n].conMat[:,idx] * wgh[idx]
            end
            ctab.sumNat[nc+i] = sum(ctab.conMat[:,nc+i])
        end
        for i=1:nations[n].ns; ctab.sumEora[i] = sum(ctab.conMat[i,:]) end

        concMat[year][n] = ctab
    end
end

function normConMat(year) # normalize concordance matrix

    global concMatNorm, natCodes, eorCodes
    concMatNorm[year] = Dict{String, conTabNorm}()
    nnc = length(natCodes[year])
    for n in collect(keys(nations))
        concMatNorm[year][n] = conTabNorm(nations[n].ns, nnc)
        for i = 1:nnc
            if concMat[year][n].sumNat[i] > 1
                for j = 1:nations[n].ns; concMatNorm[year][n].conMat[j, i] = concMat[year][n].conMat[j, i] / concMat[year][n].sumNat[i] end
            elseif concMat[year][n].sumNat[i] == 1
                for j = 1:nations[n].ns; concMatNorm[year][n].conMat[j, i] = concMat[year][n].conMat[j, i] end
            else println(n,"\tsum of ",natCodes[year][i]," is ", concMat[year][n].sumNat[i], ": concordance matrix value error")
            end
        end

        for i = 1:nnc; for j = 1:nations[n].ns; concMatNorm[year][n].sumNat[i] += concMatNorm[year][n].conMat[j, i] end end
        for i = 1:nations[n].ns; for j = 1:nnc; concMatNorm[year][n].sumEora[i] += concMatNorm[year][n].conMat[i, j] end end
    end

    cmn = Dict{String, Array{Float64,2}}()
    for n in collect(keys(concMatNorm[year])); cmn[n] = concMatNorm[year][n].conMat end
    return cmn
end

function printConMat(year, outputFile, convNat = ""; norm = false, categ = false)

    global natCodes

    f = open(outputFile, "w")
    tmpEor = sort(collect(keys(names)))

    #File printing
    print(f, "Eora/"*convNat*"\t")
    for s in natCodes[year]; print(f, "\t", s) end
    println(f, "\tSum")

    for n in tmpEor
        abb = names[n]
        for i = 1:length(nations[abb].sectors)
            if categ; print(f, abb, "\t", nations[abb].sectors[i].categ)
            else print(f, abb, "\t", nations[abb].sectors[i].code)
            end
            for j = 1:length(natCodes[year])
                if !norm; print(f, "\t", concMat[year][abb].conMat[i,j])
                elseif norm; print(f, "\t", concMatNorm[year][abb].conMat[i,j])
                end
            end
            if !norm; println(f, "\t", concMat[year][abb].sumEora[i])
            elseif norm; println(f, "\t", concMatNorm[year][abb].sumEora[i])
            end
        end
    end

    close(f)
end

function printSumNat(year, outputFile, convNat = ""; norm = false)

    global natCodes

    f = open(outputFile, "w")
    tmpEor = sort(collect(keys(names)))

    #File printing
    print(f, "Eora/"*convNat)
    for s in natCodes[year]; print(f, "\t", s) end
    println(f)

    for n in sort(collect(keys(nations)))
        print(f, n)
        if !norm; for s in concMat[year][n].sumNat; print(f, "\t", s) end
        elseif norm; for s in concMatNorm[year][n].sumNat; print(f, "\t", s) end
        end
        println(f)
    end

    close(f)
end

end
