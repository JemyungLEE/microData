module ConcMatBuilder

# Developed date: 14. Apr. 2021
# Last modified date: 5. Jul. 2021
# Subject: Build concordance matric between MRIO and HBS/CES micro-data
# Description: read sector matching information from a XLSX/TXT/CSV file and
#              build concordance matrix bewteen converting nation and Eora accounts
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX
using Statistics

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
deCodes = Array{String, 1}()        # direct emission sector code list
convSec = Dict{String, String}()    # converting nation's sectors; code, sector
concMatIe = Dict{String, conTab}()  # concordance matrix sets: abbreviation, conTab
concMatIeNorm = Dict{String, conTabNorm}()  # normalized concordance matrix sets: abbreviation, conTab
concMatDe = Array{Float64, 2}(undef, 0, 0)  # concordance matrix sets for direct emission

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

function getCommodityCodes(code_source)

    global natCodes

    if isa(code_source, Array{String, 1}); natCodes = code_source
    elseif isa(code_source, String)
        f_sep = getValueSeparator(code_source)
        f = open(code_source)
        i = findfirst(x->x=="Code", string.(strip.(split(readline(f), f_sep))))
        for l in eachline(f)
            code = string.(strip.(split(l, f_sep)))[i]
            if !(code in natCodes) push!(natCodes, code) end
        end
        close(f)
    end


end

function buildDeConcMat(nation, deCodeFile, concFile; norm = false, output = "", energy_wgh = false, de_data = [], de_year = 0)

    global natCodes, deCodes, concMatDe

    f_sep = getValueSeparator(deCodeFile)
    f = open(deCodeFile)
    i = findfirst(x->x=="DE_code", string.(strip.(split(readline(f), f_sep))))
    for l in eachline(f); push!(deCodes, string.(strip.(split(l, f_sep)))[i]) end
    close(f)

    nd, nc = length(deCodes), length(natCodes)
    concMatDe = zeros(Float64, nd, nc)
    if energy_wgh
        de_ener = de_data.de_energy[de_year]
        ene_avg = [mean(filter(x->x>0, [de_ener[n][j] for n in collect(keys(de_ener))])) for j=1:nd]
    end

    essential = ["Nation", "DE_code", "CES/HBS_code", "Weight"]
    f_sep = getValueSeparator(concFile)
    f = open(concFile)
    title = string.(strip.(split(readline(f), f_sep)))
    i = [findfirst(x->x==et, title) for et in essential]
    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        if s[i[1]] == nation
            decode, cescode, wgh = s[i[2]], s[i[3]], parse(Float64, s[i[4]])
            di, ci = findfirst(x->x==decode, deCodes), findfirst(x->x==cescode, natCodes)
            if energy_wgh; concMatDe[di, ci] += wgh * (de_ener[nation][di]>0 ? de_ener[nation][di] : ene_avg[di])
            else concMatDe[di, ci] += wgh
            end
        end
    end
    close(f)

    if norm
        for i = 1:nc
            ces_sum = sum(concMatDe[:,i])
            if ces_sum > 0; concMatDe[:, i] /= ces_sum end
        end
    end

    if length(output)>0
        mkpath(rsplit(output, '/', limit = 2)[1])
        f = open(output, "w")
        for sc in natCodes; print(f, "\t", sc) end; println(f)
        for i = 1:nd
            print(f, deCodes[i])
            for j = 1:nc; print(f, "\t", concMatDe[i, j]) end
            println(f)
        end
        close(f)
    end

    return concMatDe
end

function readIeConcMatFile(natFile, sectorFile, concFile; weight=false)

    global names, nations, convSec, natCodes, eorCodes

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
        n, c = s[1], s[2]
        push!(nations[n].sectors[c].linked, s[4])
        if length(s) > 5 && tryparse(Float64, s[6]) != nothing
            push!(nations[n].sectors[c].weight, parse(Float64, s[6]))
        else push!(nations[n].sectors[c].weight, 1.0)
        end
        if !haskey(convSec, s[4]); convSec[s[4]] = s[5] end
    end
    close(f)
end

function readXlsxData(inputFile, convNat; weight=false, nat_label = "")

    global totals, names, nations, convSec, natCodes, eorCodes
    if length(nat_label) == 0; nat_label = convNat end

    xf = XLSX.readxlsx(inputFile)

    # read all nations' abstract information
    sh = xf["Abstract"]
    nc = length(sh["A"])
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r) == 1; continue end
        nat = string(strip(r[3]))
        names[string(strip(r[2]))] = nat
        nations[nat] = nation(string(strip(r[2])), nat, r[4], r[5], string(strip(r[6])))
        totals += r[4]
    end

    # read converting nation's sectors
    sh = xf[convNat]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r) == 1; continue end
        if r[1] == nat_label; convSec[string(strip(string(r[2])))] = string(strip(string(r[3]))) end
    end
    natCodes = sort(collect(keys(convSec)))

    # read sector data
    for n in sort(collect(keys(nations)))
        sh = xf[nations[n].matchCode]
        eorCodes[n] = []
        lc = ""
        for r in XLSX.eachrow(sh)
            na = string(strip(string(r[2])))
            sc = string(strip(string(r[3])))
            ds = string(strip(string(r[4])))
            if XLSX.row_number(r) == 1
                if r[1] != "No."; println(n,": Heading error.") end
            elseif na == nations[n].matchCode[2:end]
                nations[n].sectors[sc] = sector(na, sc, ds)
                if !(sc in eorCodes[n]) push!(eorCodes[n], sc) end
                lc = sc
            elseif na == nat_label
                if convSec[sc] == ds
                    push!(nations[n].sectors[lc].linked, sc)
                    if weight; push!(nations[n].sectors[lc].weight, parse(Float64, string(r[5])))
                    else push!(nations[n].sectors[lc].weight, 1)
                    end
                else
                    println(n,"\t", r[1], "\t", na, "\t", convSec[sc], "\t", ds, "\tsectors do not match")
                end
            else
                println(n,"\t", r[1], "\t", na, "\t", nat_label, "\tsource error.")
            end
        end
    end

    close(xf)

    return nations
end

function exportLinkedSectors(outputFile, nation; mrio="Eora")

    global nations, convSec, eorCodes

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    println(f, "Nation\t"*mrio*"_code\t"*mrio*"category\t"*nation*"_code\t"*nation*"_category\tWeight")
    for n in sort(collect(keys(nations)))
        for sc in eorCodes[n]
            s = nations[n].sectors[sc]
            for i = 1:length(s.linked)
                c = s.linked[i]
                w = s.weight[i]
                println(f, n, "\t", s.code, "\t", s.categ, "\t", c, "\t", convSec[c], "\t", w)
            end
        end
    end
    close(f)
end

function buildIeConMat()  # build concordance matrix for all countries

    global nations, concMatIe, natCodes, eorCodes

    for n in collect(keys(nations))
        concMatIe[n] = conTab(nations[n].ns, length(natCodes))
        for idxEor = 1:length(eorCodes[n])
            s = nations[n].sectors[eorCodes[n][idxEor]]
            for i = 1:length(s.linked)
                l = s.linked[i]
                w = s.weight[i]
                idxNat = findfirst(x -> x==l, natCodes)
                concMatIe[n].conMat[idxEor, idxNat] += w
                concMatIe[n].sumEora[idxEor] += w
                concMatIe[n].sumNat[idxNat] += w
            end
        end
    end

    cm = Dict{String, Array{Int,2}}()
    for n in collect(keys(concMatIe)); cm[n] = concMatIe[n].conMat end
    return cm
end

function addSubstSec(substCodes, subDict, secDict)
    # rebuild concordance matrix adding substitute sectors: {[substitute code], Dict[subst.code, [sub.code]], Dict[code, sector]}

    global nations, concMatIe, natCodes, eorCodes, convSec

    nc = length(natCodes)   # without substitution codes
    append!(natCodes, substCodes)
    for c in substCodes; convSec[c] = secDict[c] end
    ntc = length(natCodes)  # with substitution codes

    for n in collect(keys(nations))
        ctab = conTab(nations[n].ns, ntc)
        ctab.conMat[:,1:nc] = concMatIe[n].conMat
        ctab.sumNat[1:nc] = concMatIe[n].sumNat

        for i=1:length(substCodes)
            subcds = subDict[substCodes[i]]
            for j=1:length(subcds); ctab.conMat[:,nc+i] += concMatIe[n].conMat[:,findfirst(x->x==subcds[j], natCodes)] end
            ctab.sumNat[nc+i] = sum(ctab.conMat[:,nc+i])
        end
        for i=1:nations[n].ns; ctab.sumEora[i] = sum(ctab.conMat[i,:]) end

        concMatIe[n] = ctab
    end
end

function normConMat() # normalize concordance matrix

    global concMatIeNorm, natCodes, eorCodes
    nnc = length(natCodes)
    for n in collect(keys(nations))
        concMatIeNorm[n] = conTabNorm(nations[n].ns, nnc)
        for i = 1:nnc
            if concMatIe[n].sumNat[i] > 1
                for j = 1:nations[n].ns; concMatIeNorm[n].conMat[j, i] = concMatIe[n].conMat[j, i] / concMatIe[n].sumNat[i] end
            elseif concMatIe[n].sumNat[i] == 1
                for j = 1:nations[n].ns; concMatIeNorm[n].conMat[j, i] = concMatIe[n].conMat[j, i] end
            else println(n,"\tsum of ",natCodes[i]," is ", concMatIe[n].sumNat[i], ": concordance matrix value error")
            end
        end

        for i = 1:nnc; for j = 1:nations[n].ns; concMatIeNorm[n].sumNat[i] += concMatIeNorm[n].conMat[j, i] end end
        for i = 1:nations[n].ns; for j = 1:nnc; concMatIeNorm[n].sumEora[i] += concMatIeNorm[n].conMat[i, j] end end
    end

    cmn = Dict{String, Array{Float64,2}}()
    for n in collect(keys(concMatIeNorm)); cmn[n] = concMatIeNorm[n].conMat end
    return cmn
end

function printConMat(outputFile, convNat = ""; norm = false, categ = false)

    global natCodes, eorCodes, nations

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
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
                if !norm; print(f, "\t", concMatIe[abb].conMat[i,j])
                elseif norm; print(f, "\t", concMatIeNorm[abb].conMat[i,j])
                end
            end
            if !norm; println(f, "\t", concMatIe[abb].sumEora[i])
            elseif norm; println(f, "\t", concMatIeNorm[abb].sumEora[i])
            end
        end
    end

    close(f)
end

function printSumNat(outputFile, convNat = ""; norm = false)

    global natCodes

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    tmpEor = sort(collect(keys(names)))

    #File printing
    print(f, "Eora/"*convNat)
    for s in natCodes; print(f, "\t", s) end
    println(f)

    for n in sort(collect(keys(nations)))
        print(f, n)
        if !norm; for s in concMatIe[n].sumNat; print(f, "\t", s) end
        elseif norm; for s in concMatIeNorm[n].sumNat; print(f, "\t", s) end
        end
        println(f)
    end

    close(f)
end

end
