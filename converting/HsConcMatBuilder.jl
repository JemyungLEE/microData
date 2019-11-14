module HsConcMatBuilder

import XLSX

# Developed date: 17. Oct. 2019
# Last modified date: 17. Oct. 2019
# Subject: Harmonized System (HS) commodity classification codes
#           and Household consumption micr-data codes concordance matrix builing
# Description: read sector matching information from a XLSX file and build concordance matrix
#              bewteen converting nation and HS codes
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

mutable struct sector       # category data
    class::String           # nation or H3
    code::String            # classification codes of each Eora classification
    categ::String           # category of Eora industry or commodity classification
    level::Int8             # whether 2, 4, or 6-digit code: just for HS codes
    parent::String          # parent code: just for HS codes
    linked::Array{String,1}  # linked converting nation's sector codes

    sector(cls::String, cod::String, cat::String, lv::Int8, par::String, lin=[]) = new(cls, cod, cat, lv, par, lin)
end

mutable struct conTab       # concordance tables
    conMat::Array{Int,2}    # concordance matrix
    sumHS::Array{Int,1}     # sums of HS sectors
    sumNat::Array{Int,1}    # sums of converting nation's sectors

    function conTab(hsSecNum, natSecNum)
        new(zeros(Int, hsSecNum, natSecNum), zeros(Int, hsSecNum), zeros(Int, natSecNum))
    end
end

mutable struct conTabNorm       # normalized concordance tables
    conMat::Array{Float32,2}    # concordance matrix
    sumHS::Array{Float32,1}     # sums of HS sectors
    sumNat::Array{Float32,1}    # sums of converting nation's sectors

    function conTabNorm(hsSecNum, natSecNum)
        new(zeros(Int, hsSecNum, natSecNum), zeros(Int, hsSecNum), zeros(Int, natSecNum))
    end
end

convSec = Dict{String, String}()     # converting nation's sectors; code, sector
hsSec = Dict{String, sector}()       # HS classification's sectors; code, sector

function readXlsxData(inputFile, convNat, HSsheet)  # input filename, converting nation, HS code sheet name

    global convSec, hsSec
    hsCode = []

    xf = XLSX.readxlsx(inputFile)

    # read converting nation's sectors
    sh = xf[convNat]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r) == 1; continue end
        convSec[string(r[1])] = r[2]
    end

    # read sector data
    sh = xf[HSsheet]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r) == 1
            if r[1] != "Source"
                println("Heading error: ",r[1],"\t", r[2],"\t", r[3],"\t", r[4],"\t", r[5])
            end
        elseif r[1] == "H3"
            hsSec[r[2]] = sector(r[1], r[2], r[3], parse(Int8,r[5]), r[4])
            push!(hsCode, r[2])
        elseif r[1] == convNat
            if convSec[string(r[2])] == r[3]
                push!(hsSec[hsCode[end]].linked, string(r[2]))
            else
                println(n,"\t", r[1], "\t", r[2], "\t", convSec[r[2]], "\t", r[3], "\tsectors do not match")
            end
        else
            println(n,"\t", r[1], "\t", r[2], "\t", r[3], "\tsource error.")
        end
    end

    close(xf)

    return hsSec, convSec
end

function buildConMat()  # build concordance matrix, row: HS, column: nation

    ct = conTab(length(keys(hsSec)), length(keys(convSec)))

    natCodes = sort(collect(keys(convSec)))
    hsCodes = sort(collect(keys(hsSec)))

    for c in hsCodes
        idxHS = findfirst(x -> x==c, hsCodes)
        for l in hsSec[c].linked
            idxNat = findfirst(x -> x==l, natCodes)
            ct.conMat[idxHS, idxNat] += 1
            ct.sumHS[idxHS] += 1
            ct.sumNat[idxNat] += 1
        end
    end

    return ct
end

function normConMat(ct) # normalize concordance matrix of ct (concordance tables), row: HS, column: nation

    nhs = length(keys(hsSec))
    ncs = length(keys(convSec))

    nct = conTabNorm(nhs, ncs)

    for i = 1:ncs
        if ct.sumNat[i] > 1
            for j = 1:nhs
                nct.conMat[j, i] = ct.conMat[j, i] / ct.sumNat[i]
            end
        elseif ct.sumNat[i] == 1
            for j = 1:nhs
                nct.conMat[j, i] = ct.conMat[j, i]
            end
        elseif ct.sumNat[i] != 0
            println(i, "\t", ct.sumNat[i], "\tconcordance matrix value error")
        end
    end

    for i = 1:ncs
        for j = 1:nhs
            nct.sumNat[i] += nct.conMat[j, i]
        end
    end

    for i = 1:nhs
        for j = 1:ncs
            nct.sumHS[i] += nct.conMat[i, j]
        end
    end

    return nct

end

function printConMat(outputFile, ct =[], convNat = "")


    f = open(outputFile, "w")
    natCodes = sort(collect(keys(convSec)))
    hsCodes = sort(collect(keys(hsSec)))

    #File printing

    print(f, "HS/"*convNat)
    for c in natCodes
        print(f, "\t", c)
    end
    println(f, "\tSum")

    for i = 1:length(hsCodes)
        print(f, hsCodes[i])
        for j = 1:length(natCodes)
            print(f, "\t", ct.conMat[i,j])
        end
        println(f, "\t", ct.sumHS[i])
    end

    print(f, "Sum")
    for i = 1:length(natCodes)
        print(f, "\t",ct.sumNat[i])
    end
    println(f)

    close(f)

end

end
