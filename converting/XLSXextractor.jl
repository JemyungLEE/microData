module XLSXextractor

import XLSX

# Developed date: 1. Oct. 2019
# Last modified date: 1. Oct. 2019
# Subject: XLSX data extractor
# Description: read sector matching information from XLSX files
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

mutable struct sector       # category data
    source::String          # 3-digit abbreviation of nation, cf) "_SD": standard 26 categcories, "_EU": EU 61 categcories
    code::Int               # classification codes of each Eora and India classification
    categ::String           # category of Eora industry or India commodity classification
    linked::Array{sector,1} # list linked converting nation's sectors

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

nc = 0     # number of countries
nations = Dict{String, nation}()
convSec = Dict{Int16, String}()    # converting nation's sectors

function readXlsxFile(inputFile, nation)

    global nc, nations, convSec

    xf = XLSX.readxlsx(inputFile)

    # read abstract information
    sh = xf["Abstract"]
    nc = length(sh["A"])
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r) == 1; continue end
        nations[r[3]] = nation(r[2], r[3], r[4], r[5], r[6])
    end

    # read converting nation's data
    sh = xf[nation]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r) == 1; continue end
        if r[1] == nation; convSec[r[2]] = r[3] end
    end

    # read each country's sector data
    for n in sort(collect(keys(nations)))
        sh = xf[nations[n].matchCode]
        for r in XLSX.eachrow(sh)
            if XLSX.row_number(r) == 1
                if r[1] != "No."; println(n,": Heading error.") end
                continue
            end
            push!(nations[n].sectors, sector(r[2], r[3], r[4]))
        end
    end


    # check read data
    #=
    for n in sort(collect(keys(nations)))
        for c in nations[n].sectors
            println(n,"\t",c.source,"\t",c.code,"\t",c.categ)
        end

    end
    =#

    close(xf)

    return nations
end


end
