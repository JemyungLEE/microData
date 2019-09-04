# Developed date: 30. Aug. 2019
# Last modified date: 4. Sep. 2019
# Subject: Eora sector classification
# Description: classify the categories of Eora sectors by nation
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

module SectorClassifier

struct sector
    index::Int
    country::String
    abb::String         # abbreviation of country name
    entity::String      # 'Industries' or 'Commodities'
    categ::String       # categeory of this sector

    function sector(str::String)
        idx, name, abb, ent, sec = split(str, '\t')
        new(parse(Int, idx), name, abb, ent, strip(sec))
    end
end

mutable struct nation
    name::String
    abb::String         # abbreviation of country name
    ns::Int16           # number of sectors
    hasComEn::Bool      # wether the nation have 'Commodities'-entity data
    isSimple::Bool      # wether a nation's sectors are classified by the basic 26 categories
    sectors::Array{String,1}

    nation(nstr::String, astr::String) = new(nstr, astr, 0, false, false, [])
end

nc = 0     # number of counties
nations = Dict{String, nation}()

function readSectorData(inputFile)
    f = open(inputFile)
    cmm = "Commodities"
    ind = "Industries"
    global nc
    global nations

    # check nations and wether they have commodity accounts
    readline(f)
    for l in eachline(f)
        s = sector(l)
        if !haskey(nations, s.abb)
            nations[s.abb] = nation(s.country, s.abb)
        elseif s.entity == cmm
                nations[s.abb].hasComEn = true
        end
    end

    # store data
    seekstart(f)
    readline(f)
    for l in eachline(f)
        s = sector(l)
        if haskey(nations, s.abb)
            n = nations[s.abb]
            if (n.hasComEn && s.entity == cmm) || (!n.hasComEn && s.entity == ind)
                push!(n.sectors, s.categ)
            elseif !n.hasComEn && s.entity == cmm
                println("Data entity error: ", l)
            end
        else
            println("Nation key(A3 abbriviation) wrong data: ", l)
        end
    end
    close(f)

    # check number of sectors by nation
    for n in keys(nations)
        nations[n].ns = length(nations[n].sectors)
        if nations[n].ns == 26
            nations[n].isSimple = true
        end
    end
    nc = length(nations)
end

function printNationData(outputFile = "")
    if isempty(outputFile)
        for key in sort(collect(keys(nations)))
            n = nations[key]
            println(key, '\t', n.name, '\t', n.abb, '\t', n.ns, '\t', n.hasComEn, '\t', n.isSimple)
        end
    else
        idx = 1
        f = open(outputFile, "w")
        println(f, "Index\tName\tAbb\tN_sectors\tHas commodity\tIs simple")
        for key in sort(collect(keys(nations)))
            n = nations[key]
            println(f, idx, '\t', n.name, '\t', n.abb, '\t', n.ns, '\t', n.hasComEn, '\t', n.isSimple)
            idx += 1
        end
        println("Written File: ", outputFile)
        close(f)
    end
end

function compareSectors(comparedFile = "", containedFile = "")
    compared = fill(false, nc, nc)
    compMat = fill(0, nc, nc)       # comparison results matrix: same or not
    contMat = fill(0, nc, nc)       # comparison results matrix: amount of co-existing sectors
    nationList = []

    col = 0
    for kcol in sort(collect(keys(nations)))
        push!(nationList, kcol)
        col += 1
        row = 0
        colSectors = sort(nations[kcol].sectors)
        for krow in sort(collect(keys(nations)))
            row += 1
            rowSectors = sort(nations[krow].sectors)
            compared[row, col] = (colSectors == rowSectors)
            if compared[row, col] && kcol != krow
                compMat[row, col] = nations[kcol].ns
            end
            for s in rowSectors
                if s in colSectors; contMat[row, col] += 1 end
            end
        end
    end

    if !isempty(comparedFile)
        f = open(comparedFile, "w")
        for n in nationList; print(f, '\t', n) end
        println(f)
        for i = 1:nc
            print(f, nationList[i], '\t')
            for j = 1:nc
                print(f, compMat[i, j], '\t')
            end
            println(f)
        end
        println("Written File: ", comparedFile)
        close(f)
    else
        println(compMat)
    end

    if !isempty(containedFile)
        f = open(containedFile, "w")
        for n in nationList; print(f, '\t', n) end
        println(f)
        for i = 1:nc
            print(f, nationList[i], '\t')
            for j = 1:nc
                if contMat[i, j] > 10 && contMat[i, j] != nations[nationList[i]].ns
                    print(f, contMat[i, j],'/', nations[nationList[i]].ns," /", nations[nationList[j]].ns, '\t')
                else print(f, "\t")
                end
            end
            println(f)
        end
        println("Written File: ", containedFile)
        close(f)
    end

end


end
