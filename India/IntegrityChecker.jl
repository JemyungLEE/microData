module IntegrityChecker

# Developed date: 10. Oct. 2019
# Last modified date: 11. Oct. 2019
# Subject: India household consumption microdata checker
# Description: check each sector's data integrity reading all the lines
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

mutable struct sector       # category data
    code::String             # index number of the category item(sector)
    desc::String            # description of the sector

    total::UInt32           # total data number of the sector
    nQt::UInt32             # number lines that conatin quantity data
    nVl::UInt32             # number lines that conatin value data

    isQt::Bool              # whether the sector has any quantity data (Kg, litre, no., ot etc.)
    isVl::Bool              # whether the sector has any monetary value (Pr. or other currency)
    allQt::Bool             # whether all the data-lines of the sector have quantity data
    allVl::Bool             # whether all the data-lines of the sector have value data

    function sector(cod, des="", tot=0, nq=0, nv=0, iqt=false, ivl=false, aqt=false, avl=false)
        new(cod, des, tot, nq, nv, iqt, ivl, aqt, avl)
    end
end

sectors = Dict{String, sector}()

function checkIntegrity(inputFile, idxCd, idxQt, idxVl, idxLs=[])
    # idx:index of, Cd:code, Qt:quantity, Vl:value, idxLs:check lines that contain at least a value in the list

    global sectors

    f = open(inputFile)

    readline(f)
    for l in eachline(f)
        tmpArray = split(l, '\t')
        code = tmpArray[idxCd]

        if length(idxLs) > 0
            chk = false
            for i in idxLs
                if length(tmpArray[i]) > 0 && parse(Float32, tmpArray[i]) != 0
                    chk = true
                end
            end
        else chk = true
        end

        if !haskey(sectors, code); sectors[code] = sector(code) end

        if chk
            s = sectors[code]
            s.total += 1

            if length(tmpArray[idxQt]) > 0 && isinteger(parse(Float32, tmpArray[idxQt]))
                s.nQt += 1
            end
            if length(tmpArray[idxVl]) > 0 && isinteger(parse(Float32, tmpArray[idxVl]))
                s.nVl += 1
            end

            if length(tmpArray[idxQt]) == 0 && length(tmpArray[idxVl]) == 0
                println(tmpArray[20],"\t",tmpArray[22],"\t",tmpArray[28],"\t",tmpArray[32],"\t",tmpArray[27],"\t",tmpArray[34])
            end
        end
    end

    for code in keys(sectors)
        s= sectors[code]
        if s.nQt > 0; s.isQt = true end
        if s.nVl > 0; s.isVl = true end
        if s.nQt == s.total; s.allQt = true end
        if s.nVl == s.total; s.allVl = true end
    end
    close(f)
end

function printIntegrity(outputFile = "")

    if length(outputFile)>0; f = open(outputFile, "w")
    else f = IOBuffer()
    end

    println(f, "Code\tTotal\tQuantities\tValues\tAll_Q\tAll_V")
    for code in sort(collect(keys(sectors)))
        s = sectors[code]
        println(f, s.code, "\t", s.total, "\t", s.nQt, "\t", s.nVl, "\t", s.allQt, "\t", s.allVl)
    end

    close(f)
end

end
