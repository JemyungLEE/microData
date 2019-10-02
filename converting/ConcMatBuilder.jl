module ConcMatBuilder

# Developed date: 26. Aug. 2019
# Last modified date: 2. Oct. 2019
# Subject: classification concordance matrix builder, India-Eora
# Description: matching India commodity and Eora industry classifications and build a concordance matrix
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

mutable struct sector
    source::String              # India: India classification, Eora: Eora classification
    code::String                # classification sector of each Eora and India classification
    categ::String               # category of Eora industry or India commodity classification
    linked::Array{sector,1}      # list of Eora-India classfication links

    function sector()
        new()
        linked = []
    end
end

mutable struct nation
    name::String
    abb::String         # abbreviation of country name
    ns::Int16           # number of sectors
    hasComEn::Bool      # wether the nation have 'Commodities'-entity data
    matchCode::String   # matching nation code of this nation. cf) "_SD": standard 26 categcories, "_EU": EU 61 categcories
    sectors::Array{String,1}

    nation(n::String, a::String, ns:Int16, has::Bool, mc::String) = new(n, a, ns, has, mc, [])
end

nc = 0     # number of counties
nations = Dict{String, nation}()


function readClasssector(inputFile, eoraTag = "Eora", nationTag = "India")
    f = open(inputFile)

    global eorCodeList = sector[]
    global indCodeList = sector[]
    global indClass = []
    global eorClass = []

    for l in eachline(f)
        c = sector()
        no, c.source, c.code, c.categ = split(l, '\t')
        if c.source == eoraTag
            push!(eorCodeList, c)
            push!(eorClass, c.code)
        elseif c.source == nationTag
            if !(c.code in indClass)
                push!(indClass, c.code)
                push!(indCodeList, c)
            else
                c = indCodeList[findfirst(x -> x==c.code, indClass)]
            end

            push!(c.linked, eorCodeList[end])
            push!(eorCodeList[end].linked, c)
        elseif c.source != "Source"
            println("tag error: ", c.source)
        end
    end
    indClass = sort(indClass)

    return eorCodeList, indCodeList, eorClass, indClass

end

function makeConMat()

    ne = length(eorClass)
    ni = length(indClass)

    global conMat = zeros(Int, ni, ne)      # concordance matrix between Eora-India commodity classifications
    global sumEor = zeros(Int, ne)
    global sumInd = zeros(Int, ni)

    for c in eorCodeList
        idxEor = findfirst(x -> x==c.code, eorClass)
        for l in c.linked
            idxInd = findfirst(x -> x==l.code, indClass)
            conMat[idxInd, idxEor] += 1
            sumEor[idxEor] += 1
            sumInd[idxInd] += 1
        end
    end

    return conMat, sumEor, sumInd
end

function printConMat(outputFile)
    f = open(outputFile, "w")
    ne = length(eorClass)
    ni = length(indClass)
    total = 0

    #File print
    print(f, "India/Eora")
    for i = 1:ne
        print(f, "\t", eorClass[i])
    end
    println(f, "\tSum")
    for i = 1:ni
        print(f, indClass[i], "\t")
        for j = 1:ne
            print(f, conMat[i, j], "\t")
        end
        println(f, sumInd[i])
    end
    print(f, "Sum\t")
    for i = 1:ne
        print(f, sumEor[i], "\t")
        total += sumEor[i]
    end
    println(f, total)
    close(f)

    #Screen print
    #=
    total = 0
    for i = 1:ne
        print("\t", eorClass[i])
    end
    println("\tSum")
    for i = 1:ni
        print(indClass[i], "\t")
        for j = 1:ne
            print(conMat[i, j], "\t")
        end
        println(sumInd[i])
    end
    print(f, "Sum\t")
    for i = 1:ne
        print(sumEor[i], "\t")
        total += sumEor[i]
    end
    println(total)
    =#
end

end
