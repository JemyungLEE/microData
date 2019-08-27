module ConcMatBuilder

# Developed date: 26. Aug. 2019
# Last modified date: 27. Aug. 2019
# Subject: classification concordance matrix builder, India-Eora
# Description: matching India commodity and Eora industry classifications and build a concordance matrix
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

mutable struct codes
    source::String              # India: India classification, Eora: Eora classification
    code::String                # classification codes of each Eora and India classification
    categ::String               # category of Eora industry or India commodity classification
    linked::Array{codes,1}      # list of Eora-India classfication links
    codes() = new()
end

function readClassCodes(inputFile)
    f = open(inputFile)

    eorTag = "Eora"
    indTag = "India"

    eorCodeList = codes[]
    indCodeList = codes[]
    indClassList = []
    eorClassList = []

    for l in eachline(f)
        c = codes()
        c.linked = []
        no, c.source, c.code, c.categ = split(l, '\t')
        if c.source == eorTag
            push!(eorCodeList, c)
            push!(eorClassList, c.code)
        elseif c.source == indTag
            if !(c.code in indClassList)
                push!(indClassList, c.code)
                push!(indCodeList, c)
            else
                c = indCodeList[findfirst(x -> x==c.code, indClassList)]
            end

            push!(c.linked, eorCodeList[end])
            push!(eorCodeList[end].linked, c)
        else
            println("tag error.")
        end
    end

    return eorCodeList, indCodeList, eorClassList, indClassList

end

function makeConcMat(eorCodeList, eorClass, indClass)

    ne = length(eorClass)
    ni = length(indClass)

    conMat = zeros(Float16, ni, ne)      # concordance matrix between Eora-India commodity classifications

    for c in eorCodeList
        idxEor = findfirst(x -> x==c.code, eorClass)
        for l in c.linked
            idxInd = findfirst(x -> x==l.code, indClass)
            conMat[idxInd, idxEor] += 1.0
        end
    end

    return conMat

end

end
