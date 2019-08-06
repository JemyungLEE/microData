module CommClassReader

# Developed date: 31. Jul. 2019
# Modified date: 6. Aug. 2019
# Subject: commodity classification reader, KR
# Description: read and store the commodity classification of Bank of Korea
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

mutable struct codes
    code::String
    desc::String
    depth::Int8
    upp::String             # upper-category
    sub::Array{String,1}    # sub-categories
    codes() = new()
end

function readCommCode(inputFile)
    f = open(inputFile)
    tag = ["Large-sized", "Medium-sized", "Small-sized", "Aggregation", "Large-Medium-Link"]

    idx = 0
    codeList = []

    for l in eachline(f)
        if l in tag
            idx = findfirst(x -> x==l, tag)
        elseif !isempty(l) && !startswith(l, "Code")
            if idx < 5
                tmpCode = codes()
                tmpCode.code, tmpCode.desc = split(l, '\t')

                if idx == 4
                    tmpCode.depth = 0
                    tmpCode.upp = "-1"
                    tmpCode.sub = ["-1"]
                elseif idx < 4
                    tmpCode.depth = idx
                    if idx == 1
                        tmpCode.upp = "-1"
                        tmpCode.sub = []
                    elseif idx == 2
                        tmpCode.sub = []
                    elseif idx == 3
                        tmpCode.upp = tmpCode.code[1:end-1]
                        tmpCode.sub = ["-1"]
                        for c in codeList
                            if c.depth == 2 && c.code == tmpCode.upp
                                push!(c.sub, tmpCode.code)
                            end
                        end
                    end
                end
                push!(codeList, tmpCode)
            elseif idx == 5
                upper, lower = split(l, '\t')
                for c in codeList
                    if c.depth == 1 && c.code == upper
                        push!(c.sub, lower)
                    elseif c.depth == 2 && c.code == lower
                        c.upp = upper
                    end
                end
            end
        end
    end
    close(f)

    for c in codeList
        println(c)
    end
end

end
