# Developed date: 31. Jul. 2019
# Modified date:
# Subject: commodity classification reader, KR
# Description: read and store the commodity classification of Bank of Korea
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

function readCommCode(inputFile)
    f = open(inputFile)

    tag = ["Index", "Large-sized", "Large-Medium-Link", "Medium-sized", "Small-sized"]

    i = 1
    for l in eachline(f)
        if i<=length(tag) && occursin(tag[i], l)
            println(i, "\t", tag[i], "\t", l)
            i += 1
        end
    end

    close(f)
end

inputFile = "/Users/leejimac/Desktop/Workspace/Julia/IOT/Comm_Class_KR.txt"
readCommCode(inputFile)
