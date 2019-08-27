# Developed date: 27. Aug. 2019
# Last modified date: 27. Aug. 2019
# Subject: classification category analyzer, India-Eora
# Description: match categories of India commodity classificartion and Eora industry classification
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
include("ConcMatBuilder.jl")
using .ConcMatBuilder
cmb = ConcMatBuilder

inputFile = "India_Eora_singleLink.txt"
#inputFile = "India_Eora_multipleLinks.txt"
inputFile = Base.source_dir()*"/"*inputFile

eorCodeList, indCodeList, eorClass, indClass = cmb.readClassCodes(inputFile)
conMat = cmb.makeConcMat(eorCodeList, eorClass, indClass)

for cls in indClass
println(cls)
end

#println(conMat)
#=
for c in eorCodeList
    print(c.source,"\t",c.code,"\t",c.categ,"\t")
    for lc in c.linked
        print(lc.code,"\t",lc.categ,"\t")
    end
    println()
end

for c in indCodeList
    print(c.source,"\t",c.code,"\t",c.categ,"\t")
    for lc in c.linked
        print(lc.code,"\t",lc.categ,"\t")
    end
    println()
end
=#
