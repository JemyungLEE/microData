# Developed date: 10.Oct. 2019
# Last modified date: 9. Jan. 2020
# Subject: India microdata analyzer
# Description: proceed data analysis process for India household consumption microdata
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()
cd(Base.source_dir())
include("IntegrityChecker.jl")
using .IntegrityChecker
ic = IntegrityChecker

#for type 1
#microdata = ["test_lv5.txt", 20, 23, 24]    # level 5
#microdata = ["test_lv6.txt", 20, 23, 24]    # level 6
#microdata = ["test_lv7.txt", 20, 21, 22]    # level 7 for 30 days and 365 days monetary
#microdata = ["test_lv8.txt", 20, 21, 21]    # level 8 for only 30 dats monetary
#microdata = ["test_lv9.txt", 20, 27, 34, [22,24,25,26,28,30,31,32,33]]    # level 9

#for type 2
#microdata = ["test_lv5.txt", 19, 22, 23]    # level 5
#microdata = ["test_lv6.txt", 19, 20, 21]    # level 6
#microdata = ["test_lv7.txt", 19, 20, 20]    # level 7 for 30 days and 365 days monetary
#microdata = ["test_lv8.txt", 19, 20, 20]    # level 8 for only 30 dats monetary
microdata = ["test_lv9.txt", 20, 22, 28, [22,24,25,26,27]]    # level 9

inputFile = Base.source_dir()*"/data/type_2/"*microdata[1]
integrityFile = replace(inputFile, ".txt" => "_integrity.txt")

if length(microdata) == 4
    ic.checkIntegrity(inputFile, microdata[2], microdata[3], microdata[4])
elseif length(microdata) == 5
    ic.checkIntegrity(inputFile, microdata[2], microdata[3], microdata[4], microdata[5])
end

ic.printIntegrity(integrityFile)
