# Developed date: 17. Oct. 2022
# Last modified date: 27. Oct. 2022
# Subject: Analysis group information
# Description: Calculate statistics, CF, or related figures of grouped information
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

clearconsole()

cd(Base.source_dir())

include("MicroDataReader.jl")
include("GroupEstimator.jl")

using .MicroDataReader
using .GroupEstimator

mdr = MicroDataReader
gre = GroupEstimator

cesYear = 2018; exchYear = cesYear
eoraYear = 2015
nation = "Indonesia"
natA3 = "IDN"
natCurr = "IDR"
hh_label = "URUT"
member_mode = false
resfile = "kor18rt_diseminasi_modified.csv"
qt_label = ["R1401","R1402","R1403","R1404","R1405","R1406","R1407","R1408","R1507","R1508","R1509","R1510A","R1510B","R1511A","R1511B"]
qt_respo = ["1","1","1","1","1","1","1","1",["1","2"],"1",["7","8","9"],"1","1",["8","9","10","11"],"1"]


filePath = Base.source_dir() * "/data/" * natA3 * "/"
indexFilePath = filePath * "index/"
microDataPath = filePath * "microdata/"
extractedPath = filePath * "extracted/"
commonIndexPath = Base.source_dir() * "/data/Common/"

natFileTag = natA3 * "_" * string(cesYear)
regInfoFile = filePath * natFileTag * "_MD_RegionInfo.txt"
cmmfile = filePath * natFileTag * "_MD_Commodities.txt"
hhsfile = filePath * natFileTag * "_MD_Households_"*natCurr*".txt"
mmsfile = filePath * natFileTag * "_MD_Members.txt"

resfile = microDataPath * resfile
sttfile = extractedPath * natA3 * "_" * string(cesYear) * "_response_stt.txt"


println("[Process]")

print(" Micro-data reading:")
print(" regions"); mdr.readPrintedRegionData(cesYear, natA3, regInfoFile)
print(", households"); mdr.readPrintedHouseholdData(cesYear, natA3, hhsfile)
if member_mode; print(", members"); mdr.readPrintedMemberData(cesYear, natA3, mmsfile) end

gre.importMicroData(mdr)
gre.readResponseData(cesYear, natA3, resfile, qst_label = qt_label, qst_res = qt_respo, hhid_label = hh_label)
gre.printResponseStatistics(cesYear, natA3, sttfile)

println(" ... completed")
