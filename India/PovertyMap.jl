module PovertyMap

# Developed date: 5. Mar. 2020
# Last modified date: 5. Mar. 2020
# Subject: Poverty ratio estimatetion and maps generation
# Description: read poverty lines and calculate poverty ratios by state. Export map making csv files.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX

stlist = Array{String, 1}()    # State code list
stname = Dict{String, String}()    # State name {State code, name}
pline = Dict{String, Tuple{Float64, Float64}}()    # {State code, Rural poverty line, Urban poverty line}
stpop = Dict{String, Tuple{Int, Int, Int}}() # State population, {State code, population{total, rural, urban}}

stsmp = Dict{String, Array{Int, 1}}() # State sample size, {State code, sample number{total, rural, urban}}
stpov = Dict{String, Array{Int, 1}}() # State poverty population, {State code, population{total, rural, urban}}

stpovrt = Dict{String, Array{Float64, 1}}()  # State poverty rate, {State code, rates{total, rural, urban}}
stpovwt = Dict{String, Array{Float64, 1}}()  # State-population weighted poverty, {State code, weighted poverty{total, rural, urban}}

gidst = Dict{String, String}()    # states' gis_codes: {state code, gis id (GIS_1)}

households = Dict{String, Any}()

function migrateData(mdr)

    global households = mdr.households

end

function applyPovertyLine(pvlineFile, outputFile="") # pvlineFile: poverty line 'csv' file

    global households
    global stlist, stname, pline, stpop, stsmp, stpov, stpovrt, stpovwt

    f = open(pvlineFile)
    readline(f)
    for l in eachline(f)
        s = split(l, ",")
        stname[s[1]] = s[2]
        pline[s[1]] = (parse(Float64, s[3]), parse(Float64, s[4]))
        stpop[s[1]] = (parse(Int, s[6]), parse(Int, s[8]), parse(Int, s[10]))
        stsmp[s[1]] = zeros(Int, 3)
        stpov[s[1]] = zeros(Int, 3)
    end
    close(f)
    stlist = sort(collect(keys(pline)))

    for h in collect(values(households))
        if h.sector == "1"; stidx=1         # rural
        elseif h.sector == "2"; stidx=2     # urban
        else println("HH sector error: not \"urban\" nor \"rural\"")
        end

        pl = pline[h.state][stidx]
        stsmp[h.state][1] += h.size
        stsmp[h.state][stidx+1] += h.size
        if h.mpceMrp < pl
            h.pov = true
            stpov[h.state][1] += h.size
            stpov[h.state][stidx+1] += h.size
        else h.pov = false
        end
    end

    for st in stlist
        stpovrt[st] = zeros(Float64, 3)
        stpovwt[st] = zeros(Float64, 3)
        for i in [2,3]
            stpovrt[st][i] = stpov[st][i] / stsmp[st][i]
            stpovwt[st][i] = stpop[st][i] * stpovrt[st][i]
            stpovwt[st][1] += stpovwt[st][i]
        end
        stpovrt[st][1] = stpovwt[st][1] / stpop[st][1]
    end

    if length(outputFile)>0
        tag=("total", "rural", "urban")
        f = open(outputFile, "w")
        print(f, "State\tName")
        for i=1:length(tag)
            print(f, "\tPopulation_",tag[i],"\tPoverty_weighted_",tag[i],"\tPoverty_rate_",tag[i])
            print(f, "\tSample_population_",tag[i],"\tSample_poverty_population_",tag[i])
        end
        println(f)
        for st in stlist
            print(f, st,"\t",stname[st])
            for i=1:length(tag)
                print(f, "\t", stpop[st][i],"\t",stpovwt[st][i],"\t",stpovrt[st][i],"\t",stsmp[st][i],"\t",stpov[st][i])
            end
            println(f)
        end
        close(f)
    end
end

function exportPovertyMap(gidIndexFile, outputFile)

    global gidst
    global stlist, stname, stpop, stpovwt

    xf = XLSX.readxlsx(gidIndexFile)
    sh = xf["IND_stat"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1 && !ismissing(r[2]); gidst[string(r[1])] = string(r[2]) end end
    close(xf)

    f = open(outputFile, "w")
    println(f,"GID_1,Population,Under_PovertyLine,Over_PovertyLine")
    for st in stlist; println(f, gidst[st],",",stpop[st][1],",",stpovwt[st][1],",",(stpop[st][1]-stpovwt[st][1])) end
    close(f)
end

end
