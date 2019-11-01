module MicroDataReader

# Developed date: 21. Oct. 2019
# Last modified date: 31. Oct. 2019
# Subject: India Household Consumer Expenditure microdata reader
# Description: read and store specific data from India microdata, integrate the consumption data from
#              different files, and export the data as DataFrames
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

mutable struct item
    code::String        # item code
    value::Float64      # total consumption value
    quantity::Float64   # total consumption quantity
    period::Int16       # during the last 30 or 365 days

    source::Int8        # 1:purchase, 2:home-grown, 3:both, 4:free collection, 5:exchange, 6:gifts/charities
    homeVal::Float64    # home produce value
    homeQt::Float64     # home produce quantity

    item(c, v=0, q=0, p=0, s=0, hv=0, hq=0) = new(c, v, q, p, s, hv, hq)
end

mutable struct member
    age::Int16
    sex::Int8   # [1]male, [2]female
    mar::Int8   # [1]never, [2]married, [3]widowed, [4]divorced/seperated
    edu::Int8   # [1]non literate, [2-4]literate [5]below primary [6]primary, [7]middle, [8]secondary
                # [10]higher secondary, [11]diploma/certificate, [12]graduate, [13]post gradutate/above
    member(a=0, s=0, m=0, e=0) = new(a, s, m, e)
end

mutable struct household
    id::String          # Household identification no.
    date::String        # Survey date
    fsu::String         # Fisrt Stage Unit serial no.
    state::String       # State code
    district::String    # District code
    sector::String      # urban or rural
    size::Int16         # household size
    mpceUrp::Float64    # monthly per capita expenditure (uniform reference period)
    mpceMrp::Float64    # monthly per capita expenditure (mixed reference period)

    members::Array{member,1}    # household member(s)
    items::Array{item,1}        # consumed items

    household(i,da="",fa="",st="",di="",sec="",si=0,mu=0,mm=0,me=[],it=[])=new(i,da,fa,st,di,sec,si,mu,mm,me,it)
end

global households = Dict{String, household}()

function readHouseholdData(hhData)
    # Read household identification data, index: [1]hhid, [2]date, [3]fsu, [4]state, [5]district, [6]sector
    f = open(hhData[1][1])
    idx = hhData[1][2]
    for l in eachline(f)
        s = split(l, '\t')
        id = s[idx[1]]
        if !haskey(households, id)
            households[id] = household(id, s[idx[2]], s[idx[3]], s[idx[4]], s[idx[5]], s[idx[6]])
        elseif haskey(households, id)
            println("Household ID duplication: ", l)
        end
    end
    close(f)

    # Read household size data, index: [1]hhid, [2]household size
    f = open(hhData[2][1])
    for l in eachline(f)
        s = split(l, '\t')
        id = s[hhData[2][2][1]]
        if haskey(households, id)
            households[id].size = parse(Int16, s[hhData[2][2][2]])
        elseif !haskey(households, id)
            println("Household ID absence: ", id)
        end
    end
    close(f)

    # Read household monthly per capita expenditure, index: [1]hhid, [2]MPCE URP, [3]MPCE MRP
    f = open(hhData[3][1])
    for l in eachline(f)
        s = split(l, '\t')
        id = s[hhData[3][2][1]]
        if haskey(households, id)
            households[id].mpceUrp = parse(Float64, s[hhData[3][2][2]])
            households[id].mpceMrp = parse(Float64, s[hhData[3][2][3]])
        elseif !haskey(households, id)
            println("Household ID absence: ", id)
        end
    end
    close(f)

    # Read household family members, index: [1]hhid, [2]age, [3]sex, [4]marriage, [5]education
    f = open(hhData[4][1])
    for l in eachline(f)
        s = split(l, '\t')
        id = s[hhData[4][2][1]]
        m = member()
        if haskey(households, id)
            if length(s[hhData[4][2][2]])>0; m.age = parse(Int16, s[hhData[4][2][2]])
            else m.age = -1 end
            if length(s[hhData[4][2][3]])>0; m.sex = parse(Int8, s[hhData[4][2][3]])
            else m.sex = -1 end
            if length(s[hhData[4][2][4]])>0; m.mar = parse(Int8, s[hhData[4][2][4]])
            else m.mar = -1 end
            if length(s[hhData[4][2][5]])>0; m.edu = parse(Int8, s[hhData[4][2][5]])
            else m.edu = -1 end

            push!(households[id].members, m)
        elseif !haskey(households, id)
            println("Household ID absence: ", id)
        end
    end
    close(f)

    return households
end

function readMicroData(mdata)
    # index: [1]id, [2]item code, [3]value, [4]quantity, [5]period,
    #        [6]source, [7]home produce value, [8]home produce quantity
    # A value of '-1' means 'no data'.

    for m in mdata
        f = open(m[1])
        for l in eachline(f)
            s = split(l, '\t')
            id = s[m[2][1]]
            if haskey(households, id)
                i = item(s[m[2][2]])
                if m[2][3]>0 && length(s[m[2][3]])>0; i.value = parse(Float64, s[m[2][3]])
                else i.value = -1 end
                if m[2][4]>0 && length(s[m[2][4]])>0; i.quantity = parse(Float64, s[m[2][4]])
                else i.quantity = -1 end
                i.period = m[2][5]
                if length(m[2])>5 && length(s[m[2][6]])>0; i.source = parse(Int8, s[m[2][6]])
                else i.source = -1 end
                if length(m[2])>6 && length(s[m[2][7]])>0; i.homeVal = parse(Float64, s[m[2][7]])
                else i.homeVal = -1 end
                if length(m[2])>7 && length(s[m[2][8]])>0; i.homeQt = parse(Float64, s[m[2][8]])
                else i.homeQt = -1 end

                push!(households[id].items, i)
            elseif !haskey(households, id)
                println("Household ID absence: ", id)
            end
        end
        close(f)
    end

    return households
end

function printHouseholdData(outputFile)
    f = open(outputFile, "w")
    sd = ["rural", "urban"]
    count = 0

    println(f, "HHID\tSurvey Date\tFSU\tState\tDistrict\tSector\tHH size\tMPCE_URP\tMPCE_MRP")
    for hhid in sort(collect(keys(households)))
        h = households[hhid]
        sector = sd[parse(Int8, h.sector)]
        println(f, h.id,"\t",h.date,"\t",h.fsu,"\t",h.state,"\t",h.district,"\t",sector,"\t",h.size,"\t",h.mpceUrp,"\t",h.mpceMrp)
        count += 1
    end
    close(f)

    println("$count households' data is printed.")
end

function printMemberData(outputFile)
    f = open(outputFile, "w")
    count = 0

    sd = ["male", "female"]
    md = ["never married", "currently married", "widowed", "divorced/seperated"]
    ed = ["non literate", "literate without formal schooling: through EGC/NFEC/AEC",
            "literate without formal schooling: through TLC", "literate without formal schooling: Others",
            "literate with formal schooling:below primary", "literate with formal schooling:primary",
            "literate with formal schooling:middle", "literate with formal schooling:secondary", "",
            "literate with formal schooling:higher secondary",
            "literate with formal schooling:diploma/certificate course",
            "literate with formal schooling:graduate",
            "literate with formal schooling:post graduate and above"]

    println(f, "HHID\tAge\tSex\tMarital status\tEducation")
    for hhid in sort(collect(keys(households)))
        h = households[hhid]
        for m in h.members
            if m.sex>0; sex = sd[m.sex]; else sex = ""; end
            if m.mar>0; mar = md[m.mar]; else mar = ""; end
            if m.edu>0; edu = ed[m.edu]; else edu = ""; end

            println(f, h.id,"\t",m.age,"\t",sex,"\t",mar,"\t",edu)
            count += 1
        end
    end
    close(f)

    println("$count members' data is printed.")
end

function printMicroData(outputFile)
    f = open(outputFile, "w")
    count = 0

    sd = ["only purchase", "only home-grown stock", "both purchase and home-grown stock",
          "only free collection", "only exchange of goods and services", "only gifts/charities", "", "", "others"]

    println(f, "HHID\tItem code\tTotal value\tTotal quantity\tConsumption period\tConsumption source\tHome produce value\tHome produce quantity")
    for hhid in sort(collect(keys(households)))
        h = households[hhid]
        for i in h.items
            if i.source>0; source = sd[i.source]; else source = ""; end

            println(f, h.id,"\t",i.code,"\t",i.value,"\t",i.quantity,"\t",i.period,"\t",source,"\t",i.homeVal,"\t",i.homeQt)
            count += 1
        end
    end
    close(f)

    println("$count items' data is printed.")
end

end
