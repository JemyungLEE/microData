module MicroDataReader

# Developed date: 21. Oct. 2019
# Last modified date: 15. Jan. 2020
# Subject: India Household Consumer Expenditure microdata reader
# Description: read and store specific data from India microdata, integrate the consumption data from
#              different files, and export the data as DataFrames
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using DataFrames

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

    totExp::Float64     # household's total one-year expenditure calculated by the item data
    totExpMrp::Float64  # household's total one-year expenditure calculated by MPCE_MRP

    regCds::Array{String, 1}   # additional region codes: [State_region, District, Stratum, Substratum_No, FOD_Sub_Region]

    household(i,da="",fa="",st="",di="",sec="",si=0,mu=0,mm=0,me=[],it=[],te=0,tem=0,rcd=[]) = new(i,da,fa,st,di,sec,si,mu,mm,me,it,te,tem,rcd)
end

global households = Dict{String, household}()
global categories = Dict{String, String}()    # expenditure category: {code, description}

function readHouseholdData(hhData, tag="")
    # Read household identification data, index: [1]hhid, [2]date, [3]fsu, [4]state, [5]district, [6]sector
    f = open(hhData[1][1])
    idx = hhData[1][2]
    for l in eachline(f)
        s = split(l, '\t')
        id = tag * s[idx[1]]
        if !haskey(households, id)
            households[id] = household(id, s[idx[2]], s[idx[3]], s[idx[4]], s[idx[5]], s[idx[6]])
            if length(hhData[1][2])>6
                households[id].regCds = [string(s[idx[7][1]]),string(s[idx[7][2]]),string(s[idx[7][3]]),string(s[idx[7][4]]),string(s[idx[7][5]])]
            end
        elseif haskey(households, id)
            println("Household ID duplication: ", l)
        end
    end
    close(f)

    # Read household size data, index: [1]hhid, [2]household size
    f = open(hhData[2][1])
    for l in eachline(f)
        s = split(l, '\t')
        id = tag * s[hhData[2][2][1]]
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
        id = tag * s[hhData[3][2][1]]
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
        id = tag * s[hhData[4][2][1]]
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

function readMicroData(mdata, tag="")
    # index: [1]id, [2]item code, [3]value, [4]quantity, [5]period,
    #        [6]source, [7]home produce value, [8]home produce quantity
    #        [9]another period, (period, start item code, end item code)
    # A value of '-1' means 'no data'.

    for m in mdata
        f = open(m[1])
        for l in eachline(f)
            s = split(l, '\t')
            id = tag * s[m[2][1]]
            if haskey(households, id)
                i = item(s[m[2][2]])
                if m[2][3]>0 && length(s[m[2][3]])>0; i.value = parse(Float64, s[m[2][3]])
                else i.value = -1 end
                if m[2][4]>0 && length(s[m[2][4]])>0; i.quantity = parse(Float64, s[m[2][4]])
                else i.quantity = -1 end
                if length(m[2])>5 && length(s[m[2][6]])>0; i.source = parse(Int8, s[m[2][6]])
                else i.source = -1 end
                if length(m[2])>6 && length(s[m[2][7]])>0; i.homeVal = parse(Float64, s[m[2][7]])
                else i.homeVal = -1 end
                if length(m[2])>7 && length(s[m[2][8]])>0; i.homeQt = parse(Float64, s[m[2][8]])
                else i.homeQt = -1 end
                if length(m[2])>8 && string(m[2][9][2]) <= i.code <= string(m[2][9][3]); i.period = m[2][9][1]
                else i.period = m[2][5] end

                push!(households[id].items, i)
            elseif !haskey(households, id)
                println("Household ID absence: ", id)
            end
        end
        close(f)
    end

    return households
end

function currencyExchange(exchangeRate)     # exchangeRate: Rupees to USD currency exchange rate
                                            # Dict[MMYY] or Dict[YY] are also applicable

    if typeof(exchangeRate) <: Number
        for hhid in sort(collect(keys(households)))
            h = households[hhid]
            for i in h.items
                i.value *= exchangeRate
                i.homeVal *= exchangeRate
            end
        end
    elseif typeof(exchangeRate) <: AbstractDict
        for hhid in sort(collect(keys(households)))
            h = households[hhid]
            if haskey(exchangeRate, h.date[3:6]); er = exchangeRate[h.date[3:6]]
            elseif haskey(exchangeRate, h.date[5:6]); er = exchangeRate[h.date[5:6]]
            else println("Exchange rate error: no exchange rate data for ", h.date[5:6], " year")
            end

            for i in h.items
                i.value *= exchangeRate
                i.homeVal *= exchangeRate
            end
        end
    end

end

function readCategory(inputFile)
    f= open(inputFile)

    readline(f)
    for l in eachline(f)
        s = split(l, '\t')
        categories[s[1]] = s[2]
    end

    close(f)
end

function makeExpenditureMatrix(outputFile = "")
    # build expenditure matrix with 365 days consumption monetary values

    row = sort(collect(keys(households)))
    col = sort(collect(keys(categories)))

    mat = zeros(Float64, length(row), length(col))
    rowErr = zeros(Int32, length(row))
    colErr = zeros(Int32, length(col))

    # make expenditure matrix, row: households, col: expenditure categories
    for r = 1:length(row)
        h = households[row[r]]
        total = 0
        for i in h.items
            if haskey(categories, i.code)
                c = findfirst(x -> x==i.code, col)
                if i.value > 0
                    val = i.value
                    if i.period != 365; val *= 365 / i.period end
                    mat[r,c] += val
                    total += val
                else rowErr[r] += 1; colErr[c] += 1
                end
            end
        end
        h.totExp = total
        h.totExpMrp = h.mpceMrp * h.size * 12
    end

    # print expenditure matrix as a file
    if length(outputFile) > 0
        f = open(outputFile, "w")

        for c in col; print(f, "\t", c) end
        print(f, "\tRow_error")
        println(f)
        for r = 1:length(row)
            print(f, row[r])
            for c = 1:length(col); print(f, "\t", mat[r,c]) end
            print(f, "\t", rowErr[r])
            println(f)
        end
        print(f, "Column error")
        for c = 1:length(col); print(f, "\t", colErr[c]) end
        print(f, "\t", sum(colErr))
        println(f)
        close(f)
    end

    return mat, row, col, rowErr, colErr
end

function convertHouseholdData(outputFile = "")
    id = String[]; dat = String[]; fsu = String[]; sta = String[]; dis = String[]
    sec = String[]; siz = Int16[]; urp = Float64[]; mrp = Float64[]
    tot = Float64[]; totMrp = Float64[]

    avg = Float16[]; mal = Int32[]; fem = Int32[]; chi = Int32[]; mid = Int32[]; old = Int32[]

    for hhid in sort(collect(keys(households)))
        h = households[hhid]

        agesum = 0
        total = 0
        male = 0
        female = 0
        child = 0
        grownup = 0
        aged = 0

        for m in h.members
            total += 1
            agesum += m.age
            if m.sex == 1; male += 1
            elseif m.sex == 2; female +=1
            end
            if 0<= m.age < 15; child += 1
            elseif 15<= m.age <65; grownup += 1
            elseif m.age >= 65; aged += 1
            end
        end

        push!(id, h.id)
        push!(dat, h.date)
        push!(fsu, h.fsu)
        push!(sta, h.state)
        push!(dis, h.district)
        push!(sec, h.sector)
        push!(siz, h.size)
        push!(urp, h.mpceUrp)
        push!(mrp, h.mpceMrp)
        push!(tot, h.totExp)
        push!(totMrp, h.totExpMrp)
        push!(avg, agesum / total)
        push!(mal, male)
        push!(fem, female)
        push!(chi, child)
        push!(mid, grownup)
        push!(old, aged)
    end

    df = DataFrame(HHID=id, Survey_date=dat, FSU=fsu, State=sta, District=dis, Sector=sec, Size=siz,
                    MPCE_Urp=urp, MPCE_Mrp=mrp, Total_Exp=tot, Total_Exp_Mrp=totMrp,
                    Avg_age=avg, Male=mal, Female=fem, Children=chi, Grownups=mid, Aged_persons=old)

    if length(outputFile) > 0
        f = open(outputFile, "w")
        println(f, df)
        close(f)
    end

    return df
end


function printHouseholdData(outputFile, addCds=false)
    f = open(outputFile, "w")
    sd = ["rural", "urban"]
    count = 0

    print(f, "HHID\tSurvey Date\tFSU\tState\tDistrict\tSector\tHH size\tMPCE_URP\tMPCE_MRP")
    if addCds; print(f, "\tState_code\tDistrict_code\tStratum\tSubstratum_No\tFOD_Sub_Region") end
    println(f)
    for hhid in sort(collect(keys(households)))
        h = households[hhid]
        sector = sd[parse(Int8, h.sector)]
        print(f, h.id,"\t",h.date,"\t",h.fsu,"\t",h.state,"\t",h.district,"\t",sector,"\t",h.size,"\t",h.mpceUrp,"\t",h.mpceMrp)
        if addCds; print(f, "\t", h.regCds[1],"\t",h.regCds[2],"\t",h.regCds[3],"\t",h.regCds[4],"\t",h.regCds[5]) end
        println(f)
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

function initVars()
    global households = Dict{String, household}()
    global categories = Dict{String, String}()
end

end
