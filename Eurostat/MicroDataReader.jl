module MicroDataReader

# Developed date: 9. Jun. 2020
# Last modified date: 5. Oct. 2020
# Subject: EU Household Budget Survey (HBS) microdata reader
# Description: read and store specific data from EU HBS microdata, integrate the consumption data from
#              different files, and export the data
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using XLSX

mutable struct member
    hhid::String    # Household id
    nation::String  # Nation code
    birthNat::Int16    # Country of Birth, 10:national, 21:EU, 22:non-EU, 99:Not specified
    citizNat::Int16    # Country of Citizenship, 10:national, 21:EU, 22:non-EU, 99:Not specified
    residNat::Int16    # Country of Residence, 10:national, 21:EU, 22:non-EU, 99:Not specified

    gender::Int8    # 1:male, 2:female, 9: not specified
    mar::Int8       # Marital status of household member
                    # 1:Never married and never in a registered partnership, 2:Married or in a registered partnership,
                    # 3:Widowed or with registered partnership that ended with death of partner (not remarried or in new registered partnership)
                    # or Divorced or with registered partnership that was legally dissolved (not remarried or in new registered partnership,
                    # 9:Not Specified
    union::Int8     # Consensual union of household member, 1:Person living in consensual union, 2:Person not living in consensual union, 9:Not specified
    relat::Int8     # Relationship with reference person, The Eurostat Definition of the Reference Person is the following:
                    # "The member that contributes most to the Household consumption budget."
                    # 1:Reference person, 2:spouse or partner, 3:child of Reference person and/or of the spouse,
                    # 4:parent of Reference person and/or of the spouse,
                    # 5:other relative, 6:no family relationship (e.g.: resident employee), 9:not specified

    edu::Int8       # Level of studies completed by the household member,
    educur::Int     # Level of studies currently being followed by the household member,
                    # 0:No formal education or below ISCED 1, 1:ISCED 1 - Primary education, 2:ISCED 2 - Lower secondary education,
                    # 3:ISCED 3 - Upper secondary education, 4:ISCED 4 - Post-secondary non-tertiary education,
                    # 5:ISCED 5 – Tertiary education first stage, 6:ISCED 6 – Tertiary education second stage, 9:not specified

    age::String     # 5 year-range classes, 0_4, 5_9, 10_14, 15_19, 20_24, …, 80_84, 85_Inf, NA=missing

    activ::Int8     # Current activity status of household member, 1:working including with employment but temporarily absent,
                    # 2:unemployed, 3:In retirement or early retirement or has given up business,
                    # 4:Pupil, student, further training, unpaid work experience, 5:Fulfilling domestic tasks,
                    # 6:Permanently disabled, 7:In compulsory military or community service,
                    # 8:not applicable (legal age to work unfulfilled), 9:not specified
    workhrs::Int8   # Hours worked, 1:Full time, 2:Part time, 8:Not applicable, 9:not specified
    worktyp::Int8   # Type of work contract for the household member, 1:permanent job/work contract of unlimited duration,
                    # 2:temporary job/work contract of limited duration, 8:not applicable (does not work), 9:not specified
    worksec::String # Economic sector in Employmen t of household member (reflecting NACE rev 2)
    worksts::Int8   # Status in employment household member, 1:employer, 2:self employed person, 3:employee, 4:unpaid family worker,
                    # 5:apprentice, 6:persons not classified by status, 8:not applicable (legal age to work unfulfilled), 9:not specified
    occup::String   # Occupation of household member (ISCO88), May not appear in file if ISCO08 is provided,
                    # Z1:Legislators, senior officials and managers, Z2:Professionals, Z3:Technicians and associate professionals,
                    # Z4:Clerks, Z5:Service workers, shop & market sales workers, Z6:Skilled agricultural and fishery workers,
                    # Z7:Craft and related trades workers, Z8:Plant and machine operators and assemblers, Z9:Elementary occupations,
                    # Z0:Armed forces, 88:Not applicable (legal age to work unfulfilled), 99:Not specified
    occup08::String # Occupation of household member (ISCO08)

    income::Float64 # Total income from all sources (net amount) corresponding to each single member of the family (in Euro)
                    # This variable does not include any household allowances

    member(id, na="") = new(id, na, 0,0,0,0,0,0,0,0,0,"",0,0,0,"",0,"","",0)
end

mutable struct household
    hhid::String        # Household identification no.
    nation::String      # Nation code
    nuts1::String       # NUTS1 code
    district::String    # District code
    sector::String      # urban or rural
    size::Int16         # household size
    weight::Float64     # Sample weight: The weights are those supplied by the Member State
    income::Float64     # Net income (total income from all sources including non-monetary components minus income taxes)
    totexp::Float64     # Total consumption expenditure (All expenditures are reported in euro)
    domexp::Float64     # Total domestic consumption expenditure (All expenditures are reported in euro)
    abrexp::Float64     # Total abroad consumption expenditure (All expenditures are reported in euro)

    popdens::Int8       # 1:Densely populated (at least 500 inhabitants/km2),  2:Intermediate (between 100 and 499 inhabitants/km2)
                        # 3:Sparsely populated (less than 100 inhabitants/km2), 9:Not specified
    eqsize::Float64     # Equivalent size (OECD scale)
    eqmodsize::Float64  # Equivalent size (modified OECD scale)

    incomes::Array{Float64,1}
        # [1]: Income in kind from employment (wages and salaries in kind). Benefits provided within the framework of paid employment (except imputed rent)
        # [2]: Income in kind from nonsalaried activities. Including withdrawals from own garden, farm or enterprise for the household's private consumption. Excluding imputed rent
        # [3]: Imputed rent. The owners' imputed rent and that of tenants living free of charge
        # [4]: Monetary net income (total monetary income from all sources minus income taxes)
    source::Int8        # 1: wages or salary, 2: income from self-employment, 3: property income, 4: pensions and retirement benefits, 5: unemployment benefit, 6: other current benefits and other income, 9: not specified

    hhtype1::Int16  # Type of Household 1 - Age limit for children set at 16 years of age
                    # 1:one adult, 2:two adults, 3:more than 2 adults, 4:one adult with dependant children, 5:two adults with dependant children, 6:more than 2 adults with dependant children, 9:other
    hhtype2::Int16  # Type of Household 2
                    # 10:One person household, 21:Lone parent with child(ren) aged less than 25, 22:Couple without child(ren) aged less than 25, 23:Couple with child(ren) aged less than 25, 24:Couple or lone parent with child(ren) aged less than 25 and other persons living in the household, 99:other type of household
    ageprof::Array{Int,1}   # age profile, Number of persons aged [1]:less than or equal to 4, [2]:from 5 to 13, [3]:from 14 to 15,
                            # [4]:from 16 to 24, [5]:from 16 to 24 who are students, [6]:from 25 to 64, [7]:more than or equal to 65
    working::Int16      # Number of persons aged 16-64 in household who are at work, '-1' denotes 'NA'
    notworking::Int16   # Number of persons aged 16-64 in household who are unemployed or are economically inactive, '-1' denotes 'NA'
    activating::Int16   # Number of members economnically active, '-1' denotes 'NA'
    occupation::String  # Socio-economic situation of the reference person
                    # Private sector, Z1:manual worker except agriculture, Z2:non manual worker except agriculture
                    # Public sector, Z3:manual worker except agriculture, Z4:non manual worker except agriculture
                    # Other, Z5:self employed person except agriculture, Z6:farmer or agricultural worker, Z7:unemployed, Z8:retired, Z9:student or in national service
                    # 10:housewife or person engaged in a non economic activity, 11:unable to work, 88:not applicable (legal age to work unattained), 99:not specified

    pov::Bool           # [true]: expends under poverty line, [false]: expends over poverty line

    members::Array{member,1}    # household member(s)
    expends::Array{Float64,1}   # consumption expenditure tables, matching with 'heCodes' and 'heDescs'
    substExp::Dict{String, Float64} # substitute code's expenditrue value: {Substitute code, expenditure}

    household(i,na) = new(i,na,"","","",0,0,0,0,0,0,0,0,0,[],0,0,0,[],0,0,0,"",false,[],[],Dict{String, Float64}())
end

global nations = Array{String, 1}()         # nation list: {Nation code}
global nationNames = Dict{String, String}() # Nation names: {Nation code, Name}

global heCats = Dict{String, String}()      # household expenditure category: {code, description}
global heCodes = Array{String, 1}()         # household expenditure item code list
global heDescs = Array{String, 1}()         # household expenditure item description list
global heCdHrr = Array{Dict{String, String}, 1}()   # household expenditure item code hierarchy: {category depth, Dict{Sub cat., Upper cat.}}
global heSubst = Array{String, 1}()         # substitute codes list
global heRplCd = Dict{String, Array{String, 1}}()   # replaced codes: {substitute code, [replaced code]}
global substCodes = Dict{Int, Dict{String, Dict{String, String}}}()    # substitute code-matching: {year, {nation, {replaced code, substitute code}}}

global hhCodes = Array{String, 1}()         # household micro-data sector code list
global hmCodes = Array{String, 1}()         # household member micro-data sector code list

global mdata = Dict{Int, Dict{String, Dict{String, household}}}()   # HBS micro-data: {year, {nation, {hhid, household}}}
global hhsList = Dict{Int, Dict{String, Array{String, 1}}}()        # household id list: {year, {nation, {hhid}}}
global expTable = Dict{Int, Dict{String, Array{Float64, 2}}}()      # household expenditure table: {year, {nation, {hhid, category}}}


function checkDepthIntegrity(year, catFiles=[], expFiles=[], fragFile="", outputFile=[]; startDepth = 1, subst = false, fixed = false)

    global mdata, heCdHrr
    fragment = Array{Dict{String, Float64}, 1}()                    # data depth fragmentation: {depth, {nation, difference}}
    fragRate = Array{Dict{String, Array{Float64, 1}}, 1}()          # data depth fragmented rate: {depth, {nation, rates}}
    integrity = Array{Dict{String, Dict{String, Float64}}, 1}()     # data depth integrity: {depth, {nation, {code, difference}}}
    exptb = Array{Dict{String, Array{Float64, 2}}, 1}()             # expenditure tables: {depth, {nation, {hhid, expenditure}}}

    nats = Array{String, 1}()                   # nation list
    codes = Array{Array{String, 1}, 1}()        # consumption codes: {depth, {code}}
    hhids = Dict{String, Array{String, 1}}()    # hhid list: {nation, {hhid}}

    nidx = Dict{String, Array{Int, 1}}()        # epxenditure index list by nation: {nation, {exp. index}}

    if length(expFiles) == length(catFiles); nd = length(expFiles)
    else println("Expenditure file number and category file number doesn't match.")
    end

    # read categories
    for i=1:nd
        cds = Array{String, 1}()
        f = open(catFiles[i])
        readline(f)
        for l in eachline(f)
            c = string(split(l,',')[1])
            if c[5:6]=="HE"; push!(cds, c) end
        end
        close(f)
        if subst
            f = open(replace(catFiles[i], ".csv"=>"_subst.csv"))
            readline(f)
            for l in eachline(f); push!(cds, string(split(l,',')[1])) end
            close(f)
        end
        push!(codes, cds)
    end

    # read expenditure data
    for i=1:nd
        f = open(expFiles[i])
        cds = string.(split(readline(f),','))[4:end]
        idx = [findfirst(x->x==c, cds) for c in codes[i]]
        push!(exptb, Dict{String, Array{Float64, 2}}())

        nh = countlines(f); seek(f, 0); readline(f)
        cnt = 0
        for l in eachline(f)
            cnt += 1
            s = string.(split(l,','))
            if parse(Int,s[1]) != year; year = parse(Int,s[1]); println("Year has changed to be ", year) end
            if i==1
                if !(s[2] in nats)
                    push!(nats, s[2])
                    hhids[s[2]] = Array{String, 1}()
                    nidx[s[2]] = Array{Int, 1}()
                end
                push!(hhids[s[2]], s[3])
                push!(nidx[s[2]], cnt)
            else
                if !(s[3] in hhids[s[2]]); println("Depth ",i," household ",s[3]," included in ",s[2]) end
                if !(cnt in nidx[s[2]]); println("Depth ",i," index ",cnt," included in ",s[2]) end
            end

            if !haskey(exptb[i], s[2]); exptb[i][s[2]] = zeros(Float64, nh, length(codes[i])) end
            exptb[i][s[2]][cnt,:] = [parse(Float64, x) for x in s[4:end][idx]]
        end

        for n in nats; exptb[i][n] = exptb[i][n][nidx[n],:] end

        close(f)
    end

    # check fragmentation
    for i = 1:nd
        depth = i + startDepth - 1
        frag = Dict{String, Float64}()
        frat = Dict{String, Array{Float64, 1}}()    # {nation, {fragmented hhs rate, average differene (within frag. hhs), std. (within)}}
        for n in nats
            nhhs = length(hhids[n])
            frag[n] = 0
            frat[n] = Array{Float64, 1}()
            nfrag = 0           # number of data fragmented households
            fravg = 0           # average difference within fragmented households
            frstd = 0           # standard deviation of differences within fragmented households
            for j= 1:nhhs
                hhexp = mdata[year][n][hhids[n][j]].domexp
                totexp = sum(exptb[depth][n][j,1:length(codes[depth])])
                diff = abs(hhexp - totexp)
                frag[n] += diff
                if diff > 0.001 * hhexp
                    nfrag += 1
                    fravg += diff
                end
            end
            if nfrag>0; fravg /= nfrag end

            for j= 1:nhhs
                hhexp = mdata[year][n][hhids[n][j]].domexp
                totexp = sum(exptb[depth][n][j,1:length(codes[depth])])
                diff = hhexp - totexp
                if diff > 0.001 * hhexp; frstd += (diff - fravg)^2 end
            end
            if nfrag>0; frstd = sqrt(frstd/nfrag) end

            frag[n] /= nhhs
            frat[n] = [nfrag/nhhs, fravg, frstd]

        end
        push!(fragment, frag)
        push!(fragRate, frat)
    end

    # print fragmented results
    f = open(fragFile, "w")
    depthTag = ["1st", "2nd", "3rd", "4th"]
    print(f, "Nation")
    for i = 1:nd
        print(f, ",Diff_",depthTag[i + startDepth - 1])
        print(f, ",FragRate_",depthTag[i + startDepth - 1])
        print(f, ",FragAvg_",depthTag[i + startDepth - 1])
        print(f, ",FragStd_",depthTag[i + startDepth - 1])
    end
    println(f)
    for n in nats
        print(f, n)
        for i = 1:nd
            print(f, ",",fragment[i + startDepth - 1][n])
            for j=1:3; print(f, ",",fragRate[i + startDepth - 1][n][j]) end
        end
        println(f)
    end
    close(f)

    # check integrity
    for i = 1:nd-1
        depth = i + startDepth - 1
        if fixed; ul = startDepth + 7 else ul = depth + 7 end     # upper code length
        ll = depth + 8      # lower code length
        push!(integrity, Dict{String, Dict{String, Float64}}())
        for n in nats
            integrity[i][n] = Dict{String, Float64}()
            if fixed; ui = 1; else ui = i end
            for j = 1:length(codes[ui])
                c = codes[ui][j]
                if length(c)==ul; integrity[i][n][c] = sum(exptb[ui][n][:,j]) end
            end
            for k = 1:length(codes[i+1])
                c = codes[i+1][k]
                integrity[i][n][c[1:ul]] -= sum(exptb[i+1][n][:,k])
            end
            for j = 1:length(codes[ui])
                c = codes[ui][j]
                if length(c)==ul; integrity[i][n][c] /= length(hhids[n]) end
            end
        end
    end

    # print results
    for i = 1:length(outputFile)
        uc = []
        if fixed; ui = 1; else ui = i end
        for n in nats
            if length(uc)==0 || uc == sort(collect(keys(integrity[i][n]))); uc = sort(collect(keys(integrity[i][n])))
            else println("Nation ", n, " has different codes in depth ", i)
            end
        end
        nc = length(uc)
        f = open(outputFile[i], "w")
        print(f, "Nation"); for n in nats; print(f, ",", n) end; println(f)
        for c in uc
            print(f, c)
            for n in nats; print(f, ",", integrity[i][n][c]) end
            println(f)
        end
        close(f)
    end
end

function readCategory(inputFile; depth=4, catFile="", inclAbr=false)

    global hhCodes, hmCodes, nationNames
    global heCats, heCodes, heDescs, heCdHrr
    codes = []
    descs = []

    xf = XLSX.readxlsx(inputFile)
    sh = xf["Nation"]
    for r in XLSX.eachrow(sh) if XLSX.row_number(r)>1 && !ismissing(r[1]) && !ismissing(r[2]); nationNames[string(r[1])] = string(r[2]) end end
    sh = xf["Consumption_categories"]
    for r in XLSX.eachrow(sh)
        if XLSX.row_number(r)>1 && !ismissing(r[1]) && !ismissing(r[2])
            push!(codes, string(r[1]))
            push!(descs, string(r[2]))
            heCats[codes[end]] = descs[end]
        end
    end
    sh = xf["HH_code"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1 && !ismissing(r[1]); push!(hhCodes, string(r[1])) end end
    sh = xf["HM_code"]
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1 && !ismissing(r[1]); push!(hmCodes, string(r[1])) end end
    close(xf)

    predpt = 0
    for i = 1:length(codes)
        curdpt = length(codes[i]) - 7
        if curdpt == depth && codes[i] != "EUR_HE00" && codes[i] != "EUR_HJ00" && (inclAbr || codes[i][5:6]=="HE")
            push!(heCodes, codes[i])
            push!(heDescs, descs[i])
        elseif predpt < depth && curdpt <= predpt && codes[i] != codes[i-1] && codes[i-1] != "EUR_HE00" && codes[i-1] != "EUR_HJ00" && (inclAbr || codes[i-1][5:6]=="HE")
            push!(heCodes, codes[i-1])
            push!(heDescs, descs[i-1])
        end
        predpt = curdpt
    end
    if predpt < depth && (inclAbr || codes[end][5:6]=="HE"); push!(heCodes, codes[end]); push!(heDescs, descs[end]) end

    for i =1:depth-1
        tmpDic = Dict{String, String}()
        uppcat = ""
        upplen = i + 7
        sublen = i + 8
        for c in codes
            if length(c) == upplen; uppcat = c
            elseif length(c) == sublen && c[1:end-1] == uppcat; tmpDic[c] = uppcat
            end
        end
        push!(heCdHrr, tmpDic)
    end
end

function readHouseholdData(year, mdataPath; visible=false, substitute=false)

    global heCodes, hhCodes, hhsList, nations, heCdHrr, heSubst, substCodes, heRplCd
    global mdata[year] = Dict{String, Dict{String, household}}()
    global hhsList[year] = Dict{String, Array{String, 1}}()

    if substitute; substCodes[year] = Dict{String, Dict{String, String}}() end

    files = []
    if isa(mdataPath, AbstractArray); files = [x*"_HBS_hh.xlsx" for x in mdataPath]
    elseif isa(mdataPath, AbstractString); for f in readdir(mdataPath); if endswith(f, "_HBS_hh.xlsx"); push!(files, f) end end
    end
    cnt = 0
    for f in files
        nc = 0
        hhidx = []
        expidx = []
        hhdata = Dict{String, household}()
        hhs = []
        if visible; cnt += 1; print("  ", cnt,"/",length(files),": ",f) end
        if substitute; sbcd = Dict{String, String}(); hrridx = Dict{String, Int}(); hrrsum = Dict{String, Float64}() end

        if isa(mdataPath, AbstractArray); nation = f[end-13:end-12]; xf = XLSX.readxlsx(f)
        elseif isa(mdataPath, AbstractString); nation = f[1:2]; xf = XLSX.readxlsx(joinpath(mdataPath, f))
        end
        sh = xf[XLSX.sheetnames(xf)[1]]
        for r in XLSX.eachrow(sh)
            if XLSX.row_number(r) == 1
                nc = XLSX.column_bounds(r)[2]
                sectors = [string(r[i]) for i=1:nc]
                hhidx = [findfirst(x->x==hc, sectors) for hc in hhCodes]
                expidx = [findfirst(x->x==he, sectors) for he in heCodes]
                if hhidx[1] == nothing; hhidx[1] = findfirst(x->x=="new_"*hhCodes[1], sectors) end
                if substitute && length(heCdHrr)>0
                    for hrr in heCdHrr; for c in collect(keys(hrr)); hrridx[c] = findfirst(x->x==c, sectors) end end
                    for c in collect(values(heCdHrr[1])); hrridx[c] = findfirst(x->x==c, sectors) end
                end
            elseif !ismissing(r[1])
                if nc != XLSX.column_bounds(r)[2]; println(year, " ", nation, " data's column size doesn't match with: ", nc) end
                d = [string(r[i]) for i=1:nc]
                if nation == "MT"; d = [if occursin("_Inf", x); replace(x,"_Inf"=>""); else x end for x in d] end  # for Malta micro-data

                if d[hhidx[2]] != nation; println(year, " ", nation, " data doesn't match with: ", XLSX.row_number(r)) end
                if d[hhidx[3]] != string(year); println(year, " ", nation, " data's year doesn't match with: ", XLSX.row_number(r)) end
                hh = household(d[hhidx[1]],d[hhidx[2]])
                push!(hhs, d[hhidx[1]])
                if d[hhidx[7]] == "NA" || d[hhidx[7]] == "missing"; inc = 0; else inc = parse(Float64, d[hhidx[7]]) end
                if d[hhidx[8]] == "NA" || d[hhidx[8]] == "missing"; domexp = 0; else domexp = parse(Float64, d[hhidx[8]]) end
                if d[hhidx[9]] == "NA" || d[hhidx[9]] == "missing"; abrexp = 0; else abrexp = parse(Float64, d[hhidx[9]]) end

                hh.nuts1, hh.size, hh.weight, hh.income = d[hhidx[4]], parse(Int16, d[hhidx[5]]), parse(Float64, d[hhidx[6]]), inc
                hh.domexp, hh.abrexp, hh.totexp = domexp, abrexp, (domexp+abrexp)
                hh.popdens, hh.eqsize, hh.eqmodsize= parse(Int8, d[hhidx[10]]), parse(Float64, d[hhidx[11]]), parse(Float64, d[hhidx[12]])
                hh.incomes = [if d[hhidx[i]] == "NA" || d[hhidx[i]] == "missing"; 0; else parse(Float64, d[hhidx[i]]) end for i=13:16]
                hh.source, hh.hhtype1, hh.hhtype2 = parse(Int8, d[hhidx[17]]), parse(Int16, d[hhidx[18]]), parse(Int16, d[hhidx[19]])
                hh.ageprof = [if d[hhidx[i]] == "missing"; 0; else parse(Int, d[hhidx[i]]) end for i=20:26]
                if d[hhidx[27]] == "NA" || d[hhidx[27]] == "missing"; wrk = -1; else wrk = parse(Int16, d[hhidx[27]]) end
                if d[hhidx[28]] == "NA" || d[hhidx[28]] == "missing"; nwrk = -1; else nwrk = parse(Int16, d[hhidx[28]]) end
                if d[hhidx[29]] == "NA" || d[hhidx[29]] == "missing"; act = -1; else act = parse(Int16, d[hhidx[29]]) end
                hh.working, hh.notworking, hh.activating, hh.occupation = wrk, nwrk, act, d[hhidx[30]]

                hh.expends = [if he == nothing; 0; elseif d[he] == "missing"; 0; else parse(Float64, d[he]) end for he in expidx]
                hh.members = []

                hhdata[d[hhidx[1]]] = hh

                if substitute
                    for cdic in heCdHrr
                        for c in collect(keys(cdic))
                            tmpcd = cdic[c]
                            tmpstr = d[hrridx[c]]
                            if !haskey(hrrsum, tmpcd); hrrsum[tmpcd] = 0 end
                            if tmpstr != "NA" && tmpstr != "missing"; hrrsum[tmpcd] = hrrsum[tmpcd]+parse(Float64, tmpstr) end
                        end
                    end
                end
            end
        end

        if substitute
            for cdic in heCdHrr
                for c in collect(values(cdic))
                    if hrrsum[c]==0
                        tmpsum = 0
                        for r in XLSX.eachrow(sh)
                            if XLSX.row_number(r)>1 && !ismissing(r[1])
                                tmpstr = string(r[hrridx[c]])
                                if tmpstr!="NA" && tmpstr!="missing"; tmpsum += parse(Float64, tmpstr) end
                            end
                        end
                        if haskey(sbcd, c); for sc in [k for (k,v) in cdic if v==c]; sbcd[sc] = sbcd[c] end
                        elseif tmpsum>0; for sc in [k for (k,v) in cdic if v==c]; sbcd[sc] = c end
                        end
                    end
                end
            end
            if length(sbcd)>0
                substCodes[year][nation] = filter((x,y)->haskey(heCdHrr[end], x), sbcd)
                for sc in collect(values(substCodes[year][nation]))
                    if !(sc in heSubst); push!(heSubst, sc) end
                    for r in XLSX.eachrow(sh)
                        if XLSX.row_number(r)>1 && !ismissing(r[1])
                            tmpstr = string(r[hrridx[sc]])
                            if tmpstr!="NA" && tmpstr!="missing"; hhdata[string(r[1])].substExp[sc] = parse(Float64, tmpstr) end
                        end
                    end
                end
            end
        end

        if visible; println(", ", length(hhdata), " households") end
        if substitute; sort!(heSubst) end

        mdata[year][nation] = hhdata
        hhsList[year][nation] = hhs
        close(xf)
    end

    nations = sort(collect(keys(mdata[year])))

    if substitute
        for sc in heSubst; heRplCd[sc] = [] end
        for scdict in collect(values(substCodes[year]))
            for sc in sort(collect(keys(scdict))); if !(sc in heRplCd[scdict[sc]]); push!(heRplCd[scdict[sc]], sc) end end
        end
    end
end

function readMemberData(year, mdataPath; visible=false)

    global hmCodes, mdata

    files = []
    if isa(mdataPath, AbstractArray); files = [x*"_HBS_hm.xlsx" for x in mdataPath]
    elseif isa(mdataPath, AbstractString); for f in readdir(mdataPath); if endswith(f, "_HBS_hm.xlsx"); push!(files, f) end end
    end
    cnt = 0
    for f in files
        nc = 0
        hmidx = []
        if visible; cnt += 1; print("  ", cnt,"/",length(files),": ",f) end
        nm = 0

        if isa(mdataPath, AbstractArray); nation = f[end-13:end-12]; xf = XLSX.readxlsx(f)
        elseif isa(mdataPath, AbstractString); nation = f[1:2]; xf = XLSX.readxlsx(joinpath(mdataPath, f))
        end
        sh = xf[XLSX.sheetnames(xf)[1]]
        for r in XLSX.eachrow(sh)
            nm += 1
            if XLSX.row_number(r) == 1
                nc = XLSX.column_bounds(r)[2]
                sectors = [string(r[i]) for i=1:nc]
                filter!(x->x!="missing",sectors)
                nc = length(sectors)
                hmidx = [findfirst(x->x==hc, sectors) for hc in hmCodes]
                if hmidx[1] == nothing; hmidx[1] = findfirst(x->x=="new_"*hmCodes[1], sectors) end
            elseif !ismissing(r[1])
                if nc != XLSX.column_bounds(r)[2]; println(year, " ", nation, " data's column size doesn't match with: ", nc) end
                d = [string(r[i]) for i=1:nc]
                if nation == "MT"; d = [if occursin("_Inf", x); replace(x,"_Inf"=>""); else x end for x in d] end  # for Malta micro-data

                if d[hmidx[2]] != nation; println(year, " ", nation, " data doesn't match with: ", XLSX.row_number(r)) end
                if d[hmidx[3]] != string(year); println(year, " ", nation, " data's year doesn't match with: ", XLSX.row_number(r)) end
                hm = member(d[hmidx[1]], d[hmidx[2]])
                if d[hmidx[21]] == "NA" || d[hmidx[21]] == "missing"; inc = 0; else inc = parse(Float64, d[hmidx[21]]) end
                if nation != "MT";  hm.birthNat, hm.citizNat, hm.residNat = parse(Int16, d[hmidx[4]]), parse(Int16, d[hmidx[5]]), parse(Int16, d[hmidx[6]])
                elseif nation == "MT"; hm.birthNat, hm.citizNat, hm.residNat = parse(Int16, split(d[hmidx[4]],"_")[1]), parse(Int16, split(d[hmidx[5]],"_")[1]), parse(Int16, split(d[hmidx[6]],"_")[1])
                end
                hm.gender, hm.mar, hm.union, hm.relat = parse(Int8, d[hmidx[7]]), parse(Int8, d[hmidx[8]]), parse(Int8, d[hmidx[9]]), parse(Int8, d[hmidx[10]])
                hm.edu, hm.educur, hm.age = parse(Int8, d[hmidx[11]]), parse(Int, d[hmidx[12]]), d[hmidx[13]]
                if nation == "MT" && d[hmidx[14]] == "3_7"; hm.activ = 2; else hm.activ = parse(Int8, d[hmidx[14]]) end
                hm.workhrs, hm.worktyp, hm.worksec = parse(Int8, d[hmidx[15]]), parse(Int8, d[hmidx[16]]), d[hmidx[17]]
                if nation == "MT"; hm.worksts = parse(Int8, split(d[hmidx[18]],"_")[1]); else hm.worksts = parse(Int8, d[hmidx[18]]) end
                hm.occup, hm.occup08, hm.income = d[hmidx[19]], d[hmidx[20]], inc

                push!(mdata[year][nation][d[hmidx[1]]].members, hm)
            end
        end
        if visible; println(", ", nm, " members") end
    end
end

function buildExpenditureMatrix(year, outputFile=""; substitute=false)
    global nations, mdata, hhsList, heCodes, heSubst, substCodes
    global expTable[year] = Dict{String, Array{Float64, 2}}()
    ne = length(heCodes)
    ns = length(heSubst)
    if substitute; nt = ne+ns; else nt = ne end

    for n in nations
        nh = length(hhsList[year][n])
        hhdata = mdata[year][n]
        hhlist = hhsList[year][n]
        etable = zeros(Float64, length(hhsList[year][n]), nt)
        for i=1:nh; etable[i,1:ne] = mdata[year][n][hhsList[year][n][i]].expends[1:ne] end

        if substitute && haskey(substCodes[year], n)>0
            for sc in collect(values(substCodes[year][n]))
                scidx = ne + findfirst(x->x==sc, heSubst)
                for i=1:nh; etable[i,scidx] = hhdata[hhlist[i]].substExp[sc] end
            end
        end

        expTable[year][n] = etable
    end

    if length(outputFile) > 0
        f = open(outputFile, "w")
        print(f, "Year,Nation,Househod")
        for hc in heCodes; print(f, ",",hc) end
        if substitute; for sc in heSubst; print(f, ",",sc) end end
        println(f)
        for n in nations
            nh = length(hhsList[year][n])
            for i = 1:nh
                print(f, year, ",", n, ",", hhsList[year][n][i])
                for j = 1:nt; print(f, ",", expTable[year][n][i,j]) end
                println(f)
            end
        end
        close(f)
    end

    return expTable
end

function makeStatistics(year, outputFile; substitute=false)

    global nations, nationNames, hhsList, heCodes, mdata

    domIdx = findall(x->x[5:6]=="HE", heCodes)
    abrIdx = findall(x->x[5:6]=="HJ", heCodes)

    f = open(outputFile, "w")
    print(f,"Year,NC,Nation,Households,Members,Inc_PerCap,Exp_PerCap,Wgh_hhs,Wgh_mm,Wgh_IncPerCap,Wgh_ExpPerCap,ExpPerHH,ExpPerEqSiz")
    println(f, ",Exp_PerHH,Exp_PerCap,ExpTotChk")
    for n in nations
        nm = incs = exps = mexps = wghhhs = wghmms = wghincs = wghexps = eqsize = expchk = 0
        nh = length(hhsList[year][n])
        for h in hhsList[year][n]
            hh = mdata[year][n][h]
            nm += hh.size
            eqsize += hh.eqsize
            incs += hh.income
            exps += hh.domexp
            wghhhs += hh.weight
            wghmms += hh.weight * hh.size
            wghincs += hh.weight * hh.income
            wghexps += hh.weight * hh.domexp
            tmpexp = sum(hh.expends[domIdx])
            if substitute; tmpexp += sum(values(hh.substExp)) end
            mexps += tmpexp         # aggregated household expenditure
            expchk += abs(hh.domexp - tmpexp)
        end
        expPerHHs = exps / nh
        expPerEqSize = exps / eqsize
        incPerCap = incs / nm
        expPerCap = exps / nm
        wghincs /= wghmms
        wghexps /= wghmms
        expchk /= nh
        mexpPerHHs = mexps / nh
        mexpPerCap = mexps / nm

        print(f, year,",",n,",",nationNames[n],",",nh,",",nm,",",incPerCap,",",expPerCap,",",wghhhs,",",wghmms,",",wghincs,",",wghexps)
        println(f,",",expPerHHs,",",expPerEqSize,",",mexpPerHHs,",",mexpPerCap,",",expchk)
    end
    close(f)
end

function readSubstCodesCSV(inputFile)

    global nations, hhsList, mdata, substCodes, heSubst, heRplCd
    year = 0

    f = open(inputFile)
    readline(f)
    for l in eachline(f)
        s = string.(split(l, ','))
        if year != parse(Int, s[1])
            year = parse(Int, s[1])
            if !haskey(substCodes, year); substCodes[year] = Dict{String, Dict{String, String}}() end
        end
        if !haskey(substCodes[year], s[2]); substCodes[year][s[2]] = Dict{String, String}() end
        if !(s[4] in heSubst); push!(heSubst, s[4]); heRplCd[s[4]] = [] end
        substCodes[year][s[2]][s[3]] = s[4]
        push!(heRplCd[s[4]], s[3])
    end
    close(f)

    sort!(heSubst)
end

function readPrintedHouseholdData(inputFile)

    global nations = Array{String, 1}()
    global hhsList = Dict{Int, Dict{String, Array{String, 1}}}()
    global mdata = Dict{Int, Dict{String, Dict{String, household}}}()

    year = 0
    nation = ""
    f = open(inputFile)
    readline(f)
    for l in eachline(f)
        s = string.(split(l, ','))
        if year != parse(Int, s[1])
            year = parse(Int, s[1])
            mdata[year] = Dict{String, Dict{String, household}}()
            hhsList[year] = Dict{String, Array{String, 1}}()
        end
        if nation != s[2]
            nation = s[2]
            mdata[year][nation] = Dict{String, household}()
            hhsList[year][nation] = Array{String, 1}()
            push!(nations, nation)
        end
        hh = household(s[3],s[2])
        hh.nuts1,hh.size,hh.weight,hh.income = s[4],parse(Int16,s[5]),parse(Float64,s[6]),parse(Float64,s[7])
        hh.totexp, hh.domexp, hh.abrexp = parse(Float64,s[8]), parse(Float64,s[9]), parse(Float64,s[10])
        hh.popdens, hh.eqsize,hh.eqmodsize = parse(Int8,s[11]),parse(Float64,s[12]),parse(Float64,s[13])
        hh.incomes = [parse(Float64, s[i]) for i=14:17]
        hh.source,hh.hhtype1,hh.hhtype2 = parse(Int8,s[18]),parse(Int16,s[19]),parse(Int16,s[20])
        hh.ageprof = [parse(Int, s[i]) for i=21:27]
        hh.working,hh.notworking,hh.activating,hh.occupation = parse(Int16,s[28]),parse(Int16,s[29]),parse(Int16,s[30]),s[31]
        mdata[year][nation][s[3]] = hh
        push!(hhsList[year][nation], s[3])
    end
    close(f)
end

function readPrintedMemberData(inputFile)

    global nations, hhsList, mdata

    f = open(inputFile)
    readline(f)
    for l in eachline(f)
        s = string.(split(l, ','))
        hm = member(s[3], s[2])
        hm.birthNat,hm.citizNat,hm.residNat,hm.gender,hm.mar,hm.union,hm.relat = parse(Int16,s[4]),parse(Int16,s[5]),parse(Int16,s[6]),parse(Int8,s[7]),parse(Int8,s[8]),parse(Int8,s[9]),parse(Int8,s[10])
        hm.edu,hm.educur,hm.age,hm.activ,hm.workhrs,hm.worktyp,hm.worksec,hm.worksts = parse(Int8,s[11]),parse(Int,s[12]),s[13],parse(Int8,s[14]),parse(Int8,s[15]),parse(Int8,s[16]),s[17],parse(Int8,s[18])
        hm.occup,hm.occup08,hm.income = s[19],s[20],parse(Float64,s[21])
        push!(mdata[parse(Int,s[1])][s[2]][s[3]].members, hm)
    end
    close(f)
end

function readPrintedExpenditureData(inputFile; substitute=false, buildHhsExp=false)

    global nations, mdata, hhsList, heCodes, heSubst, expTable

    nc = length(heCodes)
    ns = length(heSubst)
    if substitute; nt = nc+ns; else nt = nc end

    f = open(inputFile)
    readline(f)
    for l in eachline(f)
        s = string.(split(l, ','))
        year = parse(Int,s[1]); n = s[2]; hh = s[3]
        if !haskey(expTable, year); expTable[year] = Dict{String, Array{Float64, 2}}() end
        if !haskey(expTable[year], n); expTable[year][n] = zeros(Float64, length(hhsList[year][n]), nt) end
        idx = findfirst(x -> x==hh, hhsList[year][n])
        expTable[year][n][idx,1:nt] = [parse(Float64, x) for x in s[4:nt+3]]

        if buildHhsExp
            mdata[year][n][hh].expends = [parse(Float64, x) for x in s[4:nc+3]]
            if substitute && haskey(substCodes[year], n)
                for sc in collect(values(substCodes[year][n]))
                    mdata[year][n][hh].substExp[sc] = parse(Float64, s[3+nc+findfirst(x->x==sc, heSubst)])
                end
            end
        end
    end
    close(f)
end

function printCategory(year, outputFile; substitute=false)

    global heCodes, heDescs, heSubst, substCodes

    f = open(outputFile, "w")
    println(f,"Code,Description")
    for i = 1:length(heCodes); println(f,heCodes[i],",\"",heDescs[i],"\"") end
    close(f)
    if substitute
        f = open(replace(replace(outputFile, ".csv"=>"_subst.csv"), ".txt"=>"_subst.txt"), "w")
        println(f,"Code,Description")
        for sc in heSubst; println(f,sc,",\"",heCats[sc],"\"") end
        close(f)
    end

    if substitute
        f = open(replace(outputFile, "Category_"=>"SubstituteCodes_"), "w")
        println(f, "Year,Nation,Replaced_code,Substitute_code")
        for n in sort(collect(keys(substCodes[year])))
            for c in sort(collect(keys(substCodes[year][n]))); println(f, year, ",", n, ",", c, ",", substCodes[year][n][c]) end
        end
        close(f)
    end
end

function printHouseholdData(year, outputFile)

    global nations, hhsList, mdata
    f = open(outputFile, "w")
    cnt = 0

    print(f, "Year,Nation,HHID,NUTS1,HH_size,Weight")
    print(f, ",Income,Tot_exp,Dom_exp,Abr_exp,Pop_dens,Eq_size,EqMod_size")
    print(f, ",Inc_empl,Inc_nonSal,Inc_rent,Inc_monNet,Inc_source,HHtype1,HHtype2")
    println(f, ",age_0_4,age_5_13,age_14_15,age_16_24,age_16_24_stu,age_25_64,age_65_,working,notworking,activating,occupation")
    for n in nations
        for h in hhsList[year][n]
            d = mdata[year][n][h]
            print(f, year, ",", d.nation, ",", d.hhid, ",", d.nuts1, ",", d.size, ",", d.weight)
            print(f, ",", d.income, ",", d.totexp, ",", d.domexp, ",", d.abrexp, ",", d.popdens)
            print(f, ",", d.eqsize, ",", d.eqmodsize, ",", d.incomes[1], ",", d.incomes[2], ",", d.incomes[3], ",", d.incomes[4], ",", d.source)
            print(f, ",", d.hhtype1, ",", d.hhtype2, ",", d.ageprof[1], ",", d.ageprof[2], ",", d.ageprof[3], ",", d.ageprof[4], ",", d.ageprof[5], ",", d.ageprof[6], ",", d.ageprof[7])
            print(f, ",", d.working, ",", d.notworking, ",", d.activating, ",", d.occupation)
            println(f)
            cnt += 1
        end
    end
    close(f)

    println("$cnt households' data is printed.")
end

function printMemberData(year, outputFile)

    global nations, hhsList, mdata
    f = open(outputFile, "w")
    cnt = 0

    println(f, "Year,Nation,HHID,BirthNat,CitizNat,ResidNat,Gender,Mar,Union,Relat,Edu,EduCur,Age,ActSta,WorkHrs,WorkTyp,WorkSec,WorkSta,Occup,Occup08,Income")
    for n in nations
        for h in hhsList[year][n]
            for m in mdata[year][n][h].members
                print(f, year, ",", m.nation, ",", m.hhid, ",", m.birthNat, ",", m.citizNat, ",", m.residNat, ",", m.gender, ",", m.mar, ",", m.union, ",", m.relat)
                print(f, ",", m.edu, ",", m.educur, ",", m.age, ",", m.activ, ",", m.workhrs, ",", m.worktyp, ",", m.worksec, ",", m.worksts, ",", m.occup, ",", m.occup08, ",", m.income)
                println(f)
                cnt += 1
            end
        end
    end
    close(f)
    println("$cnt members' data is printed.")
end

function exchangeExpCurrency(exchangeRate; inverse=false)
    # exchangeRate: can be a file path that contains excahnge rates, a constant value of
    #               EUR to USD currency exchange rate (USD/EUR), or a set of values of Dict[MMYY] or Dict[YY]
    global nations, hhsList, mdata, expTable

    # read exchange rate from the recieved file if 'exchangeRate' is 'String'
    if typeof(exchangeRate) <: AbstractString
        erd = Dict{Int, Float64}()
        f = open(exchangeRate)
        readline(f)
        for l in eachline(f); s = split(l, '\t'); erd[parse(Int,s[1])] = parse(Float64,s[2]) end
        close(f)
        if inverse; for x in collect(keys(erd)); erd[x] = 1/erd[x] end end
        exchangeRate = erd
    end

    # if 'exchangeRate' is a Tuple of (year, constant rate)
    if typeof(exchangeRate) <: Tuple
        year, er = exchangeRate
        for n in nations; for h in hhsList[year][n]; mdata[year][n][h].expends .*= er end end
        if length(expTable) > 0; for n in nations; expTable[year][n] .*= er end end
    # if 'exchangeRate' is a set of rates
    elseif typeof(exchangeRate) <: AbstractDict
        for year in collect(keys(mdata))
            if haskey(exchangeRate, year); er = exchangeRate[year] else println(year," year exchange rate is not on the list.") end
            for n in nations; for h in hhsList[year][n]; mdata[year][n][h].expends .*= er end end
            if length(expTable) > 0; for n in nations; expTable[year][n] .*= er end end
        end
    end
end

function convertToPPP(pppRateFile; inverse=false)
    # PPP rates: can be a file path that contains PPP rates, a constant value of
    #            EUR to USD currency exchange rate (USD/EUR), or a set of values of Dict[MMYY] or Dict[YY]
    global nations, hhsList, mdata
    ppps = Dict{Int, Dict{String, Float64}}()

    # read converting rate from the recieved file if 'pppFile' is 'String'
    f = open(pppRateFile)
    readline(f)
    for l in eachline(f)
        s = split(l, '\t')
        year = parse(Int, s[1])
        if !(year in collect(keys(ppps))); ppps[year] = Dict{String, Float64}() end
        ppps[year][s[2]] = parse(Float64,s[3])
    end
    close(f)

    for year in collect(keys(mdata))
        for n in nations
            if haskey(ppps, year) && haskey(ppps[year], n)
                ppp = ppps[year][n]
                for h in hhsList[year][n]
                    hh = mdata[year][n][h]
                    hh.income /= ppp
                    hh.totexp /= ppp
                    hh.domexp /= ppp
                    hh.abrexp /= ppp
                    hh.incomes /= ppp

                    for m in hh.members; m.income /= ppp end
                end
            else println(year," year ", n, " nation's PPP is not on the list.")
            end
        end
    end
end

function initVars()
    global mdata = Dict{Int, Dict{String, Dict{String, household}}}()
end

end
