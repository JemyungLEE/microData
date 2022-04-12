module MicroDataReader

# Developed date: 9. Jun. 2020
# Last modified date: 8. Apr. 2022
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
    occup::String   # Occupation of household member (ISCO08)
                    # Z1:Legislators, senior officials and managers, Z2:Professionals, Z3:Technicians and associate professionals,
                    # Z4:Clerks, Z5:Service workers, shop & market sales workers, Z6:Skilled agricultural and fishery workers,
                    # Z7:Craft and related trades workers, Z8:Plant and machine operators and assemblers, Z9:Elementary occupations,
                    # Z0:Armed forces, 88:Not applicable (legal age to work unfulfilled), 99:Not specified
    # occup::String   # Occupation of household member (ISCO88), May not appear in file if ISCO08 is provided (previous version)

    income::Float64 # Total income from all sources (net amount) corresponding to each single member of the family (in Euro)
                    # This variable does not include any household allowances

    member(id, na="") = new(id, na, 0,0,0,0,0,0,0,0,0,"",0,0,0,"",0,"",0)
end

mutable struct household
    hhid::String        # Household identification no.
    nation::String      # Nation code
    nuts1::String       # NUTS1 code
    size::Int16         # household size
    weight::Float64     # Sample weight: The weights are those supplied by the Member State
    weight_nt::Float64  # Sample weight: calculated NUTS-level weight
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
    substCds::Dict{String, String}  # expenditrue code correspoding substitute code's: {expenditure item code, substitute code}

    household(i,na) = new(i,na,"",0,0,0,0,0,0,0,0,0,0,[],0,0,0,[],0,0,0,"",false,[],[],Dict{String, Float64}(),Dict{String, String}())
end

global year_list = Array{Int, 1}()          # year list
global nations = Array{String, 1}()         # nation list: {Nation code}
global nationNames = Dict{String, String}() # Nation names: {Nation code, Name}

global heCats = Dict{Int, Dict{String, String}}()   # household expenditure category: {year, {code, description}}
global heCodes = Dict{Int, Array{String, 1}}()      # household expenditure item code list: {year, code}
global heDescs = Dict{Int, Array{String, 1}}()      # household expenditure item description list: {year, description}
global heCdHrr = Dict{Int, Array{Dict{String, String}, 1}}()    # household expenditure item code hierarchy: {year, {category depth, Dict{Sub cat., Upper cat.}}}
global heSubHrr = Dict{Int, Array{Dict{String, Array{String, 1}}, 1}}() # household expenditure item upper-code corresponding sub-codes: {year, {category depth, Dict{Upper cat., {Sub cat.}}}}
global heSubst = Dict{Int, Array{String, 1}}()      # substitute codes list: {year, code}
global heRplCd = Dict{Int, Dict{String, Array{String, 1}}}()            # replaced codes: {year, {substitute code, [replaced code]}}
global substCodes = Dict{Int, Dict{String, Dict{String, String}}}()     # substitute code-matching: {year, {nation, {replaced code, substitute code}}}
global cpCodes = Dict{Int, Array{String, 1}}()      # Eurostat household expenditure COICOP code list: {year, code}
global crrHeCp = Dict{Int, Dict{String, String}}()  # Corresponding COICOP code in the national expenditure statistics: {year, {HE_code, COICOP_code(3digit)}}
global altCp = Dict{Int, Dict{String, String}}()    # Alternative COICOP sectors: {year, {COICOP_code(original), COICOP_code(alternative)}}

global hhCodes = Dict{Int, Array{String, 1}}()      # household micro-data sector code list: {year, {code}}
global hmCodes = Dict{Int, Array{String, 1}}()      # household member micro-data sector code list: {year, {code}}

global mdata = Dict{Int, Dict{String, Dict{String, household}}}()   # HBS micro-data: {year, {nation, {hhid, household}}}
global hhsList = Dict{Int, Dict{String, Array{String, 1}}}()        # household id list: {year, {nation, {hhid}}}
global expTable = Dict{Int, Dict{String, Array{Float64, 2}}}()      # household expenditure table: {year, {nation, {hhid, category}}}
global expTableSc = Dict{Int, Dict{String, Array{Float64, 2}}}()    # scaled household expenditure table: {year, {nation, {hhid, category}}}
global expStat = Dict{Int, Dict{String, Dict{String, Float64}}}()   # Eurostat exp. statistics: {year, {nation, {COICOP_category, expenditure}}}
global expSum = Dict{Int, Dict{String, Dict{String, Float64}}}()    # HBS exp. summation: {year, {nation, {COICOP_category, expenditure}}}
global expSumSc = Dict{Int, Dict{String, Dict{String, Float64}}}()  # Scaled HBS exp. summation: {year, {nation, {COICOP_category, expenditure}}}
global expQtSc = Dict{Int, Dict{String, Dict{String, Float64}}}()   # Scaled HBS exp. quota: {year, {nation, {COICOP_category, quota}}}

global cpi_list = Dict{Int, Dict{String, Array{String, 1}}}()       # Consumption price indexes: {year, {nation, {COICOP_category}}}
global cpis = Dict{Int, Dict{String, Dict{String, Float64}}}()      # Consumption price indexes: {year, {nation, {COICOP_category, CPI}}}
global sclRate = Dict{Int, Dict{String, Dict{String, Float64}}}()   # CPI scaling rate: {year, {nation, {HBS code, rate}}}

function getValueSeparator(file_name)
    fext = rsplit(file_name, '.', limit=2)[end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

function readCPIs(years=[], cpiFile=""; idx_sep = ',', freq="A", unit="INX_A_AVG", topLev = "EU")

    global nations, cpCodes, cpi_list, cpis
    if isa(years, Int); years = [years] end
    nats, codes = Array{String, 1}(), Array{String, 1}()

    mkpath(rsplit(cpiFile, '/', limit = 2)[1])
    f_sep = getValueSeparator(cpiFile)
    f = open(cpiFile)
    yrs = parse.(Int, string.(strip.(split(readline(f), f_sep)[2:end])))
    if length(years)>0; yi = [findfirst(x->x==y, yrs) for y in years]
    else yi = 1:length(yrs)
    end

    for i in yi
        cpi_list[yrs[i]] = Dict{String, Array{String, 1}}()
        cpis[yrs[i]] = Dict{String, Dict{String, Float64}}()
        cpi_list[yrs[i]][topLev] = Array{String, 1}()
    end

    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        ts = string.(strip.(split(s[1], idx_sep)))
        s = s[2:end]

        if length(ts)>1 && (ts[1], ts[2]) == (freq, unit) && ts[3][1:2] == "CP"
            code, nat = ts[3], ts[4]
            for i in yi
                cpi = strip(filter(x -> !(x in ['b', 'c', 'd', 'e']), s[i]))
                if tryparse(Float64, cpi) != nothing
                    if !haskey(cpi_list[yrs[i]], nat); cpi_list[yrs[i]][nat] = Array{String, 1}() end
                    if !haskey(cpis[yrs[i]], nat); cpis[yrs[i]][nat] = Dict{String, Float64}() end

                    if !(code in cpi_list[yrs[i]][nat]); push!(cpi_list[yrs[i]][nat], code) end
                    if !(code in cpi_list[yrs[i]][topLev]); push!(cpi_list[yrs[i]][topLev], code) end
                    cpis[yrs[i]][nat][code] = parse(Float64, cpi)
                end
            end
            if !(nat in nats); push!(nats, nat) end
        end
    end
    close(f)

    for n in nations; if !(n in nats); println("Nation ", n, " is not included in CPI data") end end
end

function scalingByCPI(years, std_year; codeDepth=0, topLev = "EU", subst = false)

    global nations, hhsList, heCodes, heSubst, mdata, cpCodes, cpi_list, cpis, sclRate
    if isa(years, Int); years = [years] end

    for y in years
        sclRate[y] = Dict{String, Dict{String, Float64}}()

        # determine code depth that all nations have CPIs
        cds_top = filter(x->x in cpi_list[std_year][topLev], filter(x->length(x)==4, cpi_list[y][topLev]))
        nct = length(cds_top)
        if codeDepth == 0
            cds_dep = ones(Int, nct)
            for i = 1:nct
                dpt_chk = true
                while dpt_chk && cds_dep[i] < 4
                    cds_dep[i] += 1
                    cds = filter(x->x in cpi_list[std_year][topLev], filter(x->startswith(x, cds_top[i]) && length(x)==cds_dep[i]+3, cpi_list[y][topLev]))
                    if length(cds) == 0; dpt_chk = false
                    else for n in nations, c in cds; if !(c in cpi_list[y][n]) || !(c in cpi_list[std_year][n]); dpt_chk = false end end
                    end
                    if !dpt_chk; cds_dep[i] -= 1 end
                end
            end
        elseif isa(codeDepth, Int); cds_dep = [codeDepth for i = 1:nct]
        end
        cpi_cds = Array{String, 1}()
        for i = 1:nct; append!(cpi_cds, filter(x->startswith(x, cds_top[i]) && length(x)<=cds_dep[i]+3, cpi_list[y][topLev])) end
        filter!(x->x in cpi_list[std_year][topLev], cpi_cds)

        # Expenditure scaling
        dpth = Dict(cds_top .=> cds_dep)
        hecs = [heCodes[y]]
        if subst; push!(hecs, heSubst[y]) end
        cp_cds = []
        for cs in hecs
            cp_cd = []
            for c in cs
                hc = c[7:min(length(c),dpth["CP"*c[7:8]]+7)]
                idx = findfirst(x->x[3:end]==hc, cpi_cds)
                while idx == nothing
                    hc = hc[1:end-1]
                    idx = findfirst(x->x[3:end]==hc, cpi_cds)
                end
                push!(cp_cd, cpi_cds[idx])
            end
            push!(cp_cds, cp_cd)
        end
        cp_he = cp_cds[1]
        if subst
            cp_sb = cp_cds[2]
            hesb_cp = Dict(heSubst[y] .=> cp_sb)
        end

        for n in filter(x -> haskey(mdata[y], x), nations)
            cr_he = [cpis[std_year][n][c] / cpis[y][n][c] for c in cp_he]
            sclRate[y][n] = Dict(heCodes[y] .=> cr_he)
            if subst
                cr = [cpis[std_year][n][c] / cpis[y][n][c] for c in cp_sb]
                merge!(sclRate[y][n], Dict(heSubst[y] .=> cr))
                cr_sb = Dict(cp_sb .=> cr)
            end

            hhs = mdata[y][n]
            for h in hhsList[y][n]
                hhs[h].expends .*= cr_he
                if subst; for c in collect(keys(hhs[h].substExp)); hhs[h].substExp[c] *= cr_sb[hesb_cp[c]] end end
                # if subst; for c in collect(keys(hhs[h].substExp)); hhs[h].substExp[c] *= sclRate[y][n][c] end end
            end

            if haskey(expTable, y) && haskey(expTable[y], n)
                if size(expTable[y][n], 2) == length(cp_he); expTable[y][n] .*= cr_he'
                elseif subst && size(expTable[y][n], 2) == length(cp_he)+length(cp_sb); expTable[y][n] .*= [cr_he;cr]'
                else println("CPI rate list and exp_table size do not match.")
                end
            end
            if haskey(expTableSc, y) && haskey(expTableSc[y], n)
                if size(expTableSc[y][n], 2) == length(cp_he); expTableSc[y][n] .*= cr_he'
                elseif subst && size(expTableSc[y][n], 2) == length(cp_he)+length(cp_sb); expTableSc[y][n] .*= [cr_he;cr]'
                else println("CPI rate list and exp_table size do not match.")
                end
            end
        end
    end
end

function mitigateExpGap(years, statFile; subst=false, percap=false, eqsize="none", cdrepl=false, alter=false, dist_mode = false)
    # fill the expenditure differences between national expenditure accounts (COICOP) and HBS by scaling the HBS expenditures
    # cdrepl: if a sub-section's COICOP does not have a value, then the HBS value moves to its higher-order section.
    # alter: apply alternative COICOP codes if sub-cateogry has COICOP data but not HBS data
    # dist_mode: distribute upper COICOP residuals to HBS sub-categories

    global nations, hhsList, mdata, heCodes, heSubst, expTable
    global cpCodes, crrHeCp, expTableSc, expStat, expSum, altCp
    if isa(years, Int); years = [years] end

    corrCds = Dict{Int, Dict{String, Array{String, 1}}}()
    for y in years
        expTableSc[y] = Dict{String, Array{Float64, 2}}()    # {nation, {hhid, HBS_category}}
        expStat[y] = Dict{String, Dict{String, Float64}}()   # {nation, {2 or 3-digit COICOP_code, exp.}}
        corrCds[y] = Dict{String, Array{String, 1}}()        # HBS code correspoding COICOP code: {nation, {COICOP_code; heCodes, heSubst order}}
        expSum[y] = Dict{String, Dict{String, Float64}}()    # HBS expenditure totals corresponding COICOP
        expSumSc[y] = Dict{String, Dict{String, Float64}}()  # scaled HBS expenditure total for checking
        expQtSc[y] = Dict{String, Dict{String, Float64}}()   # scaled HBS expenditure quota for checking

        for n in filter(x -> haskey(mdata[y], x), nations)
            expStat[y][n] = Dict{String, Float64}()
            expSum[y][n] = Dict{String, Float64}()
            corrCds[y][n] = Array{String, 1}()
            expSumSc[y][n] = Dict{String, Float64}()
            expQtSc[y][n] = Dict{String, Float64}()
        end
    end

    # reading Eurostat COICOP expenditure national accounts
    mkpath(rsplit(statFile, '/', limit = 2)[1])
    f = open(statFile)
    title = strip.(string.(split(readline(f), '\t')))
    yridx = Dict(years .=> [findfirst(x->x==string(y), title) for y in years])
    if percap; dt = "CP_EUR_HAB" else dt = "CP_MEUR" end
    if percap; unit = 1 else unit = 10^6 end
    for l in eachline(f)
        s = strip.(string.(split(l, '\t')))
        typ, cod, nat = strip.(string.(split(s[1], ',')))[[2,3,4]]
        if typ == dt && nat in nations && (cod[1:2] == "CP" || cod == "TOTAL")
            for y in years
                expval = tryparse(Float64, strip(replace(s[yridx[y]], ['b','d','e']=>"")))
                if expval !== nothing; expStat[y][nat][cod] = expval * unit end
            end
        end
    end
    close(f)

    for y in years
        ne = nt = length(heCodes[y])
        if subst
            ns = length(heSubst[y])
            nt += ns
        end
        nuc = 0
        nats = filter(x -> haskey(mdata[y], x), nations)

        # correspoding codes matching
        cpc = Array{String, 1}()    # HBS HE code list matching COICOP code list
        for c in heCodes[y]; push!(cpc, crrHeCp[y][c]) end
        if subst; for c in heSubst[y]; push!(cpc, crrHeCp[y][c]) end end

        for n in nats
            for c in cpc
                uc = c[1:end-1]
                if haskey(expStat[y][n], c)
                    if !cdrepl; push!(corrCds[y][n], c)
                    elseif cdrepl
                        if expStat[y][n][c] > 0; push!(corrCds[y][n], c)
                        elseif haskey(expStat[y][n], uc) && expStat[y][n][uc]>0 ; push!(corrCds[y][n], uc)
                        else push!(corrCds[y][n], "NA"); println("Core matching error: ", n, ",", c, ", ", expStat[y][n][c])
                        end
                    end
                elseif haskey(expStat[y][n], uc) && (!cdrepl || expStat[y][n][uc]>0); push!(corrCds[y][n], uc)
                else push!(corrCds[y][n], "NA"); println("Core matching error: ", n, ",", c, ", ", expStat[y][n][c])
                end
            end
        end

        # for keeping upper COICOP residuals
        if !dist_mode
            crrCpHe = Dict(v => k for (k, v) in crrHeCp[y])
            uppCodes = [crrCpHe[cp] for cp in unique([corrCds[y][n][i][1:4] for n in nats, i=1:nt])]
            uppCodes_flt = sort(filter(x -> !(x in heSubst[y]), uppCodes))
            nuc = length(uppCodes_flt)
        end

        # scaling expenditures
        for n in nats
            hhlist, etable = hhsList[y][n], expTable[y][n]
            cplist = sort(unique(corrCds[y][n]))    # HBS corresponded COICOP code list
            cpidx = Dict{String, Int}()             # COICOP code index: {COICOP_code, index in 'cplist'}
            allcplist = sort(filter(x -> x != "TOTAL", collect(keys(expStat[y][n]))))   # All COICOP code list
            allcpidx = Dict{String, Int}()          # All COICOP code index: {COICOP_code, index in 'allcplist'}

            nh, nc, nac = length(hhlist), length(cplist), length(allcplist)
            cpval = zeros(Float64, nc)              # COICOP expenditures
            hesum = zeros(Float64, nc)              # COICOP corresponding HBS expenditure summation
            scexp = zeros(Float64, nh, nt + nuc)    # scaled expenditure table
            schesum = zeros(Float64, nac)           # COICOP corresponding scaled HBS expenditure summation, for checking
            nsamp = 0                               # total sample household members

            if percap
                hhs = mdata[y][n]
                for h in hhlist
                    if eqsize=="none"; nsamp += hhs[h].size
                    elseif eqsize=="eq"; nsamp += hhs[h].eqsize
                    elseif eqsize=="eqmod"; nsamp += hhs[h].eqmodsize
                    end
                end
            end

            for i=1:nc; cpidx[cplist[i]] = i end
            for i=1:nac; allcpidx[allcplist[i]] = i end
            for i=1:nc; cpval[i] = expStat[y][n][cplist[i]] end
            for i=1:nt; hesum[cpidx[corrCds[y][n][i]]] += sum(etable[:,i]) end
            if percap; hesum ./= nsamp end

            cpsumlist = zeros(Float64, nac)         # All COICOP sectors' assigned total COICOP expenditures
            subhesum = Dict{String, Float64}()      #{CP code (2-digit), total corresponding HE sum}
            subhesumHHs = Dict{String, Array{Float64, 1}}()    # Upper COICOP sector's corresponding HBS_total by household:{CP_code(2-digit),{hhid}}

            # allocate COICOP upper-sectors' expenditure quota
            for i=1:nc
                cp = cplist[i]
                if length(cp)==5 && cpval[i]>0
                    if hesum[i]>0 || altCp[y][cp]=="RT"; cpsumlist[allcpidx[cp]] += cpval[i]
                    elseif hesum[i]==0 && length(altCp[y][cp])==5; cpsumlist[allcpidx[altCp[y][cp]]] += cpval[i]
                    end
                end
            end
            for i=1:nac
                cp = allcplist[i]
                if length(cp)==4; cpsumlist[i] += expStat[y][n][cp]
                elseif length(cp)==5; cpsumlist[allcpidx[cp[1:4]]] -= cpsumlist[i]
                end
            end
            for i=1:nac; expQtSc[y][n][allcplist[i]] = cpsumlist[i] end

            # calculate total COICOP sub-sectors' corresponding HBS expenditures by COICOP upper-sector
            for i=1:nt
                upcp = corrCds[y][n][i][1:4]
                if !haskey(subhesum, upcp)
                    subhesum[upcp] = 0
                    subhesumHHs[upcp] = zeros(Float64, nh)
                end
                subhesum[upcp] += sum(etable[:,i])
                subhesumHHs[upcp][:] += etable[:,i]
            end
            if percap
                for c in collect(keys(subhesum))
                    subhesum[c] /= nsamp
                    for j=1:nh
                        hsize = 0
                        if eqsize=="none"; hsize = mdata[y][n][hhlist[j]].size
                        elseif eqsize=="eq"; hsize = mdata[y][n][hhlist[j]].eqsize
                        elseif eqsize=="eqmod"; hsize = mdata[y][n][hhlist[j]].eqmodsize
                        end
                        subhesumHHs[c][j] /= hsize
                    end
                end
            end

            # scale HBS expenditures to match with corresponding COICOP sub-sector
            # distmode = "sqrt"
            distmode = "ln"
            nDisSamp = 0
            if alter
                hhs = mdata[y][n]
                if distmode == "sqrt"; for i=1:nh; nDisSamp += sqrt(hhs[hhlist[i]].size) end
                elseif distmode == "ln"; for i=1:nh; nDisSamp += log(hhs[hhlist[i]].size) + 1 end
                end
            end

            for i=1:nt
                cp = corrCds[y][n][i]
                cpid = cpidx[cp]
                if cp != "NA"
                    if hesum[cpid]>0
                        # scaling ratio = COICOP_expenditures/HBS_expenditures
                        if length(cp)==5
                            scrat = cpsumlist[allcpidx[cp]] / hesum[cpidx[cp]]
                            scexp[:,i] = etable[:,i] .* scrat
                        end
                    elseif !alter || altCp[y][cp]!="RT"; scexp[:,i] = etable[:,i]
                    end
                end
            end

            if alter
                # distribute a sub-sector's expenditrue to corresponding HBS sectors: if alternative code is "RT"
                rtcnt = Dict{String, Int}()
                for i=1:ne
                    cp = corrCds[y][n][i]
                    if altCp[y][cp]=="RT"; if !haskey(rtcnt, cp); rtcnt[cp] = 1 else rtcnt[cp] += 1 end end
                end
                rtrat = [if altCp[y][corrCds[y][n][i]]=="RT"; 1/rtcnt[corrCds[y][n][i]] else 0 end for i=1:ne]

                for i=1:ne
                    cp = corrCds[y][n][i]
                    cpid = cpidx[cp]
                    if cp != "NA"
                        if hesum[cpid] == 0 && altCp[y][cp]=="RT"
                            hhs = mdata[y][n]
                            scrat = rtrat[i] * cpsumlist[allcpidx[cp]] * nsamp / nDisSamp
                            if distmode == "sqrt"; for j=1:nh; scexp[j,i] += scrat * sqrt(hhs[hhlist[j]].size) end
                            elseif distmode == "ln"; for j=1:nh; scexp[j,i] += scrat * (log(hhs[hhlist[j]].size)+1) end
                            end
                        end
                    end
                end
            end

            if dist_mode
                # distribute upper-sector's expenditure to sub-sectors
                hecnt = zeros(Int, nac)
                for i=1:ne; hecnt[allcpidx[corrCds[y][n][i][1:4]]] += 1 end
                dsrat = [1/hecnt[allcpidx[corrCds[y][n][i][1:4]]] for i=1:ne]
                for i=1:ne
                    cp = corrCds[y][n][i][1:4]
                    scrat = cpsumlist[allcpidx[cp]] / subhesum[cp]
                    for j=1:nh
                        hsize = 0
                        if eqsize=="none"; hsize = mdata[y][n][hhlist[j]].size
                        elseif eqsize=="eq"; hsize = mdata[y][n][hhlist[j]].eqsize
                        elseif eqsize=="eqmod"; hsize = mdata[y][n][hhlist[j]].eqmodsize
                        end
                        scexp[j,i] += dsrat[i] * subhesumHHs[cp][j] * scrat * hsize
                    end
                end
            else
                # keep upper-sector's expenditure
                for uc in uppCodes
                    cp = crrHeCp[y][uc]
                    heidx = uc in heSubst[y] ? ne + findfirst(x -> heSubst[y][x] == uc, 1:ns) : nt + findfirst(x -> uppCodes_flt[x] == uc, 1:nuc)
                    scrat = cpsumlist[allcpidx[cp]] / subhesum[cp]
                    for j=1:nh
                        hsize = 0
                        if eqsize=="none"; hsize = mdata[y][n][hhlist[j]].size
                        elseif eqsize=="eq"; hsize = mdata[y][n][hhlist[j]].eqsize
                        elseif eqsize=="eqmod"; hsize = mdata[y][n][hhlist[j]].eqmodsize
                        end
                        scexp[j,heidx] += subhesumHHs[cp][j] * scrat * hsize
                    end
                end
            end

            for i=1:nc; expSum[y][n][cplist[i]] = hesum[i] end

            # for checking
            for i=1:nt; schesum[allcpidx[corrCds[y][n][i]]] += sum(scexp[:,i]) end
            if percap; schesum ./= nsamp end
            for i=1:nac; expSumSc[y][n][allcplist[i]] = schesum[i] end

            expTableSc[y][n] = scexp
        end

        if !dist_mode; append!(heSubst[y], uppCodes_flt) end
    end

    return expTableSc
end

function printExpStats(year, outputFile; scaled=false)

    global nations, cpCodes, expStat, expSum, expSumSc, expQtSc
    nats = filter(x -> haskey(mdata[year], x), nations)

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    print(f, "Year,Nation,COICOP_Code,Eurostat(national),HBS")
    if scaled; print(f, ",Scaled_HBS") end
    println(f)
    for n in nats
        for c in cpCodes[year]
            print(f, year,",",n,",",c,",")
            if haskey(expStat[year][n], c); print(f, expStat[year][n][c],",") else print(f, "NA,") end
            if haskey(expSum[year][n], c); print(f, expSum[year][n][c],",") else print(f, "NA,") end
            if scaled; if haskey(expSumSc[year][n], c); print(f, expSumSc[year][n][c],",") else print(f, "NA,") end end
            println(f)
        end
    end
    close(f)

    f = open(replace(outputFile, ".csv"=>"_quota.csv"), "w")
    print(f, "Year,Nation,COICOP_Code,Eurostat(COICOP),HBS")
    if scaled; print(f, ",Scaled_HBS") end
    println(f)
    for n in nats
        sume = 0
        sumq = 0
        for c in cpCodes[year]
            if haskey(expSum[year][n], c); sume += expSum[year][n][c] end
            if haskey(expQtSc[year][n], c); sumq += expQtSc[year][n][c] end
        end
        println(f, year,",",n,",Total,",expStat[year][n]["TOTAL"],",",sume,",",sumq)
        for c in cpCodes[year]
            print(f, year,",",n,",",c,",")
            if haskey(expStat[year][n], c); print(f, expStat[year][n][c],",") else print(f, "NA,") end
            if haskey(expSum[year][n], c); print(f, expSum[year][n][c],",") else print(f, "NA,") end
            if scaled; if haskey(expQtSc[year][n], c); print(f, expQtSc[year][n][c],",") else print(f, "NA,") end end
            println(f)
        end
    end
    close(f)
end

function printExpTable(year, outputFile; scaled=false, subst = false)

    global nations, hhsList, heCodes, heSubst, expTable, expTableSc

    if scaled; etab = expTableSc; else etab = expTable end
    nt = length(heCodes[year]); if subst; nt += length(heSubst[year]) end

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    print(f, "Year,Nation,Househod")
    for hc in heCodes[year]; print(f, ",",hc) end
    if subst; for sc in heSubst[year]; print(f, ",",sc) end end
    println(f)
    for n in filter(x -> haskey(mdata[year], x), nations)
        nh = length(hhsList[year][n])
        for i = 1:nh
            print(f, year, ",", n, ",", hhsList[year][n][i])
            for j = 1:nt; print(f, ",", etab[year][n][i,j]) end
            println(f)
        end
    end
    close(f)
end

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

function readCategory(years, inputFile; depth=4, catFile="", inclAbr=false, coicop=false)

    global hhCodes, hmCodes, nationNames, year_list
    global heCats, heCodes, heDescs, heCdHrr, heSubHrr, cpCodes, crrHeCp, altCp
    if isa(years, Int); years = [years] end
    append!(year_list, filter(y-> !(y in year_list), years))

    for y in years
        codes, descs = [], []
        hhCodes[y], hmCodes[y] = Array{String, 1}(), Array{String, 1}()
        heCats[y], heCodes[y], heDescs[y] = Dict{String, String}(), Array{String, 1}(), Array{String, 1}()
        heCdHrr[y], heSubHrr[y] = Array{Dict{String, String}, 1}(), Array{Dict{String, Array{String, 1}}, 1}()
        cpCodes[y], crrHeCp[y], altCp[y] = Array{String, 1}(), Dict{String, String}(), Dict{String, String}()

        xf = XLSX.readxlsx(inputFile)
        sh = xf["Nation"]
        for r in XLSX.eachrow(sh) if XLSX.row_number(r)>1 && !ismissing(r[1]) && !ismissing(r[2]); nationNames[string(r[1])] = string(r[2]) end end
        sh = xf["Consumption_categories_" * string(y)]
        for r in XLSX.eachrow(sh)
            if XLSX.row_number(r)>1 && !ismissing(r[1]) && !ismissing(r[2])
                push!(codes, string(r[1]))
                push!(descs, string(r[2]))
                heCats[y][codes[end]] = descs[end]
                if coicop
                    if !ismissing(r[3]); crrHeCp[y][codes[end]] = string(r[3])
                    elseif codes[end][5:6]=="HE"; println("No corrsponding COICOP code.")
                    end
                end
            end
        end
        cpCodes[y] = sort(filter(x->x!="TOTAL",unique(collect(values(crrHeCp[y])))))
        sh = xf["HH_code_" * string(y)]
        for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1 && !ismissing(r[1]); push!(hhCodes[y], string(r[1])) end end
        sh = xf["HM_code_" * string(y)]
        for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1 && !ismissing(r[1]); push!(hmCodes[y], string(r[1])) end end
        if coicop
            sh = xf["COICOP"]
            for r in XLSX.eachrow(sh)
                if XLSX.row_number(r)>1 && !ismissing(r[1]) && !ismissing(r[4]); altCp[y][string(r[1])] = string(r[4]) end
            end
        end
        close(xf)

        predpt = 0
        for i = 1:length(codes)
            curdpt = length(codes[i]) - 7
            if curdpt == depth && codes[i] != "EUR_HE00" && codes[i] != "EUR_HJ00" && (inclAbr || codes[i][5:6]=="HE")
                push!(heCodes[y], codes[i])
                push!(heDescs[y], descs[i])
            elseif predpt < depth && curdpt <= predpt && codes[i] != codes[i-1] && codes[i-1] != "EUR_HE00" && codes[i-1] != "EUR_HJ00" && (inclAbr || codes[i-1][5:6]=="HE")
                push!(heCodes[y], codes[i-1])
                push!(heDescs[y], descs[i-1])
            end
            predpt = curdpt
        end
        if predpt < depth && (inclAbr || codes[end][5:6]=="HE"); push!(heCodes[y], codes[end]); push!(heDescs[y], descs[end]) end

        for i =1:depth-1
            sub_upp = Dict{String, String}()
            upp_subs = Dict{String, Array{String, 1}}()
            uppcat = ""
            upplen = i + 7
            sublen = i + 8
            for c in codes
                if length(c) == upplen
                    uppcat = c
                    upp_subs[uppcat] = Array{String, 1}()
                elseif length(c) == sublen && c[1:end-1] == uppcat
                    sub_upp[c] = uppcat
                    push!(upp_subs[uppcat], c)
                end
            end
            push!(heCdHrr[y], sub_upp)
            push!(heSubHrr[y], upp_subs)
        end
    end
end

function readHouseholdData(year, mdataPath; visible=false, substitute=false, per=0.01, pea=0)
    # per, pea = permissible error_rate, absolute_error

    global heCodes, hhCodes, hhsList, nations, heCdHrr, heSubHrr, heSubst, substCodes, heRplCd
    global mdata[year] = Dict{String, Dict{String, household}}()
    global hhsList[year] = Dict{String, Array{String, 1}}()

    heSubst[year], heRplCd[year] = Array{String, 1}(), Dict{String, Array{String, 1}}()
    if substitute; substCodes[year] = Dict{String, Dict{String, String}}() end

    nullVal = ["missing", "NA"]
    nd = length(heCdHrr[year])

    files = []
    if isa(mdataPath, AbstractArray)
        hbs_ext = Dict(2010=>"_HBS_hh.xlsx", 2015=>"_MFR_hh.xlsx")
        files = [mp * hbs_ext[year] for mp in mdataPath]
    elseif isa(mdataPath, AbstractString)
        for f in readdir(mdataPath); if endswith(f, "_HBS_hh.xlsx") || endswith(f, "_MFR_hh.xlsx"); push!(files, f) end end
    end
    cnt = 0
    if substitute
        hrrcds, upp_cds, sub_cds = Array{String, 1}(), Array{Array{String, 1}, 1}(), Array{Array{String, 1}, 1}()
        for hrr in heCdHrr[year]
            hrr_cd = sort(collect(keys(hrr)))
            append!(hrrcds, hrr_cd)
            push!(sub_cds, hrr_cd)
            push!(upp_cds, sort(unique(collect(values(hrr)))))
        end
    end

    for f in files
        nc = 0
        hhidx = []
        expidx = []
        hhdata = Dict{String, household}()
        hhs = []
        if visible; cnt += 1; print("  ", cnt,"/",length(files),": ",f) end
        if substitute; subst, hrridx, hrrsum = Dict{String, String}(), Dict{String, Int}(), Dict{String, Float64}() end

        if isa(mdataPath, AbstractArray); nation = f[end-13:end-12]; xf = XLSX.readxlsx(f)
        elseif isa(mdataPath, AbstractString); nation = f[1:2]; xf = XLSX.readxlsx(joinpath(mdataPath, f))
        end
        sh = xf[XLSX.sheetnames(xf)[1]]
        for r in XLSX.eachrow(sh)
            if XLSX.row_number(r) == 1
                nc = XLSX.column_bounds(r)[2]
                sectors = [string(r[i]) for i=1:nc]
                hhidx = [findfirst(x->x==hc, sectors) for hc in hhCodes[year]]
                expidx = [findfirst(x->x==he, sectors) for he in heCodes[year]]
                if hhidx[1] == nothing; hhidx[1] = findfirst(x->x=="new_"*hhCodes[year][1], sectors) end
                if substitute && length(heCdHrr[year])>0
                    for c in hrrcds; hrridx[c] = findfirst(x->x==c, sectors) end
                    for c in collect(values(heCdHrr[year][1])); hrridx[c] = findfirst(x->x==c, sectors) end
                elseif substitute && length(heCdHrr[year])==0; println("Error: Hierarchy code list is empty.")
                end
            elseif !ismissing(r[1])
                if nc != XLSX.column_bounds(r)[2]; println(year, " ", nation, " data's column size doesn't match with: ", nc) end
                d = [string(r[i]) for i=1:nc]
                if nation == "MT"; d = [if occursin("_Inf", x); replace(x,"_Inf"=>""); else x end for x in d] end  # for Malta micro-data

                if d[hhidx[2]] != nation; println(year, " ", nation, " data doesn't match with: ", XLSX.row_number(r)) end
                if d[hhidx[3]] != string(year); println(year, " ", nation, " data's year doesn't match with: ", XLSX.row_number(r)) end
                hh = household(d[hhidx[1]],d[hhidx[2]])
                push!(hhs, d[hhidx[1]])
                if d[hhidx[7]] in nullVal; inc = 0; else inc = parse(Float64, d[hhidx[7]]) end
                if d[hhidx[8]] in nullVal; domexp = 0; else domexp = parse(Float64, d[hhidx[8]]) end
                if d[hhidx[9]] in nullVal; abrexp = 0; else abrexp = parse(Float64, d[hhidx[9]]) end

                hh.nuts1, hh.size, hh.weight, hh.income = d[hhidx[4]], parse(Int16, replace(d[hhidx[5]], "+"=>"")), parse(Float64, d[hhidx[6]]), inc
                hh.domexp, hh.abrexp, hh.totexp = domexp, abrexp, (domexp+abrexp)
                hh.popdens, hh.eqsize, hh.eqmodsize= parse(Int8, d[hhidx[10]]), parse(Float64, d[hhidx[11]]), parse(Float64, d[hhidx[12]])
                hh.incomes = [if d[hhidx[i]] in nullVal; 0; else parse(Float64, d[hhidx[i]]) end for i=13:16]
                hh.source, hh.hhtype1, hh.hhtype2 = parse(Int8, d[hhidx[17]]), parse(Int16, d[hhidx[18]]), parse(Int16, d[hhidx[19]])
                hh.ageprof = [if d[hhidx[i]] in nullVal; 0; else parse(Int, replace(d[hhidx[i]], "+"=>"")) end for i=20:26]
                if d[hhidx[27]] in nullVal; wrk = -1; else wrk = parse(Int16, replace(d[hhidx[27]], "+"=>"")) end
                if d[hhidx[28]] in nullVal; nwrk = -1; else nwrk = parse(Int16, replace(d[hhidx[28]], "+"=>"")) end
                if d[hhidx[29]] in nullVal; act = -1; else act = parse(Int16, replace(d[hhidx[29]], "+"=>"")) end
                hh.working, hh.notworking, hh.activating, hh.occupation = wrk, nwrk, act, d[hhidx[30]]

                hh.expends = [if he == nothing; 0; elseif d[he] in nullVal; 0; else parse(Float64, d[he]) end for he in expidx]
                hh.members = []

                hhdata[d[hhidx[1]]] = hh

                if substitute
                    for i = 1:nd, c in upp_cds[i]   # c = higher-category
                        hg_val = !(d[hrridx[c]] in nullVal) ? parse(Float64, d[hrridx[c]]) : 0    # higher-category's value
                        sb_cds = heSubHrr[year][i][c]   # corresponding sub-categories
                        sb_sum = sum([!(d[hrridx[sc]] in nullVal) ? parse(Float64, d[hrridx[sc]]) : 0 for sc in sb_cds])
                        diff = hg_val - sb_sum
                        if diff < -pea && (hg_val == 0 ? abs(diff)/sb_sum : abs(diff)/hg_val) > per
                            if haskey(hh.substCds, c)
                                sb_c = hh.substCds[c]
                                sb_exp = hh.substExp[sb_c] + diff
                                if sb_exp > 0; hh.substExp[sb_c] = sb_exp
                                else
                                    hh.substCds = filter(x -> last(x) != sb_c, hh.substCds)
                                    hh.substExp = filter(x -> first(x) != sb_c, hh.substExp)
                                end
                            end
                            # println("Upper-sector $c ($hg_val) is less than sub-sectors ($sb_sum)")
                        elseif diff > pea && diff / hg_val > per
                            sbst = filter(x->d[hrridx[x]] in nullVal || parse(Float64, d[hrridx[x]]) == 0, sb_cds)
                            if length(sbst) > 0; for sc in sbst; hh.substCds[sc] = c end
                            elseif length(sbst) == 0; hh.substCds[c] = c
                            end
                            hh.substExp[c] = diff
                        elseif hg_val == 0 && sb_sum == 0 && haskey(hh.substCds, c)
                            for sc in sb_cds; hh.substCds[sc] = hh.substCds[c] end
                        end
                    end
                    hh.substExp = filter(x -> last(x) > 0, hh.substExp)
                    hh.substCds = filter(x -> haskey(hh.substExp, last(x)), hh.substCds)
                end
            end
        end
        close(xf)

        if substitute
            for hh in collect(values(hhdata)), hc in filter(x->haskey(hh.substCds, x), heCodes[year])
                if !haskey(subst, hc) || length(subst[hc]) < length(hh.substCds[hc])
                    subst[hc] = hh.substCds[hc]
                end
            end
            substCodes[year][nation] = subst
            append!(heSubst[year], sort(unique(collect(values(subst)))))
        end

        mdata[year][nation] = hhdata
        hhsList[year][nation] = hhs

        if visible; println(", ", length(hhdata), " households") end
    end
    nations = sort(collect(keys(mdata[year])))

    if substitute
        sort!(unique!(heSubst[year]))
        for sc in heSubst[year]
            heRplCd[year][sc] = Array{String, 1}()
            for n in filter(x -> haskey(substCodes[year], x), nations)
                append!(heRplCd[year][sc], filter(x -> substCodes[year][n][x] == sc, collect(keys(substCodes[year][n]))))
            end
            sort!(heRplCd[year][sc])
        end
    end
end

# function readHouseholdData(year, mdataPath; visible=false, substitute=false, per=0.01, pea=0)
#     # per, pea = permissible error_rate, absolute_error
#
#     global heCodes, hhCodes, hhsList, nations, heCdHrr, heSubHrr, heSubst, substCodes, heRplCd
#     global mdata[year] = Dict{String, Dict{String, household}}()
#     global hhsList[year] = Dict{String, Array{String, 1}}()
#
#     heSubst[year], heRplCd[year] = Array{String, 1}(), Dict{String, Array{String, 1}}()
#     if substitute; substCodes[year] = Dict{String, Dict{String, String}}() end
#
#     nullVal = ["missing", "NA"]
#     nd = length(heCdHrr[year])
#
#     files = []
#     if isa(mdataPath, AbstractArray)
#         hbs_ext = Dict(2010=>"_HBS_hh.xlsx", 2015=>"_MFR_hh.xlsx")
#         files = [mp * hbs_ext[year] for mp in mdataPath]
#     elseif isa(mdataPath, AbstractString)
#         for f in readdir(mdataPath); if endswith(f, "_HBS_hh.xlsx") || endswith(f, "_MFR_hh.xlsx"); push!(files, f) end end
#     end
#     cnt = 0
#     if substitute
#         hrrcds, upp_cds, sub_cds = Array{String, 1}(), Array{Array{String, 1}, 1}(), Array{Array{String, 1}, 1}()
#         for hrr in heCdHrr[year]
#             hrr_cd = sort(collect(keys(hrr)))
#             append!(hrrcds, hrr_cd)
#             push!(sub_cds, hrr_cd)
#             push!(upp_cds, sort(unique(collect(values(hrr)))))
#         end
#     end
#
#     for f in files
#         nc = 0
#         hhidx = []
#         expidx = []
#         hhdata = Dict{String, household}()
#         hhs = []
#         if visible; cnt += 1; print("  ", cnt,"/",length(files),": ",f) end
#         if substitute; subst, hrridx, hrrsum = Dict{String, String}(), Dict{String, Int}(), Dict{String, Float64}() end
#
#         if isa(mdataPath, AbstractArray); nation = f[end-13:end-12]; xf = XLSX.readxlsx(f)
#         elseif isa(mdataPath, AbstractString); nation = f[1:2]; xf = XLSX.readxlsx(joinpath(mdataPath, f))
#         end
#         sh = xf[XLSX.sheetnames(xf)[1]]
#         for r in XLSX.eachrow(sh)
#             if XLSX.row_number(r) == 1
#                 nc = XLSX.column_bounds(r)[2]
#                 sectors = [string(r[i]) for i=1:nc]
#                 hhidx = [findfirst(x->x==hc, sectors) for hc in hhCodes[year]]
#                 expidx = [findfirst(x->x==he, sectors) for he in heCodes[year]]
#                 if hhidx[1] == nothing; hhidx[1] = findfirst(x->x=="new_"*hhCodes[year][1], sectors) end
#                 if substitute && length(heCdHrr[year])>0
#                     for c in hrrcds; hrridx[c] = findfirst(x->x==c, sectors) end
#                     for c in collect(values(heCdHrr[year][1])); hrridx[c] = findfirst(x->x==c, sectors) end
#                 elseif substitute && length(heCdHrr[year])==0; println("Error: Hierarchy code list is empty.")
#                 end
#             elseif !ismissing(r[1])
#                 if nc != XLSX.column_bounds(r)[2]; println(year, " ", nation, " data's column size doesn't match with: ", nc) end
#                 d = [string(r[i]) for i=1:nc]
#                 if nation == "MT"; d = [if occursin("_Inf", x); replace(x,"_Inf"=>""); else x end for x in d] end  # for Malta micro-data
#
#                 if d[hhidx[2]] != nation; println(year, " ", nation, " data doesn't match with: ", XLSX.row_number(r)) end
#                 if d[hhidx[3]] != string(year); println(year, " ", nation, " data's year doesn't match with: ", XLSX.row_number(r)) end
#                 hh = household(d[hhidx[1]],d[hhidx[2]])
#                 push!(hhs, d[hhidx[1]])
#                 if d[hhidx[7]] in nullVal; inc = 0; else inc = parse(Float64, d[hhidx[7]]) end
#                 if d[hhidx[8]] in nullVal; domexp = 0; else domexp = parse(Float64, d[hhidx[8]]) end
#                 if d[hhidx[9]] in nullVal; abrexp = 0; else abrexp = parse(Float64, d[hhidx[9]]) end
#
#                 hh.nuts1, hh.size, hh.weight, hh.income = d[hhidx[4]], parse(Int16, replace(d[hhidx[5]], "+"=>"")), parse(Float64, d[hhidx[6]]), inc
#                 hh.domexp, hh.abrexp, hh.totexp = domexp, abrexp, (domexp+abrexp)
#                 hh.popdens, hh.eqsize, hh.eqmodsize= parse(Int8, d[hhidx[10]]), parse(Float64, d[hhidx[11]]), parse(Float64, d[hhidx[12]])
#                 hh.incomes = [if d[hhidx[i]] in nullVal; 0; else parse(Float64, d[hhidx[i]]) end for i=13:16]
#                 hh.source, hh.hhtype1, hh.hhtype2 = parse(Int8, d[hhidx[17]]), parse(Int16, d[hhidx[18]]), parse(Int16, d[hhidx[19]])
#                 hh.ageprof = [if d[hhidx[i]] in nullVal; 0; else parse(Int, replace(d[hhidx[i]], "+"=>"")) end for i=20:26]
#                 if d[hhidx[27]] in nullVal; wrk = -1; else wrk = parse(Int16, replace(d[hhidx[27]], "+"=>"")) end
#                 if d[hhidx[28]] in nullVal; nwrk = -1; else nwrk = parse(Int16, replace(d[hhidx[28]], "+"=>"")) end
#                 if d[hhidx[29]] in nullVal; act = -1; else act = parse(Int16, replace(d[hhidx[29]], "+"=>"")) end
#                 hh.working, hh.notworking, hh.activating, hh.occupation = wrk, nwrk, act, d[hhidx[30]]
#
#                 hh.expends = [if he == nothing; 0; elseif d[he] in nullVal; 0; else parse(Float64, d[he]) end for he in expidx]
#                 hh.members = []
#
#                 hhdata[d[hhidx[1]]] = hh
#
#                 if substitute
#                     for i = 1:nd, c in upp_cds[i]   # c = higher-category
#                         hg_val = !(d[hrridx[c]] in nullVal) ? parse(Float64, d[hrridx[c]]) : 0    # higher-category's value
#                         sb_cds = heSubHrr[year][i][c]   # corresponding sub-categories
#                         sb_sum = sum([!(d[hrridx[sc]] in nullVal) ? parse(Float64, d[hrridx[sc]]) : 0 for sc in sb_cds])
#                         diff = hg_val - sb_sum
#                         if diff < -pea && (hg_val == 0 ? abs(diff)/sb_sum : abs(diff)/hg_val) > per
#                             if haskey(hh.substCds, c)
#                                 sb_c = hh.substCds[c]
#                                 sb_exp = hh.substExp[sb_c] - sb_sum
#                                 if sb_exp > 0; hh.substExp[sb_c] = sb_exp
#                                 else
#                                     hh.substCds = filter(x -> last(x) != sb_c, hh.substCds)
#                                     hh.substExp = filter(x -> first(x) != sb_c, hh.substExp)
#                                 end
#                             end
#                             # println("Upper-sector $c ($hg_val) is less than sub-sectors ($sb_sum)")
#                         elseif diff > pea && diff / hg_val > per
#                             sbst = filter(x->d[hrridx[x]] in nullVal || parse(Float64, d[hrridx[x]]) == 0, sb_cds)
#                             if length(sbst) > 0
#                                 for sc in sbst
#                                     hh.substCds[sc] = c
#                                     # if !haskey(subst, sc); subst[sc] = c end
#                                 end
#                             elseif length(sbst) == 0
#                                 hh.substCds[c] = c
#                                 # subst[c] = c
#                             end
#                             hh.substExp[c] = diff
#                         elseif hg_val == 0 && sb_sum == 0 && haskey(hh.substCds, c)
#                             for sc in sb_cds
#                                 hh.substCds[sc] = hh.substCds[c]
#                                 # subst[sc] = c
#                             end
#                         end
#                     end
#                     hh.substExp = filter(x -> last(x) > 0, hh.substExp)
#                     hh.substCds = filter(x -> haskey(hh.substExp, last(x)), hh.substCds)
#                 end
#             end
#         end
#         close(xf)
#
#         if substitute
#             for hh in collect(values(hhdata)), hc in filter(x->haskey(hh.substCds, x), heCodes[year])
#                 if !haskey(subst, hc) || length(subst[hc]) < length(hh.substCds[hc])
#                     # if haskey(subst, hc) && !(subst[hc] in heSubst[year]); push!(heSubst[year], subst[hc]) end
#                     subst[hc] = hh.substCds[hc]
#                     # if !(subst[hc] in heSubst[year]); push!(heSubst[year], subst[hc]) end
#                     # println("Substitute code of $hc collapese between ", subst[hc]," and ", hh.substCds[hc])
#                 end
#             end
#             substCodes[year][nation] = subst
#             append!(heSubst[year], sort(unique(collect(values(subst)))))
#
#             # # Remove this part if there is no problem
#             # sb_cds = collect(keys(subst))
#             # if length(subst) > 0
#                 # for i = 1:nd, c in filter(x->haskey(heSubHrr[year][i], x), sb_cds)
#                 #     for sc in filter(x->!haskey(subst, x), heSubHrr[year][i][c])
#                 #         subst[sc] = subst[c]
#                 #     end
#                 # end
#                 # substCodes[year][nation] = filter(x->first(x) in heCodes[year], subst)
#             # end
#         end
#
#         mdata[year][nation] = hhdata
#         hhsList[year][nation] = hhs
#
#         if visible; println(", ", length(hhdata), " households") end
#     end
#     nations = sort(collect(keys(mdata[year])))
#
#     if substitute
#         sort!(unique!(heSubst[year]))
#         for sc in heSubst[year]
#             heRplCd[year][sc] = Array{String, 1}()
#             for n in filter(x -> haskey(substCodes[year], x), nations)
#                 append!(heRplCd[year][sc], filter(x -> substCodes[year][n][x] == sc, collect(keys(substCodes[year][n]))))
#             end
#             sort!(heRplCd[year][sc])
#         end
#     end
# end

function readMemberData(year, mdataPath; visible=false)

    global hmCodes, mdata

    nullVal = ["missing", "NA"]

    files = []
    if isa(mdataPath, AbstractArray)
        hbs_ext = Dict(2010=>"_HBS_hm.xlsx", 2015=>"_MFR_hm.xlsx")
        files = [mp * hbs_ext[year] for mp in mdataPath]
    elseif isa(mdataPath, AbstractString)
        for f in readdir(mdataPath); if endswith(f, "_HBS_hm.xlsx")|| endswith(f, "_MFR_hm.xlsx"); push!(files, f) end end
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
                hmidx = [findfirst(x->x==hc, sectors) for hc in hmCodes[year]]
                if hmidx[1] == nothing; hmidx[1] = findfirst(x->x=="new_"*hmCodes[year][1], sectors) end
                if nation == "SI"
                    if hmidx[11] == nothing; hmidx[11] = findfirst(x->x==hmCodes[year][11]*"_Recoded_4Classes", sectors) end
                    if hmidx[12] == nothing; hmidx[12] = findfirst(x->x==hmCodes[year][12]*"_Recoded_4Classes", sectors) end
                    if hmidx[17] == nothing; hmidx[17] = findfirst(x->x==hmCodes[year][17]*"_Recoded", sectors) end
                end
            elseif !ismissing(r[1])
                if nc != XLSX.column_bounds(r)[2]; println(year, " ", nation, " data's column size doesn't match with: ", nc) end
                d = [string(r[i]) for i=1:nc]
                if nation == "MT"; d = [if occursin("_Inf", x); replace(x,"_Inf"=>""); else x end for x in d] end  # for Malta micro-data

                if d[hmidx[2]] != nation; println(year, " ", nation, " data doesn't match with: ", XLSX.row_number(r)) end
                if d[hmidx[3]] != string(year); println(year, " ", nation, " data's year doesn't match with: ", XLSX.row_number(r)) end
                hm = member(d[hmidx[1]], d[hmidx[2]])
                if nation != "MT";  hm.birthNat, hm.citizNat, hm.residNat = parse(Int16, d[hmidx[4]]), parse(Int16, d[hmidx[5]]), parse(Int16, d[hmidx[6]])
                elseif nation == "MT"; hm.birthNat, hm.citizNat, hm.residNat = parse(Int16, split(d[hmidx[4]],"_")[1]), parse(Int16, split(d[hmidx[5]],"_")[1]), parse(Int16, split(d[hmidx[6]],"_")[1])
                end
                hm.gender, hm.mar, hm.union = parse(Int8, d[hmidx[7]]), parse(Int8, d[hmidx[8]]), parse(Int8, d[hmidx[9]])
                hm.relat = d[hmidx[10]] in nullVal ? 0 : parse(Int8, d[hmidx[10]])
                hm.edu, hm.educur, hm.age = parse(Int8, d[hmidx[11]]), parse(Int, d[hmidx[12]]), d[hmidx[13]]
                if nation == "MT" && d[hmidx[14]] == "3_7"; hm.activ = 2; else hm.activ = parse(Int8, d[hmidx[14]]) end
                hm.workhrs, hm.worktyp, hm.worksec = parse(Int8, d[hmidx[15]]), parse(Int8, d[hmidx[16]]), d[hmidx[17]]
                if nation == "MT"; hm.worksts = parse(Int8, split(d[hmidx[18]],"_")[1]); else hm.worksts = parse(Int8, d[hmidx[18]]) end
                hm.occup, hm.income = d[hmidx[19]], (d[hmidx[20]] in [nullVal; "X"] ? 0 : parse(Float64, d[hmidx[20]]))

                if d[hmidx[1]] != "NA"; push!(mdata[year][nation][d[hmidx[1]]].members, hm) end
            end
        end
        if visible; println(", ", nm, " members") end
    end
end

function buildExpenditureMatrix(year; substitute=false)
    global nations, mdata, hhsList, heCodes, heSubst, substCodes
    global expTable[year] = Dict{String, Array{Float64, 2}}()
    ne = length(heCodes[year])
    ns = length(heSubst[year])
    if substitute; nt = ne+ns; else nt = ne end

    for n in filter(x -> haskey(mdata[year], x), nations)
        nh = length(hhsList[year][n])
        hhdata = mdata[year][n]
        hhlist = hhsList[year][n]
        etable = zeros(Float64, length(hhsList[year][n]), nt)
        for i=1:nh; etable[i,1:ne] = mdata[year][n][hhsList[year][n][i]].expends[1:ne] end

        if substitute
            for i=1:nh
                etable[i,ne+1:nt] = [haskey(hhdata[hhlist[i]].substExp, sc) ? hhdata[hhlist[i]].substExp[sc] : 0.0 for sc in heSubst[year]]
            end
        end

        # if substitute && haskey(substCodes[year], n)
        #     sc_list = sort(collect(values(substCodes[year][n])))
        #     sc_idxs = [ne + findfirst(x->x==sc, heSubst[year]) for sc in sc_list]
        #     for i=1:nh; etable[i,sc_idxs] = [haskey(hhdata[hhlist[i]].substExp, sc) ? hhdata[hhlist[i]].substExp[sc] : 0.0 for sc in sc_list] end
        # end

        expTable[year][n] = etable
    end

    return expTable
end

function makeStatistics(year, outputFile; substitute=false)

    global nations, nationNames, hhsList, heCodes, mdata

    domIdx = findall(x->x[5:6]=="HE", heCodes[year])
    abrIdx = findall(x->x[5:6]=="HJ", heCodes[year])

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    print(f,"Year,NC,Nation,Households,Members,Inc_PerCap,Exp_PerCap,Wgh_hhs,Wgh_mm,Wgh_IncPerCap,Wgh_ExpPerCap,ExpPerHH,ExpPerEqSiz")
    println(f, ",Exp_PerHH,Exp_PerCap,ExpTotChk,ExpAbroad_PerCap")
    for n in filter(x -> haskey(mdata[year], x), nations)
        nm = incs = exps = expabr = mexps = wghhhs = wghmms = wghincs = wghexps = wghexpabr = eqsize = expchk = 0
        nh = length(hhsList[year][n])
        for h in hhsList[year][n]
            hh = mdata[year][n][h]
            nm += hh.size
            eqsize += hh.eqsize
            incs += hh.income
            exps += hh.domexp
            expabr += hh.abrexp
            wghhhs += hh.weight
            wghmms += hh.weight * hh.size
            wghincs += hh.weight * hh.income
            wghexps += hh.weight * hh.domexp
            wghexpabr += hh.weight * hh.abrexp
            tmpexp = sum(hh.expends[domIdx])
            if substitute; tmpexp += sum(values(hh.substExp)) end
            mexps += tmpexp         # aggregated household expenditure
            expchk += abs(hh.domexp - tmpexp)
        end
        expPerHHs = exps / nh
        expPerEqSize = exps / eqsize
        incPerCap = incs / nm
        expPerCap = exps / nm
        expAbrPerCap = expabr / nm
        wghincs /= wghmms
        wghexps /= wghmms
        wghexpabr /= wghmms
        expchk /= nh
        mexpPerHHs = mexps / nh
        mexpPerCap = mexps / nm

        print(f, year,",",n,",",nationNames[n],",",nh,",",nm,",",incPerCap,",",expPerCap,",",wghhhs,",",wghmms,",",wghincs,",",wghexps)
        println(f,",",expPerHHs,",",expPerEqSize,",",mexpPerHHs,",",mexpPerCap,",",expchk,",",expAbrPerCap)
    end
    close(f)
end

function readSubstCodesCSV(year, substcdFile, linkFile)

    global hhsList, mdata, substCodes, heSubst, heRplCd

    if !haskey(heSubst, year); heSubst[year] = Array{String, 1}() end
    if !haskey(substCodes, year); substCodes[year] = Dict{String, Dict{String, String}}() end
    if !haskey(heRplCd, year); heRplCd[year] = Dict{String, Array{String, 1}}() end

    f = open(substcdFile)
    readline(f)
    for l in eachline(f)
        sc = strip(string(split(l, ',')[1]))
        if length(sc) > 0; push!(heSubst[year], sc) end
    end
    close(f)

    f = open(linkFile)
    readline(f)
    for l in eachline(f)
        s = string.(split(l, ','))
        if year == parse(Int, s[1])
            if !haskey(substCodes[year], s[2]); substCodes[year][s[2]] = Dict{String, String}() end
            if !haskey(heRplCd[year], s[4]); heRplCd[year][s[4]] = Array{String, 1}() end
            substCodes[year][s[2]][s[3]] = s[4]
            if !(s[3] in heRplCd[year][s[4]]); push!(heRplCd[year][s[4]], s[3]) end
            if !(s[4] in heSubst[year]); println("HBS Subst code list does not contain: ", s[4]) end
        else println("Substitution code file has wrong year: ", s[1])
        end
    end
    close(f)
end

function readPrintedHouseholdData(inputFile)

    global nations, hhsList, mdata

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
            if !(nation in nations); push!(nations, nation) end
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

    global hhsList, mdata

    f = open(inputFile)
    readline(f)
    for l in eachline(f)
        s = string.(split(l, ','))
        hm = member(s[3], s[2])
        hm.birthNat,hm.citizNat,hm.residNat,hm.gender,hm.mar,hm.union,hm.relat = parse(Int16,s[4]),parse(Int16,s[5]),parse(Int16,s[6]),parse(Int8,s[7]),parse(Int8,s[8]),parse(Int8,s[9]),parse(Int8,s[10])
        hm.edu,hm.educur,hm.age,hm.activ,hm.workhrs,hm.worktyp,hm.worksec,hm.worksts = parse(Int8,s[11]),parse(Int,s[12]),s[13],parse(Int8,s[14]),parse(Int8,s[15]),parse(Int8,s[16]),s[17],parse(Int8,s[18])
        hm.occup,hm.income = s[19],parse(Float64,s[20])
        push!(mdata[parse(Int,s[1])][s[2]][s[3]].members, hm)
    end
    close(f)
end

function readPrintedExpenditureData(inputFile; substitute=false, buildHhsExp=false)

    global year_list, mdata, hhsList, heCodes, heSubst, expTable

    nc = [length(heCodes[y]) for y in year_list]
    ns = [length(heSubst[y]) for y in year_list]
    if substitute; nt = nc+ns; else nt = nc end

    f = open(inputFile)
    # check sectors on the title line
    scl = string.(strip.(split(readline(f), ',')[4:end]))
    for y in filter(x -> startswith(inputFile, string(x)), year_list)
        sc_list = substitute ? [heCodes[y];heSubst[y]] : heCodes[y][:]
        if scl != sc_list; println("Sector list is not consistent with expenditure matrix, ", filter(x->!(x in sc_list), scl), ", ", filter(x->!(x in scl), sc_list)) end
        # if scl == sc_list; println("Sector list and matrix columns are consistent.") end
    end
    # read matrix
    for l in eachline(f)
        s = string.(strip.(split(l, ',')))
        year = parse(Int,s[1]); n = s[2]; hh = s[3]
        yi = findfirst(x->x==year, year_list)
        if !haskey(expTable, year); expTable[year] = Dict{String, Array{Float64, 2}}() end
        if !haskey(expTable[year], n); expTable[year][n] = zeros(Float64, length(hhsList[year][n]), nt[yi]) end
        idx = findfirst(x -> x==hh, hhsList[year][n])
        expTable[year][n][idx,1:nt[yi]] = [parse(Float64, x) for x in s[4:nt[yi]+3]]

        if buildHhsExp
            mdata[year][n][hh].expends = [parse(Float64, x) for x in s[4:nc[yi]+3]]
            if substitute
                for sc in heSubst[year]
                    mdata[year][n][hh].substExp[sc] = parse(Float64, s[3+nc[yi]+findfirst(x->x==sc, heSubst[year])])
                end
            end

            # previous version
            #
            # if substitute && haskey(substCodes[year], n)
            #     for sc in collect(values(substCodes[year][n]))
            #         mdata[year][n][hh].substExp[sc] = parse(Float64, s[3+nc[yi]+findfirst(x->x==sc, heSubst[year])])
            #     end
            # end
        end
    end
    close(f)
end

function printCategory(year, outputFile; substitute=false)

    global heCodes, heDescs, heSubst, substCodes

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    println(f,"Code,Description")
    for i = 1:length(heCodes[year]); println(f,heCodes[year][i],",\"",heDescs[year][i],"\"") end
    close(f)
    if substitute
        f = open(replace(replace(outputFile, ".csv"=>"_subst.csv"), ".txt"=>"_subst.txt"), "w")
        println(f,"Code,Description")
        for sc in heSubst[year]; println(f,sc,",\"",heCats[year][sc],"\"") end
        close(f)

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

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    cnt = 0

    print(f, "Year,Nation,HHID,NUTS1,HH_size,Weight")
    print(f, ",Income,Tot_exp,Dom_exp,Abr_exp,Pop_dens,Eq_size,EqMod_size")
    print(f, ",Inc_empl,Inc_nonSal,Inc_rent,Inc_monNet,Inc_source,HHtype1,HHtype2")
    println(f, ",age_0_4,age_5_13,age_14_15,age_16_24,age_16_24_stu,age_25_64,age_65_,working,notworking,activating,occupation")
    for n in filter(x -> haskey(mdata[year], x), nations)
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

    print(" $cnt households")
end

function printMemberData(year, outputFile)

    global nations, hhsList, mdata

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    cnt = 0

    println(f, "Year,Nation,HHID,BirthNat,CitizNat,ResidNat,Gender,Mar,Union,Relat,Edu,EduCur,Age,ActSta,WorkHrs,WorkTyp,WorkSec,WorkSta,Occup,Income")
    for n in filter(x -> haskey(mdata[year], x), nations)
        for h in hhsList[year][n]
            for m in mdata[year][n][h].members
                print(f, year, ",", m.nation, ",", m.hhid, ",", m.birthNat, ",", m.citizNat, ",", m.residNat, ",", m.gender, ",", m.mar, ",", m.union, ",", m.relat)
                print(f, ",", m.edu, ",", m.educur, ",", m.age, ",", m.activ, ",", m.workhrs, ",", m.worktyp, ",", m.worksec, ",", m.worksts, ",", m.occup, ",", m.income)
                println(f)
                cnt += 1
            end
        end
    end
    close(f)
    print(" $cnt members")
end

function printExpenditureMatrix(year, outputFile; substitute=false)

    global nations, hhsList, heCodes, heSubst, expTable
    ne = length(heCodes[year])
    ns = length(heSubst[year])
    if substitute; nt = ne+ns; else nt = ne end

    mkpath(rsplit(outputFile, '/', limit = 2)[1])
    f = open(outputFile, "w")
    print(f, "Year,Nation,Househod")
    for hc in heCodes[year]; print(f, ",",hc) end
    if substitute; for sc in heSubst[year]; print(f, ",",sc) end end
    println(f)
    for n in filter(x -> haskey(mdata[year], x), nations)
        nh = length(hhsList[year][n])
        for i = 1:nh
            print(f, year, ",", n, ",", hhsList[year][n][i])
            for j = 1:nt; print(f, ",", expTable[year][n][i,j]) end
            println(f)
        end
    end
    close(f)
end

function exchangeExpCurrency(exchangeRate; inverse=false, year=0)
    # exchangeRate: can be a file path that contains excahnge rates, a constant value of
    #               EUR to USD currency exchange rate (USD/EUR), or a set of values of Dict[MMYY] or Dict[YY]
    global nations, hhsList, mdata, expTable, expTableSc
    if year > 0; yrs = [year] else yrs = sort(collect(keys(mdata))) end

    # read exchange rate from the recieved file if 'exchangeRate' is 'String'
    if typeof(exchangeRate) <: AbstractString
        erd = Dict{Int, Float64}()
        mkpath(rsplit(exchangeRate, '/', limit = 2)[1])
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
        nats = filter(x -> haskey(hhsList[year], x), nations)
        for n in nats; for h in hhsList[year][n]; mdata[year][n][h].expends .*= er end end
        if length(expTable) > 0; for n in nats; expTable[year][n] .*= er end end
    # if 'exchangeRate' is a set of rates
    elseif typeof(exchangeRate) <: AbstractDict
        for year in yrs
            if haskey(exchangeRate, year); er = exchangeRate[year] else println(year," year exchange rate is not on the list.") end
            nats = filter(x -> haskey(hhsList[year], x), nations)
            for n in nats
                for h in hhsList[year][n]
                    mdata[year][n][h].expends .*= er
                    for sc in collect(keys(mdata[year][n][h].substExp)); mdata[year][n][h].substExp[sc] *= er end
                end
                if haskey(expTable, year) && haskey(expTable[year], n); expTable[year][n] .*= er end
                if haskey(expTableSc, year) haskey(expTableSc[year], n); expTableSc[year][n] .*= er end
            end
        end
    end
end

function convertToPPP(pppRateFile; inverse=false)
    # PPP rates: can be a file path that contains PPP rates, a constant value of
    #            EUR to USD currency exchange rate (USD/EUR), or a set of values of Dict[MMYY] or Dict[YY]
    global nations, hhsList, mdata
    ppps = Dict{Int, Dict{String, Float64}}()

    # read converting rate from the recieved file if 'pppFile' is 'String'
    mkpath(rsplit(pppRateFile, '/', limit = 2)[1])
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
        for n in filter(x -> haskey(mdata[year], x), nations)
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

function mergeNUTS(year)
    global mdata

    for n in collect(keys(mdata[year])), h in collect(keys(mdata[year][n]))
        mdata[year][n][h].nuts1 = mdata[year][n][h].nuts1[1:2]*"0"
    end
end

function checkPopDens(year)
    global mdata, hhsList
    pd_idx = Dict(1 => 1, 2 => 2, 3 => 3, 9 => 4)

    println("Nation,HHs,Dense,Inter,Sparse,No-data")
    for n in collect(keys(hhsList[year]))
        nh = length(hhsList[year][n])
        pd = zeros(Int, 4)
        for h in hhsList[year][n]; pd[pd_idx[mdata[year][n][h].popdens]] += 1 end

        print(n, ",", nh)
        for i=1:4; print(",", pd[i]) end
        println()
    end
end

function initVars(;year = [], nation = [])

    global mdata, hhsList, expTable, expTableSc
    if isa(year, Number); year = [year] end
    if isa(nation, String); nation = [nation] end

    if length(year) == 0
        if length(nation) == 0;
            mdata = Dict{Int, Dict{String, Dict{String, household}}}()
            hhsList = Dict{Int, Dict{String, Array{String, 1}}}()
            expTable = Dict{Int, Dict{String, Array{Float64, 2}}}()
            expTableSc = Dict{Int, Dict{String, Array{Float64, 2}}}()
        else
            for y in collect(keys(mdata)), n in nation
                if haskey(mdata, y); mdata[y][n] = Dict{String, household}() end
                if haskey(hhsList, y); hhsList[y][n] = Array{String, 1}() end
                if haskey(expTable, y); expTable[y][n] = Array{Float64, 2}(undef, 0, 0) end
                if haskey(expTableSc, y); expTableSc[y][n] = Array{Float64, 2}(undef, 0, 0) end
            end
        end
    else
        if length(nation) == 0;
            for y in year
                mdata[y] = Dict{String, Dict{String, household}}()
                hhsList[y] = Dict{String, Array{String, 1}}()
                expTable[y] = Dict{String, Array{Float64, 2}}()
                expTableSc[y] = Dict{String, Array{Float64, 2}}()
            end
        else
            for y in year, n in nation
                if haskey(mdata, y); mdata[y][n] = Dict{String, household}() end
                if haskey(hhsList, y); hhsList[y][n] = Array{String, 1}() end
                if haskey(expTable, y); expTable[y][n] = Array{Float64, 2}(undef, 0, 0) end
                if haskey(expTableSc, y); expTableSc[y][n] = Array{Float64, 2}(undef, 0, 0) end
            end
        end
    end
end

end
