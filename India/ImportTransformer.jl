module ImportTransformer

# Developed date: 12. Nov. 2019
# Last modified date: 22. Nov. 2019
# Subject: Import account transformation
# Description: Transform import nation's accounts to India accounts.
#              Utilize Eora MRIO, concordance tables, and Comtrade micro-data.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

import XLSX

# India data
sec = Array{String, 1}          # India products or services sectors
cons = Array{Float64, 1}        # India total consumptions by sectors
hhid = Array{String, 1}         # Household ID
hhShare = Array{Float64, 2}     # Enpenditure share by households and by sectors

# Import data
nations = Array{String, 1}                      # Nations list
sectors = Dict{String, Array{String, 1}}()      # products or services sectors
origTrade = Dict{String, Array{Float64, 1}}()   # original trade data
concMat = Dict{String, Array{Float64, 2}}()     # concordance matrix
transf = Dict{String, Array{Float64, 1}}()      # transformed trade matrix for Nation
transfHH = Dict{String, Array{Float64, 2}}()    # transformed trade matrix for Households

function analyzeExpenditures(expMat, row, col)  # row: households, col: product or service categories

    global hhid = row
    global sec = col
    global cons = zeros(Float64, length(sec))
    global hhShare = zeros(Float64, length(sec), length(hhid))

    # calculate total expenditures by sector
    for r = 1:length(sec); for c = 1:length(hhid); cons[r] += expMat[c,r] end end

    #calculate each household's share by sector
    for r = 1:length(sec)
        for c = 1:length(hhid)
            if cons[r] > 0; hhShare[r,c] = expMat[c,r] / cons[r] end
        end
    end

    return cons, hhShare
end

function transformHStoIND(nat, hsSec, tdMat, conMat)  # Nations, HS sectors, trade matrix, concordance matrix

    global nations = nat
    global sectors["hs"] = hsSec
    global origTrade
    global concMat["hs"] = conMat
    global transf

    for n = 1:length(nat)
        origTrade[nat[n]] = tdMat[:, n]

        # calculate transformed trade matrix for nation
        tr = zeros(Float64, length(sec))
        mat = zeros(Float64, length(hsSec), length(sec))
        for r = 1:length(hsSec); for c = 1:length(sec); tr[c] += conMat[r,c] * cons[c] * tdMat[r, n] end end

    #    for r = 1:length(hsSec); for c = 1:length(sec); mat[r,c] = conMat[r,c] * cons[c] end end
    #    for c = 1:length(sec); for r = 1:length(hsSec); tr[c] += mat[r,c] * tdMat[r, n] end end
        transf[nat[n]] = tr
    end

    return transf
end

function transformEORAtoIND(nat, eorSec, fdMat, conMat)
                            # Nations, Eora sectors, import final demand matrix, concordance matrix
    global nations = nat
    global sectors = eorSec
    global origTrade = fdMat
    global concMat
    global transf

    for n in nat
        cmat = conMat[n]    # row: eora index, col: india index
        eosec = eorSec[n]
        fd = fdMat[n]

        tr = zeros(Float64, length(sec))
        mat = zeros(Float64, length(eosec), length(sec))
        for r = 1:length(eosec); for c = 1:length(sec); mat[r,c] = cmat[r,c] * cons[c] end end
        for c = 1:length(sec); for r = 1:length(eosec); tr[c] += mat[r,c] * fd[r] end end

        concMat[n] = cmat
        transf[n] = tr
    end

    return transf
end

function calculateHHshare(outputFile = "", saveData = true)

    global transfHH

    # print a file if "outputFile" has a file name
    if length(outputFile) > 0
        f = open(outputFile, "w")
        print(f, "Sectors/HHID")
        for h in hhid; print(f, "\t$h") end
        println(f)
    end

    # converting process
    for n = 1:length(nations)

        # calculate transformed trade matrix for households
        trhh = zeros(Float64, length(sec), length(hhid))
        for r = 1:length(sec)
            for c = 1:length(hhid)
                trhh[r,c] = transf[nations[n]][r] * hhShare[r,c]
            end
        end

        # for memory management: it do not save matrix data if "saveData" is 'false'.
        if saveData; transfHH[nations[n]] = trhh end

        # print transformed trade matrix by households
        if length(outputFile) > 0
            for s = 1:length(sec)
                print(f, nations[n],"_",sec[s])
                for h = 1:length(hhid); print(f, "\t", trhh[s,h]) end
                println(f)
            end
        end
    end

    if length(outputFile) > 0; close(f) end

    return transfHH
end

function integrateEoraHS(matchingFile, nation) # Eora-HS code matching XLSX file, merging nation

    natMatch = Dict{String, String}()   # {HS nation name, Eora nation name}
    secMatch = Dict{String, String}()   # {India sector code, "Product" or "Service"}
    transfInt = Dict{String, Array{Float64, 1}}()   # Integrated transformed import vector
    transfEor
    transfHs

    # read matching condition
    xf = XLSX.readxlsx(matchingFile)
    sh = xf["Nation"]   # Eora-HS nation code matching
    for r in XLSX.eachrow(sh); if XLSX.row_number(r)>1 && (r[1]|r[2])!==missing; natMatch[r[2]] = r[1] end end
    sh = xf[nation]     # merging Eora-based sectors of 'nation'
    for r in XLSX.eachrow(sh); secMatch[r[1]] = r[3] end
    close(xf)

    # integrate transformed vectors
    for n in sort(collect(keys(natMatch)))
        trEora = transfEor[natMatch[n]]
        trHs = transfHs[n]
        trInt = Array{Float64, 1}()

        for c in sort(collect(keys(secMatch)))
            if secMatch[c] == "Product"; push!(trInt, trHs[findfirst(x -> x==c, sec)])
            elseif secMatch[c] == "Service"; push!(trInt, trEora[findfirst(x -> x==c, sec)])
            else println("Sector entitiy error: ", c,"\t", secMatch[c])
            end
        end

        transfInt[natMatch[n]] = trInt
    end

end

function printTransfImport(outputFile = "")
    f = open(outputFile, "w")

    for n in nations; print(f, "\t$n") end
    println(f)
    for s = 1:length(sec)
        print(f, sec[s])
        for n in nations; print(f, "\t", transf[n][s]) end
        println(f)
    end
    close(f)
end

function printTransfImportHH(outputFile = "")
    f = open(outputFile, "w")

    print(f, "Sectors/HHID")
    for h in hhid; print(f, "\t$h") end
    println(f)
    for n in nations
        for s = 1:length(sec)
            print(f, n,"_",sec[s])
            for h = 1:length(hhid); print(f, "\t", transfHH[n][s,h]) end
            println(f)
        end
    end

    close(f)
end

end
