module ImportTransformer

# Developed date: 12. Nov. 2019
# Last modified date: 15. Nov. 2019
# Subject: Import account transformation
# Description: Transform import nation's accounts to India accounts.
#              Utilize Eora MRIO, concordance tables, and Comtrade micro-data.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

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
    for r = 1:length(sec); for c = 1:length(hhid); hhShare[r,c] = expMat[c,r] / cons[r] end end

    return cons, hhShare
end

function transformHStoIND(nat, hsSec, tdMat, conMat)  # Nations, HS sectors, trade matrix, concordance matrix

    global nations = nat
    global sectors["hs"] = hsSec
    global origTrade
    global concMat["hs"] = conMat
    global transf
    global transfHH

    for n = 1:length(nat)
        origTrade[nat[n]] = tdMat[:, n]

        # calculate transformed trade matrix for nation
        tr = zeros(Float64, length(sec))
        mat = zeros(Float64, length(hsSec), length(sec))
        for r = 1:length(hsSec); for c = 1:length(sec); mat[r,c] = conMat[r,c] * cons[c] end end
        for c = 1:length(sec); for r = 1:length(hsSec); tr[c] += mat[r,c] * tdMat[r, n] end end
        transf[nat[n]] = tr

        # calculate transformed trade matrix for households
        trhh = zeros(Float64, length(sec), length(hhid))
        for r = 1:length(sec); for c = 1:length(hhid); trhh[r,c] = tr[r] * hhShare[r,c] end end
        transfHH[nat[n]] = trhh
    end
end

function transformEORAtoIND(nat, eorSec, impMat, conMat)
                            # Nations, Eora sectors, import final demand matrix, concordance matrix
    global nations = nat
    global sectors
    global origTrade
    global concMat
    global transf
    global transfHH


end

function integrateTransformedImport(nations, sec, transf, transfhh)

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
