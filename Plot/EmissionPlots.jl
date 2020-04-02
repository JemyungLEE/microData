module EmissionPlots

# Developed date: 26. Mar. 2020
# Last modified date: 1. Apr. 2020
# Subject: Plotting emission charts
# Description: Read emission data and plot violin and box charts
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using StatsPlots

hhid = Array{String, 1}()   # Household ID
sec = Array{String, 1}()    # India products or services sectors
secName = Dict{String, String}()    # sector name dictionary: {sector code, name}
relName = Array{String, 1}()    # religion name list
siz = Dict{String, Int}()       # hhid's family size: {hhid, number of members}

catList = Array{String, 1}()    # category list: sub-category of food sections
disList = Array{String, 1}()    # district list
relList = Array{String, 1}()    # religion list
expList = Array{Float64, 1}()   # expenditure-level sector list

emissionsHHs = Dict{Int16, Array{Float64, 2}}() # categozied emission by household: {year, {hhid, category}}
emissionsDis = Dict{Int16, Array{Float64, 2}}()     # categozied emission by district: {year, {district, category}}
emissionsRel = Dict{Int16, Array{Float64, 2}}()     # categozied emission by religion: {year, {religion, category}}
emissionsExp = Dict{Int16, Array{Float64, 2}}()     # categozied emission by incomes: {year, {income level, category}}


function migrateData(year, ec)
    global sec, hhid, secName, relName, siz = ec.sec, ec.hhid, ec.secName, ec.relName, ec.siz
    global catList, disList, relList, expList = ec.catList, ec.disList, ec.relList, ec.incList

    global emissionsHHs[year] = ec.emissionsHHs[year]
    global emissionsDis[year] = ec.emissionsDis[year]
    global emissionsRel[year] = ec.emissionsRel[year]
    global emissionsExp[year] = ec.emissionsInc[year]
end

function plotExpStackedBarChart(input="", output=""; reverse=false, perCap=false, disp=false)

    exp = Array{String,1}()
    val = Array{String,1}()
    emat = []

    # read data
    f = open(input)
    global catList = filter(x->!(x in ["Total","Pop.","HH.","WghPop."]),string.(split(readline(f),','))[3:end])
    nc = length(catList)
    ec = Array{Array{Float64,1},1}()
    for l in eachline(f)
        s = string.(split(l,','))
        push!(exp, s[1]); push!(val, s[2]); push!(ec, parse.(Float64, s[3:nc+2]))
    end
    close(f)
    ne = length(exp)
    emat = zeros(Float64, ne, nc)
    for i=1:ne; emat[i,:] = ec[i][:] end

    println(exp)
    println(cat)

    # plot charts
    title = "CF stacked bars"
    xtitle = "Groups by expenditure-level"
    ytitle = "Carbon emission (tCO2/year/capita)"
    p = plot(title=title, xaxis=xtitle, yaxis=ytitle, framestyle=:origin)
    p = groupedbar!(exp, emat, bar_position=:stack, bar_width=0.67, lw=0, legend=:outertopright)

    if disp; display(p) end
    if length(output)>0; savefig(p,output) end
end

function plotExpCatViolin(year, ebiData, intv; perCap=false, boxplot=false, disp=false, output="")
    # plot expenditure-level category emission violin charts
    # ebiData: emission data categorized by expenditure-level

    idxExp = ebiData[7]
    nh = length(hhid)
    nc = length(catList)
    ne = length(expList)
    elist = Array{Array{Float64,2},1}() # emission list be expenditure-level category

    # prepare emission data list
    eh = emissionsHHs[year]
    idxcnt = zeros(Int, ne)
    if perCap
        for i=1:ne
            cnt = 0
            for j=1:nh; if idxExp[j]==i; cnt += siz[hhid[j]] end end
            push!(elist, zeros(Float64, cnt, nc))
        end
        for i=1:nh
            ms=siz[hhid[i]]; idx=idxExp[i]
            for j=1:ms; idxcnt[idx] += 1; elist[idx][idxcnt[idx],:] = eh[i,:]/ms end
        end
    else
        for i=1:ne; push!(elist, zeros(Float64, count(x->x==i, idxExp), nc)) end
        for i=1:nh; ms=siz[hhid[i]]; idx=idxExp[i]; idxcnt[idx] += 1; elist[idx][idxcnt[idx],:] = eh[i,:] end
    end

    # category
    cat = []
    push!(cat, "Bottom "*string(round(intv[1]*100,digits=1))*"%")
    for i=2:ne-1; push!(cat, string(round(intv[i-1]*100,digits=1))*"-"*string(round(intv[i]*100,digits=1))*"%") end
    push!(cat, "Top "*string(round((intv[ne]-intv[ne-1])*100,digits=1))*"%")
    #=
    push!(cat, string(round(intv[1]*100,digits=1))*"%,Less "*string(round(expList[2],digits=2)))
    for i=2:ne-1; push!(cat, string(round(intv[i-1]*100,digits=1))*"-"*string(round(intv[i]*100,digits=1))*"%,From "*string(round(expList[i],digits=2))*" to less "*string(round(expList[i+1],digits=2))) end
    push!(cat, string(round((intv[ne]-intv[ne-1])*100,digits=1))*"%,"*string(round(expList[ne],digits=2)))
    =#
    # plotting
    xtitle = "Groups by expenditure-level"
    ytitle = "Carbon emission (tCo2/year/capita)"
    for i=1:nc
        p = plot(title=catList[i], xaxis=xtitle, yaxis=ytitle)
        for j=1:ne
            p = violin!([cat[j]],elist[j][:,i],leg=false)
            p = boxplot!([cat[j]],elist[j][:,i],leg=false)
        end
        if disp; display(p) end
        if length(output)>0; savefig(p,replace(output,".png"=>"_"*catList[i]*".png")) end
    end
end

end
