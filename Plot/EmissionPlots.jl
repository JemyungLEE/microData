module EmissionPlots

# Developed date: 26. Mar. 2020
# Last modified date: 8. Apr. 2020
# Subject: Plotting emission charts
# Description: Read emission data and plot violin and box charts
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using Plots
using Bootstrap
using StatsPlots
using Statistics
using LaTeXStrings

hhid = Array{String, 1}()   # Household ID
sec = Array{String, 1}()    # India products or services sectors
secName = Dict{String, String}()    # sector name dictionary: {sector code, name}
relName = Array{String, 1}()    # religion name list

dis = Dict{String, String}()    # hhid's district: {hhid, district code}
nam = Dict{String, String}()    # districts' name: {district code, district name}
siz = Dict{String, Int}()       # hhid's family size: {hhid, number of members}
wgh = Dict{String, Float64}()   # hhid's population weight: {hhid, weight}
mpce = Dict{String, Float64}()   # hhid's income: {hhid, monthly per capita expenditure (mixed reference period)}

sam = Dict{String, Tuple{Int,Int}}()    # sample population and households by districct: {district code, (population, number of households)}
pop = Dict{String, Tuple{Int,Int,Float64}}()    # population by district: {district code, (population, number of households, area(km^2))}

catList = Array{String, 1}()    # category list: sub-category of food sections
disList = Array{String, 1}()    # district list
relList = Array{String, 1}()    # religion list
expList = Array{Float64, 1}()   # expenditure-level sector list

emissionsHHs = Dict{Int16, Array{Float64, 2}}()     # categozied emission by household: {year, {hhid, category}}
emissionsDis = Dict{Int16, Array{Float64, 2}}()     # categozied emission by district: {year, {district, category}}
emissionsRel = Dict{Int16, Array{Float64, 2}}()     # categozied emission by religion: {year, {religion, category}}
emissionsExp = Dict{Int16, Array{Float64, 2}}()     # categozied emission by incomes: {year, {income level, category}}

colPalTen = Array{String, 1}()  # color palette for 10 categories
colPalTwe = Array{String, 1}()  # color palette for 20 categories

function readColorPalette(inputfile; rev=false, tran=false)
    global colPalTen, colPalTwe

    f = open(inputfile)
    readline(f)
    for i=1:10
        p = readline(f)
        push!(colPalTen, strip(replace(replace(p,"<color>"=>""),"</color>"=>"")))
    end
    readline(f); readline(f); readline(f)
    for i=1:20
        p = readline(f)
        push!(colPalTwe, strip(replace(replace(p,"<color>"=>""),"</color>"=>"")))
    end
    readline(f)
    close(f)

    if rev; colPalTen=reverse(colPalTen); colPalTwe=reverse(colPalTwe) end
    if tran; colPalTen=permutedims(colPalTen); colPalTwe=permutedims(colPalTwe) end
end

function migrateData(year, ec)
    global sec, hhid, secName, relName = ec.sec, ec.hhid, ec.secName, ec.relName
    global dis, nam, siz, wgh, mpce, sam, pop = ec.dis, ec.nam, ec.siz, ec.wgh, ec.inc, ec.sam, ec.pop
    global catList, disList, relList, expList = ec.catList, ec.disList, ec.relList, ec.incList

    global emissionsHHs[year] = ec.emissionsHHs[year]
    global emissionsDis[year] = ec.emissionsDis[year]
#    global emissionsRel[year] = ec.emissionsRel[year]
#    global emissionsExp[year] = ec.emissionsInc[year]
end



function plotCfBubbleChart(year, output=""; disp=false, dataoutput="", povline=1.9)
    # Population density: X
    # CF per capita: Y
    # Poverty ratio: color
    # Total CF: size

    global hhid, dis, siz, wgh, mpce, pop
    global disList, catList, emissionsDis
    nc = length(catList)
    nd = length(disList)
    ed = emissionsDis[year]          # emission per capita by district
    et = zeros(Float64, nd, nc)      # total emission by district
    for i=1:nd; et[i,:] = ed[i,:] * pop[disList[i]][1] end

    popds = [pop[x][1]/pop[x][3] for x in disList]      # population density
    povr = zeros(Float64, nd)

    for h in hhid
        if mpce[h]<povline
            idx = findfirst(x->x==dis[h], disList)
            povr[idx] += siz[h]
        end
    end
    for i=1:nd; povr[i] /= sam[disList[i]][1] end

    # plot charts
    pyplot()
    title = "Carbon footprint per capita and Population density"
    xtitle = L"Population density (Persons/km^2)"
    ytitle = L"Carbon footprint per capita (ton CO_2/year/capita)"
#    for i=1:nc
    for i=nc:nc
        p = scatter(popds, ed[:,i], title=catList[i]*" "*title, xlab=xtitle, xscale=:log10, ylab=ytitle, framestyle=:origin, label="")
#        p = scatter(popds, ed[:,i], markersize=et*10, title=catList[i]*" "*title, xaxis=xtitle, yaxis=ytitle, framestyle=:origin, label="")

        if disp; display(p) end
        if length(output)>0; savefig(p,replace(output,".png"=>"_"*catList[i]*".png")) end
    end

end

function plotExpStackedBarChart(input="", output=""; hhsCF="", reverse=false, perCap=false, disp=false)
    #=
    Food
    Electricity
    Gas
    Other energy
    Medical care
    Public transport
    Private transport
    Education
    Consumable goods
    Durable goods
    Other services
    =#

    exp = Array{String,1}()
    val = Array{String,1}()
    emat = []
    cis = []

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
    smat = sum(emat, dims=2)

    # estimate confidence intervals by bootstrap method
    if length(hhsCF)>0
        expval = []
        f = open(hhsCF)
        s = string.(split(readline(f),','))
        idxCF = findfirst(x->x=="Total", s)
        idxExp = findfirst(x->x=="MPCE", s)
        idxSiz = findfirst(x->x=="HH_size", s)
        for l in eachline(f)
            s = string.(split(l,','))
            push!(expval, [parse(Float64,s[idxCF])/parse(Float64,s[idxSiz]),parse(Float64,s[idxExp])])
        end
        close(f)
        ne = length(expList)
        for i=1:ne
            if i==ne; expByLev = filter(x->expList[ne]<=x[2], expval)
            else expByLev = filter(x->expList[i]<=x[2]<expList[i+1], expval)
            end
            bs = bootstrap(mean, expByLev, BasicSampling(1000))
            ci = confint(bs, BasicConfInt(0.95))
            push!(cis, ci[1])
        end
    end

    # plot charts
    title = "CF stacked bars"
    xtitle = L"Groups\,by\,expenditure-level"
    ytitle = L"Carbon\,footprint\,per\,apita\,(ton\,CO_2/year/capita)"
    err = [(x[1]-x[2],x[1],x[3]-x[2]) for x in cis]

    pyplot()
    #fnt = font(12,"Arial")
    println(exp)
    println(val)
    p = plot(xaxis=xtitle, yaxis=ytitle, framestyle=:box, grid=false)
    p = groupedbar!(exp, emat, labels=permutedims(catList), bar_position=:stack, bar_width=0.67, lw=0, legend=:inside, color=colPalTwe, fg_legend=:transparent, yticks = ([0.0,0.5,1.0,1.5]))
#        legend=:inside, titlefont=fnt, tickfont=fnt, legendfont=fnt, color=colPalTwe, fg_legend = :transparent)
#    p = bar!(exp, smat, yerr=err, label = "", bar_width=0.67, lw=0, marker=stroke(1.5,:black), color=nothing)
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
