module SamplingError

# Developed date: 26. Mar. 2020
# Last modified date: 2. Apr. 2020
# Subject: Estimate sampling errors
# Description: Proceed Bootstrap method to estimate sampling errors
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using Bootstrap
using StatsPlots
using Statistics

hhid = Array{String, 1}()   # Household ID
sec = Array{String, 1}()    # India products or services sectors
secName = Dict{String, String}()    # sector name dictionary: {sector code, name}
relName = Array{String, 1}()    # religion name list
siz = Dict{String, Int}()       # hhid's family size: {hhid, number of members}
nam = Dict{String, String}()    # districts' name: {district code, district name}

catList = Array{String, 1}()    # category list: sub-category of food sections
disList = Array{String, 1}()    # district list
relList = Array{String, 1}()    # religion list
expList = Array{Float64, 1}()   # expenditure-level sector list

emissionsHHs = Dict{Int16, Array{Float64, 2}}() # categozied emission by household: {year, {hhid, category}}

function migrateData(year, ec)
    global sec, hhid, secName, relName,  = ec.sec, ec.hhid, ec.secName, ec.relName
    global siz, nam = ec.siz, ec.nam
    global emissionsHHs[year] = ec.emissionsHHs[year]
    global catList, disList, relList, expList = ec.catList, ec.disList, ec.relList, ec.incList
end

function proceedExpenditureBootstrap(year, ebiData, intv; perCap=false, disp=false, output="")
    # plot expenditure-level category emission violin charts
    # ebiData: emission data categorized by expenditure-level

    idxExp = ebiData[7]
    nh = length(hhid)
    nc = length(catList)
    ne = length(expList)
    elist = Array{Array{Float64,2},1}() # emission list be expenditure-level category, {Exp_level, {hhid, category}}

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
        for i=1:nh; idx=idxExp[i]; idxcnt[idx] += 1; elist[idx][idxcnt[idx],:] = eh[i,:] end
    end

    # sampling
    iter = 1000
    avgs = zeros(Float64, ne, iter, nc)
    stds = zeros(Float64, ne, iter, nc)
    for i=1:ne
        nss = size(elist[i],1)
        for j=1:iter
            idx = round.(Int, rand(nss)*nss, RoundUp)
            sample = zeros(Float64, nss, nc)
            for k=1:nss; sample[k,:] = elist[i][idx[k],:] end
            avgs[i,j,:] = mean(sample, dims=1)
            stds[i,j,:] = std(sample, dims=1)
        end

    #    print("This:\t",i); for j=1:nc; print("\t",mean(avgs[i,:,j])) end; println()
    #    println("Boot:\t",i); for j=1:nc; println(bootstrap(mean,  elist[i][:,j], BasicSampling(1000))) end
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
    ytitle = "Carbon emission (tCO2/year/capita)"
    for i=1:nc
        p = plot(title=catList[i], xaxis=xtitle, yaxis=ytitle)
        for j=1:ne
            p = violin!([cat[j]],avgs[j,:,i],leg=false)
        #    p = boxplot!([cat[j]],elist[j][:,i],leg=false)
        end
        if disp; display(p) end
        if length(output)>0; savefig(p,replace(output,".png"=>"_"*catList[i]*".png")) end
    end
end

function proceedDistrictBootstrap(year, ebdData; perCap=false, disp=false, output="")
    # plot expenditure-level category emission violin charts
    # ebiData: emission data categorized by expenditure-level

    idxDis = ebdData[8]
    nh = length(hhid)
    nc = length(catList)
    nd = length(disList)
    elist = Array{Array{Float64,2},1}() # emission list be district, {district, {hhid, category}}

    # prepare emission data list
    eh = emissionsHHs[year]
    idxcnt = zeros(Int, nd)
    if perCap
        for i=1:nd
            cnt = 0
            for j=1:nh; if idxDis[j]==i; cnt += siz[hhid[j]] end end
            push!(elist, zeros(Float64, cnt, nc))
        end
        for i=1:nh
            ms=siz[hhid[i]]; idx=idxDis[i]
            for j=1:ms; idxcnt[idx] += 1; elist[idx][idxcnt[idx],:] = eh[i,:]/ms end
        end
    else
        for i=1:nd; push!(elist, zeros(Float64, count(x->x==i, idxDis), nc)) end
        for i=1:nh
            idx=idxDis[i]
            idxcnt[idx] += 1
            elist[idx][idxcnt[idx],:] = eh[i,:]
        end
    end

    # sampling
    iter = 1000
    avgs = zeros(Float64, nd, iter, nc)
    stds = zeros(Float64, nd, iter, nc)
    for i=1:nd
        nss = size(elist[i],1)
        for j=1:iter
            idx = round.(Int, rand(nss)*nss, RoundUp)
            sample = zeros(Float64, nss, nc)
            for k=1:nss; sample[k,:] = elist[i][idx[k],:] end
            avgs[i,j,:] = mean(sample, dims=1)
            stds[i,j,:] = std(sample, dims=1)
        end

    #    print("This:\t",i); for j=1:nc; print("\t",mean(avgs[i,:,j])) end; println()
    #    println("Boot:\t",i); for j=1:nc; println(bootstrap(mean,  elist[i][:,j], BasicSampling(1000))) end
    end

    # plotting
    xtitle = "Groups by expenditure-level"
    ytitle = "Carbon emission (tCo2/year/capita)"
    for i=1:nc
        p = plot(title=catList[i], xaxis=xtitle, yaxis=ytitle)
        for j=1:nd
            p = violin!([nam[disList[j]]],avgs[j,:,i],leg=false)
        #    p = boxplot!([cat[j]],elist[j][:,i],leg=false)
        end
        if disp; display(p) end
        if length(output)>0; savefig(p,replace(output,".png"=>"_"*catList[i]*".png")) end
    end
end

end
