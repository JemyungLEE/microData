module BoxChart

# Developed date: 26. Mar. 2020
# Last modified date: 27. Mar. 2020
# Subject: Plotting box-charts
# Description: Plot an overlapping box chart reading two CSV files
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

using Plots

category = Array{String, 1}()   # category list
sector = Array{String, 1}()     # sector list
xval = Array{Float64, 2}      # x-axis values: {category, sector}
yval = Array{Float64, 2}      # y-axis values: {category, sector}

function readValues(xData, yData)

    global category, sector, xval, yval, xttl, x

    f = open(xData)
    category = string.(split(readline(f),',')[2:end-3])
    xval = zeros(Float64, length(category),0)
    i = 1
    for l in eachline(f)
        l = string.(split(l,','))
        push!(sector, "\$"*l[1])
        xval = hcat(xval,map(x->parse(Float64,x)/365,l[2:end-3]))
        i += 1
    end
    close(f)
    f = open(yData)
    if category != string.(split(readline(f),',')[2:end-3]); println("Categories do not match.") end
    yval = zeros(Float64, length(category),0)
    i = 1
    for l in eachline(f)
        l = string.(split(l,','))
        if sector[i] != "\$"*l[1]; println("Sector ",i," do not match.") end
        yval = hcat(yval,map(x->parse(Float64,x),l[2:end-3]))
        i += 1
    end
    close(f)
end

function plotBoxChart(xttl="", yttl="" ; disp=false, output="", dash=false)

    global category, sector, xval, yval
    ns = length(sector)

    for i=1:length(category)
        p = plot(title=category[i], xaxis=xttl, yaxis=yttl, size=(500,500))
        plotOrder = sortperm(xval[i,:], rev=true)
        for j=1:ns
            x = xval[i,plotOrder[j]]; y = yval[i,plotOrder[j]]
            p = plot!(Shape([0,x,x,0],[0,0,y,y]),color=:white,legend=nothing,foreground_color_legend = nothing,label=sector[plotOrder[j]])
            p = annotate!(x, y, text(sector[plotOrder[j]], :black, :right, :top, 8))
        end
        if dash; p = plot!([0,xval[i,plotOrder[1]]],[0,yval[i,plotOrder[1]]], line = (:dot, 1, :black)) end
        if disp; display(p) end
        if length(output)>0; savefig(p,replace(output,".png"=>"_"*category[i]*".png")) end
    end
end

end
