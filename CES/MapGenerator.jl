module MapGenerator

# Developed date: 16. Dec. 2022
# Last modified date: 16. Dec. 2022
# Subject: Generate regional carbon footprint maps
# Description: Read a base map file and generate carbon footprint maps, as GeoJSON files.
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

include("EmissionCategorizer.jl")

using GeoJSON
using .EmissionCategorizer

js = GeoJSON
ec = EmissionCategorizer

cat_list = Array{String, 1}()   # category list

file_names = Dict{String, Dict{String, String}}()       # map filename: {CF type, {category, name}}

base_map = Dict{Int, Dict{String, Dict{}(), 1}}()       # Base map: {year, {nation A3, {base map Dict}}}
cf_maps = Dict{Int, Dict{String, Array{Dict{}(), 1}}}() # categorized CF maps: {year, {nation A3, {map Dict by category}}}

function readBaseMap(year, nation, map_file; remove = true, alter = true)
    # remove: (default)[true] remain only "geometry", "properties", and "type" features
    # alter: (default)[true] change "GIS_ID" to "KEY_CODE" and "GIS_name" to "EN_NAME", and add "fill_carbon"

    global base_map
    if !haskey(base_map, year); base_map[year] = Dict{String, Dict{}(), 1}() end

    ess_ft = ["geometry", "properties", "type"] # essential features
    ess_pr = ["GIS_name", "GIS_ID"]             # essential properties
    bs_map = js.geo2dict(js.read(read(map_file)))

    ft = bs_map["features"]
    nf = length(ft)

    if remove
        for i = 1:nf, rf in filter(x -> !(x in ess_ft), collect(keys(ft[i]))); delete!(ft[i], rf)) end
    end

    if alter
        for i = 1:nf
            fp = ft[i]["properties"]
            for rp in filter(x -> !(x in ess_pr), collect(keys(fp))); delete!(fp, rp) end
            fp["KEY_CODE"] = fp["GIS_ID"]
            fp["EN_NAME"] = fp["GIS_name"]
            fp["fill_carbon"] = ""
            delete!(fp, "GIS_ID")
            delete!(fp, "GIS_name")
        end
    end

    base_map[year][nation] = bs_map
end

function readFileNames(input_file)

    global file_names

    f_sep = getValueSeparator(inputFile)
    f = open(inputFile)

    title = string.(strip.(split(readline(f), f_sep)))
    ti, ci, fi = [findfirst(x -> x == it, title) for it in ["CF_type", "Category", "File_name"]]

    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        if !haskey(file_names, s[1]); file_names[s[1]] = Dict{String, String}() end
        file_names[s[1]][s[2]] = s[[3]]
    end

    close(f)
end

function mappingGeoJSON()

end

function printMapFiles(year, nation, )

    global cat_list, file_names, cf_maps

    nc = length(cat_list)
    maps = cf_maps[year][nation]

    for i = 1:nc
        t, c = string.(rsplit(cat_list[i], "_", limit = 2))
        f = open(file_names[t][c], "w")
        println(f, js.write(js.dict2geo(maps[i])))
        close(f)
    end
end

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ',' elseif fext == "tsv" || fext == "txt"; return '\t' end
end

end
