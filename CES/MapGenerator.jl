module MapGenerator

# Developed date: 16. Dec. 2022
# Last modified date: 18. Jan. 2023
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
map_list = Array{String, 1}()   # map category list: ex) ["overall_total", "percap_total", "percap_food", ...]
reg_list = Dict{Int, Dict{String, Array{String, 1}}}()  # Region list: {year, {nation, {gis_id}}}
reg_id = Dict{Int, Dict{String, Dict{String,String}}}() # GIS region ID: {year, {nation, {gis_id, gis_label}}}

pc_reg_cf = Dict{Int, Dict{String, Array{Int, 2}}}()    # per capita emission rank by region: {year, {nation, {region, category}}}
ov_reg_cf = Dict{Int, Dict{String, Array{Int, 2}}}()    # overall emission rank by region: {year, {nation, {region, category}}}

file_names = Dict{String, String}()             # map filename: {'CF type'_'category', name}
field_label = Dict{String, String}()            # map field label: {'CF type'_'category', label}
hex_codes = Dict{String, Array{String, 1}}()    # Map color HEX codes: {CF sort (percap/overall), {HEX code}}
base_map = Dict{Int, Dict{String, Dict{}}}()    # Base map: {year, {nation A3, {base map Dict}}}
cf_maps = Dict{Int, Dict{String, Dict{String, Dict{}}}}() # categorized CF maps: {year, {nation A3, {map Dict by category}}}

function readBaseMap(year, nation, map_file; remove_feat = true, remove_reg = true, alter = true, label_conv = true)
    # remove: (default)[true] remain only "geometry", "properties", and "type" features
    # alter: (default)[true] change "GIS_ID" to "KEY_CODE" and "GIS_name" to "EN_NAME", and add "fill_carbon"
    # label_conv: (default)[true] change "KEY_CODE"("GIS_ID")'s gis id to gis label

    global base_map, reg_id, reg_list
    if !haskey(base_map, year); base_map[year] = Dict{String, Dict{}}() end

    ess_ft = ["geometry", "properties", "type"] # essential features
    ess_pr = ["GIS_name", "GIS_ID"]             # essential properties
    bs_map = js.geo2dict(js.read(read(map_file)))

    ft = bs_map["features"]
    nf = length(ft)

    if remove_feat
        for i = 1:nf, rf in filter(x -> !(x in ess_ft), collect(keys(ft[i]))); delete!(ft[i], rf) end
    end

    if remove_reg
        filter!(x -> string(x["properties"]["GIS_ID"]) in reg_list[year][nation], ft)
        nf = length(ft)
    end

    if alter
        for i = 1:nf
            fp = ft[i]["properties"]
            for rp in filter(x -> !(x in ess_pr), collect(keys(fp))); delete!(fp, rp) end
            fp["KEY_CODE"] = string(fp["GIS_ID"])
            fp["EN_NAME"] = fp["GIS_name"]
            fp["fill_carbon"] = ""
            delete!(fp, "GIS_ID")
            delete!(fp, "GIS_name")
        end
    end

    if label_conv
        rid = reg_id[year][nation]
        id_key = (alter ? "KEY_CODE" : "GIS_ID")
        for i = 1:nf
            fp = ft[i]["properties"]
            fp[id_key] = rid[string(fp[id_key])]
        end
    end

    base_map[year][nation] = bs_map
end

function readMapInfo(input_file)

    global file_names, map_list, field_label
    title_labels = ["CF_type", "Category", "File_name", "Field_label"]

    f_sep = getValueSeparator(input_file)
    f = open(input_file)

    title = string.(strip.(split(readline(f), f_sep)))
    ti, ci, fi, li = [findfirst(x -> x == it, title) for it in title_labels]

    for l in eachline(f)
        s = string.(strip.(split(l, f_sep)))
        ml = s[ti] * "_" * s[ci]
        file_names[ml] = s[fi]
        field_label[ml] = s[li]
        push!(map_list, ml)
    end

    close(f)
end

function readColorMap(rgbFile; reverse=false)

    rgb = Array{Tuple{Int, Int, Int}, 1}()

    f = open(rgbFile)

    n = Base.parse(Int, split(readline(f), r"[ \t]", keepempty=false)[end])
    idx_l = string.(split(replace(readline(f), "#" => ""), r"[ \t]", keepempty=false))
    ri, gi, bi = [findfirst(x -> x == i, idx_l) for i in ["r", "g", "b"]]

    cnt = 0
    for l in eachline(f)
        if isdigit(l[1])
            s = split(l, r"[ \t]", keepempty=false)[[ri, gi, bi]]
            if all(['.' in x for x in s].==false); push!(rgb, Tuple([Base.parse(Int, x) for x in s]))
            else push!(rgb, Tuple([convert(Int, floor(Base.parse(Float64, x) * 255, digits = 0)) for x in s]))
            end
            cnt += 1
        end
    end
    if reverse; reverse!(rgb) end
    close(f)

    if cnt == n; return rgb
    else println("Mismatching RGB numbers: ", n, "\t", cnt)
    end
end

function convertRgbToHex(rgbs::Array{Tuple{Int, Int, Int}, 1}; mode = "")
    # mode: "percap" or "overall"
    global hex_codes

    hex_cds = Array{String, 1}()

    for rgb in rgbs
        hex = "#"
        for c in rgb
            if c < 16; hex *= "0" end
            hex *= string(c, base = 16)
        end
        push!(hex_cds, hex)
    end

    if mode in ["percap", "overall"]; hex_codes[mode] = hex_cds
    else println("Incorrect mode: ", mode)
    end
end

function importEmissionData(ec_data::Module; emission = "cf", pc_dev = true, ov_dev = false)
    # emission: "cf" (total), "ie" (embedded), or "de" (direct)
    # pc_dev: [true] categozied emission per capita deviation from mean rank by region
    # ov_dev: [false] categozied emission rank by region

    global cat_list = ec_data.cat_list
    global reg_list = ec_data.gisRegList
    global reg_id = ec_data.gisRegID
    global pc_reg_cf, ov_reg_cf

    if emission == "cf"
        if pc_dev; pc_reg_cf = ec_data.cfRegDevPcRankGIS; else pc_reg_cf = ec_data.cfRegPcRankGIS end
        if ov_dev; ov_reg_cf = ec_data.cfRegDevRankGIS; else ov_reg_cf = ec_data.cfRegRankGIS end
    elseif emission == "ie"
        if pc_dev; pc_reg_cf = ec_data.ieRegDevPcRankGIS; else pc_reg_cf = ec_data.ieRegPcRankGIS end
        if ov_dev; ov_reg_cf = ec_data.ieRegDevRankGIS; else ov_reg_cf = ec_data.ieRegRankGIS end
    elseif emission == "de"
        if pc_dev; pc_reg_cf = ec_data.deRegDevPcRankGIS; else pc_reg_cf = ec_data.deRegPcRankGIS end
        if ov_dev; ov_reg_cf = ec_data.deRegDevRankGIS; else ov_reg_cf = ec_data.deRegRankGIS end
    else println("Incorrect emission mode: ", emission)
    end
end

function mapRegionCF(year, nation; label_conv = true, blank_color = "#EDEDED", value_mode = true)
    # value_mode: [true] store "CFAC" (rank) values in the properties field

    global cat_list, map_list, reg_list, reg_id, field_label
    global pc_reg_cf, ov_reg_cf, hex_codes
    global base_map, cf_maps

    rls, rid = reg_list[year][nation], reg_id[year][nation]
    if !haskey(cf_maps, year); cf_maps[year] = Dict{String, Dict{String, Dict{}}}() end
    if !haskey(cf_maps[year], nation); cf_maps[year][nation] = Dict{String, Dict{}}() end

    for ml in map_list
        typ, cat = string.(strip.(split(ml, "_", limit = 2)))
        ci = findfirst(x -> x == cat, cat_list)
        cmap = deepcopy(base_map[year][nation])

        if !(typ in ["percap", "overall"]); println("Incorrect CF type: ", typ) end
        if value_mode
            val_field = "CFAC_" * field_label[ml]
            blank_value = 0
        end

        ft = cmap["features"]

        for i = 1:length(ft)
            fp = ft[i]["properties"]
            gid = fp["KEY_CODE"]
            ri = findfirst(x -> (label_conv ? rid[x] : x) == gid, rls)
            if ri != nothing
                if typ == "overall"; rcf = ov_reg_cf[year][nation]
                elseif typ == "percap"; rcf = pc_reg_cf[year][nation]
                end
                fp["fill_carbon"] = hex_codes[typ][rcf[ri,ci]]
                if value_mode; fp[val_field] = rcf[ri,ci] end
            else
                fp["fill_carbon"] = blank_color
                if value_mode; fp[val_field] = blank_value end
            end
        end
        cf_maps[year][nation][ml] = cmap
    end
end

function printMapFiles(year, nation, map_filepath)

    global map_list, file_names, cf_maps
    mkpath(map_filepath)

    nc = length(map_list)
    maps = cf_maps[year][nation]

    for ml in map_list
        f = open(map_filepath * file_names[ml], "w")
        println(f, js.write(js.dict2geo(maps[ml])))
        close(f)
    end
end

function getValueSeparator(file_name)
    fext = file_name[findlast(isequal('.'), file_name)+1:end]
    if fext == "csv"; return ','
    elseif fext == "tsv" || fext == "txt"; return '\t' end
end

end
