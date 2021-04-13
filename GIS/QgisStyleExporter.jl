module QgisStyleExporter

# Developed date: 13. Feb. 2020
# Last modified date: 29. May. 2020
# Subject: Export QGIS style file(s)
# Description: Make QML (QGIS style) file(s) containing 'categories' and 'symbols' data
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

nsym = 0
rgb = Array{Tuple{Int, Int, Int}, 1}()      # Symbols' RGB

function readColorMap(inputFile; reverse=false)

    global nsym = 0
    global rgb = Array{Tuple{Int, Int, Int}, 1}()

    f = open(inputFile)
    for l in eachline(f)
        if isdigit(l[1])
            s = split(l, r"[ \t]", keepempty=false)
            if all(['.' in x for x in s].==false); push!(rgb, Tuple([parse(Int, x) for x in s]))
            else push!(rgb, Tuple([convert(Int, round(parse(Float64, x)*255, digits=0)) for x in s]))
            end
        end
    end
    nsym = length(rgb)
    close(f)
    if reverse; reverse!(rgb) end
end

function makeQML(outputFile, attr::String; empty=false, labels=[], values=[], indexValue=false, labelReverse=false)  # attr=attribute field

    global nsym, rgb
    if empty; pushfirst!(rgb, (204, 204, 204)) end # set polygon style for no-data cells
    if empty; csi=1 else csi=0 end  #color starting index
    isLabel = length(labels)>0
    isValue = length(values)>0

    if !isLabel
        labels = Array{String,1}(undef,nsym+csi)
        if isValue
            nv = length(values)
            if empty; labels[1] = "No data" end
            for i=1:nv-1
                labels[i+csi] = string(round(values[i],sigdigits=3))
                if values[i]>= 0.0001
                    if round(values[i],sigdigits=3) == round(values[i],sigdigits=1); labels[i+csi] *= "00"
                    elseif round(values[i],sigdigits=3) == round(values[i],sigdigits=2); labels[i+csi] *= "0"
                    end
                elseif values[i]< 0.0001 && round(values[i],sigdigits=3)==round(values[i],sigdigits=2)
                    labels[i+csi] = replace(labels[i+csi], "e"=>"0e")
                end
            end
            labels[nv+csi] = string(round(values[end],sigdigits=3))
            if values[end]>0.0001
                if round(values[end],sigdigits=3) == round(values[end],sigdigits=1); labels[nv+csi] *= "00"
                elseif round(values[end],sigdigits=3) == round(values[end],sigdigits=2); labels[nv+csi] *= "0"
                end
            elseif values[end]< 0.0001 && round(values[end],sigdigits=3)==round(values[end],sigdigits=2)
                labels[nv+csi] = replace(labels[nv+csi], "e"=>"0e")
            end
        else for i=1:nsym+csi; labels[i] = i-csi end
        end
        labels[1+csi] *= "less "
        labels[end] *= "over "
    elseif empty; pushfirst!(labels, "No data")
    end
    if indexValue || !isValue
        values=zeros(Int, nsym+csi)
        for i=1:nsym+csi; values[i] = i-csi end
    end

    if labelReverse
        if empty; lb = []; push!(lb, labels[1]); labels = append(lb, reverse(labels))
        else reverse!(labels)
        end
    end

    f = open(outputFile, "w")

    # print head
    println(f,"<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>")
    println(f,"<qgis version=\"3.4.15-Madeira\" styleCategories=\"Symbology\">")
    println(f,"  <renderer-v2 type=\"categorizedSymbol\" enableorderby=\"0\" attr=\""*attr*"\" symbollevels=\"0\" forceraster=\"0\">")

    # print categories
    println(f,"    <categories>")
    if isLabel
        for i=1:nsym+csi
            println(f,"      <category render=\"true\" label=\"",labels[i],"\" symbol=\"",i-1,"\" value=\"",values[i],"\"/>")
        end
    else
        for i=1:nsym+csi
            if i==csi || i==csi+1 || i==nsym+csi
                println(f,"      <category render=\"true\" label=\"",labels[i],"\" symbol=\"",i-1,"\" value=\"",values[i],"\"/>")
            else println(f,"      <category render=\"true\" label=\"",labels[i-1],"-",labels[i],"\" symbol=\"",i-1,"\" value=\"",values[i],"\"/>")
            end
        end
    end
    println(f, "    </categories>")

    # print symbols
    println(f,"    <symbols>")
    for i=1:nsym+csi
        println(f,"      <symbol type=\"fill\" clip_to_extent=\"1\" name=\"",i-1,"\" force_rhr=\"0\" alpha=\"1\">")
        println(f,"        <layer pass=\"0\" class=\"SimpleFill\" enabled=\"1\" locked=\"0\">")
        println(f,"          <prop k=\"border_width_map_unit_scale\" v=\"3x:0,0,0,0,0,0\"/>")
        println(f,"          <prop k=\"color\" v=\"",rgb[i][1],",",rgb[i][2],",",rgb[i][3],",255\"/>")
        println(f,"          <prop k=\"joinstyle\" v=\"bevel\"/>")
        println(f,"          <prop k=\"offset\" v=\"0,0\"/>")
        println(f,"          <prop k=\"offset_map_unit_scale\" v=\"3x:0,0,0,0,0,0\"/>")
        println(f,"          <prop k=\"offset_unit\" v=\"MM\"/>")
        println(f,"          <prop k=\"outline_color\" v=\"110,110,110,255\"/>")
        println(f,"          <prop k=\"outline_style\" v=\"solid\"/>")
        println(f,"          <prop k=\"outline_width\" v=\"0.000001\"/>")
        println(f,"          <prop k=\"outline_width_unit\" v=\"MM\"/>")
        println(f,"          <prop k=\"style\" v=\"solid\"/>")
        println(f,"          <data_defined_properties>")
        println(f,"            <Option type=\"Map\">")
        println(f,"              <Option type=\"QString\" name=\"name\" value=\"\"/>")
        println(f,"              <Option name=\"properties\"/>")
        println(f,"              <Option type=\"QString\" name=\"type\" value=\"collection\"/>")
        println(f,"            </Option>")
        println(f,"          </data_defined_properties>")
        println(f,"        </layer>")
        println(f,"      </symbol>")
    end
    println(f,"    </symbols>")

    # print source-symbol
    println(f,"    <source-symbol>")
    println(f,"      <symbol type=\"fill\" clip_to_extent=\"1\" name=\"0\" force_rhr=\"0\" alpha=\"1\">")
    println(f,"        <layer pass=\"0\" class=\"SimpleFill\" enabled=\"1\" locked=\"0\">")
    println(f,"          <prop k=\"border_width_map_unit_scale\" v=\"3x:0,0,0,0,0,0\"/>")
    println(f,"          <prop k=\"color\" v=\"",rgb[1+csi][1],",",rgb[1+csi][2],",",rgb[1+csi][3],",255\"/>")
    println(f,"          <prop k=\"joinstyle\" v=\"bevel\"/>")
    println(f,"          <prop k=\"offset\" v=\"0,0\"/>")
    println(f,"          <prop k=\"offset_map_unit_scale\" v=\"3x:0,0,0,0,0,0\"/>")
    println(f,"          <prop k=\"offset_unit\" v=\"MM\"/>")
    println(f,"          <prop k=\"outline_color\" v=\"110,110,110,255\"/>")
    println(f,"          <prop k=\"outline_style\" v=\"solid\"/>")
    println(f,"          <prop k=\"outline_width\" v=\"0.000001\"/>")
    println(f,"          <prop k=\"outline_width_unit\" v=\"MM\"/>")
    println(f,"          <prop k=\"style\" v=\"solid\"/>")
    println(f,"          <data_defined_properties>")
    println(f,"            <Option type=\"Map\">")
    println(f,"              <Option type=\"QString\" name=\"name\" value=\"\"/>")
    println(f,"              <Option name=\"properties\"/>")
    println(f,"              <Option type=\"QString\" name=\"type\" value=\"collection\"/>")
    println(f,"            </Option>")
    println(f,"          </data_defined_properties>")
    println(f,"        </layer>")
    println(f,"      </symbol>")
    println(f,"    </source-symbol>")

    # print color-ramp type
    println(f,"    <colorramp type=\"gradient\" name=\"[source]\">")
    println(f,"      <prop k=\"color1\" v=\"",rgb[1+csi][1],",",rgb[1+csi][2],",",rgb[1+csi][3],",255\"/>")
    println(f,"      <prop k=\"color2\" v=\"",rgb[nsym][1],",",rgb[nsym][2],",",rgb[nsym][3],",255\"/>")
    println(f,"      <prop k=\"discrete\" v=\"0\"/>")
    println(f,"      <prop k=\"rampType\" v=\"gradient\"/>")
    idx = convert(Int, round(0.25*(nsym-csi), digits=0))+csi
    print(f,"      <prop k=\"stops\" v=\"0.25;",rgb[idx][1],",",rgb[idx][2],",",rgb[idx][3],",255")
    idx = convert(Int, round(0.5*(nsym-csi), digits=0))+csi
    print(f,":0.5;",rgb[idx][1],",",rgb[idx][2],",",rgb[idx][3],",255")
    idx = convert(Int, round(0.75*(nsym-csi), digits=0))+csi
    println(f,":0.75;",rgb[idx][1],",",rgb[idx][2],",",rgb[idx][3],",255\"/>")
    println(f,"    </colorramp>")

    # others
    println(f,"    <rotation/>")
    println(f,"    <sizescale/>")
    println(f,"  </renderer-v2>")
    println(f,"  <blendMode>0</blendMode>")
    println(f,"  <featureBlendMode>0</featureBlendMode>")
    println(f,"  <layerGeometryType>2</layerGeometryType>")
    println(f,"</qgis>")

    close(f)
end

end
