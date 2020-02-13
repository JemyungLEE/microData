module QgisStyleExporter

# Developed date: 13. Feb. 2020
# Last modified date: 13. Feb. 2020
# Subject: Export QGIS style file(s)
# Description: Make QML (QGIS style) file(s) containing 'categories' and 'symbols' data
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

nsym = 0
rgb = Array{Tuple{Int, Int, Int}, 1}()      # Symbols' RGB

function readColorMap(inputFile)

    global nsym, rgb

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
end

function makeQML(outputFile, attr::String)  # attr=attribute field

    f = open(outputFile, "w")

    # print head
    println(f,"<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>")
    println(f,"<qgis version=\"3.4.15-Madeira\" styleCategories=\"Symbology\">")
    println(f,"  <renderer-v2 type=\"categorizedSymbol\" enableorderby=\"0\" attr=\""*attr*"\" symbollevels=\"0\" forceraster=\"0\">")

    # print categories
    println(f,"    <categories>")
    for i=1:nsym; println(f,"      <category render=\"true\" label=\"",i,"\" symbol=\"",i-1,"\" value=\"",i,"\"/>") end
    println(f, "    </categories>")

    # print symbols
    println(f,"    <symbols>")
    for i=1:nsym
        println(f,"      <symbol type=\"fill\" clip_to_extent=\"1\" name=\"",i-1,"\" force_rhr=\"0\" alpha=\"1\">")
        println(f,"        <layer pass=\"0\" class=\"SimpleFill\" enabled=\"1\" locked=\"0\">")
        println(f,"          <prop k=\"border_width_map_unit_scale\" v=\"3x:0,0,0,0,0,0\"/>")
        println(f,"          <prop k=\"color\" v=\"",rgb[i][1],",",rgb[i][2],",",rgb[i][3],",255\"/>")
        println(f,"          <prop k=\"joinstyle\" v=\"bevel\"/>")
        println(f,"          <prop k=\"offset\" v=\"0,0\"/>")
        println(f,"          <prop k=\"offset_map_unit_scale\" v=\"3x:0,0,0,0,0,0\"/>")
        println(f,"          <prop k=\"offset_unit\" v=\"MM\"/>")
        println(f,"          <prop k=\"outline_color\" v=\"35,35,35,255\"/>")
        println(f,"          <prop k=\"outline_style\" v=\"solid\"/>")
        println(f,"          <prop k=\"outline_width\" v=\"0.05\"/>")
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
    println(f,"          <prop k=\"color\" v=\"",rgb[1][1],",",rgb[1][2],",",rgb[1][3],",255\"/>")
    println(f,"          <prop k=\"joinstyle\" v=\"bevel\"/>")
    println(f,"          <prop k=\"offset\" v=\"0,0\"/>")
    println(f,"          <prop k=\"offset_map_unit_scale\" v=\"3x:0,0,0,0,0,0\"/>")
    println(f,"          <prop k=\"offset_unit\" v=\"MM\"/>")
    println(f,"          <prop k=\"outline_color\" v=\"35,35,35,255\"/>")
    println(f,"          <prop k=\"outline_style\" v=\"solid\"/>")
    println(f,"          <prop k=\"outline_width\" v=\"0.05\"/>")
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
    println(f,"      <prop k=\"color1\" v=\"",rgb[1][1],",",rgb[1][2],",",rgb[1][3],",255\"/>")
    println(f,"      <prop k=\"color2\" v=\"",rgb[nsym][1],",",rgb[nsym][2],",",rgb[nsym][3],",255\"/>")
    println(f,"      <prop k=\"discrete\" v=\"0\"/>")
    println(f,"      <prop k=\"rampType\" v=\"gradient\"/>")
    idx = convert(Int, round(0.25*nsym, digits=0))
    print(f,"      <prop k=\"stops\" v=\"0.25;",rgb[idx][1],",",rgb[idx][2],",",rgb[idx][3],",255")
    idx = convert(Int, round(0.5*nsym, digits=0))
    print(f,":0.5;",rgb[idx][1],",",rgb[idx][2],",",rgb[idx][3],",255")
    idx = convert(Int, round(0.75*nsym, digits=0))
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
