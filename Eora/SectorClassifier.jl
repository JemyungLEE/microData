# Developed date: 30. Aug. 2019
# Last modified date: 2. Sep. 2019
# Subject: Eora sector classification
# Description: classify the categories of Eora sectors by nation
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

module SectorClassifier

struct sector
    index::Int
    country::String
    abb::String          # abbreviation of country name
    entity::String      # 'industries' or 'commodities'
    categ::String       # categeory of this sector

    function sector(str::String)
        idx, name, abb, ent, sec = split(str, '\t')
        new(parse(Int, idx), name, abb, ent, sec)
    end
end

mutable struct nation
    name::String
    abb::String         # abbreviation of country name
    ns::Int8            # number of sectors
    hasComEn::Bool      # wether the nation have 'commodities'-entity data
    isSimple::Bool      # wether a nation's sectors are classified by the basic 26 categories
    sectors::Array{String,1}

    nation(nstr::String, astr:;String) = new(nstr, astr, 0, false, false, Array{String, 1})
end

global nc::Int8     # number of counties
global nations = Dict{String, nation}()

function readSectorData(inputFile)
    # check the number of nations and store their names, abbreviation,
    #  and number of sections of each nation.

    f = open(inputFile)

    # check nations and wether they have commodity accounts
    for l in eachline(f)
        s = sector(l)
        if !haskey(nationList, s.abb)
            nationList[s.abb] = nation(s.country, s.abb)
        elseif entity == "commodities"
                nationList[s.abb].hasComEn = true
        end
    end



end
