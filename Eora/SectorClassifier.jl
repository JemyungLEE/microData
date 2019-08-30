# Developed date: 30. Aug. 2019
# Last modified date: 30. Aug. 2019
# Subject: Eora sector classification
# Description: classify the categories of Eora sectors by nation
# Developer: Jemyung Lee
# Affiliation: RIHN (Research Institute for Humanity and Nature)

module SectorClassifier

mutable struct sectors
    index::Int
    country::String
    ct::String          # abbreviation of country name
    entity::String      # 'industries' or 'commodities'
    sector::String
end

mutable struct nations
    name::String
    abb::String         # abbreviation of country name
    ns::Int8            # number of sectors
    hasComEn::Bool      # wether the nation have 'commodities'-entity data
    isSimple::Bool      # wether a nation's sectors are classified by the basic 26 categories

end

global nc::Int8     # number of counties
global sectorList = []
global nationList

functoin check

function readSectorData(inputFile)
    f = open(inputFile)



end
