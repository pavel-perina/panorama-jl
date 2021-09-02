import JSON
using  Printf

elements = JSON.Parser.parsefile("e:/dev-jl/pano/data-osm/osm-cz-sk.json ")["elements"]
@printf("Found %d elements.", size(elements)[1])

io = open("e:/dev-jl/pano/data-osm/osm-cz-sk.tsv", "w")
@printf(io,"\"Summit\"\t\"Elevation\"\t\"Latitude\"\t\"Longitude\"\n")

for i in 1:size(elements)[1]
    ele::String = get(elements[i]["tags"], "ele", "")
    if length(ele) == 0
        continue
    end
    name::String = get(elements[i]["tags"], "name", "")
    if length(name) == 0
        continue
    end
    lat = elements[i]["lat"]
    lon = elements[i]["lon"]
    @printf(io, "\"%s\"\t%0.1f\t%0.5f\t%0.5f\n", name, parse(Float64, ele), lat, lon)
end
close(io)

