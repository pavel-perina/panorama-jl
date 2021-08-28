#=
Julia script to convert data of hills in Czech Republic into useable format ("N50° 44.161 E15° 44.380" is not a nice field)

Original data were downloaded from https://www.ultratisicovky.cz/ke-stazeni/
and saved in OpenOffice as CSV file and header was manually modified.

Author: Pavel Perina
Date:   2021-08-27
=#

import CSV
import DataFrames

dfIn  = CSV.File("cz-prom100.csv",select=[4,5,12],type=String) |> DataFrames.DataFrame
dfOut = DataFrames.DataFrame(Summit = String[], Elevation = Float64[], Latitude=Float64[], Longitude=Float64[])

nRows = size(dfIn)[1]
for row in 1:nRows
#    println(row)
    mountain  = replace(dfIn[row,1], r"^\s+|\s+$|\s+(?=\s)" => "")
    elevation = parse(Float64, match(r"(\d+)", dfIn[row,2])[1])
    m = match(r"^(\s*)N(?<ndeg>\d*)° (?<nmin>([0-9]*[.])?[0-9]*) E(?<edeg>\d*)° (?<emin>([0-9]*[.])?[0-9]*)", dfIn[row,3])
    north = parse(Float64, m["ndeg"]) + parse(Float64, m["nmin"])/60.0
    east  = parse(Float64, m["edeg"]) + parse(Float64, m["emin"])/60.0
    push!(dfOut, (mountain, elevation, north, east))
end
CSV.write("data-cz-prom100.tsv",dfOut; delim='\t',quotestrings=true)