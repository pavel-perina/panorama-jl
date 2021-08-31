#=
Julia script to convert data of hills in Czech Republic into useable format ("N50° 44.161 E15° 44.380" is not a nice field)

Original data were downloaded from http://www.peaklist.org/WWlists/euro600/slovakia/Slovakia200m.html
parsed via https://pdftables.com to XLS
and saved in OpenOffice as CSV file after minor modifications

Author: Pavel Perina
Date:   2021-08-31
=#

import CSV
import DataFrames

dfIn  = CSV.File("Slovakia_P200_prominence.csv"; type=String) |> DataFrames.DataFrame
dfOut = DataFrames.DataFrame(Summit = String[], Elevation = Int[], Latitude=Float64[], Longitude=Float64[])

for row in eachrow(dfIn)
    println(row["Peak"])
    mountain  = row["Peak"]
    elevation = parse(Int,row["Height"])
    n = match(r"(\d*):(\d*):(\d*)", row["Lat"])
    north = parse(Float64, n[1]) + parse(Float64, n[2])/60.0 + parse(Float64, n[3])/3600.0
    e = match(r"(\d*):(\d*):(\d*)", row["Long"])
    east  = parse(Float64, e[1]) + parse(Float64, e[2])/60.0 + parse(Float64, e[3])/3600.0
    push!(dfOut, (mountain, elevation, north, east))
end
CSV.write("data-sk-prom200.tsv",dfOut; delim='\t',quotestrings=true)