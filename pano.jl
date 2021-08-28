using Printf, LinearAlgebra
import FixedPointNumbers: N0f8
import Images
import CSV, DataFrames

#===========================================
  ____            _   _ _   _ _
 / ___| ___  ___ | | | | |_(_) |___
| |  _ / _ \/ _ \| | | | __| | / __|
| |_| |  __/ (_) | |_| | |_| | \__ \
 \____|\___|\___/ \___/ \__|_|_|___/

=============================================#

struct PositionLLH
    lat::Float64
    lon::Float64
    height::Float64
end

struct PositionXYZ
    x::Float64
    y::Float64
    z::Float64
end

### WGS84 utils
struct Ellipsoid
    a::Float64
    b::Float64
    f::Float64
    e::Float64
    e²::Float64
end

function createWgs84Ellipsoid()
    a  = 6378137.0
    f  = 1.0/298.257223563
    b  = a*(1.0-f)
    e² = 1.0 - b*b/(a*a) # or f * (2.0-f)
    e  = sqrt(e²)
    return Ellipsoid(a,b,f,e,e²)
end

wgs84 = createWgs84Ellipsoid()

toRadians(α) = α / 180.0 * π
toDegrees(α) = α * 180.0 / π

function llh_to_xyz(ellipsoid::Ellipsoid, llh::PositionLLH)::PositionXYZ
    Φ    = toRadians(llh.lat)
    λ    = toRadians(llh.lon)
    h    = llh.height
    e²   = ellipsoid.e²
    sinΦ = sin(Φ)
    cosΦ = cos(Φ)
    sinλ = sin(λ)
    cosλ = cos(λ)
    v = ellipsoid.a / sqrt(1.0 - e² * sinΦ * sinΦ)
    x = (v+h) * cosΦ * cosλ
    y = (v+h) * cosΦ * sinλ
    z = ((1.0 - e²) * v + h) * sinΦ;
    return PositionXYZ(x,y,z)
end

# TODO: make it work for south and west
function xyz_to_llh(ellipsoid::Ellipsoid, xyz::PositionXYZ)::PositionLLH
    e²       = ellipsoid.e²
    p        = sqrt(xyz.x*xyz.x + xyz.y*xyz.y)
    lon      = atan(xyz.y / xyz.x)
    lat_init = atan(xyz.z/(p*(1.0 - e²)))
    v        = ellipsoid.a / sqrt(1.0-e²*sin(lat_init)*sin(lat_init))
    lat      = atan((xyz.z + e²*v*sin(lat_init))/p)
    height   = (p/cos(lat))-v
    return PositionLLH(toDegrees(lat), toDegrees(lon), height)
end    

# GeoUtils test
#testXYZ = llh_to_xyz(wgs84, PositionLLH(49.0, 16.0, 225.0))
#println(testXYZ)
#testLLH= xyz_to_llh(wgs84, testXYZ)
#println(testLLH)

#================================================================
 _   _      _       _     _   __  __
| | | | ___(_) __ _| |__ | |_|  \/  | __ _ _ __
| |_| |/ _ \ |/ _` | '_ \| __| |\/| |/ _` | '_ \
|  _  |  __/ | (_| | | | | |_| |  | | (_| | |_) |
|_| |_|\___|_|\__, |_| |_|\__|_|  |_|\__,_| .__/
              |___/                       |_|
===============================================================#


struct LatLonRange
    minLat::Int16
    minLon::Int16
    maxLat::Int16
    maxLon::Int16
end


getHgtFileName(lat, lon) = @sprintf "N%02dE%03d.hgt" lat lon


getHgtFilePath(lat, lon, tileDir) = @sprintf "%s\\%s" tileDir getHgtFileName(lat, lon)


function loadData(range::LatLonRange, tileDir)
    nTilesHoriz = range.maxLon - range.minLon + 1
    nTilesVert  = range.maxLat - range.minLat + 1
    nTilesTotal = nTilesHoriz * nTilesVert
    dataWidth   = nTilesHoriz*1200+1
    dataHeight  = nTilesVert*1200+1

    println(@sprintf("Requesting data for area %d°N %d°E - %d°N %d°E ... (aprox. %3.0fx%3.0f km).",
        range.minLat, range.minLon,
        range.maxLat, range.maxLon,
        (nTilesHoriz)*111.1*cos(range.minLat/180.0*π),
        (nTilesVert)*111.1
        ))
    println(@sprintf("I will read %dx%d=%d tiles, heightmap size is %dx%d (%d MB).",
        nTilesHoriz, nTilesVert, nTilesTotal,
        dataWidth, dataHeight, dataWidth*dataHeight*2/1000000
        ))

    # Note: array is indexed by row, column and starting index is 1
    data = Array{UInt16}(undef, dataHeight, dataWidth)
    progress = 0
    for lat in range.minLat:range.maxLat
        Threads.@threads for lon in range.minLon:range.maxLon
            # Print progress
            progress = progress + 1
            println(@sprintf("Loading tile %03d/%03d lat=%02d, lon=%02d", progress, nTilesTotal, lat, lon))
            # Load tile
            tile = Array{Int16}(undef, 1201, 1201)
            path = getHgtFilePath(lat, lon, tileDir)
            io = open(path, "r")
            read!(io, tile)
            close(io)
            # Fix endianity, clamp, transpose and copy to data
            tile .= ntoh.(tile)
            clamp!(tile, 0, 6000)
            tile = transpose(tile)
            rOffset = (range.maxLat-lat) * 1200
            cOffset = (lon-range.minLon) * 1200
            for r in 1:1201
                @simd for c in 1:1201
                    data[rOffset+r, cOffset+c] = tile[r, c]
                end
            end
        end
    end
    return data
end


function saveHeightMap(data)
   minValue = minimum(data)
   maxValue = maximum(data)
   norm = Gray.(data/maxValue)
   println("Saving heightmap-gray.png")
   save("heightmap-gray.png", norm)
end


# getHeight(srtmRange, heightMap, 49.142158, 16.627978) -> 192 Svratka, Svitava
function getHeight(range::LatLonRange, data::Matrix{UInt16}, lat::Float64, lon::Float64)
    r = Int64(trunc((range.maxLat+1 - lat)*1200))
    c = Int64(trunc((lon-range.minLon)*1200))
    r = clamp(r, 1, size(data)[1])
    c = clamp(c, 1, size(data)[2])
    return data[r,c]
end


makeEarthCurve(radius, distMax, distStep) = [sqrt(radius*radius-x*x)-radius for x=range(0, distMax, step=distStep)]


function makeDistMap(eye::PositionLLH, latLonRange::LatLonRange, heightMap::Matrix{UInt16})

    wgs84  = createWgs84Ellipsoid()
    eyeXYZ = llh_to_xyz(wgs84, eye)
    
    angleMin  = toRadians( 90.0)
    angleMax  = toRadians(135.0)
    angleStep = 0.0001
    distStep  = 100.0
    distMax   = 250.0*1000.0

#    vertAngleMin = toRadians(-5.0)
#    vertAngleMax = toRadians( 2.0)
# This matches reference output from panorama generator
    vertAngleMin = -0.0560
    vertAngleMax =  0.0339

    pRef   = [eyeXYZ.x, eyeXYZ.y, eyeXYZ.z]
    vZ     = [0.0, 0.0, 1.0]
    vDown  = -normalize(pRef)
    vEast  = normalize(cross(vDown,vZ))
    vNorth = normalize(cross(vEast,vDown))
# Note: 1.18 matches photo during winter inversion quite well
    diffractionModifier = 1.18
    earthRadius = sqrt(dot(pRef, pRef)) * diffractionModifier
    earthCurve  = makeEarthCurve(earthRadius, distMax, distStep)
    xMax        = Int64(trunc( (angleMax-angleMin)/angleStep ))+1
    outWidth    = xMax+1
    outHeight   = size(range(vertAngleMin, vertAngleMax, step=angleStep))[1]+1
    println(@sprintf("Earth radius is %6.1f km (diffraction x%4.2f)", earthRadius/diffractionModifier/1000.0, diffractionModifier))
    println(@sprintf("Output size is %d x %d pixels", outWidth, outHeight))
    println(@sprintf("Output resolution is %f mrad per pixel or %f pixels per degree", angleStep * 1000.0, 1/toDegrees(angleStep)))
    output      = zeros(UInt16, outHeight, outWidth)
    distances   = range(0, distMax, step=distStep)
    Threads.@threads for x in 0:xMax
        #println(@sprintf("Rendering line %d of %d", x, xMax))
        azimuth = angleMin + x*angleStep
        cosAz   = cos(azimuth)
        sinAz   = sin(azimuth)

        vertAngle     = vertAngleMin
        h0            = eye.height
        rayCastHeight = h0
        index         = 1
        direction     = vNorth*cosAz + vEast*sinAz
        for dist in distances
            #point = pRef + dist * direction;
            point = [pRef[1]+dist*direction[1], pRef[2]+dist*direction[2], pRef[3]+dist*direction[3]]
            llh = xyz_to_llh(wgs84, PositionXYZ(point[1], point[2], point[3]))
            rayCastHeight = h0 + sin(vertAngle) * dist
            terrainHeight = earthCurve[index] + getHeight(latLonRange, heightMap, llh.lat, llh.lon)
            if terrainHeight > rayCastHeight 
                newVertAngle = atan((terrainHeight-h0)/dist)
                yTop = Int64(trunc( (vertAngleMax-newVertAngle)/angleStep ))
                yBot = Int64(trunc( (vertAngleMax-   vertAngle)/angleStep ))
                v = UInt16(trunc( dist / distStep )) 
                for y in yTop:yBot 
                    output[y, x+1] = v
                end
                vertAngle = newVertAngle
            end
            index = index + 1
        end
    end

    println("SUMMIT TEST CODE ---- ")
    # TEST code
    hills = CSV.File("data-cz-prom100.tsv") |> DataFrames.DataFrame
    # convert to ours azimuth, angle above horizon and distance - project into 
    hill_to_xyz(ellipsoid::Ellipsoid, dfRow)::PositionXYZ = llh_to_xyz(ellipsoid, PositionLLH(dfRow["Latitude"], dfRow["Longitude"], dfRow["Elevation"])) 
    mLocalToWorld = hcat(vEast,vNorth,-vDown)
    mWorldToLocal = inv(mLocalToWorld)
    for hill in eachrow(hills)
        hill_world = hill_to_xyz(wgs84, hill)
        hill_local_xyz = mWorldToLocal * [hill_world.x, hill_world.y, hill_world.z]
        hill_local_xyz[3] = hill_local_xyz[3] - sqrt(dot(pRef, pRef))
        distance = sqrt(dot(hill_local_xyz, hill_local_xyz))
        azimuth  = toDegrees(atan(hill_local_xyz[1], hill_local_xyz[2])) # goes from north to east so y/x is swapped
        #print(hill["Summit"])
        if azimuth < 0
            azimuth = azimuth + 360
        end
        if (azimuth < toDegrees(angleMin) || azimuth > toDegrees(angleMax))
#            println(@sprintf(" is out of visible azimuths (%f not in interval from %f tp %f", azimuth, toDegrees(angleMin), toDegrees(angleMax)))
            continue
        end
        if (distance > distMax)
 #           println(@sprintf(" is out of range (%f is greater than %f)", distance, distMax))
            continue
        end
        println(@sprintf("%20s is possibly visible at azimuth %5.1f, distance %5.1f km", hill["Summit"], azimuth, distance/1000.0))
    end
    
    return output
end


function extractHorizon(distMap::Matrix{UInt16})
    nRows = size(distMap)[1]
    nCols = size(distMap)[2]
    #ouput = Array{Int16}(undef, nRows, nCols)
    output = zeros(N0f8, nRows, nCols)
    Threads.@threads for row in 2:nRows
        for col in 1:nCols
            diff::UInt8 = 255-clamp(abs(reinterpret(Int16,distMap[row-1,col]) - reinterpret(Int16,distMap[row,col])), 0, 255)
            output[row, col] = reinterpret(N0f8, diff)
        end
    end
    return output
end


function tryMountains(distMap::Matrix{UInt16})

end

# Info (from https://www.udeuschle.de/panoramas/makepanoramas_en.htm)
# Lat: 50.08309 Lon 17.23094 Alt(auto+10m): 1500+10
# View direction: 112.5, extension 45 left: 90, right: 135, resolution 20pix/deg
# Tilt, range, vert. exaggeration 1.2

# TODO: determine lat/long range automatically
function main()
    tileDir     = "d:\\_disk_d_old\\devel-python\\panorama\\data_srtm"
    latLonRange = LatLonRange(47, 15, 50, 21)
    eyePos      = PositionLLH(50.08309, 17.23094, 1510)

    heightMap = loadData(latLonRange, tileDir)
    #saveHeightMap(data)
    distMap   = makeDistMap(eyePos, latLonRange, heightMap)

    minValue = minimum(distMap)
    maxValue = maximum(distMap)
    println(@sprintf("min=%d max=%d", minValue, maxValue))
    println("Saving distmap-gray.png")
    Images.save("disttmap-gray.png", Images.Gray.(distMap/maxValue))

    println("Saving horizon.png....")
    hrz=extractHorizon(distMap)
    Images.save("horizon.png", Images.Gray.(hrz))

    # This can be fun: https://wiki.flightgear.org/Atmospheric_light_scattering
    # http://www.science-and-fiction.org/rendering/als.html

end

main()
