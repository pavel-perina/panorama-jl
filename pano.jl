using Printf
# dot,cross
using LinearAlgebra
using Images
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


struct SrtmRange
    minLat::Int16
    minLon::Int16
    maxLat::Int16
    maxLon::Int16
end


getHgtFileName(lat, lon) = @sprintf "N%02dE%03d.hgt" lat lon


getHgtFilePath(lat, lon, tileDir) = @sprintf "%s\\%s" tileDir getHgtFileName(lat, lon)


function loadData(range::SrtmRange, tileDir)
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
    data = Array{Int16}(undef, dataHeight, dataWidth)
    progress = 0
    for lat in range.minLat:range.maxLat
        for lon in range.minLon:range.maxLon
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
                for c in 1:1201
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

#============================
============================#

# getHeight(srtmRange, heightMap, 49.142158, 16.627978) -> 192 Svratka, Svitava
function getHeight(range::SrtmRange, data::Matrix{Int16}, lat::Float64, lon::Float64)
    r = Int64(trunc((range.maxLat+1 - lat)*1200))
    c = Int64(trunc((lon-range.minLon)*1200))
    r = clamp(r, 1, size(data)[1])
    c = clamp(c, 1, size(data)[2])
    return data[r,c]
end


function makeEarthCurve(radius, distMax, distStep)
    distances = range(0, distMax, step=distStep)
    [sqrt(radius*radius-x*x)-radius for x=distances]
end




function makeDistMap(eye::PositionLLH, srtmRange::SrtmRange, heightMap::Matrix{Int16})
    wgs84  = createWgs84Ellipsoid()
    eyeXYZ = llh_to_xyz(wgs84, eye)
    
    angleMin  = toRadians( 90.0)
    angleMax  = toRadians(135.0)
    angleStep = 0.0001
    distStep  = 100.0
    distMax   = 250.0*1000.0

    vertAngleMin = toRadians(-10.0)
    vertAngleMax = toRadians( 5.0)
    pRef   = [eyeXYZ.x, eyeXYZ.y, eyeXYZ.z]
    vZ     = [0.0, 0.0, 1.0]
    vDown  = -normalize(pRef)
    vEast  = normalize(cross(vDown,vZ))
    vNorth = normalize(cross(vEast,vDown))
# Note: 1.18 matches photo during winter inversion quite well
    diffractionModifier = 1.18
    earthRadius = sqrt(dot(pRef, pRef)) * diffractionModifier
    earthCurve = makeEarthCurve(earthRadius, distMax, distStep)
    println(@sprintf("Earth radius is %f", earthRadius))
    xMax = Int64(trunc( (angleMax-angleMin)/angleStep ))+1
    outWidth  = xMax+1
    outHeight = size(range(vertAngleMin, vertAngleMax, step=angleStep))[1]+1
    println(@sprintf("Output size is %d x %d pixels", outWidth, outHeight))
    println(@sprintf("Output resolution is %f mrad per pixel or %f pixels per degree", angleStep * 1000.0, 1/toDegrees(angleStep)))
    output = zeros(UInt16, outHeight, outWidth)
    distances = range(0, distMax, step=distStep)
    for x in 0:xMax
        azimuth = angleMin + x*angleStep
        cosAz   = cos(azimuth)
        sinAz   = sin(azimuth)

        vertAngle     = vertAngleMin
        h0            = eye.height
        rayCastHeight = h0
        index         = 1
        for dist in distances
            point = pRef + vNorth*dist*cosAz + vEast*dist*sinAz;
            llh = xyz_to_llh(wgs84, PositionXYZ(point[1], point[2], point[3]))
            rayCastHeight = h0 + sin(vertAngle) * dist
            terrainHeight = earthCurve[index] + getHeight(srtmRange, heightMap, llh.lat, llh.lon)
            #println(@sprintf("dist=%f cast=%f, terr=%f, va=%f", dist, rayCastHeight, terrainHeight, vertAngle ))
            #sleep(0.1)
            if terrainHeight > rayCastHeight 
                newVertAngle = atan((terrainHeight-h0)/dist)
                yTop = Int64(trunc( (vertAngleMax-newVertAngle)/angleStep ))
                yBot = Int64(trunc( (vertAngleMax-   vertAngle)/angleStep ))
                v = Int64(trunc( dist / distStep )) 
                #println(@sprintf("x=%d, yTop=%d, yBot=%d, val=%d, newAngle=%f", x, yTop, yBot, v, newVertAngle))
                #sleep(0.1)
                for y in yTop:yBot 
                    output[y, x+1] = v
                end
                vertAngle = newVertAngle
            end
            index = index + 1
        end
    end
    return output
end

# Info (from https://www.udeuschle.de/panoramas/makepanoramas_en.htm)
# Lat: 50.08309 Lon 17.23094 Alt(auto+10m): 1500+10
# View direction: 112.5, extension 45 left: 90, right: 135, resolution 20pix/deg
# Tilt, range, vert. exaggeration 1.2


function main()
    tileDir   = "d:\\_disk_d_old\\devel-python\\panorama\\data_srtm"
    srtmRange = SrtmRange(47, 15, 50, 21)
    #eyePos    = PositionLLH(49.5460, 18.448, 1330.0)
    eyePos = PositionLLH(50.08309, 17.23094, 1510)

    heightMap=loadData(srtmRange, tileDir)
    #saveHeightMap(data)
    output = makeDistMap(eyePos, srtmRange, heightMap)
   
    minValue = minimum(output)
    maxValue = maximum(output)  
    println(@sprintf("min=%d max=%d", minValue, maxValue))  
    norm = Gray.(output/maxValue)
    println("Saving distmap-gray.png")
    save("disttmap-gray.png", norm)     
end

main()
