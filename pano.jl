using Printf, LinearAlgebra
import FixedPointNumbers: N0f8
import Images
import CSV, DataFrames

toRadians(α::Float64) = α / 180.0 * π
toDegrees(α::Float64) = α * 180.0 / π

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


abstract type Ellipsoid
end


struct Wgs84 <: Ellipsoid
    a::Float64
    b::Float64
    f::Float64
    e::Float64
    e²::Float64
    function Wgs84()
        a  = 6378137.0
        f  = 1.0/298.257223563
        b  = a*(1.0-f)
        e² = 1.0 - b*b/(a*a) # or f * (2.0-f)
        e  = sqrt(e²)
        new(a,b,f,e,e²)
    end
end


struct SphericalEarth <: Ellipsoid
    r::Float64
    SphericalEarth() = new(6378137.0)
end


function llh_to_xyz(ellipsoid::Wgs84, llh::PositionLLH)::PositionXYZ
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
function xyz_to_llh(ellipsoid::Wgs84, xyz::PositionXYZ)::PositionLLH
    e²       = ellipsoid.e²
    p        = sqrt(xyz.x*xyz.x + xyz.y*xyz.y)
    lon      = atan(xyz.y / xyz.x)
    lat_init = atan(xyz.z/(p*(1.0 - e²)))
    v        = ellipsoid.a / sqrt(1.0-e²*sin(lat_init)*sin(lat_init))
    lat      = atan((xyz.z + e²*v*sin(lat_init))/p)
    height   = (p/cos(lat))-v
    return PositionLLH(toDegrees(lat), toDegrees(lon), height)
end


function llh_to_xyz(ellipsioid::SphericalEarth, llh::PositionLLH)::PositionXYZ
    Φ    = toRadians(llh.lat)
    λ    = toRadians(llh.lon)
    h    = llh.height
    cosΦ = cos(Φ)
    v    = ellipsioid.r + h
    x    = v * cosΦ * cos(λ)
    y    = v * cosΦ * sin(λ)
    z    = v * sin(Φ)
    return PositionXYZ(x,y,z)
end


function xyz_to_llh(ellipsoid::SphericalEarth, xyz::PositionXYZ)::PositionLLH
    v      = sqrt(xyz.x*xyz.x + xyz.y*xyz.y + xyz.z*xyz.z) 
    height = v - ellipsoid.r
    lon    = atan(xyz.y, xyz.x)    
    lat    = asin(xyz.z/v)
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
    progressLock = Threads.SpinLock()
    Threads.@threads for i in 0:(nTilesTotal-1)
        lat = range.minLat + div(i, nTilesHoriz)
        lon = range.minLon + mod(i, nTilesHoriz)
        # Print progress
        lock(progressLock) do
            progress = progress + 1
            println(@sprintf("Loading tile %03d/%03d lat=%02d, lon=%02d", progress, nTilesTotal, lat, lon))
        end
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
    return data
end


function saveHeightMap(data)
#  minValue = minimum(data)
   maxValue = maximum(data)
   norm = Images.Gray.(data/maxValue)
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

struct ViewPort
    ellipsoid::Ellipsoid
    eye::PositionLLH

    angleMin::Float64
    angleMax::Float64
    vertAngleMin::Float64
    vertAngleMax::Float64
    angleStep::Float64

    distMax::Float64
    distStep::Float64

    refractionCoef::Float64

    vUp::Vector{Float64}
    vNorth::Vector{Float64}
    vEast::Vector{Float64}

    outWidth::Int
    outHeight::Int

    function ViewPort(ellipsoid::Ellipsoid, eye::PositionLLH, azimuthMinR, azimuthMaxR, elevationMinR, elevationMaxR, angularStepR, distMaxM, refractionCoef)
        distStep = 100.0
        eyeXYZ = llh_to_xyz(ellipsoid, eye)
        pRef   = [eyeXYZ.x; eyeXYZ.y; eyeXYZ.z]
        vZ     = [0.0; 0.0; 1.0]
        vUp    = normalize(pRef)
        vEast  = normalize(cross(-vUp,vZ))
        vNorth = normalize(cross(vEast,-vUp))

        xMax        = Int64(trunc( (azimuthMaxR-azimuthMinR)/angularStepR ))+1
        outWidth    = xMax+1
        outHeight   = size(range(elevationMinR, elevationMaxR, step=angularStepR))[1]+1
    
        new(ellipsoid, eye,
            azimuthMinR, azimuthMaxR, elevationMinR, elevationMaxR, angularStepR, 
            distMaxM, distStep, 
            refractionCoef,
            vUp, vNorth, vEast,
            outWidth, outHeight)
    end
end


function eyeVec(vp::ViewPort)::Vector{Float64}
    eyeXYZ = llh_to_xyz(vp.ellipsoid, vp.eye)
    return [eyeXYZ.x; eyeXYZ.y; eyeXYZ.z]
end


function makeDistMap(vp::ViewPort, latLonRange::LatLonRange, heightMap::Matrix{UInt16})::Matrix{UInt16}
    pRef = eyeVec(vp)
    earthRadius = sqrt(dot(pRef, pRef)) * vp.refractionCoef
    earthCurve  = makeEarthCurve(earthRadius, vp.distMax, vp.distStep)
    println(@sprintf("Earth radius is %6.1f km (refraction x%4.2f)", earthRadius/vp.refractionCoef/1000.0, vp.refractionCoef))
    println(@sprintf("Output size is %d x %d pixels", vp.outWidth, vp.outHeight))
    println(@sprintf("Output resolution is %f mrad per pixel or %f pixels per degree", vp.angleStep * 1000.0, 1.0/toDegrees(vp.angleStep)))
    output      = zeros(UInt16, vp.outHeight, vp.outWidth)
    distances   = range(0.0, vp.distMax, step=vp.distStep)
    Threads.@threads for x in 0:(vp.outWidth-1)
        #println(@sprintf("Rendering line %d of %d", x, xMax))
        azimuth = vp.angleMin + x*vp.angleStep
        cosAz   = cos(azimuth)
        sinAz   = sin(azimuth)

        vertAngle     = vp.vertAngleMin
        h0            = vp.eye.height
        rayCastHeight = h0
        index         = 1
        direction     = vp.vNorth*cosAz + vp.vEast*sinAz
        point = [0.0,0.0,0.0]
        for dist in distances
            #point = pRef + dist * direction;
            #point = [pRef[1]+dist*direction[1], pRef[2]+dist*direction[2], pRef[3]+dist*direction[3]]
            point[1] = pRef[1]+dist*direction[1]
            point[2] = pRef[2]+dist*direction[2]
            point[3] = pRef[3]+dist*direction[3]
            llh = xyz_to_llh(vp.ellipsoid, PositionXYZ(point[1], point[2], point[3]))
            rayCastHeight = h0 + sin(vertAngle) * dist
            terrainHeight = earthCurve[index] + getHeight(latLonRange, heightMap, llh.lat, llh.lon)
            if terrainHeight > rayCastHeight 
                newVertAngle = atan((terrainHeight-h0)/dist)
                yTop = Int64(trunc( (vp.vertAngleMax-newVertAngle)/vp.angleStep ))
                yBot = Int64(trunc( (vp.vertAngleMax-   vertAngle)/vp.angleStep ))
                v = UInt16(trunc( dist / vp.distStep ))
                for y in yTop:yBot 
                    output[y, x+1] = v
                end
                vertAngle = newVertAngle
            end
            index = index + 1
        end
    end
    return output
    println("SUMMIT TEST CODE ---- ")
    # TEST code
    hills = CSV.File("data-cz-prom100.tsv") |> DataFrames.DataFrame
    # convert to ours azimuth, angle above horizon and distance - project into 
    hill_to_xyz(ellipsoid::Ellipsoid, dfRow)::PositionXYZ = llh_to_xyz(ellipsoid, PositionLLH(dfRow["Latitude"], dfRow["Longitude"], dfRow["Elevation"])) 
    mLocalToWorld = hcat(vEast,vNorth,-vDown)
    mWorldToLocal = inv(mLocalToWorld)
    for hill in eachrow(hills)
        hill_world = hill_to_xyz(ellipsoid, hill)
        hill_local_xyz = mWorldToLocal * [hill_world.x, hill_world.y, hill_world.z]
        hill_local_xyz[3] = hill_local_xyz[3] - sqrt(dot(pRef, pRef))
        distance = sqrt(dot(hill_local_xyz, hill_local_xyz))
        azimuth  = toDegrees(atan(hill_local_xyz[1], hill_local_xyz[2])) # goes from north to east so y/x is swapped
        #print(hill["Summit"])
        if azimuth < 0
            azimuth = azimuth + 360
        end
        if (distance > distMax || azimuth < toDegrees(angleMin) || azimuth > toDegrees(angleMax))
            continue
        end
        print(@sprintf("%20s is possibly visible at azimuth %5.1f, distance %5.1f km", hill["Summit"], azimuth, distance/1000.0))
        # FIXME: do correction for refraction (?)
        elevationAngle=atan(distance, hill_local_xyz[3])
        println(@sprintf("hill=%f, dist=%f", hill_local_xyz[3], distance))
        println(@sprintf(", pixel.x,y=%5.0f,%5.0f", (toRadians(azimuth)-angleMin)/angleStep , (vertAngleMax-elevationAngle)/angleStep ))
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

function drawHorizon()
end

function drawAzimuths()
end

function testSummit(distMap::Matrix{UInt16})
end

function drawSummits()
end

# Info (from https://www.udeuschle.de/panoramas/makepanoramas_en.htm)
# Lat: 50.08309 Lon 17.23094 Alt(auto+10m): 1500+10
# View direction: 112.5, extension 45 left: 90, right: 135, resolution 20pix/deg
# Tilt, range, vert. exaggeration 1.2

# TODO: determine lat/long range automatically
function main()
    tileDir     = "d:\\_disk_d_old\\devel-python\\panorama\\data_srtm"
    latLonRange = LatLonRange(47, 15, 50, 21)
    eye         = PositionLLH(50.08309, 17.23094, 1510)

    heightMap = loadData(latLonRange, tileDir)
    #saveHeightMap(data)
    ellipsoid = SphericalEarth()
    vp = ViewPort(ellipsoid, eye, toRadians( 90.0), toRadians(135.0), -0.0560, 0.0339, 0.0001, 250.0e3, 1.18)
    distMap   = makeDistMap(vp, latLonRange, heightMap)

    minValue = minimum(distMap)
    maxValue = maximum(distMap)
    println(@sprintf("min=%d max=%d", minValue, maxValue))
    println("Saving distmap-gray.png")
    Images.save("distmap-gray.png", Images.Gray.(distMap/maxValue))

    println("Saving horizon.png....")
    hrz=extractHorizon(distMap)
    Images.save("horizon.png", hrz)

    # This can be fun: https://wiki.flightgear.org/Atmospheric_light_scattering
    # http://www.science-and-fiction.org/rendering/als.html

end

main()
