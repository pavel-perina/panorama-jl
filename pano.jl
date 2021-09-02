using Printf, LinearAlgebra
import FixedPointNumbers: N0f8
import Images
import CSV, DataFrames
import Cairo

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


# https://www.movable-type.co.uk/scripts/latlong.html
# for SphericalEarth, but best source found
function bearing(p1::PositionLLH, p2::PositionLLH)::Float64
    φ1 = toRadians(p1.lat)
    φ2 = toRadians(p2.lat)
    Δλ = toRadians(p2.lon-p1.lon)
    y  = sin(Δλ)*cos(φ2)
    x  = cos(φ1)*sin(φ2) - sin(φ1)*cos(φ2)*cos(Δλ)
    θ  = atan(y, x)
    return mod(toDegrees(θ) + 360.0, 360.0)
end


xyz_to_vector(xyz::PositionXYZ)::Vector{Float64} = [xyz.x; xyz.y; xyz.z]



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


function loadData(range::LatLonRange, tileDir)
    getHgtFileName(lat, lon) = @sprintf "N%02dE%03d.hgt" lat lon
    getHgtFilePath(lat, lon, tileDir) = @sprintf "%s/%s" tileDir getHgtFileName(lat, lon)

    nTilesHoriz = range.maxLon - range.minLon + 1
    nTilesVert  = range.maxLat - range.minLat + 1
    nTilesTotal = nTilesHoriz * nTilesVert
    dataWidth   = nTilesHoriz*1200+1
    dataHeight  = nTilesVert*1200+1

    @printf("Requesting data for area %d°N %d°E - %d°N %d°E ... (aprox. %3.0fx%3.0f km).\n",
        range.minLat, range.minLon,
        range.maxLat, range.maxLon,
        (nTilesHoriz)*111.1*cos(range.minLat/180.0*π),
        (nTilesVert)*111.1
        )
    @printf("I will read %dx%d=%d tiles, heightmap size is %dx%d (%d MB).\n",
        nTilesHoriz, nTilesVert, nTilesTotal,
        dataWidth, dataHeight, dataWidth*dataHeight*2/1000000
        )

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
            @printf("Loading tile %03d/%03d lat=%02d, lon=%02d\n", progress, nTilesTotal, lat, lon)
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
        distStep = 50.0
        pRef   = xyz_to_vector(llh_to_xyz(ellipsoid, eye))
        vZ     = [0.0; 0.0; 1.0]
        vUp    = normalize(pRef)
        vEast  = normalize(cross(-vUp,vZ))
        vNorth = normalize(cross(vEast,-vUp))

        azimuthDiff = abs(azimuthMinR-azimuthMaxR)
        if (azimuthDiff > 2.0*π)
            throw(DomainError(azimuthDiff, "Azimuth difference should not exceed 2π"))
        end
        if (azimuthMinR > azimuthMaxR)
            azimuthMinR = azimuthMinR - 2.0*π
        end
        if (azimuthMinR < 0.0 && azimuthMaxR < 0.0)
            azimuthMinR = azimuthMinR + 2.0*π
            azimuthMaxR = azimuthMaxR + 2.0*π
        end
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
    @printf("Earth radius is %6.1f km (refraction x%4.2f)\n", earthRadius/vp.refractionCoef/1000.0, vp.refractionCoef)
    @printf("Output size is %d x %d pixels\n", vp.outWidth, vp.outHeight)
    @printf("Output resolution is %f mrad per pixel or %f pixels per degree\n", vp.angleStep * 1000.0, 1.0/toDegrees(vp.angleStep))
    output      = zeros(UInt16, vp.outHeight, vp.outWidth)
    distances   = range(0.0, vp.distMax, step=vp.distStep)
    Threads.@threads for x in 0:(vp.outWidth-1)
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
end


function extractOutlines(distMap::Matrix{UInt16})
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


function testPixel(distMap::Matrix{UInt16}, x::UInt64, y::UInt64, radius::Int64, value::UInt16, valueTolerance::UInt16)::Bool
    if (x < radius+1) || (y < radius+1) || (x+radius> size(distMap)[2]) || (y+radius> size(distMap)[1])
        return false
    end
    for row in (y-radius):(y+radius)
        for col in (x-radius):(x+radius)
            mapValue = distMap[row, col]
            if (mapValue > value) && ((mapValue - value) ≤ valueTolerance)
                return true
            end
            if (mapValue < value) && ((value - mapValue) ≤ valueTolerance)
                return true
            end
        end
    end
    return false
end


function drawSummits(vp::ViewPort, distMap::Matrix{UInt16})
    println("Loading data")
    dfFiltered = DataFrames.DataFrame(Summit = String[], Elevation = Float64[], Distance=Float64[], X=UInt64[], Y=UInt64[])

    # TODO: options
    #hillsCZ = CSV.File("data-cz-prom100.tsv") |> DataFrames.DataFrame
    #hillsSK = CSV.File("data-sk-prom200.tsv") |> DataFrames.DataFrame
    #hills = vcat(hillsCZ, hillsSK)
    hills =  CSV.File("osm-cz-sk.tsv") |> DataFrames.DataFrame
    # convert to ours azimuth, angle above horizon and distance - project into 
    hill_to_xyz(ellipsoid::Ellipsoid, dfRow)::PositionXYZ = llh_to_xyz(ellipsoid, PositionLLH(dfRow["Latitude"], dfRow["Longitude"], dfRow["Elevation"])) 
    # difference between true and seen earth curvature
    elevationDropAtDistance(distance::Float64, radius::Float64)::Float64 = sqrt(radius*radius-distance*distance)-radius
    elevationDropCompensation(distance::Float64, radius::Float64, refractionCoef::Float64)::Float64 = elevationDropAtDistance(distance, radius*refractionCoef) - elevationDropAtDistance(distance, radius)
    # TODO: make distance function
    mLocalToWorld = hcat(vp.vEast,vp.vNorth,vp.vUp)
    mWorldToLocal = inv(mLocalToWorld)
    #ground = PositionLLH(vp.eye.lat, vp.eye.lon, 0.0)
    pRef = xyz_to_vector(llh_to_xyz(vp.ellipsoid, vp.eye))
    earthRadius = sqrt(dot(pRef, pRef))
    for hill in eachrow(hills)
        hill_world = hill_to_xyz(vp.ellipsoid, hill)
        hill_local_xyz = mWorldToLocal * [hill_world.x; hill_world.y; hill_world.z]
#        hill_local_xyz[3] = hill_local_xyz[3] - earthRadius  <---- this gives weird altitutude, let's restore it in the same way as raytracer
        hill_local_xyz[3] = 0.0
        distance = sqrt(dot(hill_local_xyz, hill_local_xyz))
        if distance > vp.distMax
            continue
        end
        hill_local_xyz[3] = hill["Elevation"] + elevationDropAtDistance(distance, earthRadius * vp.refractionCoef) - vp.eye.height
        azimuth = bearing(vp.eye, PositionLLH(hill["Latitude"], hill["Longitude"], 0.0))
        if (azimuth < toDegrees(vp.angleMin) || azimuth > toDegrees(vp.angleMax)) # FIXME: test for weird angles
            continue
        end
        @printf("%20s is possibly visible at azimuth %6.2f, distance %6.2f km", hill["Summit"], azimuth, distance/1000.0)
        elevationAngle=atan(hill_local_xyz[3], distance)
#        @printf("              hill=%f+%f, dist=%f\n", hill_local_xyz[3], elevationDropCompensation(distance, earthRadius, vp.refractionCoef), distance)
#        @printf("              , pixel.x,y=%5.0f,%5.0f\n", (toRadians(azimuth)-vp.angleMin)/vp.angleStep , (vp.vertAngleMax-elevationAngle)/vp.angleStep )
        testX::UInt64 = round((toRadians(azimuth)-vp.angleMin)/vp.angleStep)
        testY::UInt64 = round((vp.vertAngleMax-elevationAngle)/vp.angleStep)
        visible::Bool = testPixel(distMap, testX, testY, 4, UInt16(trunc(distance/vp.distStep)), UInt16(5))
        @printf(", visible=%s\n", visible)
        if visible
            push!(dfFiltered, (hill["Summit"], hill["Elevation"], distance, testX, testY))
        end
    end
    println("Drawing ...")
    # Initialize
    Cairo_set_line_color(ctx::Cairo.CairoContext) = Cairo.set_source_rgb(ctx, 131/255, 148/255, 150/255)
    Cairo_set_text_color(ctx::Cairo.CairoContext) = Cairo.set_source_rgb(ctx,  38/255, 139/255, 210/255)
    function Cairo_line(ctx::Cairo.CairoContext, x1::Core.Real, y1::Core.Real, x2::Core.Real, y2::Core.Real)
        Cairo.move_to(ctx, x1, y1)
        Cairo.line_to(ctx, x2, y2)
        Cairo.stroke(ctx)
    end
    surf = Cairo.CairoARGBSurface(vp.outWidth, vp.outHeight)
    ctx  = Cairo.CairoContext(surf)
    # Background (previous image)
    bgSurf = Cairo.read_from_png("outlines.png")
    Cairo.set_source_surface(ctx, bgSurf, 0.0, 0.0)
    Cairo.paint(ctx)
    # Annotations
    Cairo.select_font_face(ctx, "Fira Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    Cairo.set_font_size(ctx, 18.0)
    Cairo.set_line_width(ctx, 1.0)
    for poi in eachrow(dfFiltered)
        Cairo_set_line_color(ctx)
        Cairo_line(ctx, poi["X"]+0.5, poi["Y"], poi["X"]+0.5, 300.0)
        Cairo_set_text_color(ctx)
        Cairo.move_to(ctx, poi["X"]+5, 300.0-5.0)
        Cairo.save(ctx)
        Cairo.rotate(ctx, toRadians(-45.0))
        Cairo.show_text(ctx, poi["Summit"])
        Cairo.show_text(ctx, @sprintf(" (%3.0f km)", poi["Distance"]/1000.0))
        Cairo.restore(ctx)
    end
    # Azimuth ticks
    azMinD::Int = Int(ceil(toDegrees(vp.angleMin)))
    azMaxD::Int = Int(floor(toDegrees(vp.angleMax)))
    for az in azMinD:azMaxD
        x = round((toRadians(Float64(az))-vp.angleMin)/vp.angleStep)+0.5
        Cairo_set_line_color(ctx)
        Cairo_line(ctx, x, 38, x, 42)
        Cairo_line(ctx, x, 63, x, 68)
        Cairo_set_text_color(ctx)
        ext = Cairo.text_extents(ctx, @sprintf("%d", az))
        Cairo.move_to(ctx, x - ext[3]/2, 58)
        Cairo.show_text(ctx, @sprintf("%d °", az))
    end
    # Horizon line
    horizY = round(vp.vertAngleMax/vp.angleStep)+0.5
    Cairo.set_source_rgb(ctx, 238/255, 232/255, 213/255)
    Cairo_line(ctx, 0.0, horizY, vp.outWidth, horizY)
    # Finalize
    Cairo.write_to_png(surf, "outline-with-annotations.png" )
end

# Info (from https://www.udeuschle.de/panoramas/makepanoramas_en.htm)
# Lat: 50.08309 Lon 17.23094 Alt(auto+10m): 1500+10
# View direction: 112.5, extension 45 left: 90, right: 135, resolution 20pix/deg
# Tilt, range, vert. exaggeration 1.2

# TODO: determine lat/long range automatically
function main()
    tileDir     = "d:/_disk_d_old/devel-python/panorama/data_srtm"
	#tileDir     = "data_srtm"
    latLonRange = LatLonRange(47, 15, 50, 21)
    eye         = PositionLLH(50.08309, 17.23094, 1510)

    heightMap = loadData(latLonRange, tileDir)
    #saveHeightMap(data)
    ellipsoid = SphericalEarth()
    #ellipsoid = Wgs84()
    vp = ViewPort(ellipsoid, eye, toRadians(90.0), toRadians(135.0), -0.0560, 0.0339, 0.0001, 250.0e3, 1.18)
    distMap   = makeDistMap(vp, latLonRange, heightMap)

    minValue = minimum(distMap)
    maxValue = maximum(distMap)
    @printf("min=%d max=%d\n", minValue, maxValue)
    println("Saving distmap-gray.png")
    Images.save("distmap-gray.png", Images.Gray.(distMap/maxValue))

    println("Extracting outlines ...")
    outlines=extractOutlines(distMap)
    println("Saving outlines.png")
    Images.save("outlines.png", outlines)

    println("Creating annotations")
    drawSummits(vp, distMap)
    println("All done")

    # This can be fun: https://wiki.flightgear.org/Atmospheric_light_scattering
    # http://www.science-and-fiction.org/rendering/als.html
end

main()

