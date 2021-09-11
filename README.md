# What's this

First ever Julia code attempt :-)

Panorama generator / voxel rendering engine that shows panorama as seen from certain location and annotates hills visible in the distance.

As of 2021-08-28 it's for personal use as it contains many hard-coded variables (heightmap path, view set to Praděd,CZ)

## Historical background

Original code for panorama rendering was written in Python and unfinished/unpolished.
I've found https://www.udeuschle.de/panoramas/makepanoramas_en.htm that was by far superior.
Program in Python run in order of tens of minutes (2.6GHz Core2Duo) and manual measuring of angle (encoded in x position) and distance (encoded in color) to find far away hills on the map was not convenient. It was enough to confirm that High Tatras are visible from Praděd and if Schneeberg is visible from few places around Brno. Nothing else.

But because I found Julia language accidentally (by finding JuliaMono font) and it stated that it's way faster than Python, I tried it and after three evenings extended it beyond original functionality. Hopefully program runs in order of seconds on historic i5-4590@3.3GHz using three cores.

Honestly, for every other Python script I wrote speed does not matter.

## TODOs

* [x] Automatic anotations (summits from database, azimuths, horizon line)
* [x] Change annotation rendering from Luxor to Cairo (it's more useful to learn)
* [x] Render annoations to the same image
* [ ] WIP: Summit database for Slovakia and Alps, some manual entries for local hills
* [x] Make it work over any azimuth range (e.g. 350° to 10°, 170° to -170°)
* [ ] See above and check/fix code that discards points of interest
* [ ] WIP: Optimizations (wgs84->sphere), remove atan
* [ ] Explore where WGS84 vs Sphere matters (cause annotation are few pixels off)
* [ ] WIP: Object oriented code (how classes work in Julia?)
* [ ] Coloring (and shading?), legend (color->distance)
* [ ] Input data interpolation
* [ ] Test for western/southern hemisphere? I don't need it.
* [ ] Sky, haze, clouds, athmospheric light scattering, stars, sun, texturing (just kidding)

## Interesting problems

* Earth is not flat - if you draw line from point A to B, it has different azimuth in both points. Imagine plane that starts in USA, goes north east and arrives to Europe going south east. This makes annotating distant hills somewhat difficult.
* Earth is not perfect sphere - hopefully this does not seem to matter as much.
* Athmospehere gets less dense with altitude due to decreasing pressure. As it gets colder, density slightly increases (for given pressure), which slightly compensates this effect. Refractive index changes a bit and ray of light that travels up at low angle is bend back towards the ground. By trial-error and comparing rendered image with photo, surrounding terrain is mapped onto sphere that is 18% bigger than Earth. Reference photo was taken at temperature inversion. It's not typical, but it helps to keep humid air, dust and smog close to the ground and it makes extreme visibility posible.
* Surprisingly not all data I've found are accurate. Czech hills with prominence over 100m, Slovak hills over 200m have sometimes position erros up to low hundreds of meters. Original SRTM data have some voids as optical measurement failed on snow covered places and this is another source of misplaced annotations (hopefully fixed SRTM data can be found and there are other data sources)

## Interesting links

* Below the horizon—the physics of extreme visual ranges (Michael Vollmer): https://www.osapublishing.org/ao/abstract.cfm?uri=ao-59-21-f11
* Online tool with very same idea https://www.udeuschle.de/panoramas/makepanoramas_en.htm
* Application with very same idea https://www.peakfinder.org
* Useful calculations on sphere https://www.movable-type.co.uk/scripts/latlong.html
* Useful calculations on spheroid https://www.movable-type.co.uk/scripts/latlong-vincenty.html
* Terrain rendering algorithm in less than 20 lines of code https://github.com/s-macke/VoxelSpace (when I wrote Python version 10 years ago I had that game and Mars demo/Time Clarke/1994 in mind)

## Sample output

```
Requesting data for area 47°N 15°E - 50°N 21°E ... (aprox. 530x444 km).
I will read 7x4=28 tiles, heightmap size is 8401x4801 (81 MB).
Loading tile 002/028 lat=47, lon=19
Loading tile 001/028 lat=47, lon=17
Loading tile 024/028 lat=50, lon=20
  ⋮ 
Loading tile 025/028 lat=50, lon=16
Loading tile 026/028 lat=50, lon=18
Earth radius is 6367.1 km (diffraction x1.18)
Output size is 7855 x 901 pixels
Output resolution is 0.100000 mrad per pixel or 174.532925 pixels per degree
HILL TEST CODE ---- 
           Lysá hora is possibly visible at azimuth 123.8, distance 106.0 km
                Smrk is possibly visible at azimuth 127.5, distance 104.0 km
             Kněhyně is possibly visible at azimuth 129.5, distance 101.6 km
              Travný is possibly visible at azimuth 121.8, distance 108.6 km
              Skalka is possibly visible at azimuth 127.0, distance  97.0 km
      Velký Javorník is possibly visible at azimuth 132.4, distance  91.1 km
  ⋮         
min=0 max=2376
Saving distmap-gray.png
Saving horizon.png
```

## Images and comparison

![](pano-20210828.png)

## Data sources

Heightmaps:
SRTM (Shuttle Radar Topography Mission) - no longer available online?
LiDar data for SK,AT https://data.opendataportal.at/dataset/dtm-europe

POIs:
Czech hills with prominence above 100m https://www.ultratisicovky.cz/products/ultrakopec/
Slovak hills with prominence above 200m http://www.peaklist.org/WWlists/euro600/slovakia/Slovakia200m.html
Open Street Maps (via Overpass API)


## Thoughts about Julia

After writing like 600 lines of code in 10 evenings. Highly subjective and I may be wrong. 
TL;DR: If you are not limited by speed in Python, stick with Python.

Pros:
* It seems really fast - compared to Python
* It's easier to handle dependencies than let's say C++ & CMake on Windows
* Many pros are shared with Python
  * Packages are easy to install and use
  * REPL as command line interface
  * Jupyter notebooks
* Basics are not that hard to learn and it's easy to set-up in VS Code
* Good tools for profiling/benchmarking (do they even exist in Python?)

Cons:
* Python is more widespread and solutions are much easier to find
* For 95% of code I wrote in Python, speed does not matter at all
* Performance is sometimes unpredicable and depends on implementation nuances
  * vector addition and multiplication may allocate memory and be very slow
  * threaded-for may allocate memory when single threaded implementation works fine and it ends up being slower
  * in general memory allocation may occur at any place where one would expect using just stack in C++
* Some features are "weird" e.g. 
  * no const keyword
  * classes are immutable by default on the other hand
  * class methods are more like function, they are not even part of class declaration
  * class can't be redefined, once you touch it (e.g. add member or change constructor), you have to restart REPL
* It's not clear how to structure larger code into files and how to work with includes/modules
* Some operators and constants are UTF-8 character. Good luck copy-pasting π, integer division operator, ...
* Very high run time when code is compiled or external library used for the first time - starting new REPL is pain
* Very slow debug (issue with running new REPL and evaluating variables takes long)
* Sometimes poor error reporting