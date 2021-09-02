#!/usr/bin/bash
#wget -O test.osm "https://overpass-api.de/api/interpreter?data=node[natural=peak](49.0,16.0,50.0,17.0);out;"
#curl --globoff -o test.osm.json "https://overpass-api.de/api/interpreter?data=[out:json];node[natural=peak](49.0,16.0,50.0,17.0);out;"
curl --globoff -o osm-cz-sk.json "https://overpass-api.de/api/interpreter?data=[out:json];node[natural=peak](47.5,12.0,51.2,22.6);out;"
