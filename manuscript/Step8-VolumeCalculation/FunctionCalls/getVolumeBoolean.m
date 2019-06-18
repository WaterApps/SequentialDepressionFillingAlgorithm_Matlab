function[volume] = getVolumeBoolean(spilloverElevation, pitId, dem, pits, cellsize)
    volume = sum((spilloverElevation - double(dem(pits == pitId & dem <= spilloverElevation))).*cellsize*cellsize);
end