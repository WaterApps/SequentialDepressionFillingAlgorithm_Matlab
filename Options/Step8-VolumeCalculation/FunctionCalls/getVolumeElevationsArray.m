function[volume, evs] = getVolumeElevationsArray(elevations1, elevations2, spilloverElevation, cellsize)
    evs = [elevations1; elevations2];
    lessThans = evs(:, 1) <= spilloverElevation;
    volume = sum((spilloverElevation - evs(lessThans, 1)).*evs(lessThans, 2)*cellsize*cellsize);
end