function[volume] = getVolumeCellIndexes(cellIndexes, spilloverElevation, dem, cellsize)
    lessThans = cellIndexes(dem(cellIndexes) <= spilloverElevation);
    volume = sum((spilloverElevation - dem(lessThans))*cellsize*cellsize);
    lessThans = [];
end