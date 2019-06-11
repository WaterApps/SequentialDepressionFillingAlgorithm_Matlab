function[volume] = getVolumeNonAdjacentParents(areaCellCount, pitCell, spilloverElevation, dem, flow_direction_parents, cellsize)
    indicesToCheck = zeros(areaCellCount, 1);
    indicesToCheck(1) = pitCell;
    j = 1;
    i = 1;
    volume = double(0);
    while i <= j
        volume = volume + ((spilloverElevation - double(dem(indicesToCheck(i))))*cellsize*cellsize);
        parents = flow_direction_parents{indicesToCheck(i)};
        lessThans = parents(dem(parents) <= spilloverElevation);
        k = sum(length(lessThans));
        indicesToCheck(j+1:j+k) = lessThans;
        j = j + k;
        i = i + 1;
    end
end