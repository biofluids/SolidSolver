function npt = no_quad_pts(nsd,nen,reduced)
if reduced
    npt = 1;
else
if nsd == 1 && nen == 2
    npt = 2;
elseif nsd == 2 && nen == 4
    npt = 4;
end
end
end
