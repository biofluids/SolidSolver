function npt = no_quad_pts(nsd,nen,reduced)
if reduced
    if nsd == 1
        if nen == 2
            npt = 2;
        end
    elseif nsd == 2
        if nen == 3
            npt = 1;
        elseif nen == 4
            npt = 1;
        end
    end
else
    if nsd == 1
        if nen == 2
            npt = 2;
        end
    elseif nsd == 2
        if nen == 3
            npt = 1;
        elseif nen == 4
            npt = 4;
        end
    end
end
end
