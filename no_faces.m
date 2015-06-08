function n=no_faces(nsd,nen)
if nsd==2
    if nen==3
        n=3;
    elseif nen==4
        n=4;
    end
elseif nsd==3
    if nen==4
        n=4;
    elseif nen==8
        n=6;
    end
end
end