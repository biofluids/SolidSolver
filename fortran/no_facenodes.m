function n=no_facenodes(nsd,nen)
if nsd==2
    n=2;
elseif nsd==3
    if nen==4
        n=3;
    elseif nen==8
        n=4;
    end
end
end