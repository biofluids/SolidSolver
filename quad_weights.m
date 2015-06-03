function weights = quad_weights(nsd,nen,npt)
weights = zeros(npt,1);
if nsd == 1 % 1D
    if npt == 2
        weights = [1.,1.];
    end
elseif nsd == 2 % 2D
    if npt == 4
        weights = [1.,1.,1.,1.];
    elseif npt == 1
        weights = 4;
    end
end
end
    