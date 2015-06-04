function N = sf(nen,nsd,xi)
N = zeros(nen, 1);
if nsd == 1 
    if nen == 2
        N(1) = 0.5*(1.+xi(1));
        N(2) = 0.5*(1.-xi(1));
    end
elseif nsd == 2
    if nen == 3
        N(1) = xi(1);
        N(2) = xi(2);
        N(3) = 1.-xi(1)-xi(2);   
    elseif nen == 4
        N(1) = 0.25*(1.-xi(1))*(1.-xi(2));
        N(2) = 0.25*(1.+xi(1))*(1.-xi(2));
        N(3) = 0.25*(1.+xi(1))*(1.+xi(2));
        N(4) = 0.25*(1.-xi(1))*(1.+xi(2));
    end
end
end
        