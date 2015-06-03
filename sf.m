function N = sf(nen,nsd,xi)
N = zeros(nen, 1);
% 1D
if nsd == 1 && nen == 2
    N(1) = 0.5*(1.+xi(1));
    N(2) = 0.5*(1.-xi(1));
% 2D rectangle
elseif nsd == 2 && nen == 4
    N(1) = 0.25*(1.-xi(1))*(1.-xi(2));
    N(2) = 0.25*(1.+xi(1))*(1.-xi(2));
    N(3) = 0.25*(1.+xi(1))*(1.+xi(2));
    N(4) = 0.25*(1.-xi(1))*(1.+xi(2));
end
end