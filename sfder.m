function dNdxi = sfder(nen,nsd,xi)
dNdxi = zeros(nen,nsd);
if nsd == 1
    if nen == 2
        dNdxi = [0.5;-0.5];
    end
elseif nsd == 2
    if nen == 3
        dNdxi = [1,0;0,1;-1,-1];
    elseif nen == 4
        dNdxi(1,1) = -0.25*(1.-xi(2));
        dNdxi(1,2) = -0.25*(1.-xi(1));
        dNdxi(2,1) = 0.25*(1.-xi(2));
        dNdxi(2,2) = -0.25*(1.+xi(1));
        dNdxi(3,1) = 0.25*(1.+xi(2));
        dNdxi(3,2) = 0.25*(1.+xi(1));
        dNdxi(4,1) = -0.25*(1.+xi(2));
        dNdxi(4,2) = 0.25*(1.-xi(1));
    end
end
end
