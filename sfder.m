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
elseif nsd==3
    if nen==4
        dNdxi(1,1) = 1.;
        dNdxi(2,2) = 1.;
        dNdxi(3,3) = 1.;
        dNdxi(4,1) = -1.;
        dNdxi(4,2) = -1.;
        dNdxi(4,3) = -1.;
    elseif nen==8
        dNdxi(1,1) = -(1.-xi(2))*(1.-xi(3))/8.;
        dNdxi(1,2) = -(1.-xi(1))*(1.-xi(3))/8.;
        dNdxi(1,3) = -(1.-xi(1))*(1.-xi(2))/8.;
        dNdxi(2,1) = (1.-xi(2))*(1.-xi(3))/8.;
        dNdxi(2,2) = -(1.+xi(1))*(1.-xi(3))/8.;
        dNdxi(2,3) = -(1.+xi(1))*(1.-xi(2))/8.;
        dNdxi(3,1) = (1.+xi(2))*(1.-xi(3))/8.;
        dNdxi(3,2) = (1.+xi(1))*(1.-xi(3))/8.;
        dNdxi(3,3) = -(1.+xi(1))*(1.+xi(2))/8.;
        dNdxi(4,1) = -(1.+xi(2))*(1.-xi(3))/8.;
        dNdxi(4,2) = (1.-xi(1))*(1.-xi(3))/8.;
        dNdxi(4,3) = -(1.-xi(1))*(1.+xi(2))/8.;
        dNdxi(5,1) = -(1.-xi(2))*(1.+xi(3))/8.;
        dNdxi(5,2) = -(1.-xi(1))*(1.+xi(3))/8.;
        dNdxi(5,3) = (1.-xi(1))*(1.-xi(2))/8.;
        dNdxi(6,1) = (1.-xi(2))*(1.+xi(3))/8.;
        dNdxi(6,2) = -(1.+xi(1))*(1.+xi(3))/8.;
        dNdxi(6,3) = (1.+xi(1))*(1.-xi(2))/8.;
        dNdxi(7,1) = (1.+xi(2))*(1.+xi(3))/8.;
        dNdxi(7,2) = (1.+xi(1))*(1.+xi(3))/8.;
        dNdxi(7,3) = (1.+xi(1))*(1.+xi(2))/8.;
        dNdxi(8,1) = -(1.+xi(2))*(1.+xi(3))/8.;
        dNdxi(8,2) = (1.-xi(1))*(1.+xi(3))/8.;
        dNdxi(8,3) = (1.-xi(1))*(1.+xi(2))/8.;
    end
end
end
