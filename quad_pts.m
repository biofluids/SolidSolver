function xi = quad_pts(nsd,nen,npt)
xi = zeros(nsd,npt);
if nsd == 1 % 1D
    if npt == 2
        xi(1,1) = -1./sqrt(3);
        xi(1,2) = -xi(1,1);
    elseif npt == 1
        xi(1,1) = 0.;       
    end
elseif nsd == 2 % 2D
    if nen == 3 % triangle
        if npt == 1
            xi = [1/3.;1/3.];
        end
    elseif nen == 4 % rectangle
        if npt == 4
            xi = [-0.5773502692,0.5773502692,-0.5773502692,0.5773502692;...
                -0.5773502692,-0.5773502692,0.5773502692,0.5773502692];
        elseif npt == 1
            xi = [0.;0.];
        end
    end
end
end