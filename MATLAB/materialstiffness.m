function C = materialstiffness(ned,nsd,B,J,materialprops)
% C_{ijkl} is returned as a 4th order tensor
C = zeros(ned,nsd,ned,nsd);
mu1 = materialprops(2); 
mu2 = materialprops(3); 
K1 = materialprops(4);
dl = [[1,0,0];[0,1,0];[0,0,1]];
temp1 = trace(B);
temp2 = 0;
for i = 1:nsd
    for j = 1:nsd
        temp2 = temp2 + B(i,j)^2;
    end
end
if nsd==2
    temp1 = temp1+1;
    temp2 = temp2+1;
end
for i = 1:ned
    for j = 1:nsd
        for k = 1:ned
            for l = 1:nsd
                if materialprops(1)==2
                    C(i,j,k,l) = mu1*( dl(i,k)*B(j,l)+B(i,l)*dl(j,k)...
                        - (2/3)*(B(i,j)*dl(k,l)+dl(i,j)*B(k,l))...
                        + (2/3)*temp1*dl(i,j)*dl(k,l)/3 )/J^(2/3) ...
                        + K1*(2*J-1)*J*dl(i,j)*dl(k,l);
                elseif materialprops(1)==3
                    C(i,j,k,l) = mu1*( dl(i,k)*B(j,l)+B(i,l)*dl(j,k)...
                        - (2/3)*(B(i,j)*dl(k,l)+dl(i,j)*B(k,l))...
                        + (2/3)*temp1*dl(i,j)*dl(k,l)/3 )/J^(2/3) ...
                        + K1*(2*J-1)*J*dl(i,j)*dl(k,l)...
                        + mu2/J^(4/3.)*(2*B(i,j)*B(k,l)...
                        + temp1*(B(l,j)*dl(i,k)+B(l,i)*dl(j,k)-4/3.*(B(i,j)*dl(l,k)+B(l,k)*dl(i,j)))...
                        - B(k,j)*B(l,i)-B(k,i)*B(l,j)+4/9.*(temp1^2-temp2)*dl(i,j)*dl(l,k));
                    for p = 1:nsd
                        C(i,j,k,l) = C(i,j,k,l)+mu2/J^(4/3.)*(-B(l,p)*(B(p,j)*dl(i,k)+B(p,i)*dl(j,k))...
                            +4/3.*(B(l,p)*B(k,p)*dl(i,j)+B(i,p)*B(j,p)*dl(l,k)));
                    end
                end
            end
        end
    end
end
end
    