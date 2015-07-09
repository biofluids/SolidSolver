function stress = Kirchhoffstress(ned,nsd,B,J,materialprops)
% tau = J*sigma
stress = zeros(ned,nsd);
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
% 2d strain, assuming the third direction is undeformable
if ned==2
    temp1 = temp1 + 1; 
    temp2 = temp2 + 1;
end
for i=1:ned
    for j=1:nsd
        stress(i,j) = mu1*(B(i,j) - temp1*dl(i,j)/3.)/J^(2./3.) ...
            + K1*J*(J-1)*dl(i,j);
        if materialprops(1)==3
            stress(i,j)=stress(i,j)+...
                mu2*(temp1*B(i,j)-temp1^2*dl(i,j)/3+temp2*dl(i,j)/3.)/J^(5/3.);
            for k=1:nsd
                stress(i,j)=stress(i,j) - mu2*B(i,k)*B(k,j)/J^(5/3.);
            end
        end
     end
end
end
