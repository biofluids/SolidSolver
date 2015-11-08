module reduced_material
	implicit none
	! this module contains a function to compute the material stiffness tensor, and a function to compute Kirchhoff stress
	
	
contains
	
	function materialstiffness(nsd,ned,B,Ja,materialprops)
		implicit none
		
		integer, intent(in) :: nsd, ned
		real(8), intent(in) :: Ja
		real(8), dimension(nsd,nsd), intent(in) :: B
		real(8), dimension(5), intent(in) :: materialprops
		real(8), dimension(ned,nsd,ned,nsd):: materialstiffness
		
		integer :: i, j, k, l, p
		real(8) :: mu1, mu2, K1, temp1, temp2
		real(8), dimension(3,3) :: dl
		
		! initialize 
		dl = reshape([1.,0.,0.,0.,1.,0.,0.,0.,1.],shape(dl))
		do i=1,ned
			do j=1,nsd
				do k=1,ned
					do l=1,nsd
						materialstiffness(i,j,k,l) = 0.
					end do
				end do
			end do
		end do
		mu1 = materialprops(2)
		mu2 = materialprops(3)
		K1 = materialprops(4)
		temp1 = 0.
		temp2 = 0.
		
		do i=1,nsd
			temp1 = temp1 + B(i,i)
		end do
		
		do i=1,nsd
			do j=1,nsd
				temp2 = temp2 + B(i,j)**2
			end do 
		end do
		
		if (nsd == 2) then
			temp1 = temp1 + 1.
			temp2 = temp2 + 1.
		end if
		
		do i=1,ned
			do j=1,nsd
				do k=1,ned
					do l=1,nsd
						if (materialprops(1) == 2.) then
							materialstiffness(i,j,k,l) = mu1*( dl(i,k)*B(j,l)+B(i,l)*dl(j,k) &
							             - (2/3.)*(B(i,j)*dl(k,l)+dl(i,j)*B(k,l)) &
							             + (2/3.)*temp1*dl(i,j)*dl(k,l)/3.)/Ja**(2/3.) &
							             + K1*(2*Ja-1.)*Ja*dl(i,j)*dl(k,l);
						else if (materialprops(1) == 3.) then
						    materialstiffness(i,j,k,l) = mu1*( dl(i,k)*B(j,l)+B(i,l)*dl(j,k) &
						                 - (2/3.)*(B(i,j)*dl(k,l)+dl(i,j)*B(k,l)) &
						                 + (2/3.)*temp1*dl(i,j)*dl(k,l)/3. )/Ja**(2/3.) &
				                         + K1*(2*Ja-1.)*Ja*dl(i,j)*dl(k,l) &
		     	                         + mu2/Ja**(4/3.)*(2*B(i,j)*B(k,l) &
						                 + temp1*(B(l,j)*dl(i,k)+B(l,i)*dl(j,k)-4/3.*(B(i,j)*dl(l,k)+B(l,k)*dl(i,j))) &
						                 - B(k,j)*B(l,i)-B(k,i)*B(l,j)+4/9.*(temp1**2-temp2)*dl(i,j)*dl(l,k));
						    do p=1,nsd
								materialstiffness(i,j,k,l) = materialstiffness(i,j,k,l)+mu2/Ja**(4/3.)*(-B(l,p)*(B(p,j)*dl(i,k)+B(p,i)*dl(j,k)) &
							             	 + 4/3.*(B(l,p)*B(k,p)*dl(i,j)+B(i,p)*B(j,p)*dl(l,k)));
						    end do
						end if
					end do
				end do
			end do
		end do
	end function materialstiffness
	
	function Kirchhoffstress(nsd,ned,B,Ja,materialprops)
		implicit none
		
		integer, intent(in) :: nsd, ned
		real(8), intent(in) :: Ja
		real(8), dimension(nsd,nsd), intent(in) :: B
		real(8), dimension(5), intent(in) :: materialprops
		real(8), dimension(ned,nsd) :: Kirchhoffstress
		
		integer :: i, j, k
		real(8) :: mu1, mu2, K1, temp1, temp2
		real(8), dimension(3,3) :: dl
		
		! initialize 
		dl = reshape([1.,0.,0.,0.,1.,0.,0.,0.,1.],shape(dl))
		mu1 = materialprops(2)
		mu2 = materialprops(3)
		K1 = materialprops(4)
		temp1 = 0.
		temp2 = 0.	
		do i=1,nsd
			temp1 = temp1 + B(i,i)
		end do
		do i=1,ned
			do j=1,nsd
				Kirchhoffstress(i,j) = 0.
			end do
		end do
		do i=1,nsd
			do j=1,nsd
				temp2 = temp2 + B(i,j)**2
			end do 
		end do
		
		if (nsd == 2) then
			temp1 = temp1 + 1.
			temp2 = temp2 + 1.
		end if
		
		do i=1,ned
			do j=1,nsd
				Kirchhoffstress(i,j) = mu1*(B(i,j) - temp1*dl(i,j)/3.)/Ja**(2./3.) + K1*Ja*(Ja-1)*dl(i,j);
				if ( materialprops(1) == 3.) then
					Kirchhoffstress(i,j) = Kirchhoffstress(i,j)+mu2*(temp1*B(i,j)-temp1**2.*dl(i,j)/3.+temp2*dl(i,j)/3.)/Ja**(5/3.);
					do k=1,nsd
						Kirchhoffstress(i,j)=Kirchhoffstress(i,j) - mu2*B(i,k)*B(k,j)/Ja**(5./3.);
					end do
				end if
			end do
		end do
	end function Kirchhoffstress	
	
end module reduced_material