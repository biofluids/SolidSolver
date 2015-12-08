module mixed_material
	
	implicit none

contains
	function materialstiffness(nsd,ned,B,Ja,pressure,materialprops)
		! Compute the material stiffness tensor
		implicit none
	
		integer, intent(in) :: nsd, ned
		real(8), intent(in) :: Ja, pressure
		real(8), dimension(nsd,nsd), intent(in) :: B
		real(8), dimension(5), intent(in) :: materialprops
		real(8), dimension(ned,nsd,ned,nsd):: materialstiffness
	
		integer :: i, j, k, l, p
		real(8) :: mu1, mu2, K1, I1, I2, c1, c2, c3
		real(8), dimension(nsd,nsd) :: Bbar, BB, eye
	
		! square matrix
		do i=1,nsd
			do j=1,nsd
				if(i==j) then
					eye(i,j) = 1.
				else
					eye(i,j) = 0.
				end if
			end do
		end do
		materialstiffness = 0.
		
		! Change B to Bbar and compute Bbar squared
		Bbar = B/Ja**(2/3.)
		BB = matmul(Bbar,Bbar)
		! I1 I2 denote I1bar I2bar actually
		I1 = 0.
		do i=1,nsd
			I1 = I1 + Bbar(i,i)
		end do
	    I2 = I1**2
		do i=1,nsd
			do j=1,nsd
				I2 = I2 - Bbar(i,j)**2
			end do 
		end do
		if (nsd == 2) then
			I1 = I1 + 1.
			I2 = I2 - 1.
		end if
		I2 = I2/2.
		
		if (dAbs(materialprops(1)-4) < 1d-4) then
			c1 = materialprops(2)
			c2 = materialprops(3)
			c3 = materialprops(5)
			do i = 1, ned
				do j = 1, nsd
					do k = 1, ned
						do l = 1, nsd
							materialstiffness(i,j,k,l) = (8*c2 + 24*c3*(I1 - 3))*Bbar(i,j)*Bbar(k,l) &
													   - (4/3.*c1 + (16/3.*I1 - 8)*c2 + 12*(I1-1)*(I1-3)*c3)*(eye(i,j)*Bbar(k,l) + eye(k,l)*Bbar(i,j)) &
													   + (4/9.*I1*c1 + (16/9.*I1**2 - 8/3.*I1)*c2 + 4*I1*(I1-1)*(I1-3)*c3)*eye(i,j)*eye(k,l) &
													   + (2/3.*I1*c1 + 4/3.*I1*(I1-3)*c2 + (2*I1*(I1-3)**2)*c3)*(eye(i,k)*eye(j,l) + eye(i,l)*eye(j,k)) &
													   + (eye(i,j)*eye(k,l) - (eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)))*Ja*pressure
						end do
					end do
				end do
			end do	
		else
			mu1 = materialprops(2)
			mu2 = materialprops(3)
			K1 = materialprops(4)
			do i=1,ned
				do j=1,nsd
					do k=1,ned
						do l=1,nsd
							materialstiffness(i,j,k,l) = mu2*(2*Bbar(i,j)*Bbar(k,l)-(Bbar(i,k)*Bbar(j,l)+Bbar(i,l)*Bbar(j,k))) &
										 - 2/3.*(mu1+2*mu2*I1)*(Bbar(k,l)*eye(i,j)+Bbar(i,j)*eye(k,l)) &
										 + 4/3.*mu2*(BB(i,j)*eye(k,l)+BB(k,l)*eye(i,j)) &
										 + 2/9.*(mu1*I1+4*mu2*I2)*eye(i,j)*eye(k,l) &
										 + 1/3.*(mu1*I1+2*mu2*I2)*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)) &
										 !+ (Ja*(2*Ja-1)*eye(i,j)*eye(k,l) - Ja*(Ja-1)*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)))*K1
										 + (eye(i,j)*eye(k,l) - (eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)))*Ja*pressure
						end do
					end do
				end do
			end do
		end if
	
		
	end function materialstiffness
	
	function Kirchhoffstress(nsd,ned,B,Ja,pressure,materialprops)
		! Compute the Kirchhoff stress
		implicit none
		
		integer, intent(in) :: nsd, ned
		real(8), intent(in) :: Ja, pressure
		real(8), dimension(nsd,nsd), intent(in) :: B
		real(8), dimension(5), intent(in) :: materialprops
		real(8), dimension(ned,nsd) :: Kirchhoffstress
		real(8), dimension(nsd,nsd) :: Bbar, BB, eye
		integer :: i, j, k
		real(8) :: mu1, mu2, K1, I1, I2, c1, c2, c3
		
		! square matrix
		do i=1,nsd
			do j=1,nsd
				if(i==j) then
					eye(i,j) = 1.
				else
					eye(i,j) = 0.
				end if
			end do
		end do
		! Change B to Bbar and compute Bbar squared
		Bbar = B/Ja**(2/3.)
		BB = matmul(Bbar,Bbar)
		! I1 I2 denote I1bar I2bar actually
		I1 = 0.
		do i=1,nsd
			I1 = I1 + Bbar(i,i)
		end do
	    I2 = I1**2
		do i=1,nsd
			do j=1,nsd
				I2 = I2 - Bbar(i,j)**2
			end do 
		end do
		if (nsd == 2) then
			I1 = I1 + 1.
			I2 = I2 - 1.
		end if
		I2 = I2/2.
		Kirchhoffstress = 0.
		
		if (dAbs(materialprops(1)-4) < 1d-4) then
			c1 = materialprops(2)
			c2 = materialprops(3)
			c3 = materialprops(5)
			do i = 1, ned
				do j = 1, nsd
					Kirchhoffstress(i,j) = (2*c1 + 4*c2*(I1 - 3) + 6*c3*(I1-3)**2)*(Bbar(i,j) - 1/3.*I1*eye(i,j)) + Ja*pressure*eye(i,j)
				end do
			end do
		else
			mu1 = materialprops(2)
			mu2 = materialprops(3)
			K1 = materialprops(4)
			do i=1,ned
				do j=1,nsd
					Kirchhoffstress(i,j) = -1/3.*(mu1*I1+2*mu2*I2)*eye(i,j) - mu2*BB(i,j) + (mu1+mu2*I1)*Bbar(i,j) + Ja*pressure*eye(i,j);
				end do
			end do
		end if
	end function Kirchhoffstress	
	
end module mixed_material