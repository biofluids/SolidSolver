module tagentstiffness
	
	
	implicit none
	
contains
	function tagent_internal(dofs)
		use read_file, only: nsd, ned, nn, coords, nel, nen, connect, materialprops
		use shapefunction
		use integration
		use material
		implicit none
		real(8), dimension(nn*ned), intent(in) :: dofs
		real(8), dimension(nn*ned,nn*ned) :: tagent_internal
		real(8), dimension(nsd,nen) :: elecoord
		real(8), dimension(ned,nen) :: eledof
		real(8), dimension(nen,nsd) :: dNdx, dNdy
		real(8), dimension(ned,nsd) :: stress
		real(8), dimension(nen*ned) :: fele
		real(8), dimension(nsd) :: xi
		real(8), dimension(nen,nsd) :: dNdxi 
		real(8), dimension(nsd,nsd) :: dxdxi, dxidx, F, Finv, B, eye
		real(8), allocatable, dimension(:,:) :: xilist, xilist1
		real(8), allocatable, dimension(:) :: weights, weights1
		integer :: ele,a,i,npt,j,row,intpt,npt1
		real(8) :: det, Ja
		real(8), dimension(nsd) :: work ! for lapack inverse
		integer, dimension(nsd) :: ipiv ! for lapack inverse
		integer :: info, n1 ! for lapack inverse
		
		external DGETRF
		external DGETRI
		n1 = nsd
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