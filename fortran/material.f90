module material
    implicit none
contains
    subroutine getStress(theta, F, stress)
        use read_file, only: nsd, materialprops, delta
        real(8), intent(in) :: theta
        real(8), dimension(nsd, nsd), intent(in) :: F
        real(8), dimension(nsd, nsd), intent(inout) :: stress
        real(8), dimension(nsd, nsd) :: Fe
        real(8) :: mu, lambda, Je

        mu = materialprops(3)
        lambda = materialprops(4)
        Fe = F/theta
        if (nsd == 2) then
            Je = (Fe(1,1)*Fe(2,2) - Fe(1,2)*Fe(2,1))
        else if (nsd == 3) then
            Je = Fe(1,1)*Fe(2,2)*Fe(3,3) - Fe(1,1)*Fe(3,2)*Fe(2,3) &
                    - Fe(1,2)*Fe(2,1)*Fe(3,3) + Fe(1,2)*Fe(2,3)*Fe(3,1) &
                    + Fe(1,3)*Fe(2,1)*Fe(3,2) - Fe(1,3)*Fe(2,2)*Fe(3,1)
        end if
        stress = (lambda*log(Je) - mu)*delta + mu*matmul(Fe, transpose(Fe))
    end subroutine getStress

    subroutine update(F, theta, Fe, Je, Be, Ce, trace_Ce, trace_Me)
        use read_file, only: nsd, delta, materialprops
        real(8), intent(in) :: theta
        real(8), intent(inout) :: Je, trace_Ce, trace_Me
        real(8), dimension(nsd, nsd), intent(in) :: F
        real(8), dimension(nsd, nsd), intent(inout) :: Fe, Be, Ce
        real(8), dimension(nsd, nsd) :: Me
        real(8) :: lambda, mu
        integer :: i, j
        Fe = F/theta
        if (nsd == 2) then
            Je = (Fe(1,1)*Fe(2,2) - Fe(1,2)*Fe(2,1))
        else if (nsd == 3) then
            Je = Fe(1,1)*Fe(2,2)*Fe(3,3) - Fe(1,1)*Fe(3,2)*Fe(2,3) &
                    - Fe(1,2)*Fe(2,1)*Fe(3,3) + Fe(1,2)*Fe(2,3)*Fe(3,1) &
                    + Fe(1,3)*Fe(2,1)*Fe(3,2) - Fe(1,3)*Fe(2,2)*Fe(3,1)
        end if
        mu = materialprops(3)
        lambda = materialprops(4)
        Ce = matmul(transpose(Fe),Fe)
        Be = matmul(Fe, transpose(Fe))
        Me = (lambda*log(Je) - mu)*delta + mu*Ce
        trace_Me = 0.0
        trace_Ce = 0.0
        do i = 1, nsd
            trace_Me = trace_Me + Me(i, i)
            trace_Ce = trace_Ce + Ce(i,i)
        end do
    end subroutine update

    subroutine growth_increment(theta_pre, theta, F, stress, mstiff)
        use read_file, only: nsd, dt, tol, materialprops, delta, step
        real(8), intent(inout) :: theta, theta_pre
        real(8), dimension(nsd, nsd), intent(in) :: F
        real(8), dimension(nsd, nsd), intent(inout) :: stress
        real(8), dimension(nsd, nsd, nsd, nsd), intent(inout) :: mstiff
        real(8), dimension(nsd, nsd) :: Fe, Ce, Be
        real(8) :: mu, lambda, Je, theta_max, tau_g, gamma, kg, trace_Me, Me_cr, trace_Ce
        real(8) :: dr1, dr2, residual, local_tangent
        integer :: i, j, k, l, cnt
        logical :: flag
        theta_max = 1.5
        tau_g = 1.0
        gamma = 2.0
        Me_cr = 1.2d4
        mu = materialprops(3)
        lambda = materialprops(4)
        residual = 1.0
        theta = theta_pre
        call update(F, theta, Fe, Je, Be, Ce, trace_Ce, trace_Me)
        cnt = 0
        flag = .false.
        ! local Newton iteration
        do while (residual > tol .and. dabs(theta_max - theta) > tol .and. trace_Me - Me_cr > tol)
            kg = 1.0/tau_g*((theta_max - theta_pre)/(theta_max - 1.0))**gamma
            residual = theta - theta_pre - kg*(trace_Me-Me_cr)*dt
            dr1 = -(9*lambda + 2*mu*trace_Ce)/theta ! partial derivative of phi w.r.t. theta
            dr2 = -gamma*kg/(theta_max - theta) ! partial derivative of k w.r.t. theta
            local_tangent = 1 - (kg*dr1 + (trace_Me-Me_cr)*dr2)*dt
            theta = theta - residual/local_tangent
            flag = .true.
            call update(F, theta, Fe, Je, Be, Ce, trace_Ce, trace_Me)
            cnt = cnt + 1
            if (theta - theta_max > tol) then
                write(*,*) 'maximum value exceeded, theta = ', theta
                stop
            else if (cnt == 500) then
                write(*,*) "max number of theta iteration reached!"
                stop
            end if
        end do
        mstiff = 0.0
        do i = 1, nsd
            do j = 1, nsd
                do k = 1, nsd
                    do l = 1, nsd
                        mstiff(i, j, k, l) = (mu - lambda*log(Je))*(delta(i, k)*delta(j, l) + delta(i, l)*delta(j, k)) &
                            + lambda*(delta(i, j)*delta(k, l))
                        if (flag) then
                            mstiff(i, j, k, l) = mstiff(i, j, k, l) - kg*dt/(local_tangent*theta)&
                                *(3*lambda*delta(i,j)+2*mu*Be(i,j))*(3*lambda*delta(k,l)+2*mu*Be(k,l))
                        end if
                    end do
                end do
            end do
        end do
        stress = (lambda*log(Je) - mu)*delta + mu*matmul(Fe, transpose(Fe))
   end subroutine

    subroutine Kirchhoffstress(nsd, intcoord, F, materialtype, materialprops, stress)
        ! Compute the Kirchhoff stress
        integer, intent(in) :: nsd, materialtype
        real(8), dimension(nsd,nsd), intent(inout) :: F
        real(8), dimension(nsd), intent(in) :: intcoord
        real(8), dimension(5), intent(in) :: materialprops
        real(8), dimension(nsd,nsd), intent(inout) :: stress
        real(8), dimension(nsd,nsd) :: B, Bbar, BB, eye, C
        integer :: i, j, k
        real(8) :: kappa, mu1, mu2, I1, I2, c1, c2, c3, Ja, I4, I6, lambda4, lambda6, kk1, kk2, R, beta
        real(8), dimension(nsd) :: a0, g0, a, g

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
        if (nsd == 2) then
            Ja = (F(1,1)*F(2,2) - F(1,2)*F(2,1))
        else if (nsd == 3) then
            Ja = F(1,1)*F(2,2)*F(3,3) - F(1,1)*F(3,2)*F(2,3) &
                    - F(1,2)*F(2,1)*F(3,3) + F(1,2)*F(2,3)*F(3,1) &
                    + F(1,3)*F(2,1)*F(3,2) - F(1,3)*F(2,2)*F(3,1)
        end if
        B = matmul(F,transpose(F))
        Bbar = B/Ja**(2/3.)
        BB = matmul(Bbar,Bbar)

        ! I1 I2 denote I1bar I2bar actually
        if (nsd == 2) then
            I1 = 1.0/Ja**(2/3.)
        else if (nsd == 3) then
            I1 = 0.0
        end if
        do i = 1, nsd
            I1 = I1 + Bbar(i,i)
        end do
        I2 = I1**2
        do i = 1, nsd
            do j = 1, nsd
                I2 = I2 - Bbar(i,j)**2
            end do
        end do
        if (nsd == 2) then
            I2 = I2 - 1.0/Ja**(4/3.)
        end if
        I2 = I2/2.

        stress = 0.

        if (materialtype == 1) then ! Mooney-Rivlin
            kappa = materialprops(2)
            mu1 = materialprops(3)
            mu2 = materialprops(4)
            do i=1,nsd
                do j=1,nsd
                    stress(i,j) = -(mu1*I1+2*mu2*I2)*eye(i,j)/3. - mu2*BB(i,j) + (mu1+mu2*I1)*Bbar(i,j) + Ja*(Ja-1)*kappa*eye(i,j)
                end do
            end do
        else if (materialtype == 2) then ! Yeoh
            kappa = materialprops(2)
            c1 = materialprops(3)
            c2 = materialprops(4)
            c3 = materialprops(5)
            do i = 1, nsd
                do j = 1, nsd
                    stress(i,j) = (2*c1 + 4*c2*(I1 - 3) + 6*c3*(I1-3)**2)*(Bbar(i,j) - 1/3.*I1*eye(i,j)) + Ja*(Ja-1)*kappa*eye(i,j)
                end do
            end do
        else if (materialtype == 3) then ! HGO
            kappa = materialprops(2)
            mu1 = materialprops(3)
            kk1 = materialprops(4)
            kk2 = materialprops(5)
            ! Material parameters are hard-coded
            if (intcoord(1)**2 + intcoord(2)**2 <= (0.97d-3)**2) then
                beta = 29*3.14159/180
                kk1 = 2.3632d3
                kk2 = 0.8393
                mu1 = 3000
            else
                beta = 62*3.14159/180
                kk1 = 0.5620d3
                kk2 = 0.7112
                mu1 = 300
            end if

            R = sqrt((intcoord(1))**2 + (intcoord(2))**2)
            if (nsd == 3) then
                a0 = [cos(beta)*intcoord(2)/R, -cos(beta)*intcoord(1)/R,  sin(beta)]
                g0 = [cos(beta)*intcoord(2)/R, -cos(beta)*intcoord(1)/R, -sin(beta)]
            else if (nsd == 2) then
                a0 = [cos(beta)*intcoord(2)/R, -cos(beta)*intcoord(1)/R]
                g0 = [cos(beta)*intcoord(2)/R, -cos(beta)*intcoord(1)/R]
            end if
            C = matmul(transpose(F),F)
            I4 = dot_product(a0,matmul(C,a0))
            if (nsd == 2) then
                I4 = I4 + (sin(beta))**2
            end if
            lambda4 = sqrt(I4)
            a = matmul(F,a0)/lambda4
            I4 = I4/Ja**(2/3.)
            I6 = dot_product(g0,matmul(C,g0))
            if (nsd == 2) then
                I6 = I6 + (sin(beta))**2
            end if
            lambda6 = sqrt(I6)
            g = matmul(F,g0)/lambda6
            I6 = I6/Ja**(2/3.)

            do i = 1, nsd
                do j = 1, nsd
                    stress(i,j) = -1/3.*mu1*I1*eye(i,j) + mu1*Bbar(i,j) + Ja*(Ja-1)*kappa*eye(i,j)
                    ! anisotropic part
                    if (lambda4 > 1.0) then
                        stress(i,j) = stress(i,j) + 2*kk1*( (I4-1)*exp(kk2*(I4-1)**2)*(a(i)*a(j) - 1/3.*eye(i,j))*I4 )
                    end if
                    if (lambda6 > 1.0) then
                        stress(i, j) = stress(i, j) + 2*kk1*( (I6-1)*exp(kk2*(I6-1)**2)*(g(i)*g(j) - 1/3.*eye(i,j))*I6 )
                    end if
                end do
            end do
        end if
    end subroutine Kirchhoffstress

    subroutine materialstiffness(nsd, intcoord, F, materialtype, materialprops, mstiff)
        ! Compute the material stiffness tensor
        integer, intent(in) :: nsd, materialtype
        real(8), dimension(nsd,nsd), intent(in) :: F
        real(8), dimension(nsd), intent(in) :: intcoord
        real(8), dimension(5), intent(in) :: materialprops
        real(8), dimension(nsd, nsd, nsd, nsd), intent(out):: mstiff

        integer :: i, j, k, l, p
        real(8) :: mu1, mu2, kappa, I1, I2, c1, c2, c3, Ja, I4, I6, lambda4, lambda6, kk1, kk2, der14, der24, der16, der26, R, beta
        real(8), dimension(nsd,nsd) :: B, Bbar, BB, eye, C
        real(8), dimension(nsd) :: a0, g0, a, g

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

        mstiff = 0.

        ! Change B to Bbar and compute Bbar squared
        if (nsd == 2) then
            Ja = (F(1,1)*F(2,2) - F(1,2)*F(2,1))
        else if (nsd == 3) then
            Ja = F(1,1)*F(2,2)*F(3,3) - F(1,1)*F(3,2)*F(2,3) &
                    - F(1,2)*F(2,1)*F(3,3) + F(1,2)*F(2,3)*F(3,1) &
                    + F(1,3)*F(2,1)*F(3,2) - F(1,3)*F(2,2)*F(3,1)
        end if
        B = matmul(F,transpose(F))
        Bbar = B/Ja**(2/3.)
        BB = matmul(Bbar,Bbar)

        ! I1 I2 denote I1bar I2bar actually
        if (nsd == 2) then
            I1 = 1.0/Ja**(2/3.)
        else if (nsd == 3) then
            I1 = 0.0
        end if
        do i = 1, nsd
            I1 = I1 + Bbar(i,i)
        end do
        I2 = I1**2
        do i = 1, nsd
            do j = 1, nsd
                I2 = I2 - Bbar(i,j)**2
            end do
        end do
        if (nsd == 2) then
            I2 = I2 - 1.0/Ja**(4/3.)
        end if
        I2 = I2/2.

        if (materialtype == 1) then ! Mooney-Rivlin
            kappa = materialprops(2)
            mu1 = materialprops(3)
            mu2 = materialprops(4)
            do i=1,nsd
                do j=1,nsd
                    do k=1,nsd
                        do l=1,nsd
                            mstiff(i,j,k,l) = mu2*(2*Bbar(i,j)*Bbar(k,l)-(Bbar(i,k)*Bbar(j,l)+Bbar(i,l)*Bbar(j,k))) &
                                            - 2/3.*(mu1+2*mu2*I1)*(Bbar(k,l)*eye(i,j)+Bbar(i,j)*eye(k,l)) &
                                            + 4/3.*mu2*(BB(i,j)*eye(k,l)+BB(k,l)*eye(i,j)) &
                                            + 2/9.*(mu1*I1+4*mu2*I2)*eye(i,j)*eye(k,l) &
                                            + 1/3.*(mu1*I1+2*mu2*I2)*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)) &
                                            + (Ja*(2*Ja-1)*eye(i,j)*eye(k,l) - Ja*(Ja-1)*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k))) &
                                            *kappa
                                            !+ (eye(i,j)*eye(k,l) - (eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)))*Ja*pressure
                        end do
                    end do
                end do
            end do
        else if (materialtype == 2) then ! Yeoh
            kappa = materialprops(2)
            c1 = materialprops(3)
            c2 = materialprops(4)
            c3 = materialprops(5)
            do i = 1, nsd
                do j = 1, nsd
                    do k = 1, nsd
                        do l = 1, nsd
                            mstiff(i,j,k,l) = (8*c2 + 24*c3*(I1 - 3))*Bbar(i,j)*Bbar(k,l) &
                                            - (4/3.*c1 + (16/3.*I1 - 8)*c2 + 12*(I1-1)*(I1-3)*c3) &
                                            *(eye(i,j)*Bbar(k,l) + eye(k,l)*Bbar(i,j)) &
                                            + (4/9.*I1*c1 + (16/9.*I1**2 - 8/3.*I1)*c2 + 4*I1*(I1-1)*(I1-3)*c3)*eye(i,j)*eye(k,l) &
                                            + (2/3.*I1*c1 + 4/3.*I1*(I1-3)*c2 + (2*I1*(I1-3)**2)*c3) &
                                            *(eye(i,k)*eye(j,l) + eye(i,l)*eye(j,k)) &
                                            + (Ja*(2*Ja-1)*eye(i,j)*eye(k,l) - Ja*(Ja-1)*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k))) &
                                            *kappa
                                            !+ (eye(i,j)*eye(k,l) - (eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)))*Ja*pressure
                        end do
                    end do
                end do
            end do
        else if (materialtype == 3) then ! HGO
            kappa = materialprops(2)
            mu1 = materialprops(3)
            kk1 = materialprops(4)
            kk2 = materialprops(5)
            ! Material parameters are hard-coded
            if (intcoord(1)**2 + intcoord(2)**2 <= (0.97d-3)**2) then
                beta = 29*3.14159/180
                kk1 = 2.3632d3
                kk2 = 0.8393
                mu1 = 3000
            else
                beta = 62*3.14159/180
                kk1 = 0.5620d3
                kk2 = 0.7112
                mu1 = 300
            end if

            R = sqrt((intcoord(1))**2 + (intcoord(2))**2)
            if (nsd == 3) then
                a0 = [cos(beta)*intcoord(2)/R, -cos(beta)*intcoord(1)/R,  sin(beta)]
                g0 = [cos(beta)*intcoord(2)/R, -cos(beta)*intcoord(1)/R, -sin(beta)]
            else if (nsd == 2) then
                a0 = [cos(beta)*intcoord(2)/R, -cos(beta)*intcoord(1)/R]
                g0 = [cos(beta)*intcoord(2)/R, -cos(beta)*intcoord(1)/R]
            end if
            C = matmul(transpose(F),F)
            I4 = dot_product(a0,matmul(C,a0))
            if (nsd == 2) then
                I4 = I4 + (sin(beta))**2
            end if
            lambda4 = sqrt(I4)
            a = matmul(F,a0)/lambda4
            I4 = I4/Ja**(2/3.)
            I6 = dot_product(g0,matmul(C,g0))
            if (nsd == 2) then
                I6 = I6 + (sin(beta))**2
            end if
            lambda6 = sqrt(I6)
            g = matmul(F,g0)/lambda6
            I6 = I6/Ja**(2/3.)

            ! 1st order derivative of psi w.r.t. I4/I6
            der14 = kk1*(I4-1.)*exp(kk2*(I4-1.)**2)
            der16 = kk1*(I6-1.)*exp(kk2*(I6-1.)**2)
            ! 2nd order derivative of psi w.r.t. I4/I6
            der24 = kk1*(1+2*kk2*(I4-1.)**2)*exp(kk2*(I4-1.)**2)
            der26 = kk1*(1+2*kk2*(I6-1.)**2)*exp(kk2*(I6-1.)**2)

            do i=1,nsd
                do j=1,nsd
                    do k=1,nsd
                        do l=1,nsd
                            mstiff(i,j,k,l) = - 2/3.*(mu1)*(Bbar(k,l)*eye(i,j)+Bbar(i,j)*eye(k,l)) &
                                            + 2/9.*(mu1*I1)*eye(i,j)*eye(k,l) &
                                            + 1/3.*(mu1*I1)*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)) &
                                            + (Ja*(2*Ja-1)*eye(i,j)*eye(k,l) - Ja*(Ja-1)*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k))) &
                                            *kappa
                                            !+ (eye(i,j)*eye(k,l) - (eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k)))*Ja*pressure
                            if (lambda4 > 1.0) then
                                ! anisotropic part 1
                                mstiff(i,j,k,l) = mstiff(i,j,k,l) + 4*I4**2*der24*a(i)*a(j)*a(k)*a(l) &
                                                - 4/3.*(I4*der24+der14)*I4*(a(i)*a(j)*eye(k,l) + eye(i,j)*a(k)*a(l)) &
                                                + (4/9.*I4**2*der24 + 4/9.*I4*der14)*eye(i,j)*eye(k,l) &
                                                + 2/3.*I4*der14*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k))
                            end if
                            if (lambda6 > 1.0) then
                                ! anisotropic part 2
                                mstiff(i,j,k,l) = mstiff(i,j,k,l) + 4*I6**2*der26*g(i)*g(j)*g(k)*g(l) &
                                                - 4/3.*(I6*der26+der16)*I6*(g(i)*g(j)*eye(k,l) + eye(i,j)*g(k)*g(l)) &
                                                + (4/9.*I6**2*der26 + 4/9.*I6*der16)*eye(i,j)*eye(k,l) &
                                                + 2/3.*I6*der16*(eye(i,k)*eye(j,l)+eye(i,l)*eye(j,k))
                            end if
                        end do
                    end do
                end do
            end do
        end if
    end subroutine materialstiffness

end module material
