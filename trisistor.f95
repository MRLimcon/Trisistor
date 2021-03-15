module trisis
    
    contains

        function geracao_polos(m) result(polos)
            use OMP_LIB
            implicit none
            
            INteger, INTENT(IN) :: m
            integer, parameter :: n = 500
            integer :: i, polos(m,m)
            integer, allocatable :: posicoes(:,:)
            real, allocatable :: r(:,:), r2(:)

            allocate(posicoes(2,n))
            allocate(r(2,n))
            allocate(r2(n))

            call RANDOM_NUMBER(r)
            call RANDOM_NUMBER(r2)
            
            !$OMP PARALLEL DO SIMD
            do i = 1, n
                if ( r2(i) <=0.531914894 ) then
                    posicoes(1,i) = (r(1,i) * (m*0.24)) + 1
                    posicoes(2,i) = (r(2,i) * m-1) + 1
                elseif ( r2(i) > 0.531914894 .and. r2(i) <= 0.765957447 ) then
                    posicoes(1,i) = (r(1,i) * (m*0.55)) + m*0.45
                    posicoes(2,i) = (r(2,i) * m*0.19) + 1
                elseif ( r2(i) > 0.765957447 ) then
                    posicoes(1,i) = (r(1,i) * (m*0.55)) + m*0.45
                    posicoes(2,i) = (r(2,i) * m*0.2) + m*0.8
                end if
            end do
            !$OMP END PARALLEL DO SIMD
            
            !$OMP PARALLEL DO SIMD
            do i = 1,n
                polos(posicoes(1,i),posicoes(2,i)) = -1
            end do
            !$OMP END PARALLEL DO SIMD
            
            deallocate(posicoes)
            allocate(posicoes(2,n))

            !$OMP PARALLEL DO SIMD
            do i = 1, n
                posicoes(1,i) = (r(1,i) * (m*0.2)) + m*0.25
                posicoes(2,i) = r(2,i) * m
            end do
            !$OMP END PARALLEL DO SIMD
            
            !$OMP PARALLEL DO SIMD
            do i = 1,n
                polos(posicoes(1,i),posicoes(2,i)) = 1
            end do
            !$OMP END PARALLEL DO SIMD

            deallocate(posicoes)
            deallocate(r)
            deallocate(r2)
            
            
            !polos(45,45:55) = -7
            polos(55:m,1) = 7
            polos(1,45:55) = -7
            polos(55:m,m) = 7
            polos(25:60,m) = 7

        end function geracao_polos

        function geracao_campo(polos,m) result(campo)
            use OMP_LIB
            implicit none
            
            INteger, INTENT(IN) :: m, polos(m,m)
            real :: c1 = 1, quad!.4418*(10**-9), c2 = 111*(10**-12)
            real :: d, comp_angular1, comp_angular2, campo(m,m,2)
            integer :: i,j,k,l

            campo = 0

            !$OMP PARALLEL DO SIMD 
            do i = 1,m
                do j = 1,m
                    do k = 1,m
                        do l = 1,m
                            quad = (k-i)**2 + (l-j)**2
                            if ( k /= i .and. l /= j .and. polos(l,k) /= 0 .and. &
                            sqrt(quad) <= 100 ) then
                                d = ((l-j)**2)+((k-i)**2)!*c2**2
                                comp_angular1 = (j-l)/sqrt(d)
                                comp_angular2 = (i-k)/sqrt(d)
                                campo(j,i,1) = campo(j,i,1) + comp_angular2*(c1*polos(l,k)/d)
                                campo(j,i,2) = campo(j,i,2) + comp_angular1*(c1*polos(l,k)/d)
                            end if
                        end do
                    end do
                end do
            end do
            !$OMP END PARALLEL DO SIMD

        end function geracao_campo
        
        function geracao_bur_elec(m) result(elecbur)
            use OMP_LIB
            implicit none
            
            INteger, INTENT(IN) :: m
            integer, parameter :: n = 1500
            integer :: i, elecbur(m,m), posicoes(2,n)
            real :: r(2,n), r2(n)

            call RANDOM_NUMBER(r)
            call RANDOM_NUMBER(r2)

            elecbur = 0

            !$OMP PARALLEL DO SIMD
            do i = 1, n
                posicoes(1:2,i) = (r(1:2,i)*m-1)+1
            end do
            !$OMP END PARALLEL DO SIMD

            !$OMP PARALLEL DO SIMD
            do i = 1, n
                if ( r2(i) >= 0.5 ) then
                    elecbur(posicoes(1,i),posicoes(2,i)) = 1
                else
                    elecbur(posicoes(1,i),posicoes(2,i)) = -1
                end if
            end do
            !$OMP END PARALLEL DO SIMD

            elecbur(46:m,21:79) = 0

        end function geracao_bur_elec

        function atualizacao_campo(elec_bur,campoinit,m) result(campo)
            use OMP_LIB
            implicit none
            
            INteger, INTENT(IN) :: m, elec_bur(m,m)
            real, intent(in) :: campoinit(m,m,2)
            real :: c1 = 1, quad!.4418*(10**-9), c2 = 111*(10**-12)
            real :: d, comp_angular1, comp_angular2, campo(m,m,2), polos(m,m)
            integer :: i,j,k,l

            polos  = elec_bur
            campo = campoinit

            !$OMP PARALLEL DO SIMD 
            do i = 1,m
                do j = 1,m
                    do k = 1,m
                        do l = 1,m
                            quad = (k-i)**2 + (l-j)**2
                            if ( k /= i .and. l /= j .and. polos(l,k) /= 0 .and. &
                            sqrt(quad) <= 100 ) then
                                d = ((l-j)**2)+((k-i)**2)!*c2**2
                                comp_angular1 = (j-l)/sqrt(d)
                                comp_angular2 = (i-k)/sqrt(d)
                                campo(j,i,1) = campo(j,i,1) + comp_angular2*(c1*polos(l,k)/d)
                                campo(j,i,2) = campo(j,i,2) + comp_angular1*(c1*polos(l,k)/d)
                            end if
                        end do
                    end do
                end do
            end do
            !$OMP END PARALLEL DO SIMD

        end function atualizacao_campo

        function atualizacao_elecbur(elec_bur,campoinit,t,m) result(elecbur_fim)
            use OMP_LIB
            implicit none
            
            INteger, INTENT(IN) :: m, elec_bur(m,m), t
            real, intent(in) :: campoinit(m,m,2)
            real :: r(m,m,2)
            integer :: i,j,k(m,m,2), elecbur_fim(m,m), c(2),l,n

            elecbur_fim = elec_bur
            k = 0
            call RANDOM_NUMBER(r)

            r = (r-0.03)
            

            !$OMP PARALLEL DO SIMD
            do i = 1, m
                do j = 1, m
                    if ( elec_bur(j,i) /= 0 .and. elecbur_fim(j,i) /= 0) then
                        k(j,i,1) = j
                        k(j,i,2) = i
                        if ((campoinit(j,i,2)*r(j,i,1))*elec_bur(j,i) > 0 .and. &
                        elecbur_fim(j+1,i) /= elec_bur(j,i) ) then
                            k(j,i,1) = j + 1
                        elseif ((campoinit(j,i,2)*r(j,i,1))*elec_bur(j,i) < 0 .and. &
                        elecbur_fim(j-1,i) /= elec_bur(j,i)) then
                            k(j,i,1) = j - 1
                        end if
                        if ((campoinit(j,i,1)*r(j,i,2))*elec_bur(j,i) > 0 .and. &
                        elecbur_fim(j,i+1) /= elec_bur(j,i)) then
                            k(j,i,2) = i + 1
                        elseif ((campoinit(j,i,1)*r(j,i,2))*elec_bur(j,i) < 0 .and. &
                        elecbur_fim(j,i-1) /= elec_bur(j,i) ) then
                            k(j,i,2) = i - 1
                        end if
                    end if
                end do
            end do
            !$OMP END PARALLEL DO SIMD

            !$OMP PARALLEL DO SIMD
            do i = 1,m
                do j = 1, m
                    if ( ((k(j,i,1) >= 45 .and. k(j,i,2) >= 20 .and. k(j,i,2) <= 80)) &
                    .or. (k(j,i,1) == j .and. k(j,i,2) == i) .or. elecbur_fim(k(j,i,1),k(j,i,2)) == elec_bur(j,i) ) then
                    
                    elseif ( k(j,i,1) >= 1 .and. k(j,i,1) <= m .and. k(j,i,2) >= 1 .and. k(j,i,2) <= m ) then

                        if ( elec_bur(j,i) == 1 ) then
                            elecbur_fim(k(j,i,1),k(j,i,2)) = elecbur_fim(k(j,i,1),k(j,i,2)) +1
                            elecbur_fim(j,i) = 0
                        elseif ( elec_bur(j,i) == -1 ) then
                            elecbur_fim(k(j,i,1),k(j,i,2)) = elecbur_fim(k(j,i,1),k(j,i,2)) -1
                            elecbur_fim(j,i) = 0
                        end if

                    end if
                end do
            end do
            !$OMP END PARALLEL DO SIMD

            !elecbur_fim(45,50) = -1
            elecbur_fim(2,49:51) = -1 
            !elecbur_fim(55:m,1) = 0
            !elecbur_fim(55:m,m) = 0

            if ( t == 1 ) then
                elecbur_fim(55:m,1:2) = 0
                elecbur_fim(55:m,m-1:m) = 0
            end if

        end function atualizacao_elecbur

        function c_densidade(elec_bur,polos,m) result(densidade)
            use OMP_LIB
            implicit none
            
            INteger, INTENT(IN) :: m, elec_bur(m,m), polos(m,m)
            integer :: i,j,k,l, densidade(m,m), calculo(m+8,m+8), polos2(m,m)

            polos2 = polos
            polos2(15:60,1) = 0
            polos2(15:60,m) = 0
            !polos2(45,45:55) = 0
            polos2(35:m,1) = 0
            polos2(1,30:(m-30)) = 0
            polos2(35:m,m) = 0

            calculo = 0
            calculo(5:m+4,5:m+4) = elec_bur+polos2!

            !$OMP PARALLEL DO SIMD
            do i = 5,m+4
                do j = 5,m+4
                    densidade(j-4,i-4) = sum(calculo(j-4:j+4,i-4:i+4))
                end do
            end do
            !$OMP END PARALLEL DO SIMD

        end function c_densidade

end module trisis