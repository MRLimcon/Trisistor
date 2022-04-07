module trisistor

    contains

        function update_field(elec_bur, campoinit, m) result(campo)
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
                            if ( k /= i .and. l /= j .and. polos(l,k) /= 0) then
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

        end function update_field

        function update_positions(elec_bur, campoinit, t, m) result(elecbur_fim)
            use OMP_LIB
            implicit none

            INteger, INTENT(IN) :: m, elec_bur(m,m), t
            real, intent(in) :: campoinit(m,m,2)
            real :: r(m,m,2)
            integer :: i,j,k(m,m,2), elecbur_fim(m,m), c(2),l,n

            elecbur_fim = elec_bur
            k = 0
            call RANDOM_NUMBER(r)

            r = (r-0.33)


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
            !elecbur_fim(55:m,1) = 0
            !elecbur_fim(55:m,m) = 0

            if ( t == 1 ) then
                !elecbur_fim(2,49:51) = -1
                elecbur_fim(m-1,83:87) = -1
            end if

        end function update_positions

end module trisistor
