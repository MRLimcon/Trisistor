module utils

  contains

    function generate_doping(m, is_activated) result(polos)
        use OMP_LIB
        implicit none

        integer, intent(in) :: m, is_activated
        integer, parameter :: n = 100
        integer :: i,j, polos(m,m)
        integer :: posicoes(m,m)
        real :: r(m,m), r2(m,m)

        call RANDOM_NUMBER(r)
        call RANDOM_NUMBER(r2)

        polos = 0

        !$OMP PARALLEL DO SIMD
        do i = 1, m
            do j = 1, m
                if ( r(i,j) < 0.06 ) then
                    if ( i < (m*0.25) .or. (i > m*0.45 .and. (j > m*0.8 .or. j < m*0.2))) then
                        polos(i,j) = -1
                    else if ( i >= (m*0.25) .and. i <= (m*0.45) ) then
                        polos(i,j) = 1
                    end if
                end if
            end do
        end do
        !$OMP END PARALLEL DO SIMD

        do i = 1, m, 2
            polos(1,i) = 4
        end do

        if (is_activated == 1) then
          do i = 1, 45, 2
              polos(i,m) = 5
          end do
        end if

        polos(90:m,m) = -3
        polos(1,0:m) = 1
        polos(m,80:m) = -4

    end function generate_doping

    function generate_field(polos,m) result(campo)
        use OMP_LIB
        implicit none

        INteger, INTENT(IN) :: m, polos(m,m)
        real :: c1 = 1!, quad!.4418*(10**-9), c2 = 111*(10**-12)
        real :: d, comp_angular1, comp_angular2, campo(m,m,2)
        integer :: i,j,k,l

        campo = 0

        !$OMP PARALLEL DO SIMD
        do i = 1,m
            do j = 1,m
                do k = 1,m
                    do l = 1,m
                        if ( k /= i .and. l /= j .and. polos(l,k) /= 0) then
                            d = (k-i)**2 + (l-j)**2
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

    end function generate_field

    function generate_electrons(m) result(elecbur)
        use OMP_LIB
        implicit none

        INteger, INTENT(IN) :: m
        integer, parameter :: n = 1000
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

    end function generate_electrons

    function calculate_density(elec_bur,polos,m) result(densidade)
        use OMP_LIB
        implicit none

        INteger, INTENT(IN) :: m, elec_bur(m,m), polos(m,m)
        integer :: i,j,k,l, densidade(m,m), calculo(m+8,m+8), polos2(m,m)

        polos2 = polos
        polos2(1:m,1) = 0
        polos2(1:m,m) = 0
        polos2(1,1:m) = 0
        polos2(m,1:m) = 0
        polos2(35:m,m) = 0

        calculo = 0
        calculo(5:m+4,5:m+4) = elec_bur + polos2!

        !$OMP PARALLEL DO SIMD
        do i = 5,m+4
            do j = 5,m+4
                densidade(j-4,i-4) = sum(calculo(j-4:j+4,i-4:i+4))
            end do
        end do
        !$OMP END PARALLEL DO SIMD

    end function calculate_density

end module utils
