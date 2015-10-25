! Adriaan Larmuseau - s0190593 - 1ste Master CW
! ## Architecture :: 
! (uname -r) Linux VPCEA2S1E 2.6.35-25-generic #44-Ubuntu SMP Fri Jan 21 17:40:44 UTC 2011 x646_64 GNU/Linux
! ## Compilers :: 
! G95 (GCC 4.1.2 (g95 0.93!) Jun 16 2010)
! ## Lib ::
! blas -- version 1.2-7
! ## Compiler Commando ::
! zie makefile
! ## Bestede tijd ::
! 1 uur

! ===  Program  ======================================================================
!         Name:  test_steps
!  Description:  test the correctness of the steps functionality
! ====================================================================================
program test_steps
use zmatrix
use ieee_arithmetic
implicit none

integer :: it,h,b,i,j,n
type(z_matrix) :: zg,zr,zs
real(kind=dp),dimension(:,:),allocatable :: g,r,s,rs,e

do  i = 3,8
    n = 2 ** i 
    b = 1
    do while (b .le. n)

        allocate(g(n,n),r(n,n),s(n,n),e(n,n)) 

        g = ieee_value(g,ieee_positive_inf)
        do j = 1,n
            g(j,j) = 0
            if(j+1 .le. n) then
                g(j+1,j) = 1
            endif
        enddo

        r = 0
        s = 0

        e = 0
        do j = 1,n 
            it = 1
            do h = j+1,n
               e(h,j) = it
               it = it +1
            enddo
        enddo

        call matrix_to_zmatrix(g,zg,b)
        call matrix_to_zmatrix(r,zr,b)
        call matrix_to_zmatrix(s,zs,b)

        call find_distance(zg,zr,steps=zs)

        call zmatrix_to_matrix(zs,rs)

        do j = 1,n 
            do h = 1,n
                if((abs(rs(h,j)-e(h,j))) .gt. 0.000001_dp )then
                    print *,"error !! for steps :: ",h,j,e(h,j),rs(h,j)
                endif            
            enddo
        enddo

        deallocate(e,g,r,s,rs,zg%matrix,zr%matrix,zs%matrix)

        b = b *2

    enddo
enddo


end program test_steps
