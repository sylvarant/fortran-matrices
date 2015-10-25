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
!         Name:  test_distance
!  Description:  test the minimal distance finding function - mostly used for
!                trouble shooting, test_steps is much better
! ====================================================================================
program test_distance
use zmatrix
use ieee_arithmetic
implicit none

integer :: i,j,bl
real(kind=dp),dimension(4,4) :: a,b,r,s,steps
real(kind=dp),dimension(:,:),allocatable :: c,d
type(z_matrix) :: za,zb,zs

a = 0   
a(1,2) = ieee_value(a(1,2),ieee_positive_inf)
a(1,3) = ieee_value(a(1,3),ieee_positive_inf)
a(1,4) = ieee_value(a(1,4),ieee_positive_inf)
a(2,1) = 1
a(2,3) = ieee_value(a(2,3),ieee_positive_inf)
a(2,4) = ieee_value(a(2,4),ieee_positive_inf)
a(3,1) = ieee_value(a(3,1),ieee_positive_inf)
a(3,2) = 1
a(3,4) = ieee_value(a(3,4),ieee_positive_inf)
a(4,1) = ieee_value(a(4,1),ieee_positive_inf)
a(4,2) = ieee_value(a(4,2),ieee_positive_inf)
a(4,3) = 1

r = 0
r(1,2) = ieee_value(r(1,2),ieee_positive_inf)
r(1,3) = ieee_value(r(1,3),ieee_positive_inf)
r(1,4) = ieee_value(r(1,4),ieee_positive_inf)
r(2,1) = 1
r(2,3) = ieee_value(r(2,3),ieee_positive_inf)
r(2,4) = ieee_value(r(2,4),ieee_positive_inf)
r(3,1) = 2
r(3,2) = 1
r(4,1) = 3
r(4,2) = 2
r(4,3) = 1
r(3,4) = ieee_value(r(3,4),ieee_positive_inf)

steps = 0
steps(2,1) = 1
steps(3,1) = 2
steps(3,2) = 1
steps(4,1) = 3
steps(4,2) = 2
steps(4,3) = 1

bl = 1

do while (bl .le. 4)
    b = 0
    s = 0

    call matrix_to_zmatrix(a,za,bl)
    call matrix_to_zmatrix(b,zb,bl)
    call matrix_to_zmatrix(s,zs,bl)

    call find_distance(za,zb,steps=zs)

    call zmatrix_to_matrix(zb,c)
    call zmatrix_to_matrix(zs,d)

    do i = 1,4
        do j = 1,4
            if((abs(r(i,j)-c(i,j))) .gt. 0.000001_dp )then
                    print *,"error !! for result :: ",i,j,c(i,j),r(i,j)
            endif
            if (steps(i,j) .ne. d(i,j)) then
                print *,"error !! for steps :: ",d(i,j),steps(i,j)
            endif
        enddo
    enddo

    bl = bl * 2
    deallocate(c,d,za%matrix,zb%matrix,zs%matrix)
enddo

end program test_distance
