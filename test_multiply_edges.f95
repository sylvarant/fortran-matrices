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
! 20 min

! ===  Program  ======================================================================
!         Name:  test_multiply_edges
!  Description:  test c = min(a+b)  functionality for the z-matrices
! ====================================================================================
program test_multiply_edges
use zmatrix
use multiplication
use ieee_arithmetic
implicit none

real(kind=dp),dimension(:,:),allocatable :: a,b,c,resc
type(z_matrix) :: za,zb,zc
integer :: i,j,n,bl

n = 1
bl = 1

do while ( n .le. 1024 )

    do while (bl .le. n)

        allocate(a(n,n),b(n,n),c(n,n))

        call random_number(a)
        call random_number(b)
        c = 0.0_DP
    
        call matrix_to_zmatrix(a,za,bl)
        call matrix_to_zmatrix(b,zb,bl)
        call matrix_to_zmatrix(c,zc,bl) ! useless

        call zmatrix_multiply_edges(za,zb,zc)
        call three_loops_e(a,b,c)

        call zmatrix_to_matrix(zc,resc)

        do i = 1,n
            do j = 1,n
                if((abs(resc(i,j)-c(i,j))) .gt. 0.000001_dp) then
                    print *," error :: ",i,j," - ",resc(i,j), "!=",c(i,j)
                endif
            enddo
        enddo
    
        deallocate(a,b,c,za%matrix,zb%matrix,zc%matrix,resc)
        bl = bl *2
    end do
    n = n * 2
end do

end program test_multiply_edges
