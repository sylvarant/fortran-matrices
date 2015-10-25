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
!         Name:  testperformance_edges 
!  Description:  test the performance of the edge multiplication computations
! ====================================================================================

program testperformance_edges
use zmatrix
use multiplication
implicit none

integer,parameter :: maxit = 20


integer :: i,it,n,bl,b
real :: t1,t2
real(kind=dp),allocatable,dimension(:,:) :: g,r
type(z_matrix) :: zg,zr
real(kind=dp),dimension(maxit) :: three
real(kind=dp),allocatable,dimension(:,:) :: edges

print *,"Give power i :: "
read *,i

n = 2 ** i
allocate(g(n,n),r(n,n),edges(maxit,i+1))
three = 0
edges = 0

do it = 1,maxit
    r = 0
    call random_number(g)
    call cpu_time(t1)
    call three_loops_e(g,g,r)
    call cpu_time(t2)
    three(it) = (t2-t1)

    b = 0
    do while (b .le. i)
        r = 0
        bl = 2 ** b
        call matrix_to_zmatrix(g,zg,bl)
        call matrix_to_zmatrix(r,zr,bl)
        call cpu_time(t1)
        call zmatrix_multiply_edges(zg,zg,zr)
        call cpu_time(t2)
        b = b + 1 
        deallocate(zg%matrix,zr%matrix)
        edges(it,b) = (t2 -t1)
    enddo
enddo

print *,three
call write_matrix(edges)

end program testperformance_edges
