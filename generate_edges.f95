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
!         Name:  generate_edges  
!  Description:  generate edge matrices to test with
! ====================================================================================
program generate_edges
use zmatrix
use ieee_arithmetic
implicit none

integer :: k,p,i,j
real(kind=dp) :: perc,x
character(len=80) :: filename
real(kind=dp),dimension(:,:),allocatable :: a

print *, "dimension : 2 ^ "
read *,k
print *, "percentage of connectivity (1000 = all nodes are connected)"
read *,p
perc = (real(p,kind=dp)/1000_DP)

call random_seed()
allocate(a(2**k,2**k))

do i = 1,2**k
    do j = 1,2**k
        if ( i .eq. j) then
            a(i,j) = 0
        else
           call random_number(x)
           if (x .le. perc) then
                call random_number(x)
                a(i,j) = x
           else
                a(i,j) = ieee_value(a(i,j),ieee_positive_inf) 
           endif
        endif
    enddo
enddo


write(filename,'(a,"_",i0,"_",i0,".data")') "G",2**k,p
open(13,file=filename)
call write_matrix(a,13)


deallocate(a)

end program generate_edges
