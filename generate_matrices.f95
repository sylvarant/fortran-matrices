! Adriaan Larmuseau - s0190593 - 1ste Master CW
! ## Architecture :: 
! (uname -r) Linux VPCEA2S1E 2.6.35-25-generic #44-Ubuntu SMP Fri Jan 21 17:40:44 UTC 2011 x86_64 GNU/Linux
! ## Compilers :: 
! G95 (GCC 4.1.2 (g95 0.93!) Jun 16 2010)
! ## Lib ::
! blas -- version 1.2-7
! ## Compiler Commando ::
! zie makefile : generate
! ## Bestede tijd ::
! 10min

! ===  Program  ======================================================================
!         Name:  generate_matrices 
!  Description:  generate matrices to test with
! ====================================================================================
program generate_matrices
use zmatrix
implicit none

integer :: i,k
character(len=80) :: filenameA,filenameB
real(kind=dp),dimension(:,:),allocatable :: a,b

print *, "give  power of 2 "
read *,i

allocate(a(2**i,2**i),b(2**i,2**i))
call random_number(a)
call random_number(b)
write(filenameA,'(a,"_",i0,".data")') "A",2**i
write(filenameB,'(a,"_",i0,".data")') "B",2**i
open(13,file=filenameA)
open(14,file=filenameB)
call write_matrix(a,13)
call write_matrix(b,14)
i=i+1
deallocate(a,b)
close(13)
close(14)

end program generate_matrices
