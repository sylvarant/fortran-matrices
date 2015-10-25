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
!         Name:  test_read
!  Description:  test the reading and writing of matrices
! ====================================================================================
program test_read
use zmatrix
implicit none

integer,parameter :: maxk = 10
integer :: i,p,q
character(len=80) :: filenameA,filenameB
real(kind=dp),dimension(:,:),allocatable :: a,b

! initialisation
i = 0

do while (i .le. maxk)
    allocate(a(2**i,2**i))
    call random_number(a)
    write(filenameA,'(a,"_",i0,".data")') "temp_",2**i
    open(13,file=filenameA)
    call write_matrix(a,13)
    close(13)
    open(13,file=filenameA)
    call read_matrix(b,13)
    do p=1,2**i
        do q = 1,2**i
            if((abs(a(p,q)-b(p,q))) .gt. 0.000001_DP) then
            print *," error :: ",p,q," - ",a(p,q), "!=",b(p,q)
            endif
        end do
    end do
    ! go to next
    i=i+1
    deallocate(a,b)
end do

end program test_read
