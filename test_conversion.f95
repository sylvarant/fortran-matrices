! Adriaan Larmuseau - s0190593 - 1ste Master CW
! ## Architecture :: 
! (uname -r) Linux VPCEA2S1E 2.6.35-25-generic #44-Ubuntu SMP Fri Jan 21 17:40:44 UTC 2011 x46_64 GNU/Linux
! ## Compilers :: 
! G95 (GCC 4,1,2 (g95 0.93!) Jun 16 2010)
! ## Lib ::
! blas -- version 1.2-7
! ## Compiler Commando ::
! zie makefile : generate
! ## Bestede tijd ::
! 20min

! ===  Program  ======================================================================
!         Name:  test_conversion
!  Description:  test the conversions of the matrices
! ====================================================================================
program test_conversion
use zmatrix
implicit none

integer :: i,j,tel
real(kind=dp),dimension(:,:),allocatable :: a,resulta
real(kind=dp),dimension(8,8) :: sha,expected
type(z_matrix) ::  za

! create number sequence matrix
tel = 0
allocate(a(8,8))
do i = 1,8
    do j = 1,8
        tel = tel+1
        a(i,j) = tel
    end do
end do

call matrix_to_zmatrix(a,za,4)
sha = (reshape(za%matrix ,(/8,8/)))

! what we expect
expected = reshape((/ 1, 2, 3, 4, 9, 10, 11, 12, 17, 18, 19, 20, 25, 26, 27, 28,& 
5, 6, 7, 8, 13, 14, 15, 16, 21, 22, 23, 24, 29, 30, 31, 32, 33, 34, 35, 36, 41, 42,& 
43, 44, 49, 50, 51, 52, 57, 58, 59, 60, 37, 38, 39, 40, 45, 46, 47, 48, 53, 54, 55, 56,&
61, 62, 63, 64 /),(/8,8/))

! what we get
do i = 1,8
    do j = 1,8
        if((abs(expected(i,j)-sha(i,j))) .gt. 0.000001) then
            print *,"Error, incorrect z sequencing !!  :: ",i,j," - ",expected(i,j), "!=",sha(i,j)
        endif
    end do
end do

call zmatrix_to_matrix(za,resulta)

! test deconversion correctness
do i = 1,8
    do j = 1,8
        if((abs(resulta(i,j)-a(i,j))) .gt. 0.000001) then
            print *,"Error, incorrect deconversion :: ",i,j," - ",resulta(i,j), "!=",a(i,j)
        endif
    end do
end do

end program test_conversion


