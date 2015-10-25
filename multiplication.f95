! Adriaan Larmuseau - s0190593 - 1ste Master CW
! ## Architecture :: 
! (uname -r) Linux VPCEA2S1E 2.6.35-25-generic #44-Ubuntu SMP Fri Jan 21 17:40:44 UTC 2011 x86_64 GNU/Linux
! ## Compilers :: 
! G95 (GCC 4.1.2 (g95 0.93!) Jun 16 2010)
! ## Lib ::
! blas -- version 1.2-7
! ## Compiler Commando ::
! zie makefile
! ## Bestede tijd ::
! 1 uur

! ===  MODULE  ======================================================================
!         Name:  multiplication
!  Description:  contains additional methods to perform matrix multiplication
! ===================================================================================

module multiplication
use zmatrix
implicit none
private

!---------------------------- 
! Public methods
!----------------------------
public :: three_loops
public :: dgemm_blas
public :: std_matmul
public :: three_loops_e
contains 

!===============================================
! 3 loops matrix multiplication - kji style as
! it was the fastest in my previous assignment
!-----------------------------------------------
! Arguments ::
!   => a -- the left side matrix
!
!   => b -- the right side matrix
!
!   => c -- the result of a * b
!===============================================
subroutine three_loops(a,b,c)
    real(kind=DP), dimension(:,:), intent(out) :: c
    real(kind=DP), dimension(:,:), intent(in)  :: a, b
    integer :: i,k,j
    c = 0.0_DP
    do k = 1, size(c,2)
        do j = 1, size(c,1)
            do i = 1,size(a,2)
                c(i,k) = c(i,k) + (a(i,j) * b(j,k))
            end do
        end do
    end do
end subroutine three_loops


!===============================================
! 3 loops matrix multiplication -for edge matrix
! not meant to be efficient, just correct
!-----------------------------------------------
! Arguments ::
!   => a -- the left side matrix
!
!   => b -- the right side matrix
!
!   => c -- the result of min(a+b)
!===============================================
subroutine three_loops_e(a,b,c)
    real(kind=DP), dimension(:,:), intent(out) :: c
    real(kind=DP), dimension(:,:), intent(in)  :: a, b

    integer :: i,k,j
    real(kind=DP),dimension(size(a,2)) :: scal

    c = 0.0_DP
    do i = 1, size(c,1)
        do j = 1, size(c,2)
            do k = 1,size(a,2)
                scal(k) = a(i,k) + b(k,j)
            enddo
                c(i,j) = minval(scal)
        enddo
    enddo
end subroutine three_loops_e



!===============================================
! c = a * b using dgemm method from BLAS
!-----------------------------------------------
! Arguments ::
!   => a -- the left side matrix
!
!   => b -- the right side matrix
!
!   => c -- the result of dgemm(a,b)
!===============================================
subroutine dgemm_blas(a,b,c)
    real(kind=DP), dimension(:,:), intent(out) :: c
    real(kind=DP), dimension(:,:), intent(in)  :: a, b
    c = 0.0_DP ! initialise c
    call dgemm('n','n',size(a,1),size(b,2),size(a,2),1.0_DP,a,&
               size(a,2),b,size(a,2),0.0_DP,c,size(a,1)) 
end subroutine dgemm_blas


!===============================================
! c = a*b using built-in matul
!-----------------------------------------------
! Arguments ::
!   => a -- the left side matrix
!
!   => b -- the right side matrix
!
!   => c -- the result of matmul(a,b)
!===============================================
subroutine std_matmul(a,b,c)
    real(kind=DP), dimension(:,:), intent(out) :: c
    real(kind=DP), dimension(:,:), intent(in)  :: a, b
    c = matmul(a,b) 
end subroutine std_matmul

end module multiplication
