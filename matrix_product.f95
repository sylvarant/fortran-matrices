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
! 3u

! ===  Program  ======================================================================
!         Name:  matrix_product
!  Description:  performs the functionality requested of matrix_product as stated in 
!  the first phase of the assignment 
! ====================================================================================
program matrix_product
use zmatrix
use multiplication
implicit none


!---------------------------- 
! Constants
!----------------------------
! maximum command line argument length
integer,parameter :: ARG_LEN = 128


!---------------------------- 
! Variables
!----------------------------
integer :: i,blocksize,method
real :: time_result,t1,t2
logical :: normal_result,matrix_a,matrix_b
character(len=ARG_LEN) :: arg
real(kind=DP),dimension(:,:),allocatable :: a,b,c
type(z_matrix) :: za,zb,zc


! initialisation
normal_result = .true.
method = 1 ! default method is z-matrix
i = 0
matrix_a = .true.
matrix_b = .true.
blocksize = 0


!---------------------------- 
! parse commandline commands
!----------------------------
do while (i .lt. command_argument_count()) 
    
    i = i+1
    call get_command_argument(i, arg)

    select case(arg)
    
        case("--3loops")
            method = 2

        case("--zmatrix")
            method = 1
        
        case("--matmul")
            method = 3
        
        case("--dgemm")
            method = 4
                
        case("--blocksize")
            call get_command_argument(i+1,arg)
            read (arg,*)  blocksize
            i = i+1

        case("--time")
            normal_result = .false.

        case("--help")
            call info(.false.)

        case("--version")
            print *, "Version 1.0 -- by Adriaan Larmuseau"
            stop

        ! read in matrix files other wise
        case default
            if(arg .eq. "-") then
                if(matrix_a) then
                    call read_matrix(a)
                    matrix_a = .false.
                else
                    if(matrix_b) then
                        call read_matrix(b)
                        matrix_b = .false.
                    else
                        print *,"ignored command :: ",arg
                    endif
                endif
            else
                open (13,file=arg,err=100)
                if(matrix_a) then
                    call read_matrix(a,13)
                    matrix_a = .false.
                else
                    if(matrix_b) then
                        call read_matrix(b,13)
                        matrix_b = .false.
                    else
                        print *,"ignored command :: ",arg
                    endif
                endif
                close(13)
            endif
            
    end select

end do

!---------------------------- 
! execute 
!----------------------------
if(allocated(a) .and. allocated(b) )then

    if(size(a,1) .ne. size(b,1)) then
        print *,"A & B do not have the same dimensions !!"
        stop 1
    endif

    allocate(c(size(a,1),size(a,1)))

    ! set blocksize
    if(blocksize .eq. 0) then
        ! default blocksize
        if(size(a,1) .le. BLOCK_SIZE) then
            blocksize = real(size(a,1)) / real(4) 
            if(blocksize .eq. 0) then
            blocksize = blocksize + 1
            endif
        else
         blocksize = BLOCK_SIZE
        endif
    elseif(blocksize .gt. size(a,1)) then
        print *,"Chosen block size ",blocksize," is larger than the matrix dim ",size(a,1) 
        stop 1
    endif

    ! parse the different possible methods
    select case(method)
    
        case(1) ! zmatrix
            call matrix_to_zmatrix(a,za,blocksize)
            call matrix_to_zmatrix(b,zb,blocksize)
            call matrix_to_zmatrix(c,zc,blocksize) ! to lazy to set by hand
            call cpu_time(t1)
            call zmatrix_multiply(za,zb,zc)
            call cpu_time(t2)
            deallocate(c)
            call zmatrix_to_matrix(zc,c) 

        case(2) ! 3loops
            call cpu_time(t1)
            call three_loops(a,b,c)
            call cpu_time(t2)

        case(3) ! matmum
            call cpu_time(t1)
            call std_matmul(a,b,c)
            call cpu_time(t2)

        case(4) ! dgemm
            call cpu_time(t1)
            call dgemm_blas(a,b,c)
            call cpu_time(t2)

    end select

    ! print result or time if requested
    if (normal_result) then
        call write_matrix(c)
    else
        time_result = (t2-t1)
        write(6,'(f12.7)') time_result
    endif

    ! free memory
    deallocate(a,b,c)

endif

! end program
stop 

! handle file opening errors
100 call info 


contains


! ===  Subroutine  ===================================================================
!         Name:  info
!  Description:  print info on how to use this program
!  Parameters:  
!       * error  -- ending in error ?       -- INPUT    -- OPTIONAL
! ====================================================================================
subroutine info(error)
    logical,optional,intent(in) :: error

    print *, "This program demands two text files for the A & B parameter and calculates the product"
    print *, "One of these can be replaced by stdio using - in it's place"
    print *, "extra options are ::"
    print *, "  −−3loops, −−zmatrix, −−matmul & −−dgemm : select multiplication method"
    print *, "  −−blocksize n : set the zmatrix blocksize "
    print *, "  −−time : print out time needed to calculate the product instead of C"
    print *, "  −−help , −−version : this menu and version number "

    if(present(error) .and. (.not. error)) then
        stop 
    else
        stop 1
    endif
end subroutine info

end program matrix_product
