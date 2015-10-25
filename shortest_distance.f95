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
!         Name:  shortest_distance
!  Description:  calculate shortest distance for matrix G
! ====================================================================================
program shortest_distance
use zmatrix
implicit none

!---------------------------- 
! Constants
!----------------------------
! maximum command line argument length
integer,parameter :: ARG_LEN = 128

!---------------------------- 
! Variables
!----------------------------
integer :: i,blocksize,upperbound
real :: time_result,t2,t1
logical :: normal_result,distance,matrix_g,time_steps
character(len=ARG_LEN) :: arg
real(kind=DP),dimension(:,:),allocatable :: g,r,s
type(z_matrix) :: zg,zr,stepsmatrix

! initialisation
i = 0
upperbound = 0
normal_result = .true.
time_steps = .false.
distance = .false.
matrix_g = .true.
blocksize = 0
time_result = 0


!---------------------------- 
! parse commandline commands
!----------------------------
do while (i .lt. command_argument_count()) 

    i = i+1
    call get_command_argument(i,arg)

    select case(arg)

        case("--blocksize")
            call get_command_argument(i+1,arg)
            read (arg,*)  blocksize
            i = i+1

        case("--maxsteps")
            call get_command_argument(i+1,arg)
            read (arg,*)  upperbound
            i = i+1

        case("--steps")
            distance = .true. 

        case("--time")
            normal_result = .false.

        case("--timesteps")
            normal_result = .false.
            time_steps = .true.

        case("--help")
            call info(.false.)

        case("--version")
            print *, "Version 1.0 -- by Adriaan Larmuseau"
            stop

        ! parse file name of g matrix
        case default
            if(arg .eq. "-") then
                if(matrix_g) then
                    call read_matrix(g)
                    matrix_g = .false.
                else
                    print *,"ignored command :: ",arg
                endif
            else
                if(matrix_g) then
                    open (13,file=arg,err=100)
                    call read_matrix(g,13)
                    close(13)
                    matrix_g = .false.
                else
                    print *,"ignored command :: ",arg
                endif
            endif

    end select

end do

!---------------------------- 
! execute
!----------------------------
if(allocated(g)) then

    ! set blocksize
    if(blocksize .eq. 0) then
        ! default blocksize
        if(size(g,1) .le. BLOCK_SIZE) then
            blocksize = real(size(g,1)) / real(4) 
            if(blocksize .eq. 0) then
            blocksize = blocksize + 1
            endif
        else
         blocksize = BLOCK_SIZE
        endif
    elseif(blocksize .gt. size(g,1)) then
        print *,"Chosen block size ",blocksize," is larger than the matrix dim ",size(g,1) 
        stop 1
    endif
    
    ! set up variables
    call matrix_to_zmatrix(g,zg,blocksize)
    zr = zg

    if(distance) then
        stepsmatrix = zr
        if(upperbound .gt. 0) then
            call cpu_time(t1)
            call find_distance(zg,zr,upperbound,stepsmatrix,time=time_result)
            call cpu_time(t2)
        else
            call cpu_time(t1)
            call find_distance(zg,zr,steps=stepsmatrix,time=time_result)
            call cpu_time(t2)
        end if
    else
        if(upperbound .gt. 0) then
            call cpu_time(t1)
            call find_distance(zg,zr,l=upperbound,time=time_result)
            call cpu_time(t2)
        else
            call cpu_time(t1)
            call find_distance(zg,zr,time=time_result)
            call cpu_time(t2)
        end if   
    end if

    ! print results
    if(normal_result) then
        if (distance) then
            call zmatrix_to_matrix(stepsmatrix,s)
            call write_matrix(s) 
        else
            call zmatrix_to_matrix(zr,r)
            call write_matrix(r)
        endif
    else
        if(time_steps) then
            write(6,*) (t2-t1)
        else
            write(6,*) time_result
        endif
    endif

endif

! stop program
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

    print *, "This program demands a matrix G as first parameter and"
    print *, "calculates the shortest distance"
    print *, "this can be replaced by stdio using - in it's place"
    print *, "extra options are ::"
    print *, "  −−blocksize n : set the zmatrix blocksize "
    print *, "  −−stepsmatrix : replaces the distances by the number of stepsmatrix"
    print *, "  −−maxstepsmatrix : sets a limit to the number of stepsmatrix"
    print *, "  −−time : print out time needed to calculate the matrix product"
    print *, "  −−help , −−version : this menu and version number "

    if(present(error) .and. (.not. error)) then
        stop 
    else
        stop 1
    endif
end subroutine info

end program shortest_distance
