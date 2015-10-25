! Adriaan Larmuseau - s0190593 - 1ste Master CW
! ## Architecture :: 
! (uname -r) Linux VPCEA2S1E 2.6.35-25-generic #44-Ubuntu SMP Fri Jan 21 17:40:44 UTC 2011 x86_64 GNU/Linux
! ## Compilers :: ! z-matrix uses double precision storing
! G95 (GCC 4.1.2 (g95 0.93!) Jun 16 2010)
! ## Lib ::
! blas -- version 1.2-7
! ## Compiler Commando ::
! zie makefile
! ## Bestede tijd ::
! 16 - 18 uur

! ===  MODULE  ======================================================================
!         Name:  zmatrix
!  Description:  contains all functionality related to zmatrices
! ===================================================================================
module zmatrix
use ieee_arithmetic
implicit none
private

!---------------------------- 
! Constants
!----------------------------
! z-matrix uses double precision storing
integer, public, parameter :: DP = selected_real_kind(15,307)
! use , as DELimiter for elements
character,public, parameter :: DEL = " "
! default format for writing the numbers
character(len=*),public,parameter :: DEF_FRMT = '(ES11.5)'
! maximum numbers of chars per line
integer,public,parameter :: LINE_W = 1000
! maximum numbers of chars for the format , as all strings are not well supported
integer,public,parameter :: DEF_FRMT_W = 50
! default blocksize chose by me
integer,public,parameter :: BLOCK_SIZE = 64


!---------------------------- 
! Z-Matrix type
!----------------------------
type z_matrix
   real(kind=DP),allocatable,dimension(:)   :: matrix
   integer                                  :: n
   integer                                  :: blockdim
end type
public z_matrix

!---------------------------- 
! Public methods
!----------------------------
public :: matrix_to_zmatrix
public :: zmatrix_to_matrix
public :: zmatrix_multiply
public :: zmatrix_multiply_edges
public :: find_distance
public :: read_matrix
public :: write_matrix


contains


! ===  Subroutine  ===================================================================
!         Name:  fillz
!  Description:  copy form dest to src row-row  zmatrix style
!  Parameters:  
!       * dest      -- the destination must be nxn      -- INPUT
!       * src       -- the source zmatrix               -- OUTPUT
!       * blocksize -- the lowest level blocksize       -- INPUT    -- OPTIONAL
! ====================================================================================
recursive subroutine fillz(dest,src,blocksize)
    integer,optional,intent(in) :: blocksize
    real(kind=DP),intent(in),dimension(:,:) :: src
    real(kind=DP),intent(out),dimension(size(src,1)**2) :: dest
    
    integer :: i,j,n,ncube
    n = size(src,1)
    ncube = n**2
    if(n .eq. blocksize) then
        ! load in the block row by row
        do j = 1,n
            do i = 1,n
                dest(j + ((i-1)*n)) = src(i,j)
            end do
        end do
    else
        ! split up in 4
        call fillz(dest( : ncube/4),src( : (n/2), : (n/2)),blocksize)
        call fillz(dest((ncube/4)+1 : (ncube/2)),src( : (n/2),(n/2)+1 : ),blocksize)
        call fillz(dest((ncube/2)+1 : (ncube/4)*3),src((n/2)+1 : , : (n/2)),blocksize) 
        call fillz(dest(((ncube/4)*3)+1 : ),src((n/2)+1 :,(n/2)+1 :),blocksize)
    endif
end subroutine fillz


! ===  Subroutine  ===================================================================
!         Name:  fill
!  Description:  copy form dest to src normal matrix style
!  Parameters:  
!       * dest      -- the destination must be nxn      -- INPUT
!       * src       -- the source zmatrix               -- OUTPUT
!       * blocksize -- the lowest level blocksize       -- INPUT    -- OPTIONAL
! ====================================================================================
recursive subroutine fill(dest,src,blocksize)
    integer,optional,intent(in) :: blocksize
    real(kind=DP),intent(in),dimension(:) :: src
    real(kind=DP),intent(out),dimension(:,:) :: dest
    
    integer :: i,j,n,ncube
    n = size(dest,1)
    ncube = size(src)

    if(n .eq. blocksize) then
        do j = 1,n
            do i = 1,n
                dest(i,j) = src(j + ((i-1) *n))
            end do
        end do
    else
        ! split up in 4
        call fill(dest( : (n/2), : (n/2)),src( : ncube/4),blocksize)
        call fill(dest( : (n/2),(n/2)+1 : ),src((ncube/4)+1 : (ncube/2)),blocksize)
        call fill(dest((n/2)+1 : , : (n/2)),src((ncube/2)+1 : (ncube/4)*3),blocksize) 
        call fill(dest((n/2)+1 :,(n/2)+1 :),src(((ncube/4)*3)+1 : ),blocksize)
    endif
end subroutine fill


! ===  Subroutine  ===================================================================
!         Name:  matrix_to_zmatrix
!  Description:  convert a matrix to a z_matrix
!  Parameters:  
!       * arg       -- the to be converted matrix       -- INPUT
!       * res       -- the resulting zmatrix            -- OUTPUT
!       * blocksize -- the lowest level blocksize       -- INPUT    -- OPTIONAL
! ====================================================================================
subroutine matrix_to_zmatrix(arg,res,blocksize) 
    type(z_matrix),intent(out) :: res
    integer,optional,intent(in) :: blocksize
    real(kind=DP),intent(in),dimension(:,:) :: arg

    ! create the zmatrix
    res%n = size(arg,1)
    if(present(blocksize)) then
        res%blockdim = blocksize
    else
        res%blockdim = 1
    endif
    allocate(res%matrix(res%n*res%n))
    call fillz(res%matrix,arg,res%blockdim) 
end subroutine matrix_to_zmatrix


! ===  Subroutine  ===================================================================
!         Name:  zmatrix_to_matrix
!  Description:  deconverts a zmatrix to a matrix
!  Parameters:  
!       * arg     -- the zmatrix to deconvert           -- INPUT
!       * res     -- the resulting matrix               -- OUTPUT
! ====================================================================================
subroutine zmatrix_to_matrix(arg,res) 
    type(z_matrix),intent(in) :: arg
    real(kind=DP),dimension(:,:),allocatable,intent(out) :: res

    ! creat the matrix
    allocate(res(arg%n,arg%n))

    call fill(res,arg%matrix,arg%blockdim)
end subroutine zmatrix_to_matrix


! ===  Subroutine  ===================================================================
!         Name:  matrix_multiply
!  Description:  for the low level matrices coming from the z-matrices (row by row !)
!  Parameters:  
!       * a     -- the A matrix                     -- INPUT
!       * b     -- the B matrix                     -- INPUT
!       * c     -- the resulting C matrix           -- INPUT/OUTPUT
!       * n     -- the block size                   -- INPUT
!       * alpha -- the alpha value (default = 1)    -- INTPUT       -- OPTIONAL
!       * beta  -- the beta value  (default = 0)    -- INPUT        -- OPTIONAL
! ====================================================================================
subroutine matrix_multiply(a,b,c,alpha,beta,n)
    integer,intent(in):: n
    real(kind=DP),intent(in),optional :: alpha,beta
    real(kind=DP),intent(in),dimension(:) :: a
    real(kind=DP),intent(in),dimension(size(a)) :: b
    real(kind=DP),intent(inout),dimension(size(a)) :: c

    integer :: row,col,i
    real(kind=DP) :: prod
    real(kind=DP) :: al,be

    !---------------------------- 
    ! handle optional paramters
    !----------------------------
    if( present(alpha) ) then
        al = alpha
    else
        al = 1.0_DP
    endif

    if (present(beta)) then
        be = beta
    else
        be = 0.0_DP;
    endif

    !---------------------------- 
    ! calculate C = al*AB + be*C
    !----------------------------
    do row=1,n
        do col=1,n
            prod = 0
            do i=1,n
                prod = prod + (b(col+((i-1) * n)) * a(i+((row -1)*n)))
            end do
            c(col+((row-1)*n)) = (al * prod) + (be * c(col+((row-1)*n)))
        end do
    end do
end subroutine matrix_multiply


! ===  Subroutine  ===================================================================
!         Name:  split
!  Description:  split up the matrix into small blocks to perform multiplication
!  Parameters:  
!       * a     -- the A matrix                     -- INPUT
!       * b     -- the B matrix                     -- INPUT
!       * c     -- the resulting C matrix           -- INPUT/OUTPUT
!       * n     -- the block size                   -- INPUT
!       * alpha -- the alpha value                  -- INTPUT      
! ====================================================================================
recursive subroutine split(a,b,c,alpha,n)
    integer,intent(in):: n
    real(kind=DP),intent(in) :: alpha
    real(kind=DP),intent(in),dimension(:) :: a
    real(kind=DP),intent(in),dimension(size(a)) :: b
    real(kind=DP),intent(inout),dimension(size(a)) :: c

    integer :: l
    l = size(a)

    if(l .eq. n**2) then
        ! multiply blocks beta = 1 !! to do c +=
        call matrix_multiply(a,b,c,alpha,1.0_DP,n) 
    else
        ! top left
        call split(a( :(l/4)), b( :(l/4)),c(:(l/4)),alpha,n) ! a11 * b11 += c11
        call split(a((l/4)+1 :(l/2)), b((l/2)+1:(l/4)*3),c(:(l/4)), alpha,n) ! a12 * b21 += c11
        ! top right
        call split(a( :(l/4)),b((l/4)+1 : (l/2)) ,c((l/4)+1:(l/2)),alpha,n) ! a11*b12 += c12
        call split(a((l/4)+1 : (l/2)),b(((l/4)*3)+1: ),c((l/4)+1:(l/2)),alpha,n) ! a12 * b22 +=c12
        ! bottom left
        call split(a((l/2)+1 : (l/4)*3 ),b(:(l/4)),c((l/2)+1:(l/4)*3),alpha,n) ! a21 * b11+=c21
        call split(a(((l/4)*3)+1 :),b((l/2)+1:(l/4)*3),c((l/2)+1:(l/4)*3),alpha,n) ! a22 * b21+=c21
        ! bottom right
        call split(a((l/2)+1 :(l/4)*3),b((l/4)+1 : (l/2)),c(((l/4)*3)+1 : ),alpha,n) ! a21 * b12+= c22
        call split(a(((l/4)*3)+1:),b(((l/4)*3)+1: ),c(((l/4)*3)+1 :),alpha,n)  ! a22 * b22+=c22
    endif
end subroutine split


! ===  Subroutine  ===================================================================
!         Name:  zmatrix_multiply
!  Description:  perform C = alpha * (Ax*B) + beta * C
!  Parameters:  
!       * a     -- the A matrix                     -- INPUT
!       * b     -- the B matrix                     -- INPUT
!       * c     -- the C matrix                     -- INPUT/OUTPUT
!       * alpha -- the alpha value (default = 1)    -- INTPUT       -- OPTIONAL
!       * beta  -- the beta value  (default = 0)    -- INPUT        -- OPTIONAL
! ====================================================================================
subroutine zmatrix_multiply(a,b,c,alpha,beta)
    real(kind=DP),intent(in),optional :: alpha,beta
    type(z_matrix),intent(in) :: a,b
    type(z_matrix),intent(inout) :: c

    real(kind=DP) :: al

    !---------------------------- 
    ! handle optional paramters
    !----------------------------
    if( present(alpha) ) then
        al = alpha
    else
        al = 1.0_DP
    endif

    if(present(beta) ) then
        c%matrix = beta * c%matrix
    else
        c%matrix = 0.0_DP * c%matrix
    endif
    
    ! call split 
    call split(a%matrix,b%matrix,c%matrix,al,a%blockdim)
end subroutine zmatrix_multiply


! ===  Subroutine  ===================================================================
!         Name:  matrix_multiply_edges
!  Description:  for the low level matrices coming from the z-matrices (row by row !)
!  Parameters:  
!       * a     -- the A matrix                     -- INPUT
!       * b     -- the B matrix                     -- INPUT
!       * c     -- the resulting C matrix           -- INPUT/OUTPUT
!       * n     -- the block size                   -- INPUT
!       * alpha -- the alpha value (default = 1)    -- INTPUT       -- OPTIONAL
!       * beta  -- the beta value  (default = 0)    -- INPUT        -- OPTIONAL
! ====================================================================================
subroutine matrix_multiply_edges(a,b,c,alpha,beta,n)
    integer,intent(in):: n
    real(kind=DP),intent(in),optional :: alpha,beta
    real(kind=DP),intent(in),dimension(:) :: a
    real(kind=DP),intent(in),dimension(size(a)) :: b
    real(kind=DP),intent(inout),dimension(size(a)) :: c

    integer :: row,col,k
    real(kind=DP):: scal_sum
    real(kind=DP) :: al,be

    !---------------------------- 
    ! handle optional paramters
    !----------------------------
    if( present(alpha) ) then
        al = alpha
    else
        al = 1.0_DP
    endif

    if (present(beta)) then
        be = beta
    else
        be = 0.0_DP;
    endif

    !---------------------------- 
    ! calculate C = al*AB + be*C
    ! 0 times infinity is a problem so we 
    ! split up ! todo al == 0
    !----------------------------
    if( be .eq. 0) then
        do row=1,n
            do col=1,n
                scal_sum = ieee_value(scal_sum,ieee_positive_inf)
                do k=1,n
                    scal_sum = min(scal_sum,(b(col+((k-1) * n)) + a(k+((row-1)*n))))
                end do
                c(col+((row-1)*n)) = (scal_sum * al)
            end do
        end do
    else
        do row=1,n
            do col=1,n
                scal_sum = ieee_value(scal_sum,ieee_positive_inf)
                do k=1,n
                    scal_sum = min(scal_sum,(b(col+((k-1) * n)) + a(k+((row-1)*n))))
                end do
                c(col+((row-1)*n)) = (scal_sum * al) + (be *c(col+((row-1)*n))) 
            end do
        end do
    endif
end subroutine matrix_multiply_edges


! ===  Subroutine  ===================================================================
!         Name:  split_edges
!  Description:  split up the matrix into small blocks to perform edge multiplication
!  Parameters:  
!       * a     -- the A matrix                     -- INPUT
!       * b     -- the B matrix                     -- INPUT
!       * c     -- the resulting C matrix           -- INPUT/OUTPUT
!       * n     -- the block size                   -- INPUT
!       * alpha -- the alpha value                  -- INTPUT    
!       * beta  -- the beta value                   -- INTPUT   
! ====================================================================================
recursive subroutine split_edges(a,b,c,alpha,beta,n)
    integer,intent(in):: n
    real(kind=DP),intent(in) :: alpha,beta
    real(kind=DP),intent(in),dimension(:) :: a
    real(kind=DP),intent(in),dimension(size(a)) :: b
    real(kind=DP),intent(inout),dimension(size(a)) :: c

    integer :: l
    real(kind=dp),dimension(size(a)) :: c2
    c2 = c
    l = size(a)
    if(l .eq. n**2) then
        ! to do c = min(a+b,c)
        call matrix_multiply_edges(a,b,c,alpha,beta,n) 
    else
        ! top left
        call split_edges(a( :(l/4)), b( :(l/4)),c(:(l/4)),alpha,beta,n) 
        call split_edges(a((l/4)+1 :(l/2)), b((l/2)+1:(l/4)*3),c2(:(l/4)),alpha,beta,n) 
        c(:(l/4)) = min(c(:(l/4)),c2(:(l/4)))

        ! top right
        call split_edges(a( :(l/4)),b((l/4)+1 : (l/2)) ,c((l/4)+1:(l/2)),alpha,beta,n) 
        call split_edges(a((l/4)+1 : (l/2)),b(((l/4)*3)+1: ),c2((l/4)+1:(l/2)),alpha,beta,n)
        c((l/4)+1:(l/2)) = min(c((l/4)+1:(l/2)),c2((l/4)+1:(l/2)))

        ! bottom left
        call split_edges(a((l/2)+1 : (l/4)*3 ),b(:(l/4)),c((l/2)+1:(l/4)*3),alpha,beta,n) 
        call split_edges(a(((l/4)*3)+1 :),b((l/2)+1:(l/4)*3),c2((l/2)+1:(l/4)*3),alpha,beta,n) 
        c((l/2)+1:(l/4)*3) = min(c((l/2)+1:(l/4)*3),c2((l/2)+1:(l/4)*3))

        ! bottom right
        call split_edges(a((l/2)+1 :(l/4)*3),b((l/4)+1 : (l/2)),c(((l/4)*3)+1 : ),alpha,beta,n)
        call split_edges(a(((l/4)*3)+1:),b(((l/4)*3)+1: ),c2(((l/4)*3)+1 :),alpha,beta,n) 
        c(((l/4)*3)+1:) = min(c(((l/4)*3)+1:),c2(((l/4)*3)+1:))
    endif
end subroutine split_edges


! ===  Subroutine  ===================================================================
!         Name:  zmatrix_multiply_edges
!  Description:  perform C = alpha * (AxB) + beta * C
!  Parameters:  
!       * a     -- the A matrix                     -- INPUT
!       * b     -- the B matrix                     -- INPUT
!       * c     -- the C matrix                     -- INPUT/OUTPUT
!       * alpha -- the alpha value (default = 1)    -- INTPUT       -- OPTIONAL
!       * beta  -- the beta value  (default = 0)    -- INPUT        -- OPTIONAL
! ====================================================================================
subroutine zmatrix_multiply_edges(a,b,c,alpha,beta)
    real(kind=DP),intent(in),optional :: alpha,beta
    type(z_matrix),intent(in) :: a,b
    type(z_matrix),intent(inout) :: c
    real(kind=DP) :: al,be

    !---------------------------- 
    ! handle optional paramters
    !----------------------------
    if( present(alpha) ) then
        al = alpha
    else
        al = 1.0_DP
    endif

    ! check beta
    if(present(beta)) then
        be = beta 
    else
        be = 0.0_DP
    endif

    ! call split 
    call split_edges(a%matrix,b%matrix,c%matrix,al,be,a%blockdim)
end subroutine zmatrix_multiply_edges


! ===  Subroutine  ===================================================================
!         Name:  find_distance
!  Description:  for a given G find the shortest distance between the nodes
!  Parameters:  
!       * zg    -- the input G matrix           -- INPUT
!       * zr    -- the resulting matrix         -- OUTPUT
!       * l     -- the maximum number of steps  -- INPUT     -- OPTIONAL
!       * steps -- a matrix with n of steps     -- OUTPUT    -- OPTIONAL
! ====================================================================================
subroutine find_distance(zg,zr,l,steps,time)
    type(z_matrix),intent(in) :: zg
    type(z_matrix),intent(in out) :: zr
    integer,intent(in),optional :: l
    real,intent(in out),optional :: time
    type(z_matrix),intent(in out),optional :: steps

    real :: t1,t2,time_res,h
    type(z_matrix) :: temp
    integer :: maximum,nsteps


    !---------------------------- 
    ! handle optional paramters
    !----------------------------
    if(present(l)) then
        maximum = (log(real(l))/log(2.0_DP)) 
    else
        maximum = (log(real(zg%n))/log(2.0_dp)) + 1
    endif

    !---------------------------- 
    ! initialisation
    !----------------------------
    nsteps = 1
    temp = zg
    time_res = 0
    if(present(steps)) then
        where((.not.( temp%matrix .eq. ieee_value(temp%matrix,ieee_positive_inf))) .and. & 
                (temp%matrix .gt. 0 )) temp%matrix = 1
    endif

    !---------------------------- 
    ! execution
    !----------------------------
    do while ((nsteps .le. maximum))
        call cpu_time(t1)
        call zmatrix_multiply_edges(temp,temp,zr)
        call cpu_time(t2)
        time_res = time_res + (t2-t1)
        nsteps = nsteps +1
        if (abs(sum(temp%matrix,MASK= temp%matrix .ne.ieee_value(temp%matrix,ieee_positive_inf)) - &
           sum(zr%matrix,MASK= zr%matrix .ne.ieee_value(zr%matrix,ieee_positive_inf))) .eq. 0 ) then
            exit
        endif
        temp = zr
    enddo

    if( present(time) .and. nsteps .gt. 1) then
        time = time_res / (real(nsteps)-real(1))
    endif

    if(present(steps)) then
        where (temp%matrix .eq. ieee_value(temp%matrix,ieee_positive_inf)) temp%matrix = 0
        steps = temp
    endif
end subroutine find_distance


! ===  Subroutine  ===================================================================
!         Name:  read_matrix
!  Description:  read in a matrix from an inputstream
!  Parameters:  
!       * res   -- the result goes in here          -- OUTPUT
!       * steps -- file descriptor (def = stdin)    -- INPUT    -- OPTIONAL
! ====================================================================================
subroutine read_matrix(res,unitn) 
    integer,optional,intent(in) :: unitn
    real(kind=DP),dimension(:,:),allocatable,intent(out) :: res

    integer :: u,n,m,i,pos2,pos1
    character(len=LINE_W) :: input


    !---------------------------- 
    ! handle optional paramters
    !----------------------------
    if( present(unitn) ) then
        u = unitn
    else
        u = 5 ! stdin
    endif


    !---------------------------- 
    ! read first line from file
    !----------------------------
    read (u,'(a)') input
    
    ! hack to fine find n x m 
    pos1 = index(input," x ")
    pos2 = index(input,"[")
    read (input((pos2+1):(pos1-1)),*) n
    pos2 = index(input,"]")
    read (input((pos1+3):(pos2-1)),*) m

    ! check given information about matrix is 2^k x 2^k
    if(n .ne. m .or. (iand(n,(n - 1)) .ne. 0) ) then
        print *,"!! not of the form 2^k x 2^k" ! warning
        stop 1
    endif

    ! allocate the necessary memory
    allocate(res(n,n))

    !---------------------------- 
    ! read matrix from file
    !----------------------------
    do i = 1,n
        read(u,*) res(i,:)
    end do
end subroutine read_matrix


! ===  Subroutine  ===================================================================
!         Name:  write_matrix
!  Description:  write a matrix to the output stream
!  Parameters:  
!       * arg   -- the matrix to write              -- INPUT
!       * unitn -- file descriptor (def = stdout)   -- INPUT    -- OPTIONAL
!       * cfmt  -- custom fmt to write              -- INPUT    -- OPTIONAL
! ====================================================================================
subroutine write_matrix(arg,unitn,cfmt)
    integer,optional,intent(in) :: unitn
    character(len=*) ,optional,intent(in) :: cfmt
    real(kind=DP),intent(in),dimension(:,:) :: arg

    character(len=DEF_FRMT_W) :: w_frm
    character(len=LINE_W) :: output
    integer :: u,i,j,pos,outputsize

    !---------------------------- 
    ! handle optional paramters
    !----------------------------
    if(present(unitn)) then
        u = unitn
    else
        u = 6 ! stdout
    endif

    if(present(cfmt)) then
        w_frm = cfmt
    else
        w_frm = DEF_FRMT
    endif


    !---------------------------- 
    ! print to file
    !----------------------------
    write (u,'("full matrix ","[",i0," x ",i0,"]:")') size(arg,1),size(arg,2)

    ! hack to figure out print size of custom format
    write(output,w_frm) atan(2.886)
    outputsize = len_trim(output) 

    ! print out every row of the matrix in lines of max thousand chars 
    do i = 1,size(arg,1)
        j = 0;pos = 1
        do
            j = j+ 1 
            write (u,w_frm,advance="no") arg(i,j) 
            pos = pos +outputsize + 1 

            if(pos .gt. (LINE_W - (outputsize +1))) then
                write(u,*)
                pos = 1
            else if ( j .eq. size(arg,2)) then
                write(u,*)
                exit ! quit loop
            else
                write (u,'(a)',advance="no") DEL
            endif
        end do
    end do
end subroutine write_matrix

end module zmatrix
