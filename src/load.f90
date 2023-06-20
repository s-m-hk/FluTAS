!
! SPDX-License-Identifier: MIT
!
module mod_load
  !
  use mpi
  use mod_common_mpi, only: ierr,myid,ipencil
  use mod_param     , only: ng
  use mod_types     , only: rp
  use decomp_2d
  use decomp_2d_io
  !
  implicit none
  !
  private
  public  :: load,load_int,load_weight,load_ns,load_scalar,loadIBM
  !
  contains
  !
  subroutine load(io,filename,n,fld)
    !
    ! reads/writes a restart file
    !
    implicit none
    !
    character(len=1), intent(in   )                            :: io
    character(len=*), intent(in   )                            :: filename
    integer         , intent(in   ), dimension(3)              :: n
    real(rp)        , intent(inout), dimension(n(1),n(2),n(3)) :: fld ! generic field to be read/written
    !
    integer(MPI_OFFSET_KIND) :: filesize,disp,good
    integer :: lenr,fh
    !
    select case(io)
    case('r')
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      lenr = storage_size(fld(1,1,1))/8
      !good = int(ng(1)*ng(2)*ng(3),MPI_OFFSET_KIND)*lenr
      good = int(lenr,MPI_OFFSET_KIND)*ng(1)*ng(2)*ng(3)
      !
      if(filesize.ne.good) then
        if(myid.eq.0) print*, ''
        if(myid.eq.0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
        if(myid.eq.0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call decomp_2d_finalize
        call MPI_FINALIZE(ierr)
        call exit
      endif
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_read_var(fh,disp,ipencil,fld )
      call MPI_FILE_CLOSE(fh,ierr)
      ! 
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_write_var(fh,disp,ipencil,fld)
      call MPI_FILE_CLOSE(fh,ierr)
      !
    end select
    !
  end subroutine load
  !
  subroutine load_int(io,filename,n,fld)
    !
    ! reads/writes a restart file
    !
    implicit none
    !
    character(len=1), intent(in   )                            :: io
    character(len=*), intent(in   )                            :: filename
    integer         , intent(in   ), dimension(3)              :: n
    integer       , intent(inout), dimension(n(1),n(2),n(3)) :: fld ! generic field to be read/written
    real(rp)      , dimension(n(1),n(2),n(3)) :: tmp
    !
    integer(MPI_OFFSET_KIND) :: filesize,disp,good
    integer :: lenr,fh
    !
    select case(io)
    case('r')
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      lenr = storage_size(tmp(1,1,1))/8
      !good = int(ng(1)*ng(2)*ng(3),MPI_OFFSET_KIND)*lenr
      good = int(lenr,MPI_OFFSET_KIND)*ng(1)*ng(2)*ng(3)
      !
      if(filesize.ne.good) then
        if(myid.eq.0) print*, ''
        if(myid.eq.0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
        if(myid.eq.0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call decomp_2d_finalize
        call MPI_FINALIZE(ierr)
        call exit
      endif
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_read_var(fh,disp,ipencil,tmp)
      call MPI_FILE_CLOSE(fh,ierr)
      fld = int(tmp)
      ! 
    case('w')
      !
      ! write
      !
      tmp = real(fld,rp)
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_write_var(fh,disp,ipencil,tmp)
      call MPI_FILE_CLOSE(fh,ierr)
      !
    end select
    !
  end subroutine load_int
  !
  subroutine load_weight(io,filename,n,fld)
    !
    ! reads/writes a restart file
    !
    implicit none
    !
    character(len=1), intent(in   )                            :: io
    character(len=*), intent(in   )                            :: filename
    integer         , intent(in   ), dimension(3)              :: n
    real(rp)        , intent(inout), dimension(n(1),n(2),n(3),7) :: fld ! generic field to be read/written
    real(rp)        , dimension(n(1),n(2),n(3)) :: tmp
    !
    integer(MPI_OFFSET_KIND) :: filesize,disp,good
    integer :: lenr,fh
    !
    select case(io)
    case('r')
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      lenr = storage_size(tmp(1,1,1))/8
      !good = int(ng(1)*ng(2)*ng(3),MPI_OFFSET_KIND)*lenr
      good = 7*int(lenr,MPI_OFFSET_KIND)*ng(1)*ng(2)*ng(3)
      !
      if(filesize.ne.good) then
        if(myid.eq.0) print*, ''
        if(myid.eq.0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
        if(myid.eq.0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call decomp_2d_finalize
        call MPI_FINALIZE(ierr)
        call exit
      endif
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_read_var(fh,disp,ipencil,tmp)
      fld(:,:,:,1) = tmp(:,:,:)
      call decomp_2d_read_var(fh,disp,ipencil,tmp)
      fld(:,:,:,2) = tmp(:,:,:)
      call decomp_2d_read_var(fh,disp,ipencil,tmp)
      fld(:,:,:,3) = tmp(:,:,:)
      call decomp_2d_read_var(fh,disp,ipencil,tmp)
      fld(:,:,:,4) = tmp(:,:,:)
      call decomp_2d_read_var(fh,disp,ipencil,tmp)
      fld(:,:,:,5) = tmp(:,:,:)
      call decomp_2d_read_var(fh,disp,ipencil,tmp)
      fld(:,:,:,6) = tmp(:,:,:)
      call decomp_2d_read_var(fh,disp,ipencil,tmp)
      fld(:,:,:,7) = tmp(:,:,:)
      call MPI_FILE_CLOSE(fh,ierr)
      ! 
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
      disp = 0_MPI_OFFSET_KIND
      tmp(:,:,:) = fld(:,:,:,1)
      call decomp_2d_write_var(fh,disp,ipencil,tmp)
      tmp(:,:,:) = fld(:,:,:,2)
      call decomp_2d_write_var(fh,disp,ipencil,tmp)
      tmp(:,:,:) = fld(:,:,:,3)
      call decomp_2d_write_var(fh,disp,ipencil,tmp)
      tmp(:,:,:) = fld(:,:,:,4)
      call decomp_2d_write_var(fh,disp,ipencil,tmp)
      tmp(:,:,:) = fld(:,:,:,5)
      call decomp_2d_write_var(fh,disp,ipencil,tmp)
      tmp(:,:,:) = fld(:,:,:,6)
      call decomp_2d_write_var(fh,disp,ipencil,tmp)
      tmp(:,:,:) = fld(:,:,:,7)
      call decomp_2d_write_var(fh,disp,ipencil,tmp)
      call MPI_FILE_CLOSE(fh,ierr)
      !
    end select
    !
  end subroutine load_weight
  !
  subroutine load_scalar(io,filename,time,istep,dto)
    !
    implicit none
    !
    character(len=1), intent(in   ) :: io
    character(len=*), intent(in   ) :: filename
    integer         , intent(inout) :: istep
    real(rp)        , intent(inout) :: time, dto
    !
    integer :: fh
    !
    select case(io)
    case('r')
      open(88,file=filename,status='old',action='read')
      read(88,*) time, dto, istep
      close(88)
    case('w')
      if(myid.eq.0) then
        open (88, file=filename)
        write(88,'(2E15.7, 1I9.8)') time, dto, istep
        close(88)
      endif
    end select
    !
  end subroutine
  !
  subroutine load_ns(io,filename,n,u,v,w,p,time,istep)
    !
    ! reads/writes a restart file
    !
    implicit none
    !
    character(len=1), intent(in   )                            :: io
    character(len=*), intent(in   )                            :: filename
    integer         , intent(in   ), dimension(3)              :: n
    real(rp)        , intent(inout), dimension(n(1),n(2),n(3)) :: u,v,w,p
    real(rp)        , intent(inout)                            :: time,istep
    !
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
    real(rp)  , dimension(2) :: fldinfo
    integer(8), dimension(3) :: ng
    integer(8) :: lenr
    integer    :: fh
    !
    select case(io)
    case('r')
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      !ng(:)   = n(:)
      !ng(1:3) = ng(1:3)*dims(1:3)
      lenr = storage_size(time)/8
      good = int(product(ng)*4+2,MPI_OFFSET_KIND)*lenr
      if(filesize.ne.good) then
        if(myid.eq.0) print*, ''
        if(myid.eq.0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
        if(myid.eq.0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call decomp_2d_finalize
        call MPI_FINALIZE(ierr)
        call exit
      endif
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_read_var(fh,disp,ipencil,u   )
      call decomp_2d_read_var(fh,disp,ipencil,v   )
      call decomp_2d_read_var(fh,disp,ipencil,w   )
      call decomp_2d_read_var(fh,disp,ipencil,p   )
      call decomp_2d_read_scalar(fh,disp,2,fldinfo)
      time  = fldinfo(1)
      istep = fldinfo(2)
      call MPI_FILE_CLOSE(fh,ierr)
      !
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_write_var(fh,disp,ipencil,u  )
      call decomp_2d_write_var(fh,disp,ipencil,v  )
      call decomp_2d_write_var(fh,disp,ipencil,w  )
      call decomp_2d_write_var(fh,disp,ipencil,p  )
      fldinfo = (/time,istep/)
      call decomp_2d_write_scalar(fh,disp,2,fldinfo)
      call MPI_FILE_CLOSE(fh,ierr)
      !
    end select
    !
  end subroutine load_ns
  !
  subroutine loadIBM(in,n,cell_phi_tag,cell_u_tag,cell_v_tag,cell_w_tag,level_set,nx_surf, &
                   ny_surf,nz_surf,nabs_surf,deltan,i_mirror,j_mirror,k_mirror, &
                   i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                   WP1,WP2)
  implicit none
  integer,intent(in):: in
  integer, intent(in), dimension(3) :: n
  real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(inout) :: cell_phi_tag,cell_u_tag,cell_v_tag,cell_w_tag
  integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(inout) :: level_set,i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
  integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(inout) :: i_mirror,j_mirror,k_mirror
  real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(inout) :: nx_surf,ny_surf,nz_surf,nabs_surf,deltan
  real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7),intent(inout) :: WP1,WP2
  integer :: fh
  integer(kind=MPI_OFFSET_KIND) :: filesize,disp
  real(rp), dimension(n(1),n(2),n(3)) :: temp
  integer:: imax,jmax,kmax,i1,j1,k1
  integer error,request,status(MPI_STATUS_SIZE)
  
  imax=n(1)
  jmax=n(2)
  kmax=n(3)
  i1=imax+1
  j1=jmax+1
  k1=kmax+1
  
  if (in.eq.0) then
    call MPI_FILE_OPEN(MPI_COMM_WORLD, './IBMRestart', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    call MPI_FILE_CLOSE(fh,error)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, './IBMRestart', &
         MPI_MODE_RDONLY, MPI_INFO_NULL,fh, error)
    disp = 0_MPI_OFFSET_KIND
  !1---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    cell_phi_tag(1:imax,1:jmax,1:kmax) = temp
  !2---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    cell_u_tag(1:imax,1:jmax,1:kmax) = temp
  !3---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    cell_v_tag(1:imax,1:jmax,1:kmax) = temp
  !4---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    cell_w_tag(1:imax,1:jmax,1:kmax) = temp
  !5---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    level_set(1:imax,1:jmax,1:kmax) = int(temp)
  !6---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    nx_surf(1:imax,1:jmax,1:kmax) = temp
  !7---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    ny_surf(1:imax,1:jmax,1:kmax) = temp
  !8---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    nz_surf(1:imax,1:jmax,1:kmax) = temp
  !9---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    nabs_surf(1:imax,1:jmax,1:kmax) = temp
  !10---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    deltan(1:imax,1:jmax,1:kmax) = temp
  !11---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    i_mirror(1:imax,1:jmax,1:kmax) = int(temp)
  !12---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    j_mirror(1:imax,1:jmax,1:kmax) = int(temp)
  !13---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    k_mirror(1:imax,1:jmax,1:kmax) = int(temp)
  !14---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    i_IP1(1:imax,1:jmax,1:kmax) = int(temp)
  !15---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    j_IP1(1:imax,1:jmax,1:kmax) = int(temp)
  !16---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    k_IP1(1:imax,1:jmax,1:kmax) = int(temp)
  !17---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    i_IP2(1:imax,1:jmax,1:kmax) = int(temp)
  !18---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    j_IP2(1:imax,1:jmax,1:kmax) = int(temp)
  !19---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    k_IP2(1:imax,1:jmax,1:kmax) = int(temp)
  !20---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP1(1:imax,1:jmax,1:kmax,1) = temp
  !21---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP1(1:imax,1:jmax,1:kmax,2) = temp
  !22---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP1(1:imax,1:jmax,1:kmax,3) = temp
  !23---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP1(1:imax,1:jmax,1:kmax,4) = temp
  !24---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP1(1:imax,1:jmax,1:kmax,5) = temp
  !25---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP1(1:imax,1:jmax,1:kmax,6) = temp
  !26---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP1(1:imax,1:jmax,1:kmax,7) = temp
  !27---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP2(1:imax,1:jmax,1:kmax,1) = temp
  !28---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP2(1:imax,1:jmax,1:kmax,2) = temp
  !29---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP2(1:imax,1:jmax,1:kmax,3) = temp
  !30---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP2(1:imax,1:jmax,1:kmax,4) = temp
  !31---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP2(1:imax,1:jmax,1:kmax,5) = temp
  !32---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP2(1:imax,1:jmax,1:kmax,6) = temp
  !33---------------------------------------------
    call decomp_2d_read_var(fh,disp,3,temp)
    WP2(1:imax,1:jmax,1:kmax,7) = temp
  
    call MPI_FILE_CLOSE(fh,error)
  endif
  !
  if (in.eq.1) then
    call MPI_FILE_OPEN(MPI_COMM_WORLD, './IBMRestart', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
  !1----------------------------------------------
    temp = cell_phi_tag(1:imax,1:jmax,1:kmax)
    call decomp_2d_write_var(fh,disp,3,temp)
  !2---------------------------------------------
    temp = cell_u_tag(1:imax,1:jmax,1:kmax)
    call decomp_2d_write_var(fh,disp,3,temp)
  !3---------------------------------------------
    temp = cell_v_tag(1:imax,1:jmax,1:kmax)
    call decomp_2d_write_var(fh,disp,3,temp)
  !4---------------------------------------------
    temp = cell_w_tag(1:imax,1:jmax,1:kmax)
    call decomp_2d_write_var(fh,disp,3,temp)
  !5---------------------------------------------
    temp = 1.0000*level_set(1:imax,1:jmax,1:kmax) 
    call decomp_2d_write_var(fh,disp,3,temp)
  !6---------------------------------------------
    temp = nx_surf(1:imax,1:jmax,1:kmax)
    call decomp_2d_write_var(fh,disp,3,temp)
  !7---------------------------------------------
    temp = ny_surf(1:imax,1:jmax,1:kmax)
    call decomp_2d_write_var(fh,disp,3,temp)
  !8---------------------------------------------
    temp = nz_surf(1:imax,1:jmax,1:kmax)
    call decomp_2d_write_var(fh,disp,3,temp)
  !9---------------------------------------------
    temp = nabs_surf(1:imax,1:jmax,1:kmax)
    call decomp_2d_write_var(fh,disp,3,temp)
  !10---------------------------------------------
    temp = deltan(1:imax,1:jmax,1:kmax)
    call decomp_2d_write_var(fh,disp,3,temp)
  !11---------------------------------------------
    temp = 1.0_rp*(i_mirror(1:imax,1:jmax,1:kmax))
    call decomp_2d_write_var(fh,disp,3,temp)
  !12---------------------------------------------
    temp = 1.0_rp*(j_mirror(1:imax,1:jmax,1:kmax))
    call decomp_2d_write_var(fh,disp,3,temp)
  !13---------------------------------------------
    temp = 1.0_rp*(k_mirror(1:imax,1:jmax,1:kmax))
    call decomp_2d_write_var(fh,disp,3,temp)
  !14---------------------------------------------
    temp = 1.0_rp*(i_IP1(1:imax,1:jmax,1:kmax))
    call decomp_2d_write_var(fh,disp,3,temp)
  !15---------------------------------------------
    temp = 1.0_rp*(j_IP1(1:imax,1:jmax,1:kmax))
    call decomp_2d_write_var(fh,disp,3,temp)
  !16---------------------------------------------
    temp = 1.0_rp*(k_IP1(1:imax,1:jmax,1:kmax))
    call decomp_2d_write_var(fh,disp,3,temp)
  !17---------------------------------------------
    temp = 1.0_rp*(i_IP2(1:imax,1:jmax,1:kmax))
    call decomp_2d_write_var(fh,disp,3,temp)
  !18---------------------------------------------
    temp = 1.0_rp*(j_IP2(1:imax,1:jmax,1:kmax))
    call decomp_2d_write_var(fh,disp,3,temp)
  !19---------------------------------------------
    temp = 1.0_rp*(k_IP2(1:imax,1:jmax,1:kmax))
    call decomp_2d_write_var(fh,disp,3,temp)
  !20---------------------------------------------
    temp = WP1(1:imax,1:jmax,1:kmax,1)
    call decomp_2d_write_var(fh,disp,3,temp)
  !21---------------------------------------------
    temp = WP1(1:imax,1:jmax,1:kmax,2)
    call decomp_2d_write_var(fh,disp,3,temp)
  !22--------------------------------------------
    temp = WP1(1:imax,1:jmax,1:kmax,3)
    call decomp_2d_write_var(fh,disp,3,temp)
  !23---------------------------------------------
    temp = WP1(1:imax,1:jmax,1:kmax,4)
    call decomp_2d_write_var(fh,disp,3,temp)
  !24---------------------------------------------
    temp = WP1(1:imax,1:jmax,1:kmax,5)
    call decomp_2d_write_var(fh,disp,3,temp)
  !25---------------------------------------------
    temp = WP1(1:imax,1:jmax,1:kmax,6)
    call decomp_2d_write_var(fh,disp,3,temp)
  !26---------------------------------------------
    temp = WP1(1:imax,1:jmax,1:kmax,7)
    call decomp_2d_write_var(fh,disp,3,temp)
  !27---------------------------------------------
    temp = WP2(1:imax,1:jmax,1:kmax,1)
    call decomp_2d_write_var(fh,disp,3,temp)
  !28---------------------------------------------
    temp = WP2(1:imax,1:jmax,1:kmax,2)
    call decomp_2d_write_var(fh,disp,3,temp)
  !29---------------------------------------------
    temp = WP2(1:imax,1:jmax,1:kmax,3)
    call decomp_2d_write_var(fh,disp,3,temp)
  !30---------------------------------------------
    temp = WP2(1:imax,1:jmax,1:kmax,4)
    call decomp_2d_write_var(fh,disp,3,temp)
  !31---------------------------------------------
    temp = WP2(1:imax,1:jmax,1:kmax,5)
    call decomp_2d_write_var(fh,disp,3,temp)
  !32---------------------------------------------
    temp = WP2(1:imax,1:jmax,1:kmax,6)
    call decomp_2d_write_var(fh,disp,3,temp)
  !33---------------------------------------------
    temp = WP2(1:imax,1:jmax,1:kmax,7)
    call decomp_2d_write_var(fh,disp,3,temp)
    call MPI_FILE_CLOSE(fh,error)
  endif
  !
end subroutine loadIBM
end module mod_load