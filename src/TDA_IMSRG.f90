program TDA_IMSRG
  use HFQDmod   !!! hartree fock solver
  use ME_general  !!! general coulomb stuff, full_ham type
  use IMSRG_tools   !!! generator stuff, white style derivs 
  use commutators   !!! general modular derivs
  use operators    !!! various observables (TDA,density) 
  use adams_ode   !!! diff-eq solver
  implicit none
  
  real(8),parameter ::  conv_criteria = 1e-7
  integer :: n,emax,m,g,ghol,gpart,i,j,k,l,a,b,rx(2,3),pr,q,nh,np,nb,ML,MS,ST
  integer :: p1,p2,h3,h4,h1,h2,p3,p4,neq,iwork(5),flag,info,js,m3,block_index
  real(8),allocatable,dimension(:,:) :: coefs
  real(8),allocatable,dimension(:) :: cur_vec,d_vec,work2,ex_en
  real(8) :: hw,eHF,e2nd,HF_E2,MBPT2,ex,ex2,time,time2,s,rel,abse,crit,trc
  real(8) :: omp_get_wtime,nrm1,nrm2,stp,ocrit,s_off,stp0,sm,sm2,sm3,sm4
  real(8) :: get_crit
  type(full_ham) :: HS,ETA , HD, w1,w2,H0
  type(full_sp_block_mat) :: TDA,OLD
  character(5) :: hwstr,nstr,emaxstr,offstr,mlstr,msstr,cutstr
  character(2) :: nhs,nps,nbs
  character(7) :: genstr
  logical :: test,check_conv
  external :: dG2

!=================================================================
 !!! gather inputs, start timer, open file for GS write

time=omp_get_wtime()
 
call getarg(1,nstr)
call getarg(2,hwstr)
call getarg(3,emaxstr)
call getarg(4,mlstr)
call getarg(5,msstr) 
call getarg(6,cutstr) 

read(nstr,'(I5)') n
read(hwstr,'(f5.2)') hw
read(emaxstr,'(I5)') emax
!read(offstr, '(f5.2)') s_off
read(mlstr,'(I5)') HS%mltarg
read(msstr,'(I5)') HS%mstarg
read(cutstr,'(I5)') HS%cutshell

m = emax*(emax+1) !basis

write(hwstr,'(f5.2)') hw

nstr = adjustl(nstr)
hwstr= adjustl(hwstr) 
emaxstr = adjustl(emaxstr)


!=================================================================

!=================================================================
!!! construct HF basis
  allocate(coefs(m,m))

  call construct_two_particle_HF_basis( n , hw , emax , HS, coefs )
  call import_hamiltonian(HS,'ham_'//trim(nstr)//'_'//trim(hwstr) &
       //'_'//trim(emaxstr))    

!!! allocate workspaces ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  HS%IMSRG3 = .false. ! use IMSRG(2) formulas
  neq=get_num_elements( HS )
  allocate(cur_vec(neq))
  allocate(d_vec(neq))
  allocate(work2(100+21*neq))
   
  HS%neq=neq
 
 ! call divide_work(HS) 
  call allocate_everything(HS,ETA) 
  call allocate_everything(HS,HD) 
  call allocate_everything(HS,w1) 
  call allocate_everything(HS,w2)
  
  ETA%herm = -1
 !call system('rm CI_spectrum.dat') 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !call run_simple_CI('n')  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
!!! set parameters for solver ~~~~~~~~~~~~~~~~~~~~~~~
  rel=1e-8       ! relative error
  abse=1e-6      ! absolute error
  flag=1         ! some stupid flag for the solver
  s=0.d0         ! inital flow parameter
  stp=0.1        ! initial step size
  crit = 1.d0    ! initial convergence test
  ocrit = 10.d0  ! previous critera
  s_off = 0.d0   ! offset of s 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 
!==================================================================  
 !  if you want to plot the evolution with "s"
  open(unit=31, &
  file='../output/spectrum_'//trim(hwstr)//'_'//&
 trim(nstr)//'_'//trim(emaxstr)//'.dat',position='append')
  
! open(unit=41, &
!  file='../output/convergence_'//trim(hwstr)//'_'//&
! trim(nstr)//'_'//trim(emaxstr)//'.dat')
!=================================================================
!!! start IM-SRG loop
  
  call build_TDAmat(HS,TDA)
 ! call make_map
  call calc_TDA(HS,TDA)
  call diagonalize_blocks(TDA) 
  call write_spec
  
!================================================================
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!================================================================
!!! DECOUPLE EXCITED STATES 
!================================================================
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!================================================================

!!!! set a few constants ~~~~~~
 call allocate_old(TDA,OLD)
 test = .true.
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 pr = 0
 do while ( test )
    
     sm =0.d0 
     call vectorize(HS,cur_vec,neq)
     call ode( dG2, neq , cur_vec, HS, s, s+stp, rel, abse, flag, work2, iwork ) 
     call repackage(HS,cur_vec,neq)
     call calc_TDA(HS,TDA)
     call diagonalize_blocks(TDA) 
     call write_spec
        
     test = check_conv(TDA,OLD,HS,block_index) 
     
      !pr = pr + 1 
      !if (pr == 3) then 
      !call run_simple_CI('n')
      !pr = 0 
      !end if 
  end do


call final_write_full
contains 
!======================================================
subroutine final_write_full 
  
 open(unit=35, &
  file='../output/final_spec_'//trim(hwstr)//'_'//&
 trim(nstr)//'_'//trim(emaxstr)//'.dat')
 
  write(35,'(f12.7,2(I3),f12.7)')  HS%E0 , 0,0 , 0.0
  do q = 1,TDA%blocks  
     do i = 1,TDA%map(q)   
         
         write(35,'(f12.7,2(I3),f12.7)') TDA%blkM(q)%Eigval(i) + HS%E0 , &
  TDA%blkM(q)%lmda , abs(TDA%blkM(q)%Eigval(i) - OLD%blkM(q)%matrix(i,1))
     
     end do 
 end do 
 close(35)
end subroutine  
!===================================================
subroutine run_simple_CI(evecs) 
  ! even less efficient
  implicit none 
  
  real(8) :: e_off
  character(10) :: sstr, offstr
  character(1) :: evecs
  integer :: ii,jj 
  open(unit=39,file = 'OBME.dat',form='unformatted') 
  
 
  do i = 1, m
     do j = i,m 
        
        sm = 0.d0 
        do ii = 1, n
           sm = sm + v_elem(i,ii,j,ii,HS) 
        end do 
        
        write(39) f_elem(i,j,HS)-sm 
    end do 
 end do 
 
 close(39)
 
 open(unit=39,file = 'TBME.dat',form='unformatted')     
      
 do i=1,m
    do j = i+1,m
       
       do k = i,m
          if (k > i) then
          do l = k+1,m
       
             write(39) v_elem(i,j,k,l,HS)
          end do 
          else
          do l = max(j,k+1),m
        
             write(39) v_elem(i,j,k,l,HS)
          end do 
          end if
       end do 
   end do 
 end do 
 close(39) 

  e_off = 0.d0 
  !!! calculate offset from un-normal ordering
  do i=1,n
     e_off = e_off + HS%fhh(i,i)
  end do 
 
  do q = 1,HS%nblock
     do i=1,HS%mat(q)%nhh
        e_off = e_off - HS%mat(q)%Vhhhh(i,i)
     end do 
  end do 
  
  e_off = HS%E0 - e_off
    
  write(sstr,'(f10.5)') s+s_off
  write(offstr, '(f10.5)') e_off
  sstr = adjustl(sstr)
  offstr = adjustl(offstr) 
  call system('./run_CI '//emaxstr//' '//nstr//' '//sstr//' '&
       //offstr//' '//evecs//' 0 0') 
  
end subroutine
!===============================================        
subroutine write_spec
  implicit none 
  
  character(3) :: levnum
  character(15) :: fmt
  real(8),dimension(n*(m-n)+2) :: levs
  real(8) :: off

  off = HS%E0 - off 
  levs(1) = s + s_off
  levs(2) = HS%E0
  i=3
  do q=1,TDA%blocks
     if ( TDA%map(q) > 0 ) then 
        levs(i:i+TDA%map(q)-1) = TDA%blkM(q)%eigval + HS%E0
        i = i + TDA%map(q)
     end if 
  end do 
  
  i = i-1
  write(levnum,'(I3)') i
  levnum=adjustl(levnum) 
  fmt =  '('//trim(levnum)//'(f12.7))' 
  
  write(31,trim(fmt)) levs(1:i)
 
end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
subroutine make_map
  implicit none 
  
  open(unit=30, &
  file='../output/map_'//trim(hwstr)//'_'//&
 trim(nstr)//'_'//trim(emaxstr)//'.dat')
 
  i=3
 
  do q=1,TDA%blocks
     
     write(30,'(3(I5))') TDA%blkM(q)%lmda,i
     i=i+TDA%map(q)
  end do 
 
  close(30)
end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end program
!================================================
!================================================
subroutine DG2(t,yp,H) 
  ! calculates the derivatives inside the solver
  use ME_general
  use IMSRG_tools
  use commutators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx
  real(8) :: yp(*)
  type(full_ham) :: H,HD,ETA,DER,w1,w2
   
!!! we need the full_ham structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 
  
  call allocate_everything( H , ETA ) !! holds the generator
  call allocate_everything( H , DER ) !! holds derivatives
  call allocate_everything( H , w1 ) !! workspace for computing
  call allocate_everything( H , w2 ) !! workspace for computing 
  call allocate_everything( H , HD ) !! diagonal part of hamiltonian 

  ETA%herm = -1

! construct the generator...
  call build_imaginary_time(ETA,H)

  m=H%msp
  n=H%nbody
  neq=H%neq
  
!! zero body derivatives
  call xcommutator_110(m-n,n,ETA,H,ex)
  call xcommutator_220(ETA,H,ex2)
  
  DER%E0 = ex+ex2
  
!! one body derivatives
  call xcommutator_111(ETA,H,DER) 
  call xcommutator_121(ETA,H,DER) 
  call xcommutator_221(ETA,H,DER,w1,w2) 

!! two body derivatives
  call xcommutator_122(ETA,H,DER)
  call xcommutator_222(ETA,H,DER,w1) 

!! re-write in a way that shampine and gordon can handle  

  call vectorize(DER,yp,neq)

end subroutine 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
logical function check_conv(TDA,OLD,H,hold) 
  !!! check if some of the spectra is converging
  !!! not all of it will or can converge
  use ME_general
  implicit none 
  
  integer :: q,d,i,passing,total,hold
  type(full_sp_block_mat) :: TDA,OLD
  type(full_ham) :: H
  real(8) :: difference, frac
  logical :: failure
  
  frac = .98

  passing = 0
  total = 0
  failure = .false.
  q = hold
  do q = 1, TDA%blocks
     d = TDA%map(q)
     
     OLD%blkM(q)%extra = -10.
    

!=== uncomment for entire block analysis ============= 
     do i = 1, d
        total = total + 1
        difference = &
        abs(TDA%blkM(q)%Eigval(i) - OLD%blkM(q)%Eigval(i))
        if (TDA%blkM(q)%Eigval(i) < -4*H%E0) failure = .true.
        if (difference < 1e-6) then
           passing = passing + 1
           OLD%blkM(q)%extra(i) = 10. 
        else if (difference < 1e-5) then 
           OLD%blkM(q)%extra(i) = 10. 
        end if 
     end do
!=======================================================
     
     OLD%blkM(q)%matrix(:,1) = OLD%blkM(q)%Eigval
     OLD%blkM(q)%Eigval = TDA%blkM(q)%Eigval
     
  end do 
  print*, passing, total 
  if (float(passing)/float(total) > frac )then 
     check_conv = .false.   
  else
     check_conv = .true. 
  end if 
  
  if (failure) then 
     print*, 'CONVERGENCE FAILED'
     check_conv = .false. 
  end if 
end function
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine allocate_old(TDA,OLD) 
  use ME_general
  implicit none
  
  integer :: q,d
  type(full_sp_block_mat) :: TDA, OLD
  
  allocate(OLD%blkM(TDA%blocks))
  do q = 1,TDA%blocks
     allocate(OLD%blkM(q)%Eigval(TDA%map(q)))
     allocate(OLD%blkM(q)%extra(TDA%map(q)))
     allocate(OLD%blkM(q)%matrix(TDA%map(q),1))
     OLD%blkM(q)%Eigval = 0.d0
  end do 
  
end subroutine     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

  
  
  
  
