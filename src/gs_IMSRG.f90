program gs_IMSRG
!=========================================
! BASIC GROUND STATE DECOUPLING 
! ========================================  
  use HFQDmod   !!! hartree fock solver
  use ME_general  !!! general coulomb stuff, full_ham type
  use IMSRG_tools   !!! generator stuff
  use commutators   !!! general modular derivs
  use operators    !!! various observables (TDA,density) 
  use adams_ode   !!! diff-eq solver
  use EOM_IMSRG !!! lanczos solver
  implicit none

  real(8),parameter ::  conv_criteria = 1e-7
  integer :: n,emax,m,g,ghol,gpart,i,j,k,l,a,b,rx(2,3),pr,q,nh,np,nb,numstates
  integer :: p1,p2,h3,h4,h1,h2,p3,p4,neq,iwork(5),flag,info,js,m3,block_index
  real(8),allocatable,dimension(:,:) :: coefs
  real(8),allocatable,dimension(:) :: cur_vec,d_vec,work2,ex_en
  real(8) :: hw,eHF,e2nd,HF_E2,MBPT2,ex,ex2,time,time2,s,rel,abse,crit,trc
  real(8) :: omp_get_wtime,nrm1,nrm2,stp,ocrit,s_off,stp0,sm,sm2,sm3,sm4
  real(8) :: get_crit,pert
  type(full_ham) :: HS,ETA , HD, w1,w2,H0
  type(full_sp_block_mat) :: TDA,OLD
  type(cc_mat) :: HCC,ETACC
  character(5) :: hwstr,nstr,emaxstr,mlstr,msstr,cutstr
  character(2) :: nhs,nps,nbs
  character(7) :: genstr
  character(9) :: pstr
  logical :: test,check_conv,run_tda
  external :: dGam
!=================================================================
 !!! gather inputs, start timer, open file for GS write

  time=omp_get_wtime()

  call getarg(1,nstr)
  call getarg(2,hwstr)
  call getarg(3,emaxstr)
  call getarg(4,mlstr)
  call getarg(5,msstr) 

  read(mlstr,'(I5)') HS%mltarg
  read(msstr,'(I5)') HS%mstarg
  
  read(nstr,'(I5)') n
  read(hwstr,'(f5.2)') hw
  read(emaxstr,'(I5)') emax
  
  m = emax*(emax+1) !basis

  write(hwstr,'(f5.2)') hw
  
  HS%hospace = hw
  
  nstr = adjustl(nstr)
  hwstr= adjustl(hwstr) 
  emaxstr = adjustl(emaxstr)
  !=================================================================


!!! construct HF basis
  allocate(coefs(m,m))
  call construct_two_particle_HF_basis( n , hw , emax , HS, coefs,HCC )
  eHF=HF_E2( HS )
  e2nd = eHF +MBPT2( HS )   
  print*, 'Hartree-Fock Energy: ', eHF
  
  ! these are cross coupled arrays 
  call copy_cc(HCC,ETACC) 
  
!!! DECOUPLE GROUND STATE
  

!!! allocate workspaces ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  neq=get_num_elements( HS )
  allocate(cur_vec(neq))
  allocate(d_vec(neq))
  allocate(work2(100+21*neq)) ! terrible work array 
  
  HS%neq=neq

  call allocate_everything(HS,ETA)  ! making copies of HS for eta and intermediates 
  call allocate_everything(HS,HD) 
  call allocate_everything(HS,w1) 
  call allocate_everything(HS,w2)
  
  ETA%herm = -1 

!!! set parameters for solver ~~~~~~~~~~~~~~~~~~~~~~~
  rel=1e-8       ! relative error
  abse=1e-6      ! absolute error
  flag=1         ! some stupid flag for the solver
  s=0.d0         ! inital flow parameter
  s_off=0.d0     ! offset in s
  stp=0.1        ! initial step size
  crit = 1.d0    ! initial convergence test
  ocrit = 10.d0  ! previous critera
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pr = 1

  do while (crit > conv_criteria) 
 
     sm =0.d0 
     crit=HS%E0 
     call vectorize(HS,cur_vec,neq)
     call ode( dGam, neq , cur_vec, HS,HCC,ETACC, s, s+stp, rel, abse, flag, work2, iwork ) 
     call repackage(HS,cur_vec,neq)     
     
     if (run_tda) then 
      call calc_TDA(HS,TDA)
      call diagonalize_blocks(TDA) 
      call write_spec
     end if 
     
     crit = abs( (crit - HS%E0 )/crit ) 
     write(*,'(I5,3(e14.6))') pr,s,HS%E0,crit
     pr = pr + 1 

  end do
    
  numstates = 5

  if (HS%Mstarg == 0 ) then 
     call calculate_excited_states( HS%Mltarg, HS%Mstarg, numstates , HS ) 
  else
     call calculate_1p_attached( HS%Mltarg, HS%Mstarg, numstates , HS ) 
     call calculate_1h_removed(  HS%Mltarg, HS%Mstarg, numstates , HS ) 
  end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!=================================================================
!================================================================
!!! write all the data to file
time2 = omp_get_wtime()

 open(unit=39,file='../output/TRAD_'//trim(nstr)//'_'//trim(hwstr) &
//'.dat',position='append')
  write(39,'(I4,I5,3(f15.9),f15.4)') emax,m,eHF,e2nd,HS%E0,log(time2-time)
  close(39)
  
  call export_hamiltonian(HS,'ham_'//trim(nstr)//'_'//trim(hwstr) &
       //'_'//trim(emaxstr)) 

contains 
!=================================================================
!=================================================================
subroutine run_simple_CI(evecs) 
  ! this writes the current Hamiltonian to file and
  ! runs CI on it. 
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
 
 open(unit=39,file = 'TBME.dat',form = 'unformatted')     
      
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
 ! print*, 'running CI'
  call system('./run_CI '//emaxstr//' '//nstr//' '//sstr//' '&
       //offstr//' '//evecs//' 0 0') 
  
end subroutine
!===============================================        
end program
!================================================
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!================================================
subroutine dGam(t,yp,H,HCC,ETACC) 
  ! calculates the derivatives inside the solver
  use ME_general
  use IMSRG_tools
  use commutators
  implicit none 
  
  real(8) :: t,ex,ex2
  integer :: i,j,neq,bytes,n,m,p,q,r,s,px,qx,rx,sx
  real(8) :: yp(*)
  type(full_ham) :: H,HD,ETA,DER,w1,w2
   type(cc_mat) :: HCC,ETACC
!!! we need the full_ham structure to compute the derivatives at max speed
!!! so we allocate a bunch of those to work in 

  call allocate_everything( H , ETA ) !! holds the generator
  call allocate_everything( H , DER ) !! holds derivatives
  call allocate_everything( H , w1 ) !! workspace for computing
  call allocate_everything( H , w2 )  !! transpose space for non-symmetric

  call calc_cc(H,HCC)
  ETA%herm = -1
  HCC%herm = 1.d0
! construct the generator...
  call build_white(ETA,H)
  call calc_cc(ETA,ETACC)

  ETACC%herm = -1.d0
  m=H%msp
  n=H%nbody
  neq=H%neq
  
!! zero body derivatives
  call commutator_110(m-n,n,ETA,H,ex)
  call commutator_220(ETA,H,ex2)
  
  DER%E0 = ex+ex2

!! one body derivatives
  call commutator_111(ETA,H,DER) 
  call commutator_121(ETA,H,DER)  
  call commutator_221(ETA,H,DER,w1,w2) 
  
!! two body derivatives
  call commutator_122(ETA,H,DER)
  call commutator_222(ETA,H,DER,w1,HCC,ETACC) 

!! re-write in a way that shampine and gordon can handle
  call vectorize(DER,yp,neq)
   
end subroutine 
