program magnus_IMSRG
  use HFQDmod   !!! hartree fock solver
  use ME_general  !!! general coulomb stuff, full_ham type
  use IMSRG_tools   !!! generator stuff, white style derivs 
  use commutators   !!! general modular derivs
  use operators
  implicit none

  real(8),parameter ::  conv_criteria = 1e-5
  integer :: n,emax,m,i,j,k,l,a,b,rx(2,3),pr,q,nh,np,nb
  integer :: p1,p2,h3,h4,h1,h2,p3,p4,neq,iwork(5),flag,info,js
  real(8),allocatable,dimension(:,:) :: coefs
  real(8),allocatable,dimension(:) :: cur_vec,d_vec,work2
  real(8) :: hw,eHF,e2nd,HF_E2,MBPT2,ex,ex2,time,time2
  real(8) :: s,odx,ody,crit,hk,stp,stp0,nrm1,nrm2,s_off,ocrit,xcrit
  real(8) :: omp_get_wtime
  type(full_ham) :: HAM,ETA,G,DG,HD,HS,AD,INT1,INT2,w1,w2,ETA0,G0,HS0,DEN,DS
  type(full_sp_block_mat) :: TDA
  character(5) :: hwstr,nstr,emaxstr,genlab
  character(7) :: genstr
  character(1) :: status
  !!! type from ME_general

!=================================================================
 !!! gather inputs, start timer, open file for GS write

time=omp_get_wtime()
 
call getarg(1,nstr)
call getarg(2,hwstr)
call getarg(3,emaxstr)

genlab = 'white'

read(nstr,'(I5)') n
read(hwstr,'(f5.2)') hw
read(emaxstr,'(I5)') emax
m = emax*(emax+1) !basis

write(hwstr,'(f5.2)') hw

nstr = adjustl(nstr)
hwstr= adjustl(hwstr) 
emaxstr = adjustl(emaxstr)

!=================================================================

!=================================================================
!!! construct HF basis
  allocate(coefs(m,m))   !!! useless

  !! HFQDmod
  call construct_two_particle_HF_basis( n , hw , emax , HAM, coefs )
 
 
  eHF=HF_E2( HAM ) ! hartree fock energy
  e2nd = eHF +MBPT2( HAM )  ! second order perturbation theory
 
!=================================================================  

!!! set parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  neq=get_num_elements( HAM )
  HAM%neq=neq
  s=0.d0    ! inital flow parameter
  crit = 1.d0  !initial convergence test
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!IMSRG_tools !~~~~~ allocation of fancy arrays ~~~~~~~~~~~~
call allocate_everything( HAM , ETA ) ! generator
call allocate_everything( HAM , HS )  ! evolved hamiltonian
call allocate_everything( HAM , HD )  ! diagonal hamiltonian 
call allocate_everything( HAM , DEN )  ! HF density
call allocate_everything( HAM , DS )  ! evolved density 
call allocate_everything( HAM , HS0 )  ! diagonal hamiltonian 
call allocate_everything( HAM , ETA0 )  ! diagonal hamiltonian 
call allocate_everything( HAM , G0 )  ! diagonal hamiltonian 
call allocate_everything( HAM , G  ) !omega param
call allocate_everything( HAM , DG   ) !omega deriv
call allocate_everything( HAM , INT1  ) ! space
call allocate_everything( HAM , INT2  ) ! space
call allocate_everything( HAM , AD  ) ! space
call allocate_everything( HAM , w1  ) ! space
call allocate_everything( HAM , w2  ) ! space
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!!~~~~~~set arrays as ANTI-HERMITIAN~~~~
   ETA%herm = -1 
   ETA0%herm = -1
   G%herm = -1 
   G0%herm = -1
   DG%herm = -1  
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!==================================================================  

 !  if you want to plot the evolution with "s"
  !open(unit=42, &
  !file='../output/E0flow_'//trim(hwstr)//'_'//&
  !trim(nstr)//'_'//trim(emaxstr)//'_'//genlab//'.dat')

!=================================================================

  !write(42,*) s,eHF

!!! start IM-SRG loop for white~~~~~~~~~~~~~~~~
  hk=10.
  stp=.5
  s_off = 0.d0 
  stp0=2.d0

 call copy_arrays(HAM,HS)


open(unit=31,file='spectrum_magnus.dat')

  call build_TDAmat(HAM,TDA)  
  call make_map
  call calc_TDA(HS,TDA)
  call diagonalize_blocks(TDA) 
  call write_spec  

 call s_evolve(build_wegner) 

goto 12
 stp = .001 
 s_off = s
 stp0 = .5d0
 s = 0.d0
 crit = 1.d0
 
 call copy_arrays(HS,HAM)
 call set_to_zero(G) 
 
 call s_evolve(build_white_phhh) 
 
 stp = .001 
 s_off = s+s_off
 stp0 = .5d0
 s = 0.d0
 crit = 1.d0
 
 call copy_arrays(HS,HAM)
 call set_to_zero(G) 
 
 call s_evolve(build_white_ppph) 
12 print*, 'fuck!'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!if (emax == 7) then 
!  call write_density(HAM,den,hw,coefs,G,w1,w2,INT1,INT2,AD,DS,s,emax)
!  end if

!================================================================
!!! write all the data to file
  status = 'A'
  if (stp < 1e-2) status='X'
  time2=omp_get_wtime()
  !open(unit=39,file='../output/IMSRG_'//trim(nstr)//'_'//trim(hwstr) &
!//'_'//genlab//'.dat',position='append')
!  write(39,'(I4,I5,3(f15.9),f15.4,I2)') emax,m,eHF,e2nd,HS%E0,time2-time
write(*,'(I4,I5,3(f15.9),f15.4,I2)') emax,m,eHF,e2nd,HS%E0,time2-time
!  close(39)
!  close(42)
!================================================================
contains 
!===================================================
subroutine write_spec
  implicit none 
  
  character(3) :: levnum
  character(15) :: fmt
  real(8),dimension(n*(m-n)+2) :: levs

  levs(1) = s + s_off
  levs(2) = HS%E0
  i=3
  ocrit = 0.d0
  do q=1,TDA%blocks
     if ( TDA%map(q) > 0 ) then 
        levs(i:i+TDA%map(q)-1) = TDA%blkM(q)%eigval + HS%E0
        ocrit = ocrit + sum((HS%E0+TDA%blkM(q)%eigval)**2)
        i = i + TDA%map(q)
     end if 
  end do 
  
  i = n*(m-n)+2
  write(levnum,'(I3)') i
  levnum=adjustl(levnum) 
  fmt =  '('//trim(levnum)//'(f12.7))' 
  
  write(31,trim(fmt)) levs
 
end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
subroutine make_map
  implicit none 
  

  open(unit=30,file='map_spec.dat')
  i=3
 
  do q=1,TDA%blocks
     
     write(30,'(3(I5))') TDA%blkM(q)%lmda,i
     i=i+TDA%map(q)
  end do 
 
  close(30)
end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!================================================================
subroutine s_evolve(build_gen) 
  
  external :: build_gen, conv_test
  
  call build_gen(ETA,HS)     !(IMSRG_tools)    ! construct generator
  nrm1 = mat_2_norm(ETA)
  xcrit = 10.
  
  do while ( crit > conv_criteria ) 
   
     ! G is omega here. 
     call copy_arrays(G,G0)
     call MAGNUS_EXPAND(DG,G,ETA,INT1,INT2,AD,w1,w2) ! get d(omega)       
     call euler_step(G,DG,s,stp) ! evolve with a stupid euler step          

     !! baker-campbell-hausdorff expansion 
     call copy_arrays(HS,HS0)
     call BCH_EXPAND(HS,G,HAM,INT1,INT2,AD,w1,w2) 
      
     call copy_arrays(ETA,ETA0)
     call build_gen(ETA, HS ) 
     nrm2 = mat_2_norm(ETA)
  
     if ( nrm1 < nrm2 ) then 

         s=s-stp
         call copy_arrays(G0,G)
         call copy_arrays(HS0,HS)
         call copy_arrays(ETA0,ETA)
         stp = stp / 2.d0  + 1e-3
         stp0 = stp0 /1.2d0 +1e-2
         cycle 
     
      else 
         
       call calc_TDA(HS,TDA)
       call diagonalize_blocks(TDA) 
       call write_spec  
        !stp = .1*stp0 + .9*stp
        crit = abs( (xcrit - ocrit) / xcrit  ) ! wimpy convergence test
        nrm1=nrm2
        xcrit = ocrit
       
      end if   
      
      call calc_TDA(HS,TDA)
      call diagonalize_blocks(TDA) 
      call write_spec  
      
      write(42,*) s, crit
     !hk = HS%E0
     
!!~~~~~~for plotting the evolution~~~~~~~~~~~~~~~~~~~~~~~
    !write(42,*) s,hk

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
  end do 

end subroutine 
!================================================================
end program
!================================================================
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!================================================================
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!================================================================
subroutine write_density(r,den,hw,coefs,G,w1,w2,INT1,INT2,AD,DS,s,emax)
  use ME_general
  implicit none 
  
  integer,parameter :: offset = 0
  type(full_ham) :: r,den,G,w1,w2,INT1,INT2,AD,DS
  integer :: n,m,emax,i
  real(8) :: x, hw,dx,sdx,s,evolved_density
  real(8),dimension(50,2) :: lx
  real(8),dimension(r%msp,r%msp) :: coefs
  character(5) :: hwstr,nstr,emaxstr,sstr
  
  n = r%nbody
  m = r%msp
  write(hwstr, '(f5.2)') hw
  write(sstr, '(f5.2)') s
  write(nstr, '(i5)' ) n
  write(emaxstr, '(i5)' ) emax

  hwstr = adjustl(hwstr)
  sstr = adjustl(sstr)
  nstr = adjustl(nstr)
  emaxstr = adjustl(emaxstr) 

  open(unit=38,file='../output/density_'//trim(hwstr)//'_'//trim(nstr)// &
       '_'//trim(emaxstr)//'_'//trim(sstr)//'.dat')
  
 
!$omp parallel do private(x,dx,i), shared(coefs,hw,lx)

  do i=offset+1,offset+50
     
     x= 0.1d0*i-0.1d0
     dx = evolved_density(r,G,x,hw,coefs) 

     lx(i-offset,1) = x
     lx(i-offset,2) = dx

  end do 

!$omp end parallel do

  do i=1,50
     write(38,*) lx(i,1),lx(i,2)
  end do 
  close(38)

end subroutine
!===============================================================
!===============================================================
real(8) function evolved_density(r,G,x,hw,coefs)
  use ME_general
  use operators
  use IMSRG_tools
  implicit none 
  
  type(full_ham) :: r,den,DS,G,INT1,INT2,AD,w1,w2
  real(8) :: x,hw, coefs(r%msp,r%msp)
 
   call allocate_everything( r , DS   ) !omega deriv
   call allocate_everything( r , INT1  ) ! space
   call allocate_everything( r , INT2  ) ! space
   call allocate_everything( r , AD  ) ! space
   call allocate_everything( r , w1  ) ! space
   call allocate_everything( r , w2  ) ! space
   call allocate_everything( r , DEN  ) ! space

  call density(r,den,x,hw,coefs)    
  call BCH_NONPAR(DS,G,DEN,INT1,INT2,AD,w1,w2) 
  
  evolved_density=DS%E0
end function 
!===============================================================
!===============================================================
subroutine BCH_EXPAND(HS,G,H,INT1,INT2,AD,w1,w2) 
  use commutators
  use ME_general
  use IMSRG_tools 
  implicit none 
  
  real(8), parameter :: conv = 1e-5 
  integer :: trunc,i,m,n
  type(full_ham) :: H , G, ETA, INT1, INT2, HS, AD,w1,w2
  real(8) ::  cof(11),ex,ex2
 
  cof = (/1.d0,0.5d0,0.166666666666666666d0, &
       0.04166666666666666d0,0.0083333333333333333d0,.001388888888888d0, &
       1.984126984d-4,2.48015873d-5,2.755731922d-6,2.755731922d-7,2.505210839d-8/) 

  n = H%nbody
  m = H%msp
  
  ! intermediates must be HERMITIAN
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  

  !! so here:  H is the current hamiltonian and we 
  !! copy it onto HS.  We copy this onto 
  !! INT2, just to make things easier for everyone.

  call copy_arrays( H , HS )  !ME_general
  call copy_arrays( HS , INT2 )
 
  do i = 2 , 12
          
     ! current value of HS is renamed INT1 
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_arrays( HS , INT1) 
     call copy_arrays( INT2 , AD ) 
     ! so to start, AD is equal to H
  
     !now: INT2 = [ G , AD ]  

! zero body commutator
     call set_to_zero(INT2)
     call xcommutator_110(m-n,n,G,AD,ex)
     call xcommutator_220(G,AD,ex2)
  
     INT2%E0 = ex+ex2

!! one body 
     call xcommutator_111(G,AD,INT2) 
     call xcommutator_121(G,AD,INT2) 
     call xcommutator_221(G,AD,INT2,w1,w2) 

     call set_to_zero(w1)
!! two body
     call xcommutator_122(G,AD,INT2)
     call xcommutator_222(G,AD,INT2,w1) 
   
     
     ! so now just add INT1 + c_n * INT2 to get current value of HS
     call add_arrays(INT1 , 1.d0 , INT2 , cof(i-1) , HS )   !ME_general
     if ( mat_2_norm(INT2)*cof(i-1)/mat_2_norm(INT1) < conv ) exit 
  end do 
end subroutine 
!===============================================================
!===============================================================
subroutine BCH_NONPAR(HS,G,H,INT1,INT2,AD,w1,w2) 
  use commutators
  use ME_general
  use IMSRG_tools 
  implicit none 
  
  real(8), parameter :: conv = 1e-5 
  integer :: trunc,i,m,n
  type(full_ham) :: H , G, ETA, INT1, INT2, HS, AD,w1,w2
  real(8) ::  cof(11),ex,ex2
 
  cof = (/1.d0,0.5d0,0.166666666666666666d0,&
       0.04166666666666666d0,0.0083333333333333333d0,.001388888888888d0, &
       1.984126984d-4,2.48015873d-5,2.755731922d-6,2.755731922d-7,2.505210839d-8/) 

  n = H%nbody
  m = H%msp
  
  ! intermediates must be HERMITIAN
  INT2%herm = 1
  INT1%herm = 1 
  AD%herm = 1
  

  !! so here:  H is the current hamiltonian and we 
  !! copy it onto HS.  We copy this onto 
  !! INT2, just to make things easier for everyone.

  call copy_arrays( H , HS )  !ME_general
  call copy_arrays( HS , INT2 )
 
  do i = 2 , 12
          
     ! current value of HS is renamed INT1 
     ! INT2 is renamed AD, for the AD parameters in BCH and magnus expansions
     ! AD refers AD_n and INT2 becomes AD_{n+1} . 
     call copy_arrays( HS , INT1) 
     call copy_arrays( INT2 , AD ) 
     ! so to start, AD is equal to H
  
     !now: INT2 = [ G , AD ]  

! zero body commutator
     call set_to_zero(INT2)
     call xcommutator_110(m-n,n,G,AD,ex)
     call xcommutator_220(G,AD,ex2)
  
     INT2%E0 = ex+ex2

!! one body 
     call xcommutator_111(G,AD,INT2) 
     call xcommutator_121(G,AD,INT2) 
     call xcommutator_221(G,AD,INT2,w1,w2) 

     call set_to_zero(w1)
!! two body
     call xcommutator_122(G,AD,INT2)
     call xcommutator_222_nonpar(G,AD,INT2,w1) 
   
     
     ! so now just add INT1 + c_n * INT2 to get current value of HS
     call add_arrays(INT1 , 1.d0 , INT2 , cof(i-1) , HS )   !ME_general
!   
     if ( mat_2_norm(INT2)*cof(i-1)/mat_2_norm(INT1) < conv ) exit 
  end do 
end subroutine 
!==================================================================
!==================================================================
subroutine MAGNUS_EXPAND(DG,G,ETA,INT1,INT2,AD,w1,w2)
  use commutators
  use ME_general
  use IMSRG_tools 
  implicit none 
  
  real(8), parameter :: conv = 1e-5
  integer :: trunc,i
  type(full_ham) :: G, ETA, INT1, INT2, DG, AD,w1,w2
  real(8) ::  cof(4)
 
  ! Intermediates are ANTI-HERMITIAN 
  INT2%herm = -1
  INT1%herm = -1 
  AD%herm = -1
  
  cof = (/-0.5d0,.0833333333d0,0.d0,-0.00138888888d0/) 
  
  
  !! same deal as BCH expansion, which is explained ad nauseam above. 
  call copy_arrays( ETA, DG )  !ME_general
  call copy_arrays( DG , INT2 )
  
  do i = 2 , 5
          
     call copy_arrays( DG , INT1) 
     call copy_arrays( INT2 , AD ) 
   
     call set_to_zero(INT2)
!! one body commutators
     call xcommutator_111(G,AD,INT2) 
     call xcommutator_121(G,AD,INT2) 
     call xcommutator_221(G,AD,INT2,w1,w2) 

     call set_to_zero(w1)
!! two body 
     call xcommutator_122(G,AD,INT2)
     call xcommutator_222(G,AD,INT2,w1) 
  
   
     call add_arrays(INT1 , 1.d0 , INT2 , cof(i-1) , DG ) !ME_general
     if ( mat_2_norm(INT2)*cof(i-1) / mat_2_norm(INT1) < conv ) exit 
  end do 
  
  
end subroutine 
!=====================================================
!=====================================================
subroutine euler_step(G,DG,s,stp)
  use ME_general
  use IMSRG_tools
  implicit none 

  integer :: i
  real(8) :: s , stp
  type(full_ham) :: G , DG

  call add_arrays(G,1.d0,DG,stp,G) !durr probably wrong. 
  s = s + stp 

end subroutine 
!================================================
!================================================
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!================================================

