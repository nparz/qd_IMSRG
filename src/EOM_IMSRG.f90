module EOM_IMSRG
  use ME_general
  use commutators
  use IMSRG_tools
  implicit none
  

contains 

subroutine calculate_excited_states( Ml, Ms, Numstates, HS ) 
  implicit none
  
  type(full_ham) :: HS 
  type(full_ham),allocatable,dimension(:) :: ladder_ops 
  integer :: Ml,Ms,Numstates,i,q
  
  allocate(ladder_ops(numstates)) 
  
  ladder_ops%herm = 1
 
  ladder_ops%Mltarg = Ml 
  ladder_ops%Mstarg = Ms
  
  call allocate_ex(HS,ladder_ops(1)) 
 
  do i = 2, Numstates
     call allocate_ex(Hs,ladder_ops(i))
  end do 
  
  print* 
  write(*,'((A55),(I2),(A4),(I2))') 'EXECUTING EOM CALCULATION'// &
       ' FOR EXCITED STATES: Ml=',Ml,' Ms=',Ms  
  print*
  call lanczos_diagonalize(HS,ladder_ops,Numstates)  
  
  print*
  write(*,'((A21),(f16.9))') 'Ground State Energy: ',HS%E0 
  print*
  print*, 'EXCITED STATE ENERGIES:'
  print*, '=================================='
  print*, '      dE             E_0 + dE'
  print*, '=================================='
  do i = 1, Numstates
     write(*,'(3(f16.9))') ladder_ops(i)%E0 ,ladder_ops(i)%E0+HS%E0,sum(ladder_ops(i)%fph**2)
  end do
     
end subroutine 


subroutine calculate_1p_Attached( Ml, Ms, Numstates, HS ) 
  implicit none
  
  type(full_ham) :: HS 
  real(8),dimension(Numstates) :: eigs,norms
  integer :: Ml,Ms,Numstates,i,q,holes,parts,dm
  integer :: nshell,nfill
  
  
  holes = HS%nbody
  parts = HS%msp-holes
  dm = parts+holes*parts*(parts-1)/2 
  
  nshell = (sqrt(4.d0*HS%msp+1.d0)-1)/2.d0
  nfill =  (sqrt(4.d0*HS%nbody+1.d0)-1)/2.d0
!  allocate(ladder_ops(numstates,dm)) 
   
 ! ladder_ops = 0.d0 
  
  print* 
  write(*,'((A65),(I2),(A4),(I2))') 'EXECUTING EOM CALCULATION'// &
       ' FOR PARTICLE ATTACHED STATES: Ml=',Ml,' Ms=',Ms  
  print*
  call lanczos_1p_Attached(HS,Numstates,dm,eigs,norms)  
  
  print*
  write(*,'((A21),(f16.9))') 'Ground State Energy: ',HS%E0 
  print*
  print*, '1p ATTACHED ENERGIES:'
  print*, '=================================='
  print*, '      dE             E_0 + dE'
  print*, '=================================='
  do i = 1, Numstates
     write(*,'(3(f16.9))') eigs(i) ,eigs(i)+HS%E0,norms(i)
  end do

  open(unit=39,file='../output/EOM_particle_attached.dat',position='append')
  write(39,'(4(I5),f8.2,4(f17.9))'   )  nshell,nfill,ML,MS,HS%hospace,HS%E0,eigs(1),eigs(1)+HS%E0,norms(1)
  close(39) 
  
end subroutine calculate_1p_Attached
 



subroutine calculate_1h_Removed( Ml, Ms, Numstates, HS ) 
  implicit none
  
  type(full_ham) :: HS 
  real(8),dimension(Numstates) :: eigs,norms
  integer :: Ml,Ms,Numstates,i,q,holes,parts,dm,nshell,nfill
   
  holes = HS%nbody
  parts = HS%msp-holes
  dm = holes+parts*holes*(holes-1)/2 
  
  nshell = (sqrt(4.d0*HS%msp+1.d0)-1)/2.d0
  nfill =  (sqrt(4.d0*HS%nbody+1.d0)-1)/2.d0
! 
  print* 
  write(*,'((A64),(I2),(A4),(I2))') 'EXECUTING EOM CALCULATION'// &
       ' FOR PARTICLE REMOVED STATES: Ml=',Ml,' Ms=',Ms  
  print*
  call lanczos_1h_Removed(HS,Numstates,dm,eigs,norms)  
  
  print*
  write(*,'((A21),(f16.9))') 'Ground State Energy: ',HS%E0 
  print*
  print*, '1p REMOVED ENERGIES:'
  print*, '=================================='
  print*, '      dE             E_0 + dE'
  print*, '=================================='
  do i = 1, Numstates
     write(*,'(3(f16.9))') eigs(i) ,eigs(i)+HS%E0,norms(i)
  end do
  open(unit=39,file='../output/EOM_particle_removed.dat',position='append')
  write(39,'(4(I5),f8.2, 4(f17.9))'   )  nshell,nfill,ML,MS,HS%hospace,HS%E0,eigs(1),eigs(1)+HS%E0,norms(1)
  close(39) 
   
end subroutine calculate_1h_Removed


subroutine LANCZOS_DIAGONALIZE(OP,Vecs,nev)    
  
  integer :: N 
  integer,intent(in) :: nev
  type(full_ham) :: op,V1,Q1,Q2,w1,w2
  type(full_ham),dimension(nev) :: Vecs
  real(8),allocatable,dimension(:) :: workl,eigs,resid,work,workD
  real(8),dimension(nev) :: d,norm_1p
  real(8),allocatable,dimension(:,:) :: V,Z
  integer :: i,j,ix,jx,lwork,info,ido,ncv,ldv,iparam(11),ipntr(11),q,II,JJ
  integer :: ishift,mxiter,nb,nconv,mode,np,lworkl,ldz,p,h,sps,tps,jp,jh
  integer :: li,lj,lb,la,si,sj,sa,sb,Ml,Ms,h1,h2,p1,p2,a,b
  real(8) ::  x,tol,y,sigma,t1,t2
  character(1) :: BMAT,HOWMNY 
  character(2) :: which
  logical :: rvec
  logical,allocatable,dimension(:) :: selct

  call allocate_ex(Op,Q1) !workspace
  call allocate_ex(Op,Q2) !workspace
  
  h = OP%Nbody !holes
  p = OP%Msp-h  !particles
 
  sps = 0 
  do ix = 1,p
     do jx = 1,h
                          
         li = Op%states(Op%eh(jx),2)
         si = Op%states(Op%eh(jx),3)
         la = Op%states(Op%ep(ix),2)
         sa = Op%states(Op%ep(ix),3)
              
         Ml = la-li
         Ms = sa-si 
         
         if ( Ml .ne. Op%Mltarg ) cycle
         if ( Ms .ne. Op%Mstarg ) cycle   
         sps = sps+1
     end do
  end do
  
     ! scalar operator
  tps = 0 
 ! do q = 1, OP%nblock
     
!     do II = 1,OP%mat(q)%npp
 !       do JJ = 1, OP%mat(q)%nhh 

  do a = 1,p
     do b = a+1,p
        do i = 1, h
           do j = i+1,h
      
              p1=op%ep(a)
              p2=op%ep(b)
              h1=op%eh(i)
              h2=op%eh(j)
              
              
              li = Op%states(h1,2)
              si = Op%states(h1,3)
              la = Op%states(p1,2)
              sa = Op%states(p1,3)
              lj = Op%states(h2,2)
              sj = Op%states(h2,3)
              lb = Op%states(p2,2)
              sb = Op%states(p2,3)
           

              Ml = la+lb-li-lj
              Ms = sa+sb-si-sj 
              if ( Ml .ne. Op%Mltarg ) cycle
              if ( Ms .ne. Op%Mstarg ) cycle   
              tps = tps+ 1

           end do 
        end do 
     end do
  end do
  !      end do
  !   end do
 ! end do
           
  N = sps + tps ! number of ph and pphh SDs 
  
  ido = 0  ! status integer is 0 at start
  BMAT = 'I' ! standard eigenvalue problem (N for generalized) 
  which = 'SA' ! compute smallest eigenvalues in magnitude ('SA') is algebraic. 
  tol = 0.0E+00 ! error tolerance? (wtf zero?) 
  info = 0
  ncv = 2*nev ! number of lanczos vectors I guess
  lworkl = ncv*(ncv+8) 
  allocate(V(N,NCV),workl(lworkl))
  LDV = N  
  ishift = 1
  mxiter = 500 
  mode = 1

  allocate(eigs(N),resid(N),work(10*N),workD(3*N)) 
  
  iparam(1) = ishift
  iparam(3) = mxiter
  iparam(7) = mode
  i = 0

  do 
     ! so V is the krylov subspace matrix that is being diagonalized
     ! it does not need to be initialized, so long as you have the other 
     ! stuff declared right, the code should know this.

     call dsaupd ( ido, bmat, N, which, nev, tol, resid, &
      ncv, v, ldv, iparam, ipntr, workd, workl, &
      lworkl, info )

     ! The actual matrix only gets multiplied with the "guess" vector in "matvec_prod" 
!     call progress_bar( i )
     i=i+1 
  
    if ( ido /= -1 .and. ido /= 1 ) then
      exit
    end if

    call matvec_prod(N,OP,Q1,Q2,w1,w2,workd(ipntr(1)), workd(ipntr(2)) ) 
    
    end do 
    write(6,*) 
    
   
  ! the ritz values are out of order right now. Need to do post
  ! processing to fix this, and get the eigenvectors
  rvec= .true. 
  howmny = 'A'
  
  allocate(selct(NCV)) 
  allocate(Z(N,NEV)) 
  ldz = N  
  call dseupd( rvec, howmny, selct, d, Z, ldv, sigma, &
      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
      iparam, ipntr, workd, workl, lworkl, info )
  
  ! right now Z contains the eigenvectors in the columns
  ! d contains the eigenvalues in the same order. 
  do i = 1, nev
     call unwrap(Z(:,i),Vecs(i),N)
     Vecs(i)%E0 = d(i)
  end do
  
end subroutine LANCZOS_DIAGONALIZE
!==================================================================
!==================================================================
subroutine LANCZOS_1p_ATTACHED(OP,nev,dm,d,norm_1p)    
  
  integer :: N ,dm
  integer,intent(in) :: nev
  type(full_ham) :: op
  real(8),dimension(nev,dm) :: Vecs
  real(8),allocatable,dimension(:) :: workl,eigs,resid,work,workD
  real(8),allocatable,dimension(:,:) :: V,Z
  real(8),dimension(nev) :: d,norm_1p
  integer :: i,j,ix,jx,lwork,info,ido,ncv,ldv,iparam(11),ipntr(11),q,II,JJ
  integer :: ishift,mxiter,nb,nconv,mode,np,lworkl,ldz,p,h,sps,tps,jp,jh
  integer :: li,lj,lb,la,si,sj,sa,sb,Ml,Ms,h1,h2,p1,p2,a,b
  real(8) ::  x,tol,y,sigma,t1,t2
  character(1) :: BMAT,HOWMNY 
  character(2) :: which
  logical :: rvec
  logical,allocatable,dimension(:) :: selct

  h = OP%Nbody !holes
  p = OP%Msp-h  !particles
  
  sps = 0 
  do ix = 1,p

     la = Op%states(Op%ep(ix),2)
     sa = Op%states(Op%ep(ix),3)
              
     Ml = la
     Ms = sa 
 
     if ( Ml .ne. Op%Mltarg ) cycle
     if ( Ms .ne. Op%Mstarg ) cycle   
     sps = sps+1
  end do

     ! scalar operator
  tps = 0 

  do a = 1,p-1
     do b = a+1,p
        do i = 1, h

           p1=op%ep(a)
           p2=op%ep(b)
           h1=op%eh(i)
           
           li = Op%states(h1,2)
           si = Op%states(h1,3)
           la = Op%states(p1,2)
           sa = Op%states(p1,3)
           lb = Op%states(p2,2)
           sb = Op%states(p2,3)
           
           
           Ml = la+lb-li
           Ms = sa+sb-si 
           if ( Ml .ne. Op%Mltarg ) cycle
           if ( Ms .ne. Op%Mstarg ) cycle   
           tps = tps+ 1

        end do
     end do
  end do
      
 
  N = sps + tps ! number of p and pph SDs 
  
  ido = 0  ! status integer is 0 at start
  BMAT = 'I' ! standard eigenvalue problem (N for generalized) 
  which = 'SA' ! compute smallest eigenvalues in magnitude ('SA') is algebraic. 
  tol = 0.0E+00 ! error tolerance? (wtf zero?) 
  info = 0
  ncv = 2*nev ! number of lanczos vectors I guess
  if (ncv > N ) STOP 'NCV greater than N'
  lworkl = ncv*(ncv+8) 
  allocate(V(N,NCV),workl(lworkl))
  LDV = N  
  ishift = 1
  mxiter = 500 
  mode = 1

  allocate(eigs(N),resid(N),work(10*N),workD(3*N)) 
  
  iparam(1) = ishift
  iparam(3) = mxiter
  iparam(7) = mode
  i = 0

  do 
     ! so V is the krylov subspace matrix that is being diagonalized
     ! it does not need to be initialized, so long as you have the other 
     ! stuff declared right, the code should know this.
     
     call dsaupd ( ido, bmat, N, which, nev, tol, resid, &
      ncv, v, ldv, iparam, ipntr, workd, workl, &
      lworkl, info )
     
     ! The actual matrix only gets multiplied with the "guess" vector in "matvec_prod" 
     call progress_bar( i )
     i=i+1 
     
     if ( ido /= -1 .and. ido /= 1 ) then
        exit
     end if

     call matvec_1p_attached(N,dm,OP,workd(ipntr(1)), workd(ipntr(2)) ) 
  
  end do
  
  write(6,*) 
 
  ! the ritz values are out of order right now. Need to do post
  ! processing to fix this, and get the eigenvectors
  rvec= .true. 
  howmny = 'A'
  
  allocate(selct(NCV)) 
  allocate(Z(N,NEV)) 
  ldz = N  
  call dseupd( rvec, howmny, selct, d, Z, ldv, sigma, &
      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
      iparam, ipntr, workd, workl, lworkl, info )
  
  ! right now Z contains the eigenvectors in the columns
  ! d contains the eigenvalues in the same order. 
 
  do i = 1, NEV 
     norm_1p(i) = sum(Z(1:sps,i)**2) 
  end do     
end subroutine LANCZOS_1p_ATTACHED
!==================================================================
!==================================================================
subroutine LANCZOS_1h_REMOVED(OP,nev,dm,d,norm_1p)    
  
  integer :: N ,dm
  integer,intent(in) :: nev
  type(full_ham) :: op
  real(8),dimension(nev,dm) :: Vecs
  real(8),allocatable,dimension(:) :: workl,eigs,resid,work,workD
  real(8),allocatable,dimension(:,:) :: V,Z
  real(8),dimension(nev) :: d,norm_1p
  integer :: i,j,ix,jx,lwork,info,ido,ncv,ldv,iparam(11),ipntr(11),q,II,JJ
  integer :: ishift,mxiter,nb,nconv,mode,np,lworkl,ldz,p,h,sps,tps,jp,jh
  integer :: li,lj,lb,la,si,sj,sa,sb,Ml,Ms,h1,h2,p1,p2,a,b
  real(8) ::  x,tol,y,sigma,t1,t2
  character(1) :: BMAT,HOWMNY
  character(2) :: which,ncvstr
  logical :: rvec
  logical,allocatable,dimension(:) :: selct

  h = OP%Nbody !holes
  p = OP%Msp-h  !particles
  
  sps = 0 
  do ix = 1,h

     la = Op%states(Op%eh(ix),2)
     sa = Op%states(Op%eh(ix),3)
              
     Ml = -la
     Ms = -sa 
 
     if ( Ml .ne. Op%Mltarg ) cycle
     if ( Ms .ne. Op%Mstarg ) cycle   
     sps = sps+1
  end do

     ! scalar operator
  tps = 0 

  do a = 1,h-1
     do b = a+1,h
        do i = 1, p

           p1=op%eh(a)
           p2=op%eh(b)
           h1=op%ep(i)
           
           li = Op%states(h1,2)
           si = Op%states(h1,3)
           la = Op%states(p1,2)
           sa = Op%states(p1,3)
           lb = Op%states(p2,2)
           sb = Op%states(p2,3)
           
           
           Ml = -la-lb+li
           Ms = -sa-sb+si 
           if ( Ml .ne. Op%Mltarg ) cycle
           if ( Ms .ne. Op%Mstarg ) cycle   
           tps = tps+ 1

        end do
     end do
  end do
      
 
  N = sps + tps ! number of p and pph SDs 

  ido = 0  ! status integer is 0 at start
  BMAT = 'I' ! standard eigenvalue problem (N for generalized) 
  which = 'SA' ! compute smallest eigenvalues in magnitude ('SA') is algebraic. 
  tol = 0.0E+00 ! error tolerance? (wtf zero?) 
  info = 0
  ncv = 2*nev ! number of lanczos vectors I guess
  if (ncv > N ) STOP 'NCV greater than N'
  lworkl = ncv*(ncv+8) 
  allocate(V(N,NCV),workl(lworkl))
  LDV = N  
  ishift = 1
  mxiter = 500 
  mode = 1

  allocate(eigs(N),resid(N),work(10*N),workD(3*N)) 
  
  iparam(1) = ishift
  iparam(3) = mxiter
  iparam(7) = mode
  i = 0
  
  write(ncvstr,'(I2)') ncv
  
  do 
     ! so V is the krylov subspace matrix that is being diagonalized
     ! it does not need to be initialized, so long as you have the other 
     ! stuff declared right, the code should know this.
     
     call dsaupd ( ido, bmat, N, which, nev, tol, resid, &
      ncv, v, ldv, iparam, ipntr, workd, workl, &
      lworkl, info )

     
     ! The actual matrix only gets multiplied with the "guess" vector in "matvec_prod" 
     
     call progress_bar( i )
     
     i=i+1 

     if ( ido /= -1 .and. ido /= 1 ) then
        exit
     end if

     call matvec_1h_removed(N,dm,OP,workd(ipntr(1)), workd(ipntr(2)) ) 
  
  end do
  write(6,*) 
  
  ! the ritz values are out of order right now. Need to do post
  ! processing to fix this, and get the eigenvectors
  rvec= .true. 
  howmny = 'A'
  
  allocate(selct(NCV)) 
  allocate(Z(N,NEV)) 
  ldz = N  
  call dseupd( rvec, howmny, selct, d, Z, ldv, sigma, &
      bmat, n, which, nev, tol, resid, ncv, v, ldv, &
      iparam, ipntr, workd, workl, lworkl, info )
  
  ! right now Z contains the eigenvectors in the columns
  ! d contains the eigenvalues in the same order. 
 
  do i = 1, NEV 
     norm_1p(i) = sum(Z(1:sps,i)**2) 
  end do     
end subroutine LANCZOS_1h_REMOVED
!===========================================================================
!===========================================================================
subroutine matvec_prod(N,OP,Q_op,Qout,w1,w2,v,w) 
  implicit none 
  
  integer :: N ,q,a,b,c,d,i,j,k,l,Jtot,p1,p2,h1,h2,ix,jx
  integer :: Ml,Ms,h,p,si,sj,sa,sb,li,la,lb,lj,II,JJ
  real(8) :: sm,ass
  type(full_ham) :: OP, Q_op ,Qout,w1,w2
  real(8),dimension(N) :: v,w 
  real(8) :: coef9,dfact0

  ! the name says mat-vec product, but really this
  ! is a commutator. 
  
  ! the commutator here is equivalent to the matrix-vector product 
  ! with only connected diagrams retained. 

  ! FIRST WE NEED TO CONVERT v TO a (SQ_OP) variable
 
  h = OP%Nbody !holes
  p = OP%Msp-h  !particles
  
  Q_op%fph = 0.d0
  Q_op%mat(1)%Vpphh = 0.d0 
  Qout%fph = 0.d0
  Qout%mat(1)%Vpphh = 0.d0 

  call unwrap(v,Q_op,N)
  ! now we have two full_ham operators which can be used with my commutator expressions. Noice. 
  
  do ix = 1,p
     do jx = 1,h
                          
         li = Op%states(Op%eh(jx),2)
         si = Op%states(Op%eh(jx),3)
         la = Op%states(Op%ep(ix),2)
         sa = Op%states(Op%ep(ix),3)
              
         Ml = la-li
         Ms = sa-si 
         
         if ( Ml .ne. Op%Mltarg ) cycle
         if ( Ms .ne. Op%Mstarg ) cycle   
         
        ! print*, ix,jx 
         Qout%fph(ix,jx) = HQ_comm_1b(ix+h,jx,Op,Q_op) 
         
     end do
  end do
  
 ! call print_matrix(Qout%fph) 
 ! stop
     ! scalar operator

  II = 1
 
  do a = 1,p-1
     do b = a+1,p
        
        JJ = 0
        do i = 1, h-1
           do j = i+1,h
      
              p1=op%ep(a)
              p2=op%ep(b)
              h1=op%eh(i)
              h2=op%eh(j)
                           
              li = Op%states(h1,2)
              si = Op%states(h1,3)
              la = Op%states(p1,2)
              sa = Op%states(p1,3)
              lj = Op%states(h2,2)
              sj = Op%states(h2,3)
              lb = Op%states(p2,2)
              sb = Op%states(p2,3)
           
              Ml = la+lb-li-lj
              Ms = sa+sb-si-sj 
              JJ = JJ + 1            
            
              if ( Ml .ne. Op%Mltarg ) cycle
              if ( Ms .ne. Op%Mstarg ) cycle   
              
              !print*, HQ_comm_2b(a+h,b+h,i,j,Op,Q_op) , HQ_comm_2b(b+h,a+h,i,j,Op,Q_op),a+h,b+h,i,j
              
              Qout%mat(1)%Vpphh(II,JJ) =  HQ_comm_2b(a+h,b+h,i,j,Op,Q_op)
              
              
           end do 
        end do 
        II = II + 1 
     end do
  end do
  
  ! Okay, now we have the "matrix vector product" So lets put it back in vector form:
  call rewrap(w,Qout,N) 
!  print*, w 
!  stop
end subroutine matvec_prod
!======================================================================================
!======================================================================================
subroutine matvec_1p_attached(N,dm,OP,v,w) 
  implicit none 
  
  integer :: N ,q,a,b,c,d,i,j,k,l,Jtot,p1,p2,h1,h2,ix,jx
  integer :: Ml,Ms,h,p,si,sj,sa,sb,li,la,lb,lj,II,JJ,dm
  real(8) :: sm
  type(full_ham) :: OP
  real(8),dimension(N) :: v,w 
  real(8) :: coef9,dfact0
  real(8),allocatable,dimension(:) :: vec,vec_out

  ! the name says mat-vec product, but really this
  ! is a commutator. 
  
  ! the commutator here is equivalent to the matrix-vector product 
  ! with only connected diagrams retained. 

  ! FIRST WE NEED TO CONVERT v TO a (SQ_OP) variable
 
  h = OP%Nbody !holes
  p = OP%Msp-h  !particles
  allocate(vec(dm),vec_out(dm)) 

  vec = 0.d0 
  call unwrap_1p(v,vec,OP,N,dm)
  ! now we have two full_ham operators which can be used with my commutator expressions. Noice. 

  do ix = 1,p                          
     
     la = Op%states(Op%ep(ix),2)
     sa = Op%states(Op%ep(ix),3)
              
     Ml = la
     Ms = sa 
         
     if ( Ml .ne. Op%Mltarg ) cycle
     if ( Ms .ne. Op%Mstarg ) cycle   
         
     sm = HQ_comm_1p(ix+h,Op,vec)
     
     call place_X1(ix+h,vec_out,OP,sm) 
          
  end do

  do a = 1,p-1
     do b = a+1,p
        do i = 1, h
      
           p1=op%ep(a)
           p2=op%ep(b)
           h1=op%eh(i)

                           
           li = Op%states(h1,2)
           si = Op%states(h1,3)
           la = Op%states(p1,2)
           sa = Op%states(p1,3)
           lb = Op%states(p2,2)
           sb = Op%states(p2,3)
           
           Ml = la+lb-li
           Ms = sa+sb-si 
            
           if ( Ml .ne. Op%Mltarg ) cycle
           if ( Ms .ne. Op%Mstarg ) cycle   
           
           sm = HQ_comm_2p1h(a+h,b+h,i,Op,vec)

          call place_X3(a+h,b+h,i,vec_out,OP,sm) 

        end do
     end do
  end do
 
  ! Okay, now we have the "matrix vector product" So lets put it back in vector form:
  call rewrap_1p(w,vec_out,OP,N,dm)
 
end subroutine matvec_1p_attached
!======================================================================================
!======================================================================================
subroutine matvec_1h_removed(N,dm,OP,v,w) 
  implicit none 
  
  integer :: N ,q,a,b,c,d,i,j,k,l,Jtot,p1,p2,h1,h2,ix,jx
  integer :: Ml,Ms,h,p,si,sj,sa,sb,li,la,lb,lj,II,JJ,dm
  real(8) :: sm
  type(full_ham) :: OP
  real(8),dimension(N) :: v,w 
  real(8) :: coef9,dfact0
  real(8),allocatable,dimension(:) :: vec,vec_out

  ! the name says mat-vec product, but really this
  ! is a commutator. 
  
  ! the commutator here is equivalent to the matrix-vector product 
  ! with only connected diagrams retained. 

  ! FIRST WE NEED TO CONVERT v TO a (SQ_OP) variable
 
  h = OP%Nbody !holes
  p = OP%Msp-h  !particles
  allocate(vec(dm),vec_out(dm)) 
 
  vec = 0.d0 
  vec_out = 0.d0
  call unwrap_1h(v,vec,OP,N,dm)
  ! now we have two full_ham operators which can be used with my commutator expressions. Noice. 

  
  do ix = 1,h                          
     
     la = Op%states(Op%eh(ix),2)
     sa = Op%states(Op%eh(ix),3)
              
     Ml = -la
     Ms = -sa 
         
     if ( Ml .ne. Op%Mltarg ) cycle
     if ( Ms .ne. Op%Mstarg ) cycle   
         
     sm = HQ_comm_1h(ix,Op,vec)
     
     call place_X1(ix,vec_out,OP,sm) 
          
  end do

  do a = 1,h-1
     do b = a+1,h
        do i = 1, p
      
           p1=op%eh(a)
           p2=op%eh(b)
           h1=op%ep(i)

                           
           li = Op%states(h1,2)
           si = Op%states(h1,3)
           la = Op%states(p1,2)
           sa = Op%states(p1,3)
           lb = Op%states(p2,2)
           sb = Op%states(p2,3)
           
           Ml = -la-lb+li
           Ms = -sa-sb+si 
            
           if ( Ml .ne. Op%Mltarg ) cycle
           if ( Ms .ne. Op%Mstarg ) cycle   
           
           sm = HQ_comm_2h1p(a,b,i+h,Op,vec)

          call place_X3(a,b,i+h,vec_out,OP,sm) 

        end do
     end do
  end do
  
  ! Okay, now we have the "matrix vector product" So lets put it back in vector form:

  call rewrap_1h(w,vec_out,OP,N,dm)
  
end subroutine matvec_1h_removed
!======================================================================================
!======================================================================================
subroutine unwrap( v, AX ,N) 
  implicit none 
  
  integer :: N ,i,j,a,b,g, II,JJ, p,h,q,IX,JX,Ml,ms
  integer :: h1,h2,p1,p2,li,lj,la,lb,si,sj,sa,sb
  real(8),dimension(N) :: v
  type(full_ham) :: AX 
  
  g = 1
  
  h = AX%Nbody !holes
  p = AX%Msp-h  !particles
 
  do ix = 1,p
     do jx = 1,h
                          
         li = AX%states(AX%eh(jx),2)
         si = AX%states(AX%eh(jx),3)
         la = AX%states(AX%ep(ix),2)
         sa = AX%states(AX%ep(ix),3)
              
         Ml = la-li
         Ms = sa-si 
         
         if ( Ml .ne. AX%Mltarg ) cycle
         if ( Ms .ne. AX%Mstarg ) cycle   

         AX%fph(IX,JX) = v(g) 
         g = g + 1
         
     end do
  end do

  do a = 1,p
     do b = a+1,p
        do i = 1, h
           do j = i+1,h
      
              p1=AX%ep(a)
              p2=AX%ep(b)
              h1=AX%eh(i)
              h2=AX%eh(j)
              
              
              li = AX%states(h1,2)
              si = AX%states(h1,3)
              la = AX%states(p1,2)
              sa = AX%states(p1,3)
              lj = AX%states(h2,2)
              sj = AX%states(h2,3)
              lb = AX%states(p2,2)
              sb = AX%states(p2,3)
           

              Ml = la+lb-li-lj
              Ms = sa+sb-si-sj 
              if ( Ml .ne. AX%Mltarg ) cycle
              if ( Ms .ne. AX%Mstarg ) cycle   

              II = fermionic_tp_index(a,b,p)
              JJ = fermionic_tp_index(i,j,h) 
              AX%mat(1)%Vpphh(II,JJ) = v(g)
              g = g + 1

           end do 
        end do 
     end do
  end do
  
end subroutine unwrap
!======================================================================================
!======================================================================================
subroutine unwrap_1p(v,vec,AX,N,dm) 
  implicit none 
  
  integer :: N ,i,j,a,b,g, II,JJ, p,h,q,IX,JX,Ml,ms
  integer :: h1,h2,p1,p2,li,lj,la,lb,si,sj,sa,sb,dm
  real(8),dimension(N) :: v
  real(8),dimension(dm) :: vec
  type(full_ham) :: AX 
  
  g = 1
  
  h = AX%Nbody !holes
  p = AX%Msp-h  !particles
 
  do ix = 1,p
     
     la = AX%states(AX%ep(ix),2)
     sa = AX%states(AX%ep(ix),3)
              
     Ml = la
     Ms = sa 
         
     if ( Ml .ne. AX%Mltarg ) cycle
     if ( Ms .ne. AX%Mstarg ) cycle   
     
     call place_X1(ix+h,vec,AX,v(g))
     g = g + 1
         
  end do
  

  do a = 1,p
     do b = a+1,p
        do i = 1, h
      
           p1=AX%ep(a)
           p2=AX%ep(b)
           h1=AX%eh(i)
           
           li = AX%states(h1,2)
           si = AX%states(h1,3)
           la = AX%states(p1,2)
           sa = AX%states(p1,3)
           lb = AX%states(p2,2)
           sb = AX%states(p2,3)
           
           
           Ml = la+lb-li
           Ms = sa+sb-si 
           if ( Ml .ne. AX%Mltarg ) cycle
           if ( Ms .ne. AX%Mstarg ) cycle   
      
           call place_X3(a+h,b+h,i,vec,AX,v(g))
           g = g + 1

        end do 
     end do
  end do
  
end subroutine unwrap_1p
!======================================================================================
!======================================================================================
subroutine unwrap_1h(v,vec,AX,N,dm) 
  implicit none 
  
  integer :: N ,i,j,a,b,g, II,JJ, p,h,q,IX,JX,Ml,ms
  integer :: h1,h2,p1,p2,li,lj,la,lb,si,sj,sa,sb,dm
  real(8),dimension(N) :: v
  real(8),dimension(dm) :: vec
  type(full_ham) :: AX 
  
  g = 1
  
  h = AX%Nbody !holes
  p = AX%Msp-h  !particles
 
  do ix = 1,h
     
     la = AX%states(AX%eh(ix),2)
     sa = AX%states(AX%eh(ix),3)
              
     Ml = -la
     Ms = -sa 
         
     if ( Ml .ne. AX%Mltarg ) cycle
     if ( Ms .ne. AX%Mstarg ) cycle   
     
     call place_X1(ix,vec,AX,v(g))
     g = g + 1
         
  end do
  

  do a = 1,h
     do b = a+1,h
        do i = 1, p
      
           p1=AX%eh(a)
           p2=AX%eh(b)
           h1=AX%ep(i)
           
           li = AX%states(h1,2)
           si = AX%states(h1,3)
           la = AX%states(p1,2)
           sa = AX%states(p1,3)
           lb = AX%states(p2,2)
           sb = AX%states(p2,3)
           
           
           Ml = -la-lb+li
           Ms = -sa-sb+si 
           if ( Ml .ne. AX%Mltarg ) cycle
           if ( Ms .ne. AX%Mstarg ) cycle   
      
           call place_X3(a,b,i+h,vec,AX,v(g))
           g = g + 1

        end do 
     end do
  end do
  
end subroutine unwrap_1h
!============================================================================  
!============================================================================
subroutine rewrap( v, AX ,N ) 
  implicit none 
  
  integer :: N ,i,j,a,b,g, II,JJ, p,h,q,IX,JX,Ml,Ms
  integer :: h1,h2,p1,p2,li,lj,la,lb,si,sj,sa,sb
  real(8),dimension(N) :: v
  type(full_ham) :: AX 
  
  g = 1
  
  h = AX%Nbody !holes
  p = AX%Msp-h  !particles
 
  do ix = 1,p
     do jx = 1,h
                          
         li = AX%states(AX%eh(jx),2)
         si = AX%states(AX%eh(jx),3)
         la = AX%states(AX%ep(ix),2)
         sa = AX%states(AX%ep(ix),3)
              
         Ml = la-li
         Ms = sa-si 
         
         if ( Ml .ne. AX%Mltarg ) cycle
         if ( Ms .ne. AX%Mstarg ) cycle   

         v(g) = AX%fph(IX,JX)
         g = g + 1
         
     end do
  end do

  do a = 1,p
     do b = a+1,p
        do i = 1, h
           do j = i+1,h
      
              p1=AX%ep(a)
              p2=AX%ep(b)
              h1=AX%eh(i)
              h2=AX%eh(j)
              
              
              li = AX%states(h1,2)
              si = AX%states(h1,3)
              la = AX%states(p1,2)
              sa = AX%states(p1,3)
              lj = AX%states(h2,2)
              sj = AX%states(h2,3)
              lb = AX%states(p2,2)
              sb = AX%states(p2,3)
           

              Ml = la+lb-li-lj
              Ms = sa+sb-si-sj 
              if ( Ml .ne. AX%Mltarg ) cycle
              if ( Ms .ne. AX%Mstarg ) cycle   

              II = fermionic_tp_index(a,b,p)
              JJ = fermionic_tp_index(i,j,h) 
              v(g) = AX%mat(1)%Vpphh(II,JJ)
              g = g + 1

           end do 
        end do 
     end do
  end do

end subroutine 
!==============================================================
!==============================================================
subroutine rewrap_1p(v,vec,AX,N,dm) 
  implicit none 
  
  integer :: N ,i,j,a,b,g, II,JJ, p,h,q,IX,JX,Ml,ms
  integer :: h1,h2,p1,p2,li,lj,la,lb,si,sj,sa,sb,dm
  real(8),dimension(N) :: v
  real(8),dimension(dm) :: vec
  type(full_ham) :: AX 
  
  g = 1
  
  h = AX%Nbody !holes
  p = AX%Msp-h  !particles
 
  do ix = 1,p
     
     la = AX%states(AX%ep(ix),2)
     sa = AX%states(AX%ep(ix),3)
              
     Ml = la
     Ms = sa 
         
     if ( Ml .ne. AX%Mltarg ) cycle
     if ( Ms .ne. AX%Mstarg ) cycle   

     v(g) = get_X1(ix+h,vec,AX)
     g = g + 1
         
  end do
  

  do a = 1,p
     do b = a+1,p
        do i = 1, h
      
           p1=AX%ep(a)
           p2=AX%ep(b)
           h1=AX%eh(i)
           
           li = AX%states(h1,2)
           si = AX%states(h1,3)
           la = AX%states(p1,2)
           sa = AX%states(p1,3)
           lb = AX%states(p2,2)
           sb = AX%states(p2,3)
           
           
           Ml = la+lb-li
           Ms = sa+sb-si 
           if ( Ml .ne. AX%Mltarg ) cycle
           if ( Ms .ne. AX%Mstarg ) cycle   
           
           v(g)=get_X3(a+h,b+h,i,vec,AX)
           g = g + 1

        end do 
     end do
  end do
  
end subroutine rewrap_1p
!==============================================================
!==============================================================
subroutine rewrap_1h(v,vec,AX,N,dm) 
  implicit none 
  
  integer :: N ,i,j,a,b,g, II,JJ, p,h,q,IX,JX,Ml,ms
  integer :: h1,h2,p1,p2,li,lj,la,lb,si,sj,sa,sb,dm
  real(8),dimension(N) :: v
  real(8),dimension(dm) :: vec
  type(full_ham) :: AX 
  
  g = 1
  
  h = AX%Nbody !holes
  p = AX%Msp-h  !particles
 
  do ix = 1,h
     
     la = AX%states(AX%eh(ix),2)
     sa = AX%states(AX%eh(ix),3)
              
     Ml = -la
     Ms = -sa 
         
     if ( Ml .ne. AX%Mltarg ) cycle
     if ( Ms .ne. AX%Mstarg ) cycle   

     v(g) = get_X1(ix,vec,AX)
     g = g + 1
         
  end do
  

  do a = 1,h
     do b = a+1,h
        do i = 1, p
      
           p1=AX%eh(a)
           p2=AX%eh(b)
           h1=AX%ep(i)
           
           li = AX%states(h1,2)
           si = AX%states(h1,3)
           la = AX%states(p1,2)
           sa = AX%states(p1,3)
           lb = AX%states(p2,2)
           sb = AX%states(p2,3)
           
           
           Ml = -la-lb+li
           Ms = -sa-sb+si 
           if ( Ml .ne. AX%Mltarg ) cycle
           if ( Ms .ne. AX%Mstarg ) cycle   
           
           v(g)=get_X3(a,b,i+h,vec,AX)
           g = g + 1

        end do 
     end do
  end do
  
end subroutine rewrap_1h
!============================================================================  
!============================================================================


subroutine progress_bar( step  )  
  implicit none
  
  integer :: step
 
  if ( step .ne. 0) then 
     ! hold the backspace key down
     write(6,'(A)',advance='no') char(8)//char(8)//char(8)//char(8)//char(8) &
       //char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8)&
       //char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8)//char(8) &
       //char(8)//char(8)//char(8)//char(8)
     flush 6
  end if 

  write(6,'((A19),(I5))',advance='no') 'Lanczos iteration: ',step
  flush 6 

end subroutine 

end module
