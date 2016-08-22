module HFQDmod
  implicit none 
  
!!! contains the hartree fock machinery, constructs the basis and builds 
!!! a tp basis

contains  
  
subroutine construct_two_particle_HF_BASIS( n , hw , emax , record ,coefs,HCC,pert) 
  use ME_general
  implicit none 
  
  !!!  ( number of particles , HO spacing , Max HO energy , size of 2part basis 
  !!!     full_ham storage array, ( mbas,mbas ) storage array ) 
  
  !!! record holds all of the matrix elements in blocks, sorted by particle-hole scheme
  
  !!! coefs stores the sp eigenvectors in columns , ordered by  scheme in arrange states
  
  real(8), parameter :: al=1.d0, bet=0.d0
  integer :: i,j,k,l,n,emax,Mbas,gbas,ghol,gpart,info,cnt,q,qt,ist,kron_del,k1,l1
  integer :: a,b,c,d,ii,jj,Nmax,i1,j1,n1,n2,n3,n4,m1,m2,m3,m4,s1,s2,s3,s4,iix,jjx
  integer,allocatable,dimension(:,:) ::  qnums
  integer,allocatable,dimension(:) :: ordering
  real(8),dimension(emax*(emax+1),emax*(emax+1)) :: coefs
  real(8),allocatable,dimension(:,:) :: T,V,H,F,den,Vstor
  real(8),allocatable,dimension(:) :: eold,ord,work,eig
  real(8) :: hw,crit,MBPT2,HF_E2,eHF,Vdir,Vexch,v_int
  real(8), optional :: pert
  character(5) :: hwstr,nstr,emaxstr
  character(35) :: fname
  type(full_ham) :: record,ints
  type(cc_mat) :: HCC 
  type(full_sp_block_mat) :: Fblock
 

!================= SP ALLOCATION BEGIN ========================= 

  !!! Mbas is the number of sp basis functions
  Mbas=emax*(emax+1)
  
  !!! allocate memory to arrays
  
  allocate(qnums(Mbas,3))
  allocate(T(Mbas,Mbas))
  allocate(V(Mbas,Mbas))
  allocate(H(Mbas,Mbas))
  allocate(F(Mbas,Mbas))
  allocate(den(Mbas,Mbas))
  !allocate(ints(Mbas,Mbas,Mbas,Mbas))
  allocate(eold(Mbas))
  allocate(eig(Mbas))
  allocate(ord(Mbas))
  allocate(work(10*Mbas))
  allocate(ordering(Mbas/2))
  
  
  !!! gbas is the number of tp basis functions
  gbas = Mbas*(Mbas-1)/2
  !!! ghol is the number of hole-hole states
  ghol = n*(n-1)/2
  !!! gpart is the number of particle-particle states
  gpart = (Mbas-n)*(Mbas-n-1)/2




!!!================ HF ITERATOR BEGIN ==========================

  !build states ordering, calculate 1-body potential 
  !(unchanged in HF iterations)
  
  call arrange_states(emax,qnums,ordering) 
  
  call calc_h0(T,hw,qnums) 

  ! figure out what the blocks are ahead of time
  ! so less storage is needed
  do i=1,mbas
     eig(i) = T(i,i) 
  end do

  call reindex(Mbas, gbas, 0, 0, gbas, qnums, emax,eig,ints)
  call convertF(T,ints,mbas,qnums,eig)
 
  write(hwstr,'(f5.2)') hw
  write(nstr, '(i5)' ) n
  write(emaxstr, '(i5)' ) emax
  
     !!!! get matrix elements
  
  hwstr=adjustl(hwstr)
  nstr=adjustl(nstr)
  emaxstr=adjustl(emaxstr)

 ! fname='../TBMEfiles/TBME_'//trim(hwstr)//'_'//trim(emaxstr)//'.dat'
  fname='../TBMEfiles/ME.bin' 
  
  !fname=trim(fname) 

  open(unit=37,file=fname,form='unformatted',access='stream')

  Nmax = (emax+1)*(emax)/2 

  allocate(vstor(Nmax**2,Nmax**2)) 
  vstor = 0.d0
  do 
     read(37,iostat=ist) n1,m1,n2,m2,n3,m3,n4,m4,Vdir 
     
     if (ist >0) stop 'failfile' 
     if (ist < 0) exit
     
     if (2*n1+abs(m1) > emax-1) cycle
     if (2*n2+abs(m2) > emax-1) cycle
     if (2*n3+abs(m3) > emax-1) cycle
     if (2*n4+abs(m4) > emax-1) cycle

     a = nm_index(n1,m1)
     b = nm_index(n2,m2)
     c = nm_index(n3,m3)
     d = nm_index(n4,m4)     

     II = ab_index(a,b,Nmax)
     JJ = ab_index(c,d,Nmax) 
     Vstor(II,JJ) = Vdir*sqrt(hw)
     Vstor(JJ,II) = Vdir*sqrt(hw)
  end do 
     
  close(37)
  
  do q=1,ints%nblock
     allocate(ints%mat(q)%Vpppp(ints%mat(q)%npp,ints%mat(q)%npp))
     do i=1,ints%mat(q)%npp
        i1 = ints%mat(q)%qnpp(i,1)
        j1 = ints%mat(q)%qnpp(i,2)
        
        n1 = qnums(i1,1)
        m1 = qnums(i1,2)
        s1 = qnums(i1,3)
        
        n2 = qnums(j1,1)
        m2 = qnums(j1,2)
        s2 = qnums(j1,3)
        
        a = nm_index(n1,m1)
        b = nm_index(n2,m2)
        II = ab_index(a,b,Nmax)
        IIx = ab_index(b,a,Nmax) 
        
        do j=1,ints%mat(q)%npp
           k1 = ints%mat(q)%qnpp(j,1)
           l1 = ints%mat(q)%qnpp(j,2)
           
           n3 = qnums(k1,1)
           m3 = qnums(k1,2)
           s3 = qnums(k1,3)
           
           n4 = qnums(l1,1)
           m4 = qnums(l1,2)
           s4 = qnums(l1,3)
           
           c = nm_index(n3,m3)
           d = nm_index(n4,m4)
           
           JJ = ab_index(c,d,Nmax)
           Vdir = Vstor(II,JJ) 
           JJx = ab_index(d,c,Nmax)
           Vexch = Vstor(II,JJx)            
           if (abs(Vexch)<1e-10)Vexch = Vstor(IIx,JJ)  
            
           ints%mat(q)%Vpppp(i,j) = &
                kron_del(s1,s3)*kron_del(s2,s4)*Vdir &
               -kron_del(s1,s4)*kron_del(s2,s3)*Vexch
           
        end do
     end do
     if (present(pert)) then 
        ints%mat(q)%Vpppp=pert*ints%mat(q)%Vpppp
     end if
  end do

  close(37)

!!! matrix elements are now available

  eig=0.d0
  crit=1.d0
  
  coefs=0.d0
  do i=1,mbas
     coefs(i,i)=1.d0
  end do
  
  call build_block_matrix(qnums,emax,Fblock) 
  
  do while (crit > 1d-6) 
     
     call calcDen(den,coefs,n,Mbas,qnums,ordering)   
     
     call construct_V(V,ints,den,mbas)
     
     H=T+V
     eold=eig
     
     !! because H is about to be destroyed
     F=H
     call sort_into_blocks(H,Fblock,mbas,qnums)
     call diagonalize_blocks(Fblock)       
     call eigvecs_to_normal(H,Fblock,mbas,eig,qnums)      
    
     crit=0.d0

     do i=1,n

        crit=crit + abs( eold(i)-eig(i) )             
                  
     end do 

     coefs=H  
 
 end do 

!====================== TP ALLOCATION BEGIN ============================

  !!! switch to a two particle basis, now that we have the HF basis. 
  !!! deallocate all of the single particle stuff to make room

   call dgemm('N','N',mbas,mbas,mbas,al,F,mbas,coefs,mbas,bet,H,mbas)
   call dgemm('T','N',mbas,mbas,mbas,al,coefs,mbas,H,mbas,bet,F,mbas)

   call reindex(Mbas, gbas, n, ghol, gpart, qnums, emax,eig,record)
 
   call new_COEFS(coefs,eig,mbas,record)
   
   call new_VMAT(ints,record,mbas)

   call convertF(F,record,mbas,qnums,eig)

   call build_cross_coupled(record,HCC,qnums)
   
   deallocate(qnums)
   do q=1,ints%nblock
      deallocate(ints%mat(q)%Vpppp)
   end do 
   deallocate(ints%mat) 
   deallocate(H)
   deallocate(V)
   deallocate(T)
   deallocate(den)
   deallocate(F)
   deallocate(eig)


end subroutine 
!==========================================
!==========================================
end module
!==========================================
!==========================================
subroutine construct_V(V,ints,den,m)
  use ME_general
  implicit none 
  
  integer :: i,j,k,l,m,a,b,c,d
  real(8),dimension(m,m) :: V,den
  type(full_ham) :: ints

  V=0.d0

  !!! for each i,k pair, we break into four seperate sums
  !!! because matrix elements <i,j|V|k,l> only exist
  !!! for  i .le. j  and k .le. l

  do i=1,m
     do k=1,m
        
        do j=1,m
           do  l= 1,m 
              

              V(i,k) = V(i,k) + den(l,j) * v_elem(j,i,l,k,ints) 

        
           end do 
        end do 
           
     end do 
  end do 
  
end subroutine 
!==============================================
real(8) function HF_E2(rec)
  use ME_general  
  !!! hartree fock energy using 
  !!! the two-particle basis
  implicit none 
  
  integer :: i,j
  real(8) :: sm,sm2
  type(full_ham) :: rec
  
  sm=0.d0
  sm2=0.d0

  do i=1,rec%nblock
     
     do j=1,rec%mat(i)%nhh
        
        sm = sm + rec%mat(i)%Vhhhh(j,j)
     end do
  end do 

  do i=1,rec%nbody
     sm2=sm2+rec%eig(rec%eh(i))
  end do 
 
  !side-effect: E0 is added to the array
  rec%E0=sm2-sm

  HF_E2=sm2-sm
  
  
end function 
!==============================================
real(8) function MBPT2(rec)
  use ME_general
  !!! second order correction in MBPT
  !!! uses the two-particle basis
  implicit none 
  
  real(8),parameter :: al=1.d0,bet=0.d0
  integer :: i,j,k,nh,np,II,JJ,N,p1,p2,h1,h2
  real(8) ::  sm,num,denom,f12
  type(full_ham) :: rec
  
  N = rec%nbody
  sm=0.d0
  do i=1,rec%nblock
     nh = rec%mat(i)%nhh
     np = rec%mat(i)%npp
     
     do II = 1,np 
        p1 = rec%mat(i)%qnpp(II,1)
        p2 = rec%mat(i)%qnpp(II,2)

        f12 = rec%fpp(p1-n,p1-n) + rec%fpp(p2-n,p2-n) 
        
        do JJ = 1, nh 
           h1 = rec%mat(i)%qnhh(JJ,1)
           h2 = rec%mat(i)%qnhh(JJ,2)

           denom = rec%fhh(h1,h1) + rec%fhh(h2,h2) - f12 
           
           sm = sm + rec%mat(i)%Vpphh(II,JJ)**2/denom 
        end do
     end do

  end do 
  
  MBPT2=sm

end function 
!============================================
