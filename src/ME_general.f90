module ME_general 
  implicit none 

!!! for SP states:
!!! STATES ORDER: ordered by subroutine arrange_states algorithm,
!!! which creates a block diagonal fock matrix

!!! ENERGY ORDER: a more natural ordering from lowest to highest energy

!!! NOTE: STATES ORDER is important because LAPACK cannot comprehend blocks
!!! unless they are ordered so that physical blocks appear in the diagonal. 

type block
   integer :: nhh,npp,nph,ntot
   !!! number of hh,pp,ph and total states
   integer,dimension(2) :: lam
   !!! ( ML, MS ) 
   integer,allocatable,dimension(:,:) :: qnhh,qnpp,qnph
   !!! tp to sp quantum number maps in STATES ORDER
   real(8),allocatable,dimension(:,:) :: CC,CB
   !!! C*C coefficent matrices, forward and backward to account for antisym.
   real(8),allocatable,dimension(:,:) :: Vpphh,Vppph,Vphhh,Vphph,Vhhhh,Vpppp
   !!! holds V elements in the order described by QN
   real(8),allocatable,dimension(:,:) :: Dpphh
   !!! same as Vpphh BUT WITH ENERGY DENOMINATORS
   real(8),allocatable,dimension(:) :: eig2,wrk
   !!! holds summed eigenvalues for tp states
end type block

type full_ham
   type(block),allocatable,dimension(:) :: mat
   !!! holds blocks in a vector array
   real(8),allocatable,dimension(:,:) :: fpp,fhh,fph,Xpp,Xhh,Xph,V3body
   !!! holds fock elements (IN ENERGY ORDER)
   real(8),allocatable,dimension(:) :: eig
   !!! holds eigenvalues ( IN STATES ORDER ) 
   integer,allocatable,dimension(:,:) :: states,threemap
   !!! maps state order number to quantum numbers
   integer,allocatable,dimension(:,:) :: i_array
   !!! maps sp indeces to tp indece and other stuff
   integer,allocatable,dimension(:) :: eh,ep,stoe
   !!! eh(i) gives the STATES ORDER position of the 
   !!! ith state in ENERGY ORDER
   real(8) :: E0,h3(220,220)
   !!! 0-body flow energy
   integer :: nblock,nbody,Msp,herm,neq 
   !!! number of blocks, bodies, sp states
   logical :: IMSRG3 
end type full_ham
!========================
type sp_block_mat
   real(8),allocatable,dimension(:,:) :: matrix!,eigvec
   real(8),allocatable,dimension(:) :: Eigval,extra
   integer,allocatable,dimension(:,:) :: labels
   integer,dimension(2) :: lmda
end type sp_block_mat
!========================
type full_sp_block_mat
   type(sp_block_mat),allocatable,dimension(:) :: blkM
   integer,allocatable,dimension(:) :: map
   integer :: blocks
end type full_sp_block_mat
!========================
contains
!========================================
subroutine write_matrix_pretty(ETA,fname) 
  implicit none 

  integer :: n,m,nh,np,nb,i,q
  character(2) :: nhs,nps,nbs
  character(*) :: fname
  type(full_ham) :: ETA
  
 ! fname = adjustl(fname) 
  open(unit=47,file = fname ) 
  
  n = ETA%nbody
  m = ETA%msp
  
  write(nhs,'(I2)') n
  write(nps,'(I2)') m-n
  nhs = adjustl(nhs)
  nps = adjustl(nps)
  
  write(47,*) 'fhh'
  write(47,*) 
  
  do i=1,n
     write(47,'('//trim(nhs)//'(f10.5))') ETA%fhh(i,:)
  end do 
  
  write(47,*)
  write(47,*) 'fpp'
  write(47,*) 
  
  do i=1,m-n
     write(47,'('//trim(nps)//'(f10.5))') ETA%fpp(i,:)
  end do 
  
  write(47,*)
  write(47,*) 'fph'
  write(47,*) 
  
  do i=1,m-n
     write(47,'('//trim(nhs)//'(f10.5))') ETA%fph(i,:)
  end do 
  
  write(47,*) 
  write(47,*)  '===================================================='
  
  do q = 1, ETA%nblock
     
     nh = ETA%mat(q)%nhh
     np = ETA%mat(q)%npp
     nb = ETA%mat(q)%nph
     write(nhs,'(I2)') nh
     write(nps,'(I2)') np
     write(nbs,'(I2)') nb
     nhs = adjustl(nhs)
     nps = adjustl(nps)
     nbs = adjustl(nbs)
     
 write(47,'(A,I2,A,I2)')'ML =',ETA%mat(q)%lam(1),' MS = ',ETA%mat(q)%lam(2)
     write(47,*)  '===================================================='
     write(47,*)
    
     if (nh > 0 ) then 
     write(47,*) 'Vhhhh'
     write(47,*) 
  
     
     do i=1,nh
        write(47,'('//trim(nhs)//'(f10.5))') ETA%mat(q)%Vhhhh(i,:)
     end do 
     
     write(47,*)
     end if 
     if (np > 0 ) then 
     write(47,*) 'Vpppp'
     write(47,*) 
     
     do i=1,np
        write(47,'('//trim(nps)//'(f10.5))') ETA%mat(q)%Vpppp(i,:)
     end do 
     
     write(47,*)
     end if 
     if (nb > 0 ) then 
     write(47,*) 'Vphph'
     write(47,*) 
     
     do i=1,nb
        write(47,'('//trim(nbs)//'(f10.5))') ETA%mat(q)%Vphph(i,:)
     end do 
     
     write(47,*)
     end if 
     if (nh*np > 0 ) then 
     write(47,*) 'Vpphh'
     write(47,*) 
     
     
     do i=1,np
        write(47,'('//trim(nhs)//'(f10.5))') ETA%mat(q)%Vpphh(i,:)
     end do 
     
     
     write(47,*)
     end if 
     if (np*nb > 0 ) then 
     write(47,*) 'Vppph'
     write(47,*) 
     
     do i=1,np
        write(47,'('//trim(nbs)//'(f10.5))') ETA%mat(q)%Vppph(i,:)
     end do 
     
     write(47,*)
     end if 
     if (nh*nb > 0 ) then 
     write(47,*) 'Vphhh'
     write(47,*) 
     
     do i=1,nb
        write(47,'('//trim(nhs)//'(f10.5))') ETA%mat(q)%Vphhh(i,:)
     end do 
     write(47,*)
     end if 
     write(47,*) '===================================================='
  end do 
  close(47) 
end subroutine 
!====================================================
subroutine export_hamiltonian(H) 
  implicit none 
  
  integer :: i,j,n,m,q,np,nh,nb
  type(full_ham) :: H

  n=H%nbody
  m=H%msp

  open(unit=20,file='hamiltonian.dat')
  
  write(20,*) H%E0
  do i=1,n
     do j=1,n
        write(20,*) H%fhh(i,j)
     end do 
  end do 
  
  do i=1,m-n
     do j=1,n
        write(20,*) H%fph(i,j)
     end do 
  end do 
  
  do i=1,m-n
     do j=1,m-n
        write(20,*) H%fpp(i,j)
     end do 
  end do 

  do q=1,H%nblock
     nh=H%mat(q)%nhh
     np=H%mat(q)%npp
     nb=H%mat(q)%nph
     
     do i=1,nh
        do j=1,nh
           write(20,*) H%mat(q)%Vhhhh(i,j)
        end do 
     end do 
     
     do i=1,np
        do j=1,np
           write(20,*) H%mat(q)%Vpppp(i,j)
        end do 
     end do 
     
     do i=1,np
        do j=1,nh
           write(20,*) H%mat(q)%Vpphh(i,j)
        end do 
     end do 
     
     do i=1,nb
        do j=1,nb
           write(20,*) H%mat(q)%Vphph(i,j)
        end do 
     end do 
     
     do i=1,nb
        do j=1,nh
           write(20,*) H%mat(q)%Vphhh(i,j)
        end do 
     end do 
     
     do i=1,np
        do j=1,nb
           write(20,*) H%mat(q)%Vppph(i,j)
        end do 
     end do 
     
   end do 

close(20)
end subroutine 
!=========================================
subroutine import_hamiltonian(H) 
  implicit none 
  
  integer :: i,j,n,m,q,np,nh,nb
  type(full_ham) :: H

  n=H%nbody
  m=H%msp

  open(unit=20,file='hamiltonian.dat')
  
  read(20,*) H%E0
  do i=1,n
     do j=1,n
        read(20,*) H%fhh(i,j)
     end do 
  end do 
  
  do i=1,m-n
     do j=1,n
        read(20,*) H%fph(i,j)
     end do 
  end do 
  
  do i=1,m-n
     do j=1,m-n
        read(20,*) H%fpp(i,j)
     end do 
  end do 

  do q=1,H%nblock
     nh=H%mat(q)%nhh
     np=H%mat(q)%npp
     nb=H%mat(q)%nph
     
     do i=1,nh
        do j=1,nh
           read(20,*) H%mat(q)%Vhhhh(i,j)
        end do 
     end do 
     
     do i=1,np
        do j=1,np
           read(20,*) H%mat(q)%Vpppp(i,j)
        end do 
     end do 
     
     do i=1,np
        do j=1,nh
           read(20,*) H%mat(q)%Vpphh(i,j)
        end do 
     end do 
     
     do i=1,nb
        do j=1,nb
           read(20,*) H%mat(q)%Vphph(i,j)
        end do 
     end do 
     
     do i=1,nb
        do j=1,nh
           read(20,*) H%mat(q)%Vphhh(i,j)
        end do 
     end do 
     
     do i=1,np
        do j=1,nb
           read(20,*) H%mat(q)%Vppph(i,j)
        end do 
     end do 
     
   end do 

close(20)
end subroutine         
!================================================
subroutine build_block_matrix(qn,emax,F) 
  implicit none
  
  integer :: m,emax,ml,ms,d,q,l,count_states
  integer,dimension(emax*(emax+1),3) :: qn
  type(full_sp_block_mat) :: F
 
  allocate(F%blkM(4*emax-2))
  allocate(F%map(4*emax-2))
  
  m=emax*(emax+1) 
  F%blocks= 4*emax-2
 
q = 1
do ms = -1,1,2 
  ml = 0 
  
  d = count_states(ml,ms,qn,m) 
  
  allocate(F%blkM(q)%matrix(d,d)) 
  allocate(F%blkM(q)%eigval(d)) 
  allocate(F%blkM(q)%extra(10*d)) 
 
  F%map(q) = d
  q=q+1

  do l = 1, emax - 1
   
    ml = -1*l
    d = count_states(ml,ms,qn,m) 
  
    allocate(F%blkM(q)%matrix(d,d)) 
    allocate(F%blkM(q)%eigval(d)) 
    allocate(F%blkM(q)%extra(10*d)) 
 
    F%map(q) = d
    q=q+1
    
    ml = l
    d = count_states(ml,ms,qn,m) 
  
    allocate(F%blkM(q)%matrix(d,d)) 
    allocate(F%blkM(q)%eigval(d)) 
    allocate(F%blkM(q)%extra(10*d)) 
 
    F%map(q) = d
    q=q+1
  end do 
end do 

end subroutine      
!=====================================================
subroutine sort_into_blocks(A,R,m) 
  !A is the matrix
  !Q is the block matrix holder
  implicit none
  
  integer :: m,i,q,j
  real(8),dimension(m,m) :: A
  type(full_sp_block_mat) :: R
  
  i=1
  do q=1,R%blocks
     R%blkM(q)%matrix = A(i:i+R%map(q)-1,i:i+R%map(q)-1)
     i=i+R%map(q)
  end do 
end subroutine 
!======================================================
subroutine eigvecs_to_normal(A,R,m,eig) 
  !A is the matrix
  !Q is the block matrix holder
  implicit none
  
  integer :: m,i,q,j
  real(8),dimension(m,m) :: A
  real(8),dimension(m) :: eig
  type(full_sp_block_mat) :: R
  
  A=0.d0
  i=1
  do q=1,R%blocks
     
     A(i:i+R%map(q)-1,i:i+R%map(q)-1) = R%blkM(q)%matrix
     eig(i:i+R%map(q)-1) = R%blkM(q)%eigval
     i=i+R%map(q)
  end do 
end subroutine 
!======================================================
subroutine diagonalize_blocks(R)
  implicit none 
  
  type(full_sp_block_mat) :: R
  integer :: a,q,info
  
  do q=1,R%blocks
     
     a=R%map(q)
     if (a == 0 ) cycle
     call dsyev('V','U',a,R%blkM(q)%matrix,a, &
          R%blkM(q)%eigval,R%blkM(q)%extra,10*a,info)
   
  end do 
end subroutine
!====================================================
subroutine add_arrays(r1,s1,r2,s2,r3)
  implicit none 
	
	! r3 = s1 * r1 + s2 * r2
  
  type(full_ham) :: r1,r2,r3
  integer :: i,j,k
  real(8) :: s1,s2

  
  r3%E0 = s1 * r1%E0 + s2 * r2%E0
  r3%fph = s1 * r1%fph + s2 * r2%fph
  r3%fpp = s1 * r1%fpp + s2 * r2%fpp
  r3%fhh = s1 * r1%fhh + s2 * r2%fhh
  
  do i=1,r1%nblock 
     
     r3%mat(i)%Vhhhh= s1 * r1%mat(i)%Vhhhh + s2 * r2%mat(i)%Vhhhh
     r3%mat(i)%Vpppp= s1 * r1%mat(i)%Vpppp + s2 * r2%mat(i)%Vpppp
     r3%mat(i)%Vpphh= s1 * r1%mat(i)%Vpphh + s2 * r2%mat(i)%Vpphh
     r3%mat(i)%Vphhh= s1 * r1%mat(i)%Vphhh + s2 * r2%mat(i)%Vphhh
     r3%mat(i)%Vppph= s1 * r1%mat(i)%Vppph + s2 * r2%mat(i)%Vppph
     r3%mat(i)%Vphph= s1 * r1%mat(i)%Vphph + s2 * r2%mat(i)%Vphph

  end do 

end subroutine 
!=====================================================
subroutine copy_arrays(r1,r2)
	! copy r1 onto r2
  implicit none 
  
  type(full_ham) :: r1,r2
  integer :: i
  
  r2%E0 = r1%E0
  r2%fph=r1%fph
  r2%fpp=r1%fpp
  r2%fhh=r1%fhh
  
  do i=1,r1%nblock 
             
     if ( r1%mat(i)%nhh> 0 ) then 
     r2%mat(i)%Vhhhh=r1%mat(i)%Vhhhh
     end if 
     
     if ( r1%mat(i)%npp > 0 ) then 
     r2%mat(i)%Vpppp=r1%mat(i)%Vpppp
     end if 
     
     if ( r1%mat(i)%npp *r1%mat(i)%nhh > 0 ) then 
     r2%mat(i)%Vpphh=r1%mat(i)%Vpphh
     end if 
     
     if ( r1%mat(i)%nph *r1%mat(i)%nhh > 0 ) then 
     r2%mat(i)%Vphhh=r1%mat(i)%Vphhh
     end if 

     if ( r1%mat(i)%npp *r1%mat(i)%nph > 0 ) then
     r2%mat(i)%Vppph=r1%mat(i)%Vppph
     end if 
     
     if ( r1%mat(i)%nph *r1%mat(i)%nph > 0 ) then 
     r2%mat(i)%Vphph=r1%mat(i)%Vphph
     end if 

  end do 

end subroutine 
!=============================================
subroutine set_to_zero(r1)
  implicit none
  
  type(full_ham) :: r1
  integer :: i,nh,np,nb
  
  r1%E0 = 0.d0
  r1%fhh =0.d0
  r1%fpp = 0.d0
  r1%fph = 0.d0
  
  do i=1,r1%nblock
     nh=r1%mat(i)%nhh
     np=r1%mat(i)%npp
     nb=r1%mat(i)%nph
     
     if (nh > 0 ) r1%mat(i)%Vhhhh = 0.d0
     if (np > 0 ) r1%mat(i)%Vpppp = 0.d0
     if (nb > 0 ) r1%mat(i)%Vphph = 0.d0
     if (nh*np > 0 ) r1%mat(i)%Vpphh = 0.d0
     if (nh*nb > 0 ) r1%mat(i)%Vphhh = 0.d0
     if (nb*np > 0 ) r1%mat(i)%Vppph = 0.d0
     
  end do 
end subroutine
!=============================================
subroutine arrange_states(maxE,quant_num,order)
  implicit none 
  !!! states are ordered by these criteria:
 
  !!! 1. lowest spin will have lowest ordering
  !!! 2. lowest abs(m) will have lowest ordering 
  !!! 3. negative m is lower than positive
  !!! 4. for same m, lowest n will have lowest ordering
  
  !!! this way we have a block diagonal hamiltonian
  !!! so LAPACK won't have a stroke... 

  integer,dimension((maxE*maxE+maxE),3) :: quant_num 
  integer,dimension((maxE*maxE+maxE)/2) :: order,en
  integer :: i,maxE,a,b,j,m,n
  
  order=0
  j=1
  quant_num(:,3)=-1
  
  a = maxE*(maxE+1)/2
  
  do m=0,maxE-1
     
    do n=0, (maxE-1-m)/2
          
          quant_num( j : j+1 , 1 ) = n
          quant_num( j , 2 ) = -m
          en(j)= 2*n+m+1
        
          j=j+1
    end do 
    
    if (m==0) cycle
    
    do n=0, (maxE-1-m)/2
      
          quant_num( j : j+1 , 1 ) = n
          quant_num( j , 2 ) = m
          en(j)= 2*n+m+1
        
          j=j+1
    end do 
    
 end do 
 
 quant_num(a+1:2*a,:) = quant_num( 1:a, : )
 
 quant_num(a+1:2*a,3) = 1
 
j=1
 

!! vector that stores the order states are filled is compiled
do while ( j .le. maxE*(maxE+1)/2 )
   i=1
 do while ( i .le. maxE*(maxE+1)/2 ) 
    
  
    if ( en( i ) == minval(en) ) then 
       
          
          order(j) = i 
          i=i+1
          j=j+1
          
          en( i-1 ) = 1000

          exit
      
       
   end if
   
   i=i+1

  end do 
end do 

end subroutine 
!===============================================
subroutine get_ints(ints,hw,qn,m) 
  implicit none 
  
  real(8) :: hw,v_int
  integer :: ii,jj,a,b,q,n,m,nh,np,nb,pre
  integer :: h1,h2,h3,h4,p1,p2,p3,p4
  integer,dimension(m,3) :: qn
  type(full_ham) :: ints
  
  n = ints%nbody
  m = ints%msp
  
  do q = 1,ints%nblock
     
     np = ints%mat(q)%npp

     allocate( ints%mat(q)%Vpppp( np , np ) )
  
! Vpppp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ints%mat(q)%Vpppp = 0.d0
     
     do ii=1,np
        do jj=ii+1,np
           
           p1=ints%mat(q)%qnpp(II,1)
           p2=ints%mat(q)%qnpp(II,2)
           p3=ints%mat(q)%qnpp(JJ,1)
           p4=ints%mat(q)%qnpp(JJ,2)
           
      ints%mat(q)%Vpppp(II,JJ)  = v_int( qn(p1,1:3) , qn(p2,1:3),   &
                  qn(p4,1:3) , qn(p3,1:3) , hw ) - &
                     v_int( qn(p1,1:3) , qn(p2,1:3),   &
                  qn(p3,1:3) , qn(p4,1:3) , hw )       
            
        end do 
     end do 
     
     ints%mat(q)%Vpppp = ints%mat(q)%Vpppp + Transpose( ints%mat(q)%Vpppp ) 
     
     do ii=1,np
        p1=ints%mat(q)%qnpp(II,1)
        p2=ints%mat(q)%qnpp(II,2)
   
        ints%mat(q)%Vpppp(II,II) = v_int( qn(p1,1:3) , qn(p2,1:3),   &
                  qn(p2,1:3) , qn(p1,1:3) , hw ) - &
                  v_int( qn(p1,1:3) , qn(p2,1:3),   &
                  qn(p1,1:3) , qn(p2,1:3) , hw )
     end do 
  
  end do 

end subroutine
!======================================================      
subroutine calc_h0(h0,hw,qn)
  ! assuming a harmonic oscillator trap
  implicit none 
  
  integer :: i,m,sz
  integer,dimension(:,:) :: qn 
  real(8),dimension(:,:) :: h0
  real(8) :: hw
  

  sz=size(qn(:,1))

  h0=0.d0
  do i=1,sz
     
     h0(i,i)=2*qn(i,1)+abs(qn(i,2))+1

  end do 
  
  h0=h0*hw

end subroutine
!===============================================    
subroutine print_matrix(matrix)
	implicit none 
	
	integer :: i,m
	real(8),dimension(:,:) :: matrix
	character(1) :: y
	character(10) :: fmt2

    m=size(matrix(1,:))

    write(y,'(i1)') m
  
    fmt2= '('//y//'(f14.8))'	
	
	print*
	do i=1,m
	   write(*,fmt2) matrix(i,:)
	end do
	print* 
	
end subroutine 
!===============================================    
subroutine write_matrix(matrix,fname)
	implicit none 
	
	integer :: i,m
	real(8),dimension(:,:) :: matrix
	character(1) :: y
	character(10) :: fmt2
        character(*) :: fname

    m=size(matrix(1,:))

    write(y,'(i1)') m
  
    fmt2= '('//y//'(f14.8))'	
 
    open(unit=36,file=fname)
	do i=1,m
	   write(36,fmt2) matrix(i,:)
	end do
    close(36)
	
end subroutine 
!==========================================
subroutine calcDen(rho,c,n,m,qn,ord)
	! density matrix
        ! n is no. of particles
        ! m is size of basis
        ! c is coefficient matrix
  implicit none 
  
  integer :: m,n,i,j,k
  real(8), parameter :: al=1.d0, bet=0.d0
  real(8), dimension(m,m) :: rho,c
  real(8), dimension(m,n) :: imt
  integer, dimension(m,3) :: qn
  integer, dimension(m/2) :: ord
  
  do i=1,n/2
     
     imt(:,i) = c(:, ord(i) )
     imt(:,i+n/2) =  c(: , ord(i)+m/2) 
     
  end do 
  
  call dgemm('N','T',m,m,n,al,imt,m,imt,m,bet,rho,m) 
  
end subroutine
!============================================ 
subroutine reorder(c,w,qn,m)
  implicit none 
 
  integer :: m,i,j,k,q,p,r
  real(8),dimension(m,m) :: c,c2,c3
  real(8),dimension(m) :: w,vv,w2
  integer,dimension(m,3) :: qn
  
  c2=c
  c3=abs(c)
  
  w2=w
  
  do i=1,m
     do j=1,m
        
        if ( c3(i,j) == maxval(c3(:,j)) ) then 
           
           c(:,i)=c2(:,j) 
           w(i)=w2(j)
           exit
        end if
     
     end do 
 end do 
         
end subroutine  
!==============================================      
subroutine reindex(m,g,n,gst,pst,qn,emax,eig,rec)
  !!! rewrites states into a tp index
  !!! does most of the production of the inital 
  !!! hamiltonian storage
  implicit none 
  
  integer :: a,b,c,m,g,n,i,j,k,qq,wq,II
  integer :: ML,MS,gst,pst,emax,p,l,s(3)
  integer,dimension(m,3) :: qn
  integer,dimension(g,2) :: qntp
  integer,dimension(3,4) :: ijs
  real(8) :: eF, eig(m),ecop(m),ec2(m)
  type(full_ham) :: rec 
  
  rec%nblock = 12 * emax - 9
  rec%nbody = n
  rec%herm = 1
  allocate( rec%i_array( g , 3 ) )
  allocate( rec%mat( rec%nblock ) ) 
  
  !!! this loop determines the fermi energy
  ecop = eig
  ec2 = eig
  i=1
  j=0
  if (n > 0) then 
  do 	
     if (ec2(i) == minval(ecop)) then 
	   
         ecop(i) = 1000.d0
	 ec2(i) = 1000.d0
	 j=j+1
  
	 if (j == n) then 
	     eF = eig(i) 
	   
             exit
	 end if
  
	 i=1
     else 
	 i = i + 1
     end if 
		
  end do 

  ! adding a bit to the fermi energy to make comparisons easier
  eF = eF + 1.d-6
  else
  eF = 0.d0
  end if 
  
  k=1
  l=1
  p=1
  
  a=1;b=1;c=1

  qq=1
  wq=1

  ijs(1,:)=(/ 1 , m/2 , m/2+1 , m /)
  ijs(2,:)=(/ 1 , m/2 , 1 , m/2 /)
  ijs(3,:)=(/ m/2+1 , m , m/2+1 , m /)
  s=(/ 0, -1, 1 /)
  
do MS=1,3
  !!! MS = 0,-2,2 
  
  do ML=0,2*emax-2
     
        do i = ijs(MS,1) , ijs(MS,2)

           do j= max(i+1,ijs(MS,3)) ,ijs(MS,4) 
              
               !! calculate general label for tp state
               II = (i-1) * M  - i*(i-1)/2 + j - i
               
             if   ( qn(i,2) + qn(j,2) == -ML ) then 
                
                
                if (eig(i) .le. eF ) then 
                   
                   if ( eig(j) .le. eF ) then 
                      
                      !print*, 'hh',i,j,-ML,s(MS),k
                      qntp(k,1)=i
                      qntp(k,2)=j
                      rec%i_array(II,1)=1
                      rec%i_array(II,2)=qq
                      k=k+1
                      rec%i_array(II,3)=k-a
                      wq=wq+1

                   else 
                      
                      !print*, 'ph',i,j,-ML,s(MS),l
                      qntp(gst+l,1)=i
                      qntp(gst+l,2)=j
                      rec%i_array(II,1)=2
                      rec%i_array(II,2)=qq
                      l=l+1
                      rec%i_array(II,3)=l-b
                      wq=wq+1

                   end if 

               else 
                  
                  if ( eig(j) .le. eF ) then 
                     
                      !print*, 'ph',i,j,-ML,s(MS),l
                      qntp(gst+l,1)=i
                      qntp(gst+l,2)=j
                      rec%i_array(II,1)=2
                      rec%i_array(II,2)=qq
                      l=l+1
                      rec%i_array(II,3)=l-b
                      wq=wq+1
                      
                  else 
                     
                      !print*, 'pp',i,j,-ML,s(MS),p
                      qntp(g-pst+p,1)=i
                      qntp(g-pst+p,2)=j
                      rec%i_array(II,1)=3
                      rec%i_array(II,2)=qq
                      p=p+1
                      rec%i_array(II,3)=p-c
                      wq=wq+1
                      
                  end if 
              end if 

          end if
 
       end do 
      end do

      
      !!! build the block in the full array.
     rec%mat(qq)%lam(1)=-ML 
     rec%mat(qq)%lam(2)=s(MS)
     
     allocate(rec%mat(qq)%qnhh(k-a,2))
     rec%mat(qq)%nhh=k-a
     allocate(rec%mat(qq)%qnph(l-b,2))
     rec%mat(qq)%nph=l-b
     allocate(rec%mat(qq)%qnpp(p-c,2))
     rec%mat(qq)%npp=p-c
     
     rec%mat(qq)%ntot= k+l+p-a-b-c

     rec%mat(qq)%qnhh=qntp(a:k-1,1:2)
     rec%mat(qq)%qnph=qntp(gst+b:gst+l-1,1:2) 
     rec%mat(qq)%qnpp=qntp(g-pst+c:g-pst+p-1,1:2)

     a=k
     b=l
     c=p
     qq=qq+1


              if ( ML==0 )  cycle 
              
      do i = ijs(MS,1) , ijs(MS,2)

       do j= max(i+1,ijs(MS,3)) ,ijs(MS,4)
          
          !! calculate general label for tp state
               II = (i-1) * M  - i*(i-1)/2 + j - i
              
          if   ( qn(i,2) + qn(j,2) == ML ) then 
               
                
                if ( eig(i) .le. eF ) then 
                   
                   if ( eig(j) .le. eF ) then 
                      !print*, 'hh',i,j,ML,s(MS),k
                      qntp(k,1)=i
                      qntp(k,2)=j
                      rec%i_array(II,1)=1
                      rec%i_array(II,2)=qq
                      k=k+1
                      rec%i_array(II,3)=k-a
                      wq=wq+1

                   else 
                      !print*, 'ph',i,j,ML,s(MS),l
                      qntp(gst+l,1)=i
                      qntp(gst+l,2)=j
                      rec%i_array(II,1)=2
                      rec%i_array(II,2)=qq
                      l=l+1
                      rec%i_array(II,3)=l-b
                      wq=wq+1

                   end if 

               else 
                  
                  if ( eig(j) .le. eF ) then 
                     
                      !print*, 'ph',i,j,ML,s(MS),l
                      qntp(gst+l,1)=i
                      qntp(gst+l,2)=j
                      rec%i_array(II,1)=2
                      rec%i_array(II,2)=qq
                      l=l+1
                      rec%i_array(II,3)=l-b
                      wq=wq+1
                 
                  else 
                     ! print*, 'pp',i,j,ML,s(MS),p
                      qntp(g-pst+p,1)=i
                      qntp(g-pst+p,2)=j
                      rec%i_array(II,1)=3
                      rec%i_array(II,2)=qq
                      p=p+1
                      rec%i_array(II,3)=p-c
                      wq=wq+1
                      
                  end if 
              end if 

          end if      
          
       end do 
     end do 
    
     !!! build the block array in the full array
     rec%mat(qq)%lam(1)=ML 
     rec%mat(qq)%lam(2)=s(MS)

     allocate(rec%mat(qq)%qnhh(k-a,2))
     rec%mat(qq)%nhh=k-a
     allocate(rec%mat(qq)%qnph(l-b,2))
     rec%mat(qq)%nph=l-b
     allocate(rec%mat(qq)%qnpp(p-c,2))
     rec%mat(qq)%npp=p-c
     
     rec%mat(qq)%ntot= k+l+p-a-b-c
     
     rec%mat(qq)%qnhh=qntp(a:k-1,1:2)
     rec%mat(qq)%qnph=qntp(gst+b:gst+l-1,1:2) 
     rec%mat(qq)%qnpp=qntp(g-pst+c:g-pst+p-1,1:2)
     
     a=k
     b=l
     c=p
     qq=qq+1

  end do 
           
end do 

end subroutine 
!==============================================
subroutine new_COEFS(c,w,m,rec)
  implicit none 
  
  integer :: m,i,j,q,nh,nt
  real(8),dimension(m,m) :: c
  integer,allocatable,dimension(:,:) :: qntp
  real(8) :: w(m)
  type(full_ham) :: rec

  do q=1,rec%nblock
     
     !! some things redefined to lower amount of % 
     nt=rec%mat(q)%ntot
     allocate(qntp( nt , 2 ))
     
     !!! allocate new coefficent matrices
     allocate(rec%mat(q)%CC(nt,nt))
     allocate(rec%mat(q)%CB(nt,nt))

     nh=rec%mat(q)%nhh+rec%mat(q)%nph
     
     qntp( 1 : rec%mat(q)%nhh , : ) = rec%mat(q)%qnhh
     qntp( rec%mat(q)%nhh+1 : nh , : ) = rec%mat(q)%qnph
     qntp( nh+1 : nt , : ) = rec%mat(q)%qnpp

  !!! new eigenvector storage
  do i=1,nt
     do j=1,nt
  
        rec%mat(q)%CC(i,j) = c ( qntp(j,1) , qntp(i,1) ) * &
                     c ( qntp(j,2) , qntp(i,2) ) 
        
        rec%mat(q)%CB(i,j) =  c ( qntp(j,2) , qntp(i,1) ) * &
                     c ( qntp(j,1) , qntp(i,2) ) 
     end do 
  end do
  
  
  !!! new eigenvalue storage
  allocate(rec%mat(q)%eig2(nt))

  do i=1,nt 
     rec%mat(q)%eig2(i) = w( qntp(i,1) ) + w( qntp(i,2) ) 
  end do 
  
  deallocate(qntp)

  end do 
 
end subroutine
!==============================================
subroutine new_VMAT(ints,rec,m) 
  ! transform matrix elements into HF basis
  implicit none 
  
  integer :: i,j,m,nt,nh,np,nb,q,a1,a2
  integer :: px,qx,rx,sx,II,JJ
  integer,allocatable,dimension(:,:) :: qntp
  real(8),allocatable,dimension(:,:) :: NEWV, HFV
  type(full_ham) :: rec,ints

  do q=1,rec%nblock
     
     !! some things redefined to lower amount of % 
     nt=rec%mat(q)%ntot
      
     if (nt == 0 ) then 
        allocate(rec%mat(q)%vhhhh(0,0))
        allocate(rec%mat(q)%vpppp(0,0))
        allocate(rec%mat(q)%vpphh(0,0))
        allocate(rec%mat(q)%vphhh(0,0))
        allocate(rec%mat(q)%vppph(0,0))
        allocate(rec%mat(q)%vphph(0,0))
        cycle
     end if 
     
     allocate(qntp( nt , 2 ))
     
     !!! allocate new coefficent matrices
     allocate(NEWV(nt,nt))
     allocate(HFV(nt,nt))

     nh=rec%mat(q)%nhh+rec%mat(q)%nph
     
     qntp( 1 : rec%mat(q)%nhh , : ) = rec%mat(q)%qnhh
     qntp( rec%mat(q)%nhh+1 : nh , : ) = rec%mat(q)%qnph
     qntp( nh+1 : nt , : ) = rec%mat(q)%qnpp

  NEWV=0.d0
  
  do i=1,nt
     do j=1,nt
                
        
          px = qntp(i,1) ; qx = qntp(i,2) 
          rx = qntp(j,1) ; sx = qntp(j,2) 
          II = (px-1) * M  - px*(px-1)/2 + qx - px
          JJ = (rx-1) * M  - rx*(rx-1)/2 + sx - rx

           a1 = ints%i_array(II,3)
           a2 = ints%i_array(JJ,3)
           
           NEWV(i,j) = ints%mat(q)%Vpppp(a1,a2) 

     end do 
  end do 
  
  call HF_Vmatrix(HFV , NEWV , rec%mat(q)%CC , rec%mat(q)%CB , nt)
  
  !!! so now we should have a matrix HFV full of all of the hf matrix elements of V
  !!! now we must split them up into their sub-categories

  allocate( rec%mat(q)%Vpppp( rec%mat(q)%npp , rec%mat(q)%npp ) )
  
  do i=1,rec%mat(q)%npp
     do j=1,rec%mat(q)%npp
        
        rec%mat(q)%Vpppp(i,j) = HFV( nh + i , nh + j )  
     
     end do 
  end do 
  
  allocate( rec%mat(q)%Vhhhh( rec%mat(q)%nhh , rec%mat(q)%nhh ) )
  
  do i=1,rec%mat(q)%nhh
     do j=1,rec%mat(q)%nhh
        
        rec%mat(q)%Vhhhh(i,j) = HFV(  i , j )  
     
     end do 
  end do 
  
  allocate( rec%mat(q)%Vppph( rec%mat(q)%npp , rec%mat(q)%nph ) )
  
  do i=1,rec%mat(q)%npp
     do j=1,rec%mat(q)%nph
        
        rec%mat(q)%Vppph(i,j) = HFV( nh + i , rec%mat(q)%nhh + j )  
     
     end do 
  end do 
  
  allocate( rec%mat(q)%Vphhh( rec%mat(q)%nph , rec%mat(q)%nhh ) )
  
  do i=1,rec%mat(q)%nph
     do j=1,rec%mat(q)%nhh
        
        rec%mat(q)%Vphhh(i,j) = HFV( rec%mat(q)%nhh + i ,  j )  
     
     end do 
  end do 
  
  allocate( rec%mat(q)%Vpphh( rec%mat(q)%npp , rec%mat(q)%nhh ) )
  allocate( rec%mat(q)%Dpphh( rec%mat(q)%npp , rec%mat(q)%nhh ) )
  
  do i=1,rec%mat(q)%npp
     do j=1,rec%mat(q)%nhh
        
        rec%mat(q)%Vpphh(i,j) = HFV( nh + i ,  j )  
        rec%mat(q)%Dpphh(i,j) = HFV( nh + i , j ) / &
             ( rec%mat(q)%eig2(j) - rec%mat(q)%eig2( nh + i) ) 
        
     end do 
  end do 
  
  allocate( rec%mat(q)%Vphph( rec%mat(q)%nph , rec%mat(q)%nph ) )
  
  do i=1,rec%mat(q)%nph
     do j=1,rec%mat(q)%nph
        
        rec%mat(q)%Vphph(i,j) = HFV( rec%mat(q)%nhh + i , rec%mat(q)%nhh + j )  
     
     end do 
  end do 
  
  deallocate(qntp)
  deallocate(NEWV)
  deallocate(HFV)

  end do 
  
end subroutine 
!==============================================
!==============================================
subroutine convertF(F,rec,m,qn,w)
  !!!! fill up f matrices and define energy ordering
  implicit none 

  integer :: m,i,j,q,n
  integer,dimension(m,3) :: qn
  real(8),dimension(m,m) :: F
  real(8),dimension(m) :: w
  real(8) :: mn
  type(full_ham) :: rec
  
  n=rec%nbody
  rec%Msp=m
  allocate( rec%states(m,3) ) 
  allocate( rec%eig(m) ) 
  allocate( rec%eh( n ) ) 
  allocate( rec%ep( m-n) )
  allocate( rec%stoe(m) )
  
  rec%states=qn
  rec%eig=w
  
  do i=1,n
     
     mn = minval(w)
     do j=1,m
        
        if ( w(j) == mn ) then 
           rec%eh(i)=j
           rec%stoe(j)=i
           w(j)=1000. 
           exit
        end if 
     end do 
     
  end do 

  do i=n+1,m
     
     mn = minval(w)
     do j=1,m
        
        if ( w(j) == mn ) then 
           rec%ep(i-n)=j
           rec%stoe(j)=i
           w(j)=1000. 
           exit
        end if 
     end do 
     
  end do 

  
  !!! redistribution of f_ij into pp,ph,hh
  allocate(rec%fhh(n,n))
  allocate(rec%fph(m-n,n)) 
  allocate(rec%fpp(m-n,m-n)) 
  
  do i=1,n
     do j=1,n
        
        rec%fhh(i,j) = F(rec%eh(i),rec%eh(j))
      
     end do 
     
     do j=1,m-n
        
        rec%fph(j,i) = F(rec%ep(j) , rec%eh(i) ) 
     end do 

  end do 
  
  do i=1,m-n
     do j=1,m-n
        
        rec%fpp(i,j) = F(rec%ep(i) ,rec%ep(j) ) 
     
     end do 
  end do 

end subroutine 
!==============================================
real(8) function v_elem(a1,a2,a3,a4,rec)
  !!! a,b,c,d should be in ENERGY ORDER
  implicit none 
  
  integer :: a(4),a1,a2,a3,a4,b,ML,MS,p,II,JJ
  integer :: k(4),i,q,w1,n,w2,pre,xx,m,g,p1,p2,row(2,3)
  type(full_ham) :: rec
  logical :: flg
  
  if (( a1 == a2) .or. (a3 == a4) ) then 
     v_elem = 0.d0
  else
  a(1)=a1
  a(2)=a2
  a(3)=a3
  a(4)=a4
  n=rec%nbody
  m = rec%msp
  pre=1

  !!! convert to STATES ordering
  do b=1,4
     if ( a(b) .le. n ) then 
        k(b)=rec%eh(a(b))
     else 
        k(b)=rec%ep(a(b)-n)
     end if 
  end do 

  !!! switch indeces if they are wrong
  if ( k(1) > k(2) ) then 
     xx= k(1)
     k(1)=k(2)
     k(2)=xx
     pre=-1
  end if 
  
  if ( k(3) > k(4) )  then 
     xx = k(3)
     k(3)=k(4) 
     k(4) = xx
     pre=-1*pre
  end if 
        
        II = (k(1)-1) * M  - k(1)*(k(1)-1)/2 + k(2) - k(1)
        JJ = (k(3)-1) * M  - k(3)*(k(3)-1)/2 + k(4) - k(3)
        
        w2= rec%i_array(JJ,3)   !!! qn index
        p2= rec%i_array(JJ,2)   !!! block index   
        w1 = rec%i_array(II,3)
        p1 = rec%i_array(II,2)
      
        !!!! decide if the matrix element exists
  if ( p1 .ne. p2 ) then 
     
     v_elem = 0.d0 
         
  else  
           
  !!! we have to decide which arrays to use, so 
  !!! we have this ridiculous mess:
     
  select case (rec%i_array(II,1))
     
     case (1) 
        
        select case (rec%i_array(JJ,1))
           
           case (1) 
              
              v_elem=rec%mat(p1)%Vhhhh(w1,w2)*pre
              
           case (2) 
              
              v_elem=rec%mat(p1)%Vphhh(w2,w1)*pre*rec%herm
              
           case (3) 
              
              v_elem=rec%mat(p1)%Vpphh(w2,w1)*pre*rec%herm
         
        end select
        
     case (2) 
        
        select case (rec%i_array(JJ,1))
           
           case (1) 
              
              v_elem=rec%mat(p1)%Vphhh(w1,w2)*pre
              
           case (2) 
              
              v_elem=rec%mat(p1)%Vphph(w1,w2)*pre
              
           case (3) 
              
              v_elem=rec%mat(p1)%Vppph(w2,w1)*pre*rec%herm
         
        end select
        
    case (3) 
        
        select case (rec%i_array(JJ,1))
           
           case (1) 
              
              v_elem=rec%mat(p1)%Vpphh(w1,w2)*pre
              
           case (2) 
              
              v_elem=rec%mat(p1)%Vppph(w1,w2)*pre
              
           case (3) 
              
              v_elem=rec%mat(p1)%Vpppp(w1,w2)*pre
                   
        end select
       
     end select

  end if 
  
  end if
end function
!==============================================
real(8) function w_elem(px,qx,rx,sx,tx,ux,rec) 
  implicit none
  
  integer :: p,q,r,s,t,u,II,JJ,m,x,pre
  integer :: px,qx,rx,sx,tx,ux
  type(full_ham) :: rec
  p = px
  q = qx
  r = rx
  s = sx
  t = tx
  u = ux
  
  if ((p == q) .or. ( p == r) .or. (q == r)&
       .or. (s == t) .or. (s == u) .or. (t == u) ) then 
     w_elem = 0.d0
  else
     
  m = rec%msp
  
  pre = 1
  if (p > q)  then 
     x = p
     p = q
     q = x
     pre = -1 
  end if 
  
  if (p > r)  then 
     x = p
     p = r
     r = x
     pre = -1*pre
  end if 

  if ( q > r)  then 
     x = q
     q = r
     r = x
     pre = -1*pre
  end if 
  
   if (s > t)  then 
     x = s
     s = t
     t = x
     pre = -1*pre
  end if 
  
  if (s > u)  then 
     x = s
     s = u
     u = x
     pre = -1*pre
  end if 
  
  if ( t > u)  then 
     x = t
     t = u
     u = x
     pre = -1*pre
  end if 
        
  II = m*(m-1)*(m-2)/6 + (m-p)*(m-p-1)/2 - (m+1-p) *(m-p) * (m-1-p) /6 - &
       (m+1-q)*(m-q)/2 + r - q
  
  JJ = m*(m-1)*(m-2)/6 + (m-s)*(m-s-1)/2 - (m+1-s) *(m-s) * (m-1-s) /6 - &
       (m+1-t)*(m-t)/2 + u - t
  
  w_elem = rec%v3body(II,JJ)*pre
  end if 
end function
!========================================
integer function get_num_elements(rec)
  implicit none 
  
  type(full_ham) :: rec
  integer :: i,j,k,l,n,m,nb,nh,np
   
 n = rec%nbody
 m = rec%Msp - n
 
 l = 1 + n*n + m*m + n*m
 
 do i=1,rec%nblock
    
    nb=rec%mat(i)%nph
    np=rec%mat(i)%npp
    nh=rec%mat(i)%nhh
    
    l = l + nh*nh + np*np + nh*np + &
         nb*nb + nb*nh + nb*np
  end do 
  
  get_num_elements = l
  
end function
!==============================================  
real(8) function f_elem(a,b,rec) 
  implicit none 
  
  integer :: a,b,n
  type(full_ham) :: rec
  logical :: f1,f2

  n=rec%nbody
  if ( a .le. n ) then 
     if ( b .le. n ) then 
        
        f_elem = rec%fhh(a,b) 
       
     else 
        
        f_elem = rec%fph(b-n,a)*rec%herm
       
     end if 

   else 

     if ( b .le. n ) then 
        
        f_elem = rec%fph(a-n,b) 
       
     else 
        
        f_elem = rec%fpp(a-n,b-n)
       
     end if 
     
   end if 

end function 
!==============================================
!=================================================
!=================================================
subroutine vectorize(rec,vout,neq)
  !!! maps full_ham to vector
  implicit none 

  integer ::  neq,n,m,m3,i,j,k,l,nb,np,nh
  real(8),dimension(neq) :: vout
  type(full_ham) :: rec

  n= rec%nbody
  m= rec%Msp
  m3 = m*(m-1)*(m-2)/6
  
  m = m-n
  
  vout(1) = rec%E0
  
  k=2
  do i=1,n
     do j=1,n
        
        vout(k) =  rec%fhh(j,i) 
        k=k+1
       
     end do 
  end do 

  do i=1,m
     do j=1,m
        
        vout(k) = rec%fpp(j,i) 
        k=k+1
        

     end do 
  end do 

  do i=1,n
     do j=1,m
        
        vout(k) = rec%fph(j,i) 
        k=k+1

     end do 
  end do 
  
  do l=1,rec%nblock
     
     nh=rec%mat(l)%nhh
     np=rec%mat(l)%npp
     nb=rec%mat(l)%nph
     
     do i=1,nh
        do j=1,nh
           
           vout(k) = rec%mat(l)%Vhhhh(j,i) 
           k=k+1
        
        end do 
     end do 
     
     do i=1,np
        do j=1,np
           
           vout(k) = rec%mat(l)%Vpppp(j,i)
           k=k+1
           
        end do 
     end do 
     
     do i=1,nh
        do j=1,np
           
           vout(k) = rec%mat(l)%Vpphh(j,i)
           k=k+1
           
        end do 
     end do 
     
     do i=1,nb
        do j=1,nb
           
           vout(k) = rec%mat(l)%Vphph(j,i) 
           k=k+1

        end do 
     end do 
     
     do i=1,nb
        do j=1,np
           
           vout(k) = rec%mat(l)%Vppph(j,i)
           k=k+1

        end do 
     end do 
     
     do i=1,nh
        do j=1,nb
           
           vout(k)=rec%mat(l)%Vphhh(j,i) 
           k=k+1
           
        end do 
     end do 
    
  end do 

if (rec%IMSRG3) then  
! three body stuff ===================
  do i = 1, m3
     do j = 1,m3
        
        vout(k) = rec%V3body(i,j)
        k = k+1
     end do 
  end do 
!=====================================
end if 
end subroutine 
!===============================================================
!===============================================================
subroutine repackage(rec,vout,neq)
  !! maps vector to full_ham
  implicit none 

  integer ::  neq,n,m,m3,i,j,k,l,nb,np,nh
  real(8),dimension(neq) :: vout
  type(full_ham) :: rec

  n= rec%nbody
  m= rec%Msp 
  m3 = m*(m-1)*(m-2)/6
  m = m-n


  rec%E0 = vout(1)
  
  k=2
  do i=1,n
     do j=1,n
        
        rec%fhh(j,i) = vout(k) 
        k=k+1

     end do 
  end do 

  do i=1,m
     do j=1,m
        
        rec%fpp(j,i) = vout(k) 
        k=k+1

     end do 
  end do 

  do i=1,n
     do j=1,m
        
        rec%fph(j,i) = vout(k) 
        k=k+1

     end do 
  end do 
  
  do l=1,rec%nblock
     
     nh=rec%mat(l)%nhh
     np=rec%mat(l)%npp
     nb=rec%mat(l)%nph
     
     do i=1,nh
        do j=1,nh
           
           rec%mat(l)%Vhhhh(j,i) = vout(k) 
           k=k+1
        
        end do 
     end do 
     
     do i=1,np
        do j=1,np
           
           rec%mat(l)%Vpppp(j,i) = vout(k)
           k=k+1
           
        end do 
     end do 
     
     do i=1,nh
        do j=1,np
           
           rec%mat(l)%Vpphh(j,i) = vout(k)
           k=k+1
           
        end do 
     end do 
     
     do i=1,nb
        do j=1,nb
           
           rec%mat(l)%Vphph(j,i) = vout(k) 
           k=k+1

        end do 
     end do 
     
     do i=1,nb
        do j=1,np
           
           rec%mat(l)%Vppph(j,i) = vout(k) 
           k=k+1

        end do 
     end do 
     
     do i=1,nh
        do j=1,nb
           
           rec%mat(l)%Vphhh(j,i) = vout(k)
           k=k+1
           
        end do 
     end do 
    
  end do 

if (rec%IMSRG3) then 
! three body stuff =================== 
  do i = 1, m3
     do j = 1,m3
        
        rec%V3body(i,j) = vout(k)
        k = k+1
     end do 
  end do 
!=====================================
end if 
end subroutine 
!===============================================================
!===============================================================
real(8) function mat_2_norm( r ) 
  implicit none 
  
  type(full_ham) :: r
  integer :: n,m,q,nh,np,nb  
  real(8) :: od
  
  od = 0.d0
  m=r%msp
  n=r%nbody
  
  
  od = od + sum(r%fph**2)+sum(r%fpp**2)+sum(r%fhh**2)

  do q=1,r%nblock

 	
	od = od + sum(r%mat(q)%Vhhhh**2) + sum(r%mat(q)%Vpppp**2)    &
	+ sum(r%mat(q)%Vpphh**2) + sum(r%mat(q)%Vphph**2)    &
	+ sum(r%mat(q)%Vphhh**2) + sum(r%mat(q)%Vppph**2)    
	
     
 end do   
  
     mat_2_norm = sqrt(od)
  
 end function
!=====================================================
end module  
!=====================================================
integer function count_states(ml,ms,qn,m)
  implicit none 
  
  integer :: ml,ms,m,a,i 
  integer,dimension(m,3) :: qn
  
  a=0
  do i = 1,m
     if ( (qn(i,2) == ml) .and. (qn(i,3) == ms ) ) then 
         a=a+1
     end if 
  end do 
 
  count_states = a
end function 
!==============================================           
subroutine HF_Vmatrix(Vhf,V,MCOF,MCOFrev,g)
  implicit none
  
  !!! transforms from HO to HF basis

  integer :: g
  real(8),parameter :: al=1.d0,bet=0.d0
  real(8),dimension(g,g) :: Vhf,V,MCOF,MCOFrev,temp,V1,V2,V3,V4
  
  call dgemm('N','T',g,g,g,al,V,g,MCOF,g,bet,temp,g)
  call dgemm('N','N',g,g,g,al,MCOF,g,temp,g,bet,V1,g)
  
  call dgemm('N','T',g,g,g,al,V,g,MCOFrev,g,bet,temp,g)
  call dgemm('N','N',g,g,g,al,MCOFrev,g,temp,g,bet,V2,g)
  
  call dgemm('N','T',g,g,g,al,V,g,MCOFrev,g,bet,temp,g)
  call dgemm('N','N',g,g,g,al,MCOF,g,temp,g,bet,V3,g)

  call dgemm('N','T',g,g,g,al,V,g,MCOF,g,bet,temp,g)
  call dgemm('N','N',g,g,g,al,MCOFrev,g,temp,g,bet,V4,g)
  
  Vhf=V1+V2-V3-V4

end subroutine 
!==============================================              
!==============================================   
real(8) function v_int(qi,qj,qk,ql,hw)
  ! Anisimovas, Matulis. J. Pys.: Condens. Matter 10, 601 (1998)
  implicit none 
  
  integer :: kron_del,j1,j2,j3,j4,l1,l2,l3,l4
  integer,dimension(4) :: gs
  integer,dimension(3) :: qi,qj,qk,ql
  real(8) :: G,LAM,sm,gamma,sm_int,hw,vout
  real(16) :: fac_over_fac,bin_coef,fac
  
  vout = 0.d0
  sm=0.d0
  !!! just a ridiculous sum
  !if ( (kron_del( qi(2) , qk(2) ) * &
   !                kron_del( qj(2) , ql(2) )  ) == 1 ) then
     
  if ( (kron_del( qi(3) , ql(3) ) * &
                   kron_del( qk(3) , qj(3) )  ) == 1 ) then
     
     do j1=0,qi(1)
        do j2=0,qj(1)
           do j3=0,qk(1)
              do j4=0,ql(1)
                 
            gs(1)=j1+j4+(abs(qi(2))+qi(2))/2+(abs(ql(2))-ql(2))/2
            gs(4)=j1+j4+(abs(qi(2))-qi(2))/2+(abs(ql(2))+ql(2))/2
            gs(2)=j2+j3+(abs(qj(2))+qj(2))/2+(abs(qk(2))-qk(2))/2
            gs(3)=j2+j3+(abs(qj(2))-qj(2))/2+(abs(qk(2))+qk(2))/2
  
                 
            G=sum(gs)
                 
     sm_int=0.d0

     do l1=0,gs(1)
        do l2=0,gs(2)
           do l3=0,gs(3)
              do l4=0,gs(4)
                 
                 LAM=l1+l2+l3+l4

                sm_int=sm_int + kron_del(l1+l2,l3+l4) * &
  (-1)**(gs(2)+gs(3)-l2-l3) * bin_coef(gs(1),l1) * bin_coef(gs(2),l2) * &
 bin_coef(gs(3),l3) * bin_coef(gs(4),l4) * gamma(1.d0+ LAM * 0.5d0) * &
 gamma((G-LAM+1.d0) * 0.5d0)
              
              end do 
            end do 
         end do 
      end do 
                

           sm=sm+(-1)**( j1+j2+j3+j4 )/ (fac(j1)*fac(j2)*fac(j3)*fac(j4)) &
 *bin_coef(qi(1)+abs(qi(2)),qi(1)-j1) *  bin_coef(qj(1)+abs(qj(2)),qj(1)-j2) * &
 bin_coef(qk(1)+abs(qk(2)),qk(1)-j3)  * bin_coef(ql(1)+abs(ql(2)),ql(1)-j4) * &
0.5d0**((G+1.d0) * 0.5d0) * sm_int
           
             end do 
          end do 
       end do 
    end do 

 vout = sqrt(hw) * sm * sqrt( fac_over_fac( qi(1), qi(1)+abs(qi(2)) ) * &
 fac_over_fac( qj(1), qj(1)+abs(qj(2)) ) *  fac_over_fac( qk(1), qk(1)+abs(qk(2)) ) * &
 fac_over_fac( ql(1), ql(1)+abs(ql(2)) ) ) 

   ! end if 
 end if 

 v_int = vout
 
end function 
!===============================================
real(8) function psi(n,m,hw,x)
  implicit none 
  
  integer :: n,m,m2
  real(8) :: hw,x,Laguerre
  real(16) :: fac_over_fac
  real(8),parameter :: Pi_inv = 0.318309886184d0
  
  m2 = abs(m) 
 
  psi = hw * sqrt(fac_over_fac(n,n+m2)*Pi_inv)*(hw*x)**m2 &
       *exp(-(hw*x)**2*0.5)*Laguerre(n,m2,(hw*x)**2) 

end function   
!===============================================
real(8) function Laguerre(n,m,x)
  implicit none 
  
  integer :: n,m,i
  real(8) :: x,qx
  real(16) :: fac,bin_coef
  
  qx = 0.d0
  
  do i = 0, n
     
     qx = qx +(-1)**i * bin_coef(n+m,n-i) * x**i / fac(i) 
  
  end do 
  
  Laguerre = qx
  
end function  
!=============================================
! SPECIAL FUNCTIONS      
!=============================================  
real(8) function gamma(x)
  !!! only works for half integers or whole integers
  implicit none 

  real(8),parameter :: sqpi=1.77245385090551602d0
  real(8) :: x
  real(16) :: fac,fac_over_fac
  integer :: x_c
  
  x_c=floor(x)
  
  if (abs(float(x_c)-x) < 1d-3) then 
     
     gamma=fac(x_c-1)
     
  else 
     
      x_c=floor(x-0.5)
  
      gamma=fac_over_fac(2*x_c,x_c)/4.d0**x_c*sqpi  
  end if 
  
end function
!===========================================
real(kind=16) function bin_coef(n,k) 
  implicit none 
  
  integer :: n,k
  real(kind=16) :: fac,fac_over_fac
  
  bin_coef=fac_over_fac(n,k)/fac(n-k)
end function 
!==========================================
real(kind=16) function fac_over_fac(a1,a2)
  implicit none 
  
  integer :: a1,a2,i
  real(kind=16) :: s1
 
  s1=1.d0
  if (a1>a2) then 
  
  do i=a2+1,a1
     s1=s1*i
  end do 
  fac_over_fac=s1
  else 

  do i=a1+1,a2
     s1=s1*i
  end do 
  
  fac_over_fac=1.d0/s1
  
  end if
end function
!===========================================
real(kind=16) function fac(n)
  implicit none 
  
  integer :: n,i
  real(kind=16) :: s1
  
  s1=1.d0
  
  do i=1,n
     s1=s1*i
  end do 
  
  fac=s1
end function
!===========================================
integer function kron_del(i,j) 
  implicit none 
  
  integer :: i,j
  
  if (i==j) then 
     kron_del=1
  else 
     kron_del=0
  end if 
end function
!===============================================
