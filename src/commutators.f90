module commutators
  use ME_general
  implicit none
  
  real(8), parameter :: al =1.d0, bet=0.d0
  !!! THIS MODULE CONTAINS COMMUTATION ROUTINES. 

  !!! [R1, R2]  =  R3 or C  
 
contains

subroutine commutator_110(m,n,r1,r2,C)
  !0-body part of the 1body-1body commutator
  implicit none 
  
  integer :: m,i,n
  type(full_ham) :: r1,r2
  real(8),dimension(m,n) :: A,B
  real(8),dimension(m,m) :: Q
  real(8) :: c
  
  A=r1%fph
  B=r2%fph
  call dgemm('N','T',m,m,n,al,A,m,B,m,bet,Q,m)
  
  c=0.d0
  do i=1,m
     c=c+Q(i,i)
  end do 
  
  c=(r1%herm-r2%herm)*c
  
end subroutine
!================================================
!================================================  
subroutine commutator_220(r1,r2,c) 
  !0-body part of the 2body-2body commutator
  implicit none 
  
  integer :: m,n,i,j,gen_flag
  type(full_ham) :: r1,r2
  real(8),allocatable,dimension(:,:) :: Q 
  real(8) :: c
  
  c=0.d0 
     
     !! this part only contributes if we commute a hermitian
     !! with an anti-hermitian 
 
  do j=1,r2%nblock
     if (r1%mat(j)%npp*r1%mat(j)%nhh==0) cycle
     
     m=r1%mat(j)%npp
     n=r1%mat(j)%nhh
     allocate(Q(m,m)) 
     
     call dgemm('N','T',m,m,n,al,r1%mat(j)%Vpphh,m,&
          r2%mat(j)%Vpphh,m,bet,Q,m)
  
     do i=1,m
        c=c+Q(i,i)*(r1%herm-r2%herm)
     end do 
     
     deallocate(Q)
     
  end do  

  !!! we aren't double counting, so mult by -2 instead of -.5
  
end subroutine 
!===================================================================
!===================================================================
subroutine commutator_121(r1,r2,r3) 
  implicit none 
  
  type(full_ham) :: r1,r2,r3
  integer :: n,m,i,a,h1,h2,p1,p2

  n=r2%nbody
  m=r2%msp-r2%nbody

  ! hh
  do h2=1,n
     do h1=1,n

        do i=1,n
           do a=1,m
        
              r3%fhh(h1,h2) = r3%fhh(h1,h2) &
                   + r1%fph(a,i) * ( r1%herm * v_elem(a+n,h1,i,h2,r2) - v_elem(i,h1,a+n,h2,r2) ) &
                   - r2%fph(a,i) * ( r2%herm * v_elem(a+n,h1,i,h2,r1) - v_elem(i,h1,a+n,h2,r1) )
               
              
           end do
        end do
        
     end do

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
   !ph 
     do p1 = 1,m   

        do i=1,n
           do a=1,m
        
              r3%fph(p1,h2) = r3%fph(p1,h2) &
                   + r1%fph(a,i) * ( r1%herm * v_elem(a+n,p1+n,i,h2,r2) - v_elem(i,p1+n,a+n,h2,r2) ) &
                   - r2%fph(a,i) * ( r2%herm * v_elem(a+n,p1+n,i,h2,r1) - v_elem(i,p1+n,a+n,h2,r1) )
              
           end do
        end do
        
     end do 
  end do 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  do p2=1,m
     do p1=1,m

        do i=1,n
           do a=1,m
        
              r3%fpp(p1,p2) = r3%fpp(p1,p2) &
                   + r1%fph(a,i) * ( r1%herm * v_elem(a+n,p1+n,i,p2+n,r2) - v_elem(i,p1+n,a+n,p2+n,r2) ) &
                   - r2%fph(a,i) * ( r2%herm * v_elem(a+n,p1+n,i,p2+n,r1) - v_elem(i,p1+n,a+n,p2+n,r1) )
              
           end do
        end do

     end do 
  end do 

end subroutine commutator_121
!===================================================================
!===================================================================
subroutine commutator_111(r1,r2,r3)
  implicit none 
  
  type(full_ham) :: r1,r2,r3
  integer :: n,m,i,j,k
  real(8),allocatable,dimension(:,:) :: Q1,Q2,Q3,Q4
  
  n=r2%nbody
  m=r2%msp-r2%nbody
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!fhh
  allocate(Q1(n,n)) 
  allocate(Q2(n,n))

  call dgemm('N','N',n,n,n,al,r1%fhh,n,r2%fhh,n,bet,Q1,n)
  call dgemm('T','N',n,n,m,al,r2%fph,m,r1%fph,m,bet,Q2,n) 
  
  r3%fhh = Q1 - r2%herm *  Q2 + &
       r1%herm * ( Transpose(Q2) - r2%herm * Transpose(Q1) )    

  deallocate(Q1)
  deallocate(Q2)
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!fpp  
  allocate(Q1(m,m)) 
  allocate(Q2(m,m))

  call dgemm('N','N',m,m,m,al,r2%fpp,m,r1%fpp,m,bet,Q1,m)
  call dgemm('N','T',m,m,n,al,r1%fph,m,r2%fph,m,bet,Q2,m) 
  
  r3%fpp = r2%herm * Q2 + r1%herm *  &
       (r2%herm * Transpose(Q1) - Transpose(Q2)) - Q1

  deallocate(Q1)
  deallocate(Q2)
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  ! fph
  
  allocate(Q1(m,n))
  allocate(Q2(m,n))
  allocate(Q3(m,n))
  allocate(Q4(m,n))
  
  call dgemm('N','N',m,n,n,al,r1%fph,m,r2%fhh,n,bet,Q1,m)
  call dgemm('N','N',m,n,n,al,r2%fph,m,r1%fhh,n,bet,Q2,m)
  call dgemm('N','N',m,n,m,al,r1%fpp,m,r2%fph,m,bet,Q3,m)
  call dgemm('N','N',m,n,m,al,r2%fpp,m,r1%fph,m,bet,Q4,m) 

  r3%fph = Q1 - Q2 + Q3 - Q4
  
  deallocate(Q1)
  deallocate(Q2)
  deallocate(Q3)
  deallocate(Q4)


end subroutine 
!=================================================
!=================================================
subroutine commutator_221(r1,r2,r3,w1,w2)
  implicit none 
  
  type(full_ham) :: r1,r2,r3,w1,w2
  integer :: nh,np,nb,i,j,k,w,n,m,p,p2,h,h2,a
  real(8),allocatable,dimension(:,:) :: Q1,Q2
  
m=r2%msp
n=r2%nbody


do w=1,r2%nblock 

  nh=r2%mat(w)%nhh
  np=r2%mat(w)%npp
  nb=r2%mat(w)%nph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !zero out intermediates 
  w1%mat(w)%Vhhhh=0.d0
  w2%mat(w)%Vhhhh=0.d0
  w1%mat(w)%Vpppp=0.d0
  w2%mat(w)%Vpppp=0.d0
  w1%mat(w)%Vpphh=0.d0
  w2%mat(w)%Vpphh=0.d0
  w1%mat(w)%Vppph=0.d0
  w2%mat(w)%Vppph=0.d0
  w1%mat(w)%Vphhh=0.d0
  w2%mat(w)%Vphhh=0.d0
  w1%mat(w)%Vphph=0.d0
  w1%mat(w)%Vphph=0.d0



! MATRIX MULTIPLY TO GET INTERMEDIATES 
  
if (np*nh .ne. 0 )  then
!fpp 
call dgemm('N','T',np,np,nh,al,r1%mat(w)%Vpphh,np, &
     r2%mat(w)%Vpphh,np,bet,w1%mat(w)%Vpppp,np)
end if 

if (np*nb .ne. 0 ) then 
!fpp
call dgemm('T','N',nb,nb,np,al,r2%mat(w)%Vppph,np, &
     r1%mat(w)%Vppph,np,bet,w1%mat(w)%Vphph,nb) 
end if 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


if (nh*np .ne. 0) then
!fhh
call dgemm('T','N',nh,nh,np,al,r2%mat(w)%Vpphh,np, &
     r1%mat(w)%Vpphh,np,bet,w1%mat(w)%Vhhhh,nh)
end if 

if (nh*nb .ne. 0) then 
!fhh
call dgemm('N','T',nb,nb,nh,al,r1%mat(w)%Vphhh,nb, &
     r2%mat(w)%Vphhh,nb,bet,w2%mat(w)%Vphph,nb) 
end if 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!fph
if (nh*nb*np .ne. 0) then 

call dgemm('T','N',nb,nh,np,al,r2%mat(w)%Vppph,np, &
     r1%mat(w)%Vpphh,np,bet,w1%mat(w)%Vphhh,nb)
call dgemm('T','N',nb,nh,np,al,r1%mat(w)%Vppph,np, &
     r2%mat(w)%Vpphh,np,bet,w2%mat(w)%Vphhh,nb)
call dgemm('N','T',np,nb,nh,al,r1%mat(w)%Vpphh,np, &
     r2%mat(w)%Vphhh,nb,bet,w1%mat(w)%Vppph,np) 
call dgemm('N','T',np,nb,nh,al,r2%mat(w)%Vpphh,np, &
     r1%mat(w)%Vphhh,nb,bet,w2%mat(w)%Vppph,np)

end if 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end do 

!fpp

do p=1,m-n
   do p2=1,m-n
      
      do i=1,n
         ! NOTE THAT WE ARE DIRECTLY ACCESSING fpp from the structure, but
         ! using a "get" function to access Vhphp.
         ! I am just looping over the indices of fpp, but I need to add "n" which is
         ! the number of holes when i feed these indices to v_elem, which does not discriminate
         ! between particles and holes. (Also, v_elem is acting on the intermediates w1 and w2) 
         r3%fpp(p,p2) = r3%fpp(p,p2) + r1%herm * &
              v_elem(i,p2+n,i,p+n,w1) - r2%herm * v_elem(i,p+n,i,p2+n,w1) 
         
      end do 
      
      do a=n+1,m
         
          r3%fpp(p,p2) = r3%fpp(p,p2) - r1%herm * &
               v_elem(a,p2+n,a,p+n,w1) + r2%herm * v_elem(a,p+n,a,p2+n,w1)
         
      end do 
   end do 
end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!fhh    

do h=1,n
   do h2=1,n
      
      do i=1,n
         
         r3%fhh(h,h2) = r3%fhh(h,h2) + r1%herm * &
              v_elem(i,h2,i,h,w1) - r2%herm * v_elem(i,h,i,h2,w1) 
         
      end do 
      
      do a=n+1,m
         
         r3%fhh(h,h2) = r3%fhh(h,h2) - r1%herm * &
              v_elem(a,h2,a,h,w2) + r2%herm * v_elem(a,h,a,h2,w2)
         
      end do 
   end do 
end do 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!fph

do p=1,m-n
   do h=1,n
      
      do i=1,n
         
         r3%fph(p,h) = r3%fph(p,h) + r1%herm * &
              v_elem(i,p+n,i,h,w2) - r2%herm * v_elem(i,p+n,i,h,w1)
         
      end do 

      do a=n+1,m
         
         r3%fph(p,h) = r3%fph(p,h) - r1%herm * &
              v_elem(a,p+n,a,h,w2) + r2%herm * v_elem(a,p+n,a,h,w1)
         
      end do 
      
    end do 
end do 

!===============================================================

end subroutine 
!===============================================================
!===============================================================
subroutine commutator_122(r1,r2,r3) 
  implicit none 
  
  integer :: i,a,nh,np,nb,q,II,JJ,m,n,pre
  integer :: h1,h2,h3,h4,p1,p2,p3,p4
  type(full_ham) :: r1,r2,r3
  real(8) :: sm

  
  n=r2%nbody
  m=r2%msp

  do q=1,r2%nblock
     
     nh=r2%mat(q)%nhh
     np=r2%mat(q)%npp
     nb=r2%mat(q)%nph
           
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vhhhh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (nh > 0 ) then

     r3%mat(q)%Vhhhh = 0.d0
     do JJ=1,nh
        do II=1,nh

           ! get bra and ket indices 
           h1=r2%mat(q)%qnhh(II,1)
           h2=r2%mat(q)%qnhh(II,2)

           h3=r2%mat(q)%qnhh(JJ,1)
           h4=r2%mat(q)%qnhh(JJ,2)
           
           
           
           do i=1,n
            
              r3%mat(q)%Vhhhh(II,JJ) = r3%mat(q)%Vhhhh(II,JJ) + &
                   r1%fhh(h1,i) * v_elem(i,h2,h3,h4,r2) + &
                   r1%fhh(h2,i) * v_elem(h1,i,h3,h4,r2) - r1%herm * ( &
                   r1%fhh(h3,i) * v_elem(h1,h2,i,h4,r2) + &
                   r1%fhh(h4,i) * v_elem(h1,h2,h3,i,r2) ) - (&
                
                   r2%fhh(h1,i) * v_elem(i,h2,h3,h4,r1) + &
                   r2%fhh(h2,i) * v_elem(h1,i,h3,h4,r1) - r2%herm * ( &
                   r2%fhh(h3,i) * v_elem(h1,h2,i,h4,r1) + &
                   r2%fhh(h4,i) * v_elem(h1,h2,h3,i,r1) ) ) 
              
           end do 
           
           do a=1,m-n
            
              r3%mat(q)%Vhhhh(II,JJ) = r3%mat(q)%Vhhhh(II,JJ) + &
                   r1%herm * (r1%fph(a,h1) * v_elem(a+n,h2,h3,h4,r2) + &
                   r1%fph(a,h2) * v_elem(h1,a+n,h3,h4,r2) ) - (&
                   r1%fph(a,h3) * v_elem(h1,h2,a+n,h4,r2) + &
                   r1%fph(a,h4) * v_elem(h1,h2,h3,a+n,r2) ) - (&
                   
                   r2%herm * (r2%fph(a,h1) * v_elem(a+n,h2,h3,h4,r1) + &
                   r2%fph(a,h2) * v_elem(h1,a+n,h3,h4,r1) ) - ( &
                   r2%fph(a,h3) * v_elem(h1,h2,a+n,h4,r1) + &
                   r2%fph(a,h4) * v_elem(h1,h2,h3,a+n,r1) ) )
                   
           end do 
  
        end do 
    end do 
end if     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vpppp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ( np > 0 ) then     
  r3%mat(q)%Vpppp = 0.d0
     do JJ=1,np
        do II=1,np
           
           p1=r2%mat(q)%qnpp(II,1)
           p2=r2%mat(q)%qnpp(II,2)
           p3=r2%mat(q)%qnpp(JJ,1)
           p4=r2%mat(q)%qnpp(JJ,2)
                      
           do i=1,n
            
              r3%mat(q)%Vpppp(II,JJ) = r3%mat(q)%Vpppp(II,JJ) + &
                   r1%fph(p1-n,i) * v_elem(i,p2,p3,p4,r2) + &
                   r1%fph(p2-n,i) * v_elem(p1,i,p3,p4,r2) - r1%herm * ( &
                   r1%fph(p3-n,i) * v_elem(p1,p2,i,p4,r2) + &
                   r1%fph(p4-n,i) * v_elem(p1,p2,p3,i,r2) ) - (&
                
                   r2%fph(p1-n,i) * v_elem(i,p2,p3,p4,r1) + &
                   r2%fph(p2-n,i) * v_elem(p1,i,p3,p4,r1) - r2%herm * ( &
                   r2%fph(p3-n,i) * v_elem(p1,p2,i,p4,r1) + &
                   r2%fph(p4-n,i) * v_elem(p1,p2,p3,i,r1) ) ) 
              
           end do 
           
           do a=1,m-n
            
              r3%mat(q)%Vpppp(II,JJ) = r3%mat(q)%Vpppp(II,JJ) + &
                   r1%herm * (r1%fpp(a,p1-n) * v_elem(a+n,p2,p3,p4,r2) + &
                   r1%fpp(a,p2-n) * v_elem(p1,a+n,p3,p4,r2) ) - (&
                   r1%fpp(a,p3-n) * v_elem(p1,p2,a+n,p4,r2) + &
                   r1%fpp(a,p4-n) * v_elem(p1,p2,p3,a+n,r2) ) - (&
                   
                   r2%herm * ( r2%fpp(a,p1-n) * v_elem(a+n,p2,p3,p4,r1) + &
                   r2%fpp(a,p2-n) * v_elem(p1,a+n,p3,p4,r1) ) - ( &
                   r2%fpp(a,p3-n) * v_elem(p1,p2,a+n,p4,r1) + &
                   r2%fpp(a,p4-n) * v_elem(p1,p2,p3,a+n,r1) ) )
                   
           end do 
           
        end do 
    end do 
end if     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vpphh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (np*nh > 0) then    
  r3%mat(q)%Vpphh = 0.d0
     do JJ=1,nh
        do II=1,np
           
           p1=r2%mat(q)%qnpp(II,1)
           p2=r2%mat(q)%qnpp(II,2)

           h3=r2%mat(q)%qnhh(JJ,1)
           h4=r2%mat(q)%qnhh(JJ,2)
           
           do i=1,n
            
              r3%mat(q)%Vpphh(II,JJ) = r3%mat(q)%Vpphh(II,JJ) + &
                   r1%fph(p1-n,i) * v_elem(i,p2,h3,h4,r2) + &
                   r1%fph(p2-n,i) * v_elem(p1,i,h3,h4,r2) - r1%herm * ( &
                   r1%fhh(h3,i) * v_elem(p1,p2,i,h4,r2) + &
                   r1%fhh(h4,i) * v_elem(p1,p2,h3,i,r2) ) - (&
                
                   r2%fph(p1-n,i) * v_elem(i,p2,h3,h4,r1) + &
                   r2%fph(p2-n,i) * v_elem(p1,i,h3,h4,r1) - r2%herm * ( &
                   r2%fhh(h3,i) * v_elem(p1,p2,i,h4,r1) + &
                   r2%fhh(h4,i) * v_elem(p1,p2,h3,i,r1) ) ) 

           end do 
           
           do a=1,m-n
            
              r3%mat(q)%Vpphh(II,JJ) = r3%mat(q)%Vpphh(II,JJ) + &
                   r1%herm * (r1%fpp(a,p1-n) * v_elem(a+n,p2,h3,h4,r2) + &
                   r1%fpp(a,p2-n) * v_elem(p1,a+n,h3,h4,r2) ) - (&
                   r1%fph(a,h3) * v_elem(p1,p2,a+n,h4,r2) + &
                   r1%fph(a,h4) * v_elem(p1,p2,h3,a+n,r2) ) - (&
                   
                   r2%herm * (r2%fpp(a,p1-n) * v_elem(a+n,p2,h3,h4,r1) + &
                   r2%fpp(a,p2-n) * v_elem(p1,a+n,h3,h4,r1) ) - ( &
                   r2%fph(a,h3) * v_elem(p1,p2,a+n,h4,r1) + &
                   r2%fph(a,h4) * v_elem(p1,p2,h3,a+n,r1) ) )
                                 
           end do 
      
        end do 
    end do
    
end if 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vphhh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (nh*nb > 0 )  then 
     r3%mat(q)%Vphhh = 0.d0
     do JJ=1,nh
        do II=1,nb

           !this is stupid, but my code assumes that particle states are all
           !indexed higher than holes, despite the convention in the storage structure
           !that particles come before holes. dumb. I multiply by "pre=-1" to fix this.  
           p1 = r2%mat(q)%qnph(II,2)
           h2 = r2%mat(q)%qnph(II,1)  
           pre = -1
           
           h3=r2%mat(q)%qnhh(JJ,1)
           h4=r2%mat(q)%qnhh(JJ,2)
           
           
           
           do i=1,n
            
              r3%mat(q)%Vphhh(II,JJ) = r3%mat(q)%Vphhh(II,JJ) + pre* ( &
                   r1%fph(p1-n,i) * v_elem(i,h2,h3,h4,r2) + &
                   r1%fhh(h2,i) * v_elem(p1,i,h3,h4,r2) - r1%herm * ( &
                   r1%fhh(h3,i) * v_elem(p1,h2,i,h4,r2) + &
                   r1%fhh(h4,i) * v_elem(p1,h2,h3,i,r2) ) - (&
                
                   r2%fph(p1-n,i) * v_elem(i,h2,h3,h4,r1) + &
                   r2%fhh(h2,i) * v_elem(p1,i,h3,h4,r1) - r2%herm * ( &
                   r2%fhh(h3,i) * v_elem(p1,h2,i,h4,r1) + &
                   r2%fhh(h4,i) * v_elem(p1,h2,h3,i,r1) ) )  )
              
           end do 
           
           do a=1,m-n
            
              r3%mat(q)%Vphhh(II,JJ) = r3%mat(q)%Vphhh(II,JJ) + pre * ( &
                   r1%herm * (r1%fpp(a,p1-n) * v_elem(a+n,h2,h3,h4,r2) + &
                   r1%fph(a,h2) * v_elem(p1,a+n,h3,h4,r2) ) - (&
                   r1%fph(a,h3) * v_elem(p1,h2,a+n,h4,r2) + &
                   r1%fph(a,h4) * v_elem(p1,h2,h3,a+n,r2) ) - (&
                   
                   r2%herm * (r2%fpp(a,p1-n) * v_elem(a+n,h2,h3,h4,r1) + &
                   r2%fph(a,h2) * v_elem(p1,a+n,h3,h4,r1) ) - ( &
                   r2%fph(a,h3) * v_elem(p1,h2,a+n,h4,r1) + &
                   r2%fph(a,h4) * v_elem(p1,h2,h3,a+n,r1) ) ) ) 
                   
           end do 

        end do 
    end do 
end if 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vppph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (np*nb > 0) then    
  r3%mat(q)%Vppph = 0.d0
     do JJ=1,nb
        do II=1,np
           
           p1=r2%mat(q)%qnpp(II,1)
           p2=r2%mat(q)%qnpp(II,2)

           p3 = r2%mat(q)%qnph(JJ,2)
           h4 = r2%mat(q)%qnph(JJ,1)  
           pre = -1 
           
           do i=1,n
            
              r3%mat(q)%Vppph(II,JJ) = r3%mat(q)%Vppph(II,JJ) + pre * ( &
                   r1%fph(p1-n,i) * v_elem(i,p2,p3,h4,r2) + &
                   r1%fph(p2-n,i) * v_elem(p1,i,p3,h4,r2) - r1%herm * ( &
                   r1%fph(p3-n,i) * v_elem(p1,p2,i,h4,r2) + &
                   r1%fhh(h4,i) * v_elem(p1,p2,p3,i,r2) ) - (&
                
                   r2%fph(p1-n,i) * v_elem(i,p2,p3,h4,r1) + &
                   r2%fph(p2-n,i) * v_elem(p1,i,p3,h4,r1) - r2%herm * ( &
                   r2%fph(p3-n,i) * v_elem(p1,p2,i,h4,r1) + &
                   r2%fhh(h4,i) * v_elem(p1,p2,p3,i,r1) ) ) )
              
           end do 
           
           do a=1,m-n
            
              r3%mat(q)%Vppph(II,JJ) = r3%mat(q)%Vppph(II,JJ) + pre * ( &
                   r1%herm * (r1%fpp(a,p1-n) * v_elem(a+n,p2,p3,h4,r2) + &
                   r1%fpp(a,p2-n) * v_elem(p1,a+n,p3,h4,r2) ) - (&
                   r1%fpp(a,p3-n) * v_elem(p1,p2,a+n,h4,r2) + &
                   r1%fph(a,h4) * v_elem(p1,p2,p3,a+n,r2) ) - (&
                   
                   r2%herm * (r2%fpp(a,p1-n) * v_elem(a+n,p2,p3,h4,r1) + &
                   r2%fpp(a,p2-n) * v_elem(p1,a+n,p3,h4,r1) ) - ( &
                   r2%fpp(a,p3-n) * v_elem(p1,p2,a+n,h4,r1) + &
                   r2%fph(a,h4) * v_elem(p1,p2,p3,a+n,r1) ) ) )
                   
           end do 
   
        end do 
    end do 
      
end if    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vphph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (nb> 0) then    
  
  r3%mat(q)%Vphph = 0.d0
     do JJ=1,nb
        do II=1,nb
                     
           p1 = r2%mat(q)%qnph(II,2)
           h2 = r2%mat(q)%qnph(II,1)  
       
           p3 = r2%mat(q)%qnph(JJ,2)
           h4 = r2%mat(q)%qnph(JJ,1)  
           pre = 1
          
           do i=1,n
            
              r3%mat(q)%Vphph(II,JJ) = r3%mat(q)%Vphph(II,JJ) + pre * ( &
                   r1%fph(p1-n,i) * v_elem(i,h2,p3,h4,r2) + &
                   r1%fhh(h2,i) * v_elem(p1,i,p3,h4,r2) - r1%herm * ( &
                   r1%fph(p3-n,i) * v_elem(p1,h2,i,h4,r2) + &
                   r1%fhh(h4,i) * v_elem(p1,h2,p3,i,r2) ) - (&
                
                   r2%fph(p1-n,i) * v_elem(i,h2,p3,h4,r1) + &
                   r2%fhh(h2,i) * v_elem(p1,i,p3,h4,r1) - r2%herm * ( &
                   r2%fph(p3-n,i) * v_elem(p1,h2,i,h4,r1) + &
                   r2%fhh(h4,i) * v_elem(p1,h2,p3,i,r1) ) ) )
              
           end do 
      
      
           do a=1,m-n
            
              r3%mat(q)%Vphph(II,JJ) = r3%mat(q)%Vphph(II,JJ) + pre * ( &
                   r1%herm * (r1%fpp(a,p1-n) * v_elem(a+n,h2,p3,h4,r2) + &
                   r1%fph(a,h2) * v_elem(p1,a+n,p3,h4,r2)) - (&
                   r1%fpp(a,p3-n) * v_elem(p1,h2,a+n,h4,r2) + &
                   r1%fph(a,h4) * v_elem(p1,h2,p3,a+n,r2) ) - (&
                   
                   r2%herm * (r2%fpp(a,p1-n) * v_elem(a+n,h2,p3,h4,r1) + &
                   r2%fph(a,h2) * v_elem(p1,a+n,p3,h4,r1) ) - ( &
                   r2%fpp(a,p3-n) * v_elem(p1,h2,a+n,h4,r1) + &
                   r2%fph(a,h4) * v_elem(p1,h2,p3,a+n,r1) ) ) )
                   
           end do 
    
 
        end do 
    end do 

end if     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end do 

end subroutine 
!============================================================================
!============================================================================
subroutine commutator_222(r1,r2,r3,w1,HCC,ETACC) 
  use ME_general
  implicit none 
  
  integer :: i,j,k,q,nh,np,nb,n,m,a,II,JJ
  integer :: pre,h1,h2,h3,h4,p1,p2,p3,p4
  integer :: r31,r24,r32,r14,r41,r23,r42,r13
  integer :: q31,q24,q32,q14,q41,q23,q42,q13
  type(full_ham) :: r1,r2,r3,w1,w2
  type(cc_mat) :: HCC,ETACC,XCC,YCC
  real(8) :: X,Y,Z 
  
  call construct_cc_ints(ETACC,HCC,XCC,YCC) ! constructs ph ladder intermediates 
  
!!! matrix mults
  
  n=r2%nbody ! electrons
  m=r2%msp   ! sp orbitals 
  
  !$omp parallel do private(nh,np,nb,pre,h1,h2,h3,h4,p1,p2,p3,p4,II,JJ,i,a,X) &
  !$omp& shared(XCC,YCC,r1,r2,r3)  
do q=1,r2%nblock
   
   nh=r2%mat(q)%nhh 
   np=r2%mat(q)%npp
   nb=r2%mat(q)%nph

   call matrix_mult222(r1,r2,r3,w1,q,nh,np,nb) ! constructs pp and hh ladder intermediates

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!!! obnoxious part
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vhhhh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


do JJ=1,nh
   do II=1,nh
      
        h1=r2%mat(q)%qnhh(II,1)
        h2=r2%mat(q)%qnhh(II,2)
        h3=r2%mat(q)%qnhh(JJ,1)
        h4=r2%mat(q)%qnhh(JJ,2)

        ! get indices for CC (ph ladder) intermediate
        q31 = HCC%map(r2%msp*(h3-1)+h1,1) 
        r31 = HCC%map(r2%msp*(h3-1)+h1,2) 
        q24 = HCC%map(r2%msp*(h2-1)+h4,1) 
        r24 = HCC%map(r2%msp*(h2-1)+h4,2)        

        q32 = HCC%map(r2%msp*(h3-1)+h2,1) 
        r32 = HCC%map(r2%msp*(h3-1)+h2,2) 
        q14 = HCC%map(r2%msp*(h1-1)+h4,1) 
        r14 = HCC%map(r2%msp*(h1-1)+h4,2)        

        q41 = HCC%map(r2%msp*(h4-1)+h1,1) 
        r41 = HCC%map(r2%msp*(h4-1)+h1,2) 
        q23 = HCC%map(r2%msp*(h2-1)+h3,1) 
        r23 = HCC%map(r2%msp*(h2-1)+h3,2)        

        q42 = HCC%map(r2%msp*(h4-1)+h2,1) 
        r42 = HCC%map(r2%msp*(h4-1)+h2,2) 
        q13 = HCC%map(r2%msp*(h1-1)+h3,1) 
        r13 = HCC%map(r2%msp*(h1-1)+h3,2)        

        ! this is the ph channel 222 derivative: 
        X=  &
             XCC%mat(q31)%X(r31,r24) - XCC%mat(q32)%X(r32,r14) + &
             XCC%mat(q42)%X(r42,r13) - XCC%mat(q41)%X(r41,r23)
        X= X - &
             YCC%mat(q42)%X(r42,r13) + YCC%mat(q41)%X(r41,r23) - &
             YCC%mat(q31)%X(r31,r24) + YCC%mat(q32)%X(r32,r14)

        r3%mat(q)%Vhhhh(II,JJ)=r3%mat(q)%Vhhhh(II,JJ)+X
   end do 
end do


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vpppp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do JJ=1,np
   do II=1,np
      
        p1=r2%mat(q)%qnpp(II,1)
        p2=r2%mat(q)%qnpp(II,2)
        p3=r2%mat(q)%qnpp(JJ,1)
        p4=r2%mat(q)%qnpp(JJ,2)

        q31 = HCC%map(r2%msp*(p3-1)+p1,1) 
        r31 = HCC%map(r2%msp*(p3-1)+p1,2) 
        q24 = HCC%map(r2%msp*(p2-1)+p4,1) 
        r24 = HCC%map(r2%msp*(p2-1)+p4,2)        

        q32 = HCC%map(r2%msp*(p3-1)+p2,1) 
        r32 = HCC%map(r2%msp*(p3-1)+p2,2) 
        q14 = HCC%map(r2%msp*(p1-1)+p4,1) 
        r14 = HCC%map(r2%msp*(p1-1)+p4,2)        

        q41 = HCC%map(r2%msp*(p4-1)+p1,1) 
        r41 = HCC%map(r2%msp*(p4-1)+p1,2) 
        q23 = HCC%map(r2%msp*(p2-1)+p3,1) 
        r23 = HCC%map(r2%msp*(p2-1)+p3,2)        

        q42 = HCC%map(r2%msp*(p4-1)+p2,1) 
        r42 = HCC%map(r2%msp*(p4-1)+p2,2) 
        q13 = HCC%map(r2%msp*(p1-1)+p3,1) 
        r13 = HCC%map(r2%msp*(p1-1)+p3,2)        

        X=  &
             XCC%mat(q31)%X(r31,r24) - XCC%mat(q32)%X(r32,r14) + &
             XCC%mat(q42)%X(r42,r13) - XCC%mat(q41)%X(r41,r23)
        X= X - &
             YCC%mat(q42)%X(r42,r13) + YCC%mat(q41)%X(r41,r23) - &
             YCC%mat(q31)%X(r31,r24) + YCC%mat(q32)%X(r32,r14)

        r3%mat(q)%Vpppp(II,JJ)=r3%mat(q)%Vpppp(II,JJ)+X
        
   end do 
end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vpphh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


do JJ=1,nh
   do II=1,np
      
        p1=r2%mat(q)%qnpp(II,1)
        p2=r2%mat(q)%qnpp(II,2)
        h3=r2%mat(q)%qnhh(JJ,1)
        h4=r2%mat(q)%qnhh(JJ,2)

        q31 = HCC%map(r2%msp*(h3-1)+p1,1) 
        r31 = HCC%map(r2%msp*(h3-1)+p1,2) 
        q24 = HCC%map(r2%msp*(p2-1)+h4,1) 
        r24 = HCC%map(r2%msp*(p2-1)+h4,2)        

        q32 = HCC%map(r2%msp*(h3-1)+p2,1) 
        r32 = HCC%map(r2%msp*(h3-1)+p2,2) 
        q14 = HCC%map(r2%msp*(p1-1)+h4,1) 
        r14 = HCC%map(r2%msp*(p1-1)+h4,2)        

        q41 = HCC%map(r2%msp*(h4-1)+p1,1) 
        r41 = HCC%map(r2%msp*(h4-1)+p1,2) 
        q23 = HCC%map(r2%msp*(p2-1)+h3,1) 
        r23 = HCC%map(r2%msp*(p2-1)+h3,2)        

        q42 = HCC%map(r2%msp*(h4-1)+p2,1) 
        r42 = HCC%map(r2%msp*(h4-1)+p2,2) 
        q13 = HCC%map(r2%msp*(p1-1)+h3,1) 
        r13 = HCC%map(r2%msp*(p1-1)+h3,2)        

        X = XCC%mat(q31)%X(r31,r24) - XCC%mat(q32)%X(r32,r14) + &
             XCC%mat(q42)%X(r42,r13) - XCC%mat(q41)%X(r41,r23)
        
        X= X - YCC%mat(q42)%X(r42,r13) + YCC%mat(q41)%X(r41,r23) - &
             YCC%mat(q31)%X(r31,r24) + YCC%mat(q32)%X(r32,r14)

        r3%mat(q)%Vpphh(II,JJ)=r3%mat(q)%Vpphh(II,JJ)+X
   end do 
end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vppph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do JJ=1,nb
   do II=1,np
      
        p1=r2%mat(q)%qnpp(II,1)
        p2=r2%mat(q)%qnpp(II,2)

        p3=r2%mat(q)%qnph(JJ,2)
        h4=r2%mat(q)%qnph(JJ,1)
        pre=-1

     
        q31 = HCC%map(r2%msp*(p3-1)+p1,1) 
        r31 = HCC%map(r2%msp*(p3-1)+p1,2) 
        q24 = HCC%map(r2%msp*(p2-1)+h4,1) 
        r24 = HCC%map(r2%msp*(p2-1)+h4,2)        

        q32 = HCC%map(r2%msp*(p3-1)+p2,1) 
        r32 = HCC%map(r2%msp*(p3-1)+p2,2) 
        q14 = HCC%map(r2%msp*(p1-1)+h4,1) 
        r14 = HCC%map(r2%msp*(p1-1)+h4,2)        

        q41 = HCC%map(r2%msp*(h4-1)+p1,1) 
        r41 = HCC%map(r2%msp*(h4-1)+p1,2) 
        q23 = HCC%map(r2%msp*(p2-1)+p3,1) 
        r23 = HCC%map(r2%msp*(p2-1)+p3,2)        

        q42 = HCC%map(r2%msp*(h4-1)+p2,1) 
        r42 = HCC%map(r2%msp*(h4-1)+p2,2) 
        q13 = HCC%map(r2%msp*(p1-1)+p3,1) 
        r13 = HCC%map(r2%msp*(p1-1)+p3,2)        

        X=  &
             XCC%mat(q31)%X(r31,r24) - XCC%mat(q32)%X(r32,r14) + &
             XCC%mat(q42)%X(r42,r13) - XCC%mat(q41)%X(r41,r23)

        X= X - &
             YCC%mat(q42)%X(r42,r13) + YCC%mat(q41)%X(r41,r23) - &
             YCC%mat(q31)%X(r31,r24) + YCC%mat(q32)%X(r32,r14)

        r3%mat(q)%Vppph(II,JJ)=r3%mat(q)%Vppph(II,JJ)+X*pre

   end do 
end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vphhh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do JJ=1,nh
   do II=1,nb
      
        p1=r2%mat(q)%qnph(II,2)
        h2=r2%mat(q)%qnph(II,1)
        pre=-1

        h3=r2%mat(q)%qnhh(JJ,1)
        h4=r2%mat(q)%qnhh(JJ,2)


        q31 = HCC%map(r2%msp*(h3-1)+p1,1) 
        r31 = HCC%map(r2%msp*(h3-1)+p1,2) 
        q24 = HCC%map(r2%msp*(h2-1)+h4,1) 
        r24 = HCC%map(r2%msp*(h2-1)+h4,2)        

        q32 = HCC%map(r2%msp*(h3-1)+h2,1) 
        r32 = HCC%map(r2%msp*(h3-1)+h2,2) 
        q14 = HCC%map(r2%msp*(p1-1)+h4,1) 
        r14 = HCC%map(r2%msp*(p1-1)+h4,2)        

        q41 = HCC%map(r2%msp*(h4-1)+p1,1) 
        r41 = HCC%map(r2%msp*(h4-1)+p1,2) 
        q23 = HCC%map(r2%msp*(h2-1)+h3,1) 
        r23 = HCC%map(r2%msp*(h2-1)+h3,2)        

        q42 = HCC%map(r2%msp*(h4-1)+h2,1) 
        r42 = HCC%map(r2%msp*(h4-1)+h2,2) 
        q13 = HCC%map(r2%msp*(p1-1)+h3,1) 
        r13 = HCC%map(r2%msp*(p1-1)+h3,2)        

        X=  &
             XCC%mat(q31)%X(r31,r24) - XCC%mat(q32)%X(r32,r14) + &
             XCC%mat(q42)%X(r42,r13) - XCC%mat(q41)%X(r41,r23)

        X= X - &
             YCC%mat(q42)%X(r42,r13) + YCC%mat(q41)%X(r41,r23) - &
             YCC%mat(q31)%X(r31,r24) + YCC%mat(q32)%X(r32,r14)

        r3%mat(q)%Vphhh(II,JJ)=r3%mat(q)%Vphhh(II,JJ)+X*pre
      
   end do 
end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vphph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do JJ=1,nb
   do II=1,nb

        p1=r2%mat(q)%qnph(II,2)
        h2=r2%mat(q)%qnph(II,1)
        p3=r2%mat(q)%qnph(JJ,2)
        h4=r2%mat(q)%qnph(JJ,1)
        
        q31 = HCC%map(r2%msp*(p3-1)+p1,1) 
        r31 = HCC%map(r2%msp*(p3-1)+p1,2) 
        q24 = HCC%map(r2%msp*(h2-1)+h4,1) 
        r24 = HCC%map(r2%msp*(h2-1)+h4,2)        

        q32 = HCC%map(r2%msp*(p3-1)+h2,1) 
        r32 = HCC%map(r2%msp*(p3-1)+h2,2) 
        q14 = HCC%map(r2%msp*(p1-1)+h4,1) 
        r14 = HCC%map(r2%msp*(p1-1)+h4,2)        

        q41 = HCC%map(r2%msp*(h4-1)+p1,1) 
        r41 = HCC%map(r2%msp*(h4-1)+p1,2) 
        q23 = HCC%map(r2%msp*(h2-1)+p3,1) 
        r23 = HCC%map(r2%msp*(h2-1)+p3,2)        

        q42 = HCC%map(r2%msp*(h4-1)+h2,1) 
        r42 = HCC%map(r2%msp*(h4-1)+h2,2) 
        q13 = HCC%map(r2%msp*(p1-1)+p3,1) 
        r13 = HCC%map(r2%msp*(p1-1)+p3,2)        

        X=  &
             XCC%mat(q31)%X(r31,r24) - XCC%mat(q32)%X(r32,r14) + &
             XCC%mat(q42)%X(r42,r13) - XCC%mat(q41)%X(r41,r23)

        X= X - &
             YCC%mat(q42)%X(r42,r13) + YCC%mat(q41)%X(r41,r23) - &
             YCC%mat(q31)%X(r31,r24) + YCC%mat(q32)%X(r32,r14)

        r3%mat(q)%Vphph(II,JJ)=r3%mat(q)%Vphph(II,JJ)+X

   end do 
end do

end do 
!$omp end parallel do 

end subroutine commutator_222
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine matrix_mult222(r1,r2,r3,w1,q,nh,np,nb)
  implicit none 
  
  integer :: q,nh,np,nb
  type(full_ham) :: r1,r2,r3,w1
  
  ! this is not written well. Sorry.
  ! part of it is redundant from comm221
  
  if (nh*np > 0 ) then 
!Vhhhh
  
   call dgemm('T','N',nh,nh,np,al,r2%mat(q)%Vpphh,np,&
        r1%mat(q)%Vpphh,np,bet,w1%mat(q)%Vhhhh,nh) 

    r3%mat(q)%Vhhhh = r3%mat(q)%Vhhhh + ( r1%herm * &
        Transpose(w1%mat(q)%Vhhhh) - r2%herm * w1%mat(q)%Vhhhh) 

!Vpppp
   call dgemm('N','T',np,np,nh,al,r1%mat(q)%Vpphh,np,&
        r2%mat(q)%Vpphh,np,bet,w1%mat(q)%Vpppp,np)

   r3%mat(q)%Vpppp = r3%mat(q)%Vpppp + ( r1%herm * &
        Transpose(w1%mat(q)%Vpppp) - r2%herm * w1%mat(q)%Vpppp) 


!Vpphh
   call dgemm('N','N',np,nh,np,al,r1%mat(q)%Vpppp,np,&
        r2%mat(q)%Vpphh,np,bet,w1%mat(q)%Vpphh,np) 

   r3%mat(q)%Vpphh = r3%mat(q)%Vpphh + w1%mat(q)%Vpphh
  
   call dgemm('N','N',np,nh,np,al,r2%mat(q)%Vpppp,np,&
        r1%mat(q)%Vpphh,np,bet,w1%mat(q)%Vpphh,np) 
   
   r3%mat(q)%Vpphh = r3%mat(q)%Vpphh - w1%mat(q)%Vpphh

   call dgemm('N','N',np,nh,nh,al,r1%mat(q)%Vpphh,np,&
        r2%mat(q)%Vhhhh,nh,bet,w1%mat(q)%Vpphh,np)
 
   r3%mat(q)%Vpphh = r3%mat(q)%Vpphh - w1%mat(q)%Vpphh 
   
   call dgemm('N','N',np,nh,nh,al,r2%mat(q)%Vpphh,np,&
        r1%mat(q)%Vhhhh,nh,bet,w1%mat(q)%Vpphh,np)
 
   r3%mat(q)%Vpphh = r3%mat(q)%Vpphh + w1%mat(q)%Vpphh

end if 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (nh > 0 ) then
!Vhhhh
   
  
   call dgemm('N','N',nh,nh,nh,al,r1%mat(q)%Vhhhh,nh,&
        r2%mat(q)%Vhhhh,nh,bet,w1%mat(q)%Vhhhh,nh)
   
   r3%mat(q)%Vhhhh = r3%mat(q)%Vhhhh + ( r1%herm * r2%herm * &
        Transpose(w1%mat(q)%Vhhhh) - w1%mat(q)%Vhhhh)
   
end if    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (np > 0 ) then
!Vpppp 
   
   call dgemm('N','N',np,np,np,al,r2%mat(q)%Vpppp,np,&
        r1%mat(q)%Vpppp,np,bet,w1%mat(q)%Vpppp,np) 

   r3%mat(q)%Vpppp = r3%mat(q)%Vpppp + ( r1%herm * r2%herm * &
        Transpose(w1%mat(q)%Vpppp) - w1%mat(q)%Vpppp) 
end if 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (nb*np*nh > 0)  then 

!Vphhh
   
   call dgemm('T','N',nb,nh,np,al,r1%mat(q)%Vppph,np,&
        r2%mat(q)%Vpphh,np,bet,w1%mat(q)%Vphhh,nb) 

   r3%mat(q)%Vphhh = r3%mat(q)%Vphhh + r1%herm * w1%mat(q)%Vphhh
  
   call dgemm('T','N',nb,nh,np,al,r2%mat(q)%Vppph,np,&
        r1%mat(q)%Vpphh,np,bet,w1%mat(q)%Vphhh,nb) 
   
   r3%mat(q)%Vphhh = r3%mat(q)%Vphhh - r2%herm * w1%mat(q)%Vphhh

!Vppph
   
   call dgemm('N','T',np,nb,nh,al,r1%mat(q)%Vpphh,np,&
        r2%mat(q)%Vphhh,nb,bet,w1%mat(q)%Vppph,np)
   
   r3%mat(q)%Vppph = r3%mat(q)%Vppph - r2%herm * w1%mat(q)%Vppph 
   
   call dgemm('N','T',np,nb,nh,al,r2%mat(q)%Vpphh,np,&
        r1%mat(q)%Vphhh,nb,bet,w1%mat(q)%Vppph,np)
 
   r3%mat(q)%Vppph = r3%mat(q)%Vppph + r1%herm * w1%mat(q)%Vppph

end if 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (nb*nh > 0 ) then 
!Vphhh
   
   call dgemm('N','N',nb,nh,nh,al,r1%mat(q)%Vphhh,nb,&
        r2%mat(q)%Vhhhh,nh,bet,w1%mat(q)%Vphhh,nb)
 
   r3%mat(q)%Vphhh = r3%mat(q)%Vphhh - w1%mat(q)%Vphhh 
   
   call dgemm('N','N',nb,nh,nh,al,r2%mat(q)%Vphhh,nb,&
        r1%mat(q)%Vhhhh,nh,bet,w1%mat(q)%Vphhh,nb)
 
   r3%mat(q)%Vphhh = r3%mat(q)%Vphhh + w1%mat(q)%Vphhh

!Vphph
   
   call dgemm('N','T',nb,nb,nh,al,r1%mat(q)%Vphhh,nb,&
        r2%mat(q)%Vphhh,nb,bet,w1%mat(q)%Vphph,nb)
  
   r3%mat(q)%Vphph = r3%mat(q)%Vphph + (r1%herm *  &
        Transpose(w1%mat(q)%Vphph) - r2%herm * w1%mat(q)%Vphph) 
end if 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (nb*np > 0 ) then 
!Vppph
   
   call dgemm('N','N',np,nb,np,al,r1%mat(q)%Vpppp,np,&
        r2%mat(q)%Vppph,np,bet,w1%mat(q)%Vppph,np) 

   r3%mat(q)%Vppph = r3%mat(q)%Vppph + w1%mat(q)%Vppph
  
   call dgemm('N','N',np,nb,np,al,r2%mat(q)%Vpppp,np,&
        r1%mat(q)%Vppph,np,bet,w1%mat(q)%Vppph,np) 
   
   r3%mat(q)%Vppph = r3%mat(q)%Vppph - w1%mat(q)%Vppph

!Vphph
  
   call dgemm('T','N',nb,nb,np,al,r2%mat(q)%Vppph,np,&
        r1%mat(q)%Vppph,np,bet,w1%mat(q)%Vphph,nb) 
   
   r3%mat(q)%Vphph = r3%mat(q)%Vphph + (r1%herm * &
       Transpose(w1%mat(q)%Vphph) - r2%herm * w1%mat(q)%Vphph) 


end if 

end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================
subroutine calc_cc(H,HCC)
  implicit none

  type(full_ham) :: H
  type(cc_mat) :: HCC
  integer :: i,j,k,l,q,r,nb


  do q = 1, HCC%nblocks

     do r = 1, HCC%mat(q)%block_r

        i = HCC%mat(q)%qnab(r,1)
        j = HCC%mat(q)%qnab(r,2)
        
        do nb = 1, HCC%mat(q)%block_nb

           k = HCC%mat(q)%qnhp(nb,1)
           l = HCC%mat(q)%qnhp(nb,2)
           
           HCC%mat(q)%X(r,nb) = v_elem(j,k,i,l,H)
        end do
     end do
  end do
end subroutine calc_cc
!==================================================
!==================================================
real(8) function HQ_comm_1b(a,i,H,Q) 
  implicit none 
  
  integer :: a,i,b,j,N,M,c,k
  type(full_ham) :: H,Q 
  real(8) :: sm 
  
  sm = 0.d0 
  
  N = H%Nbody
  M = H%Msp

!11->1  
  do b = N+1,M
     sm = sm + f_elem(a,b,H) * Qf_elem(b,i,Q)

  end do
  
  do j = 1,N
     sm = sm - Qf_elem(a,j,Q) *f_elem(j,i,H) 
  end do 
  
   do b = N+1,M
      do j = 1, N
 ! 12 -> 1         
         sm = sm + f_elem(j,b,H) * Q_elem(a,b,i,j,Q)

 ! 21 -> 1 
         sm = sm + Qf_elem(b,j,Q) * V_elem(a,j,i,b,H) 
     end do
   end do 

!22 -> 1
  do b = N+1,M
     do c = b+1,M 
        do j = 1, N
           sm = sm + Q_elem(c,b,i,j,Q) * V_elem(a,j,c,b,H)           
        end do
     end do
  end do

  do b = N+1,M
     do k = 1,N 
        do j = k+1,N
           sm = sm - Q_elem(a,b,k,j,Q) * V_elem(k,j,i,b,H)             
        end do
     end do
  end do

  HQ_comm_1b = sm

end function 
!==================================================
!==================================================
real(8) function HQ_comm_2b(a,b,i,j,H,Q) 
  implicit none 
  
  integer :: a,i,b,j,N,M,c,k,l,d
  type(full_ham) :: H,Q 
  real(8) :: sm 
  
  sm = 0.d0 
  
  N = H%Nbody
  M = H%Msp


 ! 12->2  
  do c = N+1,M
     sm = sm + f_elem(a,c,H) * Q_elem(c,b,i,j,Q) - &
          f_elem(b,c,H) * Q_elem(c,a,i,j,Q) 
     sm = sm + Qf_elem(c,i,Q) * V_elem(a,b,c,j,H) - &
          Qf_elem(c,j,Q) * V_elem(a,b,c,i,H)
  end do

  do k = 1,N
     sm = sm - Q_elem(a,b,k,j,Q) *f_elem(k,i,H)  + &
          Q_elem(a,b,k,i,Q) *f_elem(k,j,H) 
     sm = sm - V_elem(b,k,j,i,H) *Qf_elem(a,k,Q) + &
          V_elem(a,k,j,i,H) *Qf_elem(b,k,Q)
  end do 

 !22 -> 2
  
  do c = N+1,M
     do d = c+1,M
        sm = sm + V_elem(a,b,c,d,H) * Q_elem(c,d,i,j,Q)
     end do
  end do 

  do k = 1,N
     do l = k+1,N
        
        sm = sm + V_elem(k,l,i,j,H) * Q_elem(a,b,k,l,Q)
     end do 
  end do 

  do k = 1, N
     do c = N+1,M
        
        sm = sm - V_elem(b,k,c,i,H) * Q_elem(a,c,k,j,Q) 
        sm = sm + V_elem(a,k,c,i,H) * Q_elem(b,c,k,j,Q) 
        sm = sm + V_elem(b,k,c,j,H) * Q_elem(a,c,k,i,Q) 
        sm = sm - V_elem(a,k,c,j,H) * Q_elem(b,c,k,i,Q) 
        
      end do 
   end do 

  HQ_comm_2b = sm

end function 

!==================================================
!==================================================
real(8) function HQ_comm_1p(a,H,v) 
  implicit none 
  
  integer :: a,i,b,j,N,M,c,k
  type(full_ham) :: H
  real(8) :: sm 
  real(8),dimension(:) :: v
  
  N = H%nbody
  M = H%msp
  sm = 0.d0 
  
  do b = N+1,M
     sm = sm + f_elem(a,b,H)*get_X1(b,v,H)
  end do 
  
  do b =N+1,M
     do i = 1,N
        sm = sm + f_elem(i,b,H)*get_X3(a,b,i,v,H)
     end do 
  end do 
  
  do b = N+1,M
     do c = b+1,M
        do i = 1,N
           sm = sm + v_elem(a,i,b,c,H)*get_X3(b,c,i,v,H) 
        end do 
     end do
  end do
           
  HQ_comm_1p = sm

end function HQ_comm_1p
!==================================================
!==================================================
real(8) function HQ_comm_1h(i,H,v) 
  implicit none 
  
  integer :: a,i,b,j,N,M,c,k
  type(full_ham) :: H
  real(8) :: sm 
  real(8),dimension(:) :: v
  
  N = H%nbody
  M = H%msp
  sm = 0.d0 
  
  do j = 1,N
     sm = sm - f_elem(j,i,H)*get_X1(j,v,H)
  end do 
  
  do b =N+1,M
     do j = 1,N
        sm = sm + f_elem(j,b,H)*get_X3(i,j,b,v,H)
     end do 
  end do 
  
  do j = 1,N
     do k = j+1,N
        do b = N+1,M
           sm = sm + v_elem(j,k,b,i,H)*get_X3(j,k,b,v,H) 
        end do 
     end do
  end do
           
  HQ_comm_1h = sm

end function HQ_comm_1h
!==================================================
!==================================================
real(8) function HQ_comm_2p1h(a,b,i,H,v) 
  implicit none 
  
  integer :: a,i,b,j,N,M,c,k,l,d
  type(full_ham) :: H 
  real(8) :: sm 
  real(8),dimension(:) :: v
  
  N = H%nbody
  M = H%msp
  sm = 0.d0 
  
  do c = N+1,M
     sm = sm + v_elem(a,b,c,i,H)*get_X1(c,v,H) &
          + f_elem(a,c,H)*get_X3(c,b,i,v,H) &
          - f_elem(b,c,H)*get_X3(c,a,i,v,H)
  end do 
  
  do j = 1, N
     sm = sm - f_elem(j,i,H)*get_X3(a,b,j,v,H)
  end do 
  
  do c = N+1,M
     do d = c+1,M
        sm = sm + v_elem(a,b,c,d,H)*get_X3(c,d,i,v,H)
     end do 
  end do 
  
  do c = N+1,M 
     do j = 1, N
        
        sm = sm - v_elem(a,j,c,i,H)*get_X3(c,b,j,v,H)&
             + v_elem(b,j,c,i,H)*get_X3(c,a,j,v,H)
     end do 
  end do 

  HQ_comm_2p1h = sm 

end function HQ_comm_2p1h
!==================================================
!==================================================
real(8) function HQ_comm_2h1p(i,j,a,H,v) 
  implicit none 
  
  integer :: a,i,b,j,N,M,c,k,l,d
  type(full_ham) :: H 
  real(8) :: sm 
  real(8),dimension(:) :: v
  
  N = H%nbody
  M = H%msp
  sm = 0.d0 
  
  do k = 1,N
     sm = sm + v_elem(a,k,i,j,H)*get_X1(k,v,H) &
          - f_elem(k,j,H)*get_X3(i,k,a,v,H) &
          + f_elem(k,i,H)*get_X3(j,k,a,v,H)
  end do 
  
  do b = N+1, M
     sm = sm + f_elem(a,b,H)*get_X3(i,j,b,v,H)
  end do 
  
  do k = 1,N
     do l = k+1,N
        sm = sm + v_elem(k,l,i,j,H)*get_X3(k,l,a,v,H)
     end do 
  end do 
  
  do b = N+1,M 
     do k = 1, N
        
        sm = sm - v_elem(a,k,b,i,H)*get_X3(k,j,b,v,H)&
             + v_elem(a,k,b,j,H)*get_X3(k,i,b,v,H)
     end do 
  end do 

  HQ_comm_2h1p = sm 

end function 
!===============================================================
!===============================================================
real(8) function W_223(p,q,r,s,t,u,r1,r2) 
  implicit none 
  
  integer :: p,q,r,s,t,u,a,i,II,JJ,M3,m,n
  type(full_ham) :: r1,r2
  real(8) :: sm
  
  m = r1%msp
  n = r1%nbody
  
  sm = 0.d0
  do a = n+1,m
     
     sm = sm + v_elem(p,q,a,t,r1)*Q_elem(a,r,s,u,r2) - &
          v_elem(p,q,a,u,r1)*Q_elem(a,r,s,t,r2) - &
          v_elem(p,q,a,s,r1)*Q_elem(a,r,t,u,r2) - &
          v_elem(r,q,a,t,r1)*Q_elem(a,p,s,u,r2) + &
          v_elem(r,q,a,s,r1)*Q_elem(a,p,t,u,r2) + &
          v_elem(r,q,a,u,r1)*Q_elem(a,p,s,t,r2) - &
          v_elem(p,r,a,t,r1)*Q_elem(a,q,s,u,r2) + &
          v_elem(p,r,a,s,r1)*Q_elem(a,q,t,u,r2) + &
          v_elem(p,r,a,u,r1)*Q_elem(a,q,s,t,r2)
     
  end do
      
  do i = 1, n
          
     sm = sm - (v_elem(i,q,s,t,r1)*Q_elem(p,r,i,u,r2) - &
          v_elem(i,q,s,u,r1)*Q_elem(p,r,i,t,r2) - &
          v_elem(i,q,u,t,r1)*Q_elem(p,r,i,s,r2) - &
          v_elem(i,p,s,t,r1)*Q_elem(q,r,i,u,r2) + &
          v_elem(i,p,s,u,r1)*Q_elem(q,r,i,t,r2) + &
          v_elem(i,p,u,t,r1)*Q_elem(q,r,i,s,r2) - &
          v_elem(i,r,s,t,r1)*Q_elem(p,q,i,u,r2) + &
          v_elem(i,r,s,u,r1)*Q_elem(p,q,i,t,r2) + &
          v_elem(i,r,u,t,r1)*Q_elem(p,q,i,s,r2) ) 
  end do
  
  W_223 = sm 
end function W_223
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================
function EOM_PERTURBATIVE_TRIPLES(H,Q) result( DE ) 
  implicit none 
  
  integer :: ML,MS,n,m,la,lab,labc,a,b,c,i,j,k
  integer :: sa,sab,sabc,li,lij,lijk,si,sij,sijk
  type(full_ham),intent(IN) :: H,Q
  real(8) :: deltaE,fii,fjj,fkk,faa,fbb,fcc,Gabab,Gijij
  real(8) :: Gacac,Gbcbc,Gikik,Gjkjk,Giaia,Gjaja,Gkaka
  real(8) :: Gibib,Gjbjb,Gkbkb,Gicic,Gjcjc,Gkckc,W,denom,n3p3h,X
  real(8) :: DE(2)
  
  ML=H%MLtarg
  MS=H%MStarg
  m = H%msp
  n = H%nbody
 
  deltaE = 0.d0 
  n3p3h = 0.d0 
  do a = n+1,m 
     la = H%states(a,2)
     sa = H%states(a,3)
     faa = f_elem(a,a,H)
     
     do b = a+1,m
        lab = la+H%states(b,2)
        sab = sa+H%states(b,3)
        fbb = f_elem(b,b,H)
        Gabab = v_elem(a,b,a,b,H) 
        
        do c = b+1,m
           labc = lab + H%states(c,2) 
           sabc = sab + H%states(c,3) 
  !         print*, ml,ms,labc,sabc
           
 !          IF ( labc .ne. ml ) cycle
  !         IF ( sabc .ne. ms ) cycle
           

           fcc = f_elem(c,c,H) 
           Gacac = v_elem(a,c,a,c,H)
           Gbcbc = v_elem(b,c,b,c,H) 

  

           do i = 1,n 
              li = H%states(i,2)
              si = H%states(i,3)
              fii = f_elem(i,i,H)
             
              Giaia = v_elem(i,a,i,a,H) 
              Gibib = v_elem(i,b,i,b,H)
              Gicic = v_elem(i,c,i,c,H)              
             
              
              do j = i+1,n
                 lij = li+H%states(j,2)
                 sij = si+H%states(j,3)
                 fjj = f_elem(j,j,H)
                 Gijij = v_elem(i,j,i,j,H)      
                 
                 Gjaja = v_elem(j,a,j,a,H) 
                 Gjbjb = v_elem(j,b,j,b,H)
                 Gjcjc = v_elem(j,c,j,c,H)              
      
      
                 do k = j+1,n
                    lijk = lij + H%states(k,2) 
                    sijk = sij + H%states(k,3) 
     
                    IF ( labc-lijk .ne. Ml  ) cycle
                    IF ( sabc-sijk .ne. MS ) cycle

                    fkk = f_elem(k,k,H)
                    Gikik = v_elem(i,k,i,k,H)
                    Gjkjk = v_elem(j,k,j,k,H)

                    Gkaka = v_elem(k,a,k,a,H) 
                    Gkbkb = v_elem(k,b,k,b,H)
                    Gkckc = v_elem(k,c,k,c,H)              
      
                    
                    
                    W=W_223(a,b,c,i,j,k,H,Q) 
                    denom = Q%E0-(faa+fbb+fcc-fii-fjj-fkk+Gabab+&
                         Gacac+Gbcbc+Gijij+Gikik+Gjkjk-Giaia&
                         -Gibib-Gicic-Gjaja-Gjbjb-Gjcjc-Gkaka-&
                         Gkbkb-Gkckc) 
                    X = W/denom 

                    deltaE = deltaE + W*X
                    n3p3h = n3p3h + X*X 
                    
                 end do 
              end do
           end do
        end do
     end do
  end do
         
  DE(1) = deltaE
  DE(2) = n3p3h 
  
end function EOM_PERTURBATIVE_TRIPLES


real(8) function EOM_off_diagonal(Q1,H,Q2)
  implicit none 
  
  integer :: ML,MS,n,m,la,lab,labc,a,b,c,i,j,k,l,d
  integer :: sa,sab,sabc,li,lij,lijk,si,sij,sijk
  type(full_ham),intent(IN) :: H,Q1,Q2
  real(8) :: sm ,smx

  n = H%nbody
  m = H%msp
  ! 1p1h and 1p1h 

  sm = 0.d0 
  do i = 1, n 
     do a = n+1,m 
        
        do j = 1, n
           do b = n+1,m 
              
              smx = v_elem(a,j,i,b,H)  
              if (i==j) then 
                 smx = smx + f_elem(a,b,H)
              end if
              
              if (a==b) then 
                 smx = smx - f_elem(j,i,H)
              end if
              
              sm = sm + Qf_elem(a,i,Q1) * Qf_elem(b,j,Q2) * smx 
              
           end do
        end do
     end do
  end do

  ! 1p1h and 2p2h 
  
  do i = 1, n
     do j = i+1, n
        do a = n+1,m
           do b = a+1, m
              
              do k = 1, n
                 do c = n+1,m
                    
                    smx = 0.d0 
                    if (i==k) then 
                       smx = smx + v_Elem(a,b,c,j,H) 
                       if (b==c) then 
                          smx = smx+f_elem(a,j,H) 
                       end if 

                       if (a==c) then 
                          smx = smx+f_elem(b,j,H) 
                       end if                       
                    end if
                      
                    if (j==k) then 
                       smx = smx - v_elem(a,b,c,i,H)
                       if (b==c) then 
                          smx = smx + f_elem(a,i,H)
                       end if

                       if (a==c) then 
                          smx = smx + f_elem(b,i,H)
                       end if
                    end if 
                    
                    if (a==c) then 
                       smx = smx - v_Elem(k,b,i,j,H)
                    end if 
                    
                    if (b==c) then 
                       smx = smx + v_elem(k,a,i,j,H) 
                    end if 
                    
                    sm = sm + smx*(Qf_elem(c,k,Q1)*Q_Elem(a,b,i,j,Q2) + &
                         Qf_elem(c,k,Q2)*Q_Elem(a,b,i,j,Q1) ) 
                    
                 end do
              end do
           end do
        end do
     end do
  end do
  
  do i = 1, n
     do j = i+1,n
        do a = n+1,m
           do b = a+1,m
              
              
              do k = 1, n
                 do l = k+1,n
                    do c = n+1,m
                       do d = c+1,m
                          
                          smx = 0.d0 

                          if ( i == k ) then 
                      
                             if ( j == l) then                                 
                                smx = smx+ v_Elem(a,b,c,d,H)
                             end if 
                             
                             if ( b == d) then 
                                smx = smx - v_elem(a,l,c,j,H) 
                             end if 
                             
                             if ( b == c) then 
                                smx = smx + v_elem(a,l,d,j,H) 
                             end if 
                             
                             if ( a == d) then 
                                smx = smx + v_elem(b,l,c,j,H) 
                             end if 
                             
                             if ( a == c) then 
                                smx = smx - v_elem(b,l,d,j,H) 
                             end if
 
                          end if
                          
                          if ( j == k ) then 
                             
                             if ( b == d) then 
                                smx = smx + v_elem(a,l,c,i,H) 
                             end if 
                             
                             if ( b == c) then 
                                smx = smx - v_elem(a,l,d,i,H) 
                             end if 
                             
                             if ( a == d) then 
                                smx = smx - v_elem(b,l,c,i,H) 
                             end if 
                             
                             if ( a == c) then 
                                smx = smx + v_elem(b,l,d,i,H) 
                             end if 
                          end if
                             
                          
                          if (  i == l ) then 

                             if ( j == k) then                                 
                                smx = smx - v_Elem(a,b,c,d,H)
                             end if 
                             
                             if ( b == d) then 
                                smx = smx + v_elem(a,k,c,j,H) 
                             end if 
                             
                             if ( b == c) then 
                                smx = smx - v_elem(a,k,d,j,H) 
                             end if 
                             
                             if ( a == d) then 
                                smx = smx - v_elem(b,k,c,j,H) 
                             end if 
                             
                             if ( a == c) then 
                                smx = smx + v_elem(b,k,d,j,H) 
                             end if 
                          end if
  
                          if ( j == l ) then 
                             
                             if ( b == d) then 
                                smx = smx - v_elem(a,k,c,i,H) 
                             end if 
                             
                             if ( b == c) then 
                                smx = smx + v_elem(a,k,d,i,H) 
                             end if 
                             
                             if ( a == d) then 
                                smx = smx + v_elem(b,k,c,i,H) 
                             end if 
                             
                             if ( a == c) then 
                                smx = smx - v_elem(b,k,d,i,H) 
                             end if 
                          
                          end if
                          
                          if (a==c) then 
                             if (b==d) then 
                                
                                smx =smx+v_elem(i,j,k,l,H)
                          
                             end if
                          end if 
                          
                          if (a==d) then 
                             if (b==c) then 
                                
                                smx =smx-v_elem(i,j,k,l,H)
                          
                             end if
                          end if
                          
                          sm = sm + smx * Q_elem(a,b,i,j,Q1) * Q_elem(c,d,k,l,Q2) 
                          
                          
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do
  
  EOM_off_diagonal = sm
end function EOM_off_diagonal

end module commutators

