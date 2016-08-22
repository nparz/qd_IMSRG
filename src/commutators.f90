module commutators
  use ME_general
  implicit none
  
  real(8), parameter :: al =1.d0, bet=0.d0
  !!! THIS MODULE CONTAINS COMMUTATION ROUTINES. 
  !!! IF YOU ARE COMMUTING A GENERATOR WITH SOMETHING, 
  !!! USE R1 AS THE GENERATOR.
 
contains

subroutine xcommutator_110(m,n,r1,r2,C)
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
subroutine xcommutator_220(r1,r2,c) 
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
subroutine xcommutator_121(r1,r2,r3) 
  implicit none 
  
  type(full_ham) :: r1,r2,r3
  integer :: n,m,i,a,h1,h2,p1,p2

  n=r2%nbody
  m=r2%msp-r2%nbody
  
  do h2=1,n
     do h1=1,n
     
        call sub121_hh(r1,r2,r3,h1,h2,n,m) 
        
     end do

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     

     do p1 = 1,m   

        call sub121_ph(r1,r2,r3,p1,h2,n,m) 
        
     end do 
  end do 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  do p2=1,m
     do p1=1,m
        
        call sub121_pp(r1,r2,r3,p1,p2,n,m)

     end do 
  end do 

end subroutine 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine sub121_hh(r1,r2,r3,h1,h2,n,m)
  implicit none 
  
  type(full_ham) :: r1,r2,r3
  integer :: n,m,i,a,h1,h2
  
  do i=1,n
     do a=1,m
        
        r3%fhh(h1,h2) = r3%fhh(h1,h2) + r1%fph(a,i) * &
    ( r1%herm * v_elem(a+n,h1,i,h2,r2) - v_elem(i,h1,a+n,h2,r2) ) &
    - r2%fph(a,i) * ( r2%herm * v_elem(a+n,h1,i,h2,r1) - &
    v_elem(i,h1,a+n,h2,r1) )
               
        
      end do 
  end do 

end subroutine 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine sub121_ph(r1,r2,r3,p1,h2,n,m) 
  implicit none 
  
  type(full_ham) :: r1,r2,r3
  integer :: n,m,i,a,h2,p1
  
  do i=1,n
     do a=1,m
        
        r3%fph(p1,h2) = r3%fph(p1,h2) + r1%fph(a,i) * &
    ( r1%herm * v_elem(a+n,p1+n,i,h2,r2) - v_elem(i,p1+n,a+n,h2,r2) ) &
    - r2%fph(a,i) * ( r2%herm *  v_elem(a+n,p1+n,i,h2,r1) - &
    v_elem(i,p1+n,a+n,h2,r1) )
         
      end do 
  end do 

end subroutine 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine sub121_pp(r1,r2,r3,p1,p2,n,m)
  implicit none
  
  type(full_ham) :: r1,r2,r3
  integer :: n,m,i,a,p1,p2
  
  do i=1,n
     do a=1,m
        
          r3%fpp(p1,p2) = r3%fpp(p1,p2) + r1%fph(a,i) * &
    ( r1%herm * v_elem(a+n,p1+n,i,p2+n,r2) - v_elem(i,p1+n,a+n,p2+n,r2) ) &
    - r2%fph(a,i) * ( r2%herm * v_elem(a+n,p1+n,i,p2+n,r1) - &
    v_elem(i,p1+n,a+n,p2+n,r1) )
        
      end do 
  end do      

end subroutine 
!===================================================================
!===================================================================
subroutine xcommutator_111(r1,r2,r3)
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
subroutine xcommutator_221(r1,r2,r3,w1,w2)
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
!fpp

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
subroutine xcommutator_122(r1,r2,r3) 
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
           
           h1=(r2%mat(q)%qnhh(II,1))
           h2=(r2%mat(q)%qnhh(II,2))
           h3=(r2%mat(q)%qnhh(JJ,1))
           h4=(r2%mat(q)%qnhh(JJ,2))
           
           
           
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
           
           p1=(r2%mat(q)%qnpp(II,1))
           p2=(r2%mat(q)%qnpp(II,2))
           p3=(r2%mat(q)%qnpp(JJ,1))
           p4=(r2%mat(q)%qnpp(JJ,2))
                      
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
           
           p1=(r2%mat(q)%qnpp(II,1))
           p2=(r2%mat(q)%qnpp(II,2))
           h3=(r2%mat(q)%qnhh(JJ,1))
           h4=(r2%mat(q)%qnhh(JJ,2))
           
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
           
           p1=max((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
           h2=min((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
           pre=1
           if ( p1 == (r2%mat(q)%qnph(II,2)) ) pre= -1

           h3=(r2%mat(q)%qnhh(JJ,1))
           h4=(r2%mat(q)%qnhh(JJ,2))
           
           
           
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
           
           p1=(r2%mat(q)%qnpp(II,1))
           p2=(r2%mat(q)%qnpp(II,2))
           p3=max((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
           h4=min((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
           pre=1
           if ( p3 == (r2%mat(q)%qnph(JJ,2)) ) pre= -1
           
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
           
           
           p1=max((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
           h2=min((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
           pre=1
           if ( p1 == (r2%mat(q)%qnph(II,2)) ) pre= -1
           p3=max((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
           h4=min((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
           if ( p3 == (r2%mat(q)%qnph(JJ,2)) ) pre= -1*pre
           
       
          
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
subroutine xcommutator_222(r1,r2,r3,w1) 
  implicit none 
  
  integer :: i,j,k,q,nh,np,nb,n,m,a,II,JJ
  integer :: pre,h1,h2,h3,h4,p1,p2,p3,p4
  type(full_ham) :: r1,r2,r3,w1,w2
  
  
!!! matrix mults

n=r2%nbody
m=r2%msp
  

!$omp parallel do private(nh,np,nb,pre,h1,h2,h3,h4,p1,p2,p3,p4,II,JJ,i,a)  
do q=1,r2%nblock
   
   nh=r2%mat(q)%nhh
   np=r2%mat(q)%npp
   nb=r2%mat(q)%nph

   call matrix_mult222(r1,r2,r3,w1,q,nh,np,nb)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!!! obnoxious part
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vhhhh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


do JJ=1,nh
   do II=1,nh
      
        h1=(r2%mat(q)%qnhh(II,1))
        h2=(r2%mat(q)%qnhh(II,2))
        h3=(r2%mat(q)%qnhh(JJ,1))
        h4=(r2%mat(q)%qnhh(JJ,2))

        call sub_222_hhhh(r1,r2,r3,q,h1,h2,h3,h4,II,JJ,n,m) 

   end do 
end do


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vpppp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do JJ=1,np
   do II=1,np
      
        p1=(r2%mat(q)%qnpp(II,1))
        p2=(r2%mat(q)%qnpp(II,2))
        p3=(r2%mat(q)%qnpp(JJ,1))
        p4=(r2%mat(q)%qnpp(JJ,2))
           
        call sub_222_pppp(r1,r2,r3,q,p1,p2,p3,p4,II,JJ,n,m) 

   end do 
end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vpphh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


do JJ=1,nh
   do II=1,np
      
        p1=(r2%mat(q)%qnpp(II,1))
        p2=(r2%mat(q)%qnpp(II,2))
        h3=(r2%mat(q)%qnhh(JJ,1))
        h4=(r2%mat(q)%qnhh(JJ,2))

        call sub_222_pphh(r1,r2,r3,q,p1,p2,h3,h4,II,JJ,n,m) 

   end do 
end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vppph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do JJ=1,nb
   do II=1,np
      
        p1=(r2%mat(q)%qnpp(II,1))
        p2=(r2%mat(q)%qnpp(II,2))
        p3=max((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
        h4=min((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
        pre=1
        if ( p3 == (r2%mat(q)%qnph(JJ,2)) ) pre= -1
     

        call sub_222_ppph(r1,r2,r3,q,p1,p2,p3,h4,II,JJ,n,m,pre)

   end do 
end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vphhh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do JJ=1,nh
   do II=1,nb
      
        p1=max((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
        h2=min((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
        pre=1
        if ( p1 == (r2%mat(q)%qnph(II,2)) ) pre= -1
        h3=(r2%mat(q)%qnhh(JJ,1))
        h4=(r2%mat(q)%qnhh(JJ,2))
       
        call sub_222_phhh(r1,r2,r3,q,p1,h2,h3,h4,II,JJ,n,m,pre)
      
   end do 
end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vphph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do JJ=1,nb
   do II=1,nb
      
        p1=max((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
        h2=min((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
        pre=1
        if ( p1 == (r2%mat(q)%qnph(II,2)) ) pre= -1
        p3=max((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
        h4=min((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
        if ( p3 == (r2%mat(q)%qnph(JJ,2)) ) pre= -1 * pre
        
        call sub_222_phph(r1,r2,r3,q,p1,h2,p3,h4,II,JJ,n,m,pre)

   end do 
end do

end do 
!$omp end parallel do 


end subroutine
!============================================================================
!============================================================================
subroutine xcommutator_222_nonpar(r1,r2,r3,w1) 
  implicit none 
  
  integer :: i,j,k,q,nh,np,nb,n,m,a,II,JJ
  integer :: pre,h1,h2,h3,h4,p1,p2,p3,p4
  type(full_ham) :: r1,r2,r3,w1,w2
  
  
!!! matrix mults

n=r2%nbody
m=r2%msp
  
do q=1,r2%nblock
   
   nh=r2%mat(q)%nhh
   np=r2%mat(q)%npp
   nb=r2%mat(q)%nph

   call matrix_mult222(r1,r2,r3,w1,q,nh,np,nb)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!!! obnoxious part
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vhhhh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


do JJ=1,nh
   do II=1,nh
      
        h1=(r2%mat(q)%qnhh(II,1))
        h2=(r2%mat(q)%qnhh(II,2))
        h3=(r2%mat(q)%qnhh(JJ,1))
        h4=(r2%mat(q)%qnhh(JJ,2))

        call sub_222_hhhh(r1,r2,r3,q,h1,h2,h3,h4,II,JJ,n,m) 

   end do 
end do


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vpppp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do JJ=1,np
   do II=1,np
      
        p1=(r2%mat(q)%qnpp(II,1))
        p2=(r2%mat(q)%qnpp(II,2))
        p3=(r2%mat(q)%qnpp(JJ,1))
        p4=(r2%mat(q)%qnpp(JJ,2))
           
        call sub_222_pppp(r1,r2,r3,q,p1,p2,p3,p4,II,JJ,n,m) 

   end do 
end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vpphh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


do JJ=1,nh
   do II=1,np
      
        p1=(r2%mat(q)%qnpp(II,1))
        p2=(r2%mat(q)%qnpp(II,2))
        h3=(r2%mat(q)%qnhh(JJ,1))
        h4=(r2%mat(q)%qnhh(JJ,2))

        call sub_222_pphh(r1,r2,r3,q,p1,p2,h3,h4,II,JJ,n,m) 

   end do 
end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vppph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do JJ=1,nb
   do II=1,np
      
        p1=(r2%mat(q)%qnpp(II,1))
        p2=(r2%mat(q)%qnpp(II,2))
        p3=max((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
        h4=min((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
        pre=1
        if ( p3 == (r2%mat(q)%qnph(JJ,2)) ) pre= -1
     

        call sub_222_ppph(r1,r2,r3,q,p1,p2,p3,h4,II,JJ,n,m,pre)

   end do 
end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vphhh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do JJ=1,nh
   do II=1,nb
      
        p1=max((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
        h2=min((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
        pre=1
        if ( p1 == (r2%mat(q)%qnph(II,2)) ) pre= -1
        h3=(r2%mat(q)%qnhh(JJ,1))
        h4=(r2%mat(q)%qnhh(JJ,2))
       
        call sub_222_phhh(r1,r2,r3,q,p1,h2,h3,h4,II,JJ,n,m,pre)
      
   end do 
end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vphph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

do JJ=1,nb
   do II=1,nb
      
        p1=max((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
        h2=min((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
        pre=1
        if ( p1 == (r2%mat(q)%qnph(II,2)) ) pre= -1
        p3=max((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
        h4=min((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
        if ( p3 == (r2%mat(q)%qnph(JJ,2)) ) pre= -1 * pre
        
        call sub_222_phph(r1,r2,r3,q,p1,h2,p3,h4,II,JJ,n,m,pre)

   end do 
end do

end do 


end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine  matrix_mult222(r1,r2,r3,w1,q,nh,np,nb)
  implicit none 
  
  integer :: q,nh,np,nb
  type(full_ham) :: r1,r2,r3,w1
  


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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine sub_222_hhhh(r1,r2,r3,q,h1,h2,h3,h4,II,JJ,n,m) 
  implicit none 
  
  integer :: i,n,m,a,II,JJ
  integer :: h1,h2,h3,h4,q
  type(full_ham) :: r1,r2,r3
  
  do i=1,n
     do a=n+1,m

        r3%mat(q)%Vhhhh(II,JJ) = r3%mat(q)%Vhhhh(II,JJ) +  &
             v_elem(i,h1,a,h3,r1) * v_elem(a,h2,i,h4,r2) - &
             v_elem(a,h1,i,h3,r1) * v_elem(i,h2,a,h4,r2) - &
             v_elem(i,h2,a,h3,r1) * v_elem(a,h1,i,h4,r2) - &
             v_elem(i,h1,a,h4,r1) * v_elem(a,h2,i,h3,r2) + &
             v_elem(a,h2,i,h3,r1) * v_elem(i,h1,a,h4,r2) + &
             v_elem(a,h1,i,h4,r1) * v_elem(i,h2,a,h3,r2) + &
             v_elem(i,h2,a,h4,r1) * v_elem(a,h1,i,h3,r2) - &
             v_elem(a,h2,i,h4,r1) * v_elem(i,h1,a,h3,r2)
      
     end do
  end do

end subroutine 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine sub_222_pppp(r1,r2,r3,q,p1,p2,p3,p4,II,JJ,n,m) 
  implicit none 
  
  integer :: i,n,m,a,II,JJ
  integer :: p1,p2,p3,p4,q
  type(full_ham) :: r1,r2,r3


do i=1,n
   do a=n+1,m

      r3%mat(q)%Vpppp(II,JJ) = r3%mat(q)%Vpppp(II,JJ) +  &
           v_elem(i,p1,a,p3,r1) * v_elem(a,p2,i,p4,r2) - &
           v_elem(a,p1,i,p3,r1) * v_elem(i,p2,a,p4,r2) - &
           v_elem(i,p2,a,p3,r1) * v_elem(a,p1,i,p4,r2) - &
           v_elem(i,p1,a,p4,r1) * v_elem(a,p2,i,p3,r2) + &
           v_elem(a,p2,i,p3,r1) * v_elem(i,p1,a,p4,r2) + &
           v_elem(a,p1,i,p4,r1) * v_elem(i,p2,a,p3,r2) + &
           v_elem(i,p2,a,p4,r1) * v_elem(a,p1,i,p3,r2) - &
           v_elem(a,p2,i,p4,r1) * v_elem(i,p1,a,p3,r2)
      
   end do 
end do 

end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine sub_222_pphh(r1,r2,r3,q,p1,p2,h3,h4,II,JJ,n,m) 
  implicit none 
  
  integer :: i,n,m,a,II,JJ
  integer :: p1,p2,h3,h4,q
  type(full_ham) :: r1,r2,r3


           
do i=1,n
   do a=n+1,m

      r3%mat(q)%Vpphh(II,JJ) = r3%mat(q)%Vpphh(II,JJ) +  &
           v_elem(i,p1,a,h3,r1) * v_elem(a,p2,i,h4,r2) - &
           v_elem(a,p1,i,h3,r1) * v_elem(i,p2,a,h4,r2) - &
           v_elem(i,p2,a,h3,r1) * v_elem(a,p1,i,h4,r2) - &
           v_elem(i,p1,a,h4,r1) * v_elem(a,p2,i,h3,r2) + &
           v_elem(a,p2,i,h3,r1) * v_elem(i,p1,a,h4,r2) + &
           v_elem(a,p1,i,h4,r1) * v_elem(i,p2,a,h3,r2) + &
           v_elem(i,p2,a,h4,r1) * v_elem(a,p1,i,h3,r2) - &
           v_elem(a,p2,i,h4,r1) * v_elem(i,p1,a,h3,r2)
      
   end do 
end do 
      
end subroutine 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine sub_222_ppph(r1,r2,r3,q,p1,p2,p3,h4,II,JJ,n,m,pre)
  implicit none 
  
  integer :: i,n,m,a,II,JJ,q
  integer :: p1,p2,p3,h4,pre
  type(full_ham) :: r1,r2,r3
      
do i=1,n
   do a=n+1,m

      r3%mat(q)%Vppph(II,JJ) = r3%mat(q)%Vppph(II,JJ) + pre * ( &
           v_elem(i,p1,a,p3,r1) * v_elem(a,p2,i,h4,r2) - &
           v_elem(a,p1,i,p3,r1) * v_elem(i,p2,a,h4,r2) - &
           v_elem(i,p2,a,p3,r1) * v_elem(a,p1,i,h4,r2) - &
           v_elem(i,p1,a,h4,r1) * v_elem(a,p2,i,p3,r2) + &
           v_elem(a,p2,i,p3,r1) * v_elem(i,p1,a,h4,r2) + &
           v_elem(a,p1,i,h4,r1) * v_elem(i,p2,a,p3,r2) + &
           v_elem(i,p2,a,h4,r1) * v_elem(a,p1,i,p3,r2) - &
           v_elem(a,p2,i,h4,r1) * v_elem(i,p1,a,p3,r2) )
      
   end do 
end do 

end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine sub_222_phhh(r1,r2,r3,q,p1,h2,h3,h4,II,JJ,n,m,pre)
  implicit none 
  
  integer :: i,n,m,a,II,JJ,q
  integer :: p1,h2,h3,h4,pre
  type(full_ham) :: r1,r2,r3
    
do i=1,n
   do a=n+1,m

      r3%mat(q)%Vphhh(II,JJ) = r3%mat(q)%Vphhh(II,JJ) +  pre * (&
           v_elem(i,p1,a,h3,r1) * v_elem(a,h2,i,h4,r2) - &
           v_elem(a,p1,i,h3,r1) * v_elem(i,h2,a,h4,r2) - &
           v_elem(i,h2,a,h3,r1) * v_elem(a,p1,i,h4,r2) - &
           v_elem(i,p1,a,h4,r1) * v_elem(a,h2,i,h3,r2) + &
           v_elem(a,h2,i,h3,r1) * v_elem(i,p1,a,h4,r2) + &
           v_elem(a,p1,i,h4,r1) * v_elem(i,h2,a,h3,r2) + &
           v_elem(i,h2,a,h4,r1) * v_elem(a,p1,i,h3,r2) - &
           v_elem(a,h2,i,h4,r1) * v_elem(i,p1,a,h3,r2) )
      
   end do 
end do 

end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine sub_222_phph(r1,r2,r3,q,p1,h2,p3,h4,II,JJ,n,m,pre)
  implicit none 
  
  integer :: i,n,m,a,II,JJ,q
  integer :: p1,h2,p3,h4,pre
  type(full_ham) :: r1,r2,r3
      
do i=1,n
   do a=n+1,m
      
      r3%mat(q)%Vphph(II,JJ) = r3%mat(q)%Vphph(II,JJ) + pre * ( &
           v_elem(i,p1,a,p3,r1) * v_elem(a,h2,i,h4,r2) - &
           v_elem(i,h2,a,p3,r1) * v_elem(a,p1,i,h4,r2) + &
           v_elem(i,h2,a,h4,r1) * v_elem(a,p1,i,p3,r2) - &
           v_elem(i,p1,a,h4,r1) * v_elem(a,h2,i,p3,r2) - ( &
           
           v_elem(i,p1,a,p3,r2) * v_elem(a,h2,i,h4,r1) - &
           v_elem(i,h2,a,p3,r2) * v_elem(a,p1,i,h4,r1) + &
           v_elem(i,h2,a,h4,r2) * v_elem(a,p1,i,p3,r1) - &
           v_elem(i,p1,a,h4,r2) * v_elem(a,h2,i,p3,r1) ) )
      
   end do 
end do 

end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================
subroutine xcommutator_223(r1,r2,r3) 
  implicit none 
  
  integer :: p,q,r,s,t,u,a,i,II,JJ,M3,m,n
  type(full_ham) :: r1,r2,r3
  real(8) :: sm
  
  m = r1%msp
  n = r1%nbody
  m3 = m*(m-1)*(m-2)/6
 
  do ii = 1, m3
     do jj = ii,m3
        
        p = r1%threemap(ii,1)
        q = r1%threemap(ii,2)
        r = r1%threemap(ii,3)
        s = r1%threemap(jj,1)
        t = r1%threemap(jj,2)
        u = r1%threemap(jj,3)
        
        sm = 0.d0
        do a = n+1,m
           
           sm = sm + v_elem(p,q,a,t,r1)*v_elem(a,r,s,u,r2) - &
          v_elem(p,q,a,u,r1)*v_elem(a,r,s,t,r2) - &
          v_elem(p,q,a,s,r1)*v_elem(a,r,t,u,r2) - &
          v_elem(r,q,a,t,r1)*v_elem(a,p,s,u,r2) + &
          v_elem(r,q,a,s,r1)*v_elem(a,p,t,u,r2) + &
          v_elem(r,q,a,u,r1)*v_elem(a,p,s,t,r2) - &
          v_elem(p,r,a,t,r1)*v_elem(a,q,s,u,r2) + &
          v_elem(p,r,a,s,r1)*v_elem(a,q,t,u,r2) + &
          v_elem(p,r,a,u,r1)*v_elem(a,q,s,t,r2) - &
          
          (v_elem(p,q,a,t,r2)*v_elem(a,r,s,u,r1) - &
          v_elem(p,q,a,u,r2)*v_elem(a,r,s,t,r1) - &
          v_elem(p,q,a,s,r2)*v_elem(a,r,t,u,r1) - &
          v_elem(r,q,a,t,r2)*v_elem(a,p,s,u,r1) + &
          v_elem(r,q,a,s,r2)*v_elem(a,p,t,u,r1) + &
          v_elem(r,q,a,u,r2)*v_elem(a,p,s,t,r1) - &
          v_elem(p,r,a,t,r2)*v_elem(a,q,s,u,r1) + &
          v_elem(p,r,a,s,r2)*v_elem(a,q,t,u,r1) + &
          v_elem(p,r,a,u,r2)*v_elem(a,q,s,t,r1) )
            
      end do 
      
      do i = 1, n
         
         sm = sm + v_elem(i,q,s,t,r2)*v_elem(p,r,i,u,r1) - &
              v_elem(i,q,s,u,r2)*v_elem(p,r,i,t,r1) - &
              v_elem(i,q,u,t,r2)*v_elem(p,r,i,s,r1) - &
              v_elem(i,p,s,t,r2)*v_elem(q,r,i,u,r1) + &
              v_elem(i,p,s,u,r2)*v_elem(q,r,i,t,r1) + &
              v_elem(i,p,u,t,r2)*v_elem(q,r,i,s,r1) - &
              v_elem(i,r,s,t,r2)*v_elem(p,q,i,u,r1) + &
              v_elem(i,r,s,u,r2)*v_elem(p,q,i,t,r1) + &
              v_elem(i,r,u,t,r2)*v_elem(p,q,i,s,r1) - &
              
              (v_elem(i,q,s,t,r1)*v_elem(p,r,i,u,r2) - &
              v_elem(i,q,s,u,r1)*v_elem(p,r,i,t,r2) - &
              v_elem(i,q,u,t,r1)*v_elem(p,r,i,s,r2) - &
              v_elem(i,p,s,t,r1)*v_elem(q,r,i,u,r2) + &
              v_elem(i,p,s,u,r1)*v_elem(q,r,i,t,r2) + &
              v_elem(i,p,u,t,r1)*v_elem(q,r,i,s,r2) - &
              v_elem(i,r,s,t,r1)*v_elem(p,q,i,u,r2) + &
              v_elem(i,r,s,u,r1)*v_elem(p,q,i,t,r2) + &
              v_elem(i,r,u,t,r1)*v_elem(p,q,i,s,r2) ) 
       end do 
 
       r3%v3body(ii,jj) = sm
       r3%v3body(jj,ii) = sm 
       
       end do 
    end do 

    
end subroutine xcommutator_223
!==========================================
!==========================================
subroutine xcommutator_232(r1,r2,r3) 
  implicit none
  
  integer :: p1,p2,p3,p4,h1,h2,h3,h4,i,a,ii,jj
  integer :: nh,np,nb,n,m,q,b,j,pre
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

 
     do JJ=1,nh
        do II=1,nh
         
           h1=(r2%mat(q)%qnhh(II,1))
           h2=(r2%mat(q)%qnhh(II,2))
           h3=(r2%mat(q)%qnhh(JJ,1))
           h4=(r2%mat(q)%qnhh(JJ,2))
           
           sm = 0.d0 
           do i=1,n
              do a = n+1,m
                 do b = n+1,m
                    
             sm = sm + v_elem(i,h2,a,b,r1) * w_elem(h1,a,b,h3,i,h4,r2) - &
                 v_elem(i,h1,a,b,r1) * w_elem(h2,a,b,h3,i,h4,r2)   - &
                 v_elem(a,b,i,h4,r1) *w_elem(h1,i,h2,h3,a,b,r2)  + &
                 v_elem(a,b,i,h3,r1) * w_elem(h1,i,h2,h4,a,b,r2)  
                 end do 
                 
                 do j = 1,n
                    
             sm = sm + v_elem(a,h2,i,j,r1) *  w_elem(h1,i,j,h3,a,h4,r2) - &
                 v_elem(a,h1,i,j,r1) *  w_elem(h2,i,j,h3,a,h4,r2) - &    
                 v_elem(i,j,a,h4,r1) * w_elem(h1,a,h2,h3,i,j,r2) + &
                 v_elem(i,j,a,h3,r1) * w_elem(h1,a,h2,h4,i,j,r2) 
                 end do 
                 
              end do 
           end do 
                     
           r3%mat(q)%Vhhhh(ii,jj) = r3%mat(q)%Vhhhh(ii,jj) +  sm  /2.d0
        end do 
    end do 
end if     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vpppp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ( np > 0 ) then     
 
     do JJ=1,np
        do II=1,np
           
           p1=(r2%mat(q)%qnpp(II,1))
           p2=(r2%mat(q)%qnpp(II,2))
           p3=(r2%mat(q)%qnpp(JJ,1))
           p4=(r2%mat(q)%qnpp(JJ,2))
                      
           sm = 0.d0 
           do i=1,n
              do a = n+1,m
                 do b = n+1,m
                    
             sm = sm + v_elem(i,p2,a,b,r1) * w_elem(p1,a,b,p3,i,p4,r2) - &
                 v_elem(i,p1,a,b,r1) * w_elem(p2,a,b,p3,i,p4,r2)   - &
                 v_elem(a,b,i,p4,r1) *w_elem(p1,i,p2,p3,a,b,r2)  + &
                 v_elem(a,b,i,p3,r1) * w_elem(p1,i,p2,p4,a,b,r2)  
                 end do 
                 
                 do j = 1,n
                    
             sm = sm + v_elem(a,p2,i,j,r1) *  w_elem(p1,i,j,p3,a,p4,r2) - &
                 v_elem(a,p1,i,j,r1) *  w_elem(p2,i,j,p3,a,p4,r2) - &    
                 v_elem(i,j,a,p4,r1) * w_elem(p1,a,p2,p3,i,j,r2) + &
                 v_elem(i,j,a,p3,r1) * w_elem(p1,a,p2,p4,i,j,r2) 
                 end do 
                 
              end do 
           end do
           r3%mat(q)%Vpppp(ii,jj) = r3%mat(q)%Vpppp(ii,jj) +  sm /2.d0 
        end do 
    end do 
end if     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vpphh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (np*nh > 0) then    
    
     do JJ=1,nh
        do II=1,np
           
           p1=(r2%mat(q)%qnpp(II,1))
           p2=(r2%mat(q)%qnpp(II,2))
           h3=(r2%mat(q)%qnhh(JJ,1))
           h4=(r2%mat(q)%qnhh(JJ,2))
           
           sm = 0.d0 
           do i=1,n
              do a = n+1,m
                 do b = n+1,m
                    
             sm = sm + v_elem(i,p2,a,b,r1) * w_elem(p1,a,b,h3,i,h4,r2) - &
                 v_elem(i,p1,a,b,r1) * w_elem(p2,a,b,h3,i,h4,r2)   - &
                 v_elem(a,b,i,h4,r1) *w_elem(p1,i,p2,h3,a,b,r2)  + &
                 v_elem(a,b,i,h3,r1) * w_elem(p1,i,p2,h4,a,b,r2)  
                 end do 
                 
                 do j = 1,n
                    
             sm = sm + v_elem(a,p2,i,j,r1) *  w_elem(p1,i,j,h3,a,h4,r2) - &
                 v_elem(a,p1,i,j,r1) *  w_elem(p2,i,j,h3,a,h4,r2) - &    
                 v_elem(i,j,a,h4,r1) * w_elem(p1,a,p2,h3,i,j,r2) + &
                 v_elem(i,j,a,h3,r1) * w_elem(p1,a,p2,h4,i,j,r2) 
                 end do 
                 
                 
              end do 
           end do 
           r3%mat(q)%Vpphh(ii,jj) = r3%mat(q)%Vpphh(ii,jj) +  sm/2.d0  
        end do 
    end do
    
end if 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vphhh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (nh*nb > 0 )  then 
    
     do JJ=1,nh
        do II=1,nb
           
           p1=max((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
           h2=min((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
           pre=1
           if ( p1 == (r2%mat(q)%qnph(II,2)) ) pre= -1

           h3=(r2%mat(q)%qnhh(JJ,1))
           h4=(r2%mat(q)%qnhh(JJ,2))
           
           sm = 0.d0 
           do i=1,n
              do a = n+1,m
                 do b = n+1,m
                    
             sm = sm + v_elem(i,h2,a,b,r1) * w_elem(p1,a,b,h3,i,h4,r2) - &
                 v_elem(i,p1,a,b,r1) * w_elem(h2,a,b,h3,i,h4,r2)   - &
                 v_elem(a,b,i,h4,r1) *w_elem(p1,i,h2,h3,a,b,r2)  + &
                 v_elem(a,b,i,h3,r1) * w_elem(p1,i,h2,h4,a,b,r2)  
                 end do 
                 
                 do j = 1,n
                    
             sm = sm + v_elem(a,h2,i,j,r1) *  w_elem(p1,i,j,h3,a,h4,r2) - &
                 v_elem(a,p1,i,j,r1) *  w_elem(h2,i,j,h3,a,h4,r2) - &    
                 v_elem(i,j,a,h4,r1) * w_elem(p1,a,h2,h3,i,j,r2) + &
                 v_elem(i,j,a,h3,r1) * w_elem(p1,a,h2,h4,i,j,r2) 
                 end do 
                 
                 
              end do 
           end do 
           
           r3%mat(q)%Vphhh(ii,jj) = r3%mat(q)%Vphhh(ii,jj) +  sm *pre/2.d0 
        end do 
    end do 
end if 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vppph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (np*nb > 0) then    
   
     do JJ=1,nb
        do II=1,np
           
           p1=(r2%mat(q)%qnpp(II,1))
           p2=(r2%mat(q)%qnpp(II,2))
           p3=max((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
           h4=min((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
           pre=1
           if ( p3 == (r2%mat(q)%qnph(JJ,2)) ) pre= -1
   
           sm = 0.d0 
           do i=1,n
              do a = n+1,m
                do b = n+1,m
                    
             sm = sm + v_elem(i,p2,a,b,r1) * w_elem(p1,a,b,p3,i,h4,r2) - &
                 v_elem(i,p1,a,b,r1) * w_elem(p2,a,b,p3,i,h4,r2)   - &
                 v_elem(a,b,i,h4,r1) *w_elem(p1,i,p2,p3,a,b,r2)  + &
                 v_elem(a,b,i,p3,r1) * w_elem(p1,i,p2,h4,a,b,r2)  
                 end do 
                 
                 do j = 1,n
                    
             sm = sm + v_elem(a,p2,i,j,r1) *  w_elem(p1,i,j,p3,a,h4,r2) - &
                 v_elem(a,p1,i,j,r1) *  w_elem(p2,i,j,p3,a,h4,r2) - &    
                 v_elem(i,j,a,h4,r1) * w_elem(p1,a,p2,p3,i,j,r2) + &
                 v_elem(i,j,a,p3,r1) * w_elem(p1,a,p2,h4,i,j,r2) 
                 end do 
                 
                 
              end do 
           end do 
           r3%mat(q)%Vppph(ii,jj) = r3%mat(q)%Vppph(ii,jj) +  sm  * pre / 2.d0
        end do 
    end do 
      
end if    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vphph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (nb> 0) then    
   
   
     do JJ=1,nb
        do II=1,nb
           
           
           p1=max((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
           h2=min((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
           pre=1
           if ( p1 == (r2%mat(q)%qnph(II,2)) ) pre= -1
           p3=max((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
           h4=min((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
           if ( p3 == (r2%mat(q)%qnph(JJ,2)) ) pre= -1*pre
           
           sm = 0.d0 
           do i=1,n
              do a = n+1,m
                 do b = n+1,m
                    
             sm = sm + v_elem(i,h2,a,b,r1) * w_elem(p1,a,b,p3,i,h4,r2) - &
                 v_elem(i,p1,a,b,r1) * w_elem(h2,a,b,p3,i,h4,r2)   - &
                 v_elem(a,b,i,h4,r1) *w_elem(p1,i,h2,p3,a,b,r2)  + &
                 v_elem(a,b,i,p3,r1) * w_elem(p1,i,h2,h4,a,b,r2)  
                 end do 
                 
                 do j = 1,n
                    
             sm = sm + v_elem(a,h2,i,j,r1) *  w_elem(p1,i,j,p3,a,h4,r2) - &
                 v_elem(a,p1,i,j,r1) *  w_elem(h2,i,j,p3,a,h4,r2) - &    
                 v_elem(i,j,a,h4,r1) * w_elem(p1,a,h2,p3,i,j,r2) + &
                 v_elem(i,j,a,p3,r1) * w_elem(p1,a,h2,h4,i,j,r2) 
                 end do 
                 
              end do 
           end do 
           
           r3%mat(q)%Vphph(ii,jj) = r3%mat(q)%Vphph(ii,jj) +  sm  * pre / 2.d0
        end do 
    end do 

end if     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end do 

end subroutine xcommutator_232
!================================================
!================================================
subroutine xcommutator_132(r1,r2,r3)
  implicit none 
  
  integer :: p1,p2,p3,p4,h1,h2,h3,h4,pre,q,i
  integer :: nh,np,nb,a,n,m,II,JJ
  type(full_ham) :: r1,r2,r3
  real(8) :: sm
  
  n = r1%nbody
  m = r1%msp
  
    do q=1,r2%nblock
     
     nh=r2%mat(q)%nhh
     np=r2%mat(q)%npp
     nb=r2%mat(q)%nph
           
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vhhhh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (nh > 0 ) then

 
     do JJ=1,nh
        do II=1,nh
         
           h1=(r2%mat(q)%qnhh(II,1))
           h2=(r2%mat(q)%qnhh(II,2))
           h3=(r2%mat(q)%qnhh(JJ,1))
           h4=(r2%mat(q)%qnhh(JJ,2))
           
           sm = 0.d0 
           do i=1,n
              do a = n+1,m
                
          sm = sm + r1%herm*r1%fph(a-n,i) * w_elem(h1,a,h2,h3,i,h4,r2) - &
                      r1%fph(a-n,i) * w_elem(h1,i,h2,h3,a,h4,r2) 
           
              end do
           end do 
                     
           r3%mat(q)%Vhhhh(ii,jj) = r3%mat(q)%Vhhhh(ii,jj) +  sm  
        end do 
    end do 
end if     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vpppp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ( np > 0 ) then     
 
     do JJ=1,np
        do II=1,np
           
           p1=(r2%mat(q)%qnpp(II,1))
           p2=(r2%mat(q)%qnpp(II,2))
           p3=(r2%mat(q)%qnpp(JJ,1))
           p4=(r2%mat(q)%qnpp(JJ,2))
                      
           sm = 0.d0 
           do i=1,n
              do a = n+1,m
                 
                 sm = sm + r1%herm*r1%fph(a-n,i) * w_elem(p1,a,p2,p3,i,p4,r2) - &
                      r1%fph(a-n,i) * w_elem(p1,i,p2,p3,a,p4,r2) 
              end do 
           end do
           r3%mat(q)%Vpppp(ii,jj) = r3%mat(q)%Vpppp(ii,jj) +  sm  
        end do 
    end do 
end if     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vpphh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (np*nh > 0) then    
    
     do JJ=1,nh
        do II=1,np
           
           p1=(r2%mat(q)%qnpp(II,1))
           p2=(r2%mat(q)%qnpp(II,2))
           h3=(r2%mat(q)%qnhh(JJ,1))
           h4=(r2%mat(q)%qnhh(JJ,2))
           
           sm = 0.d0 
           do i=1,n
              do a = n+1,m
                 sm = sm + r1%herm*r1%fph(a-n,i) * w_elem(p1,a,p2,h3,i,h4,r2) - &
                      r1%fph(a-n,i) * w_elem(p1,i,p2,h3,a,h4,r2) 
              end do 
           end do 
           r3%mat(q)%Vpphh(ii,jj) = r3%mat(q)%Vpphh(ii,jj) +  sm  
        end do 
    end do
    
end if 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vphhh
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (nh*nb > 0 )  then 
    
     do JJ=1,nh
        do II=1,nb
           
           p1=max((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
           h2=min((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
           pre=1
           if ( p1 == (r2%mat(q)%qnph(II,2)) ) pre= -1

           h3=(r2%mat(q)%qnhh(JJ,1))
           h4=(r2%mat(q)%qnhh(JJ,2))
           
           sm = 0.d0 
           do i=1,n
              do a = n+1,m
                 
                 sm = sm + r1%herm*r1%fph(a-n,i) * w_elem(p1,a,h2,h3,i,h4,r2) - &
                      r1%fph(a-n,i) * w_elem(p1,i,h2,h3,a,h4,r2) 
                 
              end do 
           end do 
           
           r3%mat(q)%Vphhh(ii,jj) = r3%mat(q)%Vphhh(ii,jj) +  sm *pre 
        end do 
    end do 
end if 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vppph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (np*nb > 0) then    
   
     do JJ=1,nb
        do II=1,np
           
           p1=(r2%mat(q)%qnpp(II,1))
           p2=(r2%mat(q)%qnpp(II,2))
           p3=max((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
           h4=min((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
           pre=1
           if ( p3 == (r2%mat(q)%qnph(JJ,2)) ) pre= -1
   
           sm = 0.d0 
           do i=1,n
              do a = n+1,m
                 
                 sm = sm + r1%herm*r1%fph(a-n,i) * w_elem(p1,a,p2,p3,i,h4,r2) - &
                      r1%fph(a-n,i) * w_elem(p1,i,p2,p3,a,h4,r2) 
              end do 
           end do 
           r3%mat(q)%Vppph(ii,jj) = r3%mat(q)%Vppph(ii,jj) +  sm  * pre
        end do 
    end do 
      
end if    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Vphph
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (nb> 0) then    
   
   
     do JJ=1,nb
        do II=1,nb
           
           
           p1=max((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
           h2=min((r2%mat(q)%qnph(II,1)),(r2%mat(q)%qnph(II,2)))
           pre=1
           if ( p1 == (r2%mat(q)%qnph(II,2)) ) pre= -1
           p3=max((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
           h4=min((r2%mat(q)%qnph(JJ,1)),(r2%mat(q)%qnph(JJ,2)))
           if ( p3 == (r2%mat(q)%qnph(JJ,2)) ) pre= -1*pre
           
           sm = 0.d0 
           do i=1,n
              do a = n+1,m
                
                 sm = sm + r1%herm*r1%fph(a-n,i) * w_elem(p1,a,h2,p3,i,h4,r2) - &
                      r1%fph(a-n,i) * w_elem(p1,i,h2,p3,a,h4,r2) 
              end do 
           end do 
           
           r3%mat(q)%Vphph(ii,jj) = r3%mat(q)%Vphph(ii,jj) +  sm  * pre
          
        end do 
    end do 

end if     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end do 
    
end subroutine
!=======================================
!=======================================
subroutine xcommutator_231(r1,r2,r3)
  implicit none 
  
  integer :: p1,p2,h1,h2,i,j,a,b,n,m
  type(full_ham) :: r1,r2,r3
  real(8) ::  sm

  n = r1%nbody
  m = r1%msp
  ! fhh
  
  do h1 = 1, n
     do h2 = h1, n
        
        sm = 0.d0 
        do i = 1,n
           do j = i+1,n
              do a = n+1,m
                 do b = a+1 ,m
                    
             sm = sm + v_elem(i,j,a,b,r1) * w_elem(h1,a,b,h2,i,j,r2) - &
                  v_elem(a,b,i,j,r1) * w_elem(h1,i,j,h2,a,b,r2) 
                 
                end do
             end do 
          end do 
       end do 
       r3%fhh(h1,h2) = r3%fhh(h1,h2) + sm 
       r3%fhh(h2,h1) = r3%fhh(h2,h1) + sm 
    end do
  end do 

  ! fpp
  
  do p1 =  n+1,m
     do p2 = p1 , m
        
        sm = 0.d0 
        do i = 1,n
           do j = i+1,n
              do a = n+1,m
                 do b = a+1 ,m
                    
             sm = sm + v_elem(i,j,a,b,r1) * w_elem(p1,a,b,p2,i,j,r2) - &
                  v_elem(a,b,i,j,r1) * w_elem(p1,i,j,p2,a,b,r2) 
                 
                end do
             end do 
          end do 
       end do 
       r3%fpp(p1-n,p2-n) = r3%fpp(p1-n,p2-n) + sm 
       r3%fpp(p2-n,p1-n) = r3%fpp(p2-n,p1-n) + sm 
    end do
  end do 

   ! fph
  
  do p1 = n+1,m
     do h2 = 1 , n
        
        sm = 0.d0 
        do i = 1,n
           do j = i+1,n
              do a = n+1,m
                 do b = a+1 ,m
                    
             sm = sm + v_elem(i,j,a,b,r1) * w_elem(p1,a,b,h2,i,j,r2) - &
                  v_elem(a,b,i,j,r1) * w_elem(p1,i,j,h2,a,b,r2) 
                 
                end do
             end do 
          end do 
       end do 
       r3%fph(p1-n,h2) = r3%fph(p1-n,h2) + sm 
       
    end do
  end do 
  
  
end subroutine  
!======================================
!======================================
subroutine xcommutator_233(r1,r2,r3)
  implicit none 
  
  integer :: p,q,r,s,t,u,a,b,i,j,m,n,m3,II,JJ
  type(full_ham) :: r1,r2,r3
  real(8) :: sm,sm2
  
  m = r1%msp
  n = r1%nbody
  m3 = m*(m-1)*(m-2)/6
  
  do II = 1,m3
     do JJ = II,m3
           
        p = r1%threemap(II,1) 
        q = r1%threemap(II,2) 
        r = r1%threemap(II,3) 
        s = r1%threemap(JJ,1) 
        t = r1%threemap(JJ,2) 
        u = r1%threemap(JJ,3) 
        
        sm = 0.d0
        do i = 1,n
           do a = n+1,m
              
              sm = sm + v_elem(p,i,s,a,r1) * w_elem(q,a,r,t,i,u,r2) - &
                   v_elem(q,i,s,a,r1) * w_elem(p,a,r,t,i,u,r2) - &
                   v_elem(r,i,s,a,r1) * w_elem(q,a,p,t,i,u,r2) - &
                   v_elem(p,i,t,a,r1) * w_elem(q,a,r,s,i,u,r2) + &
                   v_elem(q,i,t,a,r1) * w_elem(p,a,r,s,i,u,r2) + &
                   v_elem(r,i,t,a,r1) * w_elem(q,a,p,s,i,u,r2) - &
                   v_elem(p,i,u,a,r1) * w_elem(q,a,r,t,i,s,r2) + &
                   v_elem(q,i,u,a,r1) * w_elem(p,a,r,t,i,s,r2) + &
                   v_elem(r,i,u,a,r1) * w_elem(q,a,p,t,i,s,r2) - &
                   
                   ( v_elem(q,a,t,i,r1) * w_elem(p,i,r,s,a,u,r2) - &
                   v_elem(p,a,t,i,r1) * w_elem(q,i,r,s,a,u,r2) - &
                   v_elem(r,a,t,i,r1) * w_elem(p,i,q,s,a,u,r2) - &
                   v_elem(q,a,s,i,r1) * w_elem(p,i,r,t,a,u,r2) + &
                   v_elem(p,a,s,i,r1) * w_elem(q,i,r,t,a,u,r2) + &
                   v_elem(r,a,s,i,r1) * w_elem(p,i,q,t,a,u,r2) - &
                   v_elem(q,a,u,i,r1) * w_elem(p,i,r,s,a,t,r2) + &
                   v_elem(p,a,u,i,r1) * w_elem(q,i,r,s,a,t,r2) + &
                   v_elem(r,a,u,i,r1) * w_elem(p,i,q,s,a,t,r2) ) 
              
           end do 
        end do 
        
        sm2 = 0.d0
        do a=n+1,m
           do b = n+1,m
              
              sm2 = sm2 +  v_elem(p,q,a,b,r1) * w_elem(a,b,r,s,t,u,r2) - &
                   v_elem(p,r,a,b,r1) * w_elem(a,b,q,s,t,u,r2) - &
                   v_elem(r,q,a,b,r1) * w_elem(a,b,p,s,t,u,r2)  - &
                   
                   (v_elem(a,b,s,t,r1) * w_elem(p,q,r,a,b,u,r2) - &
                   v_elem(a,b,s,u,r1) * w_elem(p,q,r,a,b,t,r2) - &
                   v_elem(a,b,u,t,r1) * w_elem(p,q,r,a,b,s,r2) ) 
           end do 
        end do 
        
        do i =1,n
           do j = 1,n
              
              sm2 = sm2 + v_elem(i,j,s,t,r1) * w_elem(p,q,r,i,j,u,r2) - &
                   v_elem(i,j,s,u,r1) * w_elem(p,q,r,i,j,t,r2) -&
                   v_elem(i,j,u,t,r1) * w_elem(p,q,r,i,j,s,r2) - &
                   
                   ( v_elem(p,q,i,j,r1) * w_elem(i,j,r,s,t,u,r2) - &
                    v_elem(p,r,i,j,r1) * w_elem(i,j,q,s,t,u,r2) - &
                    v_elem(r,q,i,j,r1) * w_elem(i,j,p,s,t,u,r2) )
           end do 
        end do 
        
        r3%V3body(II,JJ) = r3%V3body(II,JJ) + sm + 0.5d0 * sm2
        r3%V3body(JJ,II) = r3%V3body(II,JJ)
        
      end do
    end do 
                    
end subroutine 
!==========================
subroutine xcommutator_133(r1,r2,r3) 
  implicit none 
  
  integer :: p,q,r,s,t,u,a,i,II,JJ,m,n,m3
  type(full_ham) :: r1,r2,r3
  real(8) :: sm 
  
   
  m = r1%msp
  n = r1%nbody
  m3 = m*(m-1)*(m-2)/6
  
  do II = 1,m3
     do JJ = II,m3
           
        p = r1%threemap(II,1) 
        q = r1%threemap(II,2) 
        r = r1%threemap(II,3) 
        s = r1%threemap(JJ,1) 
        t = r1%threemap(JJ,2) 
        u = r1%threemap(JJ,3) 
        
        sm = 0.d0 
        do a=n+1,m
           sm = sm + f_elem(q,a,r1) * w_elem(p,a,r,s,t,u,r2) - &
                f_elem(p,a,r1) * w_elem(q,a,r,s,t,u,r2) - &
                f_elem(r,a,r1) * w_elem(p,a,q,s,t,u,r2) - &
                ( f_elem(a,t,r1) * w_elem(p,q,r,s,a,u,r2) - &
                 f_elem(a,s,r1) * w_elem(p,q,r,t,a,u,r2) - &
                  f_elem(a,u,r1) * w_elem(p,q,r,s,a,t,r2) )
        end do 
        
        do i=1,n
           sm = sm + f_elem(q,i,r1) * w_elem(p,i,r,s,t,u,r2) - &
                f_elem(p,i,r1) * w_elem(q,i,r,s,t,u,r2) - &
                f_elem(r,i,r1) * w_elem(p,i,q,s,t,u,r2) - &
                ( f_elem(i,t,r1) * w_elem(p,q,r,s,i,u,r2) - &
                f_elem(i,s,r1) * w_elem(p,q,r,t,i,u,r2) - &
                f_elem(i,u,r1) * w_elem(p,q,r,s,i,t,r2) ) 
        end do 
        
        r3%V3body(II,JJ) =r3%V3body(II,JJ) + sm
        r3%V3body(JJ,II) = r3%V3body(II,JJ)
     end do 
  end do
end subroutine
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
           
end module

