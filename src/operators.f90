module operators
use ME_general
implicit none 

contains

subroutine build_TDAmat(H,TDA)
  implicit none 
  
  type(full_ham) :: H
  type(full_sp_block_mat) :: TDA
  integer,allocatable,dimension(:,:) :: flb
  integer :: i,j,a,b,q,pre,II,JJ,kron_del,n,m,eF,emax
  integer :: lmax,li,la,si,sa,d,Ml,Ms
  
  n=H%nbody
  m=H%msp
  eF=nint((Sqrt(1.+4*n)-1)/2.)
  emax=nint((Sqrt(1.+4*m)-1)/2.)
  
  lmax=emax+ef-2
  TDA%blocks=3*(2*lmax+1)
  allocate(TDA%blkM(TDA%blocks))
  allocate(TDA%map(TDA%blocks))
  allocate(flb(n*(m-n),3)) 
 
  q=1
  !! name the blocks
  do Ms = -2,2,2
     do Ml = -lmax,lmax
        
        TDA%blkM(q)%lmda(1) = Ml
        TDA%blkM(q)%lmda(2) = Ms
        q=q+1
   
        
     end do 
  end do 
 
  TDA%map = 0
  !find the dimension
  do i=1,n
     do a=n+1,12 ! max out at 12 (valence space) 
              
         li = H%states(H%eh(i),2)
         si = H%states(H%eh(i),3)
         la = H%states(H%ep(a-n),2)
         sa = H%states(H%ep(a-n),3)
              
         Ml = la-li
         Ms = sa-si 
         
         q=Ml+lmax+1+(2*lmax+1)*(Ms+2)/2 
!!this is the ordering established in the previous loop
         
         TDA%map(q)=TDA%map(q)+1
     end do 
  end do        
  
!! allocate the blocks
  do q=1,TDA%blocks
     
     d = TDA%map(q)
     allocate(TDA%blkM(q)%matrix(d,d))
     allocate(TDA%blkM(q)%Eigval(d))
     allocate(TDA%blkM(q)%extra(10*d)) 
     allocate(TDA%blkM(q)%labels(d,2))
     TDA%blkM(q)%labels=0
  end do 
  
!! fill it with states
  do i=1,n
     do a=n+1,12
              
         li = H%states(H%eh(i),2)
         si = H%states(H%eh(i),3)
         la = H%states(H%ep(a-n),2)
         sa = H%states(H%ep(a-n),3)
              
         Ml = la-li
         Ms = sa-si 
         
         q=Ml+lmax+1+(2*lmax+1)*(Ms+2)/2          
        
         do d=1,TDA%map(q)
            if (TDA%blkM(q)%labels(d,1) == 0) then 
               TDA%blkM(q)%labels(d,1) = i
               TDA%blkM(q)%labels(d,2) = a 
               exit
            end if 
         end do 
         
     end do 
  end do             

end subroutine 
!==============================================================
!==============================================================
subroutine calc_TDA(H,TDA)
  implicit none 
  
  type(full_ham) :: H
  type(full_sp_block_mat) :: TDA
  integer :: n,q,info,i,j,a,b,II,JJ,kron_del
  
  n=H%nbody
  do q = 1,TDA%blocks
     if (TDA%map(q) > 0) then
     
         do II=1,TDA%map(q) 
            i=TDA%blkM(q)%labels(II,1) 
            a=TDA%blkM(q)%labels(II,2) 
              
            do JJ=1,TDA%map(q)
               
               
               j=TDA%blkM(q)%labels(JJ,1) 
               b=TDA%blkM(q)%labels(JJ,2) 
               
               TDA%blkM(q)%matrix(II,JJ ) = &
                   kron_del(i,j)*H%fpp(a-n,b-n) - &
                   kron_del(a,b)*H%fhh(j,i) + &
                   v_elem(a,j,i,b,H)
            end do 
        end do 
     end if
  end do 

end subroutine
!==========================================
!==========================================
subroutine calc_TDA_single_block(H,TDA,q)
  implicit none 
  
  type(full_ham) :: H
  type(full_sp_block_mat) :: TDA
  integer :: n,q,info,i,j,a,b,II,JJ,kron_del,ML,MS
  
  n=H%nbody
  
     if (TDA%map(q) > 0) then
     
         do II=1,TDA%map(q) 
            do JJ=1,TDA%map(q)
               
               i=TDA%blkM(q)%labels(II,1) 
               a=TDA%blkM(q)%labels(II,2) 
               j=TDA%blkM(q)%labels(JJ,1) 
               b=TDA%blkM(q)%labels(JJ,2) 
               
               TDA%blkM(q)%matrix(II,JJ ) = &
                   kron_del(i,j)*H%fpp(a-n,b-n) - &
                   kron_del(a,b)*H%fhh(j,i) + &
                   v_elem(a,j,i,b,H)
            end do 
        end do 
     end if

end subroutine
!==============================================================
!==============================================================
subroutine fullmat_TDA(H,V)
  implicit none 
  
  type(full_ham) :: H
  integer :: nb,info,m,n,i,j,q,b,p,a,kron_del
  real(8),allocatable,dimension(:,:) :: F
  real(8),allocatable,dimension(:) :: w
  real(8),dimension(:) :: V

  n=H%nbody
  m=H%msp
  
  nb=n*(m-n)

  allocate(f(nb,nb))
  allocate(w(10*nb))

  q=1
  do i=1,n
     do a=n+1,m
        
        p=1
        do j=1,n
           do b=n+1,m
              
              F(q,p) = v_elem(a,j,i,b,H) + &
                   kron_del(i,j) * H%fpp(a-n,b-n) - &
                   kron_del(a,b) * H%fhh(j,i) 
              p=p+1
           end do 
       end do 
       q=q+1
     end do 
 end do 
 
 call dsyev('V','U',nb,F,nb,V,w,10*nb,info)
 
end subroutine  

!==============================================
subroutine density(r1,den,x,hw,coefs) 
  !!! this calculates the density operator
  !!! it should only work for s=0
  !!! it only calculates for a single value of x (position) 
  implicit none 
  
  real(8),parameter :: al = 1.d0, bet=0.d0
  type(full_ham) :: r1,den
  integer :: i,j,k,n,m
  real(8) :: x,psi,hw,dn
  real(8),dimension(r1%msp,r1%msp) :: coefs,rho,im
  
  n=r1%nbody
  m=r1%msp

  do i=1,m
     do j=1,m
        
        if ( (r1%states(i,3) == r1%states(j,3) ) .and. &
             (r1%states(i,2) == r1%states(j,2) ) ) then 

              rho(i,j) = psi(r1%states(i,1),r1%states(i,2),hw,x) * &
                   psi(r1%states(j,1),r1%states(j,2),hw,x) 
          
        else
          
           rho(i,j) = 0.d0
        end if
    end do 
  end do 
 

  call dgemm('N','N',m,m,m,al,rho,m,coefs,m,bet,im,m)
  call dgemm('T','N',m,m,m,al,coefs,m,im,m,bet,rho,m)
 

  do i=1,n
    do j=1,n
       
       den%fhh(i,j) = rho( r1%eh(i) , r1%eh(j) )
    
    end do 
  end do 
 
  do i=1,m-n
    do j=1,m-n
       
       den%fpp(i,j) = rho( r1%ep(i) , r1%ep(j) )
       
    end do 
  end do 
 
  do i=1,m-n
    do j=1,n
       
       den%fph(i,j) = rho( r1%ep(i) , r1%eh(j) )
       
    end do 
  end do 
 
 dn = 0.d0
 
 do k=1,n
          
          dn = dn + rho(r1%eh(k),r1%eh(k))
          
 end do 
 
 den%E0 = dn
          
 
end subroutine 
!==============================================
  
end module
