implicit none 

integer,parameter :: R = 4
real(8),parameter :: hw = 1.0
integer :: E,Ml,Ms,n,states(R*(R+1),3),q,Mtot   
integer :: i,j,k,l
real(8) :: v_int

Mtot = R*(R+1)

!!! organize sp states
q = 1

do E = 1,R
   do Ml = -(E-1),E-1,2
      n = (E - 1 - abs(Ml))/2  
      do Ms = -1,1,2 
        
         states(q,1) = n
         states(q,2) = Ml
         states(q,3) = Ms 
         
         q = q + 1
      
      end do 
   end do 
end do 


open(unit=31,file='mat_elements.dat') 

do i=1,Mtot
   do j = i+1,Mtot
      do k = 1,Mtot
         do l = k+1,Mtot 
          
write(31,'(4(I5),f15.7)') i,j,k,l, v_int( states(i,1:3) , states(j,1:3),   &
                  states(l,1:3) , states(k,1:3) , hw ) - &
                     v_int( states(i,1:3) , states(j,1:3),   &
                  states(k,1:3) , states(l,1:3) , hw )  
          end do 
      end do 
   end do 
end do 
 
close(31)

end program 
!================================================
!================================================
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
!======================================
!======================================
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
