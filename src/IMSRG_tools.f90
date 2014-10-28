module IMSRG_tools
  use ME_general
  use commutators
  implicit none 
    
  contains

subroutine build_wegner(gen, H)
  implicit none 
  
  type(full_ham) ::  gen, H,HD,wkspc,wks2
  integer :: i,j,a,b,n,m,q,nh,np,nb,pre,kron_del
  integer :: p1,p2,p3,p4,h1,h2,h3,h4,II,JJ
  
  gen%herm = -1
  n = H%nbody
  m = H%msp 

  call allocate_everything(H,HD)
  call allocate_everything(H,wkspc)
  call allocate_everything(H,wks2) 
!! DEFINE H-DIAGONAL !!================

  do i=1,n
  HD%Fhh(i,i) = H%Fhh(i,i)
  end do 
  
  do i=1,m-n
  HD%Fpp(i,i) = H%Fpp(i,i)
  end do 

  do q=1,H%nblock

     do i=1,H%mat(q)%nhh
        
        HD%mat(q)%Vhhhh(i,i)=H%mat(q)%Vhhhh(i,i)
        
     end do 

     do i=1,H%mat(q)%npp
        
        HD%mat(q)%Vpppp(i,i)=H%mat(q)%Vpppp(i,i)
        
     end do 
     
     do i=1,H%mat(q)%nph
        
        HD%mat(q)%Vphph(i,i)=H%mat(q)%Vphph(i,i)
        
     end do 
    
  end do 

  call xcommutator_111(HD,H,GEN) 
  call xcommutator_121(HD,H,GEN) 
  call xcommutator_221(HD,H,GEN,wkspc,wks2) 

  call xcommutator_122(HD,H,GEN)
  call xcommutator_222(HD,H,GEN,wkspc) 

end subroutine 
!==============================================
!==============================================
subroutine build_white_mixed( gen, rec ) 
  ! removes phhh 
  implicit none
 
  
  integer :: i,j,k,q,q1,q2,a,n,m,p(2),h(2),be(2),he(2),pe(2),pp,row(2,3)
  type(full_ham) ::  gen , rec 
  real(8) :: efs, Aph,mx,rb,denom
  real(8), parameter :: dcut = 0.07
  
  gen%herm = -1
  n=rec%nbody
  m=rec%Msp 
  
  
  do a=n+1,m
     do i=1,n
       

        gen%fph(a-n,i) = rec%fph(a-n,i) / &
             ( rec%fpp(a-n,a-n) - rec%fhh(i,i) - v_elem(a,i,a,i,rec) ) 
      
        
     end do
  end do 

 
  !!! two body term !!!

  mx = 0.
  do q=1,rec%nblock
     
     if ( rec%mat(q)%nph*rec%mat(q)%nhh > 0 ) then 

     do i=1,gen%mat(q)%nph
        do j=1,gen%mat(q)%nhh
           
           h(:)=rec%mat(q)%qnhh(j,:)
          
           be(1)=max(rec%stoe(rec%mat(q)%qnph(I,1)),rec%stoe(rec%mat(q)%qnph(I,2)))
           be(2)=min(rec%stoe(rec%mat(q)%qnph(I,1)),rec%stoe(rec%mat(q)%qnph(I,2)))
           he(1)=rec%stoe(h(1))
           he(2)=rec%stoe(h(2))
           
           efs=rec%fpp(be(1)-n,be(1)-n) 
           efs=efs+rec%fhh(be(2),be(2)) 
           efs=efs-rec%fhh(he(1),he(1)) 
           efs=efs-rec%fhh(he(2),he(2)) 
           
           Aph=v_elem(he(1),he(2),he(1),he(2),rec) - &
               v_elem(he(1),be(1),he(1),be(1),rec) - &
               v_elem(he(2),be(1),he(2),be(1),rec) 
           
           denom = efs + Aph 
           if (abs(denom) > dcut) then
              gen%mat(q)%Vphhh(i,j) = rec%mat(q)%Vphhh(i,j)/ denom
           end if 
                
        end do 
     end do 
     
   end if 
 
   
     if ( rec%mat(q)%nph*rec%mat(q)%npp > 0 ) then 
         
     do i=1,gen%mat(q)%npp
        do j=1,gen%mat(q)%nph
           
           p(:)=rec%mat(q)%qnpp(i,:)
       
           be(1)=max(rec%stoe(rec%mat(q)%qnph(J,1)),rec%stoe(rec%mat(q)%qnph(J,2)))
           be(2)=min(rec%stoe(rec%mat(q)%qnph(J,1)),rec%stoe(rec%mat(q)%qnph(J,2)))
           pe(1)=rec%stoe(p(1))
           pe(2)=rec%stoe(p(2))
          
           efs=-rec%fpp(be(1)-n,be(1)-n) 
           efs=efs-rec%fhh(be(2),be(2)) 
           efs=efs+rec%fpp(pe(1)-n,pe(1)-n) 
           efs=efs+rec%fpp(pe(2)-n,pe(2)-n) 
           
          Aph=v_elem(pe(1),pe(2),pe(1),pe(2),rec) - &
               v_elem(pe(1),be(2),pe(1),be(2),rec) - &
               v_elem(pe(2),be(2),pe(2),be(2),rec) 
           
           denom = efs+Aph
           if (abs(denom) > dcut) then
              gen%mat(q)%Vppph(i,j) = rec%mat(q)%Vppph(i,j)/denom
           end if  
        
           
        end do 
      end do 
     
 
     end if

  end do 
  
end subroutine
!==============================================
!==============================================
subroutine build_imaginary_time( gen, rec ) 
  ! removes phhh 
  implicit none
 
  
  integer :: i,j,k,q,q1,q2,a,b,n,m,p(2),h(2),be(2),he(2),pe(2),pp,row(2,3)
  type(full_ham) ::  gen , rec 
  real(8) :: efs, Aph,mx,rb,denom

  gen%herm = -1
  n=rec%nbody
  m=rec%Msp 
  
  
  do a=n+1,m
     do i=1,n
       
        denom = rec%fpp(a-n,a-n) - rec%fhh(i,i) - v_elem(a,i,a,i,rec)
        gen%fph(a-n,i) = rec%fph(a-n,i) * &
      sign(1.d0 ,denom ) * abs(denom)**.0001 
      
        
     end do
  end do 

  do a = n+1,12
     do b = 12,m
        
        denom = rec%fpp(a-n,a-n) - rec%fpp(b-n,b-n) -v_elem(a,b,a,b,rec) 
        gen%fpp(a-n,b-n) = rec%fpp(a-n,b-n) * sign(1.d0,denom) * abs(denom)**.0001
        gen%fpp(b-n,a-n) = -1*gen%fpp(a-n,b-n) 
     end do 
  end do 

  !!! two body term !!!

  mx = 0.
  do q=1,rec%nblock
     
     if ( rec%mat(q)%nph*rec%mat(q)%nhh > 0 ) then 

     do i=1,gen%mat(q)%nph
        do j=1,gen%mat(q)%nhh
           
           h(:)=rec%mat(q)%qnhh(j,:)
          
           be(1)=max(rec%stoe(rec%mat(q)%qnph(I,1)),rec%stoe(rec%mat(q)%qnph(I,2)))
           be(2)=min(rec%stoe(rec%mat(q)%qnph(I,1)),rec%stoe(rec%mat(q)%qnph(I,2)))
           he(1)=rec%stoe(h(1))
           he(2)=rec%stoe(h(2))
           
           efs=rec%fpp(be(1)-n,be(1)-n) 
           efs=efs+rec%fhh(be(2),be(2)) 
           efs=efs-rec%fhh(he(1),he(1)) 
           efs=efs-rec%fhh(he(2),he(2)) 
           
           Aph=v_elem(he(1),he(2),he(1),he(2),rec) - &
               v_elem(he(1),be(1),he(1),be(1),rec) - &
               v_elem(he(2),be(1),he(2),be(1),rec) 
           
           denom = efs + Aph 
          
           gen%mat(q)%Vphhh(i,j) = rec%mat(q)%Vphhh(i,j)* sign(1.d0, denom)*abs(denom)**.0001
          
                
        end do 
     end do 
     
   end if 
 
   
     if ( rec%mat(q)%nph*rec%mat(q)%npp > 0 ) then 
         
     do i=1,gen%mat(q)%npp
        do j=1,gen%mat(q)%nph
           
           p(:)=rec%mat(q)%qnpp(i,:)
       
           be(1)=max(rec%stoe(rec%mat(q)%qnph(J,1)),rec%stoe(rec%mat(q)%qnph(J,2)))
           be(2)=min(rec%stoe(rec%mat(q)%qnph(J,1)),rec%stoe(rec%mat(q)%qnph(J,2)))
           if (be(1) > 12) cycle
           
           pe(1)=rec%stoe(p(1))
           pe(2)=rec%stoe(p(2))
          
           efs=-rec%fpp(be(1)-n,be(1)-n) 
           efs=efs-rec%fhh(be(2),be(2)) 
           efs=efs+rec%fpp(pe(1)-n,pe(1)-n) 
           efs=efs+rec%fpp(pe(2)-n,pe(2)-n) 
           
          Aph=v_elem(pe(1),pe(2),pe(1),pe(2),rec) - &
               v_elem(pe(1),be(2),pe(1),be(2),rec) - &
               v_elem(pe(2),be(2),pe(2),be(2),rec) 
           
           denom = efs+Aph
           
           gen%mat(q)%Vppph(i,j) = rec%mat(q)%Vppph(i,j)*sign(1.d0,denom)*abs(denom)**.0001
           
        end do 
      end do 
           
     end if

     
     if ( rec%mat(q)%nph > 0 ) then 
         
     do i=1,gen%mat(q)%nph
        pe(1)=max(rec%stoe(rec%mat(q)%qnph(I,1)),rec%stoe(rec%mat(q)%qnph(I,2)))
        pe(2)=min(rec%stoe(rec%mat(q)%qnph(I,1)),rec%stoe(rec%mat(q)%qnph(I,2)))
        if (pe(1) > 12) cycle
        
        do j=i+1,gen%mat(q)%nph
              
           be(1)=max(rec%stoe(rec%mat(q)%qnph(J,1)),rec%stoe(rec%mat(q)%qnph(J,2)))
           be(2)=min(rec%stoe(rec%mat(q)%qnph(J,1)),rec%stoe(rec%mat(q)%qnph(J,2)))
           
           if (be(1) < 13) cycle
           
           
           efs=-rec%fpp(be(1)-n,be(1)-n) 
           efs=efs-rec%fhh(be(2),be(2)) 
           efs=efs+rec%fpp(pe(1)-n,pe(1)-n) 
           efs=efs+rec%fhh(pe(2),pe(2)) 
                     
           denom = efs
           
           gen%mat(q)%Vphph(i,j) = rec%mat(q)%Vphph(i,j)*sign(1.d0,denom)*abs(denom)**.0001
           gen%mat(q)%Vphph(j,i) = -1*gen%mat(q)%Vphph(i,j) 
           
        end do 
      end do 
     
 
     end if
  end do 
  
end subroutine
!==============================================
!==============================================
subroutine xbuild_imaginary_time( gen, rec ) 
  ! removes phhh 
  implicit none
 
  
  integer :: i,j,k,q,q1,q2,a,n,m,p(2),h(2),be(2),he(2),pe(2),pp,row(2,3)
  type(full_ham) ::  gen , rec 
  real(8) :: efs, Aph,mx,rb,denom
  
  gen%herm = -1
  n=rec%nbody
  m=rec%Msp 
  
  gen%fph = rec%fph
             
  !!! two body term !!!

  do q=1,rec%nblock
     
     if ( rec%mat(q)%nph*rec%mat(q)%nhh > 0 ) then 
         gen%mat(q)%Vphhh = rec%mat(q)%Vphhh
     end if 
 
   
     if ( rec%mat(q)%nph*rec%mat(q)%npp > 0 ) then 
        gen%mat(q)%Vppph = rec%mat(q)%Vppph
     end if

  end do 
  
end subroutine
!==============================================
!==============================================
subroutine build_white( gen, rec ) 
  ! removes pphh 
  implicit none 

  integer :: i,j,k,q,q1,q2,a,n,m,p(2),h(2),pe(2),he(2),pp,row(2,3)
  type(full_ham) ::  gen , rec 
  real(8) :: efs, Aph,mx

  n=rec%nbody
  m=rec%Msp
  !!! 1 body term !!! 
  !!! only ph term exists

  do a=n+1,m
     do i=1,n
        
        gen%fph(a-n,i) = rec%fph(a-n,i) / &
             ( rec%fpp(a-n,a-n) - rec%fhh(i,i) - v_elem(a,i,a,i,rec) ) 
     end do
  end do 
  
  !!! two body term !!!

  do q=1,rec%nblock
     
     if ( rec%mat(q)%npp*rec%mat(q)%nhh == 0 ) cycle 

     do i=1,gen%mat(q)%npp
        do j=1,gen%mat(q)%nhh
           
           p(:)=rec%mat(q)%qnpp(i,:)
           h(:)=rec%mat(q)%qnhh(j,:)
           
           pe(1)=rec%stoe(p(1))
           pe(2)=rec%stoe(p(2))
           he(1)=rec%stoe(h(1))
           he(2)=rec%stoe(h(2))
           
           efs=rec%fpp(pe(1)-n,pe(1)-n) 
           efs=efs+rec%fpp(pe(2)-n,pe(2)-n) 
           efs=efs-rec%fhh(he(1),he(1)) 
           efs=efs-rec%fhh(he(2),he(2)) 
           
           Aph=v_elem(pe(1),pe(2),pe(1),pe(2),rec) + &
             v_elem(he(1),he(2),he(1),he(2),rec) - &
             v_elem(pe(1),he(1),pe(1),he(1),rec) - &
             v_elem(pe(1),he(2),pe(1),he(2),rec) - &
             v_elem(pe(2),he(1),pe(2),he(1),rec) - &
             v_elem(pe(2),he(2),pe(2),he(2),rec) 
           
           gen%mat(q)%Vpphh(i,j) = rec%mat(q)%Vpphh(i,j)/(efs + Aph ) 
           
        end do 
     end do 
  
  end do 
 
end subroutine
!========================================
!========================================
subroutine allocate_everything(rec,r2) 
  implicit none 
  
  type(full_ham) :: rec,r2
  integer :: i,n,m,m3,g
  
  n=rec%nbody
  m=rec%Msp
  m3 = m*(m-1)*(m-2)/6
  g = m*(m-1)/2
  r2%herm=1
  r2%nbody=rec%nbody
  r2%Msp=m
  r2%neq = rec%neq 
  allocate(r2%eh(n)) 
  allocate(r2%ep(m-n))
  allocate(r2%stoe(m)) 
  
  r2%stoe=rec%stoe
  r2%eh=rec%eh
  r2%ep=rec%ep
  
  allocate(r2%i_array(g,3))
  
  r2%i_array=rec%i_array

  allocate(r2%fph(m-n,n))
  
  r2%nblock=rec%nblock
  allocate( r2%mat( r2%nblock ) )

  allocate(r2%fpp(m-n,m-n))
  allocate(r2%fhh(n,n)) 
  r2%fpp=0.d0;r2%fph=0.d0;r2%fhh=0.d0
  
  
  do i=1,r2%nblock
     r2%mat(i)%nhh=rec%mat(i)%nhh
     r2%mat(i)%npp=rec%mat(i)%npp
     r2%mat(i)%nph=rec%mat(i)%nph
     r2%mat(i)%lam=rec%mat(i)%lam
     
     allocate(r2%mat(i)%qnhh(rec%mat(i)%nhh,2))
     allocate(r2%mat(i)%qnpp(rec%mat(i)%npp,2))
     allocate(r2%mat(i)%qnph(rec%mat(i)%nph,2))
     allocate(r2%mat(i)%vhhhh(rec%mat(i)%nhh,rec%mat(i)%nhh))
     allocate(r2%mat(i)%vpppp(rec%mat(i)%npp,rec%mat(i)%npp))
     allocate(r2%mat(i)%vphhh(rec%mat(i)%nph,rec%mat(i)%nhh))
     allocate(r2%mat(i)%vppph(rec%mat(i)%npp,rec%mat(i)%nph))
     allocate(r2%mat(i)%vpphh(rec%mat(i)%npp,rec%mat(i)%nhh))
     allocate(r2%mat(i)%vphph(rec%mat(i)%nph,rec%mat(i)%nph))

     r2%mat(i)%qnhh = rec%mat(i)%qnhh
     r2%mat(i)%qnpp = rec%mat(i)%qnpp
     r2%mat(i)%qnph = rec%mat(i)%qnph

     r2%mat(i)%Vhhhh=0.d0;r2%mat(i)%Vpppp=0.d0;r2%mat(i)%Vpphh=0.d0
     r2%mat(i)%Vppph=0.d0;r2%mat(i)%Vphhh=0.d0;r2%mat(i)%Vphph=0.d0
     !! thats a lot of allocation
  end do 
  r2%IMSRG3 = .false.
  
  if (rec%IMSRG3) then 
     r2%IMSRG3 = .true.
     allocate(r2%V3body(m3,m3))
     allocate(r2%threemap(m3,3))
     r2%threemap = rec%threemap
  end if 
  
end subroutine  
!===================================================
end module

