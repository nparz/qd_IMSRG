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
subroutine build_specific_space( gen, rec ) 
  ! decouples 1p1h from 2p2h for specific, specifed set of states 
  ! also considers only 1p1h exciations into a valence space (rec%cutshell) 

! labeling scheme ---------------------------------
! p - particle
! h - hole
! q - non-valence particle
! v - valence particle
! a - particle in specific set |(a)(i-)>
! i - hole in specific set |(a)(i-)>
!---------------------------------------------------  
  implicit none
 
  
  integer :: i,j,k,q,q1,q2,a,b,n,m,p(2),h(2),be(2),he(2),pe(2),pp,row(2,3)
  integer :: sz,pos
  type(full_ham) ::  gen , rec 
  real(8) :: efs, Aph,mx,rb,denom
  logical :: in,inSD
  
  gen%herm = -1
  n=rec%nbody
  m=rec%Msp 
   
  sz = size( rec%exlabels(:,2)) 
  

! decouple one body
  ! ph 
  do a=n+1,m
     do i=1,n
       
        denom = rec%fpp(a-n,a-n) - rec%fhh(i,i) - v_elem(a,i,a,i,rec)
        gen%fph(a-n,i) = rec%fph(a-n,i) * &
      sign(1.d0 ,denom ) * abs(denom)**.0001 
        
     end do
  end do 
  
  !qa
  do a = n+1,rec%cutshell
     do b = rec%cutshell,m
        if ( .not. in(a,rec%exlabels(:,2),pos,sz)) cycle
        denom = rec%fpp(a-n,a-n) - rec%fpp(b-n,b-n) -v_elem(a,b,a,b,rec) 
        gen%fpp(a-n,b-n) = rec%fpp(a-n,b-n) * sign(1.d0,denom) * abs(denom)**.0001
        gen%fpp(b-n,a-n) = -1*gen%fpp(a-n,b-n) 
     end do 
  end do 

  !!! two body term !!!

  mx = 0.
  do q=1,rec%nblock
     
     if ( rec%mat(q)%nph*rec%mat(q)%nhh > 0 ) then 

        ! pihh
     do i=1,gen%mat(q)%nph
        do j=1,gen%mat(q)%nhh
           
           h(:)=rec%mat(q)%qnhh(j,:)
          
           be(1)=max(rec%mat(q)%qnph(I,1),rec%mat(q)%qnph(I,2))
           be(2)=min(rec%mat(q)%qnph(I,1),rec%mat(q)%qnph(I,2))
           he(1)=h(1)
           he(2)=h(2)
           
           if (.not. in(be(2),rec%exlabels(:,1),pos,sz) ) cycle 
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
         !ppah
     do i=1,gen%mat(q)%npp
        do j=1,gen%mat(q)%nph
           
           p(:)=rec%mat(q)%qnpp(i,:)
       
           be(1)=max(rec%mat(q)%qnph(J,1),rec%mat(q)%qnph(J,2))
           be(2)=min(rec%mat(q)%qnph(J,1),rec%mat(q)%qnph(J,2))
           if (be(1) > rec%cutshell) cycle
           if (.not. in(be(1),rec%exlabels(:,2),pos,sz)) cycle
           pe(1)=p(1)
           pe(2)=p(2)
          
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
         ! qiah
     do i=1,gen%mat(q)%nph
        pe(1)=max(rec%mat(q)%qnph(I,1),rec%mat(q)%qnph(I,2))
        pe(2)=min(rec%mat(q)%qnph(I,1),rec%mat(q)%qnph(I,2))
        if (pe(1) > rec%cutshell) cycle
        
        do j=i+1,gen%mat(q)%nph
              
           be(1)=max(rec%mat(q)%qnph(J,1),rec%mat(q)%qnph(J,2))
           be(2)=min(rec%mat(q)%qnph(J,1),rec%mat(q)%qnph(J,2))
           
           if (be(1) < rec%cutshell+1) cycle
           if (.not. inSD(be(2),pe(1),rec%exlabels,pos,sz)) cycle
           
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
subroutine build_imaginary_time( gen, rec ) 
    ! decouples 1p1h from 2p2h for specific, for only
    ! 1p1h excitations into a specified valence space (rec%cutshell) 
! labeling scheme ---------------------------------
! p - particle
! h - hole
! q - non-valence particle
! v - valence particle
! a - particle in specific set |(a)(i-)>
! i - hole in specific set |(a)(i-)>
!---------------------------------------------------  
  implicit none
 
  
  integer :: i,j,k,q,q1,q2,a,b,n,m,p(2),h(2),be(2),he(2),pe(2),pp,row(2,3)
  integer :: sz,pos
  type(full_ham) ::  gen , rec 
  real(8) :: efs, Aph,mx,rb,denom
  logical :: in,inSD
  
  gen%herm = -1
  n=rec%nbody
  m=rec%Msp 
  
  ! one body term
  !ph
  do a=n+1,m
     do i=1,n
       
        denom = rec%fpp(a-n,a-n) - rec%fhh(i,i) - v_elem(a,i,a,i,rec)
        gen%fph(a-n,i) = rec%fph(a-n,i) * &
      sign(1.d0 ,denom ) * abs(denom)**.0001 
      
        
     end do
  end do 

! qv
  do a = n+1,rec%cutshell
     do b = rec%cutshell,m
        denom = rec%fpp(a-n,a-n) - rec%fpp(b-n,b-n) -v_elem(a,b,a,b,rec) 
        gen%fpp(a-n,b-n) = rec%fpp(a-n,b-n) * sign(1.d0,denom) * abs(denom)**.0001
        gen%fpp(b-n,a-n) = -1*gen%fpp(a-n,b-n) 
     end do 
  end do 

  !!! two body term !!!

  mx = 0.
  do q=1,rec%nblock
     
     if ( rec%mat(q)%nph*rec%mat(q)%nhh > 0 ) then 
        
     ! phhh    
     do i=1,gen%mat(q)%nph
        do j=1,gen%mat(q)%nhh
           
           h(:)=rec%mat(q)%qnhh(j,:)
          
           be(1)=max(rec%mat(q)%qnph(I,1),rec%mat(q)%qnph(I,2))
           be(2)=min(rec%mat(q)%qnph(I,1),rec%mat(q)%qnph(I,2))
           he(1)=h(1)
           he(2)=h(2)
           
 
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
         
     !ppvh    
     do i=1,gen%mat(q)%npp
        do j=1,gen%mat(q)%nph
           
           p(:)=rec%mat(q)%qnpp(i,:)
       
           be(1)=max(rec%mat(q)%qnph(J,1),rec%mat(q)%qnph(J,2))
           be(2)=min(rec%mat(q)%qnph(J,1),rec%mat(q)%qnph(J,2))
           if (be(1) > rec%cutshell) cycle
 
           pe(1)=p(1)
           pe(2)=p(2)
          
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
     ! qhvh     
     do i=1,gen%mat(q)%nph
        pe(1)=max(rec%mat(q)%qnph(I,1),rec%mat(q)%qnph(I,2))
        pe(2)=min(rec%mat(q)%qnph(I,1),rec%mat(q)%qnph(I,2))
        if (pe(1) > rec%cutshell) cycle
        
        do j=i+1,gen%mat(q)%nph
              
           be(1)=max(rec%mat(q)%qnph(J,1),rec%mat(q)%qnph(J,2))
           be(2)=min(rec%mat(q)%qnph(J,1),rec%mat(q)%qnph(J,2))
           
           if (be(1) < rec%cutshell+1) cycle
 
           
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
           
           pe(1)=p(1)
           pe(2)=p(2)
           he(1)=h(1)
           he(2)=h(2)
           
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
  allocate(r2%states(m,3))
  r2%states = rec%states
  r2%stoe=rec%stoe
  r2%eh=rec%eh
  r2%ep=rec%ep
  r2%mltarg = rec%mltarg
  r2%mstarg = rec%mstarg
  
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
!========================================
!========================================
subroutine allocate_ex(rec,r2) 
  implicit none 
  
  type(full_ham) :: rec,r2
  integer :: i,n,m,m3,g,pp,hh
  
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
  allocate(r2%states(m,3))
  r2%states = rec%states
  r2%stoe=rec%stoe
  r2%eh=rec%eh
  r2%ep=rec%ep
  r2%mltarg = rec%mltarg
  r2%mstarg = rec%mstarg
  
  allocate(r2%i_array(g,3))
  
  r2%i_array=rec%i_array

  allocate(r2%fph(m-n,n))
  
  r2%nblock=rec%nblock
  allocate( r2%mat( 1 ) )
  r2%fph=0.d0

  pp = (m-n)*(m-n-1)/2 
  hh = n*(n-1)/2 
  
  allocate(r2%mat(1)%Vpphh(pp,hh)) 
  r2%mat(1)%Vpphh = 0.d0 
  
  r2%IMSRG3 = .false.
  
  if (rec%IMSRG3) then 
     r2%IMSRG3 = .true.
     allocate(r2%V3body(m3,m3))
     allocate(r2%threemap(m3,3))
     r2%threemap = rec%threemap
  end if 
  
end subroutine  

end module




!=====================================================
!=====================================================
logical function in(element,list,position,sz) 
  implicit none 
  
  integer :: sz
  integer,dimension(sz) :: list
  integer :: element,position,i
  logical :: dum
 
  dum = .false.
  position = -1 
  do i = 1,size(list)
     if ( list(i) == element ) then 
        dum = .true. 
        exit
        position = i 
     end if 
  end do 

  in = dum 
end function
!=====================================================
!=====================================================
logical function inSD(i1,i2,list,position,sz) 
  implicit none 
  
  integer :: sz
  integer,dimension(sz,2) :: list
  integer :: position,i,i1,i2
  logical :: dum
  
  dum = .false. 
  position = -1 
  do i = 1,size(list(:,1)) 
     if ( list(i,1) == i1) then 
        if ( list(i,2) == i2 ) then 
           position = i 
           dum = .true. 
           exit
        end if 
     end if 
  end do 
  
  inSD = dum 
end function
