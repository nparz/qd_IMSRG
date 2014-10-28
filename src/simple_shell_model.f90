program simple_shell_model
  ! generic no core shell model code
  implicit none 
  
  integer :: m,n,R,RP,nx,M2b,i,j,info,emax,ml,ms
  real(8),allocatable,dimension(:,:) :: H,T,V,perc
  real(8),allocatable,dimension(:) :: eig,dwork,SD,SDpair
  integer,allocatable,dimension(:,:) :: qn
  integer,allocatable,dimension(:) :: ord 
  character(5) :: rstr,nstr,emaxstr,mlstr,msstr
  character(10) :: sstr,offstr
  character(1) :: yes,evecs
  real(8) :: gx,s,offset,t1,t2,v_elem,omp_get_wtime
  real(16) :: bin_coef
 
  
  yes = 'y' 
  call getarg(1,emaxstr)
  call getarg(2,nstr)
  call getarg(3,sstr)
  call getarg(4,offstr)
  call getarg(5,evecs)
  call getarg(6,mlstr)
  call getarg(7,msstr) 
  
  read(sstr,'(f10.5)') s
  read(offstr,'(f10.5)') offset
  read(nstr,'(I5)') n
  read(emaxstr,'(I5)') emax
  read(mlstr,'(I5)') ml
  read(msstr,'(I5)') ms
  
  m = emax*(emax+1)

  M2b = m*(m-1)/2
  R = nint(bin_coef(m,n))
  
  allocate(qn(emax*emax+emax,3),ord((emax*emax+emax)/2)) 
  call arrange_states(emax,qn,ord) 
  
  allocate(SD(R))
  
  call odometer(SD,m,n,R) ! calculate all possible SDs
  !if you don't want broken pairs: 
  if ((yes == 'y') .or. (yes == 'Y')) then 
     !! this drastically reduces the number of SDs

     allocate(SDpair(R)) 
     RP = 0
     call restrict_to_qns(SD,SDpair,M,N,R,RP,qn,emax,ml,ms) 
     deallocate(SD) 
     allocate(SD(RP)) 
     SD = SDpair(1:RP) 
     deallocate(SDpair) ! get rid of this
     R = RP
     ! continue as if you never had considered the broken pairs...
  end if 
  
  allocate(eig(R))
  allocate(dwork(10*R))
  allocate(H(R,R))
  allocate(T(m,m)) 
  allocate(V(M2b,M2b))
  allocate(perc(R,N+1))
  
  
  !! SD is the ordered list of Slater Determinants
  !! T is the one-body interaction
  !! V is the two-body interaction
  !! H is the hamiltonian

  call import_interaction(T,V,m) 

  !t1 = omp_get_wtime()
  call construct_H(SD,H,T,V,R,m,M2b,N)
  !t2 = omp_get_wtime()
  !print*, t2-t1
  call dsyev('V','U',R,H,R,eig,dwork,10*R,info) 
  
 ! t1 = omp_get_wtime()
  !print*, t1-t2
  
  if (evecs == 'y') then 
     call find_percentages(SD,H,N,R,M,perc)
     open(unit=37,file = 'level_percentages.dat') 
     do i =1,R
        write(37,*) perc(i,:)
     end do 
     close(37)
  end if
  
  open(unit = 38,file = 'CI_spectrum.dat',position='append') 
      write(38,*) s, offset+eig 
  close(38) 
  
end program
!==================================================
!==================================================
subroutine odometer(SD_bin,M,N,R)
  ! iteratively construct Hilbert space
  implicit none 
  
  integer :: M,N,R,i,q,j
  integer,dimension(R,N) :: SD
  real(8),dimension(R) :: SD_bin
   
  SD_bin = 0
  ! reference state
  do i=1,N
     SD(1,i) = i 
     SD_bin(1) = SD_bin(1) + 2.**(i-1) 
  end do 
  
  do i = 2, R 
     q = N 
     SD(i,:) = SD(i-1,:)
     ! ensure that the rightmost index increases if it can
     do 
        if ( SD(i-1,q) < (M-N+q) ) then 
           SD(i,q) = SD(i-1,q) + 1
           if (q < N) then 
              SD(i,q+1:N) = SD(1,q+1:N)+SD(i,q)-SD(1,q)
           end if 
           exit
        else
           q = q-1
        end if 
     
     end do
     do j=1,N
        SD_bin(i) = SD_bin(i) + 2.**(SD(i,j)-1) 
     end do
       
  end do 
       
end subroutine
!===============================================
!===============================================
subroutine import_interaction(H0,V,m) 
  implicit none 
  
  integer :: i,j,m,MX
  real(8),dimension(m,m) :: H0 !one body
  real(8),dimension(m*(m-1)/2,m*(m-1)/2) :: V !two body
  
  MX = m*(m-1)/2
  ! import one body matrix elements
  open(unit=39,file='OBME.dat',form='unformatted')
  H0 = 0.d0
  do i=1,m
     do j=i,m
        
        read(39) H0(i,j) 
        H0(j,i) = H0(i,j) 
     end do 
  end do 
  close(39)
  
  ! import two body matrix elements
  open(unit=39,file='TBME.dat',form='unformatted')
  V = 0 
  do i=1,MX
     do j=i,MX
   
        read(39) V(i,j)
        V(j,i) = V(i,j)
     end do 
  end do 
  close(39)
  
! note V is stored in two particle basis, ordered by odometer
end subroutine 
!==========================================
!==========================================
subroutine construct_H(SD,H,T,V,RX,m,M2,N)
 ! builds hamiltonian in N-body fock space
 implicit none 
  
 integer :: RX,m,M2,N,i,j,k,x,p,q,r,s,fac,qt
 integer :: b1(m),b2(m),res(m),sm(m)
 real(8) :: SD(RX),H(RX,RX), T(m,m), V(M2,M2),v_elem
 
 H = 0.d0
!$OMP PARALLEL DO PRIVATE( b1,b2,res,sm,x,qt,k,j,p,q,r,s,fac) 
 do i = 1, RX
    do j = i+1,RX
       ! use binary strings to compute H
       call to_binary(b1,SD(i),M)
       call to_binary(b2,SD(j),M)
        res = b1 - b2 
        sm = b1 + b2
     
        x = sum(abs(res))

        
       if (x > 4) cycle 
       
       if (x == 4) then 
          ! 2 body excitation

         
          ! find the indeces of the excitation
          p = 1
          do 
             if (res(p) == 1)  then 
                q = p+1
                exit
             end if 
             p = p+1 
          end do 
          
          r = 1
          do 
             if (res(r) == -1)  then 
                s = r+1
                exit
             end if 
             r = r+1 
          end do 
          
           do 
               if (res(q) == 1) exit
               q = q + 1
            end do 
            
             do 
               if (res(s) == -1) exit
               s = s + 1
            end do 
            
           
            ! calculate phase factor ( how many filled states between )
            sm(p) = 0;sm(q) = 0;sm(r) = 0;sm(s) = 0
            sm = sm / 2
            fac = (-1)**(sum(sm(p:q))+sum(sm(r:s))) 
            
            ! calculate matrix element
        
            H(i,j) = v_elem(p,q,r,s,V,M2,M)*fac
           
            H(j,i) = H(i,j)
          

       else
          ! one body excitation
         
            ! find indeces of the excitation
            p = 1
            
            do 
               if (abs(res(p)) == 1) then 
                  q = p+1 
                  exit
               end if 
               p = p+1
            end do 
            
            do 
               if (abs(res(q)) == 1) exit
               q = q + 1
            end do 
            
            ! calculate phase factor ( how many filled states between ) 
            sm(p) = 0; sm(q) = 0;sm = sm /2

            fac = (-1) **sum(sm(p:q)) 
            
    
            ! calculate matrix element
            H(i,j) = T(p,q)
            qt = 0
            k = 1

            do while (qt < N-1 )
                  H(i,j) = H(i,j) + v_elem(k,p,k,q,V,M2,M)*sm(k)   
                  qt = qt+sm(k)
                  k = k+1
            end do 

            H(i,j) = H(i,j)*fac
            H(j,i) = H(i,j)

       end if 
    end do 
 end do
!$OMP END PARALLEL DO 

! get the diagonal quickly
 do i=1,RX
    call to_binary(b1,SD(i),M)

    qt = 0
    k = 1
    do while (qt < N) 
       H(i,i) = H(i,i) + T(k,k)*b1(k)
      
       qt = qt + b1(k)
       k = k+1
    end do 
     
    qt = 0
    j = 1
    do while ( qt < n*(n-1)/2 ) 
       if (b1(j) == 1) then 
          k = j+1
          do while (sum(b1(k:m)) > 0 )              
             H(i,i) = H(i,i) + v_elem(j,k,j,k,V,M2,M)*b1(k)
             qt = qt + b1(k)
             k = k + 1
          end do 
       end if 
       j = j + 1
    end do

end do        
  
end subroutine
!==========================================
!==========================================
subroutine to_binary(A,binout,M)
  ! writes an integer in binary
  implicit none 
  
  integer :: M,i
  integer,dimension(M) :: A
  real(8) :: binout,x
  
  x = binout
  A = 0
  i = 1
  do while (x > 0)
    
     if (abs(modulo(x,2.)) <  1e-6) then 
        x = x/2 
     else 
        x = (x-1)/2 
        A(i) = 1
     end if 
     i = i+1 
  end do 

end subroutine 
!=========================================
!=========================================
real(8) function v_elem(px,qx,rx,sx,V,x,M) 
  ! grabs the TBME for pqrs
  implicit none 
  
  integer :: p,q,r,s,x,i,j,M,a,b,pre
  integer :: px,qx,rx,sx
  real(8),dimension(x,x) :: V
  real(8) :: out
  
  p = px
  q = qx
  r = rx 
  s = sx
  pre = 1
 
  ! formula only works for p < q 
  ! if this isn't true, switch them, 
  ! and note this with the prefactor
  if (p > q) then 
     a = p 
     p = q 
     q = a
     pre = -1 
  end if 
  
  if (r > s) then 
     a = r 
     r = s 
     s = a
     pre = -1 * pre 
  end if 
   
  if ((q == p ) .or. (s == r )) then 
     v_elem = 0.d0
  else 
     
  ! map from sp indeces to tp indeces
  i = (p-1) * M  - p*(p-1)/2 + q - p
  j = (r-1) * M  - r*(r-1)/2 + s - r
  
  
  v_elem = V(i,j) * pre
  end if
  
end function
!==========================================
!==========================================
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
subroutine restrict_to_qns(SD,SD_PAIR,M,N,R,R_P,qn,emax,ml,ms) 
  !!! restricts basis to filled spatial orbitals
  implicit none
  
  integer :: M,N,R,R_P,i,j,k,emax,ml,ms,mlsum,mssum,q
  real(8) :: SD(R),SD_PAIR(R)
  integer :: bin(M) 
  integer,dimension(emax*emax+emax,3) :: qn
  
  k = 1
  
  
  do i=1, R
     call to_binary(bin,SD(i),M)
     j = 0
     q = 1
     mlsum = 0
     mssum = 0 
     do while (j < N)
        if (bin(q) == 1) then
          
           mlsum = mlsum + qn(q,2) 
           mssum = mssum + qn(q,3) 
           j = j + 1
        end if
        q = q + 1
       
     end do 
     if ((mlsum == ml) .and. (mssum == ms) ) then  
        SD_pair(k) = SD(i)
        k = k + 1
   
     end if 
  end do 
  R_P = k-1
end subroutine 
!===========================================
subroutine find_percentages(SD,H,N,R,M,perc) 
  ! calculates the percentages of the states which 
  ! correspond to different levels of excitation
  implicit none
  
  integer :: R,N,M,i,j,level(R),kx
  real(8) :: SD(R),perc(R,N+1),H(R,R)
  integer,dimension(M) :: bin
  
  ! figure out which SDs are which npnh excitation levels
  do i = 1, R
     call to_binary(bin,SD(i),M)
     kx = sum(bin(1:N)) 
     level(i) = N-kx 
  end do    
  
  ! calculate percentages
  perc = 0.d0
  do i = 1, R
     do j = 1, R
        perc(i,level(j)+1) = perc(i,level(j)+1) + &
             H(j,i)**2
     end do 
  end do 
end subroutine    
!=============================================
subroutine arrange_states(maxE,q_actual,order)
  implicit none 
  !!! states are ordered by these criteria:
 
  !!! 1. lowest spin will have lowest ordering
  !!! 2. lowest abs(m) will have lowest ordering 
  !!! 3. negative m is lower than positive
  !!! 4. for same m, lowest n will have lowest ordering
  
  !!! this way we have a block diagonal hamiltonian
  !!! so LAPACK won't have a stroke... 

  integer,dimension((maxE*maxE+maxE),3) :: quant_num,q_actual 
  integer,dimension((maxE*maxE+maxE)/2) :: order,en
  integer :: i,maxE,a,b,j,m,n,R,Nprin
  
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

j = 1 

do R = 1, maxE
   do Nprin = 0, (R-1)/2
   do i = 1, maxE*(maxE+1) 
     If (quant_num(i,1) == Nprin) then 
        if (abs(quant_num(i,2)) == R -2*Nprin - 1 ) then 
           q_actual(j,:) = quant_num(i,:) 
           j = j+1
        end if
     end if
     end do 
   end do 
end do 
     
end subroutine 
