program get_TBME
  use ME_general
  implicit none
  
  real(8) :: hw,omp_get_wtime,t1,t2
  integer :: bytes,lines,n,emax,m,i,j,k,l,q,g
  integer,allocatable,dimension(:,:) ::  qnums
  integer,allocatable,dimension(:) :: ord
  real(8),allocatable,dimension(:,:) :: T
  real(8),allocatable,dimension(:) :: eig
  character(5) :: hwstr,nstr,emaxstr
  character(35) ::  fname
  type(full_ham) :: ints
  

    call getarg(1,hwstr)
    call getarg(2,emaxstr)
    
    hwstr=adjustl(hwstr)
    emaxstr=adjustl(emaxstr)
    
    read(hwstr,'(f5.2)') hw
    read(emaxstr,'(I5)') emax
    write(hwstr, '(f5.2)') hw  
    
    hwstr=adjustl(hwstr)
    m= emax*(emax+1) 
     
     allocate( qnums(m,3) )
     allocate( ord(m/2) )
     allocate(T(m,m))
     allocate(eig(m))

     call arrange_states(emax,qnums,ord) 
     call calc_h0(T,hw,qnums) 

  ! figure out what the blocks are ahead of time
  ! so less storage is needed
     do i=1,m
        eig(i) = T(i,i) 
     end do
      
     g = m*(m-1)/2
     call reindex(m, g, 0, 0, g, qnums, emax,eig,ints)
     call convertF(T,ints,m,qnums,eig)
     call get_ints(ints,hw,qnums,m)
          
     fname='../TBMEfiles/TBME_'//trim(hwstr)//'_'//trim(emaxstr)//'.dat'
     
     !print*, fname, 'submitted to TBMEfiles' 
     fname=trim(fname) 
    
     open(unit=37,file=fname)
    
     do q= 1, ints%nblock
        do j=1,ints%mat(q)%npp
           do k=1,ints%mat(q)%npp
              
              write(37,*) ints%mat(q)%Vpppp(j,k)
            
           end do 
        end do 
     end do 

     deallocate(qnums)
     deallocate(ord)
     
     close(37)

end program
