
	Implicit none
	INTEGER i,j,k,l,n_spec,n_type(20),check,count,counter
	INTEGER total_num,done,v
	character(len=4)::Junk,Atom(20)
	REAL::cutoff(20,20),Cell(3,3),numero,New_Coord(3)
	REAL::Coords(20000,3,20),randoms(3000000),alpha,beta,gamma
	REAL::dimless_s(20000,3,20),dimless_d(20000,3,20)
	REAL::dimless_d2(20000,3,20)
	REAL::unitv(20000,3,20),vec3(20000,3,20)
	REAL::dimlessi_ij(20000,3,20),real_ij(20000,3,20)
	REAL::real_s(20000,3,20),mag,tmp,dimless_ij(20000,3,20)
	REAL::real_d(20000,3,20),mbg

	open(unit=1,file='random_cell_input.dat',status='unknown')
	open(unit=2,file='random.xyz',status='unknown')



	read(1,*)Junk, n_spec					!number of species
	read(1,*)Junk,(Atom(i),i=1,n_spec)
	read(1,*)Junk,(n_type(i),i=1,n_spec)			!number of each type of species
	do i=1,n_spec
	 read(1,*)junk,(cutoff(i,j),j=i,n_spec)			!minimum distances between species types
	enddo
	do i=1,3
	 read(1,*)(Cell(i,j),j=1,3)				!the simulation cell
	enddo
!    Work out the vector angles

	

!	read seed						! random seed
	do i=1,3
	 unitv(1,i,1)=1
	enddo
	total_num=0
	do i=1,n_spec
	 total_num=total_num+n_type(i)
	enddo
	count=0
	done=0
	call random_gene(randoms)
	do i=1,n_spec
	 do j=1,n_type(i)
!	GENERATE THE COORDINTAES OF THE NEW ATOM
101	  count=count+1
	  dimless_s(j,1,i)=randoms(count+4)
	  dimless_s(j,2,i)=randoms(count+3)
	  dimless_s(j,3,i)=randoms(count+7)
! Calculate dimensionless separation vector (sij=(si-si)-INT(2(si-sj)-1))	  
	  do k=1,i-1
	   do l=1,n_type(k)
	   call vector_minus(i,j,k,l,dimless_s,dimless_s,dimless_d)
	    do v=1,3
	     vec3(1,v,1)=ANINT(dimless_d(1,v,1))
	    enddo
	   call vector_minus(1,1,1,1,dimless_d,vec3,dimless_ij)
! Convert dimensionless separation into real space separation
	   call dimless2real(dimless_ij,Cell,real_ij,1,1)
	   mag=SQRT(real_ij(1,1,1)**2+real_ij(1,2,1)**2
     .     +real_ij(1,3,1)**2) 
	   if(mag.le.cutoff(k,i))then
	    goto 101
	   endif						!check=1?
! Add the atom
	   enddo						!l=1,n_type(i)
	  enddo							!k=1,n_spec
! Test against own type

	 do l=1,j-1
	  call vector_minus(i,j,i,l,dimless_s,dimless_s,dimless_d)
            do v=1,3
	     vec3(1,v,1)=ANINT(dimless_d(1,v,1))
            enddo
           call vector_minus(1,1,1,1,dimless_d,vec3,dimless_ij)
! Convert dimensionless separation into real space separation
           call dimless2real(dimless_ij,Cell,real_ij,1,1)
	   call dimless2real(dimless_d,Cell,real_d,1,1)
           mag=SQRT(real_ij(1,1,1)**2+real_ij(1,2,1)**2+
     .     real_ij(1,3,1)**2)
           mbg=SQRT(real_d(1,1,1)**2+real_d(1,2,1)**2+
     .     real_d(1,3,1)**2)

           if(mag.le.cutoff(k,i))then
            goto 101
           endif                                                !check=1?
	 enddo 							!l=1,j	

	  call dimless2real(dimless_s,Cell,real_s,i,j)
	  Coords(j,1,i)=real_s(1,1,1)
	  Coords(j,2,i)=real_s(1,2,1)
	  Coords(j,3,i)=real_s(1,3,1)
	  done=done+1
	  print *, done
	 enddo							!j=1,n_type(i)
	enddo							!i=1,n_spec

	write(2,*)total_num
	write(2,*)"Created By Chance"
	do i=1,n_spec
	do j=1,n_type(i)
	 Write(2,*)Atom(i),(Coords(j,l,i),l=1,3)
	enddo
	enddo


	open(3,file='CONFIG',status='unknown')
	write(3,*)"Amorphous Cell Composition, Si:",n_type(2),
     .		  "N:",n_type(1)
	write(3,*)"         0         3 "
	do i=1,3
	 write(3,11)(CELL(i,j),j=1,3)
	enddo
	counter=1
	do i=1,n_spec
	do j=1,n_type(i)
	 Write(3,12)Atom(i),counter
	 Write(3,13)(Coords(j,l,i),l=1,3)
	 counter=counter+1
	enddo
	enddo
12	FORMAT(A4,12x,I5)
11	FORMAT(7x,F15.10,7x,F15.10,7x,F15.10)
13	FORMAT(4x,F15.10,4x,F15.10,4x,F15.10)
	END
!-----------------------------------------------------------------------------------------
	SUBROUTINE dimless2real(d,Cell,da,i,j)

	implicit none
	INTEGER i,j,k
	REAL::d(20000,3,20),da(20000,3,20),Cell(3,3)

        da(1,1,1)=Cell(1,1)*d(j,1,i)+Cell(2,1)*d(j,2,i)
     .  +Cell(3,1)*d(j,3,i)
	da(1,2,1)=Cell(1,2)*d(j,1,i)+Cell(2,2)*d(j,2,i)
     .  +Cell(3,2)*d(j,3,i)
	da(1,3,1)=Cell(1,3)*d(j,1,i)+Cell(2,3)*d(j,2,i)
     .  +Cell(3,3)*d(j,3,i)

	END Subroutine
!-----------------------------------------------------------------------------------------
	SUBROUTINE vector_minus(i,j,k,l,coords1,coords2,diff)

	INTEGER i,j,k,l
	REAL::coords1(20000,3,20),coords2(20000,3,20),diff(20000,3,20)


	diff(1,1,1)=coords1(j,1,i)-coords2(l,1,k)
	diff(1,2,1)=coords1(j,2,i)-coords2(l,2,k)
	diff(1,3,1)=coords1(j,3,i)-coords2(l,3,k)

	END Subroutine	

!-----------------------------------------------------------------------------------------
 	SUBROUTINE random_gene(i)
      	real rand                 ! Declare the type of the rand() function
	integer seed
	real i(3000000)
      	integer*4 timeArray(3)    ! Holds the hour, minute, and second
! In order to get a different sequence each time, we initialize the
! seed of the random number function with the sum of the current
! hour, minute, and second.
      	call itime(timeArray)     ! Get the current time
	seed=timeArray(1)+timeArray(2)+timeArray(3)
	do j=1,3000000
      	i(j) = ran()
	enddo
      	end subroutine
!-----------------------------------------------------------------------------------------
	SUBROUTINE test(New,Old,a,b,cut,check,Cell)
	INTEGER i,j,k,a,b,check
	REAL::New(3),Old(20000,3,20),cut
	REAL::x1,x2,y1,y2,z1,z2,dis,Cell(3,3)
	check=0					! set the return to initial value of zero
	x1=New(1)
	y1=New(2)
	z1=New(3)

	do i=-1,1
	do j=-1,1
	do k=-1,1
	 x2=Old(b,1,a)+i*Cell(1,1)
	 y2=Old(b,2,a)+j*Cell(2,2)
	 z2=Old(b,3,a)+k*Cell(3,3)
	 dis=SQRT((x1-x2)**2+(y2-y1)**2+(z2-z1)**2)	
	 if(dis.LT.cut)then
	  check=1
	 endif
	enddo
	enddo
	enddo	

	return
	END Subroutine
!-----------------------------------------------------------------------------------------
