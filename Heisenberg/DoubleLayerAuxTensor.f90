module DoubleLeyerAuxTensor_type
	use Tensor_type
	use Tools
	implicit none


	interface contract_physical_legs
		module procedure  contract_physical_legs1
		module procedure  contract_physical_legs2
		module procedure  contract_physical_leg1
		module procedure  contract_physical_leg2
	end interface



contains



	subroutine contract_physical_legs1(outA,AL,AR)
		type(Tensor),intent(inout)::outA(:,:)
		type(Tensor),intent(in)::AL(:,:),AR(:,:)
		integer::i,j,L1,L2
		character(len=characterlen)::legname1,legname2
		L1=size(AL,1)
		L2=size(AL,2)
		if(size(outA,1).ne.L1)then
			call writemess(' ERROR in leg of the output Tensor')
			call error_stop
		end if
		if(size(outA,2).ne.L2)then
			call writemess(' ERROR in leg of the output Tensor')
			call error_stop
		end if
		do i=1,L1
			do j=1,L2
				legname1=AR(i,j)%getName('n')
				legname2=legname1+'+'
				 outA(i,j)=contract(.Hn.AL(i,j),legname2,AR(i,j),legname1)
			end do
		end do
		call fuse_legs(outA)
		return
	end subroutine

	subroutine contract_physical_legs2(outA,A)
		type(Tensor),intent(inout)::outA(:,:)
		type(Tensor),intent(in)::A(:,:)
		integer::i,j,L1,L2
		character(len=characterlen)::legname1,legname2
		L1=size(A,1)
		L2=size(A,2)
		if(size(outA,1).ne.L1)then
			call writemess(' ERROR in leg of the output Tensor')
			call error_stop
		end if
		if(size(outA,2).ne.L2)then
			call writemess(' ERROR in leg of the output Tensor')
			call error_stop
		end if
		do i=1,L1
			do j=1,L2
				legname1=A(i,j)%getName('n')
				legname2=legname1+'+'
				 outA(i,j)=contract(A(i,j),legname1,.Hn.A(i,j),legname2)
			end do
		end do
		call fuse_legs(outA)
		return
	end subroutine

	subroutine contract_physical_leg1(outA,AL,AR)
		type(Tensor),intent(inout)::outA
		type(Tensor),intent(in)::AL,AR
		character(len=characterlen)::legname1,legname2
		legname1=AR%getName('n')
		legname2=legname1+'+'
		 outA=contract(.Hn.AL,legname2,AR,legname1)
		call fuse_leg(outA)
		return
	end subroutine
	subroutine contract_physical_leg2(outA,A)
		type(Tensor),intent(inout)::outA
		type(Tensor),intent(in)::A
		character(len=characterlen)::legname1,legname2
		legname1=A%getName('n')
		legname2=legname1+'+'
		 outA=contract(A,legname1,.Hn.A,legname2)
		call fuse_leg(outA)
		return
	end subroutine

	subroutine fuse_legs(A)
		type(Tensor),intent(inout)::A(:,:)
		integer::i,j,k,L1,L2,rule
		character(len=characterlen)::legname1,legname2
		L1=size(A,1)
		L2=size(A,2)
		do i=1,L1
			do j=1,L2
				call fuse_leg(A(i,j))
			end do
		end do
	end subroutine

	subroutine fuse_leg(A)
		type(Tensor),intent(inout)::A
		character(len=characterlen)::legname1,legname2
		type(Dimension)::newdimen,oldDimen
		integer::i
		if(A%ifName('L'))then
			legname1=A%getName('L')
			legname2=legname1+'+'
			call A%forward([legname1,legname2])
			oldDimen=A%Dimension()
			call A%fuse(1)
			newdimen=A%dim()
			do i=2,newdimen%getRank()
				call newdimen%setName(i,oldDimen%getName(i+1))
			end do
			call A%resetdim(newdimen)
			call A%setName(1,legname1)
		end if
		if(A%ifName('D'))then
			legname1=A%getName('D')
			legname2=legname1+'+'
			call A%forward([legname1,legname2])
			oldDimen=A%Dimension()
			call A%fuse(1)
			newdimen=A%dim()
			do i=2,newdimen%getRank()
				call newdimen%setName(i,oldDimen%getName(i+1))
			end do
			call A%resetdim(newdimen)
			call A%setName(1,legname1)
		end if
		if(A%ifName('R'))then
			legname1=A%getName('R')
			legname2=legname1+'+'
			call A%forward([legname1,legname2])
			oldDimen=A%Dimension()
			call A%fuse(1)
			newdimen=A%dim()
			do i=2,newdimen%getRank()
				call newdimen%setName(i,oldDimen%getName(i+1))
			end do
			call A%resetdim(newdimen)
			call A%setName(1,legname1)
		end if
		if(A%ifName('U'))then
			legname1=A%getName('U')
			legname2=legname1+'+'
			call A%forward([legname1,legname2])
			oldDimen=A%Dimension()
			call A%fuse(1)
			newdimen=A%dim()
			do i=2,newdimen%getRank()
				call newdimen%setName(i,oldDimen%getName(i+1))
			end do
			call A%resetdim(newdimen)
			call A%setName(1,legname1)
		end if
		return
	end subroutine


!        [1]R--L[2]R--L[3]R--L[4]
!         D      D      D      D
!         |      |      |      |
!         U      U      U      U
!        [1]R--L[2]R--L[3]R--L[4]   
!         D      D      D      D
!         |      |      |      |
!

	subroutine ApproximationDoubleLayer(outT,A1,A2,lenA,Dc,L,D,R,U)
		integer,intent(in)::Dc,lenA
		type(Tensor),intent(inout)::outT(lenA)
		type(Tensor),intent(inout)::A1(lenA),A2(lenA)
		character(len=1),intent(in)::L,D,R,U

		call QRAndSVDTwoLine(Dc,outT,A1,A2,L,D,R,U)
		
		return
	end subroutine



	subroutine QRTwoSiteFirst(QTensor,RTensor,A,B,L,D,R,U)
		character(len=1),intent(in)::L,D,R,U
		type(Tensor),intent(inout)::A,B,QTensor,RTensor
		type(Tensor)::T	
		character(len=50)::Aname,Bname
		T=contract(A,A%getName(D),B,B%getName(U))
		Aname=A%getName(R)
		Bname=B%getName(R)
		call T%QR(QTensor,RTensor,[Aname,Bname],.false.)
		call QTensor%setName(QTensor%getRank(),Bname)
		call RTensor%setName(1,'QR.'+L)
		call RTensor%setName(Aname,'QR.'+R+'1')
		call RTensor%setName(Bname,'QR.'+R+'2')
		return
	end subroutine

	subroutine QRTwoSite(QTensor,RTensor,A,B,L,D,R,U)
		character(len=1),intent(in)::L,D,R,U
		type(Tensor),intent(inout)::A,B,QTensor,RTensor
		type(Tensor)::T	
		character(len=50)::name1,name2,name3,name4
		name1='QR.'+R+'1'
		name2=A%getName(L)
		T=contract(RTensor,name1,A,name2)
		call T%setName(A%getName(R),name1)

		name1='QR.'+R+'2'
		name2=A%getName(D)

		name3=B%getName(L)
		name4=B%getName(U)

		call T%contract([name1,name2],B,[name3,name4])
		call T%setName(B%getName(R),name1)


		call T%QR(QTensor,RTensor,['QR.'+R+'1','QR.'+R+'2'],.false.)
		call QTensor%setName(QTensor%getRank(),B%getName(R))
		call QTensor%setName('QR.'+L,B%getName(L))
		call RTensor%setName(1,'QR.'+L)
		return
	end subroutine

	subroutine QRTwoSiteLast(QTensor,RTensor,A,B,L,D,R,U)
		character(len=1),intent(in)::L,D,R,U
		type(Tensor),intent(inout)::A,B,QTensor,RTensor
		character(len=50)::name1,name2,name3,name4
		name1='QR.'+R+'1'
		name2=A%getName(L)
		QTensor=contract(RTensor,name1,A,name2)


		name1='QR.'+R+'2'
		name2=A%getName(D)

		name3=B%getName(L)
		name4=B%getName(U)

		call QTensor%contract([name1,name2],B,[name3,name4])
		call QTensor%setName('QR.'+L,name3)

		return
	end subroutine




	subroutine SVDTwoSiteFirst(Dc,QTensor,UTensor,L,D,R,U)
		character(len=1),intent(in)::L,D,R,U
		type(Tensor),intent(inout)::QTensor,UTensor
		integer,intent(in)::Dc
		type(Tensor)::T,S
		character(len=50)::name1
		T=UTensor
		name1=T%getName(L)
		call T%SVD(QTensor,S,UTensor,[name1],Dc,.true.)
		S=eye(S)
		call S%setName(1,'SVD.s1')
		call S%setName(2,'SVD.s2')
		call UTensor%setName(1,name1)
		call QTensor%contract('SVD.U',S,'SVD.s1')
		call QTensor%setName('SVD.s2','SVD.'+R)
		call QTensor%setName(name1,'SVD.'+L)
		return
	end subroutine

	subroutine SVDTwoSite(Dc,QTensor,SVDTensor,L,D,R,U)
		character(len=1),intent(in)::L,D,R,U
		type(Tensor),intent(inout)::QTensor,SVDTensor
		integer,intent(in)::Dc
		type(Tensor)::T,S
		character(len=50)::name1,name2,name3,name4
		name1=QTensor%getName(R)
		name2='SVD.'+L
		call T%contract(QTensor,name1,SVDTensor,name2)
		call T%setName('SVD.'+R,name1)
		name3=T%getName(L)
		call T%SVD(SVDTensor,S,QTensor,[name3],Dc,.true.)
		S=eye(S)
		call S%setName(1,'SVD.s1')
		call S%setName(2,'SVD.s2')

		call QTensor%setName(1,name3)
		call SVDTensor%contract('SVD.U',S,'SVD.s1')
		call SVDTensor%setName('SVD.s2','SVD.'+R)
		call SVDTensor%setName(name3,'SVD.'+L)
		return
	end subroutine

	subroutine SVDTwoSiteLast(QTensor,SVDTensor,L,D,R,U)
		character(len=1),intent(in)::L,D,R,U
		type(Tensor),intent(inout)::QTensor,SVDTensor
		character(len=50)::name1,name2
		name1=QTensor%getName(R)
		name2='SVD.'+L
		call QTensor%contract(name1,SVDTensor,name2)
		call QTensor%setName('SVD.'+R,name1)
		return
	end subroutine




	subroutine QRAndSVDTwoLine(Dc,outQTensor,A,B,L,D,R,U)
		character(len=1),intent(in)::L,D,R,U
		type(Tensor),intent(inout)::outQTensor(:),A(:),B(:)
		integer,intent(in)::Dc
		type(Tensor)::RTensor
		integer::i,lenA
		lenA=size(A)
		call QRTwoSiteFirst(outQTensor(1),RTensor,A(1),B(1),L,D,R,U)
		do i=2,lenA-1
			call QRTwoSite(outQTensor(i),RTensor,A(i),B(i),L,D,R,U)
		end do
		call QRTwoSiteLast(outQTensor(lenA),RTensor,A(lenA),B(lenA),L,D,R,U)


		call SVDTwoSiteFirst(Dc,RTensor,outQTensor(lenA),L,D,R,U)
		do i=lenA-1,2,-1
			call SVDTwoSite(Dc,outQTensor(i),RTensor,L,D,R,U)
		end do
		call SVDTwoSiteLast(outQTensor(1),RTensor,L,D,R,U)
		return
	end subroutine

































end module