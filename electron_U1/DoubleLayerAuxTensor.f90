module DoubleLeyerAuxTensor_type
	use FermiTensor
	use SymDimension_typede
	use Dimension_typede
	use Tools
	implicit none


	interface contract_physical_legs
		module procedure  contract_physical_legs1
		module procedure  contract_physical_legs2
		module procedure  contract_physical_leg1
		module procedure  contract_physical_leg2
	end interface



contains


	subroutine reverse_all_legs(A)! reverse the fermi-arrows of all the non-physical legs for A^+
		type(fTensor),intent(inout)::A(:,:)
		integer::i,j,L1,L2
		character(len=characterlen)::legname1,legname2
		L1=size(A,1)
		L2=size(A,2)
		do i=1,L1
			do j=1,L2
				if(A(i,j)%ifName('R+'))then
					legname1=A(i,j)%getName('R+')
					legname2=A(i,j+1)%getName('L+')
					call ReverseFermiArrow(A(i,j),legname1,A(i,j+1),legname2,.true.)
				end if
				if(A(i,j)%ifName('D+'))then
					legname1=A(i,j)%getName('D+')
					legname2=A(i+1,j)%getName('U+')
					call ReverseFermiArrow(A(i,j),legname1,A(i+1,j),legname2,.true.)
				end if
			end do
		end do
	end subroutine

	subroutine reverse_legs(A)! reverse the fermi-arrows of all the non-physical legs for A^+
		type(fTensor),intent(inout)::A
		integer::i,j,L1,L2
		character(len=characterlen)::legname1
		if(A%ifName('L+'))then
			legname1=A%getName('L+')
			call A%setFermiArrow(legname1,-1*A%getFermiArrow(legname1))
		end if
		if(A%ifName('R+'))then
			legname1=A%getName('R+')
			call ReverseFermiArrow(A,legname1)
		end if
		if(A%ifName('D+'))then
			legname1=A%getName('D+')
			call ReverseFermiArrow(A,legname1)
		end if
		if(A%ifName('U+'))then
			legname1=A%getName('U+')
			call A%setFermiArrow(legname1,-1*A%getFermiArrow(legname1))
		end if
	end subroutine

	subroutine contract_physical_legs1(outA,AL,AR)
		type(fTensor),intent(inout)::outA(:,:)
		type(fTensor),intent(in)::AL(:,:),AR(:,:)
		integer::i,j,L1,L2,contractlen
		character(len=characterlen),allocatable::legname1(:),legname2(:)
		allocate(legname1(2))
		allocate(legname2(2))
		L1=size(AL,1)
		L2=size(AL,2)
		if(size(outA,1).ne.L1)then
			call writemess(' ERROR in leg of the output fTensor')
			call error_stop
		end if
		if(size(outA,2).ne.L2)then
			call writemess(' ERROR in leg of the output fTensor')
			call error_stop
		end if
		do i=1,L1
			do j=1,L2
				contractlen=0
				contractlen=contractlen+1
				legname1(contractlen)=AR(i,j)%getName('n')
				legname2(contractlen)=legname1(contractlen)+'+'
				if(AR(i,j)%ifName('TotalN'))then
					contractlen=contractlen+1
					legname1(contractlen)=AR(i,j)%getName('TotalN')
					legname2(contractlen)=legname1(contractlen)+'+'
				end if
				 outA(i,j)=contract(.Hn.AL(i,j),legname2(1:contractlen),AR(i,j),legname1(1:contractlen))
			end do
		end do
		call reverse_all_legs(outA)
		call fuse_legs(outA)
		return
	end subroutine

	subroutine contract_physical_legs2(outA,A)
		type(fTensor),intent(inout)::outA(:,:)
		type(fTensor),intent(in)::A(:,:)
		integer::i,j,L1,L2,contractlen
		character(len=characterlen),allocatable::legname1(:),legname2(:)
		L1=size(A,1)
		L2=size(A,2)
		if(size(outA,1).ne.L1)then
			call writemess(' ERROR in leg of the output fTensor')
			call error_stop
		end if
		if(size(outA,2).ne.L2)then
			call writemess(' ERROR in leg of the output fTensor')
			call error_stop
		end if
		allocate(legname1(2))
		allocate(legname2(2))
		do i=1,L1
			do j=1,L2
				contractlen=0
				contractlen=contractlen+1
				legname1(contractlen)=A(i,j)%getName('n')
				legname2(contractlen)=legname1(contractlen)+'+'
				if(A(i,j)%ifName('TotalN'))then
					contractlen=contractlen+1
					legname1(contractlen)=A(i,j)%getName('TotalN')
					legname2(contractlen)=legname1(contractlen)+'+'
				end if
				 outA(i,j)=contract(A(i,j),legname1(1:contractlen),.Hn.A(i,j),legname2(1:contractlen))
			end do
		end do
		call reverse_all_legs(outA)
		call fuse_legs(outA)
		return
	end subroutine

	subroutine contract_physical_leg1(outA,AL,AR)
		type(fTensor),intent(inout)::outA
		type(fTensor),intent(in)::AL,AR
		integer::contractlen
		character(len=characterlen),allocatable::legname1(:),legname2(:)
		allocate(legname1(2))
		allocate(legname2(2))
		contractlen=0
		contractlen=contractlen+1
		legname1(contractlen)=AR%getName('n')
		legname2(contractlen)=legname1(contractlen)+'+'
		if(AR%ifName('TotalN'))then
			contractlen=contractlen+1
			legname1(contractlen)=AR%getName('TotalN')
			legname2(contractlen)=legname1(contractlen)+'+'
		end if
		 outA=contract(.Hn.AL,legname2(1:contractlen),AR,legname1(1:contractlen))
		call reverse_legs(outA)
		call fuse_leg(outA)
		return
	end subroutine
	subroutine contract_physical_leg2(outA,A)
		type(fTensor),intent(inout)::outA
		type(fTensor),intent(in)::A
		integer::contractlen
		character(len=characterlen),allocatable::legname1(:),legname2(:)
		allocate(legname1(2))
		allocate(legname2(2))
		contractlen=0
		contractlen=contractlen+1
		legname1(contractlen)=A%getName('n')
		legname2(contractlen)=legname1(contractlen)+'+'
		if(A%ifName('TotalN'))then
			contractlen=contractlen+1
			legname1(contractlen)=A%getName('TotalN')
			legname2(contractlen)=legname1(contractlen)+'+'
		end if
		 outA=contract(A,legname1(1:contractlen),.Hn.A,legname2(1:contractlen))
		call fuse_leg(outA)
		return
	end subroutine

	subroutine fuse_legs(A)
		type(fTensor),intent(inout)::A(:,:)
		integer::i,j,k,L1,L2,rule
		character(len=characterlen)::legname1,legname2
		L1=size(A,1)
		L2=size(A,2)
		do i=1,L1
			do j=1,L2
				if(A(i,j)%ifName('L'))then
					legname1=A(i,j)%getName('L')
					legname2=legname1+'+'
					call A(i,j)%forward([legname1,legname2])
					rule=A(i,j)%getRule(1)
					A(i,j)=A(i,j)%FermiFuse(.true.,rule)
					call A(i,j)%setName(1,legname1)
				end if
				if(A(i,j)%ifName('D'))then
					legname1=A(i,j)%getName('D')
					legname2=legname1+'+'
					call A(i,j)%forward([legname1,legname2])
					rule=A(i,j)%getRule(1)
					A(i,j)=A(i,j)%FermiFuse(.true.,rule)
					call A(i,j)%setName(1,legname1)
				end if
				if(A(i,j)%ifName('R'))then
					legname1=A(i,j)%getName('R')
					legname2=legname1+'+'
					call A(i,j)%forward([legname1,legname2])
					rule=A(i,j)%getRule(1)
					A(i,j)=A(i,j)%FermiFuse(.true.,rule)
					call A(i,j)%setName(1,legname1)
				end if
				if(A(i,j)%ifName('U'))then
					legname1=A(i,j)%getName('U')
					legname2=legname1+'+'
					call A(i,j)%forward([legname1,legname2])
					rule=A(i,j)%getRule(1)
					A(i,j)=A(i,j)%FermiFuse(.true.,rule)
					call A(i,j)%setName(1,legname1)
				end if
			end do
		end do
	end subroutine

	subroutine fuse_leg(A)
		type(fTensor),intent(inout)::A
		character(len=characterlen)::legname1,legname2
		integer::rule
		if(A%ifName('L'))then
			legname1=A%getName('L')
			legname2=legname1+'+'
			call A%forward([legname1,legname2])
			rule=A%getRule(1)
			A=A%FermiFuse(.true.,rule)
			call A%setName(1,legname1)
		end if
		if(A%ifName('D'))then
			legname1=A%getName('D')
			legname2=legname1+'+'
			call A%forward([legname1,legname2])
			rule=A%getRule(1)
			A=A%FermiFuse(.true.,rule)
			call A%setName(1,legname1)
		end if
		if(A%ifName('R'))then
			legname1=A%getName('R')
			legname2=legname1+'+'
			call A%forward([legname1,legname2])
			rule=A%getRule(1)
			A=A%FermiFuse(.true.,rule)
			call A%setName(1,legname1)
		end if
		if(A%ifName('U'))then
			legname1=A%getName('U')
			legname2=legname1+'+'
			call A%forward([legname1,legname2])
			rule=A%getRule(1)
			A=A%FermiFuse(.true.,rule)
			call A%setName(1,legname1)
		end if
		return
	end subroutine

	subroutine checkFermiArrow(T,row,FermiArrowcheck)
		type(fTensor),intent(in)::T
		character*1,intent(in)::row
		integer,intent(in),optional::FermiArrowcheck
		integer::i,FermiArrow,rank
		if(row.eq.'r')then
			FermiArrow=T%getFermiArrow(1)
			if(present(FermiArrowcheck))then
				if(FermiArrow.ne.FermiArrowcheck)then
					call writemess('ERROR in Check FermiArrow,input Tensor is not orthogonal',-1)
					call T%diminfo()
					call error_stop()
				end if
			end if
			FermiArrow=-1*FermiArrow
			do i=2,T%getRank()
				if(FermiArrow.ne.T%getFermiArrow(i))then
					call writemess('ERROR in Check FermiArrow,input Tensor is not orthogonal',-1)
					call T%diminfo()
					call error_stop()
				end if
			end do
		else
			rank=T%getRank()
			FermiArrow=T%getFermiArrow(rank)
			if(present(FermiArrowcheck))then
				if(FermiArrow.ne.FermiArrowcheck)then
					call writemess('ERROR in Check FermiArrow,input Tensor is not orthogonal',-1)
					call T%diminfo()
					call error_stop()
				end if
			end if
			FermiArrow=-1*FermiArrow
			do i=1,T%getRank()-1
				if(FermiArrow.ne.T%getFermiArrow(i))then
					call writemess('ERROR in Check FermiArrow,input Tensor is not orthogonal',-1)
					call T%diminfo()
					call error_stop()
				end if
			end do
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
		type(fTensor),intent(inout)::outT(lenA)
		type(fTensor),intent(inout)::A1(lenA),A2(lenA)
		character(len=1),intent(in)::L,D,R,U
		integer::i,FermiArrows(lenA),ith,k,j
		character*50::nameA
		
		do i=1,lenA
			nameA=A2(i)%getName(D)
			ith=A2(i)%NameOrder(nameA)
			FermiArrows(i)=A2(i)%getFermiArrow(ith)
			call A2(i)%setFermiArrow(ith,1)
			call outT(i)%deallocate()
		end do
		
		
		call QRAndSVDTwoLine(Dc,outT,A1,A2,L,D,R,U)


		do i=1,lenA
			nameA=A2(i)%outTensorName(D)
			call outT(i)%setFermiArrow(nameA,FermiArrows(i))
			call A2(i)%setFermiArrow(nameA,FermiArrows(i))
		end do
		return
	end subroutine



	subroutine QRTwoSiteFirst(QTensor,RTensor,A,B,L,D,R,U)
		character(len=1),intent(in)::L,D,R,U
		type(fTensor),intent(inout)::A,B,QTensor,RTensor
		type(fTensor)::T	
		character(len=50)::Aname,Bname
		T=contract(A,A%getName(D),B,B%getName(U))
		Aname=A%getName(R)
		Bname=B%getName(R)
		call T%QR(QTensor,RTensor,[Aname,Bname],.false.)
		call QTensor%setName(QTensor%getRank(),Bname)
		call RTensor%setName(1,'QR.'+L)
		call RTensor%setName(Aname,'QR.'+R+'1')
		call RTensor%setName(Bname,'QR.'+R+'2')
		call checkFermiArrow(QTensor,'c',-1)
		return
	end subroutine

	subroutine QRTwoSite(QTensor,RTensor,A,B,L,D,R,U)
		character(len=1),intent(in)::L,D,R,U
		type(fTensor),intent(inout)::A,B,QTensor,RTensor
		type(fTensor)::T	
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
		call checkFermiArrow(QTensor,'c',-1)
		return
	end subroutine

	subroutine QRTwoSiteLast(QTensor,RTensor,A,B,L,D,R,U)
		character(len=1),intent(in)::L,D,R,U
		type(fTensor),intent(inout)::A,B,QTensor,RTensor
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
		type(fTensor),intent(inout)::QTensor,UTensor
		integer,intent(in)::Dc
		type(fTensor)::T,S
		character(len=50)::name1
		T=UTensor
		name1=T%getName(L)
		call T%SVD(QTensor,S,UTensor,[name1],Dc,.true.)
		call ReverseFermiArrow(S,UTensor)
		call UTensor%setName(1,name1)
		call QTensor%contract('SVD.U',S,'SVD.s1')
		call QTensor%setName('SVD.s2','SVD.'+R)
		call QTensor%setName(name1,'SVD.'+L)
		call checkFermiArrow(UTensor,'r',-1)
		return
	end subroutine

	subroutine SVDTwoSite(Dc,QTensor,SVDTensor,L,D,R,U)
		character(len=1),intent(in)::L,D,R,U
		type(fTensor),intent(inout)::QTensor,SVDTensor
		integer,intent(in)::Dc
		type(fTensor)::T,S
		character(len=50)::name1,name2,name3,name4
		name1=QTensor%getName(R)
		name2='SVD.'+L
		call T%contract(QTensor,name1,SVDTensor,name2)
		call T%setName('SVD.'+R,name1)
		name3=T%getName(L)
		call T%SVD(SVDTensor,S,QTensor,[name3],Dc,.true.)
		call ReverseFermiArrow(S,QTensor)
		call QTensor%setName(1,name3)
		call SVDTensor%contract('SVD.U',S,'SVD.s1')
		call SVDTensor%setName('SVD.s2','SVD.'+R)
		call SVDTensor%setName(name3,'SVD.'+L)
		call checkFermiArrow(QTensor,'r',-1)
		return
	end subroutine

	subroutine SVDTwoSiteLast(QTensor,SVDTensor,L,D,R,U)
		character(len=1),intent(in)::L,D,R,U
		type(fTensor),intent(inout)::QTensor,SVDTensor
		character(len=50)::name1,name2
		name1=QTensor%getName(R)
		name2='SVD.'+L
		call QTensor%contract(name1,SVDTensor,name2)
		call QTensor%setName('SVD.'+R,name1)
		return
	end subroutine




	subroutine QRAndSVDTwoLine(Dc,outQTensor,A,B,L,D,R,U)! make sure the Fermi-Arrow
		character(len=1),intent(in)::L,D,R,U
		type(fTensor),intent(inout)::outQTensor(:),A(:),B(:)
		integer,intent(in)::Dc
		type(fTensor)::RTensor
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