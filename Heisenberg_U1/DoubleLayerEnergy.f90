module DoubleLayerEnergy_type
	use SimpleUpdate
	use DoubleLeyerAuxTensor_type
	use SymTensor_type
	use Tensor_type
	use SymDimension_typede
	use Dimension_typede
	use Tools
	implicit none 
	type(SymTensor),allocatable::AuxUpOrLeft(:,:)
	type(SymTensor),allocatable::AuxDownOrRight(:,:)
	integer::AuxFlag=0
	!AuxFlag=0: no AuxTensor
	!AuxFlag=1: up and Down
	!AuxFlag=2: left and Right


contains


	
	function DirectContractNetwork_V(A)result(outRes)
		real*8::outRes
		type(SymTensor),intent(inout)::A(:,:)
		integer::i,j,lenA1,lenA2
		type(SymTensor)::Res
		character(len=characterLen)::name1,name2,name3,name4
		lenA1=size(A,1)
		lenA2=size(A,2)
		IF(lenA1.gt.4)then
			call writemess('ERROR in DirectContractNetwork. The input network is too large')
			call error_stop
		end if
		do j=1,lenA2
			if(j.eq.1)then
				name1=A(1,j)%getName('D')
				name2=A(2,j)%getName('U')
				call Res%contract(A(1,j),name1,A(2,j),name2)
				do i=3,lenA1
					name1=A(i-1,j)%getName('D')
					name2=A(i,j)%getName('U')
					call Res%contract(name1,A(i,j),name2)
				end do
			else
				name1=A(1,j-1)%getName('R')
				name2=A(1,j)%getName('L')
				call Res%contract(name1,A(1,j),name2)
				do i=2,lenA1
					name1=A(i-1,j)%getName('D')
					name2=A(i,j-1)%getName('R')
					name3=A(i,j)%getName('U')
					name4=A(i,j)%getName('L')
					call Res%contract([name1,name2],A(i,j),[name3,name4])
				end do
			end if
		end do
		outRes=Res
		return
	end function

	function DirectContractNetwork_H(A)result(outRes)
		real*8::outRes
		type(SymTensor),intent(inout)::A(:,:)
		integer::i,j,lenA1,lenA2
		type(SymTensor)::Res
		character(len=characterLen)::name1,name2,name3,name4
		lenA1=size(A,1)
		lenA2=size(A,2)
		IF(lenA2.gt.4)then
			call writemess('ERROR in DirectContractNetwork. The input network is too large')
			call error_stop
		end if

		do i=1,lenA1
			if(i.eq.1)then
				name1=A(i,1)%getName('R')
				name2=A(i,2)%getName('L')
				call Res%contract(A(i,1),name1,A(i,2),name2)
				do j=3,lenA2
					name1=A(i,j-1)%getName('R')
					name2=A(i,j)%getName('L')
					call Res%contract(name1,A(i,j),name2)
				end do
			else
				name1=A(i-1,1)%getName('D')
				name2=A(i,1)%getName('U')
				call Res%contract(name1,A(i,1),name2)
				do j=2,lenA2
					name1=A(i,j-1)%getName('R')
					name2=A(i-1,j)%getName('D')
					name3=A(i,j)%getName('L')
					name4=A(i,j)%getName('U')
					call Res%contract([name1,name2],A(i,j),[name3,name4])
				end do
			end if
		end do
		outRes=Res
		return
	end function




	subroutine initialUpDown(A,Dc)
		integer,intent(in)::Dc
		type(SymTensor),intent(inout)::A(:,:)
		integer::i,L1,L2
		L1=size(A,1)
		L2=size(A,2)
		AuxUpOrLeft(1,:)=A(1,:)
		call ApproximationDoubleLayer(AuxUpOrLeft(2,:),A(1,:),A(2,:),L2,Dc,'L','D','R','U')
		do i=3,L1-1
			call ApproximationDoubleLayer(AuxUpOrLeft(i,:),AuxUpOrLeft(i-1,:),A(i,:),L2,Dc,'L','D','R','U')
		end do

		AuxDownOrRight(L1,:)=A(L1,:)
		i=L1-1
		call ApproximationDoubleLayer(AuxDownOrRight(i,:),A(L1,:),A(i,:),L2,Dc,'L','U','R','D')
		do i=L1-2,2,-1
			call ApproximationDoubleLayer(AuxDownOrRight(i,:),AuxDownOrRight(i+1,:),A(i,:),L2,Dc,'L','U','R','D')
		end do
		AuxFlag=1
		return
	end subroutine



	subroutine initialLeftRight(A,Dc)
		integer,intent(in)::Dc
		type(SymTensor),intent(inout)::A(:,:)
		integer::j,L1,L2
		L1=size(A,1)
		L2=size(A,2)

		AuxUpOrLeft(:,1)=A(:,1)
		call ApproximationDoubleLayer(AuxUpOrLeft(:,2),A(:,1),A(:,2),L1,Dc,'U','R','D','L')
		do j=3,L2-1
			call ApproximationDoubleLayer(AuxUpOrLeft(:,j),AuxUpOrLeft(:,j-1),A(:,j),L1,Dc,'U','R','D','L')
		end do

		AuxDownOrRight(:,L2)=A(:,L2)
		j=L2-1
		call ApproximationDoubleLayer(AuxDownOrRight(:,j),A(:,L2),A(:,j),L1,Dc,'U','L','D','R')
		do j=L2-2,2,-1
			call ApproximationDoubleLayer(AuxDownOrRight(:,j),AuxDownOrRight(:,j+1),A(:,j),L1,Dc,'U','L','D','R')
		end do
		AuxFlag=2
		return
	end subroutine



	subroutine Phi_O_Phi_H(outA1,outA2,O,A1,A2)
		type(SymTensor),intent(in)::O,A1,A2
		type(SymTensor),intent(inout)::outA1,outA2
		type(SymTensor)::C,U,S,V
		character(len=100)::name1,name2
		name1=A1%getName('n')
		name2=A2%getName('n')
		C=contract(A1,A1%getName('R'),A2,A2%getName('L'))
		call C%contract(O,['H1.n2','H2.n2'],[name1,name2])
		call C%setName('H1.n',name1)
		call C%setName('H2.n',name2)
		name1=name1.subl.'.'
		name2=name2.subl.'.'
		call C%SVD(U,S,V,name1,name2)
		call U%contract('SVD.U',S,'SVD.s1')
		call U%setName('SVD.s2',A1%getName('R'))
		call V%setName('SVD.V',A2%getName('L'))

		CALL contract_physical_legs(outA1,A1,U)
		CALL contract_physical_legs(outA2,A2,V)
		return
	end subroutine

	subroutine Phi_O_Phi_V(outA1,outA2,O,A1,A2)
		type(SymTensor),intent(in)::O,A1,A2
		type(SymTensor),intent(inout)::outA1,outA2
		type(SymTensor)::C,U,S,V
		character(len=100)::name1,name2
		name1=A1%getName('n')
		name2=A2%getName('n')
		C=contract(A1,A1%getName('D'),A2,A2%getName('U'))
		call C%contract(O,['H1.n2','H2.n2'],[name1,name2])
		call C%setName('H1.n',name1)
		call C%setName('H2.n',name2)
		name1=name1.subl.'.'
		name2=name2.subl.'.'
		call C%SVD(U,S,V,name1,name2)
		call U%contract('SVD.U',S,'SVD.s1')
		call U%setName('SVD.s2',A1%getName('D'))
		call V%setName('SVD.V',A2%getName('U'))


		CALL contract_physical_legs(outA1,A1,U)
		CALL contract_physical_legs(outA2,A2,V)
		return
	end subroutine


	subroutine Phi_H_Phi_two_site(outA,outB,A,B,R,L,phya,phyb,H)
		type(SymTensor),intent(inout)::outA,outB
		type(SymTensor),intent(in)::A,B
		type(SymTensor),intent(in)::H
		character(len=*),intent(in)::R,L,phya,phyb
		character(len=50)::Aname,Bname,Anamephy,Bnamephy,name1,name2
		type(SymTensor)::C,U,S,V,QR(2),LQ(2)

		name1=A%getName(R)
		Anamephy=A%getName(phya)
		Aname=name1.subl.'.'
		call A%QR(QR(1),QR(2),[Anamephy,name1],.false.)
		if(QR(2)%getRank().ne.3)then
			call writemess("ERROR in QR in simple update")
			call QR(2)%diminfo()
			call error_stop
		end if
		call QR(1)%setName('QR.Q',name1)
		call QR(2)%setName('QR.R',Aname+'.'+L)
		
		
		name2=B%getName(L)
		Bnamephy=B%getName(phyb)
		Bname=name2.subl.'.'
		call B%LQ(LQ(1),LQ(2),[name2,Bnamephy])
		if(LQ(1)%getRank().ne.3)then
			call writemess("ERROR in LQ in simple update")
			call LQ(1)%diminfo()
			call error_stop
		end if
		call LQ(1)%setName('LQ.L',Bname+'.'+R)
		call LQ(2)%setName('LQ.Q',name2)
		
		
		C=contract(QR(2),name1,LQ(1),name2)



		C=contract(H,['H1.n2','H2.n2'],C,[Anamephy,Bnamephy])
		call C%setName('H1.n',Anamephy)
		call C%setName('H2.n',Bnamephy)
		
		call C%SVD(U,S,V,Aname,Bname)

		call U%contract('SVD.U',S,'SVD.s1')
		call U%setName('SVD.s2',name1)
		call U%contract(QR(1),name1,Aname+'.'+L)


		call V%setName('SVD.V',name2)
		call V%contract(Bname+'.'+R,LQ(2),name2)


		CALL contract_physical_legs(outA,A,U)
		CALL contract_physical_legs(outB,B,V)

		
		return
	end subroutine	


	function ValueH(ith,jth,H,singleA,DoubleA,norm,phya,phyb)
		real*8::ValueH
		integer,intent(in)::ith,jth
		character(len=*),intent(in)::phya,phyb
		type(SymTensor),intent(in)::H,singleA(:,:),DoubleA(:,:)
		real*8,intent(in)::norm
		type(SymTensor)::A1,A2
		integer::j,L1,L2,row,rowi
		type(SymTensor),allocatable::Network(:,:)
		if(AuxFlag.ne.1)then
			call writemess('ERROR AuxFlag='+AuxFlag)
			call writemess(' please reset AuxTensor')
			call error_stop
		end if
		L1=size(singleA,1)
		L2=size(singleA,2)
		row=1
		if(ith.gt.1)row=row+1
		if(ith.lt.L1)row=row+1
		allocate(Network(row,L2))

		call Phi_H_Phi_two_site(A1,A2,singleA(ith,jth),singleA(ith,jth+1),'R','L',phya,phyb,H)
		rowi=0
		if(ith.gt.1)then
			rowi=rowi+1
			do j=1,L2
				Network(rowi,j)=AuxUpOrLeft(ith-1,j)
			end do
		end if
		rowi=rowi+1
		do j=1,L2
			if(j.eq.jth)then
				Network(rowi,j)=A1
			else if(j.eq.(jth+1))then
				Network(rowi,j)=A2
			else
				Network(rowi,j)=DoubleA(ith,j)
			end if
		end do
		if(ith.lt.L1)then
			rowi=rowi+1
			do j=1,L2
				Network(rowi,j)=AuxDownOrRight(ith+1,j)
			end do
		end if
		ValueH= DirectContractNetwork_V(Network)
		ValueH=ValueH/norm
		return
	end function

	function ValueV(ith,jth,H,singleA,DoubleA,norm,phya,phyb)
		real*8::ValueV
		integer,intent(in)::ith,jth
		character(len=*),intent(in)::phya,phyb
		type(SymTensor),intent(in)::H,singleA(:,:),DoubleA(:,:)
		real*8,intent(in)::norm
		type(SymTensor)::A1,A2
		integer::i,L1,L2,col,coli
		type(SymTensor),allocatable::Network(:,:)
		if(AuxFlag.ne.2)then
			call writemess('ERROR AuxFlag='+AuxFlag)
			call writemess(' please reset AuxTensor')
			call error_stop
		end if
		L1=size(singleA,1)
		L2=size(singleA,2)
		col=1
		if(jth.gt.1)col=col+1
		if(jth.lt.L2)col=col+1
		allocate(Network(L1,col))
		call Phi_H_Phi_two_site(A1,A2,singleA(ith,jth),singleA(ith+1,jth),'D','U',phya,phyb,H)
		coli=0
		if(jth.gt.1)then
			coli=coli+1
			do i=1,L1
				Network(i,coli)=AuxUpOrLeft(i,jth-1)
			end do
		end if
		coli=coli+1
		do i=1,L1
			if(i.eq.ith)then
				Network(i,coli)=A1
			else if(i.eq.(ith+1))then
				Network(i,coli)=A2
			else
				Network(i,coli)=DoubleA(i,jth)
			end if
		end do
		if(jth.lt.L2)then
			coli=coli+1
			do i=1,L1
				Network(i,coli)=AuxDownOrRight(i,jth+1)
			end do
		end if
		ValueV= DirectContractNetwork_H(Network)
		ValueV=ValueV/norm
		return
	end function


	function norm2_network()
		real*8::norm2_network
		integer::L1,L2,j,ijth,i
		type(SymTensor),allocatable::Network(:,:)
		if(AuxFlag.eq.1)then
			L1=size(AuxUpOrLeft,1)
			L2=size(AuxUpOrLeft,2)
			allocate(Network(2,L2))
			ijth=L1/2
			do j=1,L2
				Network(1,j)=AuxUpOrLeft(ijth,j)
				Network(2,j)=AuxDownOrRight(ijth+1,j)
			end do
			norm2_network = DirectContractNetwork_V(Network)
			return
		end if
		if(AuxFlag.eq.2)then
			L1=size(AuxUpOrLeft,1)
			L2=size(AuxUpOrLeft,2)
			allocate(Network(L1,2))
			ijth=L2/2
			do i=1,L1
				Network(i,1)=AuxUpOrLeft(i,ijth)
				Network(i,2)=AuxDownOrRight(i,ijth+1)
			end do
			norm2_network = DirectContractNetwork_H(Network)
			return
		end if
		call writemess('ERROR in norm2, do not calculate any auxTensor yet')
		call error_stop
	end function


	function EnergyFunction(singleA,DoubleA,inH,Dc)
		real*8::EnergyFunction
		type(SymTensor)::inH
		type(SymTensor)::singleA(:,:)
		type(SymTensor)::DoubleA(:,:)
		integer::Dc
		real*8::norm2,Energy_i
		integer::i,j,L1,L2,a1,a2,totalTerm
		logical::Flag(4)

		L1=size(singleA,1)
		L2=size(singleA,2)
		totalTerm=L1*(L2-1)+L2*(L1-1)
		call writemess(' initializtion...')
		allocate(AuxUpOrLeft(L1,L2))
		allocate(AuxDownOrRight(L1,L2))
		call initialUpDown(DoubleA,Dc)
		norm2=norm2_network()
		call writemess(' initializtion...done')
		call reset_time_calculator(totalTerm,30)
		EnergyFunction=0
 		do i=1,L1
 			do j=1,L2-1
				Energy_i=ValueH(i,j,inH,singleA,DoubleA,norm2,'n','n')
				EnergyFunction=EnergyFunction+Energy_i
				call time_calculator()
 			end do
 		end do
		call initialLeftRight(DoubleA,Dc)
		do j=1,L2
			do i=1,L1-1
				Energy_i=ValueV(i,j,inH,singleA,DoubleA,norm2,'n','n')
				EnergyFunction=EnergyFunction+Energy_i
				call time_calculator()
			end do
		end do
		EnergyFunction=EnergyFunction/(L1*L2)
	end function



	function ContractEnergy(nodeA,H,Dc)
		real*8::ContractEnergy
		type(Node),intent(inout)::nodeA(:,:)
		type(SymTensor),intent(inout)::H
		integer,intent(in)::Dc
		type(SymTensor),allocatable::A(:,:),PEPS(:,:)
		integer::L1,L2
		L1=size(nodeA,1)
		L2=size(nodeA,2)
		allocate(A(L1,L2))
		allocate(PEPS(L1,L2))
		call contractRightDown(nodeA,PEPS)
		call contract_physical_legs(A,PEPS)
		ContractEnergy= EnergyFunction(PEPS,A,H,Dc)
		return
	end function

	
end module











