module SimpleUpdate
	use Tensor_type
	use initial_module
	implicit none



contains

	
	subroutine initialH(H)
		type(Tensor),intent(inout)::H
		type(Tensor)::Sx,Sy,Sz
		type(Tensor)::S1x,S2x,S1y,S2y,S1z,S2z
		CALL H%setType('real*8')
		call pauli_matrix(S1x,S1y,S1z,0.5d0)
		call pauli_matrix(S2x,S2y,S2z,0.5d0)
		call S1x%setName(1,'H1.n')
		call S1x%setName(2,'H1.n2')
		call S2x%setName(1,'H2.n')
		call S2x%setName(2,'H2.n2')
		
		call S1y%setName(1,'H1.n')
		call S1y%setName(2,'H1.n2')
		call S2y%setName(1,'H2.n')
		call S2y%setName(2,'H2.n2')
		
		call S1z%setName(1,'H1.n')
		call S1z%setName(2,'H1.n2')
		call S2z%setName(1,'H2.n')
		call S2z%setName(2,'H2.n2')
		
		H=(S1x.kron.S2x) + (S1y.kron.S2y) + (S1z.kron.S2z)
		call H%permute(['H1.n ','H2.n ','H1.n2','H2.n2'])

		return
	end subroutine

	type(Tensor) function gate(H,tau)
		type(Tensor),intent(in)::H
		real*8,intent(in)::tau
		gate=H*tau
		call gate%fuse(1,2)
		call gate%fuse(2,3)
		gate=expm(gate)
		call gate%split()
		return
	end function

	type(Tensor) function invLabmbda(Lambda)
		type(Tensor),intent(in)::Lambda
		integer::i,lenLambda,j
		type(Tensor)::temp
		call invLabmbda%setType(Lambda%getType())
		call temp%setType(Lambda%getType())
		invLabmbda=Lambda
		do i=1,invLabmbda.dim.1
			call invLabmbda%setValue([i,i],1d0/Lambda%di([i,i]))
		end do
		return
	end function

	subroutine step(A,lambda,B,R,L,expmH,CutOff)!Use QR and LQ decomposion,Can not work on MPS
		type(Tensor),intent(inout)::A,B,lambda
		type(Tensor),intent(in)::expmH
		character(len=*),intent(in)::R,L
		character(len=50)::Aname,Bname,AphyName,BphyName,AlegName,BlegName
		integer,intent(in)::CutOff
		character(len=50)::lambdaL,lambdaR,phya,phyb
		type(Tensor)::C,USV(3),QR(2),LQ(2)
		phya='n'
		phyb='n'
		if(expmH%getRank().ne.4)then
			call writemess('ERROR in stepGate, rank of the expm(tau*H), in step')
			call error_stop
		endif
		lambdaL=lambda%getName(L)
		lambdaR=lambda%getName(R)
		AlegName=A%getName(R)
		BlegName=B%getName(L)
		Aname=AlegName.subl.'.'
		Bname=BlegName.subl.'.'
		AphyName=A%getName(phya)
		BphyName=B%getName(phyb)
		call A%QR(QR(1),QR(2),[AphyName,AlegName],.false.)
		if(QR(2)%getRank().ne.3)then
			call writemess("ERROR in QR in simple update")
			call QR(2)%diminfo()
			call error_stop
		end if
		call QR(1)%setName('QR.Q',AlegName)
		call QR(2)%setName('QR.R',Aname+'.'+L)
		
		
		call C%contract(QR(2),AlegName,invLabmbda(lambda),lambdaL)
		
		call B%LQ(LQ(1),LQ(2),[BlegName,BphyName])
		if(LQ(1)%getRank().ne.3)then
			call writemess("ERROR in LQ in simple update")
			call LQ(1)%diminfo()
			call error_stop
		end if
		call LQ(1)%setName('LQ.L',Bname+'.'+R)
		call LQ(2)%setName('LQ.Q',BlegName)

		
		
		call C%contract(lambdaR,LQ(1),BlegName)
		call C%contract(expmH,['H1.n2','H2.n2'],[AphyName,BphyName])
		call C%setName('H1.n',AphyName)
		call C%setName('H2.n',BphyName)
		
		call C%SVD(USV(1),USV(2),USV(3),Aname,Bname,CutOff)
		USV(2)=eye(USV(2)/USV(2)%dmax('maxa'))
		call USV(2)%setName(1,'SVD.s1')
		call USV(2)%setName(2,'SVD.s2')

		lambda=USV(2)
		call A%contract(USV(1),'SVD.U',lambda,'SVD.s1')
		call A%contract(QR(1),Aname+'.'+R,Aname+'.'+L)
		call A%setName('SVD.s2',Aname+'.'+R)
		
		call B%contract(lambda,'SVD.s2',USV(3),'SVD.V')
		call B%contract(Bname+'.'+R,LQ(2),Bname+'.'+L)
		call B%setName('SVD.s1',Bname+'.'+L)

		call lambda%setName(1,lambdaL)
		call lambda%setName(2,lambdaR)
		return
	end subroutine	

	subroutine sampleUpdate(A,expH,CutOff)!A,expHType1,CutOff,MPIINFO)!each cpus update one row or col
		type(node),intent(inout)::A(:,:)
		type(Tensor),intent(inout)::expH
		integer,intent(inout)::CutOff
		integer::i,j,L1,L2
		L1=size(A,1)
		L2=size(A,2)
		do i=1,L1
			do j=1,L2-1
				call step(A(i,j)%site,A(i,j)%Right,A(i,j+1)%site,'R','L',expH,CutOff)
				A(i,j+1)%Left=A(i,j)%Right
			end do
		end do
		do j=1,L2
			do i=1,L1-1
				call step(A(i,j)%site,A(i,j)%Down,A(i+1,j)%site,'D','U',expH,CutOff)
				A(i+1,j)%Up=A(i,j)%Down
			end do
		end do
		return
	end subroutine
	

	type(Tensor) function invsqrt(T)
		type(Tensor),intent(in)::T
		type(Tensor)::Temp
		integer::i,j
		if(T%getRank().ne.2)then
			call writemess('ERROR in invsqrt,fPEPS.f90')
			call error_stop
		end if
		invsqrt=T
		do i=1,invsqrt.dim.1
			call invsqrt%setValue([i,i],1d0/dsqrt(T%di([i,i])))
		end do
		return
	end function
	
	subroutine contractRightDown(A,B)
		type(node),intent(in)::A(:,:)
		type(Tensor),intent(inout)::B(:,:)
		integer::i,j
		character(len=50)::Aname,newname
		do i=1,size(A,1)
			do j=1,size(A,2)
				Aname=A(i,j)%site%outTensorName(1)
				if(A(i,j)%site%NameOrder(Aname+'.L').ne.0) then
					call B(i,j)%contract(invsqrt(A(i,j)%Left),A(i,j)%Left%getName('R'),A(i,j)%site,Aname+'.L')
					call B(i,j)%setName(A(i,j)%Left%getName('L'),Aname+'.L')
				else
					B(i,j)=A(i,j)%site
				end if
				
				if(A(i,j)%site%NameOrder(Aname+'.D').ne.0) then
					call B(i,j)%contract(Aname+'.D',invsqrt(A(i,j)%Down),A(i,j)%Down%getName('U'))
					call B(i,j)%setName(A(i,j)%Down%getName('D'),Aname+'.D')
				end if
				
				if(A(i,j)%site%NameOrder(Aname+'.R').ne.0) then
					call B(i,j)%contract(Aname+'.R',invsqrt(A(i,j)%Right),A(i,j)%Right%getName('L'))
					call B(i,j)%setName(A(i,j)%Right%getName('R'),Aname+'.R')
				end if
				
				if(A(i,j)%site%NameOrder(Aname+'.U').ne.0) then
					call B(i,j)%contract(Aname+'.U',invsqrt(A(i,j)%Up),A(i,j)%Up%getName('D'))
					call B(i,j)%setName(A(i,j)%Up%getName('U'),Aname+'.U')
				end if
			end do
		end do
		return
	end subroutine





end module
