module initial_module
	use Tensor_type
	use SymTensor_type
	use SymDimension_typede
	use QuantumNumber_Type 
	use U1Tool
	use Tools
	implicit none
	type Node
		type(SymTensor)::site
		type(SymTensor)::Up
		type(SymTensor)::Down
		type(SymTensor)::Left
		type(SymTensor)::Right
	end type Node

	!The U(1) rule:
		!  The base of spin are |up>       
		!                       |down>    
		!                                 
		!
		!  in U(1) base( |Sz, degeneracy> ) as
		!                            |-0.5,1>=|down>       -----|1>
		!                            | 0.5,1>=|up>         -----|2>
		!
		!   Then the transformational matrix that transform the basis form
		! spin to the parity is
		!
		!  / 0  1 \
		!  |      |
		!  \ 1  0 /
		!
		!


contains

	subroutine spin_Operator(Sx,Sy,Sz,TensorName)
		type(SymTensor),intent(inout)::Sx,Sy,Sz
		character(len=*),intent(in)::TensorName
		type(QuanNum)::QN(2)
		call QN(1)%setQN([-0.5,0.5])
		call QN(1)%setDeg([1,1])
		call QN(1)%setRule(1)

		call QN(2)%setQN([-0.5,0.5])
		call QN(2)%setDeg([1,1])
		call QN(2)%setRule(-1)

		call Sx%allocate(QN,'real*8')
		call Sy%allocate(QN,'complex*16')
		call Sz%allocate(QN,'real*8')

		call Sx%setValue([1,2],Tensor(0.5))
		call Sx%setValue([2,1],Tensor(0.5))
		call Sx%setName(1,TensorName+'.n')
		call Sx%setName(2,TensorName+'.n2')

		call Sy%setValue([1,2],Tensor(dcmplx(0d0,0.5d0)))
		call Sy%setValue([2,1],Tensor(dcmplx(0d0,-0.5d0)))
		call Sy%setName(1,TensorName+'.n')
		call Sy%setName(2,TensorName+'.n2')

		call Sz%setValue([1,1],Tensor(-0.5))
		call Sz%setValue([2,2],Tensor(0.5))
		call Sz%setName(1,TensorName+'.n')
		call Sz%setName(2,TensorName+'.n2')
		return
	end subroutine


	!NetworkQN: 
		!  output the quantum numbers, degeneracies and the U1 symmetry rules of the legs in the tensors.
		!These date are not the only option. One can set any value of these data once there are non-zero
		!blocks in the tensor network.
		!
		!QNData(1) stores the quantum number of the legs
		!QNData(2) stores the degeneracy for each quantums
		!QNData(3) stores the symmetry rule of the legs

	function NetworkQN(ith,jth,L1,L2,Deg,leg)result(QNData)
		type(Tensor),allocatable::QNData(:)
		integer,intent(in)::ith,jth,L1,L2,Deg
		character(len=*),intent(in)::leg
		allocate(QNData(3))
		if(leg.equ.'R')then
			if(mod(jth,2).eq.0)then
				QNData(1)=[-1.,0.,1.]
				QNData(2)=[Deg,Deg,Deg]
			else
				QNData(1)=[-0.5,0.5]
				QNData(2)=[Deg,Deg]
			end if
			if(mod(ith,2).eq.0)then
				QNData(3)=1
			else
				QNData(3)=-1
			end if
		end if
		if(leg.equ.'L')then
			if(mod(jth,2).eq.0)then
				QNData(1)=[-0.5,0.5]
				QNData(2)=[Deg,Deg]
			else
				QNData(1)=[-1.,0.,1.]
				QNData(2)=[Deg,Deg,Deg]
			end if
			if(mod(ith,2).eq.0)then
				QNData(3)=-1
			else
				QNData(3)=1
			end if
		end if
		if(leg.equ.'D')then
			if((mod(L2,2).ne.0).and.(mod(ith,2).ne.0).and.(jth.eq.L2))then
				QNData(1)=[-0.5,0.5]
				QNData(2)=[Deg,Deg]
				QNData(3)=-1
			else
				QNData(1)=[-1.,0.,1.]
				QNData(2)=[Deg,Deg,Deg]
				QNData(3)=-1
			end if
		end if
		if(leg.equ.'U')then
			if((mod(L2,2).ne.0).and.(mod(ith,2).eq.0).and.(jth.eq.L2))then
				QNData(1)=[-0.5,0.5]
				QNData(2)=[Deg,Deg]
				QNData(3)=1
			else
				QNData(1)=[-1.,0.,1.]
				QNData(2)=[Deg,Deg,Deg]
				QNData(3)=1
			end if
		end if
		return
	end function

	!initial_One_Node
		!   OBC tensor network state, the sites in the left top corner, right top corner, left below corner
		!and the right below corner have 3 legs, while the tensors in the boundary have 4 legs, the other 
		!tensor have 5 legs.

	type(Node) function initial_One_Node(flag,deg,dspin,ith,jth,L1,L2)result(init)
		logical,intent(in)::flag(4)
		type(QuanNum),intent(in)::dspin
		integer,intent(in)::ith,jth,deg,L1,L2
		character(len=50),allocatable::dimname(:)
		integer,allocatable::rule(:)
		type(QuanNum),allocatable::QN(:)
		type(Tensor)::QNData(3)
		integer::i,rank,k
		rank=1
		i=1
		if(flag(1))rank=rank+1
		if(flag(2))rank=rank+1
		if(flag(3))rank=rank+1
		if(flag(4))rank=rank+1
		allocate(dimname(rank))
		allocate(QN(rank))
		allocate(rule(rank))
		if(flag(1))then
			dimname(i)='A'+ith+'_'+jth+'.L'
			QNData=NetworkQN(ith,jth,L1,L2,Deg,'L')
			rule(i)=QNData(3)
			call QN(i)%setQN(QNData(1)%si())
			call QN(i)%setDeg(QNData(2)%ii())
			call QN(i)%setRule(rule(i))
			call init%Left%allocate([QN(i),QN(i)],'real*8')
			call init%Left%setName(1,'Lambda.L')
			call init%Left%setName(2,'Lambda.R')
			call init%Left%setRule((/rule(i),-1*rule(i)/))
			call init%Left%eye()
			i=i+1
		else
			call init%Left%allocate((/1,1/),'real*8')
			call init%Left%setName("Lambda_Start")
		end if

		if(flag(2))then
			dimname(i)='A'+ith+'_'+jth+'.D'
			QNData=NetworkQN(ith,jth,L1,L2,Deg,'D')
			rule(i)=QNData(3)
			call QN(i)%setQN(QNData(1)%si())
			call QN(i)%setDeg(QNData(2)%ii())
			call QN(i)%setRule(rule(i))
			call init%Down%allocate([QN(i),QN(i)],'real*8')
			call init%Down%setName(1,'Lambda.U')
			call init%Down%setName(2,'Lambda.D')
			call init%Down%setRule((/-1*rule(i),rule(i)/))
			call init%Down%eye()
			i=i+1
		else
			call init%Down%allocate((/1,1/),'real*8')
			call init%Down%setName("Lambda_Start")
		end if

		if(flag(3))then
			dimname(i)='A'+ith+'_'+jth+'.R'
			QNData=NetworkQN(ith,jth,L1,L2,Deg,'R')
			rule(i)=QNData(3)
			call QN(i)%setQN(QNData(1)%si())
			call QN(i)%setDeg(QNData(2)%ii())
			call QN(i)%setRule(rule(i))
			call init%Right%allocate([QN(i),QN(i)],'real*8')
			call init%Right%setName(1,'Lambda.L')
			call init%Right%setName(2,'Lambda.R')
			call init%Right%setRule((/-1*rule(i),rule(i)/))
			call init%Right%eye()
			i=i+1
		else
			call init%Right%allocate((/1,1/),'real*8')
			call init%Right%setName("Lambda_Start")
		end if


		if(flag(4))then
			dimname(i)='A'+ith+'_'+jth+'.U'
			QNData=NetworkQN(ith,jth,L1,L2,Deg,'U')
			rule(i)=QNData(3)
			call QN(i)%setQN(QNData(1)%si())
			call QN(i)%setDeg(QNData(2)%ii())
			call QN(i)%setRule(rule(i))
			call init%Up%allocate([QN(i),QN(i)],'real*8')
			call init%Up%setName(1,'Lambda.U')
			call init%Up%setName(2,'Lambda.D')
			call init%Up%setRule((/rule(i),-1*rule(i)/))
			call init%Up%eye()
			i=i+1
		else
			call init%Up%allocate((/1,1/),'real*8')
			call init%Up%setName("Lambda_Start")
		end if

		dimname(i)='A'+ith+'_'+jth+'.n'
		QN(i)=dspin
		call init%site%allocate(QN,'real*8')
		do i=1,rank
			call init%site%setName(i,dimname(i))
		end do
		call init%site%SymRandom((/-1d0,1d0/))
		return
	end function


	subroutine initialPEPS(A,L1,L2,deg)
		type(Node),allocatable,intent(inout)::A(:,:)
		integer,intent(in)::L1,L2,deg
		integer::i,k,j
		type(QuanNum)::spinQN
		character(len=200)::TensorName
		logical::flag(4)
		allocate(A(L1,L2))
		call spinQN%setQN([-0.5,0.5])
		call spinQN%setDeg([1,1])
		call spinQN%setRule(1)
		do i=1,L1
			flag(4)=.true.!D4=PEPS_D
			flag(2)=.true.!D2=PEPS_D
			if(i.eq.1)then
				flag(4)=.false.
			end if
			if(i.eq.L1)then
				flag(2)=.false.!D2=0
			end if
			do j=1,L2
				flag(1)=.true.!D1=PEPS_D
				flag(3)=.true.!D3=PEPS_D
				if(j.eq.1)then
					flag(1)=.false.!D1=0
				end if
				if (j.eq.L2)then
					flag(3)=.false.!D3=0
				end if
				A(i,j)=initial_One_Node(flag,deg,spinQN,i,j,L1,L2)
				call A(i,j)%Site%SymCheck()
			end do
		end do

		return
	end subroutine
	

end module
