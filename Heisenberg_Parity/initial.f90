module initial_module
	use Tensor_type
	use SymTensor_type
	use SymDimension_typede
	use QuantumNumber_Type 
	use ParityTool
	use Tools
	implicit none
	type Node
		type(SymTensor)::site
		type(SymTensor)::Up
		type(SymTensor)::Down
		type(SymTensor)::Left
		type(SymTensor)::Right
	end type Node

	!The Parity rule:
		!  The base of spin are |up>       
		!                       |down>    
		!                                 
		! The number of spin that point to up is even, then the Parity=+1
		! The number of spin that point to up is odd,  then the Parity=-1
		!
		!  in Parity base( |qunatum number, degeneracy> ) as
		!                               |-1,1>=|up>           -----|1>
		!                               |+1,1>=|down>         -----|2>
		!
		!   Then the transformational matrix that transform the basis form
		! spin to the parity is
		!
		!  / 1  0 \
		!  |      |
		!  \ 0  1 /
		!
		!


contains

	subroutine spin_Operator(Sx,Sy,Sz,TensorName)
		type(SymTensor),intent(inout)::Sx,Sy,Sz
		character(len=*),intent(in)::TensorName
		type(QuanNum)::QN(2)
		QN(1)=ParityQuantumNumber([1,1])
		QN(2)=ParityQuantumNumber([1,1])
		call Sx%allocate(QN,'real*8')
		call Sy%allocate(QN,'complex*16')
		call Sz%allocate(QN,'real*8')
		call Sx%setValue([1,2],Tensor(0.5))
		call Sx%setValue([2,1],Tensor(0.5))
		call Sx%setName(1,TensorName+'.n')
		call Sx%setName(2,TensorName+'.n2')

		call Sy%setValue([1,2],Tensor(dcmplx(0d0,-0.5d0)))
		call Sy%setValue([2,1],Tensor(dcmplx(0d0,0.5d0)))
		call Sy%setName(1,TensorName+'.n')
		call Sy%setName(2,TensorName+'.n2')

		call Sz%setValue([1,1],Tensor(0.5))
		call Sz%setValue([2,2],Tensor(-0.5))
		call Sz%setName(1,TensorName+'.n')
		call Sz%setName(2,TensorName+'.n2')
		return
	end subroutine

	!initial_One_Node
		!   OBC tensor network state, the sites in the left top corner, right top corner, left below corner
		!and the right below corner have 3 legs, while the tensors in the boundary have 4 legs, the other 
		!tensor have 5 legs.

	type(Node) function initial_One_Node(flag,PEPSD,dspin,nameA)result(init)
		logical,intent(in)::flag(4)
		type(QuanNum),intent(in)::dspin,PEPSD
		character(len=*),intent(in)::nameA
		integer::rank
		character(len=50),allocatable::dimname(:)
		type(QuanNum),allocatable::QN(:)
		integer::i,j,k
		rank=1
		if(flag(1))rank=rank+1
		if(flag(2))rank=rank+1
		if(flag(3))rank=rank+1
		if(flag(4))rank=rank+1
		allocate(dimname(rank))
		allocate(QN(rank))
		i=1
		if(flag(1))then
			dimname(i)=nameA+'.L'
			QN(i)=PEPSD
			i=i+1
			call init%Left%allocate([PEPSD,PEPSD],'real*8')
			call init%Left%eye()
			call init%Left%setName(1,'Lambda.L')
			call init%Left%setName(2,'Lambda.R')
		else
			call init%Left%allocate((/1,1/),'real*8')
			call init%Left%setName("Lambda_Start")
		end if
		if(flag(2))then
			dimname(i)=nameA+'.D'
			QN(i)=PEPSD
			i=i+1
			call init%Down%allocate([PEPSD,PEPSD],'real*8')
			call init%Down%eye()
			call init%Down%setName(1,'Lambda.U')
			call init%Down%setName(2,'Lambda.D')
		else
			call init%Down%allocate((/1,1/),'real*8')
			call init%Down%setName("Lambda_Start")
		end if
		if(flag(3))then
			dimname(i)=nameA+'.R'
			QN(i)=PEPSD
			i=i+1
			call init%Right%allocate([PEPSD,PEPSD],'real*8')
			call init%Right%eye()
			call init%Right%setName(1,'Lambda.L')
			call init%Right%setName(2,'Lambda.R')
		else
			call init%Right%allocate((/1,1/),'real*8')
			call init%Right%setName("Lambda_Start")
		end if
		if(flag(4))then
			dimname(i)=nameA+'.U'
			QN(i)=PEPSD
			i=i+1
			call init%Up%allocate([PEPSD,PEPSD],'real*8')
			call init%Up%eye()
			call init%Up%setName(1,'Lambda.U')
			call init%Up%setName(2,'Lambda.D')
		else
			call init%Up%allocate((/1,1/),'real*8')
			call init%Up%setName("Lambda_Start")
		end if

		dimname(i)=nameA+'.n'
		QN(i)=dspin
		
		call init%site%allocate(QN,'real*8')
		call init%site%SymRandom((/-1d0,1d0/))
		do i=1,rank
			call init%site%setName(i,dimname(i))
		end do
		return
	end function

	subroutine initialPEPS(A,L1,L2,PEPS_D1,PEPS_D2)
		type(Node),allocatable,intent(inout)::A(:,:)
		integer,intent(in)::PEPS_D1,PEPS_D2,L1,L2
		type(QuanNum)::spinQN,QND
		integer::i,k,j
		logical::flag(4)
		character(len=100)::name1
		allocate(A(L1,L2))
		QND=ParityQuantumNumber([PEPS_D1,PEPS_D2])
		spinQN=ParityQuantumNumber([1,1])
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
				name1='A'+i+'_'+j!Ai_j
				A(i,j)=initial_One_Node(Flag,QND,spinQN,name1)
			end do
		end do
		return
	end subroutine
	

end module
