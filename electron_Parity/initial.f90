module initial_module
	use Tensor_type
	use FermiTensor
	use SymDimension_typede
	use QuantumNumber_Type 
	use ParityTool
	use Tools
	implicit none
	type Node
		type(fTensor)::site
		type(fTensor)::Up
		type(fTensor)::Down
		type(fTensor)::Left
		type(fTensor)::Right
	end type Node

	!The Parity rule:
		!  The base of spinless fermi are |0>       : no particle
		!                                 |1>       : there is one particle 
		!                                 
		!
		!  written in Parity base( |qunatum number, degeneracy> ) as
		!                               |-1,1>=|1>         -----|1>
		!                               |+1,1>=|0>         -----|2>
		!
		!   Then the transformational matrix that transform the basis form
		! spinless fermi to the parity is
		!
		!  / 0  1 \
		!  |      |
		!  \ 1  0 /


contains

	type(fTensor) function C_Operator(dag)result(C)
		character(len=*),intent(in)::dag
		type(QuanNum)::QN(2)
		QN(1)=ParityQuantumNumber([1,1])
		call QN(1)%setFermiArrow(1)
		QN(2)=ParityQuantumNumber([1,1])
		call QN(2)%setFermiArrow(-1)
		call C%allocate(QN,'real*8')
		if(dag.eq.'+')then
			call C%setValue((/1,2/),Tensor(1))
		else
			call C%setValue((/2,1/),Tensor(1))
		end if
		call C%setName(1,'H.n')
		call C%setName(2,'H.n2')
		return
	end function

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
		integer,allocatable::FermiArrow(:)
		type(QuanNum),allocatable::QN(:)
		integer::i,j,k
		rank=1
		if(flag(1))rank=rank+1
		if(flag(2))rank=rank+1
		if(flag(3))rank=rank+1
		if(flag(4))rank=rank+1
		allocate(dimname(rank))
		allocate(QN(rank))
		allocate(FermiArrow(rank))
		i=1
		if(flag(1))then
			dimname(i)=nameA+'.L'
			QN(i)=PEPSD
			FermiArrow(i)=1
			i=i+1
			call init%Left%allocate([PEPSD,PEPSD],'real*8')
			call init%Left%eye()
			call init%Left%setName(1,'Lambda.L')
			call init%Left%setName(2,'Lambda.R')
			call init%Left%setFermiArrow((/1,-1/))
		else
			call init%Left%allocate((/1,1/),'real*8')
			call init%Left%setName("Lambda_Start")
		end if
		if(flag(2))then
			dimname(i)=nameA+'.D'
			QN(i)=PEPSD
			FermiArrow(i)=-1
			i=i+1
			call init%Down%allocate([PEPSD,PEPSD],'real*8')
			call init%Down%eye()
			call init%Down%setName(1,'Lambda.U')
			call init%Down%setName(2,'Lambda.D')
			call init%Down%setFermiArrow((/1,-1/))
		else
			call init%Down%allocate((/1,1/),'real*8')
			call init%Down%setName("Lambda_Start")
		end if
		if(flag(3))then
			dimname(i)=nameA+'.R'
			QN(i)=PEPSD
			FermiArrow(i)=-1
			i=i+1
			call init%Right%allocate([PEPSD,PEPSD],'real*8')
			call init%Right%eye()
			call init%Right%setName(1,'Lambda.L')
			call init%Right%setName(2,'Lambda.R')
			call init%Right%setFermiArrow((/1,-1/))
		else
			call init%Right%allocate((/1,1/),'real*8')
			call init%Right%setName("Lambda_Start")
		end if
		if(flag(4))then
			dimname(i)=nameA+'.U'
			QN(i)=PEPSD
			FermiArrow(i)=1
			i=i+1
			call init%Up%allocate([PEPSD,PEPSD],'real*8')
			call init%Up%eye()
			call init%Up%setName(1,'Lambda.U')
			call init%Up%setName(2,'Lambda.D')
			call init%Up%setFermiArrow((/1,-1/))
		else
			call init%Up%allocate((/1,1/),'real*8')
			call init%Up%setName("Lambda_Start")
		end if

		dimname(i)=nameA+'.n'
		QN(i)=dspin
		
		FermiArrow(i)=1
		call init%site%allocate(QN,'real*8')
		call init%site%SymRandom((/-1d0,1d0/))
		do i=1,rank
			call init%site%setName(i,dimname(i))
		end do
		call init%site%setFermiArrow(FermiArrow)
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
		call spinQN%setFermiArrow(1)
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
