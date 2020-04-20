module initial_module
	use Tensor_type
	use FermiTensor
	use SymDimension_typede
	use QuantumNumber_Type 
	use U1Tool
	use Tools
	implicit none
	type Node
		type(fTensor)::site
		type(fTensor)::Up
		type(fTensor)::Down
		type(fTensor)::Left
		type(fTensor)::Right
	end type Node

	!The U(1) rule:
		!  The base of spinless fermi are |0>       : no particle
		!                                 |1>       : there is one particle 
		!                                 
		!
		!  written in U1 base( |particle number, degeneracy> ) as
		!                               |0,1>=|0>         -----|1>
		!                               |1,1>=|1>         -----|2>
		!
		!   Then the transformational matrix that transform the basis form
		! spinless fermi to the U(1) is
		!
		!  / 1  0 \
		!  |      |
		!  \ 0  1 /
		!

contains

	type(fTensor) function C_Operator(dag)result(C)
		character(len=*),intent(in)::dag
		type(QuanNum)::QN(2)
		QN(1)=U1quantumnumber([0.,1.],[1,1],1,1)
		QN(2)=U1quantumnumber([0.,1.],[1,1],-1,-1)
		call C%allocate(QN,'real*8')
		if(dag.eq.'+')then
			call C%setValue((/2,1/),Tensor(1))!C^+
		else
			call C%setValue((/1,2/),Tensor(1))
		end if
		call C%setName(1,'H.n')
		call C%setName(2,'H.n2')
		return
	end function

	!initial_QN: 
		!  output the quantum numbers, degeneracies and the U1 symmetry rules of the legs in the tensors.
		!These date are not the only option. One can set any value of these data once there are non-zero
		!blocks in the tensor network.
		!
		!QNData(1) stores the Quantum numbers of the horizontal legs
		!QNData(2) stores the Quantum numbers of the verticl legs
		!QNData(3) stores the location of the legs that specify the total particle number

	subroutine initial_QN(QNData,Qspin,M,N,total,deg,len_of_QN,totalob)
		type(Tensor),intent(inout)::QNData(3)
		type(QuanNum),intent(inout)::Qspin
		integer::M,N,total,deg,len_of_QN,deltak,totalob
		integer::i,j,k,sign,ith,jth,inum,num
		integer,allocatable::A(:,:),B(:,:)
		real*8::delta
		logical::goon
		call writemess('Initial Quantum Number:')
		call writemess(' Total Paritcle number='+total)
		allocate(A(M,N-1))
		allocate(B(M-1,N))
		A=0
		B=0
		k=0
		deltak=(total/(M*N))+1
		if(total.eq.(totalob*M*N))deltak=totalob*2
		sign=1
		j=0
		i=1
		goon=.true.
		do while (goon)
				
				j=j+sign
				if(j.eq.N)then
					if(i.eq.M)then
						write(*,*)"ERROR"
						call error_stop
					end if
					k=k+deltak
					B(i,j)=k
					sign=-1*sign
					i=i+1
					if(k.eq.total-deltak)then
						ith=i
						jth=j
					end if
				else if((j.eq.1).and.(i.ne.1))then
					if(i.eq.M)then
						write(*,*)"ERROR"
						call error_stop
					end if
					k=k+deltak
					A(i,j)=k*sign
					if(k.eq.total-deltak)then
						ith=i
						jth=j
					else
						k=k+deltak
						B(i,j)=k
						sign=-1*sign
						i=i+1
						if(k.eq.total-deltak)then
							ith=i
							jth=j
						else
							k=k+deltak
							A(i,j)=k*sign
							if(k.eq.total-deltak)then
								ith=i
								jth=j+1
							end if
						end if
					end if
				else 
					k=k+deltak
					A(i,j)=k*sign
					if(k.eq.total-deltak)then
						if(sign.gt.0)then
							ith=i
							jth=j+1
						else
							ith=i
							jth=j
						end if
					end if
				end if
				if(k.ge.total-deltak)then
					goon=.false.
				end if
		end do

		Qspin=U1quantumnumber([0.,1.],[1,1],1,1)
		call Qspin%setRule(1)
		call Qspin%setFermiArrow(1)
		QNData(1)=A
		QNData(2)=B
		QNData(3)=[ith,jth,Total,len_of_QN,deg]
		return
	end subroutine 

	!initial_One_Node
		!   OBC tensor network state, the sites in the left top corner, right top corner, left below corner
		!and the right below corner have 3 legs, while the tensors in the boundary have 4 legs, the other 
		!tensor have 5 legs.
		

	type(Node) function initial_One_Node(flag,QNData,dspin,ith,jth)result(init)
		logical,intent(in)::flag(4)
		type(Tensor),intent(in)::QNData(3)
		type(QuanNum),intent(in)::dspin
		integer,intent(in)::ith,jth
		character(len=50),allocatable::dimname(:)
		integer,allocatable::rule(:)
		integer,allocatable::Arrow(:)
		type(QuanNum),allocatable::QN(:)
		integer::i,rank,degi,k
		integer,allocatable::deg(:)
		real*4::rnum,TotalN,QNnumLen
		logical::TotalNFlag
		rank=1
		i=1
		QNnumLen=QNData(3)%si(4)
		degi=QNData(3)%ii(5)
		allocate(deg(2*int(QNnumLen)+1))
		TotalNFlag=.false.
		if( (ith.eq.QNData(3)%ii(1)).and.(jth.eq.QNData(3)%ii(2)))then
			rank=rank+1
			TotalNFlag=.true.
		end if
		if(flag(1))rank=rank+1
		if(flag(2))rank=rank+1
		if(flag(3))rank=rank+1
		if(flag(4))rank=rank+1
		allocate(dimname(rank))
		allocate(QN(rank))
		allocate(rule(rank))
		allocate(Arrow(rank))
		if(TotalNFlag)then
			TotalN=QNData(3)%si(3)
			dimname(i)='A'+ith+'_'+jth+'.TotalN'
			QN(i)=U1quantumnumber(TotalN)
			call QN(i)%setDeg([1])
			call QN(i)%setRule(-1)
			call QN(i)%setFermiArrow(1)
			i=i+1
		end if
		if(flag(1))then
			dimname(i)='A'+ith+'_'+jth+'.L'
			rnum=QNData(1)%si([ith,jth-1])
			if(rnum.gt.0)then
				rule(i)=1
			else
				rule(i)=-1
			end if
			Arrow(i)=1
			rnum=abs(rnum)
			deg=degi
			QN(i)=U1quantumnumber(rnum-QNnumLen,rnum+QNnumLen,deg,Rule(i),Arrow(i))
			call init%Left%allocate([QN(i),QN(i)],'real*8')
			call init%Left%setName(1,'Lambda.L')
			call init%Left%setName(2,'Lambda.R')
			call init%Left%setRule((/rule(i),-1*rule(i)/))
			call init%Left%setFermiArrow((/1,-1/))
			call init%Left%eye()
			i=i+1
		else
			call init%Left%allocate((/1,1/),'real*8')
			call init%Left%setName("Lambda_Start")
		end if

		if(flag(2))then
			dimname(i)='A'+ith+'_'+jth+'.D'
			rnum=QNData(2)%si([ith,jth])
			if(rnum.gt.0)then
				rule(i)=-1
			else
				rule(i)=1
			end if
			Arrow(i)=-1
			rnum=abs(rnum)
			deg=degi
			QN(i)=U1quantumnumber(rnum-QNnumLen,rnum+QNnumLen,deg,Rule(i),Arrow(i))
			call init%Down%allocate([QN(i),QN(i)],'real*8')
			call init%Down%setName(1,'Lambda.U')
			call init%Down%setName(2,'Lambda.D')
			call init%Down%setRule((/-1*rule(i),rule(i)/))
			call init%Down%setFermiArrow((/1,-1/))
			call init%Down%eye()
			i=i+1
		else
			call init%Down%allocate((/1,1/),'real*8')
			call init%Down%setName("Lambda_Start")
		end if

		if(flag(3))then
			dimname(i)='A'+ith+'_'+jth+'.R'
			rnum=QNData(1)%si([ith,jth])
			if(rnum.gt.0)then
				rule(i)=-1
			else
				rule(i)=1
			end if
			Arrow(i)=-1
			rnum=abs(rnum)
			deg=degi
			QN(i)=U1quantumnumber(rnum-QNnumLen,rnum+QNnumLen,deg,Rule(i),Arrow(i))
			call init%Right%allocate([QN(i),QN(i)],'real*8')
			call init%Right%setName(1,'Lambda.L')
			call init%Right%setName(2,'Lambda.R')
			call init%Right%setRule((/-1*rule(i),rule(i)/))
			call init%Right%setFermiArrow((/1,-1/))
			call init%Right%eye()
			i=i+1
		else
			call init%Right%allocate((/1,1/),'real*8')
			call init%Right%setName("Lambda_Start")
		end if


		if(flag(4))then
			dimname(i)='A'+ith+'_'+jth+'.U'
			rnum=QNData(2)%si([ith-1,jth])
			if(rnum.gt.0)then
				rule(i)=1
			else
				rule(i)=-1
			end if
			Arrow(i)=1
			rnum=abs(rnum)
			deg=degi
			QN(i)=U1quantumnumber(rnum-QNnumLen,rnum+QNnumLen,deg,Rule(i),Arrow(i))
			call init%Up%allocate([QN(i),QN(i)],'real*8')
			call init%Up%setName(1,'Lambda.U')
			call init%Up%setName(2,'Lambda.D')
			call init%Up%setRule((/rule(i),-1*rule(i)/))
			call init%Up%setFermiArrow((/1,-1/))
			call init%Up%eye()
			i=i+1
		else
			call init%Up%allocate((/1,1/),'real*8')
			call init%Up%setName("Lambda_Start")
		end if

		dimname(i)='A'+ith+'_'+jth+'.n'
		QN(i)=dspin
		call init%site%allocate(QN,'real*8')
		call init%site%SymRandom((/-1d0,1d0/))
		do i=1,rank
			call init%site%setName(i,dimname(i))
		end do
		return
	end function


	subroutine initialPEPS(A,L1,L2,initialDeg,totalN)
		type(Node),allocatable,intent(inout)::A(:,:)
		integer,intent(in)::L1,L2,initialDeg,totalN
		integer::i,k,j,QNLen
		type(Tensor)::QNData(3)
		type(QuanNum)::spinQN
		character(len=200)::TensorName
		logical::flag(4)
		allocate(A(L1,L2))
		QNLen=1
		if(totalN.lt.0)then
			call writemess('ERROR in TotalN')
			call error_stop
		end if
		call initial_QN(QNData,spinQN,L1,L2,totalN,initialDeg,QNLen,1)
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
				A(i,j)=initial_One_Node(flag,QNData,spinQN,i,j)
			end do
		end do

		return
	end subroutine

	

end module
