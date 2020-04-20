module SimpleUpdate
	use Tensor_type
	use Tools
	implicit none
	type Node
		type(Tensor)::site
		type(Tensor)::Up
		type(Tensor)::Down
		type(Tensor)::Left
		type(Tensor)::Right
	end type Node
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
		real*4,intent(in)::tau
		gate=H*tau
		call gate%fuse(1,2)
		call gate%fuse(2,3)
		gate=expm(gate)
		call gate%split()
		return
	end function

	subroutine initialPEPS(A,L1,L2,D)
		type(Node),allocatable,intent(inout)::A(:,:)
		integer,intent(in)::L1,L2,D
		integer::i,j
		allocate(A(L1,L2))
		do j=1,L2
			do i=1,L1
				call A(i,j)%site%allocate([D,D,D,D,2],'real*8')
				call A(i,j)%site%random([-1d0,1d0])
				call A(i,j)%site%setName(1,'A'+i+'_'+j+'.Left')
				call A(i,j)%site%setName(2,'A'+i+'_'+j+'.Down')
				call A(i,j)%site%setName(3,'A'+i+'_'+j+'.Right')
				call A(i,j)%site%setName(4,'A'+i+'_'+j+'.Up')
				call A(i,j)%site%setName(5,'A'+i+'_'+j+'.n')
				call A(i,j)%Up%allocate([D,D],'real*8')
				call A(i,j)%Up%eye()
				call A(i,j)%Up%setName(1,'Lambda.Leg1')
				call A(i,j)%Up%setName(2,'Lambda.Leg2')
				A(i,j)%Down=A(i,j)%Up
				A(i,j)%Left=A(i,j)%Up
				A(i,j)%Right=A(i,j)%Up
			end do
		end do
		return
	end subroutine


	subroutine step(A,Lambda,B,expmH,CutOff,AName,BName,cha1,cha2)
		type(Tensor),intent(inout)::A,Lambda,B
		type(Tensor),intent(in)::expmH
		character(len=*),intent(in)::cha1,cha2,AName,BName
		integer,intent(in)::CutOff
		type(Tensor)::C,SVD(3)
		C=contract(A,AName+'.'+cha1,Lambda%invTensor(),'Lambda.Leg1')
		C=contract(C,'Lambda.Leg2',B,BName+'.'+cha2)
		C=contract(C,[AName+'.n',BName+'.n'],expmH,['H1.n2','H2.n2'])
		call C%setName('H1.n',AName+'.n')
		call C%setName('H2.n',BName+'.n')
		SVD=C%SVDTensor(AName,BName,CutOff)
		Lambda=eye(SVD(2))/SVD(2)%smax()
		A=contract(SVD(1),'SVD.U',Lambda,1)
		call A%setName(A%getRank(),AName+'.'+cha1)
		B=contract(Lambda,2,SVD(3),'SVD.V')
		call B%setName(1,BName+'.'+cha2)
		call Lambda%setName(1,'Lambda.Leg1')
		call Lambda%setName(2,'Lambda.Leg2')
		return
	end subroutine

	subroutine sampleUpdate(A,expmH,L1,L2,CutOff)
		type(Node),intent(inout)::A(:,:)
		type(Tensor),intent(in)::expmH
		integer,intent(in)::L1,L2,CutOff
		integer::i,j
		do i=1,L1
			do j=1,L2-1
				call step(A(i,j)%Site,A(i,j)%Right,A(i,j+1)%Site,expmH,CutOff,'A'+i+'_'+j,'A'+i+'_'+(j+1),'Right','Left')
				A(i,j+1)%Left=A(i,j)%Right
			end do
			call step(A(i,L2)%Site,A(i,L2)%Right,A(i,1)%Site,expmH,CutOff,'A'+i+'_'+L2,'A'+i+'_1','Right','Left')
			A(i,1)%Left=A(i,L2)%Right
		end do
		do j=1,L2
			do i=1,L1-1
				call step(A(i,j)%Site,A(i,j)%Down,A(i+1,j)%Site,expmH,CutOff,'A'+i+'_'+j,'A'+(i+1)+'_'+j,'Down','Up')
				A(i+1,j)%Up=A(i,j)%Down
			end do
			call step(A(L1,j)%Site,A(L1,j)%Down,A(1,j)%Site,expmH,CutOff,'A'+L1+'_'+j,'A1_'+j,'Down','Up')
			A(1,j)%Up=A(L1,j)%Down
		end do
		return
	end subroutine


!************************** calculate the Energy ********************

	subroutine contractRightDownContract(A,B)
		type(node),intent(in)::A(:,:)
		type(Tensor),intent(inout)::B(:,:)
		integer::i,j
		character(len=50)::Aname
		do i=1,size(A,1)
			do j=1,size(A,2)
				B(i,j)=contract(A(i,j)%site,A(i,j)%site%getName('Right'),A(i,j)%Right%invTensor(),'Lambda.Leg1')
				call B(i,j)%setName('Lambda.Leg2',A(i,j)%site%getName('Right'))
				B(i,j)=contract(B(i,j),A(i,j)%site%getName('Down'),A(i,j)%Down%invTensor(),'Lambda.Leg1')
				call B(i,j)%setName('Lambda.Leg2',A(i,j)%site%getName('Down'))
			end do
		end do
	end subroutine


	type(Tensor) function contractNetWork(A) 
		type(Tensor),intent(in)::A(:,:)
		type(Tensor)::Res
		integer::i,j,leni,lenj,k
		character(len=50)::name1,name2,name3,name4,name5,name6,name7,name8
		leni=size(A,1)
		lenj=size(A,2)
		Res=A(1,1)
		do i=1,leni-1
			if(i.eq.1)then
				if(lenj.gt.1)then
					name1=Res%getName('Right')
					name2=A(i,2)%getName('Left')
					Res=contract(Res,name1,A(i,2),name2)
				end if
				do j=3,lenj-1
					name1=A(i,j-1)%getName('Right')
					name2=A(i,j)%getName('Left')
					Res=contract(Res,name1,A(i,j),name2)
				end do
				name1=A(i,lenj-1)%getName('Right')
				name2=A(i,1)%getName('Left')
				name3=A(i,lenj)%getName('Left')
				name4=A(i,lenj)%getName('Right')
				Res=contract(Res,[name1,name2],A(i,lenj),[name3,name4])
			else
				name1=A(i-1,1)%getName('Down')
				name2=A(i,1)%getName('Up')
				Res=contract(Res,name1,A(i,1),name2)
				do j=2,lenj-1
					name1=A(i-1,j)%getName('Down')
					name2=A(i,j-1)%getName('Right')
					name3=A(i,j)%getName('Up')
					name4=A(i,j)%getName('Left')
					Res=contract(Res,[name1,name2],A(i,j),[name3,name4])
				end do
				name1=A(i-1,lenj)%getName('Down')
				name2=A(i,lenj-1)%getName('Right')
				name3=A(i,1)%getName('Left')
				name4=A(i,lenj)%getName('Up')
				name5=A(i,lenj)%getName('Left')
				name6=A(i,lenj)%getName('Right')
				Res=contract(Res,[name1,name2,name3],A(i,lenj),[name4,name5,name6])
			end if
		end do
		i=leni
		name1=A(i-1,1)%getName('Down')
		name2=A(1,1)%getName('Up')
		name3=A(i,1)%getName('Up')
		name4=A(i,1)%getName('Down')
		Res=contract(Res,[name1,name2],A(i,1),[name3,name4])
		do j=2,lenj-1
			name1=A(i-1,j)%getName('Down')
			name2=A(i,j-1)%getName('Right')
			name3=A(1,j)%getName('Up')
			name4=A(i,j)%getName('Up')
			name5=A(i,j)%getName('Left')
			name6=A(i,j)%getName('Down')
			Res=contract(Res,[name1,name2,name3],A(i,j),[name4,name5,name6])
		end do
		name1=A(i-1,lenj)%getName('Down')
		name2=A(i,lenj-1)%getName('Right')
		name3=A(i,1)%getName('Left')
		name4=A(1,lenj)%getName('Up')
		name5=A(i,lenj)%getName('Up')
		name6=A(i,lenj)%getName('Left')
		name7=A(i,lenj)%getName('Right')
		name8=A(i,lenj)%getName('Down')
		Res=contract(Res,[name1,name2,name3,name4],A(i,lenj),[name5,name6,name7,name8])
		contractNetWork=Res
	end function


	subroutine Hij_phi(state,site1,site2,H)
		type(Tensor),intent(inout)::state
		type(Tensor),intent(in)::H
		integer,intent(in)::site1(2),site2(2)
		character*10::name1,name2
		type(Tensor)::allname
		allname=state%getName()
		name1='A'+site1(1)+'_'+site1(2)+'.n'
		name2='A'+site2(1)+'_'+site2(2)+'.n'
		state=contract(H,[3,4],state,[name1,name2])
		call state%setName(1,name1)
		call state%setName(2,name2)
		call state%permute(allname%ai())
		return
	end subroutine	

	subroutine  H_phi(H,state,L1,L2)
		type(Tensor),intent(in)::H
		type(Tensor),intent(inout)::state
		integer,intent(in)::L1,L2
		type(Tensor)::statei,sumState
		integer::i,j,k
		call sumState%empty()
		do i=1,L1
			do j=1,L2-1
				statei=state
				call Hij_phi(statei,[i,j],[i,j+1],H)
				if(sumState%getFlag())then
					sumState=sumState+statei
				else
					sumState=statei
				end if
			end do
			statei=state
			call Hij_phi(statei,[i,L2],[i,1],H)
			sumState=sumState+statei
		end do
		do i=1,L1-1
			do j=1,L2
				statei=state
				call Hij_phi(statei,[i,j],[i+1,j],H)
				sumState=sumState+statei
			end do
		end do
		do j=1,L2
			statei=state
			call Hij_phi(statei,[L1,j],[1,j],H)
			sumState=sumState+statei
		end do
		state=sumState
		return
	end subroutine	

	real*4 function Energy(A,H)
		type(Tensor),intent(in)::A(:,:)
		type(Tensor),intent(in)::H
		integer::L1,L2
		type(Tensor)::stateL,stateR
		L1=size(A,1)
		L2=size(A,2)
		stateR=contractNetWork(A) 
		stateL=stateR
		call H_phi(H,stateR,L1,L2)
		Energy=(stateL.x.stateR)/stateL%norm2()
		return
	end function

end module
program Example
	use SimpleUpdate
	integer::L1,L2,D,i,runningnum,runningTime
	real*4::tau,E
	type(Node),allocatable::A(:,:)
	type(Tensor),allocatable::B(:,:)
	type(Tensor)::H,expmH
	L1=3
	L2=3
	D=3
	tau=-0.05
	runningnum=100
	call set_seed(-111)
	call initialPEPS(A,L1,L2,D)
	call initialH(H)
	write(*,*)" "
	write(*,*)" "
	write(*,*)" "
	write(*,*)" "
	write(*,*)"**************************************************"
	write(*,*)"*  This is a example code for simple update      *"     
	write(*,*)"**************************************************"
	write(*,*)trim(L1+'*'+L2+' PBC Latticte system')
	write(*,*)'running Simple Update,the error of which will be near 10^-2'
	do runningTime=1,5
		write(*,*)'running evolution for tau=',tau
		expmH=gate(H,tau)
		do i=1,runningnum
			call sampleUpdate(A,expmH,L1,L2,D)
		end do
		tau=tau*0.7
	end do
	write(*,*)'Running the Energy,contract the whole network'
	allocate(B(L1,L2))
	call contractRightDownContract(A,B)
	E=Energy(B,H)/(L1*L2)
	write(*,*)'E=',E
	write(*,*)'The exact ground state enery is −0.441001111'
	write(*,*)'The error of Simple Update is',(E+0.441001111)/0.441001111
	write(*,*)" "
	write(*,*)" "

end