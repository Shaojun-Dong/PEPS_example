module initial_module
	use Tensor_type
	use FermiTensor
	use SymDimension_typede
	use QuantumNumber_Type 
	use U1Tool
	use Tools
	implicit none
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

	subroutine initialState(state,L1,L2,TotalN)
		integer::L1,L2,TotalN
		type(fTensor),intent(inout)::state
		integer::i,j,k
		type(quanNum)::QN(L1*L2+1)
		do i=1,L1*L2
			QN(i)=U1quantumnumber([0.,1.],[1,1],1,1)
		end do
		QN(L1*L2+1)=U1quantumnumber([real(TotalN)],[1],-1,-1)
		call state%allocate(QN,'real*8')
		k=0
		do j=1,L2
			do i=1,L1
				k=k+1
				call state%setName(k,'A'+i+'_'+j+'.n')
			end do
		end do
		call state%SymRandom([-1d0,1d0])
		return
	end subroutine
	
	subroutine initialH(H)
		type(fTensor),intent(inout)::H
		type(fTensor)::C1,C2,temp
		C1=C_operator('+')
		C2=C_operator('')
		call C1%setName('H1')
		call C2%setName('H2')
		H=C1.kron.C2
		C1=C_operator('')
		C2=C_operator('+')
		call C1%setName('H1')
		call C2%setName('H2')
		temp=C1.kron.C2
		H=H-temp
		call H%permute(['H1.n ','H2.n ','H1.n2','H2.n2'])
		call H%emptyZeroBlock()
		call H%SymCheck()
		return
	end subroutine

	subroutine Hij_phi(state,site1,site2,H)
		type(fTensor),intent(inout)::state
		type(fTensor),intent(in)::H
		integer,intent(in)::site1(2),site2(2)
		type(Tensor)::allname
		character*10::name1,name2
		name1='A'+site1(1)+'_'+site1(2)+'.n'
		name2='A'+site2(1)+'_'+site2(2)+'.n'
		allname=state%getName()
		state=contract(H,['H1.n2','H2.n2'],state,[name1,name2])
		call state%setName('H1.n',name1)
		call state%setName('H2.n',name2)
		call state%permute(allname%ai())
		return
	end subroutine
	
	subroutine Hi_phi(state,site1,H)
		type(fTensor),intent(inout)::state
		type(fTensor),intent(in)::H
		integer,intent(in)::site1(2)
		type(Tensor)::allname
		character*10::name1
		name1='A'+site1(1)+'_'+site1(2)+'.n'
		allname=state%getName()
		state=contract(H,'H.n2',state,name1)
		call state%setName('H.n',name1)
		call state%permute(allname%ai())
		return
	end subroutine

	subroutine H_phi(parameter,inoutstate)
		type(Tensor)::parameter(:)
		type(Tensor)::inoutstate
		type(fTensor)::sumstate,phi1
		type(fTensor),save::H,state
		integer::i,j,k,L1,L2
		integer::TotalN
		
		L1=parameter(1)%ii(1)
		L2=parameter(1)%ii(2)
		TotalN=parameter(1)%ii(3)
		if(.not.inoutstate%getFlag())then
			call initialState(state,L1,L2,TotalN)
			inoutstate=state
			call initialH(H)
			return
		end if
		state=inoutstate
		call sumstate%empty()
		do i=1,L1
			do j=1,L2-1
				phi1=state
				call Hij_phi(phi1,[i,j],[i,j+1],H)
				if(sumstate%getFlag())then
					sumstate=sumstate+phi1
				else
					sumstate=phi1
				end if
			end do
		end do
		do i=1,L1-1
			do j=1,L2
				phi1=state
				call Hij_phi(phi1,[i,j],[i+1,j],H)
				sumstate=sumstate+phi1
			end do
		end do
		inoutstate=sumstate
		return
	end subroutine
	
end module
