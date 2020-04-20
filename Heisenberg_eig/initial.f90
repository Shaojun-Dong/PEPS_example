module initial_module
	use Tensor_type
	use Tools
	implicit none
contains


	subroutine initialState(state,L1,L2)
		integer::L1,L2,TotalN
		type(Tensor),intent(inout)::state
		integer::i,j,k
		integer::dimen(L1*L2)
		dimen=2
		call state%allocate(dimen,'real*8')
		k=0
		do j=1,L2
			do i=1,L1
				k=k+1
				call state%setName(k,'A'+i+'_'+j+'.n')
			end do
		end do
		call state%Random([-1d0,1d0])
		return
	end subroutine
	
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

	subroutine Hij_phi(state,site1,site2,H)
		type(Tensor),intent(inout)::state
		type(Tensor),intent(in)::H
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
		type(Tensor),intent(inout)::state
		type(Tensor),intent(in)::H
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
		type(Tensor)::sumstate,phi1
		type(Tensor),save::H,state
		integer::i,j,k,L1,L2
		
		L1=parameter(1)%ii(1)
		L2=parameter(1)%ii(2)
		if(.not.inoutstate%getFlag())then
			call initialState(state,L1,L2)
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
