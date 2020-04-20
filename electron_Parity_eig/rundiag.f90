program aaaa
	use Tensor_type
	use SymDimension_typede
	use Dimension_typede
	use fermiTensor
	use eigen_value
	use initial_module
	use Tools
	use U1Tool
	use QuantumNumber_Type 
	implicit none
	type(eigenvalue)::eiger
	integer::L1,L2,num_ouput,i,ncv,j,TotalN,k,l
	integer::proId,proNum,ierr
	type(Tensor)::Res(2),parameter(2),tempState
	real*8::time1,time2
	logical::out_eigenvector
	character(len=200)::chartime
	L1=3
	L2=3
	num_ouput=3
	ncv=10
	out_eigenvector=.false.
	if((L1*L2).gt.16)then
		call writemess('The input lattice size is too large, the code will kill a lot of memory')
		call error_stop
	end if
	call eiger%set_ncv(ncv)
	call writemess(' ')
	call writemess(' ')
	call writemess('*********** Run the Diag for free electron model with Parity  ************')
	call set_fermi_symmetry('Parity')
	call writemess('L1='+L1+',L2='+L2)
	call writemess('ncv='+ncv)
	parameter(1)=[L1,L2]
	call eiger%set_ncv(ncv)
	call eiger%set_Hphi(H_phi)
	call eiger%set_print_num(10)
	call cpu_time(time1)
	Res=eiger%eig(parameter,'SR',num_ouput,out_eigenvector)
	call writemess('The some eigvalues are:')
	call writemess(dble(Res(1)/(L1*L2)))
	call cpu_time(time2)
	call system_time(time2-time1,chartime)
	call writemess('Running time'+(' '+chartime))

end 
