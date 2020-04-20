program example
	use SimpleUpdate
	use DoubleLayerEnergy_type
	use DoubleLeyerAuxTensor_type
	use Tensor_type
	use SymTensor_type
	use SymDimension_typede
	use QuantumNumber_Type 
	use ParityTool
	use Tools
	implicit none
	integer::L1,L2,D,i,runningnum,runningTime,initialDeg
	integer::Dc
	real*8::tau,E,Alltau(8)
	type(Node),allocatable::A(:,:)
	type(SymTensor),allocatable::B(:,:)
	type(SymTensor)::H,expmH
	write(*,*)" "
	write(*,*)" "
	write(*,*)" "
	write(*,*)" "
	write(*,*)"**************************************************"
	write(*,*)"*  This is a example code for simple update      *"    
	write(*,*)"**************************************************"
	call set_symmetry('U(1)')
	L1=3
	L2=4
	D=4
	initialDeg=2
	Dc=16
	Alltau=[0.1,0.05,0.01,0.008,0.005,0.001,0.0005,0.0001]
	Alltau=-1d0*Alltau
	runningnum=300
	call initialPEPS(A,L1,L2,initialDeg)
	call initialH(H)
	write(*,*)trim(L1+'*'+L2+' OBC Latticte system')
	write(*,*)trim('D='+D)
	write(*,*)'running Simple Update,the error of which will be near 10^-2'
	call reset_time_calculator(runningnum*size(Alltau),30) 
	do runningTime=1,size(Alltau)
		tau=Alltau(runningTime)
		expmH=gate(H,tau)
		do i=1,runningnum
			call sampleUpdate(A,expmH,D)
			call time_calculator()
		end do
	end do
	write(*,*)'Running the Energy,contract the whole network'
	E=ContractEnergy(A,H,Dc)
	write(*,*)"Energy=",E
contains

end
