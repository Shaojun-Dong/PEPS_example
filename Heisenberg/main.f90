program example
	use SimpleUpdate
	use DoubleLayerEnergy_type
	use DoubleLeyerAuxTensor_type
	use Tensor_type
	implicit none
	integer::L1,L2,D,i,runningnum,runningTime,initialD
	integer::Dc
	real*8::tau,E,Alltau(8)
	type(Node),allocatable::A(:,:)
	type(Tensor),allocatable::B(:,:)
	type(Tensor)::H,expmH
	write(*,*)" "
	write(*,*)" "
	write(*,*)" "
	write(*,*)" "
	write(*,*)"**************************************************"
	write(*,*)"*  This is a example code for simple update      *"    
	write(*,*)"**************************************************"
	L1=3
	L2=4
	D=4
	initialD=4
	Dc=16
	Alltau=[0.1,0.05,0.01,0.008,0.005,0.001,0.0005,0.0001]
	Alltau=-1d0*Alltau
	runningnum=300
	call initialPEPS(A,L1,L2,initialD)
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
