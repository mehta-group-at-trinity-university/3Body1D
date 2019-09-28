program SVD
implicit none
        real*8 :: n,w,m,Rmax,res,temp,overlap,tempA,tempB,epsout,&
                Shift,alpha,r0diatom,etaOVERpi,xMin,xMax,yMin,yMax,RDerivDelt,RFirst,RLast,dscale,&
                e,eMin,eMax,m1,m2,m3,mu3,DD,L
        real*8, allocatable :: nodesR(:),weightsR(:),holder(:),Hsparse(:),&
                Bsparse(:),Pmu(:),Pnu(:),eigenVals(:),ress(:),&
                alphaR(:),alphaI(:),work(:),EVecs(:,:),VL(:,:),beta(:)
        real*8, allocatable :: dXr(:,:),TmatR(:,:),weightsDMR(:,:),S(:,:),eigenVecs(:,:),&
                iPsi(:,:),jPsi(:,:),H(:,:),B(:,:)
        real*8, allocatable :: Uad(:,:,:),Psi(:,:,:),O(:,:)
        integer :: i,j,k,x,sR,sRc,row,mu,nu,inu,jmu,rCountA,numN0A,rCountB,numN0B,probSize,&
                xDim,yDim,ncv,matrixdim,HalfBandWidth,set,&
                NumStates,PsiFlag,CouplingFlag,LegPoints,order,Left,Right,Bottom,Top,xNumPoints,yNumPoints,&
                fparam(64),info,loop,m0,mm,&
                Tstart,Tend,rate
        integer, allocatable :: Hrow(:),Hcol(:),Brow(:),Bcol(:),indexOf(:,:)
        character*64 LegendreFile
        character* 64 filename

        !open(unit=1,file="nodesAndWeights10.dat")
        !open(unit=2,file="nodesAndWeights20.dat")
        !open(unit=3,file="nodesAndWeights30.dat")
        !open(unit=4,file="nodesAndWeights40.dat")
        !open(unit=5,file="nodesAndWeights50.dat")
        !open(unit=6,file="nwTEMP.dat")
        !open(unit=7,file="nodesAndWeights70.dat")
        !open(unit=8,file="nodesAndWeights80.dat")
        !open(unit=9,file="nodesAndWeights90.dat")
        !open(unit=10,file="nodesAndWeights100.dat")
        !open(unit=15,file="nodesAndWeights150.dat")
        !open(unit=20,file="nodesAndWeights200.dat")

        !open(unit=101,file="eigenVals-Vecs-100x60x60x200.dat")
        !open(unit=200,file="adiabaticPhi.dat")
        !open(unit=201,file="adiabaticPhi-200x60x60x100x100-m12.dat")
        !open(unit=300,file="eigenVals-200x60x60x100x100-m12.dat")
        !open(unit=301,file="eigenVecs-200x60x60x100x100-m12.dat")
        !open(unit=400,file="overlap100x400.dat")
        !open(unit=401,file="overlap200x60x60x100x100-m12.dat")

        do set=1,1
                !open(unit=15,file="nodesAndWeights150.dat")
                open(unit=10,file="nodesAndWeights100.dat")
                !open(unit=6,file="nodesAndWeights60.dat")
                if (set==1) filename = "eigenVals-100x70x70x200x50-m12Even.dat"
                if (set==2) filename = "eigenVals-100x80x80x200x50-m12Even.dat"
                if (set==3) filename = "eigenVals-100x90x90x200x50-m12Even.dat"
                if (set==4) filename = "eigenVals-150x80x80x200x100-m12Even.dat"
                if (set==5) filename = "eigenVals-200x80x80x200x250-m12.dat"
                if (set==6) filename = "eigenVals-200x80x80x200x300-m12.dat"
                open(unit=1,file=filename)
        
        call system_clock(Tstart)

        !Create the Gauss Lobatto DVR functions X for the R coordinate
        !if(set==1) x=15
        !if(set==2) x=10
        !if(set==3) x=6
        x=10
        read(x,*) sR
        Rmax=10d0
        sRc=sR-2
        
        !if(set==1 .OR. set==3 ) then
        !        m1=1d0
        !        m2=m1
        !        m3=m1
        !        m4=m1
        !else 
                m1=1d0
                m2=m1
                m3=m1
        !end if
        mu3=(m1*m2*m3/(m1+m2+m3))**(1.0d0/3.0d0);
        
        allocate(nodesR(sR),weightsR(sR),dXr(sRc,sR),holder(sR),weightsDMR(sR,sR),TmatR(src,src))

        !!!!!!
        !Create the Gauss Lobatto DVR functions X for the R coordinate and the Kmat matrix
        do i=1,sR
                read(x,*) n,w
                nodesR(i)=n
                weightsR(i)=w
        end do
        call rescale(sR,0d0,Rmax,nodesR)
        do i=1,sRc
                do j=1,sR
                        call calcDX(holder,nodesR,nodesR(j),i,sR,weightsR,res)
                        dXr(i,j)=res
                end do
        end do
        weightsDMR=0d0
        do i=1,sR
                do j=1,sR
                        if(i==j) then
                                weightsDMR(i,j)=weightsR(i)*nodesR(i)**2d0
                        end if
                end do
        end do
        TmatR=0d0
        TmatR = (0.5d0/mu3)*MATMUL(dXr,MATMUL(weightsDMR,TRANSPOSE(dXr)))
        !!!!!!
        
        !!!!!!
        !obtain the adiabatic eigenvalues and eigenvectors from Kelly's code
        NumStates = 10
        LegendreFile='Legendre.dat'
        LegPoints=10
        if(set==1) Shift=-7d0
        if(set==2) Shift=-5*50
        if(set==3) Shift=-5*50
        if(set==4) Shift=-5*100
        if(set==5) Shift=-5*250
        if(set==6) Shift=-5*300
        shift2 = -0.01d0
        order=5
        Left=1
        Right=0
        alpha=1d0
        if(set==1)xNumPoints = 40
        if(set==2)xNumPoints = 80
        if(set==2)xNumPoints = 90
        !xNumPoints=80
        xMin=0d0
        xMax=1.570796326794897
        if(set==1)yNumPoints = 70
        if(set==2)yNumPoints = 80
        if(set==3)yNumPoints = 90
        !yNumPoints=80
        yMin=0d0
        yMax =1.570796326794897
        RDerivDelt=0.0001d0
        DD = 0.853138729379683d0
        L = 1.0d0
        RFirst=nodesR(2)
        RLast=nodesR(sR-1)
        dscale=10000d0
        
        !ncv=2*NumStates
        xDim=xNumPoints+Order-3
        yDim=yNumPoints+Order-3+1
        MatrixDim=xDim*yDim
        HalfBandWidth=yDim*Order+Order

        probSize=sRc*numStates

        print *, NUmStates
        print *, HalfBandWidth, MatrixDim
        allocate(Uad(sRc,numStates,2),Psi(sRc,MatrixDim,numStates),S(HalfBandWidth+1,MatrixDim),&
                indexOf(sRc,NumStates),Pnu(numStates),Pmu(numStates),O(probSize,probSize),&
                iPsi(MatrixDim,numStates),jPsi(MatrixDim,numStates))
        call HHL1DHyperspherical(NumStates,0,0,LegendreFile,LegPoints,Shift,shift2,Order,Left,Right,alpha,m1,&
                m2,m3,xNumPoints,xMin,xMax,sRc,RDerivDelt,RFirst,RLast,DD,L,&
                nodesR(2:sR-1),Uad,Psi,numstates,MatrixDim,S,HalfBandWidth+1,set)
        !!!!!
              
        !S=0d0
        !do i=1,HalfBandWidth+1
        !        read(201,*) (S(i,j),j=1,MatrixDim)
        !end do
        !do i=1,sRc
        !        read(201,*) (Uad(i,j,1),j=1,numstates)
        !        do j=1,MatrixDim
        !                read(201,*) (Psi(i,j,k),k=1,numstates)
        !        end do
        !end do
        
        print *,"here"
        
        eMin=0d0
        do i=1,sRc
                e=Uad(i,1,1)
                if (e<eMin) then
                        eMin=e
                end if
        end do
        eMax=Uad(sRc,1,1)

        !!!!!!
        !creat the complete H matrix in sparse formation
        row=1
        do i=1,sRc
                do j=1,numStates
                        indexOf(i,j)=row
                        row=row+1
                end do
        end do

        !!!!!!!!
        O=0d0
        do i=1,sRc
                print *,i
                do j=1,sRc
                        do k=1,MatrixDim
                                do mu=1,numStates
                                        iPsi(k,mu)=Psi(i,k,mu)
                                        jPsi(k,mu)=Psi(j,k,mu)
                                end do
                        end do
                        call CalcOMatrix(Numstates,HalfBandWidth,MatrixDim,iPsi,jPsi,&
                                        S,O,sRc,i,j,indexOf,probSize)
                end do
        end do
        !do i=1,sRc
        !        print *, i
        !        do j=1,sRc
        !                do nu=1,numStates
        !                        do mu=1,numStates
        !                                inu=indexOf(i,nu)
        !                                jmu=indexOf(j,mu)
        !                                read(401,*)O(inu,jmu)
        !                        end do
        !                end do
        !        end do
        !end do
        !!!!!!!!        

        print *, 'Got S'

        numN0A=0
        rCountA=1
        numN0B=0
        rCountB=1
        allocate(H(probsize,probsize),B(probsize,probsize),alphaR(probsize),alphaI(probsize),beta(probsize),VL(probsize,probsize),EVecs(probsize,probsize),work(8*probsize))
        H=0d0
        B=0d0
        do x=1,2
                print *,x
                do i=1,sRc
                        do nu=1,numStates
                                do j=1,sRc
                                        do mu=1,numStates
                                                inu=indexOf(i,nu)
                                                jmu=indexOf(j,mu)
                                                tempA = 0d0
                                                tempB = 0d0
                                                !Pnu=Psi(i,nu,1:numStates)
                                                !Pmu=Psi(j,mu,1:numStates)
                                                Overlap=O(inu,jmu)
                                                if((overlap /= 0d0).and.(TmatR(i,j) /= 0d0)) then
                                                        tempA=tempA+TmatR(i,j)*overlap
                                                end if
                                                if((i==j).and.(nu==mu).and.(Uad(i,nu,1) /= 0d0)) then
                                                        tempA=tempA+Uad(i,nu,1)
                                                end if
                                                if((i==j).and.(overlap /= 0d0)) then
                                                        tempB=tempB+overlap*nodesR(i+1)**2
                                                end if
                                                if(tempA /= 0) then
                                                        !numN0A=numN0A+1
                                                        !rCountA=rCountA+1
                                                        !if(x==2) then
                                                        !        Hsparse(numN0A)=tempA
                                                        !        Hcol(numN0A)=jmu
                                                        !        Hrow(inu+1)=Hrow(inu+1)+1
                                                        !end if
                                                end if
                                                if(tempB /= 0) then
                                                        !numN0B=numN0B+1
                                                        !rCountB=rCountB+1
                                                        !if(x==2) then
                                                        !        Bsparse(numN0B)=tempB
                                                        !        Bcol(numN0B)=jmu
                                                        !        Brow(inu+1)=Brow(inu+1)+1
                                                        !end if
                                                end if
                                                if (x==2) then
                                                        H(inu,jmu)=H(inu,jmu)+TmatR(i,j)*overlap
                                                        if((i==j).and.(nu==mu))H(inu,jmu)=H(inu,jmu)+Uad(i,nu,1)
                                                end if
                                        end do
                                end do
                        end do
                end do
                if(x==1) then
                        !allocate(Hsparse(1:numN0A),Hcol(numN0A),Hrow(probSize+1))
                        !Hrow=0
                        !Hrow(1)=1
                        !numN0A=0
                        !rCountA=1
                        !allocate(Bsparse(1:numN0B),Bcol(numN0B),Brow(probSize+1))
                        !Brow=0
                        !Brow(1)=1
                        !numN0B=0
                        !rCountB=1
                end if
        end do
        !do i=2,probSize+1
        !        Hrow(i)=Hrow(i)+Hrow(i-1)
        !        Brow(i)=Brow(i)+Brow(i-1)
        !end do
        !write(*,*) "If",size(Hsparse),"=",Hrow(probSize+1)-1,"=",size(Hcol),"Then Good"
        !write(*,*) "If",size(Bsparse),"=",Brow(probSize+1)-1,"=",size(Bcol),"Then Good"
        !do i=1,numN0B
        !        if(Bsparse(i)<=0d0) then
                         !print *, "-B"
        !        end if
        !end do


        !!!!!!
        !Hamiltonian symatry check
        k=0
        i=6
        !do i=1,sRc
        !        do nu=1,numStates
        !                        do mu=1,numStates
        !                                inu=indexOf(i,nu)
        !                                jmu=indexOf(i,mu)
        !                                !if(O(inu,jmu) /= ) k=k+1
        !                                print *, O(inu,jmu)
        !                        end do
        !        end do
        !end do
        print *, "DIFFS"
        !do i=1,sRc
        !        do nu=1,numStates
        !                !do j=1,sRc
        !                        !do mu=1,numStates
        !                                inu=indexOf(i,nu)
        !                                jmu=indexOf(i+1,nu)
        !                                print *, O(inu,jmu)-O(jmu,inu)
        !                                !if (O(inu,jmu) .ne. O(jmu,inu)) k=k+1
        !                        !end do
        !                !end do
        !        end do
        !end do
        !print *, k
        !!!!!

        m0=100
        !allocate(EigenVals(m0),EigenVecs(probSize,m0),ress(m0*2)) 
        !call feastinit(fparam)
        !fparam(1)=1
        !call dfeast_scsrgv('F',probSize,Hsparse,Hrow,Hcol,Bsparse,Brow,Bcol,&
        !        fparam,epsout,loop,eMin,eMax,m0,EigenVals,EigenVecs,mm,ress,info)
        
        call dgeev('N','V',probsize,H,probsize,alphaR,alphaI,VL,probsize,EVecs,probsize,work,8*probsize,INFO)
        call system_clock(Tend,rate)
        print *,info
        print *,(Tend-Tstart)/rate
        !write(101,*) (Tend-Tstart)/rate, mm, eMin, eMax
        write(1,*) (Tend-Tstart)/rate, numstates
        
        !write(101,*) (EigenVals(j),j=1,mm)
        !do i=1,probSize
        !        write(101,*) (EigenVecs(i,j),j=1,mm)
        !end do
        do i=1,probsize
                write(1,*) alphaR(i)
                !write(301,*) (EVecs(i,j),j=1,probsize)
        end do
        
        close(1)
        !close(2)
        !close(3)
        !close(4)
        !close(5)
        !close(6)
        !close(8)
        !close(10)
        close(20)
        !close(15)
        !close(100)
        !close(101)
        !close(102)
        !close(200)
        !close(201)
        !close(300)
        !close(301)
        !close(400)
        !close(401)
        
                
        deallocate(nodesR,weightsR,dXr,holder,weightsDMR,TmatR,Uad,Psi,S,&
        H,B,alphaR,alphaI,beta,VL,EVecs,work,indexOf,Pnu,Pmu,O,iPsi,jPsi)

        end do
        
end program SVD

subroutine rescale(n,a,b,nodes)
implicit none
        integer :: n,i
        real*8 :: nodes(*)
        real*8 :: x1,xN,a,b
        x1= nodes(1)
        xN= nodes(n)
        do i=1,n
                if((xN-x1)*(b-a)==0) then
                        nodes(i)=0d0
                else 
                        nodes(i)=a+(nodes(i)-x1)/(xN-x1)*(b-a)
                end if
        end do
end subroutine rescale

!find derivitave of the nodes
subroutine CalcDX(holder,nodes,node,i,sz,weights,res)
implicit none
        integer :: i,sz,j,k
        real*8 nodes(sz),holder(sz),weights(sz)
        real*8 val,node,res

        res=0
        do k=1,sz
                if(k /= i) then
                        do j=1,sz
                                if((nodes(i)-nodes(j))==0d0) then
                                        holder(j)=0d0
                                else 
                                        holder(j)=(node-nodes(j))/(nodes(i)-nodes(j))
                                end if
                        end do
                        holder(i)=1d0
                        holder(k)=1d0
                        if((nodes(i)-nodes(k))==0) then
                                val=0d0
                        else 
                                val=PRODUCT(holder)/((nodes(i)-nodes(k)))
                        end if
                        res = res + val
                end if
        end do
        res=res*(weights(i)**(-.5d0))
end subroutine CalcDX

subroutine CalcOMatrix(NumStates,HalfBandWidth,MatrixDim,iPsi,jPsi,S,O,sRc,i,j,indexOf,probSize)
implicit none
        integer NumStates,HalfBandWidth,MatrixDIm,i,j,k,R,sRc,nu,mu,&
                indexOf(sRc,numStates),inu,jmu,probSize
        real*8 ddot
        double precision iPsi(MatrixDim,NumStates),jPsi(MatrixDim,numStates),S(HalfBandWidth+1,MatrixDim),&
                O(probSize,probSize),TempPsi(MatrixDim)

        do nu=1,numStates
                call dsbmv('U',MatrixDim,Halfbandwidth,1.0d0,S,Halfbandwidth+1,iPsi(:,nu),1,0.0d0,tempPsi,1)
                do mu=1,numStates
                        inu=indexOf(i,nu)
                        jmu=indexOf(j,mu)
                        O(inu,jmu) = ddot(MatrixDim,TempPsi,1,jPsi(1,mu),1)
         !               write(401,*) O(inu,jmu)
                end do
        end do
end subroutine CalcOMatrix
