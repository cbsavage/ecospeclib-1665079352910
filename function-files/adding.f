      SUBROUTINE adding(fourext,a,b,coefs,ncoefs,nlays,nmug,
     .                  nmat,asurf,ifa)

**********************************************************************
*                A D D I N G    M E T H O D
*          F O R    P O L A R I Z E D   L I G H T
**********************************************************************
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,j,m,nlays,nmug,nmat,ifa

      INTEGER ncoefs(nlevsMAX),iadd(nlevsMAX,0:nfouMAX),
     .        M0(nlevsMAX),M1(nlevsMAX),M2(nlevsMAX)

      INTEGER iunpp
      PARAMETER (iunpp=43)

*---------------------------------------------------------------------
*     ifa=0 : no Fresnel interface
*     ifa=1 : a Fresnel interface
*     xm : index of refraction of water with respect to air
*---------------------------------------------------------------------
      DOUBLE PRECISION asurf

      DOUBLE PRECISION xm
      PARAMETER (xm=1.34D0)

*---------------------------------------------------------------------
      DOUBLE PRECISION a(nlevsMAX),b(nlevsMAX),
     .                 xmu(nmuMAX),smf(nmuMAX),
     .                 ebbot(nmuMAX),ebtop(nmuMAX),eblam(nmuMAX),
     .                 coefs(nmatMAX,nmatMAX,0:ncoefsMAX,nlevsMAX)

      DOUBLE PRECISION Rmbot(nsupMAX,nsupMAX),Rmtop(nsupMAX,nsupMAX),
     .                 Tmbot(nsupMAX,nsupMAX),Tmtop(nsupMAX,nsupMAX),
     .                 Rmsbot(nsupMAX,nsupMAX)

      DOUBLE PRECISION Rmlam(nsupMAX,nsupMAX),Tmlam(nsupMAX,nsupMAX),
     .                 Rmslam(nsupMAX,nsupMAX),
     .                 Rmifa(nsupMAX,nsupMAX),Tmifa(nsupMAX,nsupMAX),
     .                 Rmsifa(nsupMAX,nsupMAX)

      DOUBLE PRECISION Zmplus(nsupMAX,nsupMAX),Zmmin(nsupMAX,nsupMAX)

      DOUBLE PRECISION xmut(nmuMAX),
     .                 Rf(nmuMAX,2),Tf(nmuMAX,3),
     .                 Rfs(nmuMAX,3),Tfs(nmuMAX,3)

!      CHARACTER fourext*40,fourfil*25
      CHARACTER fourext*120

      LOGICAL nextm,verbo
      verbo= .false.

*---------------------------------------------------------------------
*     Open the Fourier-output file:
*---------------------------------------------------------------------
!      fourfil='four_'//trim(fourext)

!      OPEN(unit=iunpp,file=trim(fourfil)//'.dat')
      OPEN(unit=iunpp,file=trim(fourext)//'.dat')
      WRITE(iunpp,'(I3)') nmat

*---------------------------------------------------------------------
*     Initialize the mu-values and the supermatrix factors:
*---------------------------------------------------------------------
      CALL setmu(nmug,iunpp,xmu,smf)

*---------------------------------------------------------------------
*     Calculate the bounds M0, M1, and M2 on the Fourier-index:
*---------------------------------------------------------------------
      CALL setfou(coefs,ncoefs,nlays,a,b,xmu,nmug,M0,M1,M2,iadd)

*---------------------------------------------------------------------
*     Calculate the reflection matrices of the interface:
*---------------------------------------------------------------------
      IF (ifa.EQ.1) CALL fmatri(xmu,nmug,xm,Rf,Tf,Tfs,Rfs,xmut)

*---------------------------------------------------------------------
*     Loop over the Fourier terms:
*---------------------------------------------------------------------
      m= -1
1000  m= m+1

*---------------------------------------------------------------------
*     Initialization:
*---------------------------------------------------------------------
      CALL inig(Rmbot,Rmtop,Tmbot,Tmtop,Rmlam,Tmlam,Rmsbot,Rmslam,
     .          Rmifa,Tmifa,Rmsifa)

*---------------------------------------------------------------------
*     Fill the arrays with reflection by a Lambertian surface:
*---------------------------------------------------------------------
      IF (m.EQ.0) THEN
         iadd(1,0)= 1

*        Fill Rmtop and Tmtop with the surface reflection:
         CALL layer0(asurf,smf,nmug,nmat,ebtop,Rmtop,Tmtop)

*        Shift everything to the surface arrays:
         CALL top2bot(nmat,nmug,ebtop,ebbot,Rmtop,Tmtop,
     .                Rmlam,Tmlam,Rmslam)
      ENDIF

!         write(*,*) Rmtop(1,1),Rmsbot(1,1)
*       write(*,*) 'bf:',Rmsbot(1,1),Rmbot(1,1),Rmlam(1,1),Rmtop(1,1)
*---------------------------------------------------------------------
*     Loop over the model atmosphere:
*---------------------------------------------------------------------
      DO i=1,nlays

*---------------------------------------------------------------------
*        Calculate the m-th Fourier coefficient of the scat. matrix:
*---------------------------------------------------------------------
         CALL setzm(m,i,coefs,ncoefs,xmu,nmug,nmat,Zmmin,Zmplus)
         IF (m.EQ.0) CALL renorm(Zmmin,Zmplus,nmug,nmat,xmu,smf)

*---------------------------------------------------------------------
*        Calculate the m-th Fourier coefficient of the layer:
*---------------------------------------------------------------------
         CALL layerm(m,M0,M1,M2,i,xmu,smf,nmug,nmat,coefs,ncoefs,
     .               Zmplus,Zmmin,a(i),b(i),ebtop,Rmtop,Tmtop)

!         write(*,*) Rmtop(1,1),Rmsbot(1,1)
*---------------------------------------------------------------------
*        Add the top layer to the bottom layer:
*---------------------------------------------------------------------
         IF (i.EQ.1) THEN
            CALL top2bot(nmat,nmug,ebtop,ebbot,Rmtop,Tmtop,
     .                      Rmbot,Tmbot,Rmsbot)
         ELSE
            CALL addlay(nmat,nmug,ebtop,ebbot,iadd(i,m),Rmtop,Tmtop,
     .                  Rmbot,Tmbot,Rmsbot)
         ENDIF

*         write(*,*) Rmtop(1,1),Rmsbot(1,1)

!        stop
      ENDDO
*         write(*,*) 'atm:',Rmbot(4:6,4:6)

*         stop
*---------------------------------------------------------------------
*     Add the Fresnel surface to the atmosphere:
*---------------------------------------------------------------------
      IF (ifa.EQ.1) THEN
         CALL addi(m,nmug,nmat,Rf,xmu,ebbot,Rmbot,Tmbot,Rmsbot)
      ENDIF

*       write(*,*) Rmbot(1,1),Rmlam(1,1)

*---------------------------------------------------------------------
*     Add the Lambertian reflecting surface to the atmosphere:
*     This only affects the m=0 term!
*---------------------------------------------------------------------
      IF (m.EQ.0) THEN
         CALL addlay(nmat,nmug,ebbot,eblam,iadd(1,0),Rmbot,Tmbot,
     .               Rmlam,Tmlam,Rmslam)
         CALL assign(Rmbot,Rmlam,nmat,nmug)
         CALL assign(Tmbot,Tmlam,nmat,nmug)
      ENDIF

*       write(*,*) 'aft:',Rmsbot(1,1),Rmbot(1,1),Rmlam(1,1),Rmtop(1,1)
*       write(*,*) Rmbot(1,1),Rmlam(1,1)
*       stop
 !     write(*,*) 'go to newfou'
*---------------------------------------------------------------------
*     Write the m-th Fourier coefficient to file:
*---------------------------------------------------------------------
      CALL newfou(m,Rmbot,smf,iunpp,nmat,nmug)

*---------------------------------------------------------------------
*     Check the convergence of the Fourier series:
*---------------------------------------------------------------------
      CALL endfou(m,M1,M0,nlays,xmu,nmug,nmat,nextm)

*---------------------------------------------------------------------
*     Next term of the Fourier-loop:
*---------------------------------------------------------------------
      IF (m.LT.2) nextm=.true.
      IF (nextm) GOTO 1000

*----------------------------------------------------------------------
*     Close the Fourier-output file:
*----------------------------------------------------------------------
      CLOSE(iunpp)
************************************************************************
      RETURN
      END
