      SUBROUTINE disk(fourext,nphase,interp,P11,P12,coefs0)

*-------------------------------------------------------------------------------
*     Coefs contains only the first column of S
*-------------------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER i,j,interp,nmugs,nmugs0,nphase,nphash,nmat,nfou,
     .        lcrite,lcrito,mf,ii,jj,k,ibase,nphase1

      INTEGER iunpp
      PARAMETER (iunpp=43)

      DOUBLE PRECISION step,coefs0,alpha

      DOUBLE PRECISION pi
      PARAMETER (pi=3.1415926535898D0)

      DOUBLE PRECISION xmu(nmuMAX),wmu(nmuMAX),
     .                 xmu0(nmuMAX),wmu0(nmuMAX),
     .                 coefs(0:ncoefsMAX,nmatMAX),
     .                 costh(nphaseMAX),
     .                 P11(nphaseMAX),P12(nphaseMAX),
     .                 rfou(nmatMAX*nmuMAX,nmuMAX,0:nfouMAX),
     .                 rint(nmatMAX*nmuMAX,nmuMAX,0:nfouMAX)

      CHARACTER  fourext*120 !fourext*40,fourfil*25

*-------------------------------------------------------------------------------
*     Initialize the planetary scattering matrix elements:
*-------------------------------------------------------------------------------
      DO i=1,nphaseMAX
         P11(i)= 0.D0
         P12(i)= 0.D0
      ENDDO

*-------------------------------------------------------------------------------
*     Fill the phase angle array:
*-------------------------------------------------------------------------------
      nphash= (nphase+1)/2
      step= 180.D0/DBLE(nphase-1)
      DO i=1,nphash
         alpha= 90.D0 + step*DBLE(i-1)
         costh(i)= DCOS(pi-alpha*pi/180.D0)
      ENDDO

*-------------------------------------------------------------------------------
*     Open and read the Fourier coefficients file:
*
*     Get the Gaussian integration points and the supermatrix weights:
*     NOTE: the supermatrix weights are not identical to the Gaussian weights!
*-------------------------------------------------------------------------------
!      fourfil='four_'//trim(fourext)
!      print*,'opening:',trim(fourfil)//'.dat'
!      OPEN(unit=iunpp,file=trim(fourfil)//'.dat',
!     .     status='old',err=999)

!      fourfil='four_'//trim(fourext)
      OPEN(unit=iunpp,file=trim(fourext)//'.dat',
     .     status='old',err=999)

      READ(iunpp,*) nmat
      READ(iunpp,*) nmugs0


      DO i=1,nmugs0
         READ(iunpp,*) xmu0(i),wmu0(i)
      ENDDO
      nfou=0
10    DO i=1,nmugs0
         ibase= (i-1)*nmat
         DO j=1,nmugs0
            READ(iunpp,*,END=100) mf,ii,jj,
     .                            (rfou(ibase+k,j,nfou),k=1,nmat)
         ENDDO

      ENDDO
      nfou= nfou+1
      GOTO 10

100   CONTINUE

      CLOSE(iunpp)

*-------------------------------------------------------------------------------
*     In case of interpolation: determine the new Gaussian abscissas
*-------------------------------------------------------------------------------
      IF (interp.EQ.1) THEN
         nmugs= 2*nmugs0
         CALL gauleg(2*nmuMAX,nmugs,0.D0,1.D0,xmu,wmu)
         DO i=1,nmugs
            wmu(i)= DSQRT(2.D0*xmu(i)*wmu(i))
         ENDDO
      ELSE
         nmugs= nmugs0
         DO i=1,nmugs
            xmu(i)= xmu0(i)
            wmu(i)= wmu0(i)
         ENDDO
      ENDIF

*-------------------------------------------------------------------------------
*     In case of interpolation: determine the new Fourier coefficients
*-------------------------------------------------------------------------------
      CALL rinter(interp,rfou,nfou,nmugs0,xmu0,wmu0,nmugs,xmu,wmu,nmat,
     .            rint)

*-------------------------------------------------------------------------------
*     Calculate the expansion coefficients of the planetary scattering matrix:
*-------------------------------------------------------------------------------
      CALL calcoefs(nmugs,xmu,wmu,nmat,rint,nfou,coefs,lcrite,lcrito)

*-------------------------------------------------------------------------------
*     Loop over the planetary expansion coefficients:
*-------------------------------------------------------------------------------

      nphase1= nphase+1
      CALL diskint(interp,nphase1,costh,lcrite,lcrito,coefs,P11,P12)
      coefs0=coefs(0,1)

      GOTO 1000

*-------------------------------------------------------------------------------
*     Error messages:
*-------------------------------------------------------------------------------
999   WRITE(*,*) 'Error opening Fourier coefficients file!'
      STOP

*-------------------------------------------------------------------------------
*     End of the integration:
*-------------------------------------------------------------------------------
1000  continue !WRITE(*,*)


*-------------------------------------------------------------------------------
      RETURN
      END
