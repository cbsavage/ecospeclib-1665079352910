      SUBROUTINE rdinput(nmat,asurf,naer,coefile,filenames,fourext,
     .                 nmug,nlevs,nlev0,wavfile,res,iabs,absfile,wav0,
     .                 baer,iaer,pres,temp,xx,nphase,
     .                 interp,ias,ifa,asurfile,nplanet)

************************************************************************
* DATE: July 2004
*
* xO3 and xH2O give the fraction of O3 and H2O molecules in an
* atmospheric layer
*
* AUTHOR: D. M. Stam
************************************************************************
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INTEGER iunin,i,k,j,iaer(nlevsMAX),interp,iabs,ias,ifa,
     .        nmat,nmug,nphase,naer,nlevs,nlev0,nplanet,res

      DOUBLE PRECISION wav0,asurf,
     .                 baer(nlevsMAX),hgpar(nlevsMAX,4),
     .                 pres(nlevsMAX),temp(nlevsMAX),
     .                 xx(7,nlevsMAX)

      CHARACTER*20 absfile,wavfile,coefile(ncoesMAX),
     .             asurfile
      CHARACTER filenames*120
      CHARACTER fourext*120,fourext1*3,fourext2*2,fourext3*3

*-----------------------------------------------------------------------
*     Open the input file and read the first lines:
*-----------------------------------------------------------------------
      iunin= 1
      OPEN(iunin,file='dap.in')
      DO i=1,4
         READ(iunin,*,err=999,end=998)
      ENDDO

*-----------------------------------------------------------------------
*     Read the GENERAL parameters:
*-----------------------------------------------------------------------
      nplanet=5

      READ(iunin,*,err=999,end=998) filenames
      READ(iunin,*,err=999,end=998) nmat
      READ(iunin,*,err=999,end=998) wavfile
      READ(iunin,*,err=999,end=998) res
      READ(iunin,*,err=999,end=998) iabs
      READ(iunin,*,err=999,end=998) absfile


      IF (nplanet.LT.2 .AND. nplanet.GT.4) THEN
         WRITE(*,*) 'rdinput: nplanet must be between 2 and 4!'
         STOP
      ENDIF

      IF (nmat.NE.1 .AND. nmat.NE.3 .AND. nmat.NE.4) THEN
         WRITE(*,*) 'rdinput: nmat is wrong (must be 1, 3, or 4)'
         STOP
      ENDIF

      DO i=1,4
         READ(iunin,*,err=999,end=998)
      ENDDO

*-----------------------------------------------------------------------
*     Read the SURFACE albedo:
*-----------------------------------------------------------------------
      READ(iunin,*,err=999,end=998) ias
      READ (iunin,*,err=999,end=998) asurf
      READ(iunin,*,err=999,end=998) asurfile
      READ(iunin,*,err=999,end=998) ifa

      IF (ias.NE.0 .AND. ias.NE.1) THEN
         WRITE(*,*) 'rdinput: wrong value for ias'
         STOP
      ENDIF
      IF (asurf.LT.0.D0 .OR. asurf.GT.1.D0) THEN
         WRITE(*,*) 'rdinput: wrong value for asurf'
         STOP
      ENDIF
      IF (ifa.NE.0 .AND. ifa.NE.1) THEN
         WRITE(*,*) 'rdinput: wrong value for ifa'
         STOP
      ENDIF

      DO i=1,4
         READ(iunin,*,err=999,end=998)
      ENDDO

*-----------------------------------------------------------------------
*     Read the GEOMETRIC parameters:
*-----------------------------------------------------------------------
      READ (iunin,*,err=999,end=998) nmug
      READ (iunin,*,err=999,end=998) interp
      READ (iunin,*,err=999,end=998) nphase

      IF (nmug.LT.0 .OR. (2*(nmug/2).NE.nmug)
     .              .OR. nmug.GT.nmuMAX) THEN
         WRITE(*,*) 'rdinput: wrong value for nmug'
         STOP
      ENDIF

      IF (nphase.GT.nphaseMAX) THEN
         WRITE(*,*) 'rdinput: nphase too large'
         STOP
      ENDIF

!       print*,2*nmug, nmuMAX

      IF (interp.EQ.1 .AND. (2*nmug).GT.nmuMAX) THEN
         WRITE(*,*) 'rdinput: too large value for nmug'
         STOP
      ENDIF

      DO i=1,4
         READ(iunin,*,err=999,end=998)
      ENDDO

*-----------------------------------------------------------------------
*     Read the AEROSOL parameters:
*-----------------------------------------------------------------------
      READ (iunin,*,err=999,end=998) naer
      IF (naer.EQ.0) THEN
         READ(iunin,*,err=999,end=998)
      ELSE
        DO i=1,naer
           READ(iunin,'(A16)') coefile(i)
        ENDDO
      ENDIF

      DO i=1,4
         READ(iunin,*,err=999,end=998)
      ENDDO

*-----------------------------------------------------------------------
*     Read the LAYER parameters:
*-----------------------------------------------------------------------
      READ (iunin,*,err=999,end=998) wav0
      READ (iunin,*,err=999,end=998) nlevs
      READ (iunin,*,err=999,end=998) nlev0

      IF (nlev0.GE.nlevs) THEN
         WRITE(*,*) 'rdinput: nlev0 should be smaller than nlevs!'
         STOP
      ENDIF
      IF (nlev0.LT.1) THEN
         WRITE(*,*) 'rdinput: nlev0 should be larger than zero!'
         STOP
      ENDIF
      IF (nlevs.GE.nlevsMAX) THEN
         WRITE(*,*) 'rdinput: nlevs should be smaller than ',nlevsMAX,
     .              ' !'
         STOP
      ENDIF

      DO i=1,3
         READ(iunin,*,err=999,end=998)
      ENDDO

      DO i=1,nlevs
         READ(iunin,*,err=999,end=998) j,baer(i),(hgpar(i,k),k=1,4),
     .       iaer(i),pres(i),temp(i),xx(1,i),xx(2,i),xx(3,i),xx(4,i),
     .       xx(5,i),xx(6,i),xx(7,i)
         xx(1,i)= xx(1,i)/1.D6
         xx(2,i)= xx(2,i)/1.D6
         xx(3,i)= xx(3,i)/1.D6
         xx(4,i)= xx(4,i)/1.D6
         xx(5,i)= xx(5,i)/1.D6
         xx(6,i)= xx(6,i)/1.D6
         xx(7,i)= xx(7,i)/1.D6
!         xx(8,i)= xx(8,i)/1.D6
!         xx(9,i)= xx(9,i)/1.D6
         pres(i)= pres(i)*1.D5

         IF (iaer(i).LT.0 .OR. iaer(i).GT.naer) iaer(i)=0
      ENDDO


*-----------------------------------------------------------------------
*     Close the input file:
*-----------------------------------------------------------------------
      CLOSE(iunin)

* make filename extention for the fourier files:
!  part1:  asurfile        part2_part3: cloud layer _ optical thickness

      fourext1=asurfile(1:3)

      DO i=1,nlevs
       IF (baer(i).gt.0.) THEN
        WRITE (fourext2, '(I2.2)') i
        WRITE (fourext3, '(I3.3)') INT(baer(i))
        GO TO 222
      ELSE
        WRITE (fourext2, '(I2.2)') 0
        WRITE (fourext3, '(I3.3)') INT(baer(1))
      ENDIF

      ENDDO
222    continue
      
       fourext=trim(fourext1)//'_'//trim(fourext2)//'_'//trim(fourext3)



************************************************************************
      RETURN
998   WRITE(*,*) 'rdinput: unexpected end of file encountered'
      STOP
999   WRITE(*,*) 'rdinput: ERROR reading file'
      STOP
      END
