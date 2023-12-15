      PROGRAM main

************************************************************************
*     This version is made for reflection by terrestrial-type
*     extrasolar planets
*
*     Daphne Stam, July 2004
*
*     wavels is in microns
************************************************************************
      IMPLICIT NONE

      INCLUDE 'max_incl'

      INCLUDE 'atm_incl'

*-----------------------------------------------------------------------
*     nlays: the number of atmospheric layers above the "ground"
*-----------------------------------------------------------------------
      INTEGER i,j,k,m,iw,ik,nmat,nwavels,naer,nphase,nplanet,iabs,
     .        nmug,nlevs,nlev0,nlays,ias,ifa,interp,nk,iso,iv,iwl,res

      INTEGER ncoefs(nlevsMAX),ncoefsa(nlevsMAX),iaer(nlevsMAX),
     .        ncoes(ncoesMAX),
     .        icoes(ncoesMAX,nwavelsMAX,2),icoes0(ncoesMAX,2),
     .        nkdis(nwavelsMAX)

      INTEGER iunss
      PARAMETER (iunss=22)

      INTEGER    Nbelow,Nall,Nmu,Nabove
      PARAMETER (Nbelow=2,Nall=3,Nmu=2,Nabove=1)

      DOUBLE PRECISION wavel,wav0,c1,c2,com,coa,asurf,grav,wa,
     .                 tot_msca,tot_mabs,tot_asca,tot_aabs,ww,
     .                 xa,a0,a1,facb,wav

      DOUBLE PRECISION wavels(nwavelsMAX),asurfs(nwavelsMAX),
     .                 pres(nlevsMAX),temp(nlevsMAX),
     .                 xx(7,nlevsMAX),
     .                 bmsca(nlevsMAX),bmabs(nlevsMAX),
     .                 basca(nlevsMAX),baabs(nlevsMAX),
     .                 baer(nlevsMAX),caer(nlevsMAX),
     .                 depols(nwavelsMAX),nmol(nlevsMAX)

      DOUBLE PRECISION coefs(nmatMAX,nmatMAX,0:ncoefsMAX,nlevsMAX),
     .                 coefsm(nmatMAX,nmatMAX,0:ncoefsMAX),
     .                 coefsa(nmatMAX,nmatMAX,0:ncoefsMAX,nlevsMAX)

      DOUBLE PRECISION a(nlevsMAX),b(nlevsMAX),wg(nkMAX)

      DOUBLE PRECISION abs0(nwavelsMAX,nkMAX),
     .                 abs1(nlevsMAX,nwavelsMAX,nkMAX)

      DOUBLE PRECISION coefs0,acoefs0,
     .                 P11(nphaseMAX),P12(nphaseMAX),
     .                 aP11(nphaseMAX),aP12(nphaseMAX)

      CHARACTER*20 wavfile,absfile,asurfile,
     .             coefile(ncoesMAX)

      CHARACTER*40 coes(ncoesMAX,nwavelsMAX)


      CHARACTER*5 atst(4)
*__________________________________________________________________
      INTEGER count_0, count_1, count_rate, count_max
      DOUBLE PRECISION startt,endt

      CHARACTER filenames*120,fourext*120,wavs*4,fourfs*120,fourfil*120

*__________________________________________________________________


*-----------------------------------------------------------------------
*     THE PROGRAM STARTS HERE:
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*     Read the standard input file:
*-----------------------------------------------------------------------

      CALL rdinput(nmat,asurf,naer,coefile,filenames,fourext,
     .             nmug,nlevs,nlev0,
     .             wavfile,res,iabs,absfile,wav0,baer,iaer,pres,temp,xx,
     .             nphase,interp,ias,ifa,asurfile,nplanet)

      nlays= nlevs-nlev0


*-----------------------------------------------------------------------
*     Create the output folder:                     
*-----------------------------------------------------------------------
      call system ( "mkdir " // filenames )

*-----------------------------------------------------------------------
*     Copy the dap.in input file to output folder:                                                                                                    
*-----------------------------------------------------------------------
                                                                       
      CALL system("cp " //"dap.in"//" "//filenames)

*-----------------------------------------------------------------------
*     Determine the planetary parameters:
*-----------------------------------------------------------------------
*     The gravitational constant (in m s^-2):
      IF (nplanet.EQ.2) grav= 9.15D0
      IF (nplanet.EQ.3) grav= 9.78D0
      IF (nplanet.EQ.4) grav= 3.69D0
      IF (nplanet.EQ.5) grav= 9.81D0

*     The molecular mass in the atmosphere (in kg.molecule^-1):
      IF (nplanet.EQ.2) wa= wCO2*xCO2_venus + wN2*xN2_venus
      IF (nplanet.EQ.3) wa= wO2*xO2_earth + wN2*xN2_earth
      IF (nplanet.EQ.4) wa= wCO2*xCO2_mars + wN2*xN2_mars +
     .                      wAr40*xAr40_mars
      IF (nplanet.EQ.5) wa= 28.3
      wa= wa*1.D-3/avogad


*     The mixing ratio of the main absorber:
      IF (nplanet.EQ.2) xa= xCO2_venus
      IF (nplanet.EQ.3) xa= xO2_earth
      IF (nplanet.EQ.4) xa= xCO2_mars
      IF (nplanet.EQ.5) xa= 0.0195297D0

!      WRITE(*,*) 'Isotropic scattering? (0=no,1=yes)'
!      WRITE(*,*) 'set to no!!!'
      iso=0

*-----------------------------------------------------------------------
*     Open the optical thickness output file:
*-----------------------------------------------------------------------
      OPEN(unit=iunss,file='tau.dat')
      WRITE(iunss,304)

*-----------------------------------------------------------------------
*     Read the wavelength file:
*-----------------------------------------------------------------------
      write(*,'(A16)') wavfile

      CALL rdwavfile(wavfile,wavels,nwavels)

*-----------------------------------------------------------------------
*     Read the surface albedo file:
*-----------------------------------------------------------------------

      CALL rdsurfile(ias,asurf,asurfile,wavels,nwavels,asurfs)

*-----------------------------------------------------------------------
*     Get the depolarization factor:
*-----------------------------------------------------------------------
      CALL getdepol(wavels,nwavels,depols)
!      WRITE(*,*) depols

*-----------------------------------------------------------------------
*     Get the number of molecules per layer:
*-----------------------------------------------------------------------
      CALL nmolecules(nlev0,nlays,grav,wa,pres,temp,nmol)

*-----------------------------------------------------------------------
*     Read the expansion coefficients data file:
*-----------------------------------------------------------------------
      CALL rdcoefile(coefile,naer,wav0,wavels,nwavels,coes,ncoes,
     .               icoes0,icoes)

*-----------------------------------------------------------------------
*     Get the aerosol optical properties at wav0:
*-----------------------------------------------------------------------

!       iaer=1
!       print*,'here2:',wav0


!        stop

      CALL getwav0(wav0,naer,coes,icoes0,nlev0,nlays,iaer,baer,caer)

*-----------------------------------------------------------------------
*     Read the gaseous absorption data files:
*-----------------------------------------------------------------------


      CALL rdabsfile(nplanet,nkdis,nwavels,wavels,iabs,absfile,xa,wg,
     .               abs0)

*-----------------------------------------------------------------------
*     Get the trace gas absorption:
*-----------------------------------------------------------------------
      CALL traces(nkdis,nwavels,wavels,nlays,xx,abs1)

************************************************************************
*     Loop over the wavelengths:
************************************************************************
      iv=0

      WRITE(*,*) 'Warning: loop bmsca!'
      facb=1.0D0

      fourfs=trim(fourext)


      DO iw=1,nwavels,res

         iv= iv+1
         wavel= wavels(iw)
         asurf= asurfs(iw)

         WRITE(*,*) wavel
         WRITE (wavs, '(I4.4)') INT(wavel*1000.)

!         fourext=trim(fourfs)//'_'//wavs
         fourext=trim(filenames)//'/four_'//trim(fourfs)//'_'//wavs

         coefs0= 0.D0
         DO i=1,nphaseMAX
            P11(i)= 0.D0
            P12(i)= 0.D0
         ENDDO

*-----------------------------------------------------------------------
*        Initialise some things:
*-----------------------------------------------------------------------
         CALL init(coefs,coefsm,coefsa)

*-----------------------------------------------------------------------
*        Calculate the molecular scattering optical thicknesses,
*        and the expansion coefficients:
*-----------------------------------------------------------------------
         CALL bmolecules(iw,wavel,nlays,nmol,depols,
     .                   bmsca,coefsm)

*-----------------------------------------------------------------------
*        Calculate the aerosol scattering and absorption optical
*        thicknesses, and the expansion coefficients:
*-----------------------------------------------------------------------
         CALL baerosols(iw,wavel,nlev0,nlays,naer,caer,iaer,coes,icoes,
     .                  nmat,basca,baabs,coefsa,ncoefsa)


*-----------------------------------------------------------------------
*        Calculate the combined expansion coefficients:
*-----------------------------------------------------------------------
         DO i=1,nlays
            ncoefs(i)=MAX0(ncoefsa(i),2)

            IF (iso.EQ.1) ncoefs(i)=2

            DO j=1,nmat
            DO k=1,nmat

               DO m=0,ncoefs(i)
                  com= bmsca(i)*coefsm(j,k,m)
                  coa= basca(i)*coefsa(j,k,m,i)
                  IF ((bmsca(i)+basca(i)).LT.1.D-10) THEN
                     coefs(j,k,m,i)= 0.D0
                  ELSE
                     coefs(j,k,m,i)= (com + coa)/(bmsca(i)+basca(i))
                  ENDIF
!                  write(*,*) com, coa, bmsca(i), basca(i)
               ENDDO
               IF (iso.EQ.1) THEN
                  coefs(j,k,0,i)= 1.D0
                  coefs(j,k,1,i)= 0.D0
                  coefs(j,k,2,i)= 0.D0
               ENDIF
            ENDDO
            ENDDO
         ENDDO

!         stop

*-----------------------------------------------------------------------
*        Check the number of k-distribution intervals:
*-----------------------------------------------------------------------
         nk= nkdis(iw)
         a0= abs0(iw,nk)
         a1= abs1(1,iw,nk)
         IF (a0.EQ.0.D0 .AND. a1.EQ.0.D0) THEN
            nk= 1
         ENDIF
         IF (wavel.LT.0.49D0) THEN
            nk=1
         ENDIF


         IF (nk.eq.0) nk=1

*-----------------------------------------------------------------------
*        Loop over the k-distribution intervals:
*-----------------------------------------------------------------------
         DO ik=1,nk
*-----------------------------------------------------------------------
*           Calculate the total optical thickness and
*           the single scattering albedo for the new k-coefficient:
*-----------------------------------------------------------------------
            tot_msca= 0.D0
            tot_mabs= 0.D0
            tot_asca= 0.D0
            tot_aabs= 0.D0

            DO i=1,nlevsMAX
              a(i)=0.
              b(i)=0.
            ENDDO


            DO i=1,nlays

               bmabs(i)= nmol(i)*abs0(iw,ik) + nmol(i)*abs1(i,iw,ik)

               b(i)=  bmabs(i) + bmsca(i) + basca(i) + baabs(i)

               IF (b(i).EQ.0.D0) THEN
                  a(i)= 0.D0
               ELSE
                  a(i)= (bmsca(i)+basca(i))/b(i)
               ENDIF

               tot_msca= tot_msca + bmsca(i)
               tot_asca= tot_asca + basca(i)
               tot_mabs= tot_mabs + bmabs(i)
               tot_aabs= tot_aabs + baabs(i)
!               write(*,*) i, tot_msca, tot_asca, tot_mabs, tot_aabs

            ENDDO
!            stop

            WRITE(iunss,302) wavel,tot_msca,tot_mabs,tot_asca,tot_aabs

*-----------------------------------------------------------------------
*           Call the doubling-adding routine:
*-----------------------------------------------------------------------

            CALL adding(fourext,a,b,coefs,ncoefs,nlays,nmug,nmat,
     .                  asurf,ifa)

*-----------------------------------------------------------------------
*           Integrate over the disk:
*-----------------------------------------------------------------------

            CALL disk(fourext,nphase,interp,aP11,aP12,acoefs0)

*-----------------------------------------------------------------------                                                                                                  
*           Move file with the Fourier coefficients to output folder:                                                                                                     
*-----------------------------------------------------------------------                                                                                                  

!            fourfil='four_'//trim(fourext)//'.dat'
!            CALL system("mv " //fourfil//" "//filenames)
!            OPEN(unit=23,file=trim(fourext)//'.dat')
*            fourfil='four_'//trim(fourext)
*            OPEN(unit=23,file=trim(fourfil)//'.dat')
!            OPEN(unit=23,file='four.dat')
*            CLOSE(23,status='delete')


*-----------------------------------------------------------------------
*           Add the k-distribution results:
*-----------------------------------------------------------------------
            IF (nk.EQ.1) THEN
               ww= 1.D0
            ELSE
               ww= wg(ik)
            ENDIF

            coefs0= coefs0 + ww*acoefs0
            DO i=1,2*((nphase+1)/2)
               P11(i)= P11(i) + ww*aP11(i)
               P12(i)= P12(i) + ww*aP12(i)
            ENDDO


222         CONTINUE
*-----------------------------------------------------------------------
*        End of loop over the k-distribution:
*-----------------------------------------------------------------------
         ENDDO

*-----------------------------------------------------------------------
*        Write results to the output files:
*-----------------------------------------------------------------------
         CALL wrout(filenames,iv,wavel,nwavels,nphase,coefs0,P11,P12)

*-----------------------------------------------------------------------
*     End of the loop over the wavelengths:
*-----------------------------------------------------------------------
111      CONTINUE
      ENDDO

*-----------------------------------------------------------------------
*     Close the output file:
*-----------------------------------------------------------------------
      CLOSE(iunss)

*-----------------------------------------------------------------------
*     Formats:
*-----------------------------------------------------------------------
300   FORMAT(' Calculation ',I4,', wavelength (in microns):',F12.8)
302   FORMAT(F12.8,1X,F14.8,1X,F14.8,1X,F14.8,1X,F14.8)
304   FORMAT('# wavelength     bmsca          bmabs',
     .       '          basca          baabs')

************************************************************************
100   WRITE(*,*)
      STOP 'end of main'
      END
