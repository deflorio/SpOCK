C***********************************************************************
C
C
C      PROGRAM MAGPOINT (GEOMAG DRIVER)
C
C
C***********************************************************************
C     Contact Information
C
C     Software and Model Support
C     	National Geophysical Data Center
C     	NOAA EGC/2
C     	325 Broadway
C     	Boulder, CO 80303 USA
C     	Attn: Susan McLean or Stefan Maus
C     	Phone:  (303) 497-6478 or -6522
C     	Email:  Susan.McLean@noaa.gov or Stefan.Maus@noaa.gov
C		Web: http://www.ngdc.noaa.gov/seg/WMM/
C
C     Sponsoring Government Agency
C	   National Geospatial-Intelligence Agency
C    	   PRG / CSAT, M.S. L-41
C    	   3838 Vogel Road
C    	   Arnold, MO 63010
C    	   Attn: Craig Rollins
C    	   Phone:  (314) 263-4186
C    	   Email:  Craig.M.Rollins@Nga.Mil
C
C      Original Program By:
C        Dr. John Quinn
C        FLEET PRODUCTS DIVISION, CODE N342
C        NAVAL OCEANOGRAPHIC OFFICE (NAVOCEANO)
C        STENNIS SPACE CENTER (SSC), MS 39522-5001
C
C***********************************************************************
C
C
C     PURPOSE: THIS ROUTINE WILL COMPUTE THE GEOMAGNETIC TOTAL
C              INTENSITY, DECLINATION, AND INCLINATION AND THEIR
C              CORRESPONDING GEOMAGNETIC SECULAR VARIATIONS AT A
C              SINGLE POSITION (ALT,LAT,LON) AND TIME (DECIMAL YEARS).
C
C
C***********************************************************************
C
C
C     NOTE: THE USER IS AUTOMATICALLY PROMPTED FOR INPUT.
C
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER EPOCH2, MAXDEG
C
      REAL NaN
      REAL EPOCH, PI, DTR 
      REAL warn_H_val, warn_H_strong_val
      REAL DLAT, DLON, ALTM, ALT, TIM, TIME, DIFF
      REAL DEC, DIP, TI, GV, TIME2, DEC2, DIP2, TI2
      REAL X, X2,Y, Y2, Z, Z2, H, H2
      REAL ATI, AX, AY, AZ, AH, ADEC, ADIP
      REAL DMIN, IMIN, DDEG, IDEG
C
      LOGICAL warn_H, warn_H_strong, warn_P
C
      CHARACTER*1 ANS
      CHARACTER*2 IPPRNT(8)
      CHARACTER*20 MODEPRNT(2)
      CHARACTER ANSWER, DLOND, DLATD
      CHARACTER*6 DIPD, DECD,  ADIPD, ADECD
C
      DATA IPPRNT  /'D','I','F','H','X','Y','Z','GV'/
      DATA MODEPRNT /'Main Field','Secular Variation'/
C
      ANSWER  = 'N'
      ANS = 'Y'
      PI  = 3.14159265359
      DTR = PI/180.0
      MAXDEG=12
      NaN = ALOG(-1.0)

      CALL GEOMAG(MAXDEG,EPOCH)
      EPOCH2 = EPOCH
      PRINT *, '      ' 
      PRINT 2,EPOCH2
    2 FORMAT (' Welcome to the World Magnetic Model (WMM)',
     +     I5,' Fortran Point Program',/)
      PRINT *,'            --- Version 2.0, September 2005 ---'
      PRINT *, '      '
      PRINT *,'This program estimates the strength and direction of',
     +     ' Earth''s main'
      PRINT *, 'magnetic field for a given location.'
      PRINT *, '      '
      PRINT *, '      '
      PRINT*,'Enter h for help and contact information or c to ' //
     +     'continue.'
      PRINT*, ' '
      READ 15, ANSWER

      IF (ANSWER .EQ. 'H' .OR. ANSWER .EQ. 'h') THEN
      PRINT *, '      '
      PRINT *, 'Input required is: '
      PRINT *, '     1) The location in geodetic latitude and longitude'
      PRINT *, '                 (positive for northern latitudes '
      PRINT *, '                  and eastern longitudes) '
      PRINT *, '     2) The WGS84 geodetic altitude AMSL in meters'
      PRINT *, '     3) The date of interest in years.'
      PRINT *, 'The program computes the estimated magnetic Declination'
      PRINT *, '(D sometimes called MAGVAR), Inclination (I), Total'
      PRINT *, 'Intensity (F or TI), Horizontal Intensity (H or HI),'
      PRINT *, 'Vertical Intensity (Z), and Grid Variation (GV).'
      PRINT *, 'Declination and Grid Variation are measured in units of'
      PRINT *, 'degrees and are considered positive when east of north.'
      PRINT *, 'Inclination is measured in units of degrees and is '
      PRINT *, 'considered positive when pointing down, into the Earth.'
      PRINT *, 'The WMM is reference to the WGS-84 ellipsoid and is '
      PRINT *, 'valid for 5 years after the base epoch.'
      PRINT *, '      '
      PRINT *, 'It is very important to note that a degree and order'
      PRINT *, '12 model, such as WMM, describes only the long'
      PRINT *, 'wavelength spatial magnetic fluctuations due to'
      PRINT *, 'Earth''s core.  Not included in the WMM series models'
      PRINT *, 'are intermediate and short wavelength spatial'
      PRINT *, 'fluctuations originating in Earth''s mantle and crust.'
      PRINT *, 'Consequently, isolated angular errors at various'
      PRINT *, 'positions on the surface (primarily over land, in'
      PRINT *, 'continental margins and over oceanic seamounts, ridges'
      PRINT *, 'and trenches) of several degrees may be expected.'
      PRINT *, 'Also not included in the model are temporal'
      PRINT *, 'fluctuations of magnetospheric and ionospheric origin.'
      PRINT *, 'On the days during and immediately following magnetic'
      PRINT *, 'storms, temporal fluctuations can cause substantial'
      PRINT *, 'deviations of the geomagnetic field from model values.'
      PRINT *, 'If the required declination accuracy is more stringent'
      PRINT *, 'than the WMM series of models provide, the user is'
      PRINT *, 'advised to request special (regional or local) surveys'
      PRINT *, 'be performed and models prepared. Please make requests'
      PRINT *, 'of this nature to the National Geospatial-Intelligence'
      PRINT *, 'Agency (NGA) at the address below.'
      PRINT *, '   '
      PRINT *, 'Contact Information:'
      PRINT *, '------- ----------- '
      PRINT *, '   '
      PRINT *, ' Software and Model Support'
      PRINT *, ' National Geophysical Data Center'
      PRINT *, ' NOAA EGC/2'
      PRINT *, ' 325 Broadway'
      PRINT *, ' Boulder, CO 80303 USA'
      PRINT *, ' Attn: Susan McLean or Stefan Maus'
      PRINT *, ' Phone:  (303) 497-6478 or -6522'
      PRINT *, ' Email:  Susan.McLean@noaa.gov or Stefan.Maus@noaa.gov'
      PRINT *, '       '
      PRINT *, 'Continue with program? (y or n)'
      READ 15, ANS
      ENDIF
C      
      IF (ANS .NE. 'Y' .AND. ANS .NE. 'y' ) GOTO 100
C
C
      WRITE(*,*)  '                                          '

 10   warn_H = .FALSE.
      warn_H_val = 99999.0
      warn_H_strong = .FALSE.
      warn_H_strong_val = 99999.0
      warn_P = .FALSE.
C
      WRITE(*,*)  ' ENTER LATITUDE IN DECIMAL DEGREES        '
      WRITE(*,*)  ' (North latitude positive, South negative)'
      WRITE(*,*)  ' (i.e. 25.5 for 25 degrees 30 minutes north)'
      READ(*,*) DLAT
C
      WRITE(*,*)  '                                          '
      WRITE(*,*)  ' ENTER LONGITUDE IN DECIMAL DEGREES       '
      WRITE(*,*)  '(East longitude positive, West negative)'
      WRITE(*,*)  '(i.e. -100.0 for 100 degrees west)'
      READ(*,*) DLON
C
      WRITE(*,*)  '                                          '
      PRINT *, 'ENTER ALTITUDE IN METERS ABOVE MEAN SEA LEVEL (WGS84)'
      READ  (*,*) ALTM
      ALT = ALTM/1000.
C
C
      WRITE(*,*)  '                                          '
      WRITE(*,*)  ' ENTER TIME IN DECIMAL YEAR               '
      READ *, TIME
      TIM=TIME+1.0
C
      DIFF=TIME-EPOCH
      IF ((DIFF .LT. 0. .OR. DIFF .GT. 5.)) THEN
      PRINT *, '      '
      PRINT *, ' WARNING - TIME EXTENDS BEYOND MODEL 5-YEAR LIFE SPAN'
      PRINT *, ' CONTACT NGDC FOR PRODUCT UPDATES:'
      PRINT *, '      '
      PRINT *, '  Software and Model Support'
      PRINT *, '      '
      PRINT *, '  National Geophysical Data Center'
      PRINT *, '  NOAA EGC/2'
      PRINT *, '  325 Broadway'
      PRINT *, '  Boulder, CO 80303 USA'
      PRINT *, '  Attn: Susan McLean or Stefan Maus'
      PRINT *, '  Phone:  (303) 497-6478 or -6522'
      PRINT *, '  Email:  Susan.McLean@noaa.gov or Stefan.Maus@noaa.gov'
      PRINT *, '      '
C
      PRINT *, '  EPOCH  = ',EPOCH
      PRINT *, '  TIME   = ',TIME
      PRINT *, '  Continue with program? (y or n)'
      READ 15, ANSWER
 15   FORMAT (A1)
C
      ELSE 
         ANSWER = 'Y'
      ENDIF
C
      IF (ANSWER .NE. 'Y' .AND. ANSWER .NE. 'y') GOTO 100
C
      WRITE(*,*)  '   '
      WRITE(*,*)  '   '
C
C
      CALL GEOMG1(ALT,DLAT,DLON,TIME,DEC,DIP,TI,GV,EPOCH)
C     
      TIME2 = TIME + 1.0
C
      CALL GEOMG1(ALT,DLAT,DLON,TIME2,DEC2,DIP2,TI2,GV,EPOCH)

C     COMPUTE X, Y, Z, AND H COMPONENTS OF THE MAGNETIC FIELD
      
      X  = TI  * (COS((DEC*DTR)) * COS((DIP*DTR)))
      X2 = TI2 * (COS((DEC2*DTR)) * COS((DIP2*DTR)))
      Y  = TI  * (COS((DIP*DTR)) * SIN((DEC*DTR)))
      Y2 = TI2 * (COS((DIP2*DTR)) * SIN((DEC2*DTR)))
      Z  = TI  * (SIN((DIP*DTR)))
      Z2 = TI2 * (SIN((DIP2*DTR)))
      H  = TI  * (COS((DIP*DTR)))
      H2 = TI2 * (COS((DIP2*DTR)))
C
C      COMPUTE ANNUAL RATES OF CHANGE (nT/yr)
C
C
      ATI = TI2 - TI
      AX = X2 - X
      AY = Y2 - Y
      AZ = Z2 - Z
      AH = H2 - H
C
C
C      COMPUTE ANNUAL CHANGE FOR DIP & DEC (Minutes/yr)
C
C
      ADIP = (DIP2 - DIP) * 60.
      ADEC = (DEC2 - DEC)
      IF (ADEC .GT.  180.) ADEC = ADEC - 360.
      IF (ADEC .LE. -180.) ADEC = ADEC + 360.
      ADEC = ADEC * 60.
C     
      IF (DIP .LT. 0 ) THEN
         DIPD = '( up )'
      ELSE
         DIPD = '(down)'
      ENDIF
C      
      IF (DEC .GE. 0 ) THEN
         DECD = '(east)'
      ELSE
         DECD = '(west)'
      ENDIF
C     
      IF (DLAT .LT. 0 ) THEN
         DLATD = 'S'
      ELSE
         DLATD = 'N'
      ENDIF
C     
      IF (DLON .GE. 0 ) THEN
         DLOND = 'E'
      ELSE
         DLOND = 'W'
      ENDIF
C     
      IF (ADIP .LT. 0 ) THEN
         ADIPD = '( up )'
      ELSE
         ADIPD = '(down)'
      ENDIF
C      
      IF (ADEC .GE. 0 ) THEN
         ADECD = '(east)'
      ELSE
         ADECD = '(west)'
      ENDIF
C
C     Warnings and behavior at poles
C

      if (H .LT. 100.0) THEN ! at magnetic poles
         DEC = NaN
         ADEC = NaN
         ADECD = '(void)'
         DECD = '(void)'
      END IF

      IF (H .LT. 1000.0) THEN
         warn_H = .FALSE.
         warn_H_strong = .TRUE.
         warn_H_strong_val = H
      ELSE IF (H .LT. 5000.0 .AND. (.NOT. warn_H_strong)) THEN 
         warn_H = .TRUE.
         warn_H_val = H
      END IF
            
      if (90.0-ABS(DLAT) .LE. 0.001) THEN ! at geographic poles
          X = NaN
          Y = NaN
          DEC = NaN
          AX = NaN
          AY = NaN
          ADEC = NaN
          ADECD = '(void)'
          DECD = '(void)'
          warn_P = .TRUE.
          warn_H = .FALSE.
          warn_H_strong = .FALSE.
       END IF
C
C     Convert D and I to deg and min
C
      DDEG = INT(DEC)
      DMIN = (DEC - DDEG) * 60.0
      IF (DDEG.NE.0.0) THEN
	DMIN = ABS(DMIN)
	END IF
      IDEG = INT(DIP)
      IMIN = (DIP - IDEG) * 60.0
      IMIN = ABS(IMIN)
C
C23456789.123456789.123456789.123456789.123456789.123456789.123456789.12
      WRITE(6,*) 'Results for'
      WRITE(6,41) ABS(DLAT), DLATD, ABS(DLON), DLOND
      WRITE(6,42) ALTM, TIME
      WRITE(6,*) '   Main Field                     Secular Change'
      IF (DDEG .GE. -180.0 .AND. DDEG .LE. 180) THEN
         WRITE(6,43) int(DDEG), int(DMIN), DECD, ADEC, ADECD
      ELSE
         PRINT *,'D:    NaN Deg NaN Min           dD:     NaN Min/yr'
      ENDIF
      WRITE(6,44) int(IDEG), int(IMIN), DIPD, ADIP, ADIPD
      WRITE(6,45) TI, ATI
      WRITE(6,46) H, AH
      IF (X .GE. -1.0e+6 .AND. X .LE. 1.0e+6) THEN
         WRITE(6,47) X, AX
      ELSE
         PRINT *,'X:            NaN nT            dX:     NaN nT/yr'
      END IF
      IF (Y .GE. -1.0e+6 .AND. Y .LE. 1.0e+6) THEN
         WRITE(6,48) Y, AY
      ELSE
         PRINT *,'Y:            NaN nT            dY:     NaN nT/yr'
      END IF
      WRITE(6,49) Z, AZ
      WRITE(6,*)  '                                        '
C
 41   FORMAT(/,' LATITUDE:  ',F9.2,A1,/,' LONGITUDE: ',F9.2,A1)
 42   FORMAT(  ' ALTITUDE:  ',F9.2,' METERS AMSL (WGS84)',/,   
     +         ' DATE:      ',F9.2,/)
 43   FORMAT(' D:   ',I4,' Deg ',I3,' Min ',
     +     A6,'    dD:',F9.1,' Min/yr ',A6)
 44   FORMAT(' I:   ',I4,' Deg ',I3,' Min ',
     +     A6,'    dI:',F9.1,' Min/yr ',A6)
 45   FORMAT(' F:',6X,F9.1,' nT',8X,'    dF:',F9.1,' nT/yr')
 46   FORMAT(' H:',6X,F9.1,' nT',8X,'    dH:',F9.1,' nT/yr')
 47   FORMAT(' X:',6X,F9.1,' nT',8X,'    dX:',F9.1,' nT/yr')
 48   FORMAT(' Y:',6X,F9.1,' nT',8X,'    dY:',F9.1,' nT/yr')
 49   FORMAT(' Z:',6X,F9.1,' nT',8X,'    dZ:',F9.1,' nT/yr')

      if (warn_H) WRITE(6,51) warn_H_val
      if (warn_H_strong) WRITE(6,52) warn_H_strong_val
      if (warn_P) WRITE(6,53)

 51   FORMAT('Warning: The horizontal field strength at this',
     1     ' location is only ',F6.1,' nT.',
     2     /,'         Compass readings have large uncertainties ',
     3     'in areas where H is',/,'         smaller than 5000 nT',/)
 52   FORMAT('Warning: The horizontal field strength at this',
     1     ' location is only',F6.1,' nT',
     2     /,'         Compass readings have VERY LARGE ',
     3     'uncertainties in areas where H is',/,
     4     '         smaller than 1000 nT',/)
 53   FORMAT('Warning: Location is at geographic pole where X,',
     5     'Y, and D are undefined',/)
C     
C     
      WRITE(*,*) 'DO YOU NEED MORE POINT DATA ?'
      READ(*,'(A)') ANSWER
      IF (ANSWER .EQ. 'y' .OR. ANSWER .EQ. 'Y') GO TO 10
 100  CONTINUE
      WRITE (*,*) 'Press RETURN to exit WMM Point Calculator'
      READ (*,*)
C     
      STOP
      END
