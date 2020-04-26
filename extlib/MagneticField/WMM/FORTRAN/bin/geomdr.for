C***********************************************************************
C
C
C     PROGRAM GEOMDR (GEOMAG DRIVER)
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
C     PURPOSE:  THIS ROUTINE PRODUCES AN UNFORMATTED GRID AND AN
C               ASCII-FORMATTED FILE OF THE GRIDDED VALUES OF
C               ANY ONE OF THE FOLLOWING GEOMAGNETIC PARAMETERS:
C
C                     DECLINATION (DEC)
C                     INCLINATION (DIP)
C                     TOTAL INTENSITY (TI)
C                     HORIZONTAL INTENSITY (H)
C                     NORTH COMPONENT (X)
C                     EAST COMPONENT (Y)
C                     VERTICALLY DOWN COMPONENT (Z)
C                     GRID VARIATION (GV) (ARCTIC/ANTARCTIC ONLY)
C
C               THE COMPUTATION OF THESE PARAMETERS IS BASED ON
C               A SPHERICAL HARMONIC MODEL OF THE EARTH'S MAGNETIC
C               POTENTIAL SUCH AS WMM-2000.
C
C               ONE GRIDDED FILE CALLED MAG.GRD IS UNFORMATTED AND
C               CONSISTS OF ONE HEADER RECORD CONTAINING THE
C               GEOGRAPHIC LIMITS AND GRID SPACING (GEODETIC
C               COORDINATES) OF THE DATA FOLLOWED BY A RECTANGULAR
C               ARRAY OF GRIDDED DATA WRITTEN ONE LATITUDE BAND AT
C               A TIME BEGINNING WITH THE LOWEST LATITUDE BAND.
C
C               ONE ASCII FILE, CALLED MAG.DAT, CONSISTS OF THREE
C               COLUMNS OF DATA FORMATTED AS (2(F8.2,5X),F11.3).
C               THE DATA IN COLUMNS 1, 2 AND 3 ARE GEODETIC LATITUDE,
C               GEODETIC LONGITUDE, AND MAGNETIC VALUE RESPECTIVELY.
C
C
C*********************************************************************
C
C
C     REFERENCES:
C
C       JOHN M. QUINN, DAVID J. KERRIDGE AND DAVID R. BARRACLOUGH,
C            WORLD MAGNETIC CHARTS FOR 1985 - SPHERICAL HARMONIC
C            MODELS OF THE GEOMAGNETIC FIELD AND ITS SECULAR
C            VARIATION, GEOPHYS. J. R. ASTR. SOC. (1986) 87,
C            PP 1143-1157
C
C       DEFENSE MAPPING AGENCY TECHNICAL REPORT, TR 8350.2:
C            DEPARTMENT OF DEFENSE WORLD GEODETIC SYSTEM 1984,
C            SEPT. 30 (1987)
C
C       JOHN M. QUINN, RACHEL J. COLEMAN, MICHAEL R. PECK, AND
C            STEPHEN E. LAUBER; THE JOINT US/UK 1990 EPOCH
C            WORLD MAGNETIC MODEL, TECHNICAL REPORT NO. 304,
C            NAVAL OCEANOGRAPHIC OFFICE (1991)
C
C       JOHN M. QUINN, RACHEL J. COLEMAN, DONALD L. SHIEL, AND
C            JOHN M. NIGRO; THE JOINT US/UK 1995 EPOCH WORLD
C            MAGNETIC MODEL, TECHNICAL REPORT NO. 314, NAVAL
C            OCEANOGRAPHIC OFFICE (1995)
C
C       SUSAN MACMILLAN, DAVID R. BARRACLOUGH, JOHN M. QUINN, AND
C            RACHEL J. COLEMAN; THE 1995 REVISION OF THE JOINT US/UK
C            GEOMAGNETIC FIELD MODELS - I. SECULAR VARIATION,
C            JOURNAL OF GEOMAGNETISM AND GEOELECTRICITY, VOL. 49,
C            PP. 229 - 243 (1997)
C
C       JOHN M. QUINN, RACHEL J. COLEMAN, SUSAN MACMILLAN AND DAVID
C            R. BARRACLOUGH; THE 1995 REVISION OF THE JOINT US/UK
C            GEOMAGNETIC FIELD MODELS. II: MAIN FIELD. JOURNAL OF
C            GEOMAGNETISM AND GEOELECTRICITY, VOL. 49, PP. 245 - 261
C            (1997)
C
C
C***********************************************************************
C
C
C     PARAMETER DESCRIPTIONS:
C
C       ALT    - WGS84 GEODETIC ALTITUDE AMSL (KM)
C       GLAT   - WGS84 GEODETIC LATITUDE (DEG.)
C       GLON   - GEODETIC LONGITUDE (DEG.)
C       IPARAM - GEOMAGNETIC PARAMETER SELECTION FLAG
C                1 = DECLINATION (DEG.)                       SV (MIN/YR)
C                2 = INCLINATION (DEG.)                       SV (MIN/YR)
C                3 = TOTAL INTENSITY (nT)                     SV ( nT/YR)
C                4 = HORIZONTAL INTENSITY (nT)                SV ( nT/YR)
C                5 = NORTH COMPONENT (nT)                     SV ( nT/YR)
C                6 = EAST COMPONEENT (nT)                     SV ( nT/YR)
C                7 = VERTICALLY DOWN COMPONENT (nT)           SV ( nT/YR)
C                8 = GRID VARIATION (DEG.) (ARCTIC/ANTARCTIC) SV (MIN/YR)
C       MODE   - FIELD TYPE
C                1 = MAIN FIELD (MF)
C                2 = SECULAR VARIATION (SV)
C       GLONE  - EASTERN CHART LIMIT (DEG.)
C       GLONW  - WESTERN CHART LIMIT (DEG.)
C       GLATN  - NORTHERN CHART LIMIT (DEG.)
C       GLATS  - SOUTHERN CHART LIMIT (DEG.)
C       GRID   - DISTANCE BETWEEN GRID POINTS (MIN. OF ARC)
C       NPTS   - NUMBER OF LATITUDE BANDS (ROWS) IN GRIDDED CHART
C       MPTS   - NUMBER OF LONGITUDE BANDS (COLUMNS) IN GRIDDED CHART
C
C
C**********************************************************************
C
      IMPLICIT NONE
C
      INTEGER MAX_MPTS
      PARAMETER ( MAX_MPTS=360*60 ) !global 1 minute resolution

      INTEGER EPOCH2, MAXDEG, IPARAM, MODE, IPRINT
      INTEGER NPTS, MPTS, N, M
      REAL PI, DTR, EPOCH
      REAL GLAT, GLON, DEC, DIP, TI, GV, H, X, Y, Z
      REAL ALT, GLATS, GLATN, GLONW, GLONE, GRID, TIM, TIME, DT
      REAL MAG(MAX_MPTS)

      CHARACTER*2 IPPRNT(8)
      CHARACTER*20 MODEPRNT(2)
      CHARACTER*1 GLATSD, GLATND, GLONWD, GLONED, ANSWER, ANS
      DATA IPPRNT /'D','I','F','H','X','Y','Z','GV'/
      DATA MODEPRNT /'Main Field','Secular Variation'/

C
C
      MAXDEG=12
      PI=3.14159265359
      DTR=PI/180.0
      ANSWER = 'N'
      ANS = 'Y'

      CALL GEOMAG(MAXDEG,EPOCH)
	
      EPOCH2 = EPOCH

      PRINT *, '      ' 
      PRINT 5,EPOCH2
    5 FORMAT('Welcome to the World Magnetic Model (WMM)',I5,
     +     ' Fortran Grid Program',/)
      PRINT *,'            --- Version 2.0, September 2005 ---'
      PRINT *, '      '

      PRINT *,'This program estimates the strength and direction of',
     +     ' Earth''s main'
      PRINT *, 'magnetic field for an area.'
      PRINT *, '    '
      PRINT*, '    '
      PRINT*, '    '
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
      PRINT *, 'It is very important to note that the output of this'
      PRINT *, 'program is given by 4 decimals ONLY to generate smooth'
      PRINT *, 'isoline maps. This resolution does not reflect the'
      PRINT *, 'actual accuracy of the model. A degree and order'
      PRINT *, '12 model, such as the WMM, describes only the long'
      PRINT *, 'wavelength spatial magnetic fluctuations due to sources'
      PRINT *, 'in Earth''s core. Not included in the WMM series models'
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

      PRINT *, 'Contact Information'

      PRINT *, ' Software and Model Support'
      PRINT *, ' National Geophysical Data Center'
      PRINT *, ' NOAA EGC/2'
      PRINT *, ' 325 Broadway'
      PRINT *, ' Boulder, CO 80303 USA'
      PRINT *, ' Attn: Susan McLean or Stefan Maus'
      PRINT *, ' Phone:  (303) 497-6478 or -6522'
      PRINT *, ' Email:  Susan.McLean@noaa.gov or Stefan.Maus@noaa.gov'
      PRINT *, '       '
      PRINT *, ' Sponsoring Government Agency'
      PRINT *, ' National Geospatial-Intelligence Agency'
      PRINT *, ' PRG / CSAT, M.S. L-41'
      PRINT *, ' 3838 Vogel Road'
      PRINT *, ' Arnold, MO 63010'
      PRINT *, ' Attn: Craig Rollins'
      PRINT *, ' Phone:  (314) 263-4186'
      PRINT *, ' Email:  Craig.M.Rollins@Nga.Mil'
      PRINT *, '       '
      PRINT *, 'Continue with program? (y or n)'
      READ 15, ANS
      ENDIF
      
      IF (ANS .NE. 'Y' .AND. ANS .NE. 'y' ) GOTO 100

      
C
C
      PRINT *, '      '
      PRINT *, '      '
      PRINT *, '      '
      PRINT *, 'ENTER ALTITUDE IN METERS ABOVE MEAN SEA LEVEL (WGS84)'
      READ *, ALT
      ALT = ALT / 1000.0 !CONVERT TO KM
C
C
      PRINT *, '      '
      PRINT *, 'ENTER CHART LIMITS'
      PRINT *, '(North latitude positive, South latitude negative;'
      PRINT *, ' East longitude positive, West longitude negative i.e.'
      PRINT *, ' 32,34,-118,-116 for San Diego, CA)'
      PRINT *, '      '
      PRINT *, '      '
      PRINT *, 'LATS,LATN,LONW,LONE'
      READ *, GLATS,GLATN,GLONW,GLONE
C
C
	
      PRINT *, '      '
      PRINT *, 'ENTER TIME OF CHART IN YEARS'
      READ *, TIME
      TIM=TIME+1.0

      DT=TIME-EPOCH
      IF ((DT .LT. 0. .OR. DT .GT. 5.)) THEN
      PRINT *, '      '
      PRINT *, 'WARNING - TIME EXTENDS BEYOND MODEL 5-YEAR LIFE SPAN'
      PRINT *, 'CONTACT NGDC FOR PRODUCT UPDATES:'
      PRINT *, '      '
      PRINT *, ' Software and Model Support'
      PRINT *, '      '
      PRINT *, '  National Geophysical Data Center'
      PRINT *, '  NOAA EGC/2'
      PRINT *, '  325 Broadway'
      PRINT *, '  Boulder, CO 80303 USA'
      PRINT *, '  Attn: Susan McLean or Stefan Maus'
      PRINT *, '  Phone:  (303) 497-6478 or -6522'
      PRINT *, '  Email:  Susan.McLean@noaa.gov or Stefan.Maus@noaa.gov'
      PRINT *, '      '
      PRINT *, '      '
      PRINT *, '      '
      PRINT *, '   Sponsoring Government Agency'
      PRINT *, '      '
      PRINT *, '   National Geospatial-Intelligence Agency'
      PRINT *, '   PRG / CSAT, M.S. L-41'
      PRINT *, '   3838 Vogel Road'
      PRINT *, '   Arnold, MO 63010'
      PRINT *, '   Attn: Craig Rollins'
      PRINT *, '   Phone:  (314) 263-4186'
      PRINT *, '   Email:  Craig.M.Rollins@Nga.Mil'
C
      PRINT *, '   EPOCH  = ',EPOCH
      PRINT *, '   TIME   = ',TIME
      PRINT *, '   Continue with program? (y or n)'
      READ 15, ANSWER
 15   FORMAT (A1)
C
      ELSE 
         ANSWER = 'Y'
      ENDIF
C
      IF (ANSWER .NE. 'Y' .AND. ANSWER .NE. 'y') THEN 
         GOTO 100
      ENDIF
C
 42   PRINT *, 'ENTER GRID SPACING (MIN. OF ARC)'
      READ *, GRID
C
C prevent array index value out of range
C
      If ((GLONE-GLONW)*60./GRID+1.GT.MAX_MPTS) THEN
         print *,"Increase grid spacing or increase"
         print *,"program parameter MAX_MPTS"
         goto 42
      end if
C
C
      PRINT *, '      '
      PRINT *, 'ENTER DESIRED MAGNETIC PARAMETER (1-8)'
      PRINT *, '      '
      PRINT *, '      1 = DECLINATION (D)'
      PRINT *, '      2 = INCLINATION (I)'
      PRINT *, '      3 = TOTAL INTENSITY (F)'
      PRINT *, '      4 = HORIZONTAL INTENSITY (H)'
      PRINT *, '      5 = NORTH COMPONENT (X)'
      PRINT *, '      6 = EAST COMPONENT (Y)'
      PRINT *, '      7 = VERTICALLY DOWN COMPONENT (Z)'
      PRINT *, '      8 = GRID VARIATION (GV)'
      READ *, IPARAM
C
C
      PRINT *, '      '
      PRINT *, 'ENTER DESIRED FIELD TYPE'
      PRINT *, ' 1 = MAIN FIELD'
      PRINT *, ' 2 = SECULAR VARIATION'
      READ  *, MODE
C
C
      PRINT *, '      '
      PRINT *, 'ENTER PRINT MODE'
      PRINT *, '      '
      PRINT *, '      0 = NO'
      PRINT *, '      1 = YES'
      READ  *, IPRINT
C
C
      OPEN(UNIT=14,FILE='MAG.GRD',STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(UNIT=15,FILE='MAG.DAT',STATUS='UNKNOWN')
      OPEN(UNIT=16,FILE='MAG.HDR',STATUS='UNKNOWN')
C
      IF (GLATS .GT. 0) THEN
         GLATSD = 'N'
      ELSE
         GLATSD = 'S'
      ENDIF
C     
      IF (GLATN .GT. 0) THEN
         GLATND = 'N'
      ELSE
         GLATND = 'S'
      ENDIF
C     
      IF (GLONW .GT. 0) THEN
         GLONWD = 'E'
      ELSE
         GLONWD = 'W'
      ENDIF
C     
      IF (GLONE .GT. 0) THEN
         GLONED = 'E'
      ELSE
         GLONED = 'W'
      ENDIF
C     
      IF (MODE .EQ. 1) THEN
         PRINT 50,IPPRNT(IPARAM),MODEPRNT(MODE),TIME
C        
      ELSE
C        
         PRINT 60,IPPRNT(IPARAM),MODEPRNT(MODE),TIME
C        
      ENDIF
      PRINT 70,GLATS,GLATSD,GLATN,GLATND
      PRINT 80,GLONW,GLONWD,GLONE,GLONED
C     
      IF (IPRINT .EQ. 1) THEN
         WRITE (16,50)IPPRNT(IPARAM),MODEPRNT(MODE),TIME
         WRITE (16,70)GLATS,GLATSD,GLATN,GLATND
         WRITE (16,80)GLONW,GLONWD,GLONE,GLONED
      ENDIF
C     
C
 50   FORMAT('Results For Element ',A,' (',A10,')', F7.1)
 60   FORMAT('Results For Element ',A,' (',A17,')',F7.1)
 70   FORMAT('Latitude Range:  ',F7.2,' (',A1,')  to ',F7.2,' (',A1,')') 
 80   FORMAT('Longitude Range: ',F7.2,' (',A1,')  to ',F7.2,' (',A1,')') 
C     
C     
      NPTS=(GLATN-GLATS)*60./GRID+1
      MPTS=(GLONE-GLONW)*60./GRID+1
C     
C     
      WRITE(14) GLONW,GLONE,GLATS,GLATN,GRID,MPTS,NPTS
C     
C     
C     INITIALIZE GEOMAG ROUTINE
C     
C     
C     GENERATE GRID VIA CONSECUTIVE CALLS TO THE PROCESSING
C     MODULE OF GEOMAG (ENTRY POINT GEOMG1)
C     
C     
      DO 30 N=1,NPTS
         GLAT=GLATS+(N-1)*GRID/60.
         DO 20 M=1,MPTS
            GLON=GLONW+(M-1)*GRID/60.
C           
C           
            CALL GEOMG1(ALT,GLAT,GLON,TIME,DEC,DIP,TI,GV,EPOCH)
C           
C           
            H=TI*COS(DIP*DTR)
            X=H*COS(DEC*DTR)
            Y=H*SIN(DEC*DTR)
            Z=TI*SIN(DIP*DTR)
C           
C           
            IF (IPARAM .EQ. 1) MAG(M)=DEC
            IF (IPARAM .EQ. 2) MAG(M)=DIP
            IF (IPARAM .EQ. 3) MAG(M)=TI
            IF (IPARAM .EQ. 4) MAG(M)=H
            IF (IPARAM .EQ. 5) MAG(M)=X
            IF (IPARAM .EQ. 6) MAG(M)=Y
            IF (IPARAM .EQ. 7) MAG(M)=Z
            IF (IPARAM .EQ. 8) MAG(M)=GV
C           
C           
            
            
            IF (MODE .EQ. 2) THEN
C              
C
               CALL GEOMG1(ALT,GLAT,GLON,TIM,DEC,DIP,TI,GV,EPOCH)
C
C
               H=TI*COS(DIP*DTR)
               X=H*COS(DEC*DTR)
               Y=H*SIN(DEC*DTR)
               Z=TI*SIN(DIP*DTR)
C
C
               IF (IPARAM .EQ. 1) MAG(M)=(DEC-MAG(M))
               IF (IPARAM .EQ. 2) MAG(M)=(DIP-MAG(M))
               IF (IPARAM .EQ. 3) MAG(M)=TI-MAG(M)
               IF (IPARAM .EQ. 4) MAG(M)=H-MAG(M)
               IF (IPARAM .EQ. 5) MAG(M)=X-MAG(M)
               IF (IPARAM .EQ. 6) MAG(M)=Y-MAG(M)
               IF (IPARAM .EQ. 7) MAG(M)=Z-MAG(M)
               IF (IPARAM .EQ. 8) MAG(M)=GV-MAG(M)
C              
C          The following lines prevent declinations and
C          grid variations out of range and scale
C
               IF (IPARAM.EQ.1.OR.IPARAM.EQ.2.OR.IPARAM.EQ.8) THEN
                  IF (MAG(M) .GT.  180.) MAG(M) = MAG(M) - 360.
                  IF (MAG(M) .LE. -180.) MAG(M) = MAG(M) + 360.
C     convert the time derivatives of D, I and GV to min/year
                  MAG(M) = MAG(M) * 60. 
               END IF
C
C
            ENDIF
C
C
            IF (IPRINT .EQ. 1) PRINT 90, GLAT,GLON,MAG(M)
            WRITE(15,90) GLAT,GLON,MAG(M)
 90         FORMAT(2(F8.2,2X),F11.4)
C           
C
 20      CONTINUE
         WRITE(14) (MAG(M),M=1,MPTS)
 30   CONTINUE
C     
C
      IF (IPRINT .EQ. 1) THEN
         WRITE (*,*) 'Press RETURN to exit WMM Grid Calculator'
         READ (*,*)
      END IF
C
 100  CONTINUE
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      PRINT *, 'End of Program'
      STOP
      END
