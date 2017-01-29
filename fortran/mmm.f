      SUBROUTINE  MMM (SKY, NSKY, HIBAD, READNS, SKYMN, SKYMED, 
     .     SKYMOD, SIGMA, SKEW)
C
C=======================================================================
C
C               Official DAO version:  1988 July 1
C
C This version of MMM (modified by PBS 1984.IV.10ff) assumes that
C the sky brightnesses in the one-dimensional array SKY are already
C sorted on entering this routine, and that pixels outside the "bad"
C limits have already been eliminated.
C
C This particular version of MMM also takes cognizance of the fact that,
C pixels falling below the LOBAD threshold already having been 
C eliminated, the contaminated sky pixels values overwhelmingly display
C POSITIVE departures from the true value.
C
C If for some reason it is impossible to obtain the mode of the sky
C distribution, this will be flagged by setting SIGMA = -1.0.
C
C Arguments
C
C     SKY (INPUT) is a real vector containing actual sorted sky values.
C    NSKY (INPUT) is the number of defined elements in SKY.
C  SKYMOD (OUTPUT) is the estimated mode of the sky values.
C   SIGMA (OUTPUT) is the computed standard deviation of the peak in
C         the sky histogram.
C    SKEW (OUTPUT) is the computed skewness of the peak in the sky
C         histogram.
C
C=======================================================================
C
      IMPLICIT NONE
      INTEGER NSKY
      REAL SKY(NSKY)
C
      DOUBLE PRECISION DSQRT, DBLE
      REAL ALOG10, AMIN1, AMAX1
C
      DOUBLE PRECISION SUM,SUMSQ
      REAL CUT, CUT1, CUT2, DELTA, SKYMID, SKYMED, SKYMN, SKYMOD
      REAL SIGMA, SKEW, R, SIGN, HIBAD, CENTER, SIDE, READNS
      REAL DMOD, OLD, CLAMP
      INTEGER I, J, K, L, M
      INTEGER MINIMM, MAXIMM, NITER, ISTEP, MAXIT, MINSKY, JSTEP
      LOGICAL REDO
      DATA MAXIT / 30 /, MINSKY / 20 /
C
C-----------------------------------------------------------------------
C
C SECTION 1
C
      IF (NSKY .LE. 0) THEN
         GO TO 9900
      END IF
      SKYMID=0.5*(SKY((NSKY+1)/2)+SKY(NSKY/2+1))
C
C SKYMID is the median value for the whole ensemble of sky pixels.
C Notice that if NSKY is odd, then (NSKY+1)/2 and (NSKY/2)+1 will be the
C same number, whereas if NSKY is even, they will be different numbers.
C This same trick will be used again later.
C
C Initialize the variables for accumulating the mean and standard
C deviation, and initialize the rejection limits.
C
      SUM=0.D0
      SUMSQ=0.D0
      CUT1=AMIN1(SKYMID-SKY(1), SKY(NSKY)-SKYMID, HIBAD-SKYMID)
C
C For the first pass we will consider only pixels in a symmetric 
C interval of brightness values about the median value.  This exploits
C the assumption that all the bad pixels are already rejected from the
C lower end of the brightness range.
C
      CUT2=SKYMID + CUT1
      CUT1=SKYMID - CUT1
C
      MINIMM=0
      DO 1010 I=1,NSKY
         IF (SKY(I) .LT. CUT1) THEN
            MINIMM=I
            GO TO 1010
         END IF
         IF (SKY(I) .GT. CUT2) GO TO 1020
         DELTA=SKY(I)-SKYMID
         SUM=SUM+DELTA
         SUMSQ=SUMSQ+DELTA**2
         MAXIMM=I
 1010 CONTINUE
C
C Henceforth in this subroutine, MINIMM will point to the highest value
C rejected at the lower end of the vector, and MAXIMM will point to the
C highest value accepted at the upper end of the vector.
C MAXIMM-MINIMM is the number of pixels within the acceptance range.
C
C Compute mean and sigma (from the first pass).
C
 1020 CONTINUE
      SKYMED=0.5*(SKY((MINIMM+MAXIMM+1)/2)+SKY((MINIMM+MAXIMM)/2+1))
      SKYMN=SUM/DBLE(MAXIMM-MINIMM)
      SIGMA=DSQRT(SUMSQ/DBLE(MAXIMM-MINIMM)-SKYMN**2)
      SKYMN=SKYMN+SKYMID
C
C The middle sky value, SKYMID, was subtracted off up above and added 
C back in down here to reduce the truncation error in the computation 
C of SIGMA.
C Note that this definition of SIGMA is incorrect by a factor of
C SQRT [NSKY/(NSKY-1.)], but for all but pathological cases (where none
C of this can be trusted anyway), it's close enough.
C
      SKYMOD=SKYMN
      IF (SKYMED .LT. SKYMN) SKYMOD=3.*SKYMED-2.*SKYMN
C
C If the mean is less than the mode, that means the contamination is
C slight, and the mean value is what we really want.  Note that this
C introduces a slight bias toward underestimating the sky when
C the scatter in the sky is caused by random fluctuations rather than
C by contamination, but I think this bias is negligible compared to the
C problem of contamination.
C
C-----------------------------------------------------------------------
C
C SECTION 2
C
C Rejection and recomputation loop:
C
      NITER=0
      OLD = 0.
      CLAMP = 1.
 2000 NITER=NITER+1
      IF ((NITER .GT. MAXIT) .OR. (MAXIMM-MINIMM .LT. MINSKY)) THEN
         GO TO 9900
      END IF
C
C Compute Chauvenet rejection criterion.
C         
      R=ALOG10(FLOAT(MAXIMM-MINIMM))
      R=AMAX1(2., (-.1042*R+1.1695)*R+.8895)
C
C Compute rejection limits (symmetric about the current mode).
C
      CUT=R*SIGMA+0.5*ABS(SKYMN-SKYMOD)
      CUT=AMAX1(1.5,CUT)
      CUT1=SKYMOD-CUT
      CUT2=SKYMOD+CUT
C
C Recompute mean and sigma by adding and/or subtracting sky values
C at both ends of the interval of acceptable values.
C
C At each end of the interval, ISTEP will show the direction we have to 
C step through the vector to go from the old partition to the new one.
C Pixels are added or subtracted depending upon whether the limit is 
C moving toward or away from the mode.
C
      REDO=.FALSE.
C
C Is CUT1 above or below the minimum currently-accepted value?
C
      ISTEP=INT(SIGN(1.0001, CUT1-SKY(MINIMM+1)))
      JSTEP=(ISTEP+1)/2
C
C If ISTEP = +1, JSTEP = 1.  If ISTEP = -1, JSTEP=0.  If ISTEP = +1, 
C then we know that at least one pixel must be deleted at the low end.
C
      IF (ISTEP .GT. 0) GO TO 2120
 2100 IF ((ISTEP .LT. 0) .AND. (MINIMM .LE. 0)) GO TO 2150
C
C Quit when SKY(MINIMM) < CUT1 <= SKY(MINIMM+1)
C
      IF ((SKY(MINIMM) .LE. CUT1) .AND. (SKY(MINIMM+1) .GE. CUT1))
     .     GO TO 2150
C
C If ISTEP is positive, subtract out the sky value at MINIMM+1; if 
C ISTEP is negative, add in the sky value at MINIMM.
C
 2120 CONTINUE
      DELTA=SKY(MINIMM+JSTEP)-SKYMID
      SUM=SUM-REAL(ISTEP)*DELTA
      SUMSQ=SUMSQ-REAL(ISTEP)*DELTA**2
      MINIMM=MINIMM+ISTEP
      REDO=.TRUE.                                 ! A change has occured
      GO TO 2100
C
 2150 CONTINUE
C
C Is CUT2 above or below the current maximum?
C
      ISTEP=INT(SIGN(1.0001, CUT2-SKY(MAXIMM)))
      JSTEP=(ISTEP+1)/2
C
C If ISTEP = +1, JSTEP = 1.  If ISTEP = -1, JSTEP=0.  If ISTEP = -1, 
C then we know that we must subtract at least one pixel from the high 
C end.
C
      IF (ISTEP .LT. 0) GO TO 2220
 2200 IF ((ISTEP .GT. 0) .AND. (MAXIMM .GE. NSKY)) GO TO 2250
C
C Quit when SKY(MAXIMM) <= CUT2 < SKY(MAXIMM+1)
C
      IF ((SKY(MAXIMM) .LE. CUT2) .AND. (SKY(MAXIMM+1) .GE. CUT2))
     .     GO TO 2250
C
C If ISTEP is positive, add in the sky value at MAXIMM+1; if ISTEP is 
C negative, subtract off the sky value at MAXIMM.
C
 2220 DELTA=SKY(MAXIMM+JSTEP)-SKYMID
      SUM=SUM+REAL(ISTEP)*DELTA
      SUMSQ=SUMSQ+REAL(ISTEP)*DELTA**2
      MAXIMM=MAXIMM+ISTEP
      REDO=.TRUE.                                 ! A change has occured
      GO TO 2200
C
 2250 CONTINUE
C
C Compute mean and sigma (from this pass).
C
      SKYMN=SUM/DBLE(MAXIMM-MINIMM)
      SIGMA=DSQRT(SUMSQ/DBLE(MAXIMM-MINIMM)-SKYMN**2)
      SKYMN=SKYMN+SKYMID
C
C Obtain the median.  To first approximation, the median would be the
C value of the sky in the middle pixel in the sorted data (if the
C total number is odd) or the mean of the two pixels straddling
C the middle (if the total number of pixels is even).
C
C     SKYMED=0.5*(SKY((MINIMM+MAXIMM+1)/2)+SKY((MINIMM+MAXIMM)/2+1))
C
C However, this is not good enough.  If you look at the estimator for
C the mode, you will note that a tiny change in the list of sky pixels,
C just sufficient to alter the median value of the sky brightness by
C one unit, will change the estimator of the mode by three units.  We
C really want something more robust than this.  As a first attempt
C at a more robust median estimator, I estimated the median
C of the distribution by the mean of the central twenty percent of sky
C values.  This involved considerable care to make sure you get
C a perfectly symmetric sample of pixels about the median, whether
C there is an even or an odd number of pixels within the acceptance
C interval.  However, even this is not good enough when you have 
C quantized intensities and detectors with low enough read noise that 
C a large fraction of the sky pixels have the same (quantized) intensity:  
C a given fixed percentage on either side of the median does not 
C necessarily give a good estimator of how the values differ to the two 
C sides. So I add another criterion:  if either limiting value is not
C sufficiently different from the median, expand the range by one pixel 
C on either side until we reach a significantly different intensity value.
C In a normal distribution, the central 20% corresponds to +/- 0.25 sigma,
C so expand the range until both limiting values differ from the central
C value by at least 0.25 times the read noise.
C
      SKYMED=0.0
      CENTER = REAL(MINIMM+1 + MAXIMM)/2.
      SIDE = REAL(NINT(0.2*REAL(MAXIMM-MINIMM)))/2. + 0.25
      J = NINT(CENTER-SIDE)
      K = NINT(CENTER+SIDE)
      L = NINT(CENTER-0.25)
      M = NINT(CENTER+0.25)
      R = 0.25*READNS
 2305 CONTINUE
      IF ((J .GT. 1) .AND. (K .LT. NSKY) .AND. (
     .     (SKY(L)-SKY(J) .LT. R) .OR. (SKY(K)-SKY(M) .LT. R) )) THEN
         J = J-1
         K = K+1
         GO TO 2305
      END IF
C
      DO 2310 I=J,K
 2310 SKYMED=SKYMED+SKY(I)
C
      SKYMED=SKYMED/REAL(K-J+1)
      IF (SKYMED .LT. SKYMN) THEN
         DMOD=3.*SKYMED-2.*SKYMN - SKYMOD
      ELSE
         DMOD=SKYMN - SKYMOD
      END IF
C
C If the mean is less than the mode, that means the contamination is
C slight, and the mean value is what we really want.  Note that this
C introduces a slight bias toward underestimating the sky when
C the scatter in the sky is caused by random fluctuations rather than
C by contamination, but I think this bias is negligible compared to the
C problem of contamination.
C
C If the limits have not yet stopped moving, try again.
C
      IF (DMOD*OLD .LT. 0.) CLAMP = 0.5*CLAMP
      SKYMOD = SKYMOD + CLAMP*DMOD
      OLD = DMOD
      IF (REDO) GO TO 2000
C
C-----------------------------------------------------------------------
C
C Normal return.
C
      SKEW=(SKYMN-SKYMOD)/AMAX1(1., SIGMA)
      NSKY=MAXIMM-MINIMM
      RETURN
C
C-----------------------------------------------------------------------
C
C An error condition has been detected.
C
 9900 SIGMA=-1.0
      SKEW=0.0
      RETURN
C
      END!
