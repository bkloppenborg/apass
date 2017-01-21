      FUNCTION NBLANK (string)
c
c find first blank character
c
      INTEGER
     *       nblank, i
      CHARACTER*1
     *       string(*)
c
      nblank = 0
      DO i=1,256
        nblank = nblank + 1
        IF (string(nblank).eq.' ') GOTO 100
      ENDDO
100   continue
      nblank = nblank - 1
      RETURN
      END
