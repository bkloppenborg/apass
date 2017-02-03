#!/bin/python

# FRED files have the following format:
## STANDARD MAGNITUDES ONLY
## FILCON ver 3.0
## RA (J2000)    DEC    CCDX      CCDY     Flags      HJD   Airmass  Set     Group    Field      Filt   Mag     Error    dmag    sys night
## 102.5140180   0.2598700  1828.670    16.950 0 0 56295.536090 2.133    1         92 0020110040     8  10.7481  0.0050  0.0460    31 56295

apass_col_names = ['ra', 'dec', 'ccdx', 'ccdy', 'flag1', 'flag2', 'hjd', 'avexx', 'kset', 'group', 'star', 'filter_id', 'xmag1', 'xerr1', 'dmag', 'sys', 'night']
apass_col_types = ['float64', 'float64', 'float32', 'float32', 'bool', 'bool', 'float32', 'float32', 'int32', 'int32', 'int32', 'uint8', 'float32', 'float32', 'float32', 'int32', 'int32']
