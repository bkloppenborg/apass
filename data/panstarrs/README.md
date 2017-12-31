PANSTARRS Data
======

The data in this directory were downloaded from the PANSTARRS survey
data release 1 on 2017-12-30 using the MAST Query / CasJobs query tool:

http://mastweb.stsci.edu/ps1casjobs/home.aspx

Documentation on PANSTARRS DR1 implies saturation around 12-14 mag
depending on filter. As such, our queries are magnitude limited
to 12 < M_x < 16 to best match SRO expectations.

The following table indicates the field numbers and corresponding query
regions:

| Field | RA_min | RA_max | DEC_min | DEC_max |
|=======|========|========|=========|=========|
|   308 | 177.01 | 180.37 |   -1.51 |    1.52 |
|   309 | 278.35 | 281.65 |   -1.52 |    1.53 |
|   310 |  11.60 |  14.89 |   -1.52 |    1.56 |
|   311 | 102.34 | 105.65 |   -1.74 |    1.55 |

## Example Query

select o.objID, o.nDetections, o.raMean, o.decMean,
     m.gMeanPSFMag, m.gMeanPSFMagErr, m.rMeanPSFMag, m.rMeanPSFMagErr,
     m.iMeanPSFMag, m.iMeanPSFMagErr, m.zMeanPSFMag, m.zMeanPSFMagErr 
into mydb.sro_311 
from ObjectThin o
inner join MeanObject m on o.objID=m.objID
where
    o.nDetections>1
    and o.raMean between 102.34 and 105.65
    and o.decMean between -1.74 and 1.55
    and m.gMeanPSFMag between 12.0 and 16.0
    and m.rMeanPSFMag between 12.0 and 16.0
    and m.iMeanPSFMag between 12.0 and 16.0
    and m.zMeanPSFMag between 12.0 and 16.0
