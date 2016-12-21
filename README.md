# pskyline

## Description
A program which implements a spatial p-skyline, i.e. one that receives as INPUT: 
a) a set S of points with imprecise coordinates, 
b) a threshold probability p which is the minimum non dominance probability 
of every point which will be in the final solution;
and OUTPUTs the p-skyline subset R from S, i.e. the set of all points in S
with a Pareto efficience probability at least equal to p

## Arguments
### mandatory
tab: a data.frame with each line representing a candidate to be evaluated by pSkyline
latCol: column number for the latitudes of the spatial
lonCol: column number for the longitudes of the spatial
meanError: the mean error in METERS in the coordinates
sdError: the standard error in METERS in the coordinates

### recommended but not mandatory 
p:      the threshold probability required for a spatial data object to be in the 
        skyline. By default 0.5
refLat, refLon: arrays of reference coordinates of points from which the distance of 
         each dataset spatial object must be computed. If NULL is passed 
         (default), the coordinates of the spatial object in the FIRST ROW are 
         used as reference coordinates.
degreeNotation: TRUE if the coordinates are in degree notation (default). 
         FALSE if they are in radian
pdfSpatial: pdf option for error in the coordinates (default = "normal")

## Details
The algorithm proceeds by deriving new columns to each reference point. The 
distances of each spatial object i is computed to each reference point and these 
values are stored in the derived columns. THOSE derived values are used after to 
compute the spatial p-skyline under the criterium of minimization. Thus the goal is 
to find those spatial points which are the nearest ones with respect to the 
references points.
----------
The pdf options are: "normal" (default), "exponencial" and "chisquare". These 
names also may be provides as "norm", "expo", "exp" and "chisq".