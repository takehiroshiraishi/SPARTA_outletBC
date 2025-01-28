# SPARTA_outletBC

## outletv input parameters
## outletvgrad input parameters
```
int directionv;         // velocity direction x:0 y:1 z:2
int directiongrad;      // direction in which velocity has a gradient x:0 y:1 z:2
double grad;            // gradient
double refvel;          // reference velocity at the zeropoint
double width;           // width of region where gradient exists
double maxwidth;        // width of the system

double margin = (maxwidth - width) * 0.5; // width of the region without gradient
```
