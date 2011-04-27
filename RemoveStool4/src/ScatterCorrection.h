#ifndef SCATTER_CORRECTION_H
#define SCATTER_CORRECTION_H

#define FEATURES 1
#define FILTER  1
#define SCALE 12
#define HIST_THRESHOLD 0.95
#define tolerance 13.0

typedef unsigned char SSlice[512][512]; 

typedef int Int4;
int  mean_density_value;
int  background = 0;
float *scale_map;

int pcol,prow,pslice;
int *sphere_no_points;
short sphere_points[100][100][3];
double scatter_mask[2*FILTER+1][2*FILTER+1], mask_total;
unsigned short *data_scale16; 

unsigned char *scale_image;

float filterP = 0.07;
int scaleP = 2;

#endif  // SCATTER_CORRECTION_H