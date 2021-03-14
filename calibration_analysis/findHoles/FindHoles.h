#ifndef _FINDHOLES_
#define _FINDHOLES_

#include <time.h>
#include <stdio.h>
//#include <math.h>
#include <cmath>

//#define W 700
//#define H 550

#define NP 256

struct HoleStruct {
  double xCoord, yCoord; //pixels
  double low_x, high_x, low_y, high_y;
  float bboxArea;
  int rowLabel, colLabel; 
  int nHits;    //number of photon hits for the hole
  int nPixels;  //number of pixels in the hole
  int maxHoleI;
  bool ignore;
};


class FindHoles  
{
 private:
   double *array1, *array2, *array3, *array4, *manip_array, *gaussBlur;
   int *visited;
   int nHoles;
   int fWidth, fHeight;
 
   double *holesx, *holesy;
   int *blobSize;
   double *save_thresh;

 public:
  FindHoles(int width, int height);
  virtual ~FindHoles();
  double *outline;
  int nHoleCut;
  
 public:

   void avgFilter(int win_rad, int filter);
   void adThresh(double factor, int nbox);
   void countHoles();
   int getHoleSize(int x, int y, int label, bool recursive);
   void calcCentroids(int MAXBBOXAREA, int MINBBOXAREA, int MINHIT, int MINPIXELS, bool gaussCentre);
   void loadArray(double *array, int arr_numb);
   void getArray(double *array, int arr_numb);
   void getManipArray(double *array);
   void setManipArray(int arr_numb);
   void getOutline(double *array);
   void getHoles(double *xholes, double *yholes, int *size);
   int goodHoles;
  
};

#endif
