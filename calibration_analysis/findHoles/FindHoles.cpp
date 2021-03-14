#include <iostream>
#include <stdlib.h>


#include "FindHoles.h"

#define cwdth 500
#define chght 300

using namespace std;

///////////////////////////////////////////////////////////////////////
//                      CTOR                                         //
///////////////////////////////////////////////////////////////////////  
FindHoles::FindHoles(int width, int height)
{
	fWidth = width;
        fHeight = height;
        nHoles = 0;
        array1 = new double[fWidth*fHeight];
        array2 = new double[fWidth*fHeight];
        array3 = new double[fWidth*fHeight];
        array4 = new double[fWidth*fHeight];
        gaussBlur = new double[fWidth*fHeight];
        save_thresh = new double[fWidth*fHeight];
        visited = new int[fWidth*fHeight];
        holesx = NULL;
        holesy = NULL;
        blobSize = NULL;
        nHoleCut = 0;
        goodHoles = 0;
        outline = new double[fWidth*fHeight];
}

///////////////////////////////////////////////////////////////////////
//                      DTOR                                         //
///////////////////////////////////////////////////////////////////////  
FindHoles::~FindHoles() {
  delete [] array1;
  delete [] array2;
  delete [] array3;
  delete [] array4;
  delete [] gaussBlur;
  delete [] save_thresh;
  delete [] visited;
  delete [] outline;
}

void FindHoles::loadArray(double *array, int arr_numb){
  if(arr_numb<0 || arr_numb>3) cout << "Number must be 0-3" << endl;
  for(int i=0; i<fWidth; i++){
    for(int j=0; j<fHeight; j++){
      if(arr_numb==0) array1[i+j*fWidth] = array[i+j*fWidth];
      if(arr_numb==1) array2[i+j*fWidth] = array[i+j*fWidth];
      if(arr_numb==2) array3[i+j*fWidth] = array[i+j*fWidth];
      if(arr_numb==3) array4[i+j*fWidth] = array[i+j*fWidth];
      outline[i+j*fWidth] = 0;
    }
  }
}

void FindHoles::getArray(double *array, int arr_numb){
  if(arr_numb<0 || arr_numb>3) cout << "Number must be 0-3" << endl;
  for(int i=0; i<fWidth; i++){
    for(int j=0; j<fHeight; j++){
      if(arr_numb==0) array[i+j*fWidth] = array1[i+j*fWidth] ;
      if(arr_numb==1) array[i+j*fWidth] = array2[i+j*fWidth];
      if(arr_numb==2) array[i+j*fWidth] = array3[i+j*fWidth];
      if(arr_numb==3) array[i+j*fWidth] = array4[i+j*fWidth];
    }
  }
}

void FindHoles::getManipArray(double *array){
  for(int i=0; i<fWidth*fHeight; i++){
    array[i] = manip_array[i];
    
  }
//   for(int i=0; i<fWidth; i++){
//     for(int j=0; j<fHeight; j++){
//       array[i+j*fWidth] = manip_array[i+j*fWidth];
//     }
//   }
}

void FindHoles::getOutline(double *array){
  for(int i=0; i<fWidth*fHeight; i++){
    array[i] = outline[i];
  }
  
}


void FindHoles::setManipArray(int arr_numb){
  if(arr_numb<0 || arr_numb>3) cout << "Number must be 0-3" << endl;
  if(arr_numb==0) manip_array = array1;
  if(arr_numb==1) manip_array = array2;
  if(arr_numb==2) manip_array = array3;
  if(arr_numb==3) manip_array = array4;
}


void FindHoles::avgFilter(int win_rad, int filter){
  cout << "Applying average filter" << endl;
  int gauss_win_rad = 15;
  double *tmp_array = new double[fWidth*fHeight];
  double gaussKernel[4*(gauss_win_rad+1)*(gauss_win_rad+1)];
  double sigma = 10;

  for(int i=0; i<=2*gauss_win_rad; i++){
      for(int j=0; j<=2*gauss_win_rad; j++){
          double x2 = pow(i-gauss_win_rad, 2) + pow(j-gauss_win_rad, 2);
          gaussKernel[j+i*(2*gauss_win_rad + 1)] = exp(-(x2)/(2*sigma*sigma));
      }
  }
  for(int i=0; i<fWidth; i++){
    for(int j=0; j<fHeight; j++){ 

      double gaussAvg = 0;
      double weight = 0;
      int kernelPos = 0;
      for(int l=-1*gauss_win_rad; l<=gauss_win_rad; l++){
        for(int k=-1*gauss_win_rad; k<=gauss_win_rad; k++){
          int bin_x = i+l;
          int bin_y = j+k;
          if(bin_x>=0 && bin_x<fWidth && bin_y>=0 && bin_y<fHeight){
            int iter = bin_x+bin_y*fWidth;
            kernelPos = k+gauss_win_rad + ((l + gauss_win_rad)*(gauss_win_rad*2 + 1));
            gaussAvg += gaussKernel[kernelPos]*manip_array[iter];
            if(gaussKernel[kernelPos] > 1 || manip_array < 0){
                cout<<"error: "<<gaussKernel[kernelPos]<<endl;
                cout<<kernelPos<<endl;
            }
            
            weight += gaussKernel[kernelPos];
          }
        }
      }
      gaussBlur[i + j*fWidth] = gaussAvg/weight;

    }
  }

  for(int i=0; i<fWidth; i++){
    for(int j=0; j<fHeight; j++){ 
      double avg = 0.;
      int bins = 0;
      for(int l=-1*win_rad; l<=win_rad; l++){
        for(int k=-1*win_rad; k<=win_rad; k++){
          int bin_x = i+l;
          int bin_y = j+k;
          if(bin_x>=0 && bin_x<fWidth && bin_y>=0 && bin_y<fHeight){
            int iter = bin_x+bin_y*fWidth;
            avg = avg+manip_array[iter];
            bins++;
          }
        }
      }
      avg = avg/(double)bins; 
      if(avg<filter){avg = 0.1;}
      tmp_array[i+j*fWidth] = avg;
    }   
 } 
 for(int i=0; i<fWidth*fHeight; i++){
   manip_array[i] = tmp_array[i];
 }
 delete [] tmp_array;

}
         
void FindHoles::adThresh(double factor, int nbox)
{
  cout << "Performing adaptive thresholding..." << endl;
  int i,j; //counters

  //partition image into n x n blocks
  //const int n = 10; //HARDCODE 
  const int n = nbox;

  //number of pixels per block row/column
  int blockRows = fWidth/n;
  int blockCols = fHeight/n;

  int br,bc;               //block row/col label (counters)
  double maxIblock[n][n];  //highest intensity for block (n,n)
  double minIblock[n][n];  //lowest intensity for block (n,n)
  double threshold[n][n];  //cutoff intensity for a pixel to be considered part of a hole
  //double x_center[n][n];
  //double y_center[n][n];
  double rangeIblock[n][n];//range of intensities for block (n,n)
  double avgRange=0;       //average range of intensities over each block
  double totalI;           //total intensity/hits
  double maxI = 254; 

  

  /*calculate max and min intensity of each block*/

  //for each block
  for (br=0; br<n; br++)
    for (bc=0; bc<n; bc++) {
      
      totalI = 0;  //initialize total intensity for current block

      //for each pixel in current block
      for (i=br*blockRows; i<(br+1)*blockRows; i++)
      for(j=bc*blockCols; j<(bc+1)*blockCols; j++) { 

	  //initialize comparison point to be at the corner of current block
	  if (i==br*blockRows && j==bc*blockCols) { 
	    maxIblock[br][bc] = manip_array[i+fWidth*j];  
	    minIblock[br][bc] = manip_array[i+fWidth*j];  
	  }

	  //find the max and min intensity of the current block
	  if (manip_array[i+fWidth*j] > maxIblock[br][bc]) 
		maxIblock[br][bc] = manip_array[i+fWidth*j];  
	  if (manip_array[i+fWidth*j] < minIblock[br][bc]) 
		minIblock[br][bc] = manip_array[i+fWidth*j];

	  //sum the intensities over every pixel in the current block
	  totalI += manip_array[i+fWidth*j];
	}

      //impose upper limit (maxI) on maximum intensity of current block
      maxIblock[br][bc] = (maxIblock[br][bc] > maxI ? maxI : maxIblock[br][bc]);
      
      //calculate range of intensities within current block
      rangeIblock[br][bc] = maxIblock[br][bc] - minIblock[br][bc];

      //sum ranges over all blocks (to be averaged later)
      avgRange += rangeIblock[br][bc];

      //define the threshold value for current block
      //(HARDCODE THESE FACTORS AND FUNCTION, look at threshed .txt image)
      threshold[br][bc] = totalI/(factor*blockRows*blockCols);

    }

  //calculate average range of intensities over each block
  avgRange = avgRange/(n*n);


  /*perform adaptive thresholding*/

  //for each block
  /*for (br=0; br<n; br++)
    for (bc=0; bc<n; bc++) {

      for(i=br*blockRows; i<(br+1)*blockRows; i++) 
	for(j=bc*blockCols; j<(bc+1)*blockCols; j++){	  
         	 
          //Take weighted average of surrounding blocks 
	  if(manip_array[i+fWidth*j] < threshold[br][bc]) 
             manip_array[i+fWidth*j]=0.1; 
            

	  
	}
    }*/
   for(int i=0; i<fWidth; i++){
     for(int j=0; j<fHeight; j++){
       double thresh = 0.;
       double weight = 0.;
       for (br=0; br<n; br++){
         for (bc=0; bc<n; bc++){
           //distance to bloc
           double dist = sqrt( pow((double)i-(double)fWidth/(double)n*((double)br+0.5),2)+
                        pow((double)j-(double)fHeight/(double)n*((double)bc+0.5),2) );
           if(dist<0.0001) dist = 0.0001;
           thresh += threshold[br][bc]/dist/dist;
           weight += 1./dist/dist;
         }
       }
       save_thresh[i+fWidth*j] = thresh/weight;
       if(manip_array[i+fWidth*j] <  thresh/weight)
         manip_array[i+fWidth*j]=0.1; 
     }
  }  


}

void FindHoles::countHoles() {
  cout << "Count Holes" << endl;

  // First clear out the visited array. 
  // Mark every non-zero element that it finds by setting the 
  // corresponding element of the array to "label".  Once a square 
  // has been marked as visited, it will stay marked until all the
  // holes have been counted.  This will prevent the same hole from 
  // being counted more than once.

  int label = 1;
  for (int x = 0; x < fWidth; x++)
    for (int y = 0; y < fHeight; y++)
      visited[x+y*fWidth]=0;

  // For each position in the grid, call getHoleSize() to get the
  // size of the hole at that position.  If the size is not zero, 
  // count a blob.  Note that if we come to a position that was part
  // of a previously counted blob, getHoleSize() will return 0 and
  // the hole will not be counted again.
  for (int x = 0; x < fWidth; x++)
    for (int y = 0; y < fHeight; y++) {
      if (getHoleSize(x,y,label, 0) > 0) 
	label++;
    }  

  nHoles = label;
  
  cout << "Number of possible holes: " << nHoles << endl;
}  

int FindHoles::getHoleSize(int x, int y, int label, bool recursive) {

  // Counts the squares in the blob at position (x,y) in the
  // grid.  Squares are only counted if they are filled and
  // unvisited.  If this routine is called for a position that
  // has been visited, the return value will be zero.
  //cout << "test " << x << " " << y << endl;
  if (x < 0 || x >= fWidth || y < 0 || y >= fHeight) {
    // This position is not in the array, so there is
    // no blob at this position.  Return a blob size of zero.
    return 0;
  }
  
  if (manip_array[x+fWidth*y]<1. && recursive == true){
    outline[x+fWidth*y] = 1;
 //   return 0;

  }
    
  if (manip_array[x+fWidth*y]<1. || visited [x+fWidth*y] != 0) return 0;
  // Mark the pixel as visited so that we won't count it again during the
  // following recursive calls.
  // cout << "test label " << label << endl;
  visited[x+fWidth*y] = label;                        

  int size = 1;  // Count the pixel at this position

  //count the pixels that are connected to this pixel horizontally or vertically
  //cout << "test size 1" << endl;
  size += getHoleSize(x-1,y,label, 1);
  //cout << "test size 2" << endl;
  size += getHoleSize(x+1,y,label, 1);
  //cout << "test size 3" << endl;
  size += getHoleSize(x,y-1,label, 1);
  //cout << "test size 4" << endl;
  size += getHoleSize(x,y+1,label, 1);
  return size;
}

void FindHoles::calcCentroids(int MAXBBOXAREA, int MINBBOXAREA, int MINHIT, int MINPIXELS, bool gaussCentre=false) {
  cout << "Calculate Centroids" << endl;
  
  /*const int MAXBBOXAREA=10000; //maximum area of bounding box of a hole
  const int MINBBOXAREA=20; //minimum area of bounding box of a hole
  const int MINHIT=200;  //minimum intensity of a hole
  const int MINPIXELS=20;*/

  //Grid boundaries (measured after viewing image)
  const int xGridLow=0;
  const int xGridHigh=fWidth;
  const int yGridLow=0;
  const int yGridHigh=fHeight;

  //Center Triangle Boudaries (measured after viewing image)
  //for hole cutting in classImage::calcCentroids()
  const int xTriLow=288;
  const int xTriHigh=334; 
  // Computes centroid and bbox for each hole, and creates new array 
  // with centroid positions marked (rounded to nearest pixel).
  
  HoleStruct *allHole = new HoleStruct[nHoles];

  int r, c, nHole;
  int i=0,j;
  int bboxArea,xwidth,ywidth; 
  float rowsum, colsum, temp;
    if(gaussCentre){

        cout<<"gaussian centre calculation"<<endl;
    }
  for (nHole=0; nHole<nHoles; nHole++) {  

    //if (nHole%500==0) cout << "Processing Hole: " << nHole << endl;

    //initialize data for the nHole'th hole
    allHole[nHole].xCoord = 0.;
    allHole[nHole].yCoord = 0.;
    allHole[nHole].nHits = 0;
    allHole[nHole].nPixels = 0; 
    allHole[nHole].low_x = fWidth;  
    allHole[nHole].high_x = -1; 
    allHole[nHole].low_y = fHeight;  
    allHole[nHole].high_y = -1;
    allHole[nHole].maxHoleI = 0;
    allHole[nHole].ignore = false;
    rowsum = 0.; 
    colsum = 0.;
    
    /*if the background intensity is higher than the fiducial hole 
      intensity values, then we subtract from each hole the largest 
      intensity value in the hole, and consider the absolute value of
      each pixel.*/
    
    for (r=1; r<fWidth; r++)
      for (c=1; c<fHeight; c++) {
	if (visited[r+fWidth*c] == nHole+1 && visited[r+fWidth*c] > allHole[nHole].maxHoleI) 
	  allHole[nHole].maxHoleI=visited[r+fWidth*c];
      }
    
    /* scan array and add all (weighted) pixels with label nHole 
       (nHits= Sum (over pixels) of Intensities) */
    /*for (r=1; r<fWidth; r++) {
      for (c=1; c<fHeight; c++) {
	if (visited[r+fWidth*c] == nHole+1) { 
	  allHole[nHole].nPixels++;
	  allHole[nHole].nHits += abs(manip_array[r-1+fWidth*(c-1)]);
          //allHole[nHole].xCoord += manip_array[r-1+fWidth*(c-1)]*(double)(r-1);
          //allHole[nHole].yCoord += manip_array[r-1+fWidth*(c-1)]*(double)(c-1);
	  rowsum += r*abs(manip_array[r-1+fWidth*(c-1)]);  
	  colsum += c*abs(manip_array[r-1+fWidth*(c-1)]);
	  
	  if(r+1 < allHole[nHole].low_x)  allHole[nHole].low_x = r+1;  
	  if(r+1 > allHole[nHole].high_x) allHole[nHole].high_x = r+1;
	  if(c+1 < allHole[nHole].low_y)  allHole[nHole].low_y = c+1;  
	  if(c+1 > allHole[nHole].high_y) allHole[nHole].high_y = c+1;
	  
	  //if((r==1)|| (c==1)|| (r == xRes) || (c == yRes)) allHole[nHole].ignore = true;
	}
      }
    }*/

    double maxADC = 0.;
    double maxXval = 0.;
    double maxYval = 0.;
    double avg_iter = 0.;
    for (r=0; r<fWidth; r++) {
      for (c=0; c<fHeight; c++) {
        if (visited[r+fWidth*c] == nHole+1) { 
            allHole[nHole].nPixels++;
            allHole[nHole].nHits += abs(manip_array[r+fWidth*(c)]-save_thresh[r+fWidth*(c)]);
                //allHole[nHole].xCoord += manip_array[r-1+fWidth*(c-1)]*(double)(r-1);
                //allHole[nHole].yCoord += manip_array[r-1+fWidth*(c-1)]*(double)(c-1);
            rowsum += ((double)r+0.5)*(double)abs(manip_array[r+fWidth*(c)]-save_thresh[r+fWidth*(c)]);  
            colsum += ((double)c+0.5)*(double)abs(manip_array[r+fWidth*(c)]-save_thresh[r+fWidth*(c)]);
            
            if(gaussCentre){
                if(gaussBlur[r+fWidth*(c)]>maxADC){
                    maxADC = (double)gaussBlur[r+fWidth*(c)];
                }
            }
            else{
                if(manip_array[r+fWidth*(c)]>maxADC)
                maxADC = (double)manip_array[r+fWidth*(c)];
            }

            
            if(r-1 < allHole[nHole].low_x && r-1>=0)  allHole[nHole].low_x = r-1;  
            if(r+1 > allHole[nHole].high_x && r+1<fWidth) allHole[nHole].high_x = r+1;
            if(c-1 < allHole[nHole].low_y && c-1>=0)  allHole[nHole].low_y = c-1;  
            if(c+1 > allHole[nHole].high_y && c+1<fHeight) allHole[nHole].high_y = c+1;
            
            //if((r==1)|| (c==1)|| (r == xRes) || (c == yRes)) allHole[nHole].ignore = true;
        }
      }
    }

     for (r=0; r<fWidth; r++) {
      for (c=0; c<fHeight; c++) {
        if (visited[r+fWidth*c] == nHole+1) {
            if(gaussCentre){
                if(gaussBlur[r+fWidth*(c)]>0.99*maxADC){
                    maxXval = (avg_iter*maxXval+(double)r+0.5)/(avg_iter+1.);
                    maxYval = (avg_iter*maxYval+(double)c+0.5)/(avg_iter+1.);
                    avg_iter +=1.;

                    
                }
            }
            else{
                if(manip_array[r+fWidth*(c)]>0.99*maxADC){
                    maxXval = (avg_iter*maxXval+(double)r+0.5)/(avg_iter+1.);
                    maxYval = (avg_iter*maxYval+(double)c+0.5)/(avg_iter+1.);
                    avg_iter +=1.;
                }
                
            }
        }
      }
    }   
            

    //allHole[nHole].xCoord /= (double)allHole[nHole].nHits;
    //allHole[nHole].yCoord /= (double)allHole[nHole].nHits;

    /*calculate centroid*/
    /*if(allHole[nHole].nHits) allHole[nHole].xCoord = rowsum / allHole[nHole].nHits;
    else {allHole[nHole].xCoord = 0;}
    if(allHole[nHole].nHits) allHole[nHole].yCoord = colsum / allHole[nHole].nHits;
    else {allHole[nHole].yCoord = 0;}*/
    if(allHole[nHole].nHits) allHole[nHole].xCoord = maxXval;
    else {allHole[nHole].xCoord = 0;}
    if(allHole[nHole].nHits) allHole[nHole].yCoord = maxYval;
    else {allHole[nHole].yCoord = 0;}

    //boundary of grid based on actual image (../include/holePars.h)
    if ( allHole[nHole].xCoord < xGridLow || allHole[nHole].xCoord > xGridHigh || 
	 allHole[nHole].yCoord < yGridLow || allHole[nHole].yCoord > yGridHigh )
      allHole[nHole].ignore = true;      
    
    /*If the area of the bounding box is larger than MAXBOXAREA, the hole 
      is too large to be a fiducial hole, and we ignore this hole.  
      Also ignore hole if nHits is less than MINHIT (this ignores the sporadic 
      background "fake" holes)
    */
    allHole[nHole].bboxArea = (allHole[nHole].high_x - allHole[nHole].low_x)*(allHole[nHole].high_y - allHole[nHole].low_y); 
    
    if ( (allHole[nHole].bboxArea > MAXBBOXAREA || allHole[nHole].nHits < MINHIT 
	  || allHole[nHole].nPixels < MINPIXELS) && allHole[nHole].ignore==false ) {
      
      allHole[nHole].ignore = true;
            
      //cout << "Par. Ignored "; printHole(allHole,nHole);      
    }

    /*
    //ignore holes that are too close together (except center triangle)
    if ( allHole[nHole].ignore==false ) {
      for (j=nHole-1; j>=0; j--) {
	if ( Hole[j].ignore==false &&
	     abs(allHole[nHole].xCoord - Hole[j].xCoord) < xGridSpace/2 && //HARDCODE factor
	     abs(allHole[nHole].yCoord - Hole[j].yCoord) < yGridSpace/1.817197531 &&
	     (allHole[nHole].xCoord < xTriLow || allHole[nHole].xCoord > xTriHigh) ) {
	  allHole[nHole].ignore=true;
	  Hole[j].ignore=true; 
	}
      }
    }
    */
    
  } /* end of for nHole loop */
  
  //count number of actual holes 
  for (nHole=0; nHole<nHoles; nHole++) {  
    if (allHole[nHole].ignore==false) nHoleCut++;
  }
  
  //allocate array to store actual grid holes
  HoleStruct *cutHole = new HoleStruct[nHoleCut];
  //fill array with actual holes
  if(holesx) 
    delete [] holesx;
  holesx = new double[nHoles];
  if(holesy) 
    delete [] holesy;
  holesy = new double[nHoles];
  if(blobSize) 
    delete [] blobSize;
  blobSize = new int[nHoles];
  
  goodHoles = 0;
  for (nHole=0; nHole<nHoles; nHole++)
    if (allHole[nHole].ignore==false) {
      

      cutHole[i] = allHole[nHole];
      if(cutHole[i].xCoord<30 || cutHole[i].xCoord>454) continue;
      if(cutHole[i].yCoord<30 || cutHole[i].yCoord>674) continue;
      
      //cout << "Real " << cutHole[i].xCoord << " " << cutHole[i].yCoord << endl;
      holesx[i] = cutHole[i].xCoord;
      holesy[i] = cutHole[i].yCoord;
      blobSize[i] = cutHole[i].nPixels;
      i++;
      goodHoles++;
    }
  delete [] allHole;
}

void FindHoles::getHoles(double *xholes, double *yholes, int *size){
  for(int i=0; i<goodHoles; i++){
    xholes[i] = holesx[i];
    yholes[i] = holesy[i];
    size[i] = blobSize[i];
  }
}




