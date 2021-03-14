#include "FindHoles.h"
#include "differencePlotter.h"
#include "TString.h"
#include "TImage.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TH2.h"
#include "TStyle.h"
#include "TPaveStats.h"

#include <unistd.h>
#include <iostream>
#include <fstream>

int main(int argc, char** argv ){
  // get image
  
  std::vector<std::vector<int>> blobSizes;
  std::vector<std::vector<double>> holx;
  std::vector<std::vector<double>> holy;
//   std::vector<TString> names;
  blobSizes.resize(argc-1);
  holx.resize(argc-1);
  holy.resize(argc-1); 
  

  
  
  TFile *infile = new TFile("../calib.root", "read");
  for(int picID=1; picID<argc; picID++){
//     names[picID-1] = argv[picID];
      
    
    TH2D *img_histo = (TH2D*)infile->Get(argv[picID]);
    // store image in array
    int w = img_histo->GetNbinsX();
    int h = img_histo->GetNbinsY();

    double *img_array = new double[w*h];
    //removing any sub-zero values from the image and storing this in an array
 
    for(int i=0; i<w; i++){
        for(int j=0; j<h; j++){
            int index = w*j+i;
            double content = img_histo->GetBinContent(i+1,j+1);
            if( content < 0 ) img_array[index] = 0;
            else img_array[index] = content;
        }
    }

    // find holes
    FindHoles *fh = new FindHoles(w,h);
    fh->loadArray(img_array,0);
    fh->setManipArray(0);
    delete[] img_array;

    int avgFilterWin = 2;
    int avgFilterThresh = 17;
    int adThreshBlocks = 10; 
    double adThreshThresh = .45;
    fh->avgFilter(avgFilterWin, avgFilterThresh);
    double filteredImg[w*h];
    fh->getManipArray(filteredImg);
    fh->adThresh(adThreshThresh, adThreshBlocks);

    fh->countHoles();
    double filteredOutline [w*h];
    fh->getOutline(filteredOutline); 
    fh->calcCentroids(10000,200,500,100, false);
    double *Xholes = new double[fh->goodHoles];
    double *Yholes = new double[fh->goodHoles];
    int *size = new int[fh->goodHoles];
    fh->getHoles(Xholes, Yholes, size);
    int nholes = fh->goodHoles;

    // match holes with control points 

    double *Xcontrol = new double[nholes];
    double *Ycontrol = new double[nholes];
    std::ifstream fc;
    std::string controlFileName("input/controlpoints_8.txt");
    fc.open( controlFileName.c_str() );

    for(int i=0; i<nholes; i++){
        fc>>Xcontrol[i];
        fc>>Ycontrol[i];
        Ycontrol[i] = h - Ycontrol[i];
        Xcontrol[i] += 0.;
        Ycontrol[i] += 0.;
    }
    fc.close();

    int *matched = new int[nholes];
    for(int i=0;i<nholes;i++) matched[i]=-99;

    for(int i=0; i<fh->goodHoles; i++){
        //Start with large distance between hole and match
        double min_dist = 99999999.;
        int match_index = 99;
        //Iterate through control points
        for(int j=0; j<nholes; j++){
            //Calculate distance between hole and control point
            const double dx = Xcontrol[j]-Xholes[i];
            const double dy = Ycontrol[j]-Yholes[i];
            const double distance = sqrt( dx*dx + dy*dy );

            //If there is already a match for this hole, check for a better one
            if(matched[j]>=0){
                const double dx_old = Xcontrol[j]-Xholes[matched[j]];
                const double dy_old = Ycontrol[j]-Yholes[matched[j]];
                double distance_old = sqrt( dx_old*dx_old + dy_old*dy_old );
                if(distance>distance_old) continue;
            }
            //If smallest distance, record distance and control point index
            if(distance < min_dist) {
                match_index = j;
                min_dist = distance;
            }
        }
        //Set the matched value for the control point to the found hole
        matched[match_index] = i;
    }

    // store holes
    for(int k=0; k<nholes; k++){
        blobSizes[picID-1].push_back(size[k]);
        holx[picID-1].push_back(Xholes[k]);
        holy[picID-1].push_back(Yholes[k]);
    }
    
    char outfname[strlen(argv[picID])];
    strcpy(outfname, argv[picID]);
    strcat(outfname, "_foundholes.txt");
    std::string outf(outfname);
    std::ofstream fholes( outf.c_str() );
    fholes << "0\t" << img_histo->Integral() << std::endl;
    for(int i=0; i<nholes; i++){
        //fholes << i+1 << "\t";
        fholes << i+1 << "\t" << Xholes[i] << "\t" << Yholes[i] << std::endl;
        //else fholes << -99 << "\t" << -99 << std::endl;
    }
    fholes.close();

    // draw image with found and control points
    TCanvas *c = new TCanvas("image","image",2000,2000);
    img_histo->SetStats(0);
    img_histo->SetTitle(argv[picID]);
    img_histo->GetXaxis()->SetTitle("Pixel X");
    img_histo->GetYaxis()->SetTitle("Pixel Y");
    img_histo->SetTitle("Calibration Image");
    gStyle->SetPalette(57);
    img_histo->SetMinimum(-1);
    img_histo->Draw("colz");
    TGraph *gh = new TGraph(fh->goodHoles,Xholes,Yholes);
    gh->SetMarkerSize(5);
    gh->SetMarkerColor(kRed);
//     gh->Draw("* same");
    TGraph *gc = new TGraph(nholes,Xcontrol,Ycontrol);
//    gc->SetMarkerStyle(kCircle);
//    gc->Draw("p same");
//    c->SetGrid(1,1);
    
    strcpy(outfname, argv[picID]);
    strcat(outfname, "_output.png");
    c->Print(outfname);

    delete c;
    // clean


    TCanvas *canv = new TCanvas("image", "image", 2000, 2000);
    TH2D *h2 = new TH2D("h2", "",  49, 231, 280, 70, 510, 580);
    
    int iter = 0;
    for(int y=0; y<h; y++){
        for(int x=0; x<w; x++){
            
           if(x>231 && x<=280 && y>510 && y<=580){
                h2->SetBinContent(x-231, y-510, filteredImg[y*w+x]);
                iter++;
                
            }
            
            


//             if(filteredOutline[i] < 1){
//                 h2->SetBinContent(i, filteredImg[i]);
//             }
//             else{
//                 h2->SetBinContent(i, filteredImg[i]);
//             // h2->SetBinContent(i, 100);
//             }
        }
    }
    h2->SetTitle("Raw Image Spot");
    h2->SetStats(0);
    gStyle->SetPalette(57);
    gStyle->SetNumberContours(255);
    h2->SetMinimum(-1);
    h2->GetXaxis()->SetTitle("Pixel x");
    h2->GetYaxis()->SetTitle("Pixel y");
    h2->GetYaxis()->SetTitleOffset(1.5);
    h2->Draw("colz");
    canv->SetRightMargin(0.12);
    canv->SetLeftMargin(0.12);
    
//     uncommenting this will cause a seg fault
//     it does however produce the plot before crashing
     TGraph *gh1 = new TGraph(fh->goodHoles,Xholes,Yholes);

     gh1->SetMarkerSize(20);
     gh1->SetMarkerColor(kBlack);
     gh1->SetMarkerStyle(39);
    // gh1->Draw("* same");

    canv->Print("flattenedImg.png");

    delete h2;
    delete canv;
    delete[] matched;
    delete[] Xholes;
    delete[] Yholes;
    delete[] Xcontrol;
    delete[] Ycontrol;
    delete img_histo;
    delete fh;
    delete gh;
    delete gc;
    
  }
  int nHoles = 28;
//   int comp[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  
  int comp[4][2] = {{0, 1}, {0, 2}, {0, 3}, {0, 4}};
 TString title[4] = {"s000 vs s001",
                    "s000 vs s002",
                    "s000 vs s004",
                    "s000 vs s005"};
                    
 TString filename[4] = {"2d_diff_s000_vs_s001.png",
                        "2d_diff_s000_vs_s002.png",
                        "2d_diff_s000_vs_s004.png",
                        "2d_diff_s000_vs_s005.png"};
  
  
//   TString title[6] = {"Parallel Windows vs Non-Parallel Windows",
//                 "Parallel from Windows vs Parallel from Linux",
//                 "Parallel from Windows vs Non-Parallel Linux",
//                 "Non-Parallel Windows vs Parallel from Linux",
//                 "Non-Parallel Windows vs Non-Parallel Linux",
//                 "Parallel from Linux vs Non-Parallel Linux"};
//                 
//   TString filename[6] = {"2d_diff_ParW_vs_W.png",
//                         "2d_diff_ParW_vs_ParL.png",
//                         "2d_diff_ParW_vs_L.png",
//                         "2d_diff_W_vs_ParL.png",
//                         "2d_diff_W_vs_L.png",
//                         "2d_diff_ParL_vs_L.png"};
//                         
  if(argc>2){
    TH2D *img_histo1 = (TH2D*)infile->Get(argv[1]);
    TH2D *img_histo2 = (TH2D*)infile->Get(argv[2]);

    int w = img_histo1->GetNbinsX();
    int h = img_histo1->GetNbinsY();

    double *img_array = new double[w*h];

    for(int i=0; i<w; i++){
        for(int j=0; j<h; j++){
            int index=j*w+i;
            img_array[index] = img_histo1->GetBinContent(i+1, j+1);
            img_array[index] -= img_histo2->GetBinContent(i+1, j+1);
        }
    }


    delete img_histo1;
    delete img_histo2;
    
    
    TCanvas *canv = new TCanvas("image", "image", 2000, 2000);
    TH2D *diff_histo = new TH2D("h2", "",  w, 0, w, h, 0, h);
    
    for(int i=0; i<w-1; i++){
        for(int j=0; j<h-1; j++){
            int index=j*w+i;
            double val = img_array[index] + img_array[index+1] + img_array[index+w] + img_array[index + w +1];
            diff_histo->SetBinContent(i+1, j+1, val);
        }
    }
    
    diff_histo->Draw("colz");
    canv->Print("Difference_Histo.png");
    
    delete canv, diff_histo, img_array;
    
    DiffPlot *pl = new DiffPlot(argc-1, nHoles);
    pl->ingest(holx, holy);
    pl->IDHoles();
    pl->pltDevFromAvg();
    pl->pltDevRelative();
    pl->pltAbsDev();
    pl->pltAbsPos();
    pl->pltCentres();
    delete pl;
      
  }


  return 0;
}
