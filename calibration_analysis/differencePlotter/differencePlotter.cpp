#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "TString.h"
#include "TImage.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TH2.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TMultiGraph.h"
#include "TLegend.h"


#include "differencePlotter.h"


using namespace std;

DiffPlot::DiffPlot(int noI, int noH){
    nImages = noI;
    nHoles = noH;
    holesX.resize(nImages);
    holesY.resize(nImages);
    order.resize(nImages);
    meanHoleX.resize(nHoles);
    meanHoleY.resize(nHoles);
    
    for(int i=0; i<nImages; i++){
        holesX[i].resize(nHoles);
        holesY[i].resize(nHoles);
        order[i].resize(nHoles);
        for(int j=0; j<nHoles; j++){
            holesX[i][j] = -999;
            holesY[i][j] = -999;
        }
    }
}


DiffPlot::~DiffPlot(){
}

void DiffPlot::ingest(std::vector<std::vector<double>> &X, std::vector<std::vector<double>> &Y){
    
    for(int i=0; i<nImages; i++){
        for(int j=0; j<X[i].size(); j++){
           holesX[i][j] = X[i][j];
           holesY[i][j] = Y[i][j];
        }
    }
}


void DiffPlot::IDHoles(){
    for(int i=0; i<nImages; i++){
        ID(i);
    }
}

void DiffPlot::ID(int picID){
    /*
    * compares the positions of the centres in all images
    * and identifies the same hole as being the closest hole to
    * another. The first image is used as the reference for all 
    * other images hence the recovered order is only relative to
    * the first image not an absolute "row 1 column 2" style reference
    */
    for(int i=0; i<nHoles; i++){
        double minDis = 99999999;
        for(int j=0; j<nHoles; j++){
            double dx = holesX[0][i]-holesX[picID][j];
            double dy = holesY[0][i]-holesY[picID][j];
            double dist = dx*dx + dy*dy;
            if(dist<minDis){
                order[picID][i] = j;
                minDis = dist;
            }
        }
        if(minDis > 100){
            std::cout<<"distance exceeded\t"<<minDis<<std::endl; 
        }
    }

    
}

void DiffPlot::pltDevFromAvg(){
    
    calcAvgPos();

    TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
    
    TPad *centre_pad = new TPad("Center_pad", "center_pad", 0.2, 0.2, 1.0, 1.0);
    centre_pad->Draw();
    
    TPad * left_pad = new TPad("left_pad", "left_pad", 0.0, 0.2, 0.2, 1.0);
    left_pad->Draw();
    
    TPad * bottom_pad = new TPad("bottom_pad", "bottom_pad", 0.2, 0.0, 1.0, 0.2);
    bottom_pad-> Draw();
    
    TH2D *h2 = new TH2D("h2", "", 23, -5.25, 5.25, 23, -5.25, 5.25);
    
    double dx, dy;
    
    for(int i=0; i<nImages; i++){
        for(int j=0; j<nHoles; j++){
            dx = holesX[i][order[i][j]] - meanHoleX[order[0][j]];
            dy = holesY[i][order[i][j]] - meanHoleY[order[0][j]];
            h2->Fill(dx, dy);
        }
    }
    TH1D *projh2X = h2->ProjectionX();
    TH1D *projh2Y = h2->ProjectionY();
    
        //central 2d histogram
    
    centre_pad->cd();  
    auto gStyle = new TStyle("Default", "Default Style");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(56);
    h2->SetTitle("Deviation from average spot centre");

    h2->GetXaxis()->SetTitle("X Difference /Pixels");
    h2->GetYaxis()->SetTitle("Y Difference /Pixels");
    gStyle->SetOptStat(11);
//    h2->SetMaximum(10);
    h2->Draw("colz");
    centre_pad->Update();
    TPaveStats *st = (TPaveStats*)h2->FindObject("stats");
    st->SetX1NDC(0.7);
    st->SetX2NDC(0.9);
    st->SetY1NDC(0.7);
    st->SetY2NDC(0.9);
    
    delete gStyle;
    //bottom projected histogram
    
    bottom_pad->cd();
    projh2X->GetYaxis()->SetLabelSize(0.1);
    projh2X->GetXaxis()->SetLabelSize(0.1);
    //  projh2X->GetXaxis()->SetTitle("X Difference /Pixels");
    projh2X->SetStats(false);
    projh2X->Draw("bar");
    
    //left projected histogram
    
    left_pad->cd();
    projh2Y->GetYaxis()->SetLabelSize(0.1);
    projh2Y->GetXaxis()->SetLabelSize(0.1);
    //  projh2Y->GetYaxis()->SetTitle("Y Difference /Pixels");
    
    projh2Y->SetStats(false);
    projh2Y->Draw("hbar");
    
    c1->Print("diffFromAvg.png"); 
    
    delete h2;
    delete c1;

}

void DiffPlot::pltDevRelative(){

    TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
    
    TPad *centre_pad = new TPad("Center_pad", "center_pad", 0.2, 0.2, 1.0, 1.0);
    centre_pad->Draw();
    
    TPad * left_pad = new TPad("left_pad", "left_pad", 0.0, 0.2, 0.2, 1.0);
    left_pad->Draw();
    
    TPad * bottom_pad = new TPad("bottom_pad", "bottom_pad", 0.2, 0.0, 1.0, 0.2);
    bottom_pad-> Draw();
    
    TH2D *h2 = new TH2D("", "", 41, -10.25, 10.25, 41, -10.25, 10.25);
    
    double dx, dy;
    
    for(int i=0; i<nImages; i++){
        for(int k=0; k<i; k++){
            for(int j=0; j<nHoles; j++){
                dx = holesX[i][order[i][j]] - holesX[k][order[k][j]];
                dy = holesY[i][order[i][j]] - holesY[k][order[k][j]];
                h2->Fill(dx, dy);

            }
        }
    }
    TH1D *projh2X = h2->ProjectionX();
    TH1D *projh2Y = h2->ProjectionY();
    
        //central 2d histogram
    
    centre_pad->cd();  
    auto gStyle = new TStyle("Default", "Default Style");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(56);
    h2->SetTitle("Spot Position Change with Splitter");

    h2->GetXaxis()->SetTitle("X Difference /Pixels");
    h2->GetYaxis()->SetTitle("Y Difference /Pixels");
    gStyle->SetOptStat(11);
//    h2->SetMaximum(10);
    h2->Draw("colz");
    centre_pad->Update();
    TPaveStats *st = (TPaveStats*)h2->FindObject("stats");
    st->SetX1NDC(0.7);
    st->SetX2NDC(0.9);
    st->SetY1NDC(0.7);
    st->SetY2NDC(0.9);
    st->SetOptStat(1000001110);
    st->Draw();
    delete gStyle;
    //bottom projected histogram
    
    bottom_pad->cd();
    projh2X->GetYaxis()->SetLabelSize(0.1);
    projh2X->GetXaxis()->SetLabelSize(0.1);
    //  projh2X->GetXaxis()->SetTitle("X Difference /Pixels");
    projh2X->SetStats(false);
    projh2X->Draw("bar");
    
    //left projected histogram
    
    left_pad->cd();
    projh2Y->GetYaxis()->SetLabelSize(0.1);
    projh2Y->GetXaxis()->SetLabelSize(0.1);
    //  projh2Y->GetYaxis()->SetTitle("Y Difference /Pixels");
    
    projh2Y->SetStats(false);
    projh2Y->Draw("hbar");
    
    c1->Print("diffFromOther.pdf"); 
    
    delete h2;
    delete c1;    
    
}


void DiffPlot::pltAbsDev(){
    
    TCanvas *c1 = new TCanvas("c1", "c1", 900, 300);
    TH1D *h1 = new TH1D("h1", "Abs Dev", 20, 0.0, 10.0);
    
    
    
    double dx, dy;
    
    for(int i=0; i<nImages; i++){
        for(int k=0; k<i; k++){
            for(int j=0; j<nHoles; j++){
                dx = holesX[i][order[i][j]] - holesX[k][order[k][j]];
                dy = holesY[i][order[i][j]] - holesY[k][order[k][j]];
                h1->Fill(sqrt(dx*dx + dy*dy));
            }
        }
    }
    
    h1->Draw("bar");
    c1->Print("absDiff.png");
    
    delete h1;
    delete c1;
    
}

void DiffPlot::pltAbsPos(){
    TCanvas *c1 = new TCanvas("c1", "c1", 5000, 5000);
    TH2D *h2 = new TH2D("h2", "Abs Dev", 484, 0, 484, 704, 0, 704);

    
    for(int i=0; i<nImages; i++){
        for(int j=0; j<nHoles; j++){
            h2->Fill(holesX[i][order[i][j]], holesY[i][order[i][j]]);
        }
    }
    
    h2->Draw("colz");
    c1->Print("absPos.png");
    
    delete h2;
    delete c1;
    
}

void DiffPlot::calcAvgPos(){
    for(int i=0; i<nHoles; i++){
        meanHoleX[i] = 0;
        meanHoleY[i] = 0;
        for(int j=0; j<nImages; j++){
            int k = order[j][i];
            meanHoleX[i] += holesX[j][k];
            meanHoleY[i] += holesY[j][k];
        }
        meanHoleX[i] /= nImages;
       meanHoleY[i] /= nImages;
    }
}

void DiffPlot::pltCentres(){
    double x1[nHoles-1];
    double x2[nHoles];
    double y1[nHoles-1];
    double y2[nHoles];

    TCanvas *c1 = new TCanvas("c1", "c1", 2000, 2000);

    for(int i=0; i<nHoles-1; i++){
        
        x1[i] = holesX[0][i];
        y1[i] = holesY[0][i];

    
    }
    for(int i=0; i<nHoles; i++){
        x2[i] = holesX[1][i];
        y2[i] = holesY[1][i];
    }
    
    TGraph *g1 = new TGraph(nHoles-1, x1, y1);
    g1->SetMarkerStyle(2);
    g1->SetMarkerSize(11);

    g1->Draw("ap");
    g1->SetTitle("Reconstructed Spot Centre Positions");
    g1->GetXaxis()->SetTitle("X Pixel");
    g1->GetYaxis()->SetTitle("Y Pixel");

    TGraph *g2 = new TGraph(nHoles, x2, y2);
    g2->SetMarkerStyle(5);
    g2->SetMarkerSize(11);
    g2->SetMarkerColor(2);
//     g2->SetLineWidth(10);
    g2->Draw("p same");
    TLegend *legend = new TLegend(0.7, 0.75, 0.95, 0.9);
    legend->AddEntry(g1, "Current DAQ", "p");
    legend->AddEntry(g2, "Current DAQ with splitter", "p");
     
    legend->Draw();

    c1->Print("CentreDiff.pdf");
    

    delete g1;
    delete g2;
    delete c1;
  //  delete legend;
}
