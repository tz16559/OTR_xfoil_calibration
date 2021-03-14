void convert_to_root(TString filename);



void smooth(int n, double *y){
  double out[n];
  out[0] = 0.5*(y[0]+y[1]);
  out[n-1] = 0.5*(y[n-2]+y[n-1]);
  double kernel[3] = {0.25, 0.5, 0.25};
  for(int i=1; i<n-1; i++){
    out[i] = 0;
    for(int j=0; j<3; j++){
      out[i] += kernel[j]*y[i-1+j];
    }
  }
  for(int i=0; i<n; i++){
    y[i] = out[i];
  }
}

void convert_to_root(){
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetOptTitle(1);

//  TString dirname("johanPGM");
  TString dirname("calib_2021_02_25_pgm/filament1_calib/");
  //TString dirname("OTR2_20180531");
  TFile * f = new TFile("calib.root","recreate");
  TSystemDirectory dir(dirname,dirname);
  TList *files = dir.GetListOfFiles();
  if(files){
    TSystemFile * file;
    TString fname;
    TIter next(files);
    while((file=(TSystemFile*)next())){
      fname = file->GetName();
      cout<<fname<<endl;
      if(!file->IsDirectory() && fname.EndsWith("pgm")){
        cout << fname.Data() << endl;
        convert_to_root(dirname+"/"+fname);
      }
    }
  }
  f->Write();
  f->Close();
}

void flatten_x(TString filename, bool bSmooth=false){
  cout << "input: " << filename << endl;

  TImage* raw = TImage::Open(filename.Data());
  UInt_t* array = raw->GetArgbArray();
  int w = raw->GetWidth();
  int h = raw->GetHeight();
  cout << "width x height " << w << " " << h << endl;

  int sPos[2] = {113, 225};
  int size[2] = {80, 8};
  int ePos[2] = {sPos[0]+size[0], sPos[1]+size[1]};
  int counter=0;
  double xgraph[size[1]];
  double bars[size[1]];
  
  for(int y=sPos[1]; y<ePos[1]; y++){
    double rowSum=0;
    for(int x=sPos[0]; x<ePos[0]; x++){
      int index = w*y +x;
      rowSum += ( array[index] & 0xFF0000 ) >> 16;
    }
    
    bars[counter] = rowSum;
    xgraph[counter] = y;
    counter++;
    
  }
  if(bSmooth){
    smooth(size[1], bars);
  }
  TString outfile(filename);
  outfile.ReplaceAll(".pgm",".png");

  TCanvas * c1d = new TCanvas("c1d","c1d",700,800);
  TGraph *gr1 = new TGraph(size[1], xgraph, bars);
  gr1->SetFillColor(40);

  //c1d->GetXaxis()->SetTitle("Pixel Y");
  gr1->Draw("AB");
  outfile.ReplaceAll(".png","_flatx.png");
  c1d->Print(outfile.Data());
  delete c1d;
  
  
}

void flatten_y(TString filename, bool bSmooth=false){
  cout << "input: " << filename << endl;

  TImage* raw = TImage::Open(filename.Data());
  UInt_t* array = raw->GetArgbArray();
  int w = raw->GetWidth();
  int h = raw->GetHeight();
  cout << "width x height " << w << " " << h << endl;

  int sPos[2] = {123, 225};
  int size[2] = {63, 5};
  int ePos[2] = {sPos[0]+size[0], sPos[1]+size[1]};
  int counter=0;
  double xgraph[size[0]];
  double bars[size[0]];
  
  for(int x=sPos[0]; x<ePos[0]; x++){
    double colSum=0;
    for(int y=sPos[1]; y<ePos[1]; y++){
      int index = w*y +x;
      colSum += ( array[index] & 0xFF0000 ) >> 16;
    }
    
    bars[counter] = colSum;
    xgraph[counter] = x;
    counter++;
    
  }
  if(bSmooth){
    smooth(size[0], bars);
  }
  TString outfile(filename);
  outfile.ReplaceAll(".pgm",".png");

  TCanvas * c1d = new TCanvas("c1d","c1d",700,800);
  TGraph *gr1 = new TGraph(size[0], xgraph, bars);
  gr1->SetFillColor(40);

  //c1d->GetXaxis()->SetTitle("Pixel Y");
  gr1->Draw("AB");
  outfile.ReplaceAll(".png","_flaty.png");
  c1d->Print(outfile.Data());
  delete c1d;

}

void overlay_light_curves(TString filename1, TString filename2, TString filename3){
    
    TImage* raw1 = TImage::Open(filename1.Data());
    UInt_t* array1 = raw1->GetArgbArray();
    int w = raw1->GetWidth();
    int h = raw1->GetHeight();
    
    TImage* raw2 = TImage::Open(filename2.Data());
    UInt_t* array2 = raw2->GetArgbArray();
    
    TImage* raw3 = TImage::Open(filename3.Data());
    UInt_t* array3 = raw3->GetArgbArray();
    
    int xrange[2] = {10, 110};
    int yrange[2] = {10, 110}; 

    double_t background1 = 0;
    double_t background2 = 0;
    double_t background3 = 0;
    
    double_t counter = 0;
    for(int i=0; i<w; i++){
        for(int j=0; j<h; j++){
            int index = w*j+i;
            if(i>=xrange[0] && i<xrange[1] && j>=yrange[0] && j<yrange[1]){
            counter+=1;
            background1 += ( array1[index] & 0xFF0000 ) >> 16;
            background2 += ( array2[index] & 0xFF0000 ) >> 16;
            background3 += ( array3[index] & 0xFF0000 ) >> 16;
            
            }
        }
    }
    background1 /= counter;
    background2 /= counter;
    background3 /= counter;
    
    
    TH1D * histo1 = new TH1D("Old System", "OTR Image Light Curves", 256, 0, 256);
    TH1D * histo2 = new TH1D("Parallel Running", "Parallel Running", 256, 0, 256);
    TH1D * histo3 = new TH1D("New System", "New System", 256, 0, 256);
    
    
    background1=0;
    background2=0;
    background3=0;
    
    
    for(int i=0; i<w; i++){
        for(int j=0; j<h; j++){
            int index = w*j+i;
            double_t px1 = ( array1[index] & 0xFF0000 ) >> 16;
            double_t px2 = ( array2[index] & 0xFF0000 ) >> 16;
            double_t px3 = ( array3[index] & 0xFF0000 ) >> 16;
            histo1->Fill(px1 - background1);
            histo2->Fill(px2 - background2);
            histo3->Fill(px3 - background3);
            
        }
    }
    
    
    TCanvas * c1d = new TCanvas("c1d","c1d",1000,1000);
    histo1->GetXaxis()->SetTitle("ADC");
    histo1->GetYaxis()->SetTitle("Number of Pixels");
    //c1d->SetLogy(1);
    gStyle->SetOptStat(0);
    c1d->SetLeftMargin(0.15);
    c1d->SetRightMargin(0.05);
    histo1->GetYaxis()->SetRangeUser(0,25000);
    histo1->GetXaxis()->SetRangeUser(0,100);
    //histo1->GetYaxis()->SetRange(0, 13000);
    //histo1d->SetMaximum(100000);
    histo1->SetLineStyle(3);
    histo1->SetLineWidth(3);
    histo1->Draw();
    
    histo2->SetLineColor(kRed);
    histo2->SetLineWidth(2);
    histo2->Draw("same");
    
    histo3->SetLineColor(kBlack);
    histo3->SetLineWidth(2);
  //  histo3->SetLineStyle(3);
    histo3->Draw("same");
    
    TLegend *legend = new TLegend(0.6, 0.6, 0.95, 0.9);
    legend->AddEntry(histo1, "Old System", "l");
    legend->AddEntry(histo2, "Parallel Running", "l");
    legend->AddEntry(histo3, "New System", "l");
     
    legend->Draw();
    
    c1d->Print("Old_New_light_curves.png");
    delete c1d, histo1, histo2, histo3, legend;
    
    
}

void convert_to_root(TString filename){

  cout << "input: " << filename << endl;

  TImage* raw = TImage::Open(filename.Data());
  UInt_t* array = raw->GetArgbArray();
  int w = raw->GetWidth();
  int h = raw->GetHeight();
  int xrange[2] = {10, 110};
  int yrange[2] = {10, 110};
  double offset = 0;//19.391;
  cout << "width x height " << w << " " << h << endl;

  TString histoname(filename);
  histoname.ReplaceAll(".pgm","");
  histoname.Remove(0,histoname.Last('/')+1);
  //TH2D * histo = new TH2D(histoname,histoname,w,0,w,h,0,h);
  TH2D * histo = new TH2D(histoname,histoname,h,0,h,w,0,w);
  TH1D * histo1d = new TH1D(histoname+"_1d","Current DAQ",255,0,255);
  
  double_t background = 0;
  double_t counter = 0;
  for(int i=0; i<w; i++){
    for(int j=0; j<h; j++){
      int index = w*j+i;
      if(i>=xrange[0] && i<xrange[1] && j>=yrange[0] && j<yrange[1]){
        counter+=1;
        background += ( array[index] & 0xFF0000 ) >> 16;
      }
    }
  }
  background /= counter;
  
  background=0;
  
  for(int i=0; i<w; i++){
    for(int j=0; j<h; j++){
      int index = w*j+i;
      //array contains ARGB 8 bits each 
      //shift it by 16 to get R
      double_t content = ( array[index] & 0xFF0000 ) >> 16;
      //if( content < threshold ) continue;
      //histo->SetBinContent(i+1,j+1,content);
      histo->SetBinContent(h-j,w-i,content-background+offset);
      histo1d->Fill(content-background+offset);
    }
  }

  TString outfile(filename);
  outfile.ReplaceAll(".pgm",".png");

  TCanvas * c = new TCanvas("c","c",1000,1000);
  histo->GetZaxis()->SetRangeUser(-25,150);
  histo->GetXaxis()->SetNdivisions(505);
  histo->GetXaxis()->SetTitle("pixel X");
  histo->GetYaxis()->SetTitle("pixel Y");
  gStyle->SetNumberContours(255);
  gStyle->SetPalette(53);
  gStyle->SetOptStat(0);
  histo->Draw("colz");
  c->Print(outfile.Data());
  delete c;

  TCanvas * c1d = new TCanvas("c1d","c1d",700,700);
  histo1d->GetXaxis()->SetTitle("ADC");
  histo1d->GetYaxis()->SetTitle("Number of Pixels");
  c1d->SetLogy(1);
  
  
  //begin presentation style (Large stats box)
  
//   TPaveStats *ptstats = new TPaveStats(0.55,0.7,0.85,0.9,"brNDC");
//   ptstats->SetName("stats");
//   ptstats->SetBorderSize(1);
//   ptstats->SetFillColor(0);
//   ptstats->SetTextAlign(12);
//   ptstats->SetTextFont(42);
// 
//   char buffer1[20] = "";
//   char buffer2[20] = "Mean     ";
//   snprintf(buffer1, sizeof buffer1, "%.2f", histo1d->GetMean());
// 
//   strcat(buffer2, buffer1);
// //  const char *buffer3 = buffer2;
//   //     strcat(buffer1, "Mean = ");
//   //     strcat(buffer1, buffer2);
//   TText *ptstats_LaTex = ptstats->AddText(buffer2);
//   
//   char buffer3[20] = "";
//   char buffer4[20] = "Std Dev     ";
//   snprintf(buffer3, sizeof buffer3, "%.2f", histo1d->GetRMS());
// 
//   strcat(buffer4, buffer3);
// 
// 
//   ptstats_LaTex = ptstats->AddText(buffer4);
//   ptstats->SetOptStat(1100);
//   ptstats->SetOptFit(0);
//   ptstats->Draw();
//   histo1d->GetListOfFunctions()->Add(ptstats);
//   ptstats->SetParent(histo1d);
//   
  //end presentation style
  
  histo1d->Draw();

//uncomment gStyle->SetOptStat below for report mode (small stats box)
  gStyle->SetOptStat(1100);
  gPad->Update();
  TPaveStats *st = (TPaveStats*)histo1d->FindObject("stats");
  st->SetX1NDC(0.7);
  st->SetX2NDC(0.9);
  
  //histo1d->GetYaxis()->SetRange(1, 100000);
  //histo1d->SetMaximum(100000);

  outfile.ReplaceAll(".png","_1d.pdf");
  c1d->Print(outfile.Data());
  delete c1d;
  
  

  

}
