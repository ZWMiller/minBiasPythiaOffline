{

  TFile* f;
  char fileName[1000];

  TString fileList[16] = {
    "FFOutput/sysChange_etaContShiftMinus_11_12_FIT.root",
    "FFOutput/sysChange_etaContShiftPlus_11_12_FIT.root",
    "FFOutput/sysChange_etaContStatMinus_11_12_FIT.root",
    "FFOutput/sysChange_etaContStatPlus_11_12_FIT.root",
    "FFOutput/sysChange_gammaContStatMinus_11_12_FIT.root",
    "FFOutput/sysChange_gammaContStatPlus_11_12_FIT.root",
    "FFOutput/sysChange_piContStatMinus_11_12_FIT.root",
    "FFOutput/sysChange_piContStatPlus_11_12_FIT.root",
    "FFOutput/sysChange_gammaPiContShiftMinus_11_12_FIT.root",
    "FFOutput/sysChange_gammaPiContShiftPlus_11_12_FIT.root",
    "FFOutput/sysChange_fitNormPlus_11_12_FIT.root",
    "FFOutput/sysChange_fitNormMinus_11_12_FIT.root",
    "FFOutput/sysChange_pileupPlus_11_12_FIT.root",
    "FFOutput/sysChange_pileupMinus_11_12_FIT.root",
    "FFOutput/sysChange_totalRecoPlus_11_12_FIT.root",
    "FFOutput/sysChange_totalRecoMinus_11_12_FIT.root",

  };

  const int numPtBin = 14;
  int n;
  double x[16][numPtBin],y[16][numPtBin];
  double sysP[numPtBin]={0.}, sysM[numPtBin]={0.}, pT[numPtBin]={0.};
  TGraphAsymmErrors* grEr;
  for(int i = 0; i < 16; i++)
  {
    f = new TFile(fileList[i]);
    if(!f->IsOpen())
      continue;

    cout << "Opened " << fileList[i] << endl;

    TGraphErrors* gr = (TGraphErrors*)f->Get("sysChange");
    n = gr->GetN(); // Get plot array dimension
    for(int ii = 0; ii < n; ii++)
    {
      gr->GetPoint(ii,x[i][ii],y[i][ii]); // get point ii (from array), store in x,y
      cout << "pt: " << x[i][ii] << " %: " << y[i][ii] << endl;
      pT[ii] = x[i][ii];
      if(y[i][ii] > 0) // Add up the squares of all positive contrib in each pt
        sysP[ii] += y[i][ii]*y[i][ii];
      if(y[i][ii] < 0) // Add up the squares of all negative contrib in each pt
        sysM[ii] += y[i][ii]*y[i][ii];
    }
  

  }

  for(int ii = 0; ii < n; ii++) // After all runs, loop through pT, to square root
  {
    sysP[ii] = sqrt(sysP[ii]); // sqrt(x*x + y*y + ...)
    sysM[ii] = sqrt(sysM[ii]);
    cout << "pT: " << pT[ii] << " Sys%(+/-): " << sysP[ii] << " / " << sysM[ii] << endl; 
  }
  
  

}

