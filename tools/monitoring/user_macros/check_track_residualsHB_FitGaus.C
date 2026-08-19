//const int NStaves[NLayer] = { 12, 16, 20, 24, 30, 42, 48 };


void fitter(double scale, TH2D* histo, TH1D* h_mean, TF1* fgaus, TString tag = ""){
   //fitter
   int nbinsx = histo->GetNbinsX();
   for(int ibin = 0; ibin < nbinsx; ibin++){
      fgaus->SetParameter(0,0);
      fgaus->SetParameter(1,0);
      fgaus->SetParameter(2,0);
      fgaus->SetParError(0,0);
      fgaus->SetParError(1,0);
      fgaus->SetParError(2,0);
      auto proj_histo_ = (TH1D*) histo->ProjectionY(Form("proj_histo%s_%d",(const char*)tag,ibin),1+ibin,1+ibin);

      double minX = proj_histo_->GetBinCenter(proj_histo_->FindFirstBinAbove(0.5*proj_histo_->GetMaximum()));
      double maxX = proj_histo_->GetBinCenter(proj_histo_->FindLastBinAbove(0.5*proj_histo_->GetMaximum()));
      double rangeX = 5*0.5*(std::abs(minX) + std::abs(maxX));
      proj_histo_->Fit(fgaus,"rQ","",-rangeX,+rangeX);
      proj_histo_->Fit(fgaus,"rQ","",-rangeX,+rangeX);

      double parA[3] = {fgaus->GetParameter(0),fgaus->GetParameter(1),fgaus->GetParameter(2)};
      double errA[3] = {fgaus->GetParError(0),fgaus->GetParError(1),fgaus->GetParError(2)};

      proj_histo_->Fit(fgaus,"rQ","", fgaus->GetParameter(1) - scale*fgaus->GetParameter(2), fgaus->GetParameter(1) + scale*fgaus->GetParameter(2));

      double parB[3] = {fgaus->GetParameter(0),fgaus->GetParameter(1),fgaus->GetParameter(2)};
      double errB[3] = {fgaus->GetParError(0),fgaus->GetParError(1),fgaus->GetParError(2)};

      if(std::abs(parA[1]) < std::abs(parB[1])) h_mean->SetBinContent(1+ibin, parA[1]);
      else h_mean->SetBinContent(1+ibin, parB[1]);
      h_mean->SetBinError(1+ibin, errA[1]);
   }
}


void fitterStaveToPhi(double scale, TH2D* histo, TH1D* h_mean, TF1* fgaus, TString tag = ""){
   //fitter
   int nbinsx = histo->GetNbinsX();
   for(int ibin = 0; ibin < nbinsx; ibin++){

      int jbin = ibin < 0.5*nbinsx ? ibin + 0.5*nbinsx : ibin - 0.5*nbinsx;

      fgaus->SetParameter(0,0);
      fgaus->SetParameter(1,0);
      fgaus->SetParameter(2,0);
      fgaus->SetParError(0,0);
      fgaus->SetParError(1,0);
      fgaus->SetParError(2,0);
      auto proj_histo_ = (TH1D*) histo->ProjectionY(Form("proj_histo%s_%d",(const char*)tag,ibin),1+ibin,1+ibin);

      double minX = proj_histo_->GetBinCenter(proj_histo_->FindFirstBinAbove(0.5*proj_histo_->GetMaximum()));
      double maxX = proj_histo_->GetBinCenter(proj_histo_->FindLastBinAbove(0.5*proj_histo_->GetMaximum()));
      double rangeX = 5*0.5*(std::abs(minX) + std::abs(maxX));
      proj_histo_->Fit(fgaus,"rQ","",-rangeX,+rangeX);
      proj_histo_->Fit(fgaus,"rQ","",-rangeX,+rangeX);

      double parA[3] = {fgaus->GetParameter(0),fgaus->GetParameter(1),fgaus->GetParameter(2)};
      double errA[3] = {fgaus->GetParError(0),fgaus->GetParError(1),fgaus->GetParError(2)};

      proj_histo_->Fit(fgaus,"rQ","", fgaus->GetParameter(1) - scale*fgaus->GetParameter(2), fgaus->GetParameter(1) + scale*fgaus->GetParameter(2));

      double parB[3] = {fgaus->GetParameter(0),fgaus->GetParameter(1),fgaus->GetParameter(2)};
      double errB[3] = {fgaus->GetParError(0),fgaus->GetParError(1),fgaus->GetParError(2)};

      if(std::abs(parA[1]) < std::abs(parB[1])) h_mean->SetBinContent(1+jbin, parA[1]);
      else h_mean->SetBinContent(1+jbin, parB[1]);
      h_mean->SetBinError(1+jbin, errA[1]);
   }
}

void check_track_residualsHB_FitGaus(int color = 1, double scale = 1.0){

TF1* fit_gaus = new TF1("fit_gaus", "gaus", -5000, 5000);

//detail phi
TH2D* hds1_vs_phi_0 = new TH2D("hds1_vs_phi_0","hds1_vs_phi_0",100,-TMath::Pi(),TMath::Pi(),1000,-100,100);
TH2D* hds1_vs_phi_1 = new TH2D("hds1_vs_phi_1","hds1_vs_phi_1",100,-TMath::Pi(),TMath::Pi(),1000,-100,100);
TH2D* hds1_vs_phi_2 = new TH2D("hds1_vs_phi_2","hds1_vs_phi_2",100,-TMath::Pi(),TMath::Pi(),1000,-100,100);
TH2D* hds1_vs_phi_3 = new TH2D("hds1_vs_phi_3","hds1_vs_phi_3",24,0,24,300,-1500,1500);
TH2D* hds1_vs_phi_4 = new TH2D("hds1_vs_phi_4","hds1_vs_phi_4",30,0,30,300,-1500,1500);
TH2D* hds1_vs_phi_5 = new TH2D("hds1_vs_phi_5","hds1_vs_phi_5",42,0,42,300,-3000,3000);
TH2D* hds1_vs_phi_6 = new TH2D("hds1_vs_phi_6","hds1_vs_phi_6",48,0,48,300,-3000,3000);

TH2D* hds2_vs_phi_0 = new TH2D("hds2_vs_phi_0","hds2_vs_phi_0",100,-TMath::Pi(),TMath::Pi(),1000,-100,100);
TH2D* hds2_vs_phi_1 = new TH2D("hds2_vs_phi_1","hds2_vs_phi_1",100,-TMath::Pi(),TMath::Pi(),1000,-100,100);
TH2D* hds2_vs_phi_2 = new TH2D("hds2_vs_phi_2","hds2_vs_phi_2",100,-TMath::Pi(),TMath::Pi(),1000,-100,100);
TH2D* hds2_vs_phi_3 = new TH2D("hds2_vs_phi_3","hds2_vs_phi_3",24,0,24,300,-1500,1500);
TH2D* hds2_vs_phi_4 = new TH2D("hds2_vs_phi_4","hds2_vs_phi_4",30,0,30,300,-1500,1500);
TH2D* hds2_vs_phi_5 = new TH2D("hds2_vs_phi_5","hds2_vs_phi_5",42,0,42,300,-3000,3000);
TH2D* hds2_vs_phi_6 = new TH2D("hds2_vs_phi_6","hds2_vs_phi_6",48,0,48,300,-3000,3000);

ResMonitor->Draw("fds1[0]*1e+4:TMath::ATan2(fgY[0],fgX[0])>>hds1_vs_phi_0","fds1[0]>-900");
ResMonitor->Draw("fds1[1]*1e+4:TMath::ATan2(fgY[1],fgX[1])>>hds1_vs_phi_1","fds1[1]>-900");
ResMonitor->Draw("fds1[2]*1e+4:TMath::ATan2(fgY[2],fgX[2])>>hds1_vs_phi_2","fds1[2]>-900");
ResMonitor->Draw("fds1[3]*1e+4:fstv[3]>>hds1_vs_phi_3","fds1[3]>-900");
ResMonitor->Draw("fds1[4]*1e+4:fstv[4]>>hds1_vs_phi_4","fds1[4]>-900");
ResMonitor->Draw("fds1[5]*1e+4:fstv[5]>>hds1_vs_phi_5","fds1[5]>-900");
ResMonitor->Draw("fds1[6]*1e+4:fstv[6]>>hds1_vs_phi_6","fds1[6]>-900");

ResMonitor->Draw("fds2[0]*1e+4:TMath::ATan2(fgY[0],fgX[0])>>hds2_vs_phi_0","fds2[0]>-900");
ResMonitor->Draw("fds2[1]*1e+4:TMath::ATan2(fgY[1],fgX[1])>>hds2_vs_phi_1","fds2[1]>-900");
ResMonitor->Draw("fds2[2]*1e+4:TMath::ATan2(fgY[2],fgX[2])>>hds2_vs_phi_2","fds2[2]>-900");
ResMonitor->Draw("fds2[3]*1e+4:fstv[3]>>hds2_vs_phi_3","fds2[3]>-900");
ResMonitor->Draw("fds2[4]*1e+4:fstv[4]>>hds2_vs_phi_4","fds2[4]>-900");
ResMonitor->Draw("fds2[5]*1e+4:fstv[5]>>hds2_vs_phi_5","fds2[5]>-900");
ResMonitor->Draw("fds2[6]*1e+4:fstv[6]>>hds2_vs_phi_6","fds2[6]>-900");

//detail z
TH2D* hds1_vs_z_0 = new TH2D("hds1_vs_z_0","hds1_vs_z_0",30,-15,15,1000,-100,100);
TH2D* hds1_vs_z_1 = new TH2D("hds1_vs_z_1","hds1_vs_z_1",30,-15,15,1000,-100,100);
TH2D* hds1_vs_z_2 = new TH2D("hds1_vs_z_2","hds1_vs_z_2",30,-15,15,1000,-100,100);
TH2D* hds1_vs_z_3 = new TH2D("hds1_vs_z_3","hds1_vs_z_3",30,-45,45,300,-1500,1500);
TH2D* hds1_vs_z_4 = new TH2D("hds1_vs_z_4","hds1_vs_z_4",30,-45,45,300,-1500,1500);
TH2D* hds1_vs_z_5 = new TH2D("hds1_vs_z_5","hds1_vs_z_5",30,-60,60,300,-3000,3000);
TH2D* hds1_vs_z_6 = new TH2D("hds1_vs_z_6","hds1_vs_z_6",30,-60,60,300,-3000,3000);

TH2D* hds2_vs_z_0 = new TH2D("hds2_vs_z_0","hds2_vs_z_0",30,-15,15,1000,-100,100);
TH2D* hds2_vs_z_1 = new TH2D("hds2_vs_z_1","hds2_vs_z_1",30,-15,15,1000,-100,100);
TH2D* hds2_vs_z_2 = new TH2D("hds2_vs_z_2","hds2_vs_z_2",30,-15,15,1000,-100,100);
TH2D* hds2_vs_z_3 = new TH2D("hds2_vs_z_3","hds2_vs_z_3",30,-45,45,300,-1500,1500);
TH2D* hds2_vs_z_4 = new TH2D("hds2_vs_z_4","hds2_vs_z_4",30,-45,45,300,-1500,1500);
TH2D* hds2_vs_z_5 = new TH2D("hds2_vs_z_5","hds2_vs_z_5",30,-60,60,300,-3000,3000);
TH2D* hds2_vs_z_6 = new TH2D("hds2_vs_z_6","hds2_vs_z_6",30,-60,60,300,-3000,3000);

TH2D* hds1_vs_z_0_HB0 = new TH2D("hds1_vs_z_0_HB0","hds1_vs_z_0_HB0",30,-15,15,1000,-100,100);
TH2D* hds1_vs_z_1_HB0 = new TH2D("hds1_vs_z_1_HB0","hds1_vs_z_1_HB0",30,-15,15,1000,-100,100);
TH2D* hds1_vs_z_2_HB0 = new TH2D("hds1_vs_z_2_HB0","hds1_vs_z_2_HB0",30,-15,15,1000,-100,100);
TH2D* hds1_vs_z_3_HB0 = new TH2D("hds1_vs_z_3_HB0","hds1_vs_z_3_HB0",30,-45,45,300,-1500,1500);
TH2D* hds1_vs_z_4_HB0 = new TH2D("hds1_vs_z_4_HB0","hds1_vs_z_4_HB0",30,-45,45,300,-1500,1500);
TH2D* hds1_vs_z_5_HB0 = new TH2D("hds1_vs_z_5_HB0","hds1_vs_z_5_HB0",30,-60,60,300,-3000,3000);
TH2D* hds1_vs_z_6_HB0 = new TH2D("hds1_vs_z_6_HB0","hds1_vs_z_6_HB0",30,-60,60,300,-3000,3000);

TH2D* hds2_vs_z_0_HB0 = new TH2D("hds2_vs_z_0_HB0","hds2_vs_z_0_HB0",30,-15,15,1000,-100,100);
TH2D* hds2_vs_z_1_HB0 = new TH2D("hds2_vs_z_1_HB0","hds2_vs_z_1_HB0",30,-15,15,1000,-100,100);
TH2D* hds2_vs_z_2_HB0 = new TH2D("hds2_vs_z_2_HB0","hds2_vs_z_2_HB0",30,-15,15,1000,-100,100);
TH2D* hds2_vs_z_3_HB0 = new TH2D("hds2_vs_z_3_HB0","hds2_vs_z_3_HB0",30,-45,45,300,-1500,1500);
TH2D* hds2_vs_z_4_HB0 = new TH2D("hds2_vs_z_4_HB0","hds2_vs_z_4_HB0",30,-45,45,300,-1500,1500);
TH2D* hds2_vs_z_5_HB0 = new TH2D("hds2_vs_z_5_HB0","hds2_vs_z_5_HB0",30,-60,60,300,-3000,3000);
TH2D* hds2_vs_z_6_HB0 = new TH2D("hds2_vs_z_6_HB0","hds2_vs_z_6_HB0",30,-60,60,300,-3000,3000);

TH2D* hds1_vs_z_0_HB1 = new TH2D("hds1_vs_z_0_HB1","hds1_vs_z_0_HB1",30,-15,15,1000,-100,100);
TH2D* hds1_vs_z_1_HB1 = new TH2D("hds1_vs_z_1_HB1","hds1_vs_z_1_HB1",30,-15,15,1000,-100,100);
TH2D* hds1_vs_z_2_HB1 = new TH2D("hds1_vs_z_2_HB1","hds1_vs_z_2_HB1",30,-15,15,1000,-100,100);
TH2D* hds1_vs_z_3_HB1 = new TH2D("hds1_vs_z_3_HB1","hds1_vs_z_3_HB1",30,-45,45,300,-1500,1500);
TH2D* hds1_vs_z_4_HB1 = new TH2D("hds1_vs_z_4_HB1","hds1_vs_z_4_HB1",30,-45,45,300,-1500,1500);
TH2D* hds1_vs_z_5_HB1 = new TH2D("hds1_vs_z_5_HB1","hds1_vs_z_5_HB1",30,-60,60,300,-3000,3000);
TH2D* hds1_vs_z_6_HB1 = new TH2D("hds1_vs_z_6_HB1","hds1_vs_z_6_HB1",30,-60,60,300,-3000,3000);

TH2D* hds2_vs_z_0_HB1 = new TH2D("hds2_vs_z_0_HB1","hds2_vs_z_0_HB1",30,-15,15,1000,-100,100);
TH2D* hds2_vs_z_1_HB1 = new TH2D("hds2_vs_z_1_HB1","hds2_vs_z_1_HB1",30,-15,15,1000,-100,100);
TH2D* hds2_vs_z_2_HB1 = new TH2D("hds2_vs_z_2_HB1","hds2_vs_z_2_HB1",30,-15,15,1000,-100,100);
TH2D* hds2_vs_z_3_HB1 = new TH2D("hds2_vs_z_3_HB1","hds2_vs_z_3_HB1",30,-45,45,300,-1500,1500);
TH2D* hds2_vs_z_4_HB1 = new TH2D("hds2_vs_z_4_HB1","hds2_vs_z_4_HB1",30,-45,45,300,-1500,1500);
TH2D* hds2_vs_z_5_HB1 = new TH2D("hds2_vs_z_5_HB1","hds2_vs_z_5_HB1",30,-60,60,300,-3000,3000);
TH2D* hds2_vs_z_6_HB1 = new TH2D("hds2_vs_z_6_HB1","hds2_vs_z_6_HB1",30,-60,60,300,-3000,3000);

ResMonitor->Draw("fds1[0]*1e+4:fgZ[0]>>hds1_vs_z_0","fds1[0]>-900 ");
ResMonitor->Draw("fds1[1]*1e+4:fgZ[1]>>hds1_vs_z_1","fds1[1]>-900 ");
ResMonitor->Draw("fds1[2]*1e+4:fgZ[2]>>hds1_vs_z_2","fds1[2]>-900 ");
ResMonitor->Draw("fds1[3]*1e+4:fgZ[3]>>hds1_vs_z_3","fds1[3]>-900 ");
ResMonitor->Draw("fds1[4]*1e+4:fgZ[4]>>hds1_vs_z_4","fds1[4]>-900 ");
ResMonitor->Draw("fds1[5]*1e+4:fgZ[5]>>hds1_vs_z_5","fds1[5]>-900 ");
ResMonitor->Draw("fds1[6]*1e+4:fgZ[6]>>hds1_vs_z_6","fds1[6]>-900 ");

ResMonitor->Draw("fds2[0]*1e+4:fgZ[0]>>hds2_vs_z_0","fds2[0]>-900 ");
ResMonitor->Draw("fds2[1]*1e+4:fgZ[1]>>hds2_vs_z_1","fds2[1]>-900 ");
ResMonitor->Draw("fds2[2]*1e+4:fgZ[2]>>hds2_vs_z_2","fds2[2]>-900 ");
ResMonitor->Draw("fds2[3]*1e+4:fgZ[3]>>hds2_vs_z_3","fds2[3]>-900 ");
ResMonitor->Draw("fds2[4]*1e+4:fgZ[4]>>hds2_vs_z_4","fds2[4]>-900 ");
ResMonitor->Draw("fds2[5]*1e+4:fgZ[5]>>hds2_vs_z_5","fds2[5]>-900 ");
ResMonitor->Draw("fds2[6]*1e+4:fgZ[6]>>hds2_vs_z_6","fds2[6]>-900 ");

ResMonitor->Draw("fds1[0]*1e+4:fgZ[0]>>hds1_vs_z_0_HB0","fds1[0]>-900 && TMath::ATan2(fgY[0],fgX[0])>0");
ResMonitor->Draw("fds1[1]*1e+4:fgZ[1]>>hds1_vs_z_1_HB0","fds1[1]>-900 && TMath::ATan2(fgY[1],fgX[1])>0");
ResMonitor->Draw("fds1[2]*1e+4:fgZ[2]>>hds1_vs_z_2_HB0","fds1[2]>-900 && TMath::ATan2(fgY[2],fgX[2])>0");
ResMonitor->Draw("fds1[3]*1e+4:fgZ[3]>>hds1_vs_z_3_HB0","fds1[3]>-900 && TMath::ATan2(fgY[3],fgX[3])>0");
ResMonitor->Draw("fds1[4]*1e+4:fgZ[4]>>hds1_vs_z_4_HB0","fds1[4]>-900 && TMath::ATan2(fgY[4],fgX[4])>0");
ResMonitor->Draw("fds1[5]*1e+4:fgZ[5]>>hds1_vs_z_5_HB0","fds1[5]>-900 && TMath::ATan2(fgY[5],fgX[5])>0");
ResMonitor->Draw("fds1[6]*1e+4:fgZ[6]>>hds1_vs_z_6_HB0","fds1[6]>-900 && TMath::ATan2(fgY[6],fgX[6])>0");

ResMonitor->Draw("fds2[0]*1e+4:fgZ[0]>>hds2_vs_z_0_HB0","fds2[0]>-900 && TMath::ATan2(fgY[0],fgX[0])>0");
ResMonitor->Draw("fds2[1]*1e+4:fgZ[1]>>hds2_vs_z_1_HB0","fds2[1]>-900 && TMath::ATan2(fgY[1],fgX[1])>0");
ResMonitor->Draw("fds2[2]*1e+4:fgZ[2]>>hds2_vs_z_2_HB0","fds2[2]>-900 && TMath::ATan2(fgY[2],fgX[2])>0");
ResMonitor->Draw("fds2[3]*1e+4:fgZ[3]>>hds2_vs_z_3_HB0","fds2[3]>-900 && TMath::ATan2(fgY[3],fgX[3])>0");
ResMonitor->Draw("fds2[4]*1e+4:fgZ[4]>>hds2_vs_z_4_HB0","fds2[4]>-900 && TMath::ATan2(fgY[4],fgX[4])>0");
ResMonitor->Draw("fds2[5]*1e+4:fgZ[5]>>hds2_vs_z_5_HB0","fds2[5]>-900 && TMath::ATan2(fgY[5],fgX[5])>0");
ResMonitor->Draw("fds2[6]*1e+4:fgZ[6]>>hds2_vs_z_6_HB0","fds2[6]>-900 && TMath::ATan2(fgY[6],fgX[6])>0");

ResMonitor->Draw("fds1[0]*1e+4:fgZ[0]>>hds1_vs_z_0_HB1","fds1[0]>-900 && TMath::ATan2(fgY[0],fgX[0])<0");
ResMonitor->Draw("fds1[1]*1e+4:fgZ[1]>>hds1_vs_z_1_HB1","fds1[1]>-900 && TMath::ATan2(fgY[1],fgX[1])<0");
ResMonitor->Draw("fds1[2]*1e+4:fgZ[2]>>hds1_vs_z_2_HB1","fds1[2]>-900 && TMath::ATan2(fgY[2],fgX[2])<0");
ResMonitor->Draw("fds1[3]*1e+4:fgZ[3]>>hds1_vs_z_3_HB1","fds1[3]>-900 && TMath::ATan2(fgY[3],fgX[3])<0");
ResMonitor->Draw("fds1[4]*1e+4:fgZ[4]>>hds1_vs_z_4_HB1","fds1[4]>-900 && TMath::ATan2(fgY[4],fgX[4])<0");
ResMonitor->Draw("fds1[5]*1e+4:fgZ[5]>>hds1_vs_z_5_HB1","fds1[5]>-900 && TMath::ATan2(fgY[5],fgX[5])<0");
ResMonitor->Draw("fds1[6]*1e+4:fgZ[6]>>hds1_vs_z_6_HB1","fds1[6]>-900 && TMath::ATan2(fgY[6],fgX[6])<0");

ResMonitor->Draw("fds2[0]*1e+4:fgZ[0]>>hds2_vs_z_0_HB1","fds2[0]>-900 && TMath::ATan2(fgY[0],fgX[0])<0");
ResMonitor->Draw("fds2[1]*1e+4:fgZ[1]>>hds2_vs_z_1_HB1","fds2[1]>-900 && TMath::ATan2(fgY[1],fgX[1])<0");
ResMonitor->Draw("fds2[2]*1e+4:fgZ[2]>>hds2_vs_z_2_HB1","fds2[2]>-900 && TMath::ATan2(fgY[2],fgX[2])<0");
ResMonitor->Draw("fds2[3]*1e+4:fgZ[3]>>hds2_vs_z_3_HB1","fds2[3]>-900 && TMath::ATan2(fgY[3],fgX[3])<0");
ResMonitor->Draw("fds2[4]*1e+4:fgZ[4]>>hds2_vs_z_4_HB1","fds2[4]>-900 && TMath::ATan2(fgY[4],fgX[4])<0");
ResMonitor->Draw("fds2[5]*1e+4:fgZ[5]>>hds2_vs_z_5_HB1","fds2[5]>-900 && TMath::ATan2(fgY[5],fgX[5])<0");
ResMonitor->Draw("fds2[6]*1e+4:fgZ[6]>>hds2_vs_z_6_HB1","fds2[6]>-900 && TMath::ATan2(fgY[6],fgX[6])<0");

TH1D* mhds1_vs_phi_0 = new TH1D("mhds1_vs_phi_0", "ds1 vs phi [0]; phi; ds1(#mum)", hds1_vs_phi_0->GetNbinsX(), hds1_vs_phi_0->GetXaxis()->GetXmin(), hds1_vs_phi_0->GetXaxis()->GetXmax()); 
mhds1_vs_phi_0->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds1_vs_phi_1 = new TH1D("mhds1_vs_phi_1", "ds1 vs phi [1]; phi; ds1(#mum)", hds1_vs_phi_1->GetNbinsX(), hds1_vs_phi_1->GetXaxis()->GetXmin(), hds1_vs_phi_1->GetXaxis()->GetXmax()); 
mhds1_vs_phi_1->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds1_vs_phi_2 = new TH1D("mhds1_vs_phi_2", "ds1 vs phi [2]; phi; ds1(#mum)", hds1_vs_phi_2->GetNbinsX(), hds1_vs_phi_2->GetXaxis()->GetXmin(), hds1_vs_phi_2->GetXaxis()->GetXmax()); 
mhds1_vs_phi_2->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds1_vs_phi_3 = new TH1D("mhds1_vs_phi_3", "ds1 vs phi [3]; phi; ds1(#mum)", hds1_vs_phi_3->GetNbinsX(), -TMath::Pi(), TMath::Pi()); 
mhds1_vs_phi_3->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds1_vs_phi_4 = new TH1D("mhds1_vs_phi_4", "ds1 vs phi [4]; phi; ds1(#mum)", hds1_vs_phi_4->GetNbinsX(), -TMath::Pi(), TMath::Pi()); 
mhds1_vs_phi_4->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds1_vs_phi_5 = new TH1D("mhds1_vs_phi_5", "ds1 vs phi [5]; phi; ds1(#mum)", hds1_vs_phi_5->GetNbinsX(), -TMath::Pi(), TMath::Pi()); 
mhds1_vs_phi_5->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds1_vs_phi_6 = new TH1D("mhds1_vs_phi_6", "ds1 vs phi [6]; phi; ds1(#mum)", hds1_vs_phi_6->GetNbinsX(), -TMath::Pi(), TMath::Pi()); 
mhds1_vs_phi_6->GetYaxis()->SetRangeUser(-200,200); 

TH1D* mhds2_vs_phi_0 = new TH1D("mhds2_vs_phi_0", "ds2 vs phi [0]; phi; ds2(#mum)", hds2_vs_phi_0->GetNbinsX(), hds2_vs_phi_0->GetXaxis()->GetXmin(), hds2_vs_phi_0->GetXaxis()->GetXmax()); 
mhds2_vs_phi_0->GetYaxis()->SetRangeUser(-20,20);
TH1D* mhds2_vs_phi_1 = new TH1D("mhds2_vs_phi_1", "ds2 vs phi [1]; phi; ds2(#mum)", hds2_vs_phi_1->GetNbinsX(), hds2_vs_phi_1->GetXaxis()->GetXmin(), hds2_vs_phi_1->GetXaxis()->GetXmax()); 
mhds2_vs_phi_1->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds2_vs_phi_2 = new TH1D("mhds2_vs_phi_2", "ds2 vs phi [2]; phi; ds2(#mum)", hds2_vs_phi_2->GetNbinsX(), hds2_vs_phi_2->GetXaxis()->GetXmin(), hds2_vs_phi_2->GetXaxis()->GetXmax()); 
mhds2_vs_phi_2->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds2_vs_phi_3 = new TH1D("mhds2_vs_phi_3", "ds2 vs phi [3]; phi; ds2(#mum)", hds2_vs_phi_3->GetNbinsX(), -TMath::Pi(), TMath::Pi()); 
mhds2_vs_phi_3->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds2_vs_phi_4 = new TH1D("mhds2_vs_phi_4", "ds2 vs phi [4]; phi; ds2(#mum)", hds2_vs_phi_4->GetNbinsX(), -TMath::Pi(), TMath::Pi()); 
mhds2_vs_phi_4->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds2_vs_phi_5 = new TH1D("mhds2_vs_phi_5", "ds2 vs phi [5]; phi; ds2(#mum)", hds2_vs_phi_5->GetNbinsX(), -TMath::Pi(), TMath::Pi()); 
mhds2_vs_phi_5->GetYaxis()->SetRangeUser(-400,400); 
TH1D* mhds2_vs_phi_6 = new TH1D("mhds2_vs_phi_6", "ds2 vs phi [6]; phi; ds2(#mum)", hds2_vs_phi_6->GetNbinsX(), -TMath::Pi(), TMath::Pi()); 
mhds2_vs_phi_6->GetYaxis()->SetRangeUser(-400,400); 

TH1D* mhds1_vs_z_0 = new TH1D("mhds1_vs_z_0", "ds1 vs z [0]; z(cm); ds1(#mum)", hds1_vs_z_0->GetNbinsX(), hds1_vs_z_0->GetXaxis()->GetXmin(), hds1_vs_z_0->GetXaxis()->GetXmax()); 
mhds1_vs_z_0->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds1_vs_z_1 = new TH1D("mhds1_vs_z_1", "ds1 vs z [0]; z(cm); ds1(#mum)", hds1_vs_z_1->GetNbinsX(), hds1_vs_z_1->GetXaxis()->GetXmin(), hds1_vs_z_1->GetXaxis()->GetXmax()); 
mhds1_vs_z_1->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds1_vs_z_2 = new TH1D("mhds1_vs_z_2", "ds1 vs z [0]; z(cm); ds1(#mum)", hds1_vs_z_2->GetNbinsX(), hds1_vs_z_2->GetXaxis()->GetXmin(), hds1_vs_z_2->GetXaxis()->GetXmax()); 
mhds1_vs_z_2->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds1_vs_z_3 = new TH1D("mhds1_vs_z_3", "ds1 vs z [0]; z(cm); ds1(#mum)", hds1_vs_z_3->GetNbinsX(), hds1_vs_z_3->GetXaxis()->GetXmin(), hds1_vs_z_3->GetXaxis()->GetXmax()); 
mhds1_vs_z_3->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds1_vs_z_4 = new TH1D("mhds1_vs_z_4", "ds1 vs z [0]; z(cm); ds1(#mum)", hds1_vs_z_4->GetNbinsX(), hds1_vs_z_4->GetXaxis()->GetXmin(), hds1_vs_z_4->GetXaxis()->GetXmax()); 
mhds1_vs_z_4->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds1_vs_z_5 = new TH1D("mhds1_vs_z_5", "ds1 vs z [0]; z(cm); ds1(#mum)", hds1_vs_z_5->GetNbinsX(), hds1_vs_z_5->GetXaxis()->GetXmin(), hds1_vs_z_5->GetXaxis()->GetXmax()); 
mhds1_vs_z_5->GetYaxis()->SetRangeUser(-400,400); 
TH1D* mhds1_vs_z_6 = new TH1D("mhds1_vs_z_6", "ds1 vs z [0]; z(cm); ds1(#mum)", hds1_vs_z_6->GetNbinsX(), hds1_vs_z_6->GetXaxis()->GetXmin(), hds1_vs_z_6->GetXaxis()->GetXmax()); 
mhds1_vs_z_6->GetYaxis()->SetRangeUser(-400,400); 

TH1D* mhds2_vs_z_0 = new TH1D("mhds2_vs_z_0", "ds2 vs z [0]; z(cm); ds2(#mum)", hds2_vs_z_0->GetNbinsX(), hds2_vs_z_0->GetXaxis()->GetXmin(), hds2_vs_z_0->GetXaxis()->GetXmax()); 
mhds2_vs_z_0->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds2_vs_z_1 = new TH1D("mhds2_vs_z_1", "ds2 vs z [1]; z(cm); ds2(#mum)", hds2_vs_z_1->GetNbinsX(), hds2_vs_z_1->GetXaxis()->GetXmin(), hds2_vs_z_1->GetXaxis()->GetXmax()); 
mhds2_vs_z_1->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds2_vs_z_2 = new TH1D("mhds2_vs_z_2", "ds2 vs z [2]; z(cm); ds2(#mum)", hds2_vs_z_2->GetNbinsX(), hds2_vs_z_2->GetXaxis()->GetXmin(), hds2_vs_z_2->GetXaxis()->GetXmax()); 
mhds2_vs_z_2->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds2_vs_z_3 = new TH1D("mhds2_vs_z_3", "ds2 vs z [3]; z(cm); ds2(#mum)", hds2_vs_z_3->GetNbinsX(), hds2_vs_z_3->GetXaxis()->GetXmin(), hds2_vs_z_3->GetXaxis()->GetXmax()); 
mhds2_vs_z_3->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds2_vs_z_4 = new TH1D("mhds2_vs_z_4", "ds2 vs z [4]; z(cm); ds2(#mum)", hds2_vs_z_4->GetNbinsX(), hds2_vs_z_4->GetXaxis()->GetXmin(), hds2_vs_z_4->GetXaxis()->GetXmax()); 
mhds2_vs_z_4->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds2_vs_z_5 = new TH1D("mhds2_vs_z_5", "ds2 vs z [5]; z(cm); ds2(#mum)", hds2_vs_z_5->GetNbinsX(), hds2_vs_z_5->GetXaxis()->GetXmin(), hds2_vs_z_5->GetXaxis()->GetXmax()); 
mhds2_vs_z_5->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds2_vs_z_6 = new TH1D("mhds2_vs_z_6", "ds2 vs z [6]; z(cm); ds2(#mum)", hds2_vs_z_6->GetNbinsX(), hds2_vs_z_6->GetXaxis()->GetXmin(), hds2_vs_z_6->GetXaxis()->GetXmax()); 
mhds2_vs_z_6->GetYaxis()->SetRangeUser(-200,200); 

TH1D* mhds1_vs_z_0_HB0 
= new TH1D("mhds1_vs_z_0_HB0", "ds1 vs z [0]; z(cm); ds1(#mum)", hds1_vs_z_0_HB0->GetNbinsX(), hds1_vs_z_0_HB0->GetXaxis()->GetXmin(), hds1_vs_z_0_HB0->GetXaxis()->GetXmax()); 
mhds1_vs_z_0_HB0->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds1_vs_z_1_HB0 
= new TH1D("mhds1_vs_z_1_HB0", "ds1 vs z [1]; z(cm); ds1(#mum)", hds1_vs_z_1_HB0->GetNbinsX(), hds1_vs_z_1_HB0->GetXaxis()->GetXmin(), hds1_vs_z_1_HB0->GetXaxis()->GetXmax()); 
mhds1_vs_z_1_HB0->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds1_vs_z_2_HB0 
= new TH1D("mhds1_vs_z_2_HB0", "ds1 vs z [2]; z(cm); ds1(#mum)", hds1_vs_z_2_HB0->GetNbinsX(), hds1_vs_z_2_HB0->GetXaxis()->GetXmin(), hds1_vs_z_2_HB0->GetXaxis()->GetXmax()); 
mhds1_vs_z_2_HB0->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds1_vs_z_3_HB0 
= new TH1D("mhds1_vs_z_3_HB0", "ds1 vs z [3]; z(cm); ds1(#mum)", hds1_vs_z_3_HB0->GetNbinsX(), hds1_vs_z_3_HB0->GetXaxis()->GetXmin(), hds1_vs_z_3_HB0->GetXaxis()->GetXmax()); 
mhds1_vs_z_3_HB0->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds1_vs_z_4_HB0 
= new TH1D("mhds1_vs_z_4_HB0", "ds1 vs z [4]; z(cm); ds1(#mum)", hds1_vs_z_4_HB0->GetNbinsX(), hds1_vs_z_4_HB0->GetXaxis()->GetXmin(), hds1_vs_z_4_HB0->GetXaxis()->GetXmax()); 
mhds1_vs_z_4_HB0->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds1_vs_z_5_HB0 
= new TH1D("mhds1_vs_z_5_HB0", "ds1 vs z [5]; z(cm); ds1(#mum)", hds1_vs_z_5_HB0->GetNbinsX(), hds1_vs_z_5_HB0->GetXaxis()->GetXmin(), hds1_vs_z_5_HB0->GetXaxis()->GetXmax()); 
mhds1_vs_z_5_HB0->GetYaxis()->SetRangeUser(-400,400); 
TH1D* mhds1_vs_z_6_HB0 
= new TH1D("mhds1_vs_z_6_HB0", "ds1 vs z [6]; z(cm); ds1(#mum)", hds1_vs_z_6_HB0->GetNbinsX(), hds1_vs_z_6_HB0->GetXaxis()->GetXmin(), hds1_vs_z_6_HB0->GetXaxis()->GetXmax()); 
mhds1_vs_z_6_HB0->GetYaxis()->SetRangeUser(-400,400); 

TH1D* mhds2_vs_z_0_HB0 
= new TH1D("mhds2_vs_z_0_HB0", "ds2 vs z [0]; z(cm); ds2(#mum)", hds2_vs_z_0_HB0->GetNbinsX(), hds2_vs_z_0_HB0->GetXaxis()->GetXmin(), hds2_vs_z_0_HB0->GetXaxis()->GetXmax()); 
mhds2_vs_z_0_HB0->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds2_vs_z_1_HB0 
= new TH1D("mhds2_vs_z_1_HB0", "ds2 vs z [1]; z(cm); ds2(#mum)", hds2_vs_z_1_HB0->GetNbinsX(), hds2_vs_z_1_HB0->GetXaxis()->GetXmin(), hds2_vs_z_1_HB0->GetXaxis()->GetXmax()); 
mhds2_vs_z_1_HB0->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds2_vs_z_2_HB0 
= new TH1D("mhds2_vs_z_2_HB0", "ds2 vs z [2]; z(cm); ds2(#mum)", hds2_vs_z_2_HB0->GetNbinsX(), hds2_vs_z_2_HB0->GetXaxis()->GetXmin(), hds2_vs_z_2_HB0->GetXaxis()->GetXmax()); 
mhds2_vs_z_2_HB0->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds2_vs_z_3_HB0 
= new TH1D("mhds2_vs_z_3_HB0", "ds2 vs z [3]; z(cm); ds2(#mum)", hds2_vs_z_3_HB0->GetNbinsX(), hds2_vs_z_3_HB0->GetXaxis()->GetXmin(), hds2_vs_z_3_HB0->GetXaxis()->GetXmax()); 
mhds2_vs_z_3_HB0->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds2_vs_z_4_HB0 
= new TH1D("mhds2_vs_z_4_HB0", "ds2 vs z [4]; z(cm); ds2(#mum)", hds2_vs_z_4_HB0->GetNbinsX(), hds2_vs_z_4_HB0->GetXaxis()->GetXmin(), hds2_vs_z_4_HB0->GetXaxis()->GetXmax()); 
mhds2_vs_z_4_HB0->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds2_vs_z_5_HB0 
= new TH1D("mhds2_vs_z_5_HB0", "ds2 vs z [5]; z(cm); ds2(#mum)", hds2_vs_z_5_HB0->GetNbinsX(), hds2_vs_z_5_HB0->GetXaxis()->GetXmin(), hds2_vs_z_5_HB0->GetXaxis()->GetXmax()); 
mhds2_vs_z_5_HB0->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds2_vs_z_6_HB0 
= new TH1D("mhds2_vs_z_6_HB0", "ds2 vs z [6]; z(cm); ds2(#mum)", hds2_vs_z_6_HB0->GetNbinsX(), hds2_vs_z_6_HB0->GetXaxis()->GetXmin(), hds2_vs_z_6_HB0->GetXaxis()->GetXmax()); 
mhds2_vs_z_6_HB0->GetYaxis()->SetRangeUser(-200,200); 

TH1D* mhds1_vs_z_0_HB1 
= new TH1D("mhds1_vs_z_0_HB1", "ds1 vs z [0]; z(cm); ds1(#mum)", hds1_vs_z_0_HB1->GetNbinsX(), hds1_vs_z_0_HB1->GetXaxis()->GetXmin(), hds1_vs_z_0_HB1->GetXaxis()->GetXmax()); 
mhds1_vs_z_0_HB1->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds1_vs_z_1_HB1 
= new TH1D("mhds1_vs_z_1_HB1", "ds1 vs z [1]; z(cm); ds1(#mum)", hds1_vs_z_1_HB1->GetNbinsX(), hds1_vs_z_1_HB1->GetXaxis()->GetXmin(), hds1_vs_z_1_HB1->GetXaxis()->GetXmax()); 
mhds1_vs_z_1_HB1->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds1_vs_z_2_HB1 
= new TH1D("mhds1_vs_z_2_HB1", "ds1 vs z [2]; z(cm); ds1(#mum)", hds1_vs_z_2_HB1->GetNbinsX(), hds1_vs_z_2_HB1->GetXaxis()->GetXmin(), hds1_vs_z_2_HB1->GetXaxis()->GetXmax()); 
mhds1_vs_z_2_HB1->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds1_vs_z_3_HB1 
= new TH1D("mhds1_vs_z_3_HB1", "ds1 vs z [3]; z(cm); ds1(#mum)", hds1_vs_z_3_HB1->GetNbinsX(), hds1_vs_z_3_HB1->GetXaxis()->GetXmin(), hds1_vs_z_3_HB1->GetXaxis()->GetXmax()); 
mhds1_vs_z_3_HB1->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds1_vs_z_4_HB1 
= new TH1D("mhds1_vs_z_4_HB1", "ds1 vs z [4]; z(cm); ds1(#mum)", hds1_vs_z_4_HB1->GetNbinsX(), hds1_vs_z_4_HB1->GetXaxis()->GetXmin(), hds1_vs_z_4_HB1->GetXaxis()->GetXmax()); 
mhds1_vs_z_4_HB1->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds1_vs_z_5_HB1 
= new TH1D("mhds1_vs_z_5_HB1", "ds1 vs z [5]; z(cm); ds1(#mum)", hds1_vs_z_5_HB1->GetNbinsX(), hds1_vs_z_5_HB1->GetXaxis()->GetXmin(), hds1_vs_z_5_HB1->GetXaxis()->GetXmax()); 
mhds1_vs_z_5_HB1->GetYaxis()->SetRangeUser(-400,400); 
TH1D* mhds1_vs_z_6_HB1 
= new TH1D("mhds1_vs_z_6_HB1", "ds1 vs z [6]; z(cm); ds1(#mum)", hds1_vs_z_6_HB1->GetNbinsX(), hds1_vs_z_6_HB1->GetXaxis()->GetXmin(), hds1_vs_z_6_HB1->GetXaxis()->GetXmax()); 
mhds1_vs_z_6_HB1->GetYaxis()->SetRangeUser(-400,400); 

TH1D* mhds2_vs_z_0_HB1 
= new TH1D("mhds2_vs_z_0_HB1", "ds2 vs z [0]; z(cm); ds2(#mum)", hds2_vs_z_0_HB1->GetNbinsX(), hds2_vs_z_0_HB1->GetXaxis()->GetXmin(), hds2_vs_z_0_HB1->GetXaxis()->GetXmax()); 
mhds2_vs_z_0_HB1->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds2_vs_z_1_HB1 
= new TH1D("mhds2_vs_z_1_HB1", "ds2 vs z [1]; z(cm); ds2(#mum)", hds2_vs_z_1_HB1->GetNbinsX(), hds2_vs_z_1_HB1->GetXaxis()->GetXmin(), hds2_vs_z_1_HB1->GetXaxis()->GetXmax()); 
mhds2_vs_z_1_HB1->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds2_vs_z_2_HB1 
= new TH1D("mhds2_vs_z_2_HB1", "ds2 vs z [2]; z(cm); ds2(#mum)", hds2_vs_z_2_HB1->GetNbinsX(), hds2_vs_z_2_HB1->GetXaxis()->GetXmin(), hds2_vs_z_2_HB1->GetXaxis()->GetXmax()); 
mhds2_vs_z_2_HB1->GetYaxis()->SetRangeUser(-20,20); 
TH1D* mhds2_vs_z_3_HB1 
= new TH1D("mhds2_vs_z_3_HB1", "ds2 vs z [3]; z(cm); ds2(#mum)", hds2_vs_z_3_HB1->GetNbinsX(), hds2_vs_z_3_HB1->GetXaxis()->GetXmin(), hds2_vs_z_3_HB1->GetXaxis()->GetXmax()); 
mhds2_vs_z_3_HB1->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds2_vs_z_4_HB1 
= new TH1D("mhds2_vs_z_4_HB1", "ds2 vs z [4]; z(cm); ds2(#mum)", hds2_vs_z_4_HB1->GetNbinsX(), hds2_vs_z_4_HB1->GetXaxis()->GetXmin(), hds2_vs_z_4_HB1->GetXaxis()->GetXmax()); 
mhds2_vs_z_4_HB1->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds2_vs_z_5_HB1 
= new TH1D("mhds2_vs_z_5_HB1", "ds2 vs z [5]; z(cm); ds2(#mum)", hds2_vs_z_5_HB1->GetNbinsX(), hds2_vs_z_5_HB1->GetXaxis()->GetXmin(), hds2_vs_z_5_HB1->GetXaxis()->GetXmax()); 
mhds2_vs_z_5_HB1->GetYaxis()->SetRangeUser(-200,200); 
TH1D* mhds2_vs_z_6_HB1 
= new TH1D("mhds2_vs_z_6_HB1", "ds2 vs z [6]; z(cm); ds2(#mum)", hds2_vs_z_6_HB1->GetNbinsX(), hds2_vs_z_6_HB1->GetXaxis()->GetXmin(), hds2_vs_z_6_HB1->GetXaxis()->GetXmax()); 
mhds2_vs_z_6_HB1->GetYaxis()->SetRangeUser(-200,200); 

mhds1_vs_phi_0->SetMarkerColor(color); 
mhds1_vs_phi_1->SetMarkerColor(color);
mhds1_vs_phi_2->SetMarkerColor(color);
mhds1_vs_phi_3->SetMarkerColor(color);
mhds1_vs_phi_4->SetMarkerColor(color);
mhds1_vs_phi_5->SetMarkerColor(color);
mhds1_vs_phi_6->SetMarkerColor(color);

mhds2_vs_phi_0->SetMarkerColor(color);
mhds2_vs_phi_1->SetMarkerColor(color);
mhds2_vs_phi_2->SetMarkerColor(color);
mhds2_vs_phi_3->SetMarkerColor(color);
mhds2_vs_phi_4->SetMarkerColor(color);
mhds2_vs_phi_5->SetMarkerColor(color);
mhds2_vs_phi_6->SetMarkerColor(color);

mhds1_vs_z_0->SetMarkerColor(color);
mhds1_vs_z_1->SetMarkerColor(color);
mhds1_vs_z_2->SetMarkerColor(color);
mhds1_vs_z_3->SetMarkerColor(color);
mhds1_vs_z_4->SetMarkerColor(color);
mhds1_vs_z_5->SetMarkerColor(color);
mhds1_vs_z_6->SetMarkerColor(color);

mhds2_vs_z_0->SetMarkerColor(color);
mhds2_vs_z_1->SetMarkerColor(color);
mhds2_vs_z_2->SetMarkerColor(color);
mhds2_vs_z_3->SetMarkerColor(color);
mhds2_vs_z_4->SetMarkerColor(color);
mhds2_vs_z_5->SetMarkerColor(color);
mhds2_vs_z_6->SetMarkerColor(color);

mhds1_vs_phi_0->SetLineColor(color);
mhds1_vs_phi_1->SetLineColor(color);
mhds1_vs_phi_2->SetLineColor(color);
mhds1_vs_phi_3->SetLineColor(color);
mhds1_vs_phi_4->SetLineColor(color);
mhds1_vs_phi_5->SetLineColor(color);
mhds1_vs_phi_6->SetLineColor(color);

mhds2_vs_phi_0->SetLineColor(color);
mhds2_vs_phi_1->SetLineColor(color);
mhds2_vs_phi_2->SetLineColor(color);
mhds2_vs_phi_3->SetLineColor(color);
mhds2_vs_phi_4->SetLineColor(color);
mhds2_vs_phi_5->SetLineColor(color);
mhds2_vs_phi_6->SetLineColor(color);

mhds1_vs_z_0->SetLineColor(color);
mhds1_vs_z_1->SetLineColor(color);
mhds1_vs_z_2->SetLineColor(color);
mhds1_vs_z_3->SetLineColor(color);
mhds1_vs_z_4->SetLineColor(color);
mhds1_vs_z_5->SetLineColor(color);
mhds1_vs_z_6->SetLineColor(color);

mhds2_vs_z_0->SetLineColor(color);
mhds2_vs_z_1->SetLineColor(color);
mhds2_vs_z_2->SetLineColor(color);
mhds2_vs_z_3->SetLineColor(color);
mhds2_vs_z_4->SetLineColor(color);
mhds2_vs_z_5->SetLineColor(color);
mhds2_vs_z_6->SetLineColor(color);

//HB
mhds1_vs_z_0_HB0->SetMarkerColor(color);
mhds1_vs_z_1_HB0->SetMarkerColor(color);
mhds1_vs_z_2_HB0->SetMarkerColor(color);
mhds1_vs_z_3_HB0->SetMarkerColor(color);
mhds1_vs_z_4_HB0->SetMarkerColor(color);
mhds1_vs_z_5_HB0->SetMarkerColor(color);
mhds1_vs_z_6_HB0->SetMarkerColor(color);

mhds2_vs_z_0_HB0->SetMarkerColor(color);
mhds2_vs_z_1_HB0->SetMarkerColor(color);
mhds2_vs_z_2_HB0->SetMarkerColor(color);
mhds2_vs_z_3_HB0->SetMarkerColor(color);
mhds2_vs_z_4_HB0->SetMarkerColor(color);
mhds2_vs_z_5_HB0->SetMarkerColor(color);
mhds2_vs_z_6_HB0->SetMarkerColor(color);

mhds1_vs_z_0_HB1->SetMarkerColor(color);
mhds1_vs_z_1_HB1->SetMarkerColor(color);
mhds1_vs_z_2_HB1->SetMarkerColor(color);
mhds1_vs_z_3_HB1->SetMarkerColor(color);
mhds1_vs_z_4_HB1->SetMarkerColor(color);
mhds1_vs_z_5_HB1->SetMarkerColor(color);
mhds1_vs_z_6_HB1->SetMarkerColor(color);

mhds2_vs_z_0_HB1->SetMarkerColor(color);
mhds2_vs_z_1_HB1->SetMarkerColor(color);
mhds2_vs_z_2_HB1->SetMarkerColor(color);
mhds2_vs_z_3_HB1->SetMarkerColor(color);
mhds2_vs_z_4_HB1->SetMarkerColor(color);
mhds2_vs_z_5_HB1->SetMarkerColor(color);
mhds2_vs_z_6_HB1->SetMarkerColor(color);

mhds1_vs_z_0_HB0->SetLineColor(color);
mhds1_vs_z_1_HB0->SetLineColor(color);
mhds1_vs_z_2_HB0->SetLineColor(color);
mhds1_vs_z_3_HB0->SetLineColor(color);
mhds1_vs_z_4_HB0->SetLineColor(color);
mhds1_vs_z_5_HB0->SetLineColor(color);
mhds1_vs_z_6_HB0->SetLineColor(color);

mhds2_vs_z_0_HB0->SetLineColor(color);
mhds2_vs_z_1_HB0->SetLineColor(color);
mhds2_vs_z_2_HB0->SetLineColor(color);
mhds2_vs_z_3_HB0->SetLineColor(color);
mhds2_vs_z_4_HB0->SetLineColor(color);
mhds2_vs_z_5_HB0->SetLineColor(color);
mhds2_vs_z_6_HB0->SetLineColor(color);

mhds1_vs_z_0_HB1->SetLineColor(color);
mhds1_vs_z_1_HB1->SetLineColor(color);
mhds1_vs_z_2_HB1->SetLineColor(color);
mhds1_vs_z_3_HB1->SetLineColor(color);
mhds1_vs_z_4_HB1->SetLineColor(color);
mhds1_vs_z_5_HB1->SetLineColor(color);
mhds1_vs_z_6_HB1->SetLineColor(color);

mhds2_vs_z_0_HB1->SetLineColor(color);
mhds2_vs_z_1_HB1->SetLineColor(color);
mhds2_vs_z_2_HB1->SetLineColor(color);
mhds2_vs_z_3_HB1->SetLineColor(color);
mhds2_vs_z_4_HB1->SetLineColor(color);
mhds2_vs_z_5_HB1->SetLineColor(color);
mhds2_vs_z_6_HB1->SetLineColor(color);

mhds1_vs_z_0_HB0->SetMarkerStyle(22);
mhds1_vs_z_1_HB0->SetMarkerStyle(22);
mhds1_vs_z_2_HB0->SetMarkerStyle(22);
mhds1_vs_z_3_HB0->SetMarkerStyle(22);
mhds1_vs_z_4_HB0->SetMarkerStyle(22);
mhds1_vs_z_5_HB0->SetMarkerStyle(22);
mhds1_vs_z_6_HB0->SetMarkerStyle(22);

mhds2_vs_z_0_HB0->SetMarkerStyle(22);
mhds2_vs_z_1_HB0->SetMarkerStyle(22);
mhds2_vs_z_2_HB0->SetMarkerStyle(22);
mhds2_vs_z_3_HB0->SetMarkerStyle(22);
mhds2_vs_z_4_HB0->SetMarkerStyle(22);
mhds2_vs_z_5_HB0->SetMarkerStyle(22);
mhds2_vs_z_6_HB0->SetMarkerStyle(22);

mhds1_vs_z_0_HB1->SetMarkerStyle(23);
mhds1_vs_z_1_HB1->SetMarkerStyle(23);
mhds1_vs_z_2_HB1->SetMarkerStyle(23);
mhds1_vs_z_3_HB1->SetMarkerStyle(23);
mhds1_vs_z_4_HB1->SetMarkerStyle(23);
mhds1_vs_z_5_HB1->SetMarkerStyle(23);
mhds1_vs_z_6_HB1->SetMarkerStyle(23);

mhds2_vs_z_0_HB1->SetMarkerStyle(23);
mhds2_vs_z_1_HB1->SetMarkerStyle(23);
mhds2_vs_z_2_HB1->SetMarkerStyle(23);
mhds2_vs_z_3_HB1->SetMarkerStyle(23);
mhds2_vs_z_4_HB1->SetMarkerStyle(23);
mhds2_vs_z_5_HB1->SetMarkerStyle(23);
mhds2_vs_z_6_HB1->SetMarkerStyle(23);

//FitGaus
fitter(scale, hds1_vs_phi_0, mhds1_vs_phi_0, fit_gaus, "ds1_vs_phi_0");
fitter(scale, hds1_vs_phi_1, mhds1_vs_phi_1, fit_gaus, "ds1_vs_phi_1");
fitter(scale, hds1_vs_phi_2, mhds1_vs_phi_2, fit_gaus, "ds1_vs_phi_2");
fitterStaveToPhi(scale, hds1_vs_phi_3, mhds1_vs_phi_3, fit_gaus, "ds1_vs_phi_3");
fitterStaveToPhi(scale, hds1_vs_phi_4, mhds1_vs_phi_4, fit_gaus, "ds1_vs_phi_4");
fitterStaveToPhi(scale, hds1_vs_phi_5, mhds1_vs_phi_5, fit_gaus, "ds1_vs_phi_5");
fitterStaveToPhi(scale, hds1_vs_phi_6, mhds1_vs_phi_6, fit_gaus, "ds1_vs_phi_6");

fitter(scale, hds2_vs_phi_0, mhds2_vs_phi_0, fit_gaus, "ds2_vs_phi_0");
fitter(scale, hds2_vs_phi_1, mhds2_vs_phi_1, fit_gaus, "ds2_vs_phi_1");
fitter(scale, hds2_vs_phi_2, mhds2_vs_phi_2, fit_gaus, "ds2_vs_phi_2");
fitterStaveToPhi(scale, hds2_vs_phi_3, mhds2_vs_phi_3, fit_gaus, "ds2_vs_phi_3");
fitterStaveToPhi(scale, hds2_vs_phi_4, mhds2_vs_phi_4, fit_gaus, "ds2_vs_phi_4");
fitterStaveToPhi(scale, hds2_vs_phi_5, mhds2_vs_phi_5, fit_gaus, "ds2_vs_phi_5");
fitterStaveToPhi(scale, hds2_vs_phi_6, mhds2_vs_phi_6, fit_gaus, "ds2_vs_phi_6");

fitter(scale, hds1_vs_z_0, mhds1_vs_z_0, fit_gaus, "ds1_vs_z_0");
fitter(scale, hds1_vs_z_1, mhds1_vs_z_1, fit_gaus, "ds1_vs_z_1");
fitter(scale, hds1_vs_z_2, mhds1_vs_z_2, fit_gaus, "ds1_vs_z_2");
fitter(scale, hds1_vs_z_3, mhds1_vs_z_3, fit_gaus, "ds1_vs_z_3");
fitter(scale, hds1_vs_z_4, mhds1_vs_z_4, fit_gaus, "ds1_vs_z_4");
fitter(scale, hds1_vs_z_5, mhds1_vs_z_5, fit_gaus, "ds1_vs_z_5");
fitter(scale, hds1_vs_z_6, mhds1_vs_z_6, fit_gaus, "ds1_vs_z_6");

fitter(scale, hds2_vs_z_0, mhds2_vs_z_0, fit_gaus, "ds2_vs_z_0");
fitter(scale, hds2_vs_z_1, mhds2_vs_z_1, fit_gaus, "ds2_vs_z_1");
fitter(scale, hds2_vs_z_2, mhds2_vs_z_2, fit_gaus, "ds2_vs_z_2");
fitter(scale, hds2_vs_z_3, mhds2_vs_z_3, fit_gaus, "ds2_vs_z_3");
fitter(scale, hds2_vs_z_4, mhds2_vs_z_4, fit_gaus, "ds2_vs_z_4");
fitter(scale, hds2_vs_z_5, mhds2_vs_z_5, fit_gaus, "ds2_vs_z_5");
fitter(scale, hds2_vs_z_6, mhds2_vs_z_6, fit_gaus, "ds2_vs_z_6");

//HB
fitter(scale, hds1_vs_z_0_HB0, mhds1_vs_z_0_HB0, fit_gaus, "ds1_vs_z_0_HB0");
fitter(scale, hds1_vs_z_1_HB0, mhds1_vs_z_1_HB0, fit_gaus, "ds1_vs_z_1_HB0");
fitter(scale, hds1_vs_z_2_HB0, mhds1_vs_z_2_HB0, fit_gaus, "ds1_vs_z_2_HB0");
fitter(scale, hds1_vs_z_3_HB0, mhds1_vs_z_3_HB0, fit_gaus, "ds1_vs_z_3_HB0");
fitter(scale, hds1_vs_z_4_HB0, mhds1_vs_z_4_HB0, fit_gaus, "ds1_vs_z_4_HB0");
fitter(scale, hds1_vs_z_5_HB0, mhds1_vs_z_5_HB0, fit_gaus, "ds1_vs_z_5_HB0");
fitter(scale, hds1_vs_z_6_HB0, mhds1_vs_z_6_HB0, fit_gaus, "ds1_vs_z_6_HB0");

fitter(scale, hds2_vs_z_0_HB0, mhds2_vs_z_0_HB0, fit_gaus, "ds2_vs_z_0_HB0");
fitter(scale, hds2_vs_z_1_HB0, mhds2_vs_z_1_HB0, fit_gaus, "ds2_vs_z_1_HB0");
fitter(scale, hds2_vs_z_2_HB0, mhds2_vs_z_2_HB0, fit_gaus, "ds2_vs_z_2_HB0");
fitter(scale, hds2_vs_z_3_HB0, mhds2_vs_z_3_HB0, fit_gaus, "ds2_vs_z_3_HB0");
fitter(scale, hds2_vs_z_4_HB0, mhds2_vs_z_4_HB0, fit_gaus, "ds2_vs_z_4_HB0");
fitter(scale, hds2_vs_z_5_HB0, mhds2_vs_z_5_HB0, fit_gaus, "ds2_vs_z_5_HB0");
fitter(scale, hds2_vs_z_6_HB0, mhds2_vs_z_6_HB0, fit_gaus, "ds2_vs_z_6_HB0");

fitter(scale, hds1_vs_z_0_HB1, mhds1_vs_z_0_HB1, fit_gaus, "ds1_vs_z_0_HB1");
fitter(scale, hds1_vs_z_1_HB1, mhds1_vs_z_1_HB1, fit_gaus, "ds1_vs_z_1_HB1");
fitter(scale, hds1_vs_z_2_HB1, mhds1_vs_z_2_HB1, fit_gaus, "ds1_vs_z_2_HB1");
fitter(scale, hds1_vs_z_3_HB1, mhds1_vs_z_3_HB1, fit_gaus, "ds1_vs_z_3_HB1");
fitter(scale, hds1_vs_z_4_HB1, mhds1_vs_z_4_HB1, fit_gaus, "ds1_vs_z_4_HB1");
fitter(scale, hds1_vs_z_5_HB1, mhds1_vs_z_5_HB1, fit_gaus, "ds1_vs_z_5_HB1");
fitter(scale, hds1_vs_z_6_HB1, mhds1_vs_z_6_HB1, fit_gaus, "ds1_vs_z_6_HB1");

fitter(scale, hds2_vs_z_0_HB1, mhds2_vs_z_0_HB1, fit_gaus, "ds2_vs_z_0_HB1");
fitter(scale, hds2_vs_z_1_HB1, mhds2_vs_z_1_HB1, fit_gaus, "ds2_vs_z_1_HB1");
fitter(scale, hds2_vs_z_2_HB1, mhds2_vs_z_2_HB1, fit_gaus, "ds2_vs_z_2_HB1");
fitter(scale, hds2_vs_z_3_HB1, mhds2_vs_z_3_HB1, fit_gaus, "ds2_vs_z_3_HB1");
fitter(scale, hds2_vs_z_4_HB1, mhds2_vs_z_4_HB1, fit_gaus, "ds2_vs_z_4_HB1");
fitter(scale, hds2_vs_z_5_HB1, mhds2_vs_z_5_HB1, fit_gaus, "ds2_vs_z_5_HB1");
fitter(scale, hds2_vs_z_6_HB1, mhds2_vs_z_6_HB1, fit_gaus, "ds2_vs_z_6_HB1");

TCanvas* canvas_ds1_vs_phi = new TCanvas("canvas_ds1_vs_phi","canvas_ds1_vs_phi",1600,800);
canvas_ds1_vs_phi->Divide(4,2);
canvas_ds1_vs_phi->cd(1); gStyle->SetOptStat(0); mhds1_vs_phi_0->Draw();
canvas_ds1_vs_phi->cd(2); gStyle->SetOptStat(0); mhds1_vs_phi_1->Draw();
canvas_ds1_vs_phi->cd(3); gStyle->SetOptStat(0); mhds1_vs_phi_2->Draw();
canvas_ds1_vs_phi->cd(5); gStyle->SetOptStat(0); mhds1_vs_phi_3->Draw();
canvas_ds1_vs_phi->cd(6); gStyle->SetOptStat(0); mhds1_vs_phi_4->Draw();
canvas_ds1_vs_phi->cd(7); gStyle->SetOptStat(0); mhds1_vs_phi_5->Draw();
canvas_ds1_vs_phi->cd(8); gStyle->SetOptStat(0); mhds1_vs_phi_6->Draw();

TCanvas* canvas_ds2_vs_phi = new TCanvas("canvas_ds2_vs_phi","canvas_ds2_vs_phi",1600,800);
canvas_ds2_vs_phi->Divide(4,2);
canvas_ds2_vs_phi->cd(1); gStyle->SetOptStat(0); mhds2_vs_phi_0->Draw();
canvas_ds2_vs_phi->cd(2); gStyle->SetOptStat(0); mhds2_vs_phi_1->Draw();
canvas_ds2_vs_phi->cd(3); gStyle->SetOptStat(0); mhds2_vs_phi_2->Draw();
canvas_ds2_vs_phi->cd(5); gStyle->SetOptStat(0); mhds2_vs_phi_3->Draw();
canvas_ds2_vs_phi->cd(6); gStyle->SetOptStat(0); mhds2_vs_phi_4->Draw();
canvas_ds2_vs_phi->cd(7); gStyle->SetOptStat(0); mhds2_vs_phi_5->Draw();
canvas_ds2_vs_phi->cd(8); gStyle->SetOptStat(0); mhds2_vs_phi_6->Draw();

TCanvas* canvas_ds1_vs_z = new TCanvas("canvas_ds1_vs_z","canvas_ds1_vs_z",1600,800);
canvas_ds1_vs_z->Divide(4,2);
canvas_ds1_vs_z->cd(1); gStyle->SetOptStat(0); mhds1_vs_z_0->Draw(); mhds1_vs_z_0_HB0->Draw("same"); mhds1_vs_z_0_HB1->Draw("same");
canvas_ds1_vs_z->cd(2); gStyle->SetOptStat(0); mhds1_vs_z_1->Draw(); mhds1_vs_z_1_HB0->Draw("same"); mhds1_vs_z_1_HB1->Draw("same");
canvas_ds1_vs_z->cd(3); gStyle->SetOptStat(0); mhds1_vs_z_2->Draw(); mhds1_vs_z_2_HB0->Draw("same"); mhds1_vs_z_2_HB1->Draw("same");
canvas_ds1_vs_z->cd(5); gStyle->SetOptStat(0); mhds1_vs_z_3->Draw(); mhds1_vs_z_3_HB0->Draw("same"); mhds1_vs_z_3_HB1->Draw("same");
canvas_ds1_vs_z->cd(6); gStyle->SetOptStat(0); mhds1_vs_z_4->Draw(); mhds1_vs_z_4_HB0->Draw("same"); mhds1_vs_z_4_HB1->Draw("same");
canvas_ds1_vs_z->cd(7); gStyle->SetOptStat(0); mhds1_vs_z_5->Draw(); mhds1_vs_z_5_HB0->Draw("same"); mhds1_vs_z_5_HB1->Draw("same");
canvas_ds1_vs_z->cd(8); gStyle->SetOptStat(0); mhds1_vs_z_6->Draw(); mhds1_vs_z_6_HB0->Draw("same"); mhds1_vs_z_6_HB1->Draw("same");

TCanvas* canvas_ds2_vs_z = new TCanvas("canvas_ds2_vs_z","canvas_ds2_vs_z",1600,800);
canvas_ds2_vs_z->Divide(4,2);
canvas_ds2_vs_z->cd(1); gStyle->SetOptStat(0); mhds2_vs_z_0->Draw(); mhds2_vs_z_0_HB0->Draw("same"); mhds2_vs_z_0_HB1->Draw("same");
canvas_ds2_vs_z->cd(2); gStyle->SetOptStat(0); mhds2_vs_z_1->Draw(); mhds2_vs_z_1_HB0->Draw("same"); mhds2_vs_z_1_HB1->Draw("same");
canvas_ds2_vs_z->cd(3); gStyle->SetOptStat(0); mhds2_vs_z_2->Draw(); mhds2_vs_z_2_HB0->Draw("same"); mhds2_vs_z_2_HB1->Draw("same");
canvas_ds2_vs_z->cd(5); gStyle->SetOptStat(0); mhds2_vs_z_3->Draw(); mhds2_vs_z_3_HB0->Draw("same"); mhds2_vs_z_3_HB1->Draw("same");
canvas_ds2_vs_z->cd(6); gStyle->SetOptStat(0); mhds2_vs_z_4->Draw(); mhds2_vs_z_4_HB0->Draw("same"); mhds2_vs_z_4_HB1->Draw("same");
canvas_ds2_vs_z->cd(7); gStyle->SetOptStat(0); mhds2_vs_z_5->Draw(); mhds2_vs_z_5_HB0->Draw("same"); mhds2_vs_z_5_HB1->Draw("same");
canvas_ds2_vs_z->cd(8); gStyle->SetOptStat(0); mhds2_vs_z_6->Draw(); mhds2_vs_z_6_HB0->Draw("same"); mhds2_vs_z_6_HB1->Draw("same");
/**/
};
