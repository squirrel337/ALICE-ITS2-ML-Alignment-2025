#define cm_to_mu 1.0e+4
void check_vertex_eta_plots_trkvtxer_color(int color=-1, int markerstyle=21, double scale = 1.0){
   new TCanvas();
   const int nbinsx = 20;

   TH2D* h_event_Vx_eta = new TH2D("h_event_Vx_eta","Vx(cm) vs eta;eta;x(cm)", nbinsx, -1.0, +1.0, 2000,-0.1,0.1);
   TH2D* h_event_Vy_eta = new TH2D("h_event_Vy_eta","Vy(cm) vs eta;eta;y(cm)", nbinsx, -1.0, +1.0, 2000,-0.1,0.1);
   TH2D* h_event_Vz_eta = new TH2D("h_event_Vz_eta","Vz(cm) vs eta;eta;z(cm)", nbinsx, -1.0, +1.0, 1500,-15,15);

   TrkVtxer->Draw("vtxevtX:eta>>h_event_Vx_eta");
   TrkVtxer->Draw("vtxevtY:eta>>h_event_Vy_eta");
   TrkVtxer->Draw("vtxevtZ:eta>>h_event_Vz_eta");

   TH2D* hdcaY_eta   = new TH2D("hdcaY_eta",
                                "dcaY vs eta;eta;dcaY(#mum)",
                                nbinsx, -1.0, +1.0, 500, -0.05*cm_to_mu, 0.05*cm_to_mu);
   TProfile* pf_width_dcaY = new TProfile("pf_width_dcaY","pf_width_dcaY",nbinsx,-1.0,1.0,-0.02*cm_to_mu,0.02*cm_to_mu,"s");
   TProfile* pf_widthE_dcaY = new TProfile("pf_widthE_dcaY","pf_widthE_dcaY",nbinsx,-1.0,1.0,-0.02*cm_to_mu,0.02*cm_to_mu,"g");
   TF1* fgaus_dcaY = new TF1("fgaus_dcaY","gaus",-0.01*cm_to_mu,0.01*cm_to_mu);

   TH2D* hdcaZ_eta   = new TH2D("hdcaZ_eta",
                                "dcaZ vs eta;eta;dcaZ(#mum)", 
                                nbinsx, -1.0, +1.0, 500, -0.05*cm_to_mu, 0.05*cm_to_mu);
   TH2D* hdcaYZ_eta  = new TH2D("hdcaYZ_eta",
                                "#sqrt{dcaY^{2}+dcaZ^{2}} vs eta;eta;#sqrt{dcaY^{2}+dcaZ^{2}}(#mum)", 
                                nbinsx, -1.0, +1.0,1000,0*cm_to_mu,0.1*cm_to_mu);
   TProfile* pf_width_dcaZ = new TProfile("pf_width_dcaZ","pf_width_dcaZ",nbinsx,-1.0,1.0,-0.02*cm_to_mu,0.02*cm_to_mu,"s");
   TProfile* pf_widthE_dcaZ = new TProfile("pf_widthE_dcaZ","pf_widthE_dcaZ",nbinsx,-1.0,1.0,-0.02*cm_to_mu,0.02*cm_to_mu,"g");
   TF1* fgaus_dcaZ = new TF1("fgaus_dcaZ","gaus",-0.01*cm_to_mu,0.01*cm_to_mu);

   TrkVtxer->Draw(Form("fip[0]*%g:eta>>hdcaY_eta",cm_to_mu),"");
   TrkVtxer->Draw(Form("fip[0]*%g:eta>>pf_width_dcaY",cm_to_mu),"");
   TrkVtxer->Draw(Form("fip[0]*%g:eta>>pf_widthE_dcaY",cm_to_mu),"");

   TrkVtxer->Draw(Form("fip[1]*%g:eta>>hdcaZ_eta",cm_to_mu),"");
   TrkVtxer->Draw(Form("fip[1]*%g:eta>>pf_width_dcaZ",cm_to_mu),"");
   TrkVtxer->Draw(Form("fip[1]*%g:eta>>pf_widthE_dcaZ",cm_to_mu),"");

   TrkVtxer->Draw(Form("TMath::Sqrt(fip[0]*fip[0]+fip[1]*fip[1])*%g:eta>>hdcaYZ_eta",cm_to_mu),"");

   //fitter
   TH1D* h_mean_dcaY = new TH1D("h_mean_dcaY","(mean) dcaY vs #eta;#eta;dcaY(#mum)",nbinsx,-1.0,1.0);
   TH1D* h_width_dcaY = new TH1D("h_width_dcaY","(width) dcaY vs #eta;#eta;#sigma(#mum)",nbinsx,-1.0,1.0);
   for(int ibin = 1; ibin < nbinsx-1; ibin++){
      fgaus_dcaY->SetParameter(0,0);
      fgaus_dcaY->SetParameter(1,0);
      fgaus_dcaY->SetParameter(2,0);
      fgaus_dcaY->SetParError(0,0);
      fgaus_dcaY->SetParError(1,0);
      fgaus_dcaY->SetParError(2,0);
      auto proj_hdcaY_eta = (TH1D*) hdcaY_eta->ProjectionY(Form("proj_hdcaY_eta_%d",ibin),1+ibin,1+ibin);

      double minX = proj_hdcaY_eta->GetBinCenter(proj_hdcaY_eta->FindFirstBinAbove(0.5*proj_hdcaY_eta->GetMaximum()));
      double maxX = proj_hdcaY_eta->GetBinCenter(proj_hdcaY_eta->FindLastBinAbove(0.5*proj_hdcaY_eta->GetMaximum()));
      double rangeX = 5*0.5*(std::abs(minX) + std::abs(maxX));
      proj_hdcaY_eta->Fit(fgaus_dcaY,"rQ","",-rangeX,+rangeX);
      proj_hdcaY_eta->Fit(fgaus_dcaY,"rQ","",-rangeX,+rangeX);

      double parA[3] = {fgaus_dcaY->GetParameter(0),fgaus_dcaY->GetParameter(1),fgaus_dcaY->GetParameter(2)};
      double errA[3] = {fgaus_dcaY->GetParError(0),fgaus_dcaY->GetParError(1),fgaus_dcaY->GetParError(2)};

      proj_hdcaY_eta->Fit(fgaus_dcaY,"rQ","", fgaus_dcaY->GetParameter(1) - scale*fgaus_dcaY->GetParameter(2), fgaus_dcaY->GetParameter(1) + scale*fgaus_dcaY->GetParameter(2));

      double parB[3] = {fgaus_dcaY->GetParameter(0),fgaus_dcaY->GetParameter(1),fgaus_dcaY->GetParameter(2)};
      double errB[3] = {fgaus_dcaY->GetParError(0),fgaus_dcaY->GetParError(1),fgaus_dcaY->GetParError(2)};

      if(std::abs(parA[1]) < std::abs(parB[1])) h_mean_dcaY->SetBinContent(1+ibin, parA[1]);
      else h_mean_dcaY->SetBinContent(1+ibin, parB[1]);
      h_mean_dcaY->SetBinError(1+ibin, errA[1]);

      h_width_dcaY->SetBinContent(1+ibin, parA[2]);
      h_width_dcaY->SetBinError(1+ibin, errA[2]);

      //h_mean_dcaY->SetBinContent(1+ibin, fgaus_dcaY->GetParameter(1));
      //h_mean_dcaY->SetBinError(1+ibin, fgaus_dcaY->GetParError(1));

      //h_width_dcaY->SetBinContent(1+ibin, fgaus_dcaY->GetParameter(2));
      //h_width_dcaY->SetBinError(1+ibin, fgaus_dcaY->GetParError(2));
   }

   TH1D* h_mean_dcaZ = new TH1D("h_mean_dcaZ","(mean) dcaZ vs #eta;#eta;dcaZ(#mum)",nbinsx,-1.0,1.0);
   TH1D* h_width_dcaZ = new TH1D("h_width_dcaZ","(width) dcaZ vs #eta;#eta;#sigma(#mum)",nbinsx,-1.0,1.0);
   for(int ibin = 1; ibin < nbinsx-1; ibin++){
      fgaus_dcaZ->SetParameter(0,0);
      fgaus_dcaZ->SetParameter(1,0);
      fgaus_dcaZ->SetParameter(2,0);
      fgaus_dcaZ->SetParError(0,0);
      fgaus_dcaZ->SetParError(1,0);
      fgaus_dcaZ->SetParError(2,0);
      auto proj_hdcaZ_eta = (TH1D*) hdcaZ_eta->ProjectionY(Form("proj_hdcaZ_eta_%d",ibin),1+ibin,1+ibin);

      double minX = proj_hdcaZ_eta->GetBinCenter(proj_hdcaZ_eta->FindFirstBinAbove(0.5*proj_hdcaZ_eta->GetMaximum()));
      double maxX = proj_hdcaZ_eta->GetBinCenter(proj_hdcaZ_eta->FindLastBinAbove(0.5*proj_hdcaZ_eta->GetMaximum()));
      double rangeX = 5*0.5*(std::abs(minX) + std::abs(maxX));
      proj_hdcaZ_eta->Fit(fgaus_dcaZ,"rQ","",-rangeX,+rangeX);
      proj_hdcaZ_eta->Fit(fgaus_dcaZ,"rQ","",-rangeX,+rangeX);

      double parA[3] = {fgaus_dcaZ->GetParameter(0),fgaus_dcaZ->GetParameter(1),fgaus_dcaZ->GetParameter(2)};
      double errA[3] = {fgaus_dcaZ->GetParError(0),fgaus_dcaZ->GetParError(1),fgaus_dcaZ->GetParError(2)};

      proj_hdcaZ_eta->Fit(fgaus_dcaZ,"rQ","", fgaus_dcaZ->GetParameter(1) - scale*fgaus_dcaZ->GetParameter(2), fgaus_dcaZ->GetParameter(1) + scale*fgaus_dcaZ->GetParameter(2));

      double parB[3] = {fgaus_dcaZ->GetParameter(0),fgaus_dcaZ->GetParameter(1),fgaus_dcaZ->GetParameter(2)};
      double errB[3] = {fgaus_dcaZ->GetParError(0),fgaus_dcaZ->GetParError(1),fgaus_dcaZ->GetParError(2)};

      if(std::abs(parA[1]) < std::abs(parB[1])) h_mean_dcaZ->SetBinContent(1+ibin, parA[1]);
      else h_mean_dcaZ->SetBinContent(1+ibin, parB[1]);
      h_mean_dcaZ->SetBinError(1+ibin, errA[1]);

      h_width_dcaZ->SetBinContent(1+ibin, parA[2]);
      h_width_dcaZ->SetBinError(1+ibin, errA[2]);

      //h_mean_dcaZ->SetBinContent(1+ibin, fgaus_dcaZ->GetParameter(1));
      //h_mean_dcaZ->SetBinError(1+ibin, fgaus_dcaZ->GetParError(1));

      //h_width_dcaZ->SetBinContent(1+ibin, fgaus_dcaZ->GetParameter(2));
      //h_width_dcaZ->SetBinError(1+ibin, fgaus_dcaZ->GetParError(2));
   }
/*
   auto canvas_vertex = new TCanvas("canvas_vertex","canvas_vertex",1600,800);
   canvas_vertex->Divide(3,2,1e-3,1e-3);
   for(int ic = 0; ic < 6; ic++){
      canvas_vertex->cd(1+ic);
      gStyle->SetOptStat(0);
      gPad->SetMargin(0.12,0.12,0.1,0.1);
   }
   canvas_vertex->cd(1);
   h_event_Vx_eta->Draw("colz");

   canvas_vertex->cd(2);
   h_event_Vy_eta->Draw("colz");

   canvas_vertex->cd(3);
   h_event_Vz_eta->Draw("colz");

   canvas_vertex->cd(4);
   hdcaY_eta->Draw("colz");

   canvas_vertex->cd(5);
   hdcaZ_eta->Draw("colz");

   canvas_vertex->cd(6);
   hdcaYZ_eta->Draw("colz");
*/
   auto canvas_dca = new TCanvas("canvas_dca_eta","canvas_dca_eta",1600,1000);
   canvas_dca->Divide(2,2);
   for(int ic = 0; ic < 4; ic++){
      canvas_dca->cd(1+ic);
      gPad->SetMargin(0.15,0.02,0.12,0.1);
   }
   canvas_dca->cd(1);
   if(color>=0){
      h_width_dcaY->SetLineWidth(3);
      h_width_dcaY->SetLineColor(color);
      h_width_dcaY->SetMarkerColor(color);
   }
   h_width_dcaY->GetYaxis()->SetRangeUser(0,100e-4*cm_to_mu);
   gStyle->SetOptStat(0);
   h_width_dcaY->GetYaxis()->SetLabelSize(0.06);
   h_width_dcaY->GetYaxis()->SetTitleSize(0.06);
   h_width_dcaY->GetYaxis()->SetTitleOffset(1.2);
   h_width_dcaY->GetYaxis()->SetNdivisions(505);
   h_width_dcaY->GetXaxis()->SetLabelSize(0.06);
   h_width_dcaY->GetXaxis()->SetTitleSize(0.06);
   h_width_dcaY->GetXaxis()->SetTitleOffset(0.8);
   h_width_dcaY->SetMarkerStyle(markerstyle);
   h_width_dcaY->Draw();
   canvas_dca->cd(3);
   if(color>=0){
      h_mean_dcaY->SetLineWidth(3);
      h_mean_dcaY->SetLineColor(color);
      h_mean_dcaY->SetMarkerColor(color);
   }
   h_mean_dcaY->GetYaxis()->SetRangeUser(-10e-4*cm_to_mu, 10e-4*cm_to_mu);
   gStyle->SetOptStat(0);
   h_mean_dcaY->GetYaxis()->SetLabelSize(0.06);
   h_mean_dcaY->GetYaxis()->SetTitleSize(0.06);
   h_mean_dcaY->GetYaxis()->SetTitleOffset(1.2);
   h_mean_dcaY->GetYaxis()->SetNdivisions(506);
   h_mean_dcaY->GetXaxis()->SetLabelSize(0.06);
   h_mean_dcaY->GetXaxis()->SetTitleSize(0.06);
   h_mean_dcaY->GetXaxis()->SetTitleOffset(0.8);
   h_mean_dcaY->SetMarkerStyle(markerstyle);
   h_mean_dcaY->Draw();

   canvas_dca->cd(2);
   if(color>=0){
      h_width_dcaZ->SetLineWidth(3);
      h_width_dcaZ->SetLineColor(color);
      h_width_dcaZ->SetMarkerColor(color);
   }
   h_width_dcaZ->GetYaxis()->SetRangeUser(0,100e-4*cm_to_mu);
   gStyle->SetOptStat(0);
   h_width_dcaZ->GetYaxis()->SetLabelSize(0.06);
   h_width_dcaZ->GetYaxis()->SetTitleSize(0.06);
   h_width_dcaZ->GetYaxis()->SetTitleOffset(1.2);
   h_width_dcaZ->GetYaxis()->SetNdivisions(505);
   h_width_dcaZ->GetXaxis()->SetLabelSize(0.06);
   h_width_dcaZ->GetXaxis()->SetTitleSize(0.06);
   h_width_dcaZ->GetXaxis()->SetTitleOffset(0.8);
   h_width_dcaZ->SetMarkerStyle(markerstyle);
   h_width_dcaZ->Draw();
   canvas_dca->cd(4);
   if(color>=0){
      h_mean_dcaZ->SetLineWidth(3);
      h_mean_dcaZ->SetLineColor(color);
      h_mean_dcaZ->SetMarkerColor(color);
   }
   h_mean_dcaZ->GetYaxis()->SetRangeUser(-10e-4*cm_to_mu, 10e-4*cm_to_mu);
   gStyle->SetOptStat(0);
   h_mean_dcaZ->GetYaxis()->SetLabelSize(0.06);
   h_mean_dcaZ->GetYaxis()->SetTitleSize(0.06);
   h_mean_dcaZ->GetYaxis()->SetTitleOffset(1.2);
   h_mean_dcaZ->GetYaxis()->SetNdivisions(506);
   h_mean_dcaZ->GetXaxis()->SetLabelSize(0.06);
   h_mean_dcaZ->GetXaxis()->SetTitleSize(0.06);
   h_mean_dcaZ->GetXaxis()->SetTitleOffset(0.8);
   h_mean_dcaZ->SetMarkerStyle(markerstyle);
   h_mean_dcaZ->Draw();

}
