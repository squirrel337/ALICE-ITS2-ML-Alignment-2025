#define cm_to_mu 1.0e+4
void check_vertex_phi_plots_trkvtxer_color(int color=-1, int markerstyle=21, double scale = 1.0){
   new TCanvas();
   const int nbinsx = 120;

   TH2D* h_event_Vx_phi = new TH2D("h_event_Vx_phi","Vx(cm) vs phi(rad);phi(rad);x(cm)", nbinsx, -TMath::Pi(), +TMath::Pi(), 2000,-0.1,0.1);
   TH2D* h_event_Vy_phi = new TH2D("h_event_Vy_phi","Vy(cm) vs phi(rad);phi(rad);y(cm)", nbinsx, -TMath::Pi(), +TMath::Pi(), 2000,-0.1,0.1);
   TH2D* h_event_Vz_phi = new TH2D("h_event_Vz_phi","Vz(cm) vs phi(rad);phi(rad);z(cm)", nbinsx, -TMath::Pi(), +TMath::Pi(), 1500,-15,15);

   TrkVtxer->Draw("vtxevtX:TMath::ASin(fyP[2])+fyAlpha>>h_event_Vx_phi");
   TrkVtxer->Draw("vtxevtY:TMath::ASin(fyP[2])+fyAlpha>>h_event_Vy_phi");
   TrkVtxer->Draw("vtxevtZ:TMath::ASin(fyP[2])+fyAlpha>>h_event_Vz_phi");

   TH2D* hdcaY_phi   = new TH2D("hdcaY_phi",
                                "dcaY vs phi(rad);phi(rad);dcaY(#mum)",
                                nbinsx, -TMath::Pi(), +TMath::Pi(), 500, -0.05*cm_to_mu, 0.05*cm_to_mu);
   TProfile* pf_width_dcaY = new TProfile("pf_width_dcaY","pf_width_dcaY",nbinsx,-TMath::Pi(),TMath::Pi(),-0.02*cm_to_mu,0.02*cm_to_mu,"s");
   TProfile* pf_widthE_dcaY = new TProfile("pf_widthE_dcaY","pf_widthE_dcaY",nbinsx,-TMath::Pi(),TMath::Pi(),-0.02*cm_to_mu,0.02*cm_to_mu,"g");
   TF1* fgaus_dcaY = new TF1("fgaus_dcaY","gaus",-0.01*cm_to_mu,0.01*cm_to_mu);

   TH2D* hdcaZ_phi   = new TH2D("hdcaZ_phi",
                                "dcaZ vs phi(rad);phi(rad);dcaZ(#mum)", 
                                nbinsx, -TMath::Pi(), +TMath::Pi(), 500, -0.05*cm_to_mu, 0.05*cm_to_mu);
   TH2D* hdcaYZ_phi  = new TH2D("hdcaYZ_phi",
                                "#sqrt{dcaY^{2}+dcaZ^{2}} vs phi(rad);phi(rad);#sqrt{dcaY^{2}+dcaZ^{2}}(#mum)", 
                                nbinsx, -TMath::Pi(), +TMath::Pi(),1000,0*cm_to_mu,0.1*cm_to_mu);
   TProfile* pf_width_dcaZ = new TProfile("pf_width_dcaZ","pf_width_dcaZ",nbinsx,-TMath::Pi(),TMath::Pi(),-0.02*cm_to_mu,0.02*cm_to_mu,"s");
   TProfile* pf_widthE_dcaZ = new TProfile("pf_widthE_dcaZ","pf_widthE_dcaZ",nbinsx,-TMath::Pi(),TMath::Pi(),-0.02*cm_to_mu,0.02*cm_to_mu,"g");
   TF1* fgaus_dcaZ = new TF1("fgaus_dcaZ","gaus",-0.01*cm_to_mu,0.01*cm_to_mu);

   TrkVtxer->Draw(Form("fip[0]*%g:TMath::ASin(fyP[2])+fyAlpha>>hdcaY_phi",cm_to_mu),"");
   TrkVtxer->Draw(Form("fip[0]*%g:TMath::ASin(fyP[2])+fyAlpha>>pf_width_dcaY",cm_to_mu),"");
   TrkVtxer->Draw(Form("fip[0]*%g:TMath::ASin(fyP[2])+fyAlpha>>pf_widthE_dcaY",cm_to_mu),"");

   TrkVtxer->Draw(Form("fip[1]*%g:TMath::ASin(fyP[2])+fyAlpha>>hdcaZ_phi",cm_to_mu),"");
   TrkVtxer->Draw(Form("fip[1]*%g:TMath::ASin(fyP[2])+fyAlpha>>pf_width_dcaZ",cm_to_mu),"");
   TrkVtxer->Draw(Form("fip[1]*%g:TMath::ASin(fyP[2])+fyAlpha>>pf_widthE_dcaZ",cm_to_mu),"");

   TrkVtxer->Draw(Form("TMath::Sqrt(fip[0]*fip[0]+fip[1]*fip[1])*%g:TMath::ASin(fyP[2])+fyAlpha>>hdcaYZ_phi",cm_to_mu),"");

   //fitter
   TH1D* h_mean_dcaY = new TH1D("h_mean_dcaY","(mean) dcaY vs #phi;#phi(rad);dcaY(#mum)",nbinsx,-TMath::Pi(),TMath::Pi());
   TH1D* h_width_dcaY = new TH1D("h_width_dcaY","(width) dcaY vs #phi;#phi(rad);#sigma(#mum)",nbinsx,-TMath::Pi(),TMath::Pi());
   for(int ibin = 0; ibin < nbinsx; ibin++){
      fgaus_dcaY->SetParameter(0,0);
      fgaus_dcaY->SetParameter(1,0);
      fgaus_dcaY->SetParameter(2,0);
      fgaus_dcaY->SetParError(0,0);
      fgaus_dcaY->SetParError(1,0);
      fgaus_dcaY->SetParError(2,0);
      auto proj_hdcaY_phi = (TH1D*) hdcaY_phi->ProjectionY(Form("proj_hdcaY_phi_%d",ibin),1+ibin,1+ibin);

      double minX = proj_hdcaY_phi->GetBinCenter(proj_hdcaY_phi->FindFirstBinAbove(0.5*proj_hdcaY_phi->GetMaximum()));
      double maxX = proj_hdcaY_phi->GetBinCenter(proj_hdcaY_phi->FindLastBinAbove(0.5*proj_hdcaY_phi->GetMaximum()));
      double rangeX = 5*0.5*(std::abs(minX) + std::abs(maxX));
      proj_hdcaY_phi->Fit(fgaus_dcaY,"rQ","",-rangeX,+rangeX);
      proj_hdcaY_phi->Fit(fgaus_dcaY,"rQ","",-rangeX,+rangeX);

      double parA[3] = {fgaus_dcaY->GetParameter(0),fgaus_dcaY->GetParameter(1),fgaus_dcaY->GetParameter(2)};
      double errA[3] = {fgaus_dcaY->GetParError(0),fgaus_dcaY->GetParError(1),fgaus_dcaY->GetParError(2)};

      proj_hdcaY_phi->Fit(fgaus_dcaY,"rQ","", fgaus_dcaY->GetParameter(1) - scale*fgaus_dcaY->GetParameter(2), fgaus_dcaY->GetParameter(1) + scale*fgaus_dcaY->GetParameter(2));

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

   TH1D* h_mean_dcaZ = new TH1D("h_mean_dcaZ","(mean) dcaZ vs #phi;#phi(rad);dcaZ(#mum)",nbinsx,-TMath::Pi(),TMath::Pi());
   TH1D* h_width_dcaZ = new TH1D("h_width_dcaZ","(width) dcaZ vs #phi;#phi(rad);#sigma(#mum)",nbinsx,-TMath::Pi(),TMath::Pi());
   for(int ibin = 0; ibin < nbinsx; ibin++){
      fgaus_dcaZ->SetParameter(0,0);
      fgaus_dcaZ->SetParameter(1,0);
      fgaus_dcaZ->SetParameter(2,0);
      fgaus_dcaZ->SetParError(0,0);
      fgaus_dcaZ->SetParError(1,0);
      fgaus_dcaZ->SetParError(2,0);
      auto proj_hdcaZ_phi = (TH1D*) hdcaZ_phi->ProjectionY(Form("proj_hdcaZ_phi_%d",ibin),1+ibin,1+ibin);

      double minX = proj_hdcaZ_phi->GetBinCenter(proj_hdcaZ_phi->FindFirstBinAbove(0.5*proj_hdcaZ_phi->GetMaximum()));
      double maxX = proj_hdcaZ_phi->GetBinCenter(proj_hdcaZ_phi->FindLastBinAbove(0.5*proj_hdcaZ_phi->GetMaximum()));
      double rangeX = 5*0.5*(std::abs(minX) + std::abs(maxX));
      proj_hdcaZ_phi->Fit(fgaus_dcaZ,"rQ","",-rangeX,+rangeX);
      proj_hdcaZ_phi->Fit(fgaus_dcaZ,"rQ","",-rangeX,+rangeX);

      double parA[3] = {fgaus_dcaZ->GetParameter(0),fgaus_dcaZ->GetParameter(1),fgaus_dcaZ->GetParameter(2)};
      double errA[3] = {fgaus_dcaZ->GetParError(0),fgaus_dcaZ->GetParError(1),fgaus_dcaZ->GetParError(2)};

      proj_hdcaZ_phi->Fit(fgaus_dcaZ,"rQ","", fgaus_dcaZ->GetParameter(1) - scale*fgaus_dcaZ->GetParameter(2), fgaus_dcaZ->GetParameter(1) + scale*fgaus_dcaZ->GetParameter(2));

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
   h_event_Vx_phi->Draw("colz");

   canvas_vertex->cd(2);
   h_event_Vy_phi->Draw("colz");

   canvas_vertex->cd(3);
   h_event_Vz_phi->Draw("colz");

   canvas_vertex->cd(4);
   hdcaY_phi->Draw("colz");

   canvas_vertex->cd(5);
   hdcaZ_phi->Draw("colz");

   canvas_vertex->cd(6);
   hdcaYZ_phi->Draw("colz");
*/

   auto canvas_dca = new TCanvas("canvas_dca_phi","canvas_dca_phi",1600,1000);
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
   h_mean_dcaY->GetYaxis()->SetRangeUser(-30e-4*cm_to_mu, 30e-4*cm_to_mu);
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
   h_mean_dcaZ->GetYaxis()->SetRangeUser(-30e-4*cm_to_mu, 30e-4*cm_to_mu);
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
