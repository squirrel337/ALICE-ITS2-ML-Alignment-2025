#include <fstream>
// Monitoring for the ITS2 alignment test run: residuals, cost, DCA, vertex.
// Mirrors what the 2024 its2-o2-decoupling environment reports, rendered as PNG
// because this ROOT build has no GIF encoder.
#include "../../Ymlp/inc/YMultiLayerPerceptron.h"

const char* LAY[8]={"L0","L1","L2","L3","L4","L5","L6","VTX"};

struct Acc { TH1D *ds1[8],*ds2[8],*chi1[8],*chi2[8]; TH1D *chitot,*ipR,*ipZ,*pipR,*pipZ,*pt;
             TH2D *vtxXY; TH1D *dvx,*dvy,*dvz; long ntr=0; };

void book(Acc&a,const char*tag){
  for(int l=0;l<8;l++){
    a.ds1[l]=new TH1D(Form("ds1_%s_%d",tag,l),"",200,-200,200);
    a.ds2[l]=new TH1D(Form("ds2_%s_%d",tag,l),"",200,-200,200);
    a.chi1[l]=new TH1D(Form("chi1_%s_%d",tag,l),"",200,-10,10);
    a.chi2[l]=new TH1D(Form("chi2_%s_%d",tag,l),"",200,-10,10);
  }
  a.chitot=new TH1D(Form("chitot_%s",tag),"",200,0,150);
  a.ipR=new TH1D(Form("ipR_%s",tag),"",200,-0.05,0.05);
  a.ipZ=new TH1D(Form("ipZ_%s",tag),"",200,-0.05,0.05);
  a.pipR=new TH1D(Form("pipR_%s",tag),"",200,-0.05,0.05);
  a.pipZ=new TH1D(Form("pipZ_%s",tag),"",200,-0.05,0.05);
  a.pt=new TH1D(Form("pt_%s",tag),"",100,0,10);
  a.vtxXY=new TH2D(Form("vtxXY_%s",tag),"",120,-0.15,0.05,120,-0.15,0.05);
  a.dvx=new TH1D(Form("dvx_%s",tag),"",200,-0.05,0.05);
  a.dvy=new TH1D(Form("dvy_%s",tag),"",200,-0.05,0.05);
  a.dvz=new TH1D(Form("dvz_%s",tag),"",200,-0.5,0.5);
}

bool fill(Acc&a,const char* fn){
  TFile* f=TFile::Open(fn); if(!f||f->IsZombie()) return false;
  TTree* t=(TTree*)f->Get("ResMonitor"); if(!t){f->Close();return false;}
  YResidualMonitor* m=new YResidualMonitor(); t->SetBranchAddress("quicklook",&m);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); a.ntr++;
    for(int l=0;l<8;l++){
      if(m->fchipID[l]<0 && l<7) continue;
      if(std::abs(m->fds1[l])<1) a.ds1[l]->Fill(m->fds1[l]*1e4);   // cm -> um
      if(std::abs(m->fds2[l])<1) a.ds2[l]->Fill(m->fds2[l]*1e4);
      a.chi1[l]->Fill(m->chi2ls1[l]); a.chi2[l]->Fill(m->chi2ls2[l]);
    }
    a.chitot->Fill(m->chi2tot);
    a.ipR->Fill(m->fip[0]); a.ipZ->Fill(m->fip[1]);
    a.pipR->Fill(m->p*m->fip[0]); a.pipZ->Fill(m->p*m->fip[1]);
    a.pt->Fill(m->pT);
    a.vtxXY->Fill(m->vtxevtX,m->vtxevtY);
    a.dvx->Fill(m->vtxevtX-m->vtxX); a.dvy->Fill(m->vtxevtY-m->vtxY); a.dvz->Fill(m->vtxevtZ-m->vtxZ);
  }
  f->Close(); return true;
}

void style(TH1* h,int col,const char* xt,const char* yt="tracks"){
  h->SetLineColor(col); h->SetLineWidth(2); h->SetStats(0);
  h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle(yt);
  h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
}

void quicklook(const char* dir="MLPTrain_Step901/Residual", const char* out="tools/monitoring/plots", int FINAL=4){
   gSystem->mkdir(out, kTRUE);
  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.14); gStyle->SetPadBottomMargin(0.14);
  Acc b,e; book(b,"base"); book(e,"ep0");
  bool okb=false, oke=false;
  for(const char* s : {"_Train",""}){
    okb |= fill(b,Form("%s/Residual_Monitor_Epoch_At_-1%s.root",dir,s));
    oke |= fill(e,Form("%s/Residual_Monitor_Epoch_At_%d%s.root",dir,FINAL,s));
  }
  printf("\nbaseline(epoch -1) tracks=%ld   after(epoch %d) tracks=%ld\n",b.ntr,FINAL,e.ntr);
  if(!okb){printf("no quicklook file under %s\n",dir);return;}

  // ---- numeric summary ----
  printf("\n============ RESIDUALS [um] : epoch -1  ->  epoch %d ============\n",FINAL);
  printf("%-5s | %17s %17s | %17s %17s\n","lay","s1 mean","s1 rms","s2 mean","s2 rms");
  for(int l=0;l<8;l++)
    printf("%-5s | %8.2f %8.2f %8.2f %8.2f | %8.2f %8.2f %8.2f %8.2f\n",LAY[l],
           b.ds1[l]->GetMean(),e.ds1[l]->GetMean(),b.ds1[l]->GetRMS(),e.ds1[l]->GetRMS(),
           b.ds2[l]->GetMean(),e.ds2[l]->GetMean(),b.ds2[l]->GetRMS(),e.ds2[l]->GetRMS());
  printf("\n%-22s %12s %12s %10s\n","quantity","epoch -1",Form("epoch %d",FINAL),"delta");
  auto row=[&](const char* n,double x,double y,const char* u){
     printf("%-22s %12.4f %12.4f %+10.4f  %s\n",n,x,y,y-x,u); };
  row("chi2tot mean",b.chitot->GetMean(),e.chitot->GetMean(),"");
  row("DCA_r rms",1e4*b.ipR->GetRMS(),1e4*e.ipR->GetRMS(),"um");
  row("DCA_z rms",1e4*b.ipZ->GetRMS(),1e4*e.ipZ->GetRMS(),"um");
  row("p*DCA_r rms",1e4*b.pipR->GetRMS(),1e4*e.pipR->GetRMS(),"um GeV/c");
  row("vtx dz mean",1e4*b.dvz->GetMean(),1e4*e.dvz->GetMean(),"um");
  printf("\n================ COST / DCA ================\n");
  printf("  chi2tot   mean = %8.3f   median = %8.3f\n",b.chitot->GetMean(),
         [&]{double q,p=0.5;b.chitot->GetQuantiles(1,&q,&p);return q;}());
  printf("  DCA_r     mean = %8.1f um  rms = %8.1f um\n",1e4*b.ipR->GetMean(),1e4*b.ipR->GetRMS());
  printf("  DCA_z     mean = %8.1f um  rms = %8.1f um\n",1e4*b.ipZ->GetMean(),1e4*b.ipZ->GetRMS());
  printf("  p*DCA_r   rms  = %8.1f um*GeV/c  (quality gate uses /40um)\n",1e4*b.pipR->GetRMS());
  printf("  <pT>           = %8.3f GeV/c\n",b.pt->GetMean());
  printf("\n  vertex (estimated - reco):  dx=%+.1f um  dy=%+.1f um  dz=%+.1f um\n",
         1e4*b.dvx->GetMean(),1e4*b.dvy->GetMean(),1e4*b.dvz->GetMean());
  printf("                        rms:  %.1f / %.1f / %.1f um\n",
         1e4*b.dvx->GetRMS(),1e4*b.dvy->GetRMS(),1e4*b.dvz->GetRMS());

  // ---- plots ----
  int C[8]={kAzure+2,kAzure-3,kTeal+3,kOrange+7,kOrange-3,kRed+1,kMagenta+2,kGray+2};
  TCanvas* c1=new TCanvas("c1","residual",1400,700); c1->Divide(4,2);
  for(int l=0;l<8;l++){ c1->cd(l+1); gPad->SetLogy();
    style(b.ds1[l],C[l],Form("%s  #Deltas1 [#mum]",LAY[l]));
    b.ds1[l]->Draw("hist"); e.ds1[l]->SetLineColor(kBlack); e.ds1[l]->SetLineStyle(2); if(oke)e.ds1[l]->Draw("hist same");
  }
  c1->SaveAs(Form("%s/residual_s1.png",out));
  TCanvas* c2=new TCanvas("c2","residual2",1400,700); c2->Divide(4,2);
  for(int l=0;l<8;l++){ c2->cd(l+1); gPad->SetLogy();
    style(b.ds2[l],C[l],Form("%s  #Deltas2 [#mum]",LAY[l]));
    b.ds2[l]->Draw("hist"); e.ds2[l]->SetLineColor(kBlack); e.ds2[l]->SetLineStyle(2); if(oke)e.ds2[l]->Draw("hist same");
  }
  c2->SaveAs(Form("%s/residual_s2.png",out));
  TCanvas* c3=new TCanvas("c3","cost",1400,450); c3->Divide(3,1);
  c3->cd(1); gPad->SetLogy(); style(b.chitot,kAzure+2,"track #chi^{2}_{tot}");
  b.chitot->Draw("hist"); if(oke){e.chitot->SetLineColor(kRed+1);e.chitot->SetLineStyle(2);e.chitot->Draw("hist same");}
  c3->cd(2); gPad->SetLogy(); style(b.chi1[0],kAzure+2,"#chi (s1), all layers");
  for(int l=0;l<7;l++){b.chi1[l]->SetLineColor(C[l]);b.chi1[l]->SetStats(0);b.chi1[l]->Draw(l?"hist same":"hist");}
  c3->cd(3); style(b.pt,kTeal+3,"p_{T} [GeV/c]"); gPad->SetLogy(); b.pt->Draw("hist");
  c3->SaveAs(Form("%s/cost.png",out));
  TCanvas* c4=new TCanvas("c4","dca",1400,450); c4->Divide(3,1);
  c4->cd(1); style(b.ipR,kAzure+2,"DCA_{r} [cm]"); b.ipR->Draw("hist");
  if(oke){e.ipR->SetLineColor(kRed+1);e.ipR->SetLineStyle(2);e.ipR->Draw("hist same");}
  c4->cd(2); style(b.ipZ,kOrange+7,"DCA_{z} [cm]"); b.ipZ->Draw("hist");
  if(oke){e.ipZ->SetLineColor(kRed+1);e.ipZ->SetLineStyle(2);e.ipZ->Draw("hist same");}
  c4->cd(3); style(b.pipR,kTeal+3,"p #upoint DCA_{r} [cm GeV/c]"); b.pipR->Draw("hist");
  c4->SaveAs(Form("%s/dca.png",out));
  TCanvas* c5=new TCanvas("c5","vertex",1400,450); c5->Divide(3,1);
  c5->cd(1); b.vtxXY->SetStats(0); b.vtxXY->GetXaxis()->SetTitle("vtx_{x} [cm]");
  b.vtxXY->GetYaxis()->SetTitle("vtx_{y} [cm]"); b.vtxXY->Draw("colz");
  c5->cd(2); style(b.dvx,kAzure+2,"vtx^{est}_{x} - vtx^{reco}_{x} [cm]"); b.dvx->Draw("hist");
  c5->cd(3); style(b.dvz,kMagenta+2,"vtx^{est}_{z} - vtx^{reco}_{z} [cm]"); b.dvz->Draw("hist");
  c5->SaveAs(Form("%s/vertex.png",out));
  // ---- Cost vs Epoch ----
  {
    std::vector<double> ep,lr,te;
    std::ifstream in("tools/monitoring/cost_epoch.txt");
    double a,b2,c2; while(in>>a>>b2>>c2){ ep.push_back(a); lr.push_back(b2); te.push_back(c2); }
    if(ep.size()>1){
      TCanvas* c6=new TCanvas("c6","costepoch",900,620);
      c6->SetGridy(); gPad->SetLeftMargin(0.15);
      TGraph* g1=new TGraph(ep.size(),&ep[0],&lr[0]);
      TGraph* g2=new TGraph(ep.size(),&ep[0],&te[0]);
      double lo=1e9,hi=-1e9;
      for(double v:lr){lo=std::min(lo,v);hi=std::max(hi,v);} for(double v:te){lo=std::min(lo,v);hi=std::max(hi,v);}
      double pad=0.15*(hi-lo);
      g1->SetTitle(""); g1->GetXaxis()->SetTitle("epoch"); g1->GetYaxis()->SetTitle("cost / ndf");
      g1->GetYaxis()->SetRangeUser(lo-pad,hi+pad);
      g1->GetXaxis()->SetTitleSize(0.05); g1->GetYaxis()->SetTitleSize(0.05);
      g1->GetXaxis()->SetLabelSize(0.045); g1->GetYaxis()->SetLabelSize(0.045);
      g1->GetYaxis()->SetTitleOffset(1.5); g1->GetYaxis()->SetMaxDigits(5);
      g1->SetLineColor(kAzure+2); g1->SetLineWidth(3); g1->SetMarkerColor(kAzure+2);
      g1->SetMarkerStyle(20); g1->SetMarkerSize(1.4);
      g2->SetLineColor(kOrange+7); g2->SetLineWidth(3); g2->SetMarkerColor(kOrange+7);
      g2->SetMarkerStyle(21); g2->SetMarkerSize(1.4);
      g1->Draw("ALP"); g2->Draw("LP same");
      TLegend* lg=new TLegend(0.55,0.74,0.88,0.88);
      lg->AddEntry(g1,"training","lp"); lg->AddEntry(g2,"test","lp");
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042); lg->Draw();
      c6->SaveAs(Form("%s/cost_epoch.png",out));
      printf("\n=== Cost vs Epoch ===\n%8s %14s %14s\n","epoch","training","test");
      for(size_t i=0;i<ep.size();i++) printf("%8.0f %14.6f %14.6f\n",ep[i],lr[i],te[i]);
      printf("  training %+.4f%%   test %+.4f%%  (epoch -1 -> %d)\n",
             100*(lr.back()/lr.front()-1),100*(te.back()/te.front()-1),FINAL);
    }
  }
  printf("\nplots -> %s/{residual_s1,residual_s2,cost,dca,vertex,cost_epoch}.png\n",out);
}
