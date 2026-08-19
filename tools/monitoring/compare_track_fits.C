#include <fstream>
// ITS2 alignment monitoring, taken from the cluster-only track fit.
//
// The module fits each track twice. TrackerFit() uses the clusters alone; its result is
// what TrackVertexQualityEstimator() reports into the TrkVtxer tree. GetCost_Beam_CircleFit()
// then refits with the primary vertex appended as an eighth measured point, because the
// alignment gradient needs that constraint; its result goes into the ResMonitor tree.
// The impact parameter and the residuals must come from the first fit -- in the second the
// vertex is part of the fit, so its "DCA" is a constrained residual pinned to Sigma_MEAS[7].
//
// Everything below therefore reads TrkVtxer. ResMonitor is read only to show the contrast.
#include "../../Ymlp/inc/YMultiLayerPerceptron.h"

struct Fit { double mean, sigma, emean, esigma; bool ok; };

Fit gfit(TH1* h, int maxIter = 8, double tol = 1e-3)
{
   Fit r{0,0,0,0,false};
   if (!h || h->GetEntries() < 50) return r;
   int pb = h->GetMaximumBin();
   double pk = h->GetBinContent(pb), half = 0.5*pk;
   int lo = pb, hi = pb;
   while (lo > 1              && h->GetBinContent(lo) > half) --lo;
   while (hi < h->GetNbinsX() && h->GetBinContent(hi) > half) ++hi;
   double x1 = h->GetBinCenter(lo), x2 = h->GetBinCenter(hi);
   if (!(x2 > x1)) { x1 = h->GetMean()-h->GetRMS(); x2 = h->GetMean()+h->GetRMS(); }
   TF1 f("f","gaus",x1,x2);
   f.SetParameters(pk, h->GetBinCenter(pb), 0.5*(x2-x1)/1.177);
   if (h->Fit(&f,"QNR") != 0) return r;
   double prev = f.GetParameter(2);
   for (int it = 0; it < maxIter; ++it) {
      double m = f.GetParameter(1), s = std::fabs(f.GetParameter(2));
      if (s <= 0) break;
      f.SetRange(m-2*s, m+2*s);
      if (h->Fit(&f,"QNR") != 0) break;
      double now = std::fabs(f.GetParameter(2));
      if (std::fabs(now-prev) < tol*std::max(now,1e-12)) { prev = now; break; }
      prev = now;
   }
   r.mean=f.GetParameter(1); r.sigma=std::fabs(f.GetParameter(2));
   r.emean=f.GetParError(1); r.esigma=f.GetParError(2); r.ok=true;
   return r;
}

struct Diff {
   TH1D** hy; TH1D** hz; int n; double lo, hi; TString name, unit;
   void book(const char* nm, const char* un, int nb, double a, double b, double rng){
      name=nm; unit=un; n=nb; lo=a; hi=b;
      hy=new TH1D*[n]; hz=new TH1D*[n];
      for(int i=0;i<n;i++){
         hy[i]=new TH1D(Form("dy_%s_%d",nm,i),"",120,-rng,rng);
         hz[i]=new TH1D(Form("dz_%s_%d",nm,i),"",120,-rng,rng);
      }
   }
   int bin(double v) const { if(v<lo||v>=hi) return -1; return int((v-lo)/(hi-lo)*n); }
   void fill(double v,double dy,double dz){ int b=bin(v); if(b<0) return; hy[b]->Fill(dy); hz[b]->Fill(dz); }
};

void draw(Diff& d, const char* out)
{
   std::vector<double> x, my, sy, mz, sz, ey, ez, esy, esz;
   for(int i=0;i<d.n;i++){
      Fit fy=gfit(d.hy[i]), fz=gfit(d.hz[i]);
      if(!fy.ok||!fz.ok) continue;
      x.push_back(d.lo+(i+0.5)*(d.hi-d.lo)/d.n);
      my.push_back(1e4*fy.mean); ey.push_back(1e4*fy.emean); sy.push_back(1e4*fy.sigma); esy.push_back(1e4*fy.esigma);
      mz.push_back(1e4*fz.mean); ez.push_back(1e4*fz.emean); sz.push_back(1e4*fz.sigma); esz.push_back(1e4*fz.esigma);
   }
   if(x.size()<2) return;
   TCanvas* c=new TCanvas(Form("c_%s",d.name.Data()),"",1200,520); c->Divide(2,1);
   auto mk=[&](std::vector<double>& v,std::vector<double>& e,int col,int mst){
      TGraphErrors* g=new TGraphErrors(x.size(),&x[0],&v[0],nullptr,&e[0]);
      g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mst);
      g->SetMarkerSize(1.2); g->SetLineWidth(2); return g; };
   c->cd(1); gPad->SetGridy(); gPad->SetLeftMargin(0.15);
   TGraphErrors* g1=mk(my,ey,kAzure+2,20); TGraphErrors* g2=mk(mz,ez,kOrange+7,21);
   TMultiGraph* m1=new TMultiGraph(); m1->Add(g1,"lp"); m1->Add(g2,"lp"); m1->Draw("A");
   m1->GetXaxis()->SetTitle(Form("%s %s",d.name.Data(),d.unit.Data()));
   m1->GetYaxis()->SetTitle("DCA mean [#mum]"); m1->GetYaxis()->SetTitleOffset(1.5);
   TLegend* l1=new TLegend(0.55,0.76,0.9,0.9); l1->AddEntry(g1,"DCA_{y} transverse","lp");
   l1->AddEntry(g2,"DCA_{z} longitudinal","lp"); l1->SetBorderSize(0); l1->SetFillStyle(0); l1->Draw();
   c->cd(2); gPad->SetGridy(); gPad->SetLeftMargin(0.15);
   TGraphErrors* g3=mk(sy,esy,kAzure+2,20); TGraphErrors* g4=mk(sz,esz,kOrange+7,21);
   TMultiGraph* m2=new TMultiGraph(); m2->Add(g3,"lp"); m2->Add(g4,"lp"); m2->Draw("A");
   m2->GetXaxis()->SetTitle(Form("%s %s",d.name.Data(),d.unit.Data()));
   m2->GetYaxis()->SetTitle("DCA width [#mum]"); m2->GetYaxis()->SetTitleOffset(1.5);
   m2->SetMinimum(0.);
   TLegend* l2=new TLegend(0.55,0.76,0.9,0.9); l2->AddEntry(g3,"DCA_{y} transverse","lp");
   l2->AddEntry(g4,"DCA_{z} longitudinal","lp"); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->Draw();
   c->SaveAs(Form("%s/dca_clusterfit_vs_%s.png",out,d.name.Data()));
   printf("\n  %-6s %10s %10s %10s %10s\n",d.name.Data(),"mean_y","width_y","mean_z","width_z");
   for(size_t i=0;i<x.size();i++)
      printf("  %6.2f %10.2f %10.2f %10.2f %10.2f\n",x[i],my[i],sy[i],mz[i],sz[i]);
}

void compare_track_fits(const char* dir="MLPTrain_Step901/Residual", const char* out="tools/monitoring/plots", int FINAL=4)
{
   gSystem->mkdir(out, kTRUE);
   gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   gStyle->SetPadLeftMargin(0.15); gStyle->SetPadBottomMargin(0.14);

   Diff dpt, det, dz, dph;
   dpt.book("pT",  "[GeV/c]", 12, 0.0,  6.0, 0.03);   // 0.5 GeV/c per bin
   det.book("eta", "",        20, -1.0, 1.0, 0.03);
   dz .book("zvtx","[cm]",    60, -15.0,15.0,0.03);
   dph.book("phi", "[rad]",   12, 0.0,  2*TMath::Pi(), 0.03);

   // residuals per layer, from the cluster-only fit; index 7 is the vertex point
   TH1D *t1[8],*t2[8],*t1f[8],*t2f[8];      // TrkVtxer  : cluster-only fit
   TH1D *r1[8],*r2[8];                      // ResMonitor: fit including the vertex, for contrast
   // one range per layer, about 4 sigma: a single range cannot serve L0 (2 um) and L6 (415 um).
   // A wider window than the core clips the fit; a much wider one under-samples it.
   const double RNG[8] = {10., 45., 60., 500., 600., 1200., 1500., 100.};
   for(int l=0;l<8;l++){
      double rng = RNG[l];
      t1[l] =new TH1D(Form("t1_%d",l), "",60,-rng,rng);
      t2[l] =new TH1D(Form("t2_%d",l), "",60,-rng,rng);
      t1f[l]=new TH1D(Form("t1f_%d",l),"",60,-rng,rng);
      t2f[l]=new TH1D(Form("t2f_%d",l),"",60,-rng,rng);
      r1[l] =new TH1D(Form("q1_%d",l), "",60,-rng,rng);
      r2[l] =new TH1D(Form("q2_%d",l), "",60,-rng,rng);
   }
   TH1D* hipy=new TH1D("hipy","",120,-0.03,0.03), *hipz=new TH1D("hipz","",120,-0.03,0.03);
   TH1D* hdvz=new TH1D("hdvz","",120,-0.4,0.1);

   auto scanTrk=[&](int epoch, bool base){
      for(const char* sfx : {"_Train","_Test"}){
         TFile* f=TFile::Open(Form("%s/Residual_Monitor_Epoch_At_%d%s.root",dir,epoch,sfx));
         if(!f||f->IsZombie()) continue;
         TTree* t=(TTree*)f->Get("TrkVtxer"); if(!t){f->Close();continue;}
         YResidualMonitor* m=new YResidualMonitor(); bool used=false;
         t->SetBranchAddress("trkvtx",&m); t->SetBranchAddress("used",&used);
         for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
            if(!used) continue;                       // tracks the vertexer kept
            for(int l=0;l<8;l++){
               if(l<7 && m->fchipID[l]<0) continue;
               if(m->fds1[l] > -900.) (base?t1[l]:t1f[l])->Fill(1e4*m->fds1[l]);
               if(m->fds2[l] > -900.) (base?t2[l]:t2f[l])->Fill(1e4*m->fds2[l]);
            }
            if(!base) continue;
            double dy=m->fip[0], dzz=m->fip[1];
            hipy->Fill(dy); hipz->Fill(dzz);
            double phi=m->phi; if(phi<0) phi+=2*TMath::Pi();
            dpt.fill(m->pT,dy,dzz);  det.fill(m->eta,dy,dzz);
            dz .fill(m->vtxZ,dy,dzz);dph.fill(phi,dy,dzz);
            hdvz->Fill(m->vtxevtZ-m->vtxZ);
         }
         f->Close();
      }
   };
   auto scanRes=[&](int epoch){
      for(const char* sfx : {"_Train","_Test"}){
         TFile* f=TFile::Open(Form("%s/Residual_Monitor_Epoch_At_%d%s.root",dir,epoch,sfx));
         if(!f||f->IsZombie()) continue;
         TTree* t=(TTree*)f->Get("ResMonitor"); if(!t){f->Close();continue;}
         YResidualMonitor* m=new YResidualMonitor(); t->SetBranchAddress("monitor",&m);
         for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
            for(int l=0;l<8;l++){
               if(l<7 && m->fchipID[l]<0) continue;
               if(m->fds1[l] > -900.) r1[l]->Fill(1e4*m->fds1[l]);
               if(m->fds2[l] > -900.) r2[l]->Fill(1e4*m->fds2[l]);
            }
         }
         f->Close();
      }
   };
   scanTrk(-1,true); scanTrk(FINAL,false); scanRes(-1);

   const char* LAY[8]={"L0","L1","L2","L3","L4","L5","L6","VTX"};
   // s2 is the transverse residual, s1 the longitudinal one: at the vertex point the module
   // writes Residual_s1 = dGZ and Residual_s2 = +-sqrt(dGX^2+dGY^2). Transverse comes first.
   printf("\n=== RESIDUALS from the cluster-only fit (TrkVtxer) [um] ===\n");
   printf("%-5s | %9s %9s %9s %9s | %9s %9s %9s %9s\n","lay",
          "m2(-1)",Form("m2(%d)",FINAL),"w2(-1)",Form("w2(%d)",FINAL),
          "m1(-1)",Form("m1(%d)",FINAL),"w1(-1)",Form("w1(%d)",FINAL));
   printf("        <---------- transverse ----------> <--------- longitudinal --------->\n");
   for(int l=0;l<8;l++){
      Fit a=gfit(t2[l]),b=gfit(t2f[l]),c=gfit(t1[l]),d=gfit(t1f[l]);
      printf("%-5s | %9.2f %9.2f %9.2f %9.2f | %9.2f %9.2f %9.2f %9.2f\n",LAY[l],
             a.mean,b.mean,a.sigma,b.sigma,c.mean,d.mean,c.sigma,d.sigma);
   }
   printf("\n=== width comparison at epoch -1: which fit the residual came from [um] ===\n");
   printf("%-5s | %11s %11s | %11s %11s\n","lay",
          "trv cluster","trv +vertex","lng cluster","lng +vertex");
   for(int l=0;l<8;l++){
      Fit a=gfit(t2[l]),b=gfit(r2[l]),c=gfit(t1[l]),d=gfit(r1[l]);
      printf("%-5s | %11.2f %11.2f | %11.2f %11.2f\n",LAY[l],a.sigma,b.sigma,c.sigma,d.sigma);
   }
   Fit fy=gfit(hipy), fz=gfit(hipz), fv=gfit(hdvz);
   printf("\n=== DCA from the cluster-only fit, against the chi2 vertex ===\n");
   printf("  DCA_y mean %+.2f um  width %.2f um\n  DCA_z mean %+.2f um  width %.2f um\n",
          1e4*fy.mean,1e4*fy.sigma,1e4*fz.mean,1e4*fz.sigma);
   printf("\n=== dvz = v_est,z - v_reco,z ===\n  fit mean %+.1f um   width %.1f um\n",
          1e4*fv.mean,1e4*fv.sigma);

   auto rescan=[&](TH1D** A, TH1D** B, const char* nm, const char* ttl){
      TCanvas* c=new TCanvas(nm,"",1300,620); c->Divide(4,2);
      for(int l=0;l<8;l++){ c->cd(l+1); gPad->SetLeftMargin(0.17);
         A[l]->SetLineColor(kAzure+2); A[l]->SetLineWidth(2); A[l]->SetStats(0);
         A[l]->GetXaxis()->SetTitle(Form("%s  %s [#mum]",LAY[l],ttl));
         A[l]->GetYaxis()->SetTitle("tracks"); A[l]->Draw("hist");
         B[l]->SetLineColor(kRed+1); B[l]->SetLineStyle(2); B[l]->SetLineWidth(2); B[l]->Draw("hist same");
         Fit ff=gfit(A[l]);
         if(ff.ok){ TF1* g=new TF1(Form("g_%s_%d",nm,l),"gaus",ff.mean-2*ff.sigma,ff.mean+2*ff.sigma);
            g->SetParameters(A[l]->GetMaximum(),ff.mean,ff.sigma); g->SetLineColor(kBlack); g->SetLineWidth(1); g->Draw("same"); }
      }
      c->SaveAs(Form("%s/%s.png",out,nm));
   };
   rescan(t2,t2f,"residual_XY_distributions","#Deltas2 (transverse)");
   rescan(t1,t1f,"residual_Z_distributions","#Deltas1 (longitudinal)");

   draw(dpt,out); draw(det,out); draw(dz,out); draw(dph,out);
   printf("\nplots -> %s/{residual_Z_distributions,residual_XY_distributions,dca_clusterfit_vs_pT,dca_clusterfit_vs_eta,dca_clusterfit_vs_zvtx,dca_clusterfit_vs_phi}.png\n",out);
}
