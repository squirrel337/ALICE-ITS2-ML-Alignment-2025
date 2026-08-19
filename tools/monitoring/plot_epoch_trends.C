#include <fstream>
// ITS2 alignment monitoring, epoch by epoch.
//
// Everything is read from the TrkVtxer tree, i.e. from TrackerFit -- the fit that uses the
// clusters alone. The other tree, ResMonitor, comes from the refit that appends the primary
// vertex as an eighth measured point to drive the alignment gradient; its residuals are pulled
// off the clusters by that constraint and its "DCA" is the residual of a constrained point.
//
// Layout produced:
// s1 and s2 are the two sensor-plane residual coordinates. At the vertex point the code writes
//   Residual_s1 = proj_GZc - meas_GZc                      -> longitudinal
//   Residual_s2 = +-sqrt(dGX^2 + dGY^2)                    -> transverse
// so s2 is the transverse one and pairs with DCA_y, s1 is longitudinal and pairs with DCA_z.
// Transverse is shown first throughout.
//
//   residual_XY_by_epoch.png / _s1   distributions per layer, every epoch overlaid
//   residual_trend_vs_epoch.png          per-layer mean and width against epoch, absolute and relative
//   dca_trend_vs_epoch.png          DCA mean and width against epoch
//   dca_width_vs_pT_by_epoch.png         DCA width against pT, one curve per epoch
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

static const int   NL       = 8;
static const char* LAY[NL]  = {"L0","L1","L2","L3","L4","L5","L6","VTX"};
static const int   LCOL[NL] = {kAzure+2,kAzure-3,kTeal+3,kSpring-6,kOrange-3,kOrange+7,kRed+1,kGray+2};
// half-range per layer, in um: about 4 sigma, so 60 bins sample the core at ~0.15 sigma.
// A single range for all seven layers cannot work -- L0 is 1.85 um wide and L6 is 335.
static const double RNG[NL]  = {10., 45., 60., 500., 600., 1200., 1500., 100.};
static const char* AXT = "#Deltas2 (transverse)";
static const char* AXL = "#Deltas1 (longitudinal)";
static const int   NPT      = 12;                 // 0.5 GeV/c steps up to 6
static const double PTLO = 0.0, PTHI = 6.0;

// epoch colours: cool for the baseline, warm for the last pass
static int ecol(int ie, int ne){
   if (ne < 2) return kBlack;
   double t = double(ie)/(ne-1);
   return TColor::GetColor(float(0.10+0.80*t), float(0.25+0.20*(1-t)), float(0.75-0.65*t));
}

// draw a set of per-epoch or per-layer graphs into the current pad
static void multi(std::vector<TGraphErrors*>& g, const char* xt, const char* yt,
                  const std::vector<TString>& key, bool logy=false, double x0=0.60,
                  bool useErrForRange=true)
{
   gPad->SetGridy(); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14);
   if (logy) gPad->SetLogy();
   TMultiGraph* mg = new TMultiGraph();
   for (auto* gr : g) mg->Add(gr,"lp");
   mg->Draw("A");
   // leave room for the legend, and keep a single wild error bar from setting the scale
   double lo=1e30, hi=-1e30;
   for (auto* gr : g) for (int i=0;i<gr->GetN();++i) {
      double v=gr->GetY()[i], e=useErrForRange?gr->GetEY()[i]:0.;
      if (logy && v<=0) continue;
      lo=std::min(lo,v-e); hi=std::max(hi,v+e);
   }
   if (hi>lo) {
      if (logy) mg->GetYaxis()->SetRangeUser(0.6*lo, 3.0*hi);
      else      { double pad=0.28*(hi-lo); mg->GetYaxis()->SetRangeUser(lo-0.4*pad, hi+2.0*pad); }
   }
   mg->GetXaxis()->SetTitle(xt); mg->GetYaxis()->SetTitle(yt);
   mg->GetYaxis()->SetTitleOffset(1.5);
   mg->GetXaxis()->SetTitleSize(0.05); mg->GetYaxis()->SetTitleSize(0.05);
   mg->GetXaxis()->SetLabelSize(0.045); mg->GetYaxis()->SetLabelSize(0.045);
   TLegend* l = new TLegend(x0,0.62,0.94,0.90);
   l->SetNColumns(g.size()>4?2:1);
   for (size_t i=0;i<g.size();++i) l->AddEntry(g[i], key[i], "lp");
   l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.038); l->Draw();
   gPad->Modified();
}

void plot_epoch_trends(const char* dir="MLPTrain_Step901/Residual", const char* out="tools/monitoring/plots")
{
   gSystem->mkdir(out, kTRUE);
   gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

   // which epochs are on disk
   std::vector<int> EP;
   for (int e = -1; e <= 64; ++e) {
      TString p = Form("%s/Residual_Monitor_Epoch_At_%d_Train.root",dir,e);
      if (!gSystem->AccessPathName(p)) EP.push_back(e);
   }
   const int NE = EP.size();
   if (NE == 0) { printf("no Residual_Monitor_Epoch_At_*_Train.root under %s\n",dir); return; }
   printf("epochs found:"); for(int e:EP) printf(" %d",e); printf("\n");

   // histograms, one set per epoch
   std::vector<std::vector<TH1D*>> h1(NE), h2(NE);          // residuals per layer
   std::vector<TH1D*> hy(NE), hz(NE);                       // DCA integrated
   std::vector<std::vector<TH1D*>> py(NE), pz(NE);          // DCA in pT bins
   for (int ie=0; ie<NE; ++ie) {
      h1[ie].resize(NL); h2[ie].resize(NL);
      for (int l=0;l<NL;l++){
         h1[ie][l]=new TH1D(Form("h1_%d_%d",ie,l),"",60,-RNG[l],RNG[l]);
         h2[ie][l]=new TH1D(Form("h2_%d_%d",ie,l),"",60,-RNG[l],RNG[l]);
      }
      hy[ie]=new TH1D(Form("hy_%d",ie),"",120,-0.03,0.03);
      hz[ie]=new TH1D(Form("hz_%d",ie),"",120,-0.03,0.03);
      py[ie].resize(NPT); pz[ie].resize(NPT);
      for (int b=0;b<NPT;b++){
         py[ie][b]=new TH1D(Form("py_%d_%d",ie,b),"",120,-0.03,0.03);
         pz[ie][b]=new TH1D(Form("pz_%d_%d",ie,b),"",120,-0.03,0.03);
      }
   }

   for (int ie=0; ie<NE; ++ie) {
      for (const char* sfx : {"_Train","_Test"}) {
         TFile* f=TFile::Open(Form("%s/Residual_Monitor_Epoch_At_%d%s.root",dir,EP[ie],sfx));
         if(!f||f->IsZombie()) continue;
         TTree* t=(TTree*)f->Get("TrkVtxer"); if(!t){f->Close();continue;}
         YResidualMonitor* m=new YResidualMonitor(); bool used=false;
         t->SetBranchAddress("trkvtx",&m); t->SetBranchAddress("used",&used);
         for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
            if(!used) continue;
            for(int l=0;l<NL;l++){
               if(l<7 && m->fchipID[l]<0) continue;
               if(m->fds1[l] > -900.) h1[ie][l]->Fill(1e4*m->fds1[l]);
               if(m->fds2[l] > -900.) h2[ie][l]->Fill(1e4*m->fds2[l]);
            }
            hy[ie]->Fill(m->fip[0]); hz[ie]->Fill(m->fip[1]);
            if(m->pT>=PTLO && m->pT<PTHI){
               int b = int((m->pT-PTLO)/(PTHI-PTLO)*NPT);
               py[ie][b]->Fill(m->fip[0]); pz[ie][b]->Fill(m->fip[1]);
            }
         }
         f->Close();
      }
   }

   // ---- fit everything up front -------------------------------------------------
   std::vector<std::vector<Fit>> f1(NE), f2(NE);
   for (int ie=0;ie<NE;++ie){ f1[ie].resize(NL); f2[ie].resize(NL);
      for(int l=0;l<NL;l++){ f1[ie][l]=gfit(h1[ie][l]); f2[ie][l]=gfit(h2[ie][l]); } }
   std::vector<Fit> fy(NE), fz(NE);
   for (int ie=0;ie<NE;++ie){ fy[ie]=gfit(hy[ie]); fz[ie]=gfit(hz[ie]); }

   std::vector<double> xe(NE); for(int ie=0;ie<NE;++ie) xe[ie]=EP[ie];

   // ---- tables ------------------------------------------------------------------
   printf("\n=== RESIDUAL WIDTH by epoch [um] (cluster-only fit) ===\n%-5s","lay");
   for(int ie=0;ie<NE;++ie) printf(" %9s",Form("s2 e%d",EP[ie]));
   for(int ie=0;ie<NE;++ie) printf(" %9s",Form("s1 e%d",EP[ie]));
   printf("     (s2 = transverse, s1 = longitudinal)\n");
   for(int l=0;l<NL-1;l++){
      if(!f1[0][l].ok && !f2[0][l].ok) continue;
      printf("%-5s",LAY[l]);
      for(int ie=0;ie<NE;++ie) printf(" %9.2f",f2[ie][l].sigma);
      for(int ie=0;ie<NE;++ie) printf(" %9.2f",f1[ie][l].sigma);
      printf("\n");
   }
   printf("\n=== RESIDUAL MEAN by epoch [um] ===\n%-5s","lay");
   for(int ie=0;ie<NE;++ie) printf(" %9s",Form("s2 e%d",EP[ie]));
   for(int ie=0;ie<NE;++ie) printf(" %9s",Form("s1 e%d",EP[ie]));
   printf("     (s2 = transverse, s1 = longitudinal)\n");
   for(int l=0;l<NL-1;l++){
      if(!f1[0][l].ok && !f2[0][l].ok) continue;
      printf("%-5s",LAY[l]);
      for(int ie=0;ie<NE;++ie) printf(" %9.2f",f2[ie][l].mean);
      for(int ie=0;ie<NE;++ie) printf(" %9.2f",f1[ie][l].mean);
      printf("\n");
   }
   printf("\n=== change from epoch %d to %d ===\n%-5s %12s %12s %14s %14s\n",
          EP[0],EP[NE-1],"lay","w2 trv [%]","w1 lng [%]","m2 trv [um]","m1 lng [um]");
   for(int l=0;l<NL-1;l++){
      if(!f1[0][l].ok || f1[0][l].sigma<=0) continue;
      printf("%-5s %12.3f %12.3f %+14.2f %+14.2f\n",LAY[l],
             100*(f2[NE-1][l].sigma/f2[0][l].sigma-1),
             100*(f1[NE-1][l].sigma/f1[0][l].sigma-1),
             f2[NE-1][l].mean-f2[0][l].mean, f1[NE-1][l].mean-f1[0][l].mean);
   }
   printf("\n=== DCA by epoch [um] ===\n%8s %10s %10s %10s %10s\n",
          "epoch","mean_y","width_y","mean_z","width_z");
   for(int ie=0;ie<NE;++ie)
      printf("%8d %10.2f %10.2f %10.2f %10.2f\n",EP[ie],
             1e4*fy[ie].mean,1e4*fy[ie].sigma,1e4*fz[ie].mean,1e4*fz[ie].sigma);

   printf("\n=== DCA_y width vs pT by epoch [um] ===\n%12s","pT [GeV/c]");
   for(int ie=0;ie<NE;++ie) printf(" %9s",Form("e%d",EP[ie]));
   printf("\n");
   for(int b=0;b<NPT;b++){
      Fit chk=gfit(py[0][b]); if(!chk.ok) continue;
      printf("%12.2f",PTLO+(b+0.5)*(PTHI-PTLO)/NPT);
      for(int ie=0;ie<NE;++ie){ Fit q=gfit(py[ie][b]); printf(" %9.2f",q.ok?1e4*q.sigma:0.0); }
      printf("\n");
   }

   // ---- residual distributions, epochs overlaid ---------------------------------
   auto overlay=[&](std::vector<std::vector<TH1D*>>& H, std::vector<TH1D*>& D,
                    const char* nm, const char* ttl, const char* dttl){
      TCanvas* c=new TCanvas(nm,"",1500,700); c->Divide(4,2);
      for(int l=0;l<NL;l++){
         if(l==7){   // no vertex point in this fit: show the impact parameter instead
            c->cd(8); gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.15);
            double t2=0; for(int ie=0;ie<NE;++ie) t2=std::max(t2,D[ie]->GetMaximum());
            for(int ie=0;ie<NE;++ie){
               TH1D* h=D[ie];
               h->SetLineColor(ecol(ie,NE)); h->SetLineWidth(ie==NE-1?3:2); h->SetStats(0);
               h->SetMaximum(1.15*t2);
               h->GetXaxis()->SetTitle(Form("%s [cm]",dttl));
               h->GetYaxis()->SetTitle("tracks");
               h->GetXaxis()->SetTitleSize(0.055); h->GetYaxis()->SetTitleSize(0.055);
               h->GetXaxis()->SetLabelSize(0.05);  h->GetYaxis()->SetLabelSize(0.05);
               h->GetYaxis()->SetTitleOffset(1.5); h->GetYaxis()->SetMaxDigits(3);
               h->Draw(ie?"hist same":"hist");
            }
            continue;
         }
         c->cd(l+1); gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.15);
         double top=0; for(int ie=0;ie<NE;++ie) top=std::max(top,H[ie][l]->GetMaximum());
         if(top<=0) continue;
         for(int ie=0;ie<NE;++ie){
            TH1D* h=H[ie][l];
            h->SetLineColor(ecol(ie,NE)); h->SetLineWidth(ie==NE-1?3:2); h->SetStats(0);
            h->SetMaximum(1.15*top);
            h->GetXaxis()->SetTitle(Form("%s  %s [#mum]",LAY[l],ttl));
            h->GetYaxis()->SetTitle("tracks");
            h->GetXaxis()->SetTitleSize(0.055); h->GetYaxis()->SetTitleSize(0.055);
            h->GetXaxis()->SetLabelSize(0.05);  h->GetYaxis()->SetLabelSize(0.05);
            h->GetYaxis()->SetTitleOffset(1.5); h->GetYaxis()->SetMaxDigits(3);
            h->Draw(ie?"hist same":"hist");
         }
         if(l==0){
            TLegend* lg=new TLegend(0.62,0.55,0.95,0.90);
            for(int ie=0;ie<NE;++ie) lg->AddEntry(H[ie][l],Form("epoch %d",EP[ie]),"l");
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.05); lg->Draw();
         }
      }
      c->SaveAs(Form("%s/%s.png",out,nm));
   };
   overlay(h2,hy,"residual_XY_by_epoch",AXT,"DCA_{y} (transverse)");
   overlay(h1,hz,"residual_Z_by_epoch",AXL,"DCA_{z} (longitudinal)");

   // ---- residual mean / width against epoch -------------------------------------
   {
      TCanvas* c=new TCanvas("cre","",1500,900); c->Divide(3,2);
      auto mkset=[&](std::vector<std::vector<Fit>>& F, bool width, bool rel,
                     std::vector<TGraphErrors*>& g, std::vector<TString>& key){
         for(int l=0;l<NL-1;l++){
            if(!F[0][l].ok || (width && F[0][l].sigma<=0)) continue;
            std::vector<double> v(NE), e(NE);
            for(int ie=0;ie<NE;++ie){
               double a = width ? F[ie][l].sigma : F[ie][l].mean;
               double d = width ? F[ie][l].esigma: F[ie][l].emean;
               if(rel){ double b0=F[0][l].sigma; v[ie]=100*(a/b0-1); e[ie]=100*d/b0; }
               else   { v[ie]=a; e[ie]=d; }
            }
            TGraphErrors* gr=new TGraphErrors(NE,&xe[0],&v[0],nullptr,&e[0]);
            gr->SetLineColor(LCOL[l]); gr->SetMarkerColor(LCOL[l]);
            gr->SetMarkerStyle(20+(l%5)); gr->SetMarkerSize(1.1); gr->SetLineWidth(2);
            g.push_back(gr); key.push_back(LAY[l]);
         }
      };
      struct Panel { int pad; bool s2; bool width; bool rel; const char* yt; bool logy; };
      Panel P[6] = {   // transverse on top, longitudinal below
         {1,true, true, false,"#Deltas2 (transverse) width [#mum]",true},
         {2,true, true, true, "#Deltas2 (transverse) width change [%]",false},
         {3,true, false,false,"#Deltas2 (transverse) mean [#mum]",false},
         {4,false,true, false,"#Deltas1 (longitudinal) width [#mum]",true},
         {5,false,true, true, "#Deltas1 (longitudinal) width change [%]",false},
         {6,false,false,false,"#Deltas1 (longitudinal) mean [#mum]",false}};
      for(auto& p : P){
         c->cd(p.pad);
         std::vector<TGraphErrors*> g; std::vector<TString> key;
         mkset(p.s2?f2:f1, p.width, p.rel, g, key);
         if(g.empty()) continue;
         multi(g,"epoch",p.yt,key,p.logy,0.72);
      }
      c->SaveAs(Form("%s/residual_trend_vs_epoch.png",out));
   }

   // ---- DCA against epoch --------------------------------------------------------
   {
      TCanvas* c=new TCanvas("cde","",1300,560); c->Divide(2,1);
      auto g1=[&](std::vector<Fit>& F, bool width, int col, int mst){
         std::vector<double> v(NE), e(NE);
         for(int ie=0;ie<NE;++ie){
            v[ie]=1e4*(width?F[ie].sigma:F[ie].mean);
            e[ie]=1e4*(width?F[ie].esigma:F[ie].emean);
         }
         TGraphErrors* gr=new TGraphErrors(NE,&xe[0],&v[0],nullptr,&e[0]);
         gr->SetLineColor(col); gr->SetMarkerColor(col); gr->SetMarkerStyle(mst);
         gr->SetMarkerSize(1.3); gr->SetLineWidth(3); return gr; };
      c->cd(1);
      { std::vector<TGraphErrors*> g{g1(fy,false,kAzure+2,20),g1(fz,false,kOrange+7,21)};
        std::vector<TString> k{"DCA_{y} transverse","DCA_{z} longitudinal"};
        multi(g,"epoch","DCA mean [#mum]",k); }
      c->cd(2);
      { std::vector<TGraphErrors*> g{g1(fy,true,kAzure+2,20),g1(fz,true,kOrange+7,21)};
        std::vector<TString> k{"DCA_{y} transverse","DCA_{z} longitudinal"};
        multi(g,"epoch","DCA width [#mum]",k); }
      c->SaveAs(Form("%s/dca_trend_vs_epoch.png",out));
   }

   // ---- DCA width vs pT, one curve per epoch -------------------------------------
   {
      TCanvas* c=new TCanvas("cdp","",1300,560); c->Divide(2,1);
      auto panel=[&](std::vector<std::vector<TH1D*>>& H, const char* yt){
         std::vector<TGraphErrors*> g; std::vector<TString> key;
         for(int ie=0;ie<NE;++ie){
            std::vector<double> x,v,e;
            for(int b=0;b<NPT;b++){ Fit q=gfit(H[ie][b]); if(!q.ok) continue;
               x.push_back(PTLO+(b+0.5)*(PTHI-PTLO)/NPT);
               v.push_back(1e4*q.sigma); e.push_back(1e4*q.esigma); }
            if(x.size()<2) continue;
            TGraphErrors* gr=new TGraphErrors(x.size(),&x[0],&v[0],nullptr,&e[0]);
            gr->SetLineColor(ecol(ie,NE)); gr->SetMarkerColor(ecol(ie,NE));
            gr->SetMarkerStyle(20); gr->SetMarkerSize(1.1); gr->SetLineWidth(ie==NE-1?3:2);
            g.push_back(gr); key.push_back(Form("epoch %d",EP[ie]));
         }
         if(!g.empty()) multi(g,"p_{T} [GeV/c]",yt,key,false,0.66);
      };
      c->cd(1); panel(py,"DCA_{y} (transverse) width [#mum]");
      c->cd(2); panel(pz,"DCA_{z} (longitudinal) width [#mum]");
      c->SaveAs(Form("%s/dca_width_vs_pT_by_epoch.png",out));
   }

   printf("\nplots -> %s/{residual_XY_by_epoch,residual_Z_by_epoch,"
          "residual_trend_vs_epoch,dca_trend_vs_epoch,dca_width_vs_pT_by_epoch}.png\n",out);
}
