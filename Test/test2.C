#include "TH1F.h"
#include "TRandom2.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"

double low_tgpol = 0.78;
double high_tgpol = 0.85;

int target_pol = 1; // define target polarization

double mean_beampol = 0.6;
double sigma_beampol = 0.1;


double G_v = 0.5;
double Sigma_v = 0.3;

int tot_events = 10000;

TRandom *fRandom;

double GetTG_pol(int event,int pol) {

  double tg_pol = double(pol) * fRandom->Uniform(low_tgpol,high_tgpol);
  
  return tg_pol;
}


double GetBeam_pol(int event,int pol) {

  double beam_pol = double(pol) * fRandom->Gaus(mean_beampol, sigma_beampol);
  return beam_pol;
}

Double_t dsigma_domega_pol(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = par[0]*(1+par[1]*par[2] * cos(2*xx) + par[3] * par[2] *par[4] * sin(2*xx));  
   // par[0] = cross section unpolarized
   // par[1] = Sigma
   // par[2] = Beam Polarization
   // par[3] = G
   // par[4] = Target Polarization
   return f;
}


void Run(){
  fRandom = new TRandom2(0);
  double beam_pol, tg_pol;
  TFile *file_out = new TFile("output_test_pol_10000.root","RECREATE");

  TF1 *para_func = new TF1("para_func",dsigma_domega_pol,-TMath::Pi(),TMath::Pi(),5);
  TF1 *perp_func = new TF1("para_func",dsigma_domega_pol,-TMath::Pi(),TMath::Pi(),5);

  TH1F *para_h = new TH1F("para_h","PARA histogram; #phi",100,-TMath::Pi(),TMath::Pi());
  TH1F *perp_h = new TH1F("perp_h","PERP histogram; #phi",100,-TMath::Pi(),TMath::Pi());

  TH1F *para2_h = new TH1F("para2_h","PARA histogram weighted pol; #phi",100,-TMath::Pi(),TMath::Pi());
  TH1F *perp2_h = new TH1F("perp2_h","PERP histogram weighted pol; #phi",100,-TMath::Pi(),TMath::Pi());

  TH1F *amo_h = new TH1F("amo_h","AMO histogram; #phi",100,-TMath::Pi(),TMath::Pi());

  double val_phi;

  for (int j=0; j<2; j++) { // beam polarization
    for (int i =0 ; i< tot_events; i++) {
      if (i % 100 == 0) printf("Polarization %d at event %d \n",(j*2-1),i);

      if (i<int(tot_events * 0.8) ) {
	amo_h->Fill(fRandom->Uniform(-TMath::Pi(),TMath::Pi()));
      }
      beam_pol = GetBeam_pol(i,(j*2)-1);
      //      printf("Got Beam Pol \n");
      tg_pol = GetTG_pol(i,target_pol);
      //     printf("Got TG pol \n");
      if (j==0) {
	para_func->SetParameters(10.,Sigma_v,beam_pol,G_v,tg_pol);
	val_phi = para_func->GetRandom();
	para_h->Fill(val_phi);
	para2_h->Fill(val_phi,TMath::Abs(1/beam_pol));
      }

      if (j==1) {
	perp_func->SetParameters(10.,Sigma_v,beam_pol,G_v,tg_pol);
	val_phi = perp_func->GetRandom();
	perp_h->Fill(val_phi);
	perp2_h->Fill(val_phi,TMath::Abs(1/beam_pol));
      }

    }
    
  }
  para_h->Write();
  perp_h->Write();
  para2_h->Write();
  perp2_h->Write();
  amo_h->Write();

  file_out->Close();


}
