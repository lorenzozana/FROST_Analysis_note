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
  TFile *file_out = new TFile("output_test_10000.root","RECREATE");

  TF1 *para_func = new TF1("para_func",dsigma_domega_pol,-TMath::Pi(),TMath::Pi(),5);
  TF1 *perp_func = new TF1("para_func",dsigma_domega_pol,-TMath::Pi(),TMath::Pi(),5);

  TH1F *para_h = new TH1F("para_h","PARA histogram; #phi",100,-TMath::Pi(),TMath::Pi());
  TH1F *perp_h = new TH1F("perp_h","PERP histogram; #phi",100,-TMath::Pi(),TMath::Pi());


  for (int j=0; j<2; j++) { // beam polarization
    for (int i =0 ; i< tot_events; i++) {
      if (i % 100 == 0) printf("Polarization %d at event %d \n",(j*2-1),i);

      beam_pol = GetBeam_pol(i,(j*2)-1);
      //      printf("Got Beam Pol \n");
      tg_pol = GetTG_pol(i,target_pol);
      //     printf("Got TG pol \n");
      if (j==0) {
	para_func->SetParameters(10.,Sigma_v,beam_pol,G_v,tg_pol);
	para_h->Fill(para_func->GetRandom());
      }

      if (j==1) {
	perp_func->SetParameters(10.,Sigma_v,beam_pol,G_v,tg_pol);
	perp_h->Fill(perp_func->GetRandom());
      }

    }
    
  }
  para_h->Write();
  perp_h->Write();

  file_out->Close();


}
