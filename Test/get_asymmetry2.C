{
  double low_tgpol = 0.78;
  double high_tgpol = 0.85;
  double avg_tgpol = (low_tgpol + high_tgpol) /2;

  int target_pol = 1; // define target polarization
 
  double mean_beampol = 0.6;
  double sigma_beampol = 0.1;


  TFile *file_in = TFile::Open("output_test_pol_1000.root");
  
  TH1F *para_h = file_in->Get("para_h");
  TH1F *perp_h = file_in->Get("perp_h");
  
  TH1F *para2_h = file_in->Get("para2_h");
  TH1F *perp2_h = file_in->Get("perp2_h");
  
  para_h->Rebin(2);
  perp_h->Rebin(2);
  para2_h->Rebin(2);
  perp2_h->Rebin(2);

  TCanvas *c1 = new TCanvas();
  TH1F *asym = perp_h->GetAsymmetry(para_h);
  TH1F *asym2 = perp2_h->GetAsymmetry(para2_h);
  TF1 *cos2phi_sin2phi = new TF1("cos2phi_sin2phi","([0]+([1]*cos((2*(x-[2]))))+([3]*sin((2*(x-[2])))))",-TMath::Pi(),TMath::Pi());
  cos2phi_sin2phi->FixParameter(2,0);
  asym->Fit("cos2phi_sin2phi","QB");
  asym->Draw("E");
  double G_pol=((cos2phi_sin2phi->GetParameter(3))/(avg_tgpol * mean_beampol));
  double G_pol_e = ((cos2phi_sin2phi->GetParError(3))/(avg_tgpol * mean_beampol));
  printf ("G= %.3f +/- %.3f \n",G_pol,G_pol_e);

  TCanvas *c2 = new TCanvas();

  asym2->Fit("cos2phi_sin2phi","QB");
  asym2->Draw("E");
  double G_pol2=((cos2phi_sin2phi->GetParameter(3))/(avg_tgpol));
  double G_pol2_e = ((cos2phi_sin2phi->GetParError(3))/(avg_tgpol));
  printf ("weighted polarization G= %.3f +/- %.3f \n",G_pol2,G_pol2_e);


}
