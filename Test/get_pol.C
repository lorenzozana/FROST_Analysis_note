{
  Float_t x,y,z;
  Int_t nlines = 0;
  TFile *f = new TFile("poltarget.root","RECREATE");
  TTree *PolTargetTree = new TTree("PolTargetTree", "PolTargetTree");
  PolTargetTree->ReadFile("g9a_allruns_polariz.dat", "Run_number/I:Target_pol/F:e_Target_pol_st/F:e_Target_pol_sys/F");
  PolTargetTree->Draw("Target_pol:Run_number","","");
  PolTargetTree->Write();






}
