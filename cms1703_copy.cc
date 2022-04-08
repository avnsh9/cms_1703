#include "cms1703_test.h"
#include "cmath"
#include "fastjet/tools/Filter.hh"

// AUTHOR: Avinash Verma
//  EMAIL: avinash.verma@students.iiserpune.ac.in
void Cms1703_test::initialize()
{
  setAnalysisName("cms1703_test");
  setInformation(""
                 "# Search for dark matter produced with an energetic\n"
                 "# jet or a\n"
                 "# √\n"
                 "# hadronically decaying W or Z boson at\n"
                 "");
  setLuminosity(12.9 * units::INVFB);
  bookSignalRegions("MONOJ200_230;MONOJ230_260;MONOJ260_290;MONOJ290_320;MONOJ320_350;MONOJ350_390;MONOJ390_430;MONOJ430_470;MONOJ470_510;MONOJ510_550;MONOJ550_590;MONOJ590_640;MONOJ640_690;MONOJ690_740;MONOJ740_790;MONOJ790_840;MONOJ840_900;MONOJ900_960;MONOJ960_1020;MONOJ1020_1090;MONOJ1090_1160;MONOJ1160;MONOV250_300;MONOV300_350;MONOV350_400;MONOV400_500;MONOV500_600;MONOV600_750;MONOV750");
  // You can also book cutflow regions with bookCutflowRegions("CR1;CR2;..."). Note that the regions are
  //  always ordered alphabetically in the cutflow output files.

  // You should initialize any declared variables here
}

void Cms1703_test::analyze()
{
  // Your eventwise analysis code goes here
  // The following objects are always defined unless they are 'ignored' above. They form std::vector objects of the respective Delphes class type (except for Etmiss which is a single object)
  // All std::vector members and etmiss have the common properties PT, Eta, Phi and P4() with the latter giving access to the full ROOT TLorentzVector.
  // Within a std::vector, all members are ordered with highest pt coming first.

  // electronsLoose, electronsMedium, electronsTight   are list of electrons that passed respective efficiency and reconstruction cuts
  // muonsCombinedPlus, muonsCombined                  as above for muons
  // photonsMedium                                     as above for photons
  // jets are all reconstructed jets                   as above for jets. Note that electrons are most certainly also reconstructed as a jet -> overlap removal do avoid double counting necessary!
  // tracks, towers                                    calorimeter and tracker information. Usually not needed.
  // missingET                                         rec missing ET EXCLUDING muons.

  // Here is a couple of useful functions and lines:
  //------------Phase Space Cuts (defined for jets, electronsXYZ, muonsXYZ, photonsXYZ)
  // jets = filterPhaseSpace(jets, 20., -2.8, 2.8)  // The vector 'jets' only contains jets with pt >= 20 GeV and -2.8 < eta < 2.8. This function is applicable to other particles too (electronsMedium, ... ).
  // jets = overlapRemoval(jets, electronsLoose, 0.2) Removes all jets for which there exists any electron in 'electronsLoose' with deltaR < 0.2.
  // jets = overlapRemovel(jets, 0.2) If two jets overlap within deltaR < 0.2, only the harder jet is stored.

  //------------Isolation Checks (defined for electronsXYZ, muonsXYZ, photonsXYZ
  //------------        For each object, if the user entered N isolation conditions, they can be
  //------------        checked individually be the second argument (running from 0 to N-1).
  // electronsMedium = filterIsolation(electronsMedium, 0)            Removes electrons that do not pass the first isolation condition entered into the AnalysisManager by the user
  // std::vector<int> flags; flags.push_back(0); flags.push_back(2);
  // electronsMedium = filterIsolation(electronsMedium, flags)        Same as above, but both the first and the third condition have to be fulfilled
  // electronsMedium = filterIsolation(electronsMedium)               Same as above, but all conditions have to be fulfilled.

  //-----------Flavour Tag Checks (defined for jets only)
  //----------          Tau tags "loose", "medium" or "tight" can be individually checked if the user enabled tau tagging in the AM.
  //----------          For b-tags, if N working points have been defined, the ith condition can be tested by putting i-1 into the second argument (if there is only one, the argument can be omitted)
  // if checkTauTag(jets[0], "tight") leadingJetIsTagged = True;
  // if checkBTag(jets[0], 0) leadingJetIsBTagged = True;

  //-----------Auxiliary Information
  // - Always ensure that you don't access vectors out of bounds. E.g. 'if(jets[1]->PT > 150)' should rather be if (jets.size() > 1 && jets[1]->PT > 150).
  // - Use rand()/(RAND_MAX+1.) for random numbers between 0 and 1. The random seed is determined from system time or by the RandomSeed parameter in CheckMATE.
  // - The 'return' statement will end this function for the current event and hence should be called whenever the current event is to be vetoed.
  // - Many advanced kinematical functions like mT2 are implemented. Check the manual for more information.
  // - If you need output to be stored in other files than the cutflow/signal files we provide, check the manual for how to do this conveniently.

  // missingET->addMuons(muonsCombined);  // Adds muons to missing ET. This should almost always be done which is why this line is not commented out.
  ++n;

  std::vector<Electron *> electrons_veto;
  std::vector<Muon *> muons_veto;
  std::vector<Photon *> photons_veto;
  std::vector<Jet *> jets_trigger;
  std::vector<Jet *> jets_signal;
  // cout<<"missingET is : "<<missingET->P4().Pt()<<endl;
  // cout<<"missingET is PT : "<<missingET->PT<<endl;
  // std::vector<GenParticle*> true_particles;
  // construct tau
  std::vector<Jet *> taus;
  for (int i = 0; i < jets.size(); i++)
  {
    // if (checkTauTag(jets[i], "tight")) {
    // if (checkTauTag(jets[i], "medium")) {
    if (checkTauTag(jets[i], "loose") and fabs(jets[i]->Charge) == 1)
    {
      taus.push_back(jets[i]);
    }
  }
  taus = filterPhaseSpace(taus, 18., -2.3, 2.3);
  // b-jet construction
  std::vector<Jet *> bjets;
  for (int i = 0; i < jets.size(); i++)
  {
    if (checkBTag(jets[i]))
    {
      bjets.push_back(jets[i]);
    }
  }
  bjets = filterPhaseSpace(bjets, 15., -2.4, 2.4);
  electrons_veto = filterPhaseSpace(electrons, 10., -2.5, 2.5);
  muons_veto = filterPhaseSpace(muons, 10., -2.4, 2.4);
  photons_veto = filterPhaseSpace(photons, 15., -2.5, 2.5);
  jets_trigger = filterPhaseSpace(jets, 20., -5, 5);
  jets_signal = filterPhaseSpace(jets, 20., -2.5, 2.5);

  


  // print muons status
  // if (muons.size()>0){
  //   for (int i = 0; i < muons.size(); i++)
  //   {
  //     cout<<"muon "<<i<<" is : "<<muons[i]->Status<<endl;
  //   }
  // }
  // for (int i = 0; i < muons.size(); i++){

  // }


  std::vector<fastjet::PseudoJet> particles;
  for (int i = 0; i < true_particles.size(); i++)
  {
    if (true_particles[i]->Status == 1){
      int charge=true_particles[i]->Charge;
      // particles.push_back( fastjet::PseudoJet( true_particles[i]->Px, true_particles[i]->Py, true_particles[i]->Pz, true_particles[i]->E ) );
      particles.push_back(fastjet::PseudoJet(true_particles[i]->P4().Px(), true_particles[i]->P4().Py(), true_particles[i]->P4().Pz(), true_particles[i]->P4().E()));
      particles.back().set_user_index(charge);}
      // print true particle charge
      // cout<<"true particle "<<i<<" is : "<<true_particles[i]->Charge<<endl;}
  }
  fastjet::JetDefinition jetAK8(fastjet::antikt_algorithm, 0.4);
  fastjet::ClusterSequence AK8jet(particles, jetAK8);
  std::vector<fastjet::PseudoJet> AK4 = sorted_by_pt(AK8jet.inclusive_jets());

  //filtering AK4
  std::vector<fastjet::PseudoJet> AK4_filtered;
  if (AK4.size()>0){
    for (int i = 0; i < AK4.size(); i++)
    {
      if (AK4[i].pt()>100 && fabs(AK4[i].eta())<2.4){
        AK4_filtered.push_back(AK4[i]);
      }
    }
  }
  // print out the size of ak4 filtered if greater than zero
  if (AK4_filtered.size()>0){
    // cout<<"AK4_filtered size is : "<<AK4_filtered.size()<<endl;
    std::vector<fastjet::PseudoJet> constituents = AK4_filtered[0].constituents();
    // print user index value in AK4_filtered jet constituents
    for (int i = 0; i < constituents.size(); i++)
    {
      cout<<"constituent "<<i<<" is : "<<constituents[i].user_index()<<endl;
    }
    

  }





  


  // cout << "AK8 jets size is : " << AK8.size() << endl;
  // cout << "jet size inbuilt : " << jets.size() << endl;
  // cout << "true particle size inbuilt : " << true_particles.size() << endl;
  // if (jets.size() > 0)
  // {
  //   cout << "jet pt : " << jets[0]->PT << endl;
  // }
  // if (AK8.size() > 0)
  // {
  //   cout << "AK8 jet pt : " << AK8[0].pt() << endl;
  // }
  // for (int i = 0; i < true_particles.size(); i++)
  // {
  //   if (true_particles[i]->PID == 22)
  //   {
  //     cout << "photon particles status is : " << true_particles[i]->Status << endl;
  //   }
  //   else if (true_particles[i]->PID == 11){
  //     cout<<"electron particles status is : "<<true_particles[i]->Status<<endl;
  //   }
  //   else if (true_particles[i]->PID == 13){
  //     cout<<"muon particles status is : "<<true_particles[i]->Status<<endl;
  //   }
  // }

  // fastjet::JetDefinition AKT4(fastjet::antikt_algorithm,0.4);
  // vector<fastjet::PseudoJet> akt4_jets = AKT4(particles);

  // if (AK8.size()>0){
  // cout<<"AK8 jets size is : "<<AK8[0].pt()<<endl;}
  // cout<<"jet size is : "<<jets[0]->PT<<endl;

  // hadronic recoil
  float totmomentumx = 0;
  float totmomentumy = 0;
  float totmomentumz = 0;
  float totenergy = 0;

  float hadronic_recoil = 0;

  if (true_particles.size() > 0)
  {
    for (int i = 0; i < true_particles.size(); i++)
    {
      
      
        totmomentumx += true_particles[i]->P4().Px();
        totmomentumy += true_particles[i]->P4().Py();
        totmomentumz += true_particles[i]->P4().Pz();
        totenergy += true_particles[i]->P4().E();
      
    }
    for (int i = 0; i < true_particles.size(); i++)
    {
      if (((true_particles[i]->PID >= 11 && true_particles[i]->PID <= 18) || true_particles[i]->PID == 22) && true_particles[i]->Status == 1)
      {
        totmomentumx -= true_particles[i]->P4().Px();
        totmomentumy -= true_particles[i]->P4().Py();
        totmomentumz -= true_particles[i]->P4().Pz();
        totenergy -= true_particles[i]->P4().E();
      }
    }
  }
  float hadronic_recoil_temp = totmomentumx * totmomentumx + totmomentumy * totmomentumy;
  hadronic_recoil = sqrtf(hadronic_recoil_temp);


  //sum of leptons and photons transverse momentum for missingET

  // float sum_leptons_ptx = 0;
  // float sum_leptons_pty = 0;
  // if (electrons.size() > 0)
  // {
  //   for (int i = 0; i < electrons.size(); i++)
  //   {
  //     sum_leptons_ptx += electrons[i]->P4().Px();
  //     sum_leptons_pty += electrons[i]->P4().Py();
  //   }
  // }
  // if (photons.size() > 0)
  // {
  //   for (int i = 0; i < muons.size(); i++)
  //   {
  //     sum_leptons_ptx += photons[i]->P4().Px();
  //     sum_leptons_pty += photons[i]->P4().Py();
  //   }
  // }
  
  // float hadx = -missingET->P4().Px() - sum_leptons_ptx;
  // float hady = -missingET->P4().Py() - sum_leptons_pty;
  // float had_recoil = sqrtf(hadx * hadx + hady * hady);
  // cout<<"hadronic recoil is : "<<hadronic_recoil<<endl;
  // cout<<had_recoil<<" "<<hadronic_recoil<<endl;

  // sum of muons pt
  // float sum_muons_pt_using_trprtcl = 0;
  // for (int i = 0; i < true_particles.size(); i++)
  // {
  //   if (true_particles[i]->Status == 1)
  //   {
  //     if (true_particles[i]->PID == 13 || true_particles[i]->PID == -13)
  //     {
  //       sum_muons_pt_using_trprtcl += true_particles[i]->PT;
  //     }
      
  //   }
  // }
  // float sum_muons_pt_using_given_muon = 0;
  // for (int i = 0; i < muons.size(); i++)
  // {
  //   sum_muons_pt_using_given_muon += muons[i]->PT;
  // }
  // cout<<"sum of muons pt using trprtcl is : "<<sum_muons_pt_using_trprtcl<<"using given muon is : "<<sum_muons_pt_using_given_muon<<endl;

  // std::vector<fastjet::PseudoJet> particles;
  // particles.push_back( fastjet::PseudoJet(totmomentumx, totmomentumy, totmomentumz, totenergy) );
  // float momx = 0;
  // float momy = 0;
  // float reco_pt = 0;
  // float elecx = 0;
  // float elecy = 0;
  // float muonx = 0;
  // float muony = 0;
  // float photx = 0;
  // float photy = 0;
  // for (int i = 0; i < true_particles.size(); i++)
  // {
  //   momx = momx + true_particles[i]->P4().Px();
  //   momy = momy + true_particles[i]->P4().Py();
  // }
  // if (electrons.size() > 0)
  // {
  //   for (int i = 0; i < electrons.size(); i++)
  //   {
  //     elecx = elecx + electrons[i]->P4().Px();
  //     elecy = elecy + electrons[i]->P4().Py();
  //   }
  // }
  // if (muons.size() > 0)
  // {
  //   for (int i = 0; i < muons.size(); i++)
  //   {
  //     muonx = muonx + muons[i]->P4().Px();
  //     muony = muony + muons[i]->P4().Py();
  //   }
  // }
  // if (photons.size() > 0)
  // {
  //   for (int i = 0; i < photons.size(); i++)
  //   {
  //     photx = photx + photons[i]->P4().Px();
  //     photy = photy + photons[i]->P4().Py();
  //   }
  // }
  // float lepx = elecx + muonx + photx;
  // float lepy = elecy + muony + photy;
  // float recox = momx - lepx;
  // float recoy = momy - lepy;
  // reco_pt = sqrt((recox * recox) + (recoy * recoy));

  // cout<<reco_pt<<endl;

  // if (photons_veto.size() != 0 || electrons_veto.size() != 0 || taus.size() != 0 || bjets.size() != 0)
  // {

  //   // cout << "event no is :" << n << endl;
  //   // cout << "missingET is :" << missingET->PT << endl;
  //   // cout << "true particles size :" << true_particles.size() << endl;
  //   // cout << "electron size :" << electrons_veto.size() << endl;
  //   // cout << "photons size :" << photons_veto.size() << endl;
  //   // cout << "taus size :" << taus.size() << endl;
  //   // cout << "bjets size :" << bjets.size() << endl;
  // }
  // // cout << "hadronic recoil is :" << hadronic_recoil << endl;
  // if (muons_veto.size() == 2)
  // {
  //   // cout << "muon size :" << muons_veto.size() << endl;
  //   // cout << "muon charge: " << muons_veto[0]->Charge << endl;
  //   // cout << "muon charge: " << muons_veto[1]->Charge << endl;
  //   // cout << "muon charge: " << muons_veto[1]->PT;
  //   // cout << muons_veto[1]->PT << endl;
  //   // cout << "invarien mass is :" << (muons_veto[0]->P4() + muons_veto[1]->P4()).M() << endl;
  //   // cout << "muon charge multiply: " << muons_veto[0]->Charge * muons_veto[1]->Charge << endl;
  //   if (muons_veto[0]->Charge * muons_veto[1]->Charge > 0)
  //     // cout << "yes" << endl;
  // }

  /////veto//////

  if (photons_veto.size() != 0 || electrons_veto.size() != 0 || taus.size() != 0 || bjets.size() != 0)
    return;

  // if (had_recoil < 200.0)
  //   return;

  if (muons_veto.size() != 2)
    return;
  if (muons_veto[0]->Charge * muons_veto[1]->Charge > 0)
    return;

  if (muons_veto[0]->PT < 20.0 && muons_veto[1]->PT < 20.0)
    return;

  if (jets_signal.size() == 0)
    return;

  if (jets_signal[0]->PT < 100.0)
    return;

  //sum of pt of muons

  float muonx = 0;
  float muony = 0;
  for (int i = 0; i < muons_veto.size(); i++)
  {
    muonx = muonx + muons_veto[i]->P4().Px();
    muony = muony + muons_veto[i]->P4().Py();
  }
  float muon_pt=sqrtf((muonx * muonx) + (muony * muony));

  float sum_leptons_ptx = 0;
  float sum_leptons_pty = 0;
  if (electrons.size() > 0)
  {
    for (int i = 0; i < electrons.size(); i++)
    {
      sum_leptons_ptx += electrons[i]->P4().Px();
      sum_leptons_pty += electrons[i]->P4().Py();
    }
  }
  if (photons.size() > 0)
  {
    for (int i = 0; i < muons.size(); i++)
    {
      sum_leptons_ptx += photons[i]->P4().Px();
      sum_leptons_pty += photons[i]->P4().Py();
    }
  }
  
  float hadx = -missingET->P4().Px() - sum_leptons_ptx;
  float hady = -missingET->P4().Py() - sum_leptons_pty;
  float had_recoil = sqrtf(hadx * hadx + hady * hady);

  // cout<<muon_pt<<" "<<hadronic_recoil<<endl;


  // //jets constituents


  // if ((muons_veto[0]->P4() + muons_veto[1]->P4()).M() < 60.0 || (muons_veto[0]->P4() + muons_veto[1]->P4()).M() > 120.0)
  //   return;

  // std::vector<Jet *> jet30 = filterPhaseSpace(jets, 30., -2.5, 2.5);
  // int i=jet30.size();

  // if (i>4){
  //   i=4;
  // }
  // for (int j=0;j<i;j++){
  //   if (fabs(particles[0].phi()-jet30[j]->Phi) < 0.5){
  //     return;
  //   }
  // }
  // cout<<"particles[0].phi() is : "<<particles[0].phi()<<endl;
  // cout<<"jets_signal[0]->Phi is : "<<jets_signal[0]->Phi<<endl;

  // for (int j = 0; j < jets_signal[0]->Constituents.GetEntriesFast(); ++j)
  // {
  //   TObject *object = jets_signal[0]->Constituents.At(j);
  //   if(object == 0) continue;
  //   else if(object->IsA() == GenParticle::Class()){
  //     int pid= ((GenParticle*) object)->PID;
  //     // cout<<"particle pid is :"<<pid<<endl;

  //   }

  // }

  // cout<<"number of charged constituents: "<<jets_signal[0]->NCharged<<endl;
  // cout<<"number of neutral constituents: "<<jets_signal[0]->NNeutrals<<endl;

  //////signal regions //////////
  ///////////////////////////////
  ////////////////////////////////
  if (had_recoil > 200.0 && had_recoil < 230.0)
    countSignalEvent("MONOJ200_230");
  if (had_recoil > 230.0 && had_recoil < 260.0)
    countSignalEvent("MONOJ230_260");
  if (had_recoil > 260.0 && had_recoil < 290.0)
    countSignalEvent("MONOJ260_290");
  if (had_recoil > 290.0 && had_recoil < 320.0)
    countSignalEvent("MONOJ290_320");
  if (had_recoil > 320.0 && had_recoil < 350.0)
    countSignalEvent("MONOJ320_350");
  if (had_recoil > 350.0 && had_recoil < 390.0)
    countSignalEvent("MONOJ350_390");

  if (had_recoil > 390.0 && had_recoil < 430.0)
    countSignalEvent("MONOJ390_430");

  if (had_recoil > 430.0 && had_recoil < 470.0)
    countSignalEvent("MONOJ430_470");

  if (had_recoil > 470.0 && had_recoil < 510.0)
    countSignalEvent("MONOJ470_510");

  if (had_recoil > 510.0 && had_recoil < 550.0)
    countSignalEvent("MONOJ510_550");

  if (had_recoil < 590.0 && had_recoil > 550.0)
    countSignalEvent("MONOJ550_590");

  if (had_recoil < 640.0 && had_recoil > 590.0)
    countSignalEvent("MONOJ590_640");

  if (had_recoil < 690.0 && had_recoil > 640.0)
    countSignalEvent("MONOJ640_690");

  if (had_recoil < 740.0 && had_recoil > 690.0)
    countSignalEvent("MONOJ690_740");

  if (had_recoil < 790.0 && had_recoil > 740.0)
    countSignalEvent("MONOJ740_790");

  if (had_recoil < 840.0 && had_recoil > 790.0)
    countSignalEvent("MONOJ790_840");

  if (had_recoil < 900.0 && had_recoil > 840.0)
    countSignalEvent("MONOJ840_900");

  if (had_recoil < 960.0 && had_recoil > 900.0)
    countSignalEvent("MONOJ900_960");

  if (had_recoil > 960.0 && had_recoil < 1020.0)
    countSignalEvent("MONOJ960_1020");

  if (had_recoil < 1090.0 && had_recoil > 1020.0)
    countSignalEvent("MONOJ1020_1090");

  if (had_recoil > 1090.0 && had_recoil < 1160.0)
    countSignalEvent("MONOJ1090_1160");

  if (had_recoil > 1160.0)
    countSignalEvent("MONOJ1160");
}

void Cms1703_test::finalize()
{
  // Whatever should be done after the run goes here
}
