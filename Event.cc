//
//  Event Class
//
//  This is the class in which all the hard work is done, objects of other 
//  classes are contained here and pre-selections, selections and calculations
//  are done by the member functions of this class.
//


#define Event_cxx

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "Event.h"
#include "AnalysisMuonEfficiencyScaleFactors.h" 
#include "MuonTriggerEfficiency.h"
#include "MuonIsolationEfficiency.h"
#include "PtReweighting.h"
//using namespace std;


//--- Constructors
Event::Event(){

  m_muon_n  = 0;
  m_vxp_n   = 0;
  m_trigger = 0;

  m_k_factor_weight           = 1;
  m_lumi_weight               = 1; 
  m_vertex_weight             = 1; 
  m_mcevent_weight            = 1;
  m_mc_muon_eff_weight        = 1;
  m_pileup_weight             = 1;
  m_trigger_efficiency_weight = 1;
  m_isolation_efficiency_weight = 1;
  m_vertex_position_weight      = 1;
  
  m_kfactor_mode  = 0;
  m_mc_sim_period = 0;

  m_make_Z_selection  = false;
  m_make_DY_selection = false;

  m_lar_error = 0;

  m_muon1_reco_sf = 0;
  m_muon2_reco_sf = 0;

  m_muon1_reco_sf_error = 0;
  m_muon2_reco_sf_error = 0;

  m_mc_pt_reweight = 0;


}


Event::Event(int muons_n, int vxp_n, bool trigger){

  m_muon_n  = muons_n;
  m_vxp_n   = vxp_n;
  m_trigger = trigger;
 
  m_k_factor_weight             = 1;
  m_lumi_weight                 = 1; 
  m_vertex_weight               = 1; 
  m_mcevent_weight              = 1;
  m_mc_muon_eff_weight          = 1;
  m_pileup_weight               = 1;
  m_trigger_efficiency_weight   = 1;
  m_isolation_efficiency_weight = 1;
  m_vertex_position_weight      = 1;

  m_kfactor_mode  = 0;
  m_mc_sim_period = 0;
  
  m_make_Z_selection  = false;
  m_make_DY_selection = false;

  m_lar_error = 0;

  m_muon1_reco_sf = 0;
  m_muon2_reco_sf = 0;

  m_muon1_reco_sf_error = 0;
  m_muon2_reco_sf_error = 0;
}


//--- Destructor

Event::~Event() {}

bool Event::IsVerbose() const { return m_is_verbose; }


///--- Getter Functions

int   Event::KFactorMode()              const { return m_kfactor_mode;              }
int   Event::NumberOfMuons()            const { return m_muon_n;                    }
int   Event::NumberOfPrimaryVertices()  const { return m_vxp_n;                     }
int   Event::NumberOfTrks()             const { return m_trk_n;                     }
bool  Event::TriggerFired()             const { return m_trigger;                   }
int   Event::EventNumber()              const { return m_event_number;              }
int   Event::RunNumber()                const { return m_run_number;                }
int   Event::MCSimPeriod()              const { return m_mc_sim_period;             }
float Event::AverageIntPerXing()        const { return m_av_int_per_xing;           }

int   Event::NumberOfTruthMuons()       const { return m_mc_muons_truth.size();     }
int   Event::NumberOfMCBlockParticles() const { return m_mc_block_particles.size(); }

unsigned int Event::LArError()          const { return m_lar_error;                 }

Particle       Event::Muon1()        const { return m_muon1;        }
Particle       Event::Muon2()        const { return m_muon2;        }
Particle       Event::Muon()         const { return m_muon;         }
Particle       Event::AntiMuon()     const { return m_antimuon;     }
SimpleParticle Event::Propagator()   const { return m_propagator;   }
Vertex         Event::ZerothVertex() const { return m_zerothvertex; }

std::vector<MCMuonTruth> Event::MCMuonsTruth() const { return m_mc_muons_truth;     }
std::vector<MCBlock> Event::MCBlockParticles() const { return m_mc_block_particles; }


double Event::KFactorWeight() const { 
  if( UseKFactorWeighting == false ){
    return 1.0;
  }else{
    return m_k_factor_weight;
  }
}

double Event::LumiWeight()    const { 
  if( UseLumiWeighting == false){
    return 1.0;
  }else{
    return m_lumi_weight;   
  }
}

double Event::MCEventWeight() const { 
  if(UseMCEventWeighting == false ){
    return 1;
  }else{
    return m_mcevent_weight; }
}

double Event::VertexWeight()  const { 
  if( UseVertexWeighting == false){
    return 1.0;
  }else{
    return m_vertex_weight;
   } 
}

double Event::MuonMCEfficiencyWeight() const {
  if( UseMuonEfficiencyWeighting == false){
    return 1.0;
  }else{
    return m_mc_muon_eff_weight;
  }
}

double Event::Muon1RecoSf() const {
  return m_muon1_reco_sf;
}

double Event::Muon2RecoSf() const {
  return m_muon2_reco_sf;
}

double Event::Muon1RecoSfError() const {
  return m_muon1_reco_sf_error;
}

double Event::Muon2RecoSfError() const {
  return m_muon2_reco_sf_error;
}

double Event::PileupWeight() const {
  if( UsePileupWeighting == false ){
    return 1.0;
  }else{
    return m_pileup_weight;
  }
}

double Event::TriggerEfficiencyWeight() const {
  if( UseTriggerEfficiencyWeighting == false ){ 
    return 1.0; 
  }else{ 
    return m_trigger_efficiency_weight; 
  }
}

double Event::IsolationEfficiencyWeight() const {
  if( UseIsolationEfficiencyWeighting == false ){ 
    return 1.0; 
  }else{ 
    return m_isolation_efficiency_weight; 
  }
}

double Event::VertexPositionWeight() const {
  if( UseVertexPositionWeighting == false ){
    return 1.0;
  }else{
    return m_vertex_position_weight;
  }

}

double Event::PtReweight() const {
  
  if( UsePtReweighting == false ){
    return 1.0;
  }else{
    return m_mc_pt_reweight;
  }

}

double Event::EventWeight()   const {
  if( IsMCEvent == false ){
    return 1;
  }else{
    return this->KFactorWeight() * this->LumiWeight() * this->VertexWeight() * this->MCEventWeight() * this->MuonMCEfficiencyWeight() * this->PileupWeight() * this->TriggerEfficiencyWeight() * this->IsolationEfficiencyWeight() * this->VertexPositionWeight() * this->PtReweight();
  }
}


void Event::SetIsVerbose( const bool verbose )
{
  m_is_verbose = verbose; 
}

void Event::PrintWeights() const {
  
  std::cout<<"----------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<"WEIGHTS FOR EVENT "<<this->EventNumber()<<":"<<std::endl;

  std::cout<<""<<std::endl;

  std::cout<<"\t Total:      "<<"\t\t"<<                                  this->EventWeight()             <<std::endl;
  std::cout<<"\t MCEvent:    "<< UseMCEventWeighting             <<"\t"<< this->MCEventWeight()    <<       std::endl;
  std::cout<<"\t kFactor:    "<< UseKFactorWeighting             <<"\t"<< this->KFactorWeight()  <<       std::endl;
  std::cout<<"\t Lumi:       "<< UseLumiWeighting                <<"\t"<< this->LumiWeight()<<       std::endl;
  std::cout<<"\t Vertex:     "<< UseVertexWeighting              <<"\t"<< this->VertexWeight()  <<       std::endl;
  std::cout<<"\t Reco:       "<< UseMuonEfficiencyWeighting      <<"\t"<< this->MuonMCEfficiencyWeight()  <<std::endl;
  std::cout<<"\t Pileup:     "<< UsePileupWeighting              <<"\t"<< this->PileupWeight()            <<std::endl;
  std::cout<<"\t Trigger:    "<< UseTriggerEfficiencyWeighting   <<"\t"<< this->TriggerEfficiencyWeight() <<std::endl;
  std::cout<<"\t Vertex:     "<< UseVertexPositionWeighting      <<"\t"<< this->VertexPositionWeight()    <<std::endl;   
  std::cout<<"\t Iso:        "<< UseIsolationEfficiencyWeighting <<"\t"<< m_isolation_efficiency_weight   <<std::endl;
  std::cout<<"\t PtReweight: "<<  UsePtReweighting               <<"\t"<< this->PtReweight()              <<std::endl;

  std::cout<<"----------------------------------------------------------------------------------------------"<<std::endl;
}


void Event::PrintEventInfo() const{

  std::cout<<"----------------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<"EVENT INFO FOR "<<this->EventNumber()<<":"<<std::endl;
  std::cout<<"nMuons = "<< this->NumberOfMuons() <<endl;
  std::cout<<"Muon1 = ";
  this->Muon1().Print();
  std::cout<<"Pt = "<< this->Muon1().Pt() <<endl; 
  std::cout<<"Muon2 = ";
  this->Muon2().Print();
  std::cout<<"Pt = "<< this->Muon2().Pt() <<endl;   
  std::cout<<"----------------------------------------------------------------------------------------------"<<std::endl;


}


void Event::CheckSystematicVariations(){

  int trueCount = 0;

  for(std::map< TString, bool>::iterator iter = m_syst_variations.begin(); iter != m_syst_variations.end(); iter++)
    {
      if( iter->second == true ){ trueCount++ ; }      
    }
  
  if( trueCount > 1 ){
    cout<<" Event::CheckSystematicVariations() Error: Too many systematic variations turned on."<<endl;
    cout<<" Aborting "<<endl;
    abort();
  }

  return;
}

std::map<TString, bool> Event::SystematicVariations() const {
  return m_syst_variations;
}

void Event::SetSystematicVariations(const std::map< TString, bool>& vars ){
  m_syst_variations = vars;
}

void  Event::SetKFactorMode( const int mode ){
  m_kfactor_mode = mode;
}



float  Event::MCPdfScale() const { return m_mc_pdf_scale; }
float  Event::MCPdf1()     const { return m_mc_pdf1;      }
float  Event::MCPdf2()     const { return m_mc_pdf2;      }
float  Event::MCPdfX1()    const { return m_mc_pdf_x1;    }
float  Event::MCPdfX2()    const { return m_mc_pdf_x2;    }
int    Event::MCPdfId1()   const { return m_mc_pdf_id1;   }
int    Event::MCPdfId2()   const { return m_mc_pdf_id2;   }


std::vector<Vertex>   Event::Vertices()       const { return m_vertices;        }
std::vector<Trk>      Event::Trks()           const { return m_trks;            }
std::vector<Particle> Event::Muons()          const { return m_muons;           }
std::vector<Particle> Event::MuonCandidates() const { return m_muon_candidates; }

std::vector<double> Event::PrimaryVertexDistributionWeights() const { return m_primary_vertex_distribution_weights; }

bool Event::MakeZSelection(){  return m_make_Z_selection;  }
bool Event::MakeDYSelection(){ return m_make_DY_selection; }


double Event::MuonDeltaPhi() const { 

  double delta_phi = this->Muon1().LorentzVector().DeltaPhi(this->Muon2().LorentzVector());
  return delta_phi; 
}


double Event::MuonDeltaEta() const {

  double etaA;
  double etaB;

  if( this->Muon1().Eta() > this->Muon2().Eta() ){
    etaA = this->Muon1().Eta();
    etaB = this->Muon2().Eta();
  }else{
    etaA = this->Muon2().Eta();
    etaB = this->Muon1().Eta();
  }

  double delta_eta = etaA - etaB;

  return delta_eta;
}


double Event::MuonDeltaR() const {

  //  double deltaR = TMath::Sqrt( pow( this->MuonDeltaEta() , 2) + pow( this->MuonDeltaPhi() , 2) );
  
  double deltaR =  this->Muon1().LorentzVector().DeltaR(this->Muon2().LorentzVector());
  return deltaR;
}

double Event::MuonIDDeltaR() const {

  double deltaR =  this->Muon1().IDVector3().DeltaR(this->Muon2().IDVector3());
  return deltaR;
}

double Event::MuonDeltaZ0() const{
  
  float z1;
  float z2;
  if( this->Muon1().Z0() > this->Muon2().Z0() ){
    z1 = this->Muon1().Z0();
    z2 = this->Muon2().Z0();
  }else{    
    z1 = this->Muon2().Z0();
    z2 = this->Muon1().Z0();
  }

  float delta_z = z1 - z2;

  return delta_z;
}



//-- This function sets Muon1(), .Muon2()
void Event::TwoHighestPtMuonsRegardlessOfCharge( const std::vector<Particle> candidate_muons){

  
  Particle muon1(0.0, 0.0, 0.0, 0.0, 0.0);
  Particle muon2(0.0, 0.0, 0.0, 0.0, 0.0);
  
  for( int i=0; i<static_cast<int>(candidate_muons.size()); i++)
    {
     
      if( candidate_muons.at(i).Pt() > muon1.Pt() ){
	
	muon2 = muon1; 
	muon1 = candidate_muons.at(i);
	
      }else if( candidate_muons.at(i).Pt() > muon2.Pt() ){
	muon2 = candidate_muons.at(i);
      }
 
    } //end of for loop
  //  cout<<"++++++++++++++++++++++++++++++++++++++"<<endl;

  this->SetMuon1( muon1 );
  this->SetMuon2( muon2 );    

  if( this->IsVerbose() ){
    cout<<"Muon1 pt: "<< this->Muon1().Pt() <<" "<< muon1.Pt() <<endl;
    cout<<"Muon2 pt: "<< this->Muon2().Pt() <<" "<< muon2.Pt() <<endl;
  }
     

} //end of func





bool Event::EventHasPositiveAndNegativeMuons( const std::vector<Particle> candidate_muons ){
  
  bool eventHasPositiveAndNegativeMuons = false;

  int positive_count = 0;
  int negative_count = 0;
  
  for( int i=0; i<static_cast<int>(candidate_muons.size()); i++)
    {
      if( candidate_muons.at(i).Charge() > 0 ){positive_count++;}
      if( candidate_muons.at(i).Charge() < 0 ){negative_count++;}
    }

  if( positive_count !=0 && negative_count != 0){
    eventHasPositiveAndNegativeMuons = true;
  }
       
  return eventHasPositiveAndNegativeMuons;

}




bool Event::BothMuonsPassEtRatioCut( const int cone, const double et_ratio_cut  ){

  bool event_passes = false;
  if( cone == 20 ){    
    if( (this->Muon1().EtRatio20() < et_ratio_cut) && (this->Muon2().EtRatio20() < et_ratio_cut) ){ event_passes = true; }
  }
  if( cone == 30 ){    
    if( (this->Muon1().EtRatio30() < et_ratio_cut) && (this->Muon2().EtRatio30() < et_ratio_cut) ){ event_passes = true; }
  }  
  if( cone == 40 ){    
    if( (this->Muon1().EtRatio40() < et_ratio_cut) && (this->Muon2().EtRatio40() < et_ratio_cut) ){ event_passes = true; }
  }

  return event_passes;
}


bool Event::BothMuonsPassPtRatioCut( const int cone, const double pt_ratio_cut  ){

  bool event_passes = false;

  if( cone == 20 ){    
    if( (this->Muon1().PtRatio20() < pt_ratio_cut) && (this->Muon2().PtRatio20() < pt_ratio_cut) ){ event_passes = true; }
  }
  if( cone == 30 ){    
    if( (this->Muon1().PtRatio30() < pt_ratio_cut) && (this->Muon2().PtRatio30() < pt_ratio_cut) ){ event_passes = true; }
  }  
  if( cone == 40 ){    
    if( (this->Muon1().PtRatio40() < pt_ratio_cut) && (this->Muon2().PtRatio40() < pt_ratio_cut) ){ event_passes = true; }
  }
 

  return event_passes;
}

bool Event::BothMuonsFailEtRatioCut( const int cone, const double et_ratio_cut  ){

  bool event_passes = false;
  if( cone == 20 ){    
    if( (this->Muon1().EtRatio20() >= et_ratio_cut) && (this->Muon2().EtRatio20() >= et_ratio_cut) ){ event_passes = true; }
  }
  if( cone == 30 ){    
    if( (this->Muon1().EtRatio30() >= et_ratio_cut) && (this->Muon2().EtRatio30() >= et_ratio_cut) ){ event_passes = true; }
  }  
  if( cone == 40 ){    
    if( (this->Muon1().EtRatio40() >= et_ratio_cut) && (this->Muon2().EtRatio40() >= et_ratio_cut) ){ event_passes = true; }
  }

  return event_passes;
}


bool Event::BothMuonsFailPtRatioCut( const int cone, const double pt_ratio_cut  ){

  bool event_passes = false;

  if( cone == 20 ){    
    if( (this->Muon1().PtRatio20() >= pt_ratio_cut) && (this->Muon2().PtRatio20() >= pt_ratio_cut) ){ event_passes = true; }
  }
  if( cone == 30 ){    
    if( (this->Muon1().PtRatio30() >= pt_ratio_cut) && (this->Muon2().PtRatio30() >= pt_ratio_cut) ){ event_passes = true; }
  }  
  if( cone == 40 ){    
   if( (this->Muon1().PtRatio40() >= pt_ratio_cut) && (this->Muon2().PtRatio40() >= pt_ratio_cut) ){ event_passes = true; }
  }
 

  return event_passes;
}


bool Event::BothMuonsPassEtaCut( const double eta_cut  ){

  bool event_passes = false;

  if( (fabs(this->Muon1().Eta()) < eta_cut) && (fabs(this->Muon2().Eta()) < eta_cut) ){ event_passes = true; }

  return event_passes;
}

bool Event::BothMuonsPassPtCut( const double pt_cut ){
  bool event_passes = false;
  if( (this->Muon1().Pt() > pt_cut) && (this->Muon2().Pt() > pt_cut) ){ event_passes = true; }
  return event_passes;
}

bool Event::BothMuonsPassPtCut( const double pt1_cut, const double pt2_cut ){
  bool event_passes = false;
  if( (this->Muon1().Pt() > pt1_cut) && (this->Muon2().Pt() > pt2_cut) ){ event_passes = true; }
  return event_passes;
}




bool Event::BothMuonsPassD0Cut( const double d0_cut  ){

  bool event_passes = false;

  if( (fabs(this->Muon1().D0()) < d0_cut) && (fabs(this->Muon2().D0()) < d0_cut) ){ event_passes = true; }

  return event_passes;
}

bool Event::BothMuonsPassD0SignificanceCut( const double d0_sig_cut  ){

  bool event_passes = false;

  if( (fabs(this->Muon1().IDD0Significance()) < d0_sig_cut) && (fabs(this->Muon2().IDD0Significance()) < d0_sig_cut) ){ event_passes = true; }

  return event_passes;
}



bool Event::BothMuonsPassZ0Cut( const double z0_cut  ){

  bool event_passes = false;

  if( (fabs(this->Muon1().Z0()) < z0_cut) && (fabs(this->Muon2().Z0()) < z0_cut) ){ event_passes = true; }
 
  return event_passes;
}


bool Event::PropagatorPassesRapidityCut( const double rap_cut  ){

  bool event_passes = false;

  if( fabs(this->Propagator().Rapidity()) < rap_cut ){ event_passes = true; }

  return event_passes;
}


bool Event::BothMuonsPassPtSumCut( const int cone, const double pt_cut ){

 bool event_passes = false;

  if( cone == 20 ){    
    if( (this->Muon1().PtCone20() < pt_cut) && (this->Muon2().PtCone20() < pt_cut) ){ event_passes = true; }
  }
  if( cone == 30 ){    
    if( (this->Muon1().PtCone30() < pt_cut) && (this->Muon2().PtCone30() < pt_cut) ){ event_passes = true; }
  }  
  if( cone == 40 ){    
    if( (this->Muon1().PtCone40() < pt_cut) && (this->Muon2().PtCone40() < pt_cut) ){ event_passes = true; }
  }

  return event_passes;

}

bool Event::BothMuonsPassEtSumCut( const int cone, const double et_cut ){

 bool event_passes = false;

  if( cone == 20 ){    
    if( (this->Muon1().EtCone20() < et_cut) && (this->Muon2().EtCone20() < et_cut) ){ event_passes = true; }
  }
  if( cone == 30 ){    
    if( (this->Muon1().EtCone30() < et_cut) && (this->Muon2().EtCone30() < et_cut) ){ event_passes = true; }
  }  
  if( cone == 40 ){    
    if( (this->Muon1().EtCone40() < et_cut) && (this->Muon2().EtCone40() < et_cut) ){ event_passes = true; }
  }

  return event_passes;

}

bool Event::BothMuonsFailPtSumCut( const int cone, const double pt_cut ){

 bool event_passes = false;

  if( cone == 20 ){    
    if( (this->Muon1().PtCone20() >= pt_cut) && (this->Muon2().PtCone20() >= pt_cut) ){ event_passes = true; }
  }
  if( cone == 30 ){    
    if( (this->Muon1().PtCone30() >= pt_cut) && (this->Muon2().PtCone30() >= pt_cut) ){ event_passes = true; }
  }  
  if( cone == 40 ){    
    if( (this->Muon1().PtCone40() >= pt_cut) && (this->Muon2().PtCone40() >= pt_cut) ){ event_passes = true; }
  }

  return event_passes;

}


bool Event::BothMuonsFailEtSumCut( const int cone, const double et_cut ){

 bool event_passes = false;

  if( cone == 20 ){    
    if( (this->Muon1().EtCone20() >= et_cut) && (this->Muon2().EtCone20() >= et_cut) ){ event_passes = true; }
  }
  if( cone == 30 ){    
    if( (this->Muon1().EtCone30() >= et_cut) && (this->Muon2().EtCone30() >= et_cut) ){ event_passes = true; }
  }  
  if( cone == 40 ){    
    if( (this->Muon1().EtCone40() >= et_cut) && (this->Muon2().EtCone40() >= et_cut) ){ event_passes = true; }
  }

  return event_passes;

}



bool Event::MuonsPassDeltaZ0Cut( const double maxDeltaZ0 ){

  bool event_passes = false;

  //if( fabs(this->MuonDeltaZ0()) < maxDeltaZ0 ){ event_passes = true; }

  /// Elisa uses this method:
  if( fabs( this->Muon1().Z0() - this->Muon2().Z0()) < maxDeltaZ0 ) { event_passes = true; }
 
  return event_passes;
}


///--- Setter Functions

void Event::SetEventNumber( const int evt ){
  m_event_number = evt;
}

void Event::SetRunNumber( const int run ){
  m_run_number = run;
}


void Event::SetMCSimPeriod( const int period ){
  m_mc_sim_period = period;
}

void Event::SetAverageIntPerXing( const float mu ){
  m_av_int_per_xing = mu;
}


void Event::SetMakeDYSelection( const bool DY_sel ){

  m_make_DY_selection = DY_sel;

  if(m_make_DY_selection == true){
    m_make_Z_selection = false;
  }else{
    m_make_Z_selection = true;
  }

}

void Event::SetMakeZSelection(const bool Z_sel){

  m_make_Z_selection = Z_sel;

  if(m_make_Z_selection == true){
    m_make_DY_selection = false;
  }else{
    m_make_DY_selection = true;
  }

}



void Event::SetNumberOfMuons( const int num_muons ){
  m_muon_n = num_muons;
}

void Event::SetNumberOfPrimaryVertices( const int num_pvtx ){
  m_vxp_n = num_pvtx;
}

void Event::SetNumberOfTrks( const int num_trks ){
  m_trk_n = num_trks;
}


void Event::SetIfTriggerFired( const bool trigger ){
  m_trigger = trigger;
} 

void Event::SetLArError( const unsigned int error ){
  m_lar_error = error;
}


void Event::SetMuon1( const Particle muon1 ){
  m_muon1 = muon1;
}

void Event::SetMuon2( const Particle muon2 ){
  m_muon2 = muon2;
}

void Event::SetMuon( const Particle muon ){
  m_muon = muon;
}

void Event::SetAntiMuon( const Particle antimuon ){
  m_antimuon = antimuon;
}


void Event::SetMuonAndAntiMuon(){

  if( this->Muon1().Charge() < 0 ){
    this->SetMuon( this->Muon1() );
    this->SetAntiMuon( this->Muon2() );
  }else{
    this->SetMuon( this->Muon2() );
    this->SetAntiMuon( this->Muon1() );
  }

}


void Event::SetPropagator( const SimpleParticle propagator){
  m_propagator = propagator;
}


void Event::SetDYMuonCandidatesRegardlessOfCharge( const std::vector<Particle>& muon_candidates ){
  m_muon_candidates_regardless_of_charge = muon_candidates;
  this->TwoHighestPtMuonsRegardlessOfCharge(muon_candidates);
}

void Event::SetMuons( const std::vector<Particle> &mu ){
  m_muons = mu;
}

void Event::SetMCMuonsTruth( const std::vector<MCMuonTruth>& muons ){
  m_mc_muons_truth = muons;
}

void Event::SetMCBlockParticles( const std::vector<MCBlock>& particles ){
  m_mc_block_particles = particles;
} 


void Event::SetZerothVertex()
{
  m_zerothvertex = this->Vertices().at(0);
}


void  Event::SetVertices( const std::vector<Vertex>& vertices ){
  m_vertices = vertices;
}

void Event::SetTrks( const std::vector<Trk> trks ){
  m_trks = trks;
}


void Event::CorrectForOverLappingCones(){


//- the tmpMuon*() Particles are uses as Event::Muon1() Event::Muon2() are consts

//--- Look at Cone20 Variables -----
  if( this->MuonIDDeltaR() < 0.2 ){
    
    float adjustedMuon1PtCone20 = this->Muon1().PtCone20() - this->Muon2().IDVector3().Pt();
    float adjustedMuon2PtCone20 = this->Muon2().PtCone20() - this->Muon1().IDVector3().Pt();

    Particle tmpMuon1 = this->Muon1();
    tmpMuon1.SetPtCone20( adjustedMuon1PtCone20 ); 
    this->SetMuon1( tmpMuon1 );
    
    Particle tmpMuon2 = this->Muon2();
    tmpMuon2.SetPtCone20( adjustedMuon2PtCone20 ); 
    this->SetMuon2( tmpMuon2 );

  }

//--- Look at Cone30 Variables -----
  if( this->MuonIDDeltaR() < 0.3 ){
    
    float adjustedMuon1PtCone30 = this->Muon1().PtCone30() - this->Muon2().IDVector3().Pt();
    float adjustedMuon2PtCone30 = this->Muon2().PtCone30() - this->Muon1().IDVector3().Pt();

    Particle tmpMuon1 = this->Muon1();
    tmpMuon1.SetPtCone30( adjustedMuon1PtCone30 ); 
    this->SetMuon1( tmpMuon1 );
      
    Particle tmpMuon2 = this->Muon2();
    tmpMuon2.SetPtCone30( adjustedMuon2PtCone30 ); 
    this->SetMuon2( tmpMuon2 );
  }


//--- Look at Cone40 Variables -----
  if( this->MuonIDDeltaR() < 0.4 ){
    
    float adjustedMuon1PtCone40 = this->Muon1().PtCone40() - this->Muon2().IDVector3().Pt();
    float adjustedMuon2PtCone40 = this->Muon2().PtCone40() - this->Muon1().IDVector3().Pt();
    

    Particle tmpMuon1 = this->Muon1();
    tmpMuon1.SetPtCone40( adjustedMuon1PtCone40 ); 
    this->SetMuon1( tmpMuon1 );
        
    Particle tmpMuon2 = this->Muon2();
    tmpMuon2.SetPtCone40( adjustedMuon2PtCone40 ); 
    this->SetMuon2( tmpMuon2 );

  }


}


void Event::CheckVertexTrackIndexRelation(){

//-- This is an error checking fucntion. Call it if you fear that then indices between vertices and trks have become out of sync
//-- The pt_sum calcualated should be equal (or very close to) Vertex::SumPt(). 


  for( int nvtx=0; nvtx<this->NumberOfPrimaryVertices(); nvtx++)
    {
      
      float pt_sum =0;
      
      for( unsigned int i=0; i<this->Vertices().at(nvtx).VertexTrkIndex().size(); i++)
	{
	  int index =  this->Vertices().at(nvtx).VertexTrkIndex().at(i);
	  
	  pt_sum += this->Trks().at(index).Pt();	  
	}	

      
//-- Grep for "Mismatch" in the log file 
      if (this->IsVerbose() ){
	cout<<"Vertex: "<< nvtx <<endl;
	cout <<" Final pt_sum = "<< pt_sum <<endl;
	cout <<" vxp_sumPt = "<< this->Vertices().at(nvtx).SumPt() <<endl;
	cout <<" vxp_sumPt - pt_sum = "<< this->Vertices().at(nvtx).SumPt() - pt_sum <<endl;
	if( fabs( this->Vertices().at(nvtx).SumPt() - pt_sum) > 0.1 ) { cout <<"Mismatch "<<endl; }
	cout<<"---"<<endl;
      }
      
    }//end of vertex loop                    


}



void Event::SetKFactorWeight(const int sampleID){

  bool drellYanSample = false;
  if( (sampleID == 113700) || (sampleID ==113701) || (sampleID == 108319) || (sampleID == 108321 ) ){ 
    this->SetKFactorMode(0);
    drellYanSample = true; 
  }
  if( (sampleID == 113711 ) || (sampleID == 113712 ) ){
    this->SetKFactorMode(1);
    drellYanSample = true; 
  }

  if( (UseKFactorWeighting == false) || (drellYanSample == false) ){
    m_k_factor_weight = 1;
  }else{

    double mass     = this->MCPdfScale();
    double rapidity = TMath::Log( (this->MCPdfX1()*7000)/mass ) ;
    

    //    std::cout<<"Rap/M: "<<rapidity<<"\t"<<mass<<std::endl;

    m_k_factor_weight = this->kGridfn( mass, rapidity, this->KFactorMode() );
  }

}

void Event::SetLumiWeight(const double lumi, const double nEvts, const double crossSection){

  if (UseLumiWeighting == false ){
    m_lumi_weight = 1;
  }else{


    double weight = ( (lumi * crossSection )/nEvts );

    m_lumi_weight = weight;


  }


}




void Event::SetLumiWeight(const double lumi, const double nEvts){

  if (UseLumiWeighting == false ){
    m_lumi_weight = 1;
  }else{

//-- Get Monte Carlo Cross-Section from MCSamples.txt
//-- The luminosity of the data and number of MC events comes from the pileuop reweighting tool
//-- Add any Monte Carlo Samples to the file not listed already.
    
    int    a; 
    double b, c;
    double crossSection = -1;
    double filter_eff   = -1;

    std::string filename = "MC11Samples.txt"; 
 
    std::ifstream myfile; 
 
    myfile.open ( filename.c_str() );
  
    if ( myfile.fail() ){
      std::cout<< "Unable to open MC11Samples.txt"<<std::endl;
    }else{
	while ( !myfile.eof() )
	  {
	
	    myfile >> a >> b >> c ;
	  
	    if( a == this->RunNumber() ){
	      crossSection = b;
	      filter_eff   = c;
	      break;
	    }
	
	  }
	myfile.close();
    }
    
    if ( crossSection == -1 ){ 
      std::cout << "Monte Carlo Sample not in MCSamples.txt" <<std::endl;
    }
    

    double weight = ( (lumi * crossSection )/nEvts )*filter_eff;

    m_lumi_weight = weight;

  }


}


void Event::SetVertexWeight( const TH1D* h_vertex_weights ){

  if( UseVertexWeighting == false ){
    m_vertex_weight = 1;
  }else{
        
    std::vector<double> PriVtxWeights;
        
    for( int bin=1; bin<h_vertex_weights->GetNbinsX(); bin++)
      {
	PriVtxWeights.push_back( h_vertex_weights->GetBinContent(bin) );
	//	std::cout<<bin-1 <<" "<<h_vertex_weights->GetBinContent(bin)<<std::endl;
      }   
    
    
    m_primary_vertex_distribution_weights = PriVtxWeights;
   
    m_vertex_weight = PriVtxWeights.at( this->NumberOfPrimaryVertices() );    
  }
}


void Event::SetPileupWeight( const double weight ){
  
  if( UsePileupWeighting == false){
    m_pileup_weight = 1.0;
  }else{
    m_pileup_weight = weight; 
  }

}

void Event::SetPtReweight(const double propagatorPt, map<TString,TH1D*> h_hist){
    
  if( UsePtReweighting == true){ 
    //    m_mc_pt_reweight = this->CalculatePtReweight( propagatorPt );
    m_mc_pt_reweight = PtReweighting::Weight( propagatorPt, h_hist );
  }else{
    m_mc_pt_reweight = 1.0;
  }
}

double Event::CalculatePtReweight(const double propagatorPt){
 
  double weight = 1.0;

  
  if( (propagatorPt >= 0) && (propagatorPt < 3 ) ){
    weight = 1.36289;
  }
  
  if( (propagatorPt >= 3) && (propagatorPt < 6 ) ){
    weight = 1.12028;
  }
  
  if( (propagatorPt >= 6) && (propagatorPt < 9 ) ){
    weight = 0.992827;
  }
  
  if( (propagatorPt >= 9) && (propagatorPt < 12 ) ){
    weight = 0.945103;
  }
  
  if( (propagatorPt >= 12) && (propagatorPt < 15 ) ){
    weight = 0.933788;
  }
  
  if( (propagatorPt >= 15) && (propagatorPt < 18 ) ){
    weight = 0.907847;
  }
  
  if( (propagatorPt >= 18) && (propagatorPt < 21 ) ){
    weight = 0.935696;
  }
  
  if( (propagatorPt >= 21) && (propagatorPt < 24 ) ){
    weight = 0.903103;
  }
  
  if( (propagatorPt >= 24) && (propagatorPt < 27 ) ){
    weight = 0.903833;
  }
  
  if( (propagatorPt >= 27) && (propagatorPt < 30 ) ){
    weight = 0.966665;
  }
  
  if( (propagatorPt >= 30) && (propagatorPt < 35 ) ){
    weight = 0.891085;
  }
  
  if( (propagatorPt >= 35) && (propagatorPt < 40 ) ){
    weight = 0.985489;
  }
  
  if( (propagatorPt >= 40) && (propagatorPt < 45 ) ){
    weight = 0.973642;
  }
  
  if( (propagatorPt >= 45) && (propagatorPt < 50 ) ){
    weight = 1.03534;
  }
  
  if( (propagatorPt >= 50) && (propagatorPt < 55 ) ){
    weight = 1.01307;
  }
  
  if( (propagatorPt >= 55) && (propagatorPt < 60 ) ){
    weight = 1.12935;
  }
  
  if( propagatorPt >= 60 ){
    weight = 1.0;
  }
  
  return weight;
}



void Event::SetMCEventWeight( const double weight ){

  if( UseMCEventWeighting == false ){
    m_mcevent_weight = 1;
  }else{
    m_mcevent_weight = weight;
  }

}


void Event::SetMuonMCEfficiencyWeight(Analysis::AnalysisMuonEfficiencyScaleFactors* sf){


  if( UseMuonEfficiencyWeighting == false){
    m_mc_muon_eff_weight = 1;
  }else{

 
    double muon1_sf = sf->scaleFactor(this->Muon1().LorentzVector());
    double muon2_sf = sf->scaleFactor(this->Muon2().LorentzVector());
    
    double muon1_sf_stat_err = sf->scaleFactorUncertainty(this->Muon1().LorentzVector());
    double muon1_sf_syst_err = sf->scaleFactorSystematicUncertainty(this->Muon1().LorentzVector());
     
    double muon2_sf_stat_err = sf->scaleFactorUncertainty(this->Muon2().LorentzVector());
    double muon2_sf_syst_err = sf->scaleFactorSystematicUncertainty(this->Muon2().LorentzVector());
    
    m_muon1_reco_sf_error = sqrt( pow(muon1_sf_stat_err,2) + pow(muon1_sf_syst_err,2) );
    m_muon2_reco_sf_error = sqrt( pow(muon2_sf_stat_err,2) + pow(muon2_sf_syst_err,2) );

    if( m_syst_variations.find("recoUp")->second == true ){
      m_muon1_reco_sf = muon1_sf + m_muon1_reco_sf_error;
      m_muon2_reco_sf = muon2_sf + m_muon2_reco_sf_error;
    }else if( m_syst_variations.find("recoDown")->second == true ){
      m_muon1_reco_sf = muon1_sf - m_muon1_reco_sf_error;
      m_muon2_reco_sf = muon2_sf - m_muon2_reco_sf_error;
    }else if( m_syst_variations.find("recoStatUp")->second == true ){
      m_muon1_reco_sf = muon1_sf + muon1_sf_stat_err;
      m_muon2_reco_sf = muon2_sf + muon2_sf_stat_err;
    }else if( m_syst_variations.find("recoStatDown")->second == true ){
      m_muon1_reco_sf = muon1_sf - muon1_sf_stat_err;
      m_muon2_reco_sf = muon2_sf - muon2_sf_stat_err;
    }else if( m_syst_variations.find("recoSystUp")->second == true ){
      m_muon1_reco_sf = muon1_sf + muon1_sf_syst_err;
      m_muon2_reco_sf = muon2_sf + muon2_sf_syst_err;
   }else if( m_syst_variations.find("recoSystDown")->second == true ){
      m_muon1_reco_sf = muon1_sf - muon1_sf_syst_err;
      m_muon2_reco_sf = muon2_sf - muon2_sf_syst_err;
    }else{
      m_muon1_reco_sf = muon1_sf;
      m_muon2_reco_sf = muon2_sf;
    }
    
    m_mc_muon_eff_weight = m_muon1_reco_sf * m_muon2_reco_sf; 

  }
  
  return;
}

void Event::SetTriggerEfficiencyWeight( const double sf ){
  if( !UseTriggerEfficiencyWeighting ){
    m_trigger_efficiency_weight = 1.0;
  }else{
    m_trigger_efficiency_weight = sf;
  }
  return;
}


void Event::SetTriggerEfficiencyWeight( const std::map< TString, TH2F* > h_hist ){

  if( UseTriggerEfficiencyWeighting == false ){
    m_trigger_efficiency_weight = 1 ;
  }else{
    
    double muon1_sf = MuonTriggerEfficiency::ScaleFactor( this->Muon1(), h_hist, this->SystematicVariations() );
    double muon2_sf = MuonTriggerEfficiency::ScaleFactor( this->Muon2(), h_hist, this->SystematicVariations() );

    m_trigger_efficiency_weight = muon1_sf * muon2_sf;
  }
}

void Event::SetIsolationEfficiencyWeight( int sampleID, const std::map< TString, TH1D* > h_hist ){

  bool signalSample = false;
  if( (sampleID == 113700) || (sampleID ==113701) || (sampleID == 108319) || (sampleID == 106047) || (sampleID == 108321) ||  (sampleID == 113711) || (sampleID == 113712)  ){ 
    signalSample = true; 
  }

  if( (UseIsolationEfficiencyWeighting == false) || (signalSample == false) ){
    m_isolation_efficiency_weight = 1.0 ;
  }else{
    
    double muon1_sf = MuonIsolationEfficiency::ScaleFactor( this->Muon1(), h_hist, this->SystematicVariations() );
    double muon2_sf = MuonIsolationEfficiency::ScaleFactor( this->Muon2(), h_hist, this->SystematicVariations() );

    m_isolation_efficiency_weight = muon1_sf * muon2_sf;
  }

  return;
}


void Event::SetVertexPositionWeight( const float &weight ){ m_vertex_position_weight = weight; }




void Event::SetMCPdfScale(const float scale){ m_mc_pdf_scale = scale; }
void Event::SetMCPdf1(    const float pdf1 ){ m_mc_pdf1      = pdf1;  }
void Event::SetMCPdf2(    const float pdf2 ){ m_mc_pdf2      = pdf2;  }
void Event::SetMCPdfX1(   const float x1   ){ m_mc_pdf_x1    = x1;    }
void Event::SetMCPdfX2(   const float x2   ){ m_mc_pdf_x2    = x2;    }
void Event::SetMCPdfId1(  const int   id1  ){ m_mc_pdf_id1   = id1;   }
void Event::SetMCPdfId2(  const int   id2  ){ m_mc_pdf_id2   = id2;   }



bool Event::PassesMuonPreSelection(const int nMuons, const double ptMin, const double d0Max){
  
  bool passes_cut = false;
  
  int pass_count =0;
  for( int i=0; i<static_cast<int>(this->Muons().size()); i++)
    {	
      if( (this->Muons().at(i).Pt() > ptMin)  && (fabs(this->Muons().at(i).D0()) < d0Max) ){ pass_count++;}
    }
  
  if( pass_count >= nMuons ){ passes_cut = true;}
  
  return passes_cut;
     
}


//-- I think this function is no longer of use....
void Event::PreSelectionD0Cut( std::vector<Particle>& muon_candidates, double d0_cut ){
   
  std::vector<Particle> muons_that_pass;
 
  for( int npart=0; npart < static_cast<int>(muon_candidates.size()); npart++)
    {      
      if( muon_candidates.at(npart).D0() >= d0_cut ){  muons_that_pass.push_back( muon_candidates.at(npart) ); }
    }

  muon_candidates = muons_that_pass;
  
}




bool Event::PassesPVWithAssociatedTracksCut( const std::vector<Vertex> vertices, int vtx_num, int trk_num ){

  bool event_passes = false;

  int passed_vtx = 0;
  for(int vtx=0; vtx<static_cast<int>(vertices.size()); vtx++)
    {
      if( (vertices.at(vtx).NumberOfTrks() > trk_num) && (vertices.at(vtx).Type() == 1) ){ passed_vtx++; }
    }
    
  if( passed_vtx >= vtx_num ){ event_passes = true; }
  
  return event_passes;
  
}


bool Event::BothMuonsFromSameVertex( const double chi2_cut ){


  //-- Only interested in the 0th vertex...

  bool muon1matched = false;
  bool muon2matched = false;

  // Error check. All Vertices().at(0) should be type 1
  if( this->Vertices().at(0).Type() != 1 ){ cout<<"ERROR: Vertex Type = "<< this->Vertices().at(0).Type()<<endl;} 

  if( (this->Vertices().at(0).Chi2() < chi2_cut) || (chi2_cut == -1) ){
    for(int trk=0; trk<this->Vertices().at(0).NumberOfTrks(); trk++)
      {
	int index = this->Vertices().at(0).VertexTrkIndex().at(trk);
	if( this->Muon1().IDD0() == this->Trks().at(index).D0() ){ muon1matched = true; }
	if( this->Muon2().IDD0() == this->Trks().at(index).D0() ){ muon2matched = true; }
      }
  }

  bool  bothmuonsfromvertex = false;
  if( (muon1matched == true) && (muon2matched == true) ){ bothmuonsfromvertex = true;}

  return bothmuonsfromvertex;

}

bool Event::BothMuonsPassMSPtCut( const double cut )
{
  bool event_passes = false;
  if( (this->Muon1().MSVector3().Pt() > cut) && (this->Muon2().MSVector3().Pt() > cut) ){ event_passes = true; }
  return event_passes;
}

bool Event::BothMuonsPassMSIDPtDiffCut( const double cut )
{
  bool event_passes = false;
  bool muon1_passes = false;
  bool muon2_passes = false;

  if( fabs(this->Muon1().MSVector3().Pt() - this->Muon1().IDVector3().Pt() ) < cut ){ muon1_passes = true;}
  if( fabs(this->Muon2().MSVector3().Pt() - this->Muon2().IDVector3().Pt() ) < cut ){ muon2_passes = true; }

  if( (muon1_passes == true) && (muon2_passes == true) ){ event_passes = true; }

  return event_passes;
}


bool Event::BothMuonsPassZ0ZVtxDiffCut( const double cut )
{
  bool event_passes = false;
  bool muon1_pass   = false;
  bool muon2_pass   = false;

  if( fabs( this->Muon1().Z0() - this->ZerothVertex().Z() ) < cut ){ muon1_pass = true; }
  if( fabs( this->Muon2().Z0() - this->ZerothVertex().Z() ) < cut ){ muon2_pass = true; }
  
  if( muon1_pass == true && muon2_pass == true ){ event_passes = true; }
  
  return event_passes;
}


bool Event::MuonPassesRel16MCPTrackQualityCuts(const Particle muon){

  bool muon_passes = false;

  //-- Track Quality --------------------

  bool passes_blayer_cut = false;
  if( (!muon.ExpectBLayerHit())||( muon.NumberOfBLayerHits()>0) ){ passes_blayer_cut = true;}

  bool passes_pixel_cut  = false;
  if( ( muon.NumberOfPixelHits() + muon.NumberOfPixelDeadSensors() )>1 ){ passes_pixel_cut = true;}

  bool passes_sct_cut = false;
  if( (muon.NumberOfSCTHits() + muon.NumberOfSCTDeadSensors() ) >=6){ passes_sct_cut = true; }


  bool passes_hole_cut = false;
  if( (muon.NumberOfPixelHoles() + muon.NumberOfSCTHoles() ) < 2){ passes_hole_cut = true;}

  bool passes_trt_cut = false;
  int n_trt = muon.NumberOfTRTHits() + muon.NumberOfTRTOutliers() ;

  //Case 1
  if( fabs( muon.IDVector3().Eta() ) < 1.9 ){
    if( n_trt > 5 && muon.NumberOfTRTOutliers() < 0.9*n_trt ){ passes_trt_cut = true; }
  }
  //Case 2
  if( fabs( muon.IDVector3().Eta() ) >= 1.9 ){
    if( n_trt > 5 ){
      if( muon.NumberOfTRTOutliers() < 0.9*n_trt ){ passes_trt_cut = true;}
    }else{
      passes_trt_cut = true;
    }
  }

  if( (passes_blayer_cut == true) && (passes_pixel_cut == true) && (passes_sct_cut == true) && (passes_hole_cut == true) && (passes_trt_cut == true) ){
    muon_passes = true;
  }

  return muon_passes;
}


bool Event::BothMuonsPassRel16MCPTrackQualityCuts(){

  bool event_passes = false;
  if( (this->MuonPassesRel16MCPTrackQualityCuts(this->Muon1()) == true) && (this->MuonPassesRel16MCPTrackQualityCuts(this->Muon2()) == true) ){
    event_passes = true;
  }

  return event_passes;
}

bool Event::MuonPassesRel17MCPTrackQualityCuts(const Particle muon){

  bool muon_passes = false;

  //-- Track Quality --------------------

  bool passes_blayer_cut = false;
  if( (!muon.ExpectBLayerHit())||( muon.NumberOfBLayerHits()>0) ){ passes_blayer_cut = true;}

  bool passes_pixel_cut  = false;
  if( ( muon.NumberOfPixelHits() + muon.NumberOfPixelDeadSensors() )>1 ){ passes_pixel_cut = true;}

  bool passes_sct_cut = false;
  if( (muon.NumberOfSCTHits() + muon.NumberOfSCTDeadSensors() ) >5 ){ passes_sct_cut = true; }


  bool passes_hole_cut = false;
  if( (muon.NumberOfPixelHoles() + muon.NumberOfSCTHoles() ) < 3){ passes_hole_cut = true;}

  bool passes_trt_cut = false;
  int n_trt = muon.NumberOfTRTHits() + muon.NumberOfTRTOutliers() ;

  //Case 1
  if( fabs( muon.IDVector3().Eta() ) < 1.9 ){
    if( n_trt > 5 && muon.NumberOfTRTOutliers() < 0.9*n_trt ){ passes_trt_cut = true; }
  }
  //Case 2
  if( fabs( muon.IDVector3().Eta() ) >= 1.9 ){
    if( n_trt > 5 ){
      if( muon.NumberOfTRTOutliers() < 0.9*n_trt ){ passes_trt_cut = true;}
    }else{
      passes_trt_cut = true;
    }
  }

  if( (passes_blayer_cut == true) && (passes_pixel_cut == true) && (passes_sct_cut == true) && (passes_hole_cut == true) && (passes_trt_cut == true) ){
    muon_passes = true;
  }

  return muon_passes;
}


bool Event::BothMuonsPassRel17MCPTrackQualityCuts(){

  bool event_passes = false;

  if( (this->Muon1().PassRel17MCPTrackQuality() == true) && (this->Muon2().PassRel17MCPTrackQuality() == true) ){
    event_passes = true;
  }

  return event_passes;
}

bool Event::MuonInLArHole(){
  
  // LAr Hole = -0.1<eta<1.5 -0.1 and -0.9<phi<-0.5
  
  bool pass = false;

  if( (this->Muon1().Eta()> -0.1) && (this->Muon1().Eta() < 1.5) && (this->Muon1().Phi() > -0.9) && (this->Muon1().Phi() < -0.5 ) ){ pass = true; }
  if( (this->Muon2().Eta()> -0.1) && (this->Muon2().Eta() < 1.5) && (this->Muon2().Phi() > -0.9) && (this->Muon2().Phi() < -0.5 ) ){ pass = true; }
 
  return pass;

}


bool Event::WantPtSmearedMuons() const { return  UsePtSmearedMuonsInAnalysis; }


void Event::UseMCPtSmearedMuons(){

// This function takes the smeared pt and passes it back to the "main block" of muon
// variables. This will change the TLorentzVector private member variable so will also change
// variables linked to the pt.

  if( UsePtSmearedMuonsInAnalysis == true ){
    
    Particle tmpMuon1 = this->Muon1();
    Particle tmpMuon2 = this->Muon2();

    tmpMuon1.SetPtEtaPhiMass( this->Muon1().MCRecSmearedPt(), this->Muon1().Eta(), this->Muon1().Phi(), this->Muon1().M() );
    tmpMuon2.SetPtEtaPhiMass( this->Muon2().MCRecSmearedPt(), this->Muon2().Eta(), this->Muon2().Phi(), this->Muon2().M() );
    
    this->SetMuon1( tmpMuon1 );
    this->SetMuon2( tmpMuon2 );

  }

}


bool Event::PassesEventSelection(){

  int NumEventSelectionsTurnedOn    = 0;
  int NumberOfEventSelectionsPassed = 0;

  if( (PreSelection_n_vertices !=-1)  && (PreSelection_n_tracks_to_vertex !=-1) ){
    NumEventSelectionsTurnedOn++;
    if(  this->PassesPVWithAssociatedTracksCut( this->Vertices(), PreSelection_n_vertices, PreSelection_n_tracks_to_vertex ) == true ){
     NumberOfEventSelectionsPassed++;
    }
  }

  if( PreSelection_n_leptons !=-1 ){
    NumEventSelectionsTurnedOn++;
    if(  this->NumberOfMuons() >= PreSelection_n_leptons ){
      NumberOfEventSelectionsPassed++;
    }
  }


  if( PreSelection_require_trigger == true ){
    NumEventSelectionsTurnedOn++;
    if( this->TriggerFired() == true ){ NumberOfEventSelectionsPassed++; }
  }


   if( PreSelection_max_vtx_z !=-1 ){
     NumEventSelectionsTurnedOn++;
     if(  fabs(this->ZerothVertex().Z()) < PreSelection_max_vtx_z ){ NumberOfEventSelectionsPassed++; }
   }
  

   if( PreSelection_reject_high_lumi_mc == true ){
     NumEventSelectionsTurnedOn++;
     if( this->MCSimPeriod() != 185761 ){ NumberOfEventSelectionsPassed++; }
   }

   if( PreSelection_reject_lar_errors == true ){
     NumEventSelectionsTurnedOn++;
     if( this->LArError() == 0 ){ NumberOfEventSelectionsPassed++; }
   }



//-- Check if all of the individual pre-selection cuts are passed:
//-- The trigger and good run list requirements are checked before the code for the event.

  bool event_passes = false;
  if( NumEventSelectionsTurnedOn == NumberOfEventSelectionsPassed ){
    event_passes = true;
    if( this->IsVerbose() ){std::cout<<"EVENT PASSES EVENT SELECTION"<<std::endl;}
  }
  
  return event_passes;

}


bool Event::PassesPreSelection(){
  
  int NumPreSelectionsTurnedOn    = 0;
  int NumberOfPreSelectionsPassed = 0;


  if( PreSelection_lepton_max_eta != -1 ){
    NumPreSelectionsTurnedOn++;
    if( BothMuonsPassEtaCut( PreSelection_lepton_max_eta ) == true ){ NumberOfPreSelectionsPassed++; }
  }

  if( PreSelection_lepton_max_deltaz0 != -1 ){
    NumPreSelectionsTurnedOn++;
    if( MuonsPassDeltaZ0Cut( PreSelection_lepton_max_deltaz0 ) == true ){ NumberOfPreSelectionsPassed++;}
  }


  if( PreSelection_muons_from_same_vertex == true ){
    NumPreSelectionsTurnedOn++;
    if( this->BothMuonsFromSameVertex( PreSelection_vxp_chi2_cut ) == true ){ NumberOfPreSelectionsPassed++; }
  }


  if( PreSelection_lepton_min_pt  != -1 ){
    NumPreSelectionsTurnedOn++;
    if( this->BothMuonsPassPtCut( PreSelection_lepton_min_pt ) == true ){ NumberOfPreSelectionsPassed++; }
  }

  if( PreSelection_lepton_max_d0  != -1 ){
    NumPreSelectionsTurnedOn++;
    if( this->BothMuonsPassD0Cut( PreSelection_lepton_max_d0 ) == true ){ NumberOfPreSelectionsPassed++; }
  }
			     
  if( PreSelection_min_MS_pt != -1 ){
    NumPreSelectionsTurnedOn++;
    if( this->BothMuonsPassMSPtCut( PreSelection_min_MS_pt ) == true ){ 
      NumberOfPreSelectionsPassed++; 
    }
  }

  if( PreSelection_min_MS_ID_pt_diff != -1 ){
    NumPreSelectionsTurnedOn++;
    if( this->BothMuonsPassMSIDPtDiffCut( PreSelection_min_MS_ID_pt_diff ) == true ){
      NumberOfPreSelectionsPassed++; 
    }
  }

  if( PreSelection_min_z0_zvtx_diff != -1 ){
    NumPreSelectionsTurnedOn++;
    if( this->BothMuonsPassZ0ZVtxDiffCut( PreSelection_min_z0_zvtx_diff ) == true ){
      NumberOfPreSelectionsPassed++; 
    }
  }

  
  if( PreSelection_reject_lar_hole == true ){
    NumPreSelectionsTurnedOn++;
    if( this->MuonInLArHole() == false ){ 
      NumberOfPreSelectionsPassed++;
    }
  }



//-- Check if all of the individual pre-selection cuts are passed:
//-- The trigger and good run list requirements are checked before the code for the event.
  bool event_passes = false;
 if( NumPreSelectionsTurnedOn == NumberOfPreSelectionsPassed ){
    event_passes = true;
    if ( this->IsVerbose() ){std::cout<<"EVENT PASSES PRESELECTION"<<std::endl;}
  }

 

  
  return event_passes;

}



bool Event::PassesSelectionWithoutIso(){


  int NumberOfSelectionsTurnedOn = 0;
  int NumberOfSelectionsPassed   = 0;

  if( Selection_n_leptons != -1 ){
    NumberOfSelectionsTurnedOn++;
    if( this->NumberOfMuons() >= Selection_n_leptons ){ NumberOfSelectionsPassed++;}    
  }

  if( Selection_require_opposite_charge_leptons == true ){
    NumberOfSelectionsTurnedOn++;
    if( (this->Muon1().Charge()*this->Muon2().Charge()) < 0 ){ NumberOfSelectionsPassed++;}
  }

  if( Selection_muon1_min_pt != -1 && Selection_muon2_min_pt != -1 ){
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsPassPtCut(Selection_muon1_min_pt, Selection_muon2_min_pt ) == true){ NumberOfSelectionsPassed++;}
  }
  
  if( Selection_lepton_max_eta != -1 ){   
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsPassEtaCut(Selection_lepton_max_eta) == true){ NumberOfSelectionsPassed++;}
  }

  if( Selection_lepton_max_d0 != -1 ){
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsPassD0Cut(Selection_lepton_max_d0) == true ){ NumberOfSelectionsPassed++;}
  }

  if( Selection_lepton_max_z0 != -1 ){
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsPassZ0Cut(Selection_lepton_max_z0) == true ){ NumberOfSelectionsPassed++;}
  }
  
  if(  Selection_lepton_max_d0_sig != -1 ){
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsPassD0SignificanceCut( Selection_lepton_max_d0_sig ) ){ NumberOfSelectionsPassed++; }
  }
  
  if( Selection_min_mass != -1 ){
    NumberOfSelectionsTurnedOn++;
    if( this->Propagator().M() > Selection_min_mass){NumberOfSelectionsPassed++;}
  }

  if( Selection_max_mass != -1 ){
    NumberOfSelectionsTurnedOn++;
    if( this->Propagator().M() < Selection_max_mass ){ NumberOfSelectionsPassed++; }
  }


  if( Selection_require_MCPRel16_track_quality_cuts == true ){
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsPassRel16MCPTrackQualityCuts() == true ){ NumberOfSelectionsPassed++; }
  }

  if( Selection_require_MCPRel17_track_quality_cuts == true ){
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsPassRel17MCPTrackQualityCuts() == true ){ NumberOfSelectionsPassed++; }
  }


  if( Selection_require_trigger == true ){
    NumberOfSelectionsTurnedOn++;
    if( this->TriggerFired() == true ){ NumberOfSelectionsPassed++; }
  }


  if( Selection_max_prop_rap != -1 ){
    NumberOfSelectionsTurnedOn++;
    if( this->PropagatorPassesRapidityCut(Selection_max_prop_rap) == true ){ NumberOfSelectionsPassed++; }
  }


  bool event_passes = false;  

//--- Test if event passes all the selections
  if ( NumberOfSelectionsTurnedOn == NumberOfSelectionsPassed ){ 

    if( this->IsVerbose() ){std::cout<<"EVENT PASSES SELECTION WITHOUT ISO"<<std::endl;}
    event_passes = true;
  }


  return event_passes;

}
 
bool Event::PassesIsolationSelection(){
  
  int NumberOfSelectionsTurnedOn = 0;
  int NumberOfSelectionsPassed   = 0;
  
//- Look at the selection without isolation  <------ This is now already passed before calling this function.
//  NumberOfSelectionsTurnedOn++;
//  if( this->PassesSelectionWithoutIso() ){ NumberOfSelectionsPassed++;}
  

  if( (Selection_iso_et_ratio == true) && (Selection_et_cone != -1) && (Selection_et_cone_ratio_max != -1) ){    
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsPassEtRatioCut(Selection_et_cone, Selection_et_cone_ratio_max) == true ){NumberOfSelectionsPassed++;}
  }
  
  
  if( (Selection_iso_pt_ratio == true) && (Selection_pt_cone != -1) && (Selection_pt_cone_ratio_max != -1) ){
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsPassPtRatioCut(Selection_pt_cone, Selection_pt_cone_ratio_max) == true ){NumberOfSelectionsPassed++;}
  }


  if( (Selection_iso_et_sum == true) && (Selection_et_cone != -1) && (Selection_et_sum_max != -1)  ){
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsPassEtSumCut( Selection_et_cone, Selection_et_sum_max ) == true ){ NumberOfSelectionsPassed++;}
  }

  if( (Selection_iso_pt_sum == true) && (Selection_pt_cone != -1) && (Selection_pt_sum_max != -1) ) {
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsPassPtSumCut( Selection_pt_cone, Selection_pt_sum_max ) == true ){ NumberOfSelectionsPassed++;}
  }
  
  bool event_passes = false;  

//--- Test if event passes all the selections
  if ( NumberOfSelectionsTurnedOn == NumberOfSelectionsPassed ){ 
    
    if( this->IsVerbose() ){std::cout<<"EVENT PASSES SELECTION"<<std::endl;}
    event_passes = true;
  }


  return event_passes;


}





bool Event::PassesNonIsolatedSelection(){


  int NumberOfSelectionsTurnedOn = 0;
  int NumberOfSelectionsPassed   = 0;

//- Look at the selection without isolation  <--- this is now already passed before calling this function
//  NumberOfSelectionsTurnedOn++;
//  if( this->PassesSelectionWithoutIso() == true ){ NumberOfSelectionsPassed++;}

  if( (Selection_iso_et_ratio == true) && (Selection_et_cone != -1) && (Selection_et_cone_ratio_max != -1) ){    
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsFailEtRatioCut(Selection_et_cone, Selection_et_cone_ratio_max) == true ){NumberOfSelectionsPassed++;}
  }


  if( (Selection_iso_pt_ratio == true) && (Selection_pt_cone != -1) && (Selection_pt_cone_ratio_max != -1) ){
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsFailPtRatioCut(Selection_pt_cone, Selection_pt_cone_ratio_max) == true ){NumberOfSelectionsPassed++;}
  }


  if( (Selection_iso_et_sum == true) && (Selection_et_cone != -1) && (Selection_et_sum_max != -1)  ){
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsFailEtSumCut( Selection_et_cone, Selection_et_sum_max ) == true ){ NumberOfSelectionsPassed++;}
  }

  if( (Selection_iso_pt_sum == true) && (Selection_pt_cone != -1) && (Selection_pt_sum_max != -1) ) {
    NumberOfSelectionsTurnedOn++;
    if( this->BothMuonsFailPtSumCut( Selection_pt_cone, Selection_pt_sum_max ) == true  ){ NumberOfSelectionsPassed++;}
  }


  bool event_passes = false;  

//--- Test if event passes all the selections
  if ( NumberOfSelectionsTurnedOn == NumberOfSelectionsPassed ){ 

    if( this->IsVerbose() ){ std::cout<<"EVENT PASSES SELECTION"<<std::endl;}
    event_passes = true;
  }


  return event_passes;

}


void Event::CreateEventSelectionCutFlow( TH1D* h_hist, const double weight ){

  // IF YOU EDIT THIS CHANGE THE BIN LABELS IN HISTOGRAMS.CPPP

  h_hist->Fill(0.0, weight );
  if( this->TriggerFired() == true ){
    h_hist->Fill(1.0, weight);
    if( this->LArError() == 0 ){ 
      h_hist->Fill(2.0, weight);
      if( this->NumberOfMuons() >= PreSelection_n_leptons ){
	h_hist->Fill(3.0, weight);
	if( this->PassesPVWithAssociatedTracksCut( this->Vertices(), PreSelection_n_vertices, PreSelection_n_tracks_to_vertex ) == true ){
	  h_hist->Fill(4.0, weight);
	}
      }
    }
  }
  return; 
}



void Event::CreatePreSelectionCutFlow( TH1D* h_hist, const double weight ){

  // IF YOU EDIT THIS CHANGE THE BIN LABELS IN HISTOGRAMS.CPPP

  h_hist->Fill(0.0, weight );
  if( this->BothMuonsPassPtCut( 6 ) && this->BothMuonsPassEtaCut( 2.5 ) ){
    h_hist->Fill(1.0, this->EventWeight() );
    if( this->BothMuonsFromSameVertex( PreSelection_vxp_chi2_cut ) == true ){ 
      h_hist->Fill(2.0, this->EventWeight() );
    }
  }
  return; 
}




void Event::CreateSelectionCutFlow( TH1D* h_hist, const double weight){

  //
  // IF YOU EDIT THIS CHANGE THE BIN LABELS IN HISTOGRAMS.CPPP
  //

  h_hist->Fill(0.0, weight);
  if( this->TriggerFired() == true ){
    h_hist->Fill(1.0, weight);
    
    if( this->BothMuonsPassRel17MCPTrackQualityCuts() == true ){
      h_hist->Fill(2.0, weight);
      
      if( (this->BothMuonsPassPtCut(Selection_muon1_min_pt, Selection_muon2_min_pt)==true) && (this->BothMuonsPassEtaCut(Selection_lepton_max_eta)==true) ){
	h_hist->Fill(3.0, weight);
	
	if( (this->Muon1().Charge()*this->Muon2().Charge()) < 0 ){ 
	  h_hist->Fill(4.0, weight);
	  
	  if( this->BothMuonsPassPtRatioCut(Selection_pt_cone, Selection_pt_cone_ratio_max)== true ){
	    h_hist->Fill(5.0, weight);
	    
	    if( (this->Propagator().M() > Selection_min_mass) && (this->Propagator().M() < Selection_max_mass) ){
	      h_hist->Fill(6.0, weight);      
	      
	      if( fabs(this->Propagator().Rapidity()) < Selection_max_prop_rap ){
		h_hist->Fill(7.0, weight);      
		
	      }
	    }
	  }
	}else{
	  h_hist->Fill(8.0, weight);		    
	  if( this->BothMuonsPassPtRatioCut(Selection_pt_cone, Selection_pt_cone_ratio_max)== true ){
	    h_hist->Fill(9.0, weight);		    
	    if( (this->Propagator().M() > Selection_min_mass) && (this->Propagator().M() < Selection_max_mass)  ){
	      h_hist->Fill(10.0, weight);		    
	    }///
	  }
	}//end of else Q1*Q2 <0
      }
    }
  }
  return;
}

void Event::CreateFullSelectionCutFlow( TH1D* h_hist, const double weight){

  //
  // IF YOU EDIT THIS CHANGE THE BIN LABELS IN HISTOGRAMS.CPPP
  //
  h_hist->Fill(0.0, weight );
  if( this->LArError() == 0 ){ 
    h_hist->Fill(1.0, weight);
    if( this->NumberOfMuons() >= PreSelection_n_leptons ){
      h_hist->Fill(2.0, weight);
      if( this->PassesPVWithAssociatedTracksCut( this->Vertices(), PreSelection_n_vertices, PreSelection_n_tracks_to_vertex ) == true ){
	h_hist->Fill(3.0, weight);
	
	if( this->BothMuonsFromSameVertex( PreSelection_vxp_chi2_cut ) == true ){ 
	  h_hist->Fill(4.0, weight);
	  
	  h_hist->Fill(5.0, weight);
	  if( this->TriggerFired() == true ){
	    h_hist->Fill(6.0, weight);
	    
	    if( this->BothMuonsPassRel17MCPTrackQualityCuts() == true ){
	      h_hist->Fill(7.0, weight);
	      
	      if( (this->BothMuonsPassPtCut(Selection_muon1_min_pt, Selection_muon2_min_pt)==true) && (this->BothMuonsPassEtaCut(Selection_lepton_max_eta)==true) ){
		h_hist->Fill(8.0, weight);
		
		if( (this->Muon1().Charge()*this->Muon2().Charge()) < 0 ){ 
		  h_hist->Fill(9.0, weight);
		  
		  if( this->BothMuonsPassPtRatioCut(Selection_pt_cone, Selection_pt_cone_ratio_max)== true ){
		    h_hist->Fill(10.0, weight);
		    
		    if( (this->Propagator().M() > Selection_min_mass) && (this->Propagator().M() < Selection_max_mass) ){
		      h_hist->Fill(11.0, weight);      
		      
		      if( fabs(this->Propagator().Rapidity()) < Selection_max_prop_rap ){
			h_hist->Fill(12.0, weight);      
			
		      }
		    }
		  }		  
		}else{
		  h_hist->Fill(13.0, weight);		    
		  if( this->BothMuonsPassPtRatioCut(Selection_pt_cone, Selection_pt_cone_ratio_max)== true ){
		    h_hist->Fill(14.0, weight);		    
		    if( (this->Propagator().M() > Selection_min_mass) && (this->Propagator().M() < Selection_max_mass)  ){
		      h_hist->Fill(15.0, weight);		    
		    }///
		  }
		}//end of else Q1*Q2 <0
	      }
	    }
	  }
	}
      }
    }
  }
  return;
}





void Event::CreateZEventSelectionCutFlow( TH1D* h_hist, const double weight ){

  // IF YOU EDIT THIS CHANGE THE BIN LABELS IN HISTOGRAMS.CPP
  
  h_hist->Fill( 0.0, weight );
  if( this->NumberOfMuons() >= PreSelection_n_leptons ){
    h_hist->Fill( 1.0, weight );
    if( this->LArError() == 0 ){
      h_hist->Fill( 2.0, weight );
      if( this->TriggerFired() == true ){
	h_hist->Fill( 3.0, weight );
	if( this->PassesPVWithAssociatedTracksCut( this->Vertices(), PreSelection_n_vertices, PreSelection_n_tracks_to_vertex ) == true ){
	  h_hist->Fill( 4.0, weight );
	} //vertex
      } //trigger
    }// LAr Error
  }// nMuons
  
}






void Event::CreateZPreSelectionCutFlow( TH1D* h_hist, const double weight ){

  // IF YOU EDIT THIS CHANGE THE BIN LABELS IN HISTOGRAMS.CPP

  h_hist->Fill(0.0, weight );
  if( this->MuonInLArHole() == false ){
    h_hist->Fill(1.0, weight );     
  }//LAr Hole Cut
  
}




void Event::CreateZSelectionCutFlow( TH1D* h_hist, const double weight){

  //
  // IF YOU EDIT THIS CHANGE THE BIN LABELS IN HISTOGRAMS.CPP
  //
  
  h_hist->Fill(0.0, weight);
  if( this->BothMuonsPassZ0Cut(Selection_lepton_max_z0 ) ){
    h_hist->Fill(1.0, weight);
    if( this->BothMuonsPassD0SignificanceCut( Selection_lepton_max_d0_sig ) ){
      h_hist->Fill(2.0, weight);
      if( this->BothMuonsPassRel17MCPTrackQualityCuts() == true ){
	h_hist->Fill(3.0, weight);
	if( (this->BothMuonsPassPtCut(Selection_muon1_min_pt, Selection_muon2_min_pt) == true) && (this->BothMuonsPassEtaCut(Selection_lepton_max_eta)==true) ){
	  h_hist->Fill(4.0, weight);
	  if( this->BothMuonsPassPtRatioCut(Selection_pt_cone, Selection_pt_cone_ratio_max) == true ) {
	    h_hist->Fill(5.0, weight);
	    if( (this->Muon1().Charge()*this->Muon2().Charge()) < 0 ){ 
	      h_hist->Fill(6.0, weight);
	      if( (this->Propagator().M() > Selection_min_mass) && (this->Propagator().M() < Selection_max_mass) ){
		h_hist->Fill(7.0, weight);      
	      } 
	    }
	  }
	}
      }
    }
  }
}


void Event::CreateCutList( TProfile* p_cutlist )
{

//-- Cut List:
//-- Add to this list as cuts are added!
//--
//
//
// 0: Vertex Selection
// 1: PreSel N leptons
// 2: PreSel Trigger
// 3: PreSel zvtx
// 4: PreSel muon eta
// 5: PreSel muon pt
// 6: PreSel muon deltaz0
// 7: PreSel muons from same vertex
// 8: PreSel muons max d0
// 9: PreSel muon min ms pt
//10: PreSel muon ms-pt id-pt diff 
//11: PreSel min z0/zvx diff
//12: 
//13:
//14: 
//15: Sel N muons
//16: Sel opp charge
//17: Sel min pt
//18: Sel max eta
//19: Sel max d0
//20: Sel max z0
//21: Sel et ratio
//22: Sel pt ratio
//23: Sel et sum
//24: Sel pt sum
//25: Sel min mass
//26: Sel max mass
//27: Sel track quality
//28: Sel trigger
//29: Sel prop rap

  if( p_cutlist->GetEntries() == 0 ){

    TString cuts[30] = { "Vertex Selection", "PreSel N leptons", "PreSel Trigger", "|z_{vtx}|", "PreSel muon eta", "PreSel muon pt", " PreSel muon deltaz0", "PreSel muons from same vertex","PreSel muons max d0", "PreSel muon min ms pt", "PreSel muon ms-pt id-pt diff","PreSel min z0/zvx diff","", "","","Sel N muons","Sel opp charge", "Sel min pt", "Sel max eta", "Sel max d0", "Sel max z0", "Sel et ratio", "Sel pt ratio", "Sel et sum", "Sel pt sum", "Sel min mass", "Sel max mass", "Sel track quality", "Sel trigger","Sel Prop Rap"};


    for( int bin = 0; bin<p_cutlist->GetNbinsX(); bin++)
      {
	p_cutlist->GetXaxis()->SetBinLabel( bin+1, cuts[bin] );
      }



//- Event Selection
    if( PreSelection_n_vertices !=-1  && PreSelection_n_tracks_to_vertex !=-1 ){                              p_cutlist->Fill( 0.0, 1.0 ); }
    if( PreSelection_n_leptons !=-1 && PreSelection_lepton_min_pt !=-1 && PreSelection_lepton_max_d0 != -1 ){ p_cutlist->Fill( 1, 1.0 ); }
    if( PreSelection_require_trigger == true ){                                                               p_cutlist->Fill( 2, 1.0 ); }
    if( PreSelection_max_vtx_z !=-1 ){                                                                        p_cutlist->Fill( 3, 1.0 ); }
//- PreSelection
    if( PreSelection_lepton_max_eta != -1 ){           p_cutlist ->Fill( 4,  1.0 ); }
    if( PreSelection_lepton_min_pt  != -1 ){           p_cutlist ->Fill( 5,  1.0 ); }
    if( PreSelection_lepton_max_deltaz0 != -1 ){       p_cutlist ->Fill( 6,  1.0 ); }
    if( PreSelection_muons_from_same_vertex == true ){ p_cutlist ->Fill( 7,  1.0 ); }
    if( PreSelection_lepton_max_d0  != -1 ){           p_cutlist ->Fill( 8,  1.0 ); }
    if( PreSelection_min_MS_pt != -1 ){                p_cutlist ->Fill( 9,  1.0 ); }
    if( PreSelection_min_MS_ID_pt_diff != -1 ){        p_cutlist ->Fill( 10, 1.0 ); }
    if( PreSelection_min_z0_zvtx_diff != -1 ){         p_cutlist ->Fill( 11, 1.0 ); }
//- Selection
    if( Selection_n_leptons != -1 ){                          p_cutlist ->Fill( 15,  1.0 ); }
    if( Selection_require_opposite_charge_leptons == true ){  p_cutlist ->Fill( 16,  1.0 ); }
    if( Selection_muon1_min_pt != -1 ){                      p_cutlist ->Fill( 17,  1.0 ); }
    if( Selection_lepton_max_eta != -1 ){                     p_cutlist ->Fill( 18,  1.0 ); }
    if( Selection_lepton_max_d0 != -1 ){                      p_cutlist ->Fill( 19,  1.0 ); }
    if( Selection_lepton_max_z0 != -1 ){                      p_cutlist ->Fill( 10,  1.0 ); }
    if( Selection_iso_et_ratio == true && Selection_et_cone != -1 && Selection_et_cone_ratio_max != -1 ){ p_cutlist ->Fill( 21,  1.0 ); }
    if( Selection_iso_pt_ratio == true && Selection_pt_cone != -1 && Selection_pt_cone_ratio_max != -1 ){ p_cutlist ->Fill( 22,  1.0 ); }
    if( Selection_iso_et_sum == true && Selection_et_cone != -1 && Selection_et_sum_max != -1  ){         p_cutlist ->Fill( 23,  1.0 ); }
    if( Selection_iso_pt_sum == true && Selection_pt_cone != -1 && Selection_pt_sum_max != -1) {          p_cutlist ->Fill( 24,  1.0 ); }

    if( Selection_min_mass != -1 ){                               p_cutlist ->Fill( 25,  1.0 ); }
    if( Selection_max_mass != -1 ){                               p_cutlist ->Fill( 26,  1.0 ); }
    if( Selection_require_MCPRel16_track_quality_cuts == true ){  p_cutlist ->Fill( 27,  1.0 ); }
    if( Selection_require_trigger == true ){                      p_cutlist ->Fill( 28,  1.0 ); }
    if( Selection_max_prop_rap !=-1 ){                            p_cutlist ->Fill( 29,  1.0 ); }


  }

  


}
