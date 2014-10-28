//  Particle Class
//
//
//  Jack Goddard, February 2011
//


#include <iostream>
#include "Particle.h"

using namespace std;


//-- Constructors:

Particle::Particle(){
  m_particle     = TLorentzVector( 0, 0, 0, 0);
  m_charge       = 100000;
  m_ptcone20     = 0.0;
  m_ptcone30     = 0.0;
  m_ptcone40     = 0.0;
  m_etcone20     = 0.0;
  m_etcone30     = 0.0;
  m_etcone40     = 0.0;
  m_unsmeared_pt = 0.0;
  m_smeared_charge_flip = 0;
  m_id_unsmeared_pt = 0.0;
  m_id_smeared_charge_flip = 0;
}

Particle::Particle( const TLorentzVector particle, const float charge){
  m_particle     = particle;
  m_charge       = charge;
  m_ptcone20     = 0.0;
  m_ptcone30     = 0.0;
  m_ptcone40     = 0.0;
  m_etcone20     = 0.0;
  m_etcone30     = 0.0;
  m_etcone40     = 0.0;
  m_unsmeared_pt = 0.0;
  m_smeared_charge_flip = 0;
  m_id_unsmeared_pt = 0.0;
  m_id_smeared_charge_flip = 0;

}

Particle::Particle( const float part_px, const float part_py, const float part_pz, const float part_E, const float charge){
  m_particle = TLorentzVector( part_px, part_py, part_pz, part_E );
  m_charge   = charge;
  m_ptcone20 = 0.0;
  m_ptcone30 = 0.0;
  m_ptcone40 = 0.0;
  m_etcone20 = 0.0;
  m_etcone30 = 0.0;
  m_etcone40 = 0.0;
  m_unsmeared_pt = 0.0;
  m_smeared_charge_flip = 0;
  m_id_unsmeared_pt = 0.0;
  m_id_smeared_charge_flip = 0;
}


//-- Destructor

Particle::~Particle() {}




//-- Getter Functions

TLorentzVector Particle::LorentzVector() const {return m_particle; }

float Particle::PtCone20() const { return m_ptcone20; }
float Particle::PtCone30() const { return m_ptcone30; }
float Particle::PtCone40() const { return m_ptcone40; }

float Particle::EtCone20() const { return m_etcone20; }
float Particle::EtCone30() const { return m_etcone30; }
float Particle::EtCone40() const { return m_etcone40; }
float Particle::EtCore()   const { return m_et_core;  }

float Particle::Pt()        const { return m_particle.Pt();       }
float Particle::Eta()       const { return m_particle.Eta();      }
float Particle::Phi()       const { return m_particle.Phi();      }
float Particle::Rapidity()  const { return m_particle.Rapidity(); }
float Particle::M()         const { return m_particle.M();        } 
float Particle::Mass()      const { return m_particle.M();        }

float Particle::Px() const { return m_particle.Px(); }
float Particle::Py() const { return m_particle.Py(); }
float Particle::Pz() const { return m_particle.Pz(); }
float Particle::E()  const { return m_particle.E();  }
float Particle::Et() const { return m_particle.Et(); }

float Particle::Charge() const { return m_charge; }

float Particle::D0()         const { return m_d0; }
float Particle::Z0()         const { return m_z0; }
float Particle::D0Error()    const { return m_d0_err; }
float Particle::Z0Error()    const { return m_z0_err; }
float Particle::CovD0Z0()    const { return m_cov_d0z0; }
float Particle::TrkFitChi2() const { return m_trk_fit_chi2; }

int Particle::Author() const {return m_author; }


bool Particle::PassesLoose()  const { return m_passesloose; }
bool Particle::PassesMedium() const { return m_passesmedium; }
bool Particle::PassesTight()  const { return m_passestight; }


int Particle::NumberOfBLayerHits() const { return m_num_hits_blayer; }
int Particle::NumberOfPixelHits()  const { return m_num_hits_pixel;  }
int Particle::NumberOfSCTHits()    const { return m_num_hits_sct;    }
int Particle::NumberOfTRTHits()    const { return m_num_hits_trt;    }

int Particle::ExpectBLayerHit()          const { return m_expect_b_layer_hit;    }
int Particle::NumberOfSCTDeadSensors()   const { return m_num_sct_dead_sensors;  }
int Particle::NumberOfSCTHoles()         const { return m_num_sct_holes;         }
int Particle::NumberOfPixelDeadSensors() const { return m_num_pixel_dead_sensors;}
int Particle::NumberOfPixelHoles()       const { return m_num_pixel_holes;       }
int Particle::NumberOfTRTOutliers()      const { return m_num_trt_outliers;      }

float Particle::IDCharge() const { return m_id_charge;         }
float Particle::IDQoverP() const { return m_id_q_over_p;       }
float Particle::IDTheta()  const { return m_id_vector.Theta(); }
float Particle::IDPhi()    const { return m_id_vector.Phi();   }
float Particle::IDD0()     const { return m_id_d0;             }
float Particle::IDZ0()     const { return m_id_z0;             }
float Particle::IDD0Error()const { return m_id_d0_err;         }
float Particle::IDZ0Error()const { return m_id_z0_err;         }

float Particle::IDUnSmearedPt()const { return m_id_unsmeared_pt; }
float Particle::IDD0Significance() const { 
  float sig = this->IDD0()/this->IDD0Error(); 
  return sig;
}

float Particle::IDZ0Significance() const { 
  float sig = this->IDZ0()/this->IDZ0Error(); 
  return sig;
}


TVector3 Particle::IDVector3() const {
  //  TVector3 vec;
  //vec.SetMagThetaPhi( fabs( 1.0/this->IDQoverP() ), this->IDTheta(), this->IDPhi() );
  return m_id_vector;
}

float Particle::CBQoverP() const { return m_cb_q_over_p; }
float Particle::CBTheta()  const { return m_cb_theta;    }
float Particle::CBPhi()    const { return m_cb_phi;      }
float Particle::CBD0()     const { return m_cb_d0;       }
float Particle::CBZ0()     const { return m_cb_z0;       }

TVector3 Particle::CBVector3() const { 
  TVector3 vec;
  vec.SetMagThetaPhi( fabs(1.0/this->CBQoverP()), this->CBTheta(), this->CBPhi() );
  return vec;
}

float Particle::MSQoverP() const { return m_ms_q_over_p; }
float Particle::MSTheta()  const { return m_ms_theta;    }
float Particle::MSPhi()    const { return m_ms_phi;      }
float Particle::MSD0()     const { return m_ms_d0;       }
float Particle::MSZ0()     const { return m_ms_z0;       } 

TVector3 Particle::MSVector3() const { 
  TVector3 vec;
  vec.SetMagThetaPhi( fabs(1.0/this->MSQoverP() ), this->MSTheta(), this->MSPhi() );
  return vec;
}


float Particle::MEQoverP() const { return m_me_q_over_p; }
float Particle::METheta()  const { return m_me_theta;    }
float Particle::MEPhi()    const { return m_me_phi;      }
float Particle::MED0()     const { return m_me_d0;       }
float Particle::MEZ0()     const { return m_me_z0;       }

TVector3 Particle::MEVector3() const { 
  TVector3 vec;
  vec.SetMagThetaPhi( fabs(1.0/this->MEQoverP()), this->METheta(), this->MEPhi() );
  return vec;
}


float Particle::IEQoverP() const { return m_ie_q_over_p; }
float Particle::IETheta()  const { return m_ie_theta;    }
float Particle::IEPhi()    const { return m_ie_phi;      }
float Particle::IED0()     const { return m_ie_d0;       }
float Particle::IEZ0()     const { return m_ie_z0;       }

TVector3 Particle::IEVector3() const { 
  TVector3 vec;
  vec.SetMagThetaPhi( fabs(1.0/this->IEQoverP()), this->IETheta(), this->IEPhi() );
  return vec;
}

double Particle::SmearedChargeFlip() const{ return m_smeared_charge_flip;  }
float Particle::UnSmearedPt()       const{ return m_unsmeared_pt;         }
float Particle::MCRecSmearedPt()    const{ return m_mc_rec_smeared_pt;    }
float Particle::MCRecIESmearedPt()  const{ return m_mc_rec_ie_smeared_pt; }
float Particle::MCRecMESmearedPt()  const{ return m_mc_rec_me_smeared_pt; }

 
float Particle::PtRatio20() const {
 
  float ptratioVal = this->PtCone20() / this->Pt();
  
  return ptratioVal;
}

float Particle::PtRatio30() const {
 
  float ptratioVal = this->PtCone30() / this->Pt();
  
  return ptratioVal;
 }

float Particle::PtRatio40() const {
 
  float ptratioVal = this->PtCone40() / this->Pt();
  
  return ptratioVal;
}



float Particle::EtRatio20() const {
 
  float etratioVal = this->EtCone20() / this->Pt();
  return etratioVal;
}

float Particle::EtRatio30() const {
  
  float etratioVal = this->EtCone30() / this->Pt();
  return etratioVal;
}

float Particle::EtRatio40() const {
 
  float etratioVal = this->EtCone40() / this->Pt();
  return etratioVal;
}


int   Particle::MCTruthBarCode()       const { return m_mc_truth_barcode;        }
int   Particle::MCTruthMatched()       const { return m_mc_truth_matched;        }
float Particle::MCTruthE()             const { return m_mc_truth_e;              }
float Particle::MCTruthDeltaR()        const { return m_mc_truth_delta_r;        }
float Particle::MCTruthEta()           const { return m_mc_truth_eta;            }
float Particle::MCTruthPhi()           const { return m_mc_truth_phi;            }
float Particle::MCTruthPt()            const { return m_mc_truth_pt;             }
int   Particle::MCTruthStatus()        const { return m_mc_truth_status;         }  
int   Particle::MCTruthType()          const { return m_mc_truth_type;           }
int   Particle::MCTruthMotherType()    const { return m_mc_truth_mother_type;    }
int   Particle::MCTruthMotherBarCode() const { return m_mc_truth_mother_barcode; }





bool Particle::PassRel17MCPTrackQuality() const {

  bool muon_passes = false;

  //-- Track Quality --------------------

  bool passes_blayer_cut = false;
  if( (!this->ExpectBLayerHit())||( this->NumberOfBLayerHits()>0) ){ passes_blayer_cut = true;}

  bool passes_pixel_cut  = false;
  if( ( this->NumberOfPixelHits() + this->NumberOfPixelDeadSensors() )>1 ){ passes_pixel_cut = true;}

  bool passes_sct_cut = false;
  if( (this->NumberOfSCTHits() + this->NumberOfSCTDeadSensors() ) >5 ){ passes_sct_cut = true; }


  bool passes_hole_cut = false;
  if( (this->NumberOfPixelHoles() + this->NumberOfSCTHoles() ) < 3){ passes_hole_cut = true;}

  bool passes_trt_cut = false;
  int n_trt = this->NumberOfTRTHits() + this->NumberOfTRTOutliers() ;

  //Case 1
  if( fabs( this->IDVector3().Eta() ) < 1.9 ){
    if( n_trt > 5 && this->NumberOfTRTOutliers() < 0.9*n_trt ){ passes_trt_cut = true; }
  }
  //Case 2
  if( fabs( this->IDVector3().Eta() ) >= 1.9 ){
    if( n_trt > 5 ){
      if( this->NumberOfTRTOutliers() < 0.9*n_trt ){ passes_trt_cut = true;}
    }else{
      passes_trt_cut = true;
    }
  }

  if( (passes_blayer_cut == true) && (passes_pixel_cut == true) && (passes_sct_cut == true) && (passes_hole_cut == true) && (passes_trt_cut == true) ){
    muon_passes = true;
  }

  return muon_passes;
}






//--- Setter Functions  

void Particle::SetCharge(const float charge){
  m_charge = charge;
}

void Particle::SetLorentzVector(const TLorentzVector particleVec){
  m_particle = particleVec;
}

void Particle::SetPtEtaPhiMass(const float pt, const float eta, const float phi, const float m){
  m_particle.SetPtEtaPhiM( pt, eta, phi, m );
}

void Particle::SetPtCone20( const float ptcone ){
  m_ptcone20 = ptcone;
}

void Particle::SetPtCone30( const float ptcone ){
  m_ptcone30 = ptcone;
}

void Particle::SetPtCone40( const float ptcone ){
  m_ptcone40 = ptcone;
}

void Particle::SetEtCone20( const float etcone ){
  m_etcone20 = etcone;
}

void Particle::SetEtCone30( const float etcone ){
  m_etcone30 = etcone;
}

void Particle::SetEtCone40( const float etcone ){
  m_etcone40 = etcone;
}

void Particle::SetEtCore( const float etcore ){
  m_et_core = etcore;
}

void Particle::SetD0( const float D0 ){
  m_d0 = D0;
}

void Particle::SetZ0( const float Z0 ){
  m_z0 = Z0;
}

void Particle::SetD0Error( const float d0_err ){
  m_d0_err = d0_err;
}

void Particle::SetZ0Error( const float z0_err ){
  m_z0_err = z0_err;
}


void Particle::SetCovD0Z0( const float cov ){
  m_cov_d0z0 = cov;
}

void Particle::SetTrkFitChi2( const float chi2 ){
  m_trk_fit_chi2 = chi2;
}

void Particle::SetPassesLoose( const int loose ){
 
  if( loose > 0 ){ m_passesloose = true;}
  else{ m_passesloose = false; }
 
}

void Particle::SetPassesMedium( const int medium ){

  if( medium > 0 ){ m_passesmedium = true;}
  else{ m_passesmedium = false; }

}

void Particle::SetPassesTight( const int tight ){
  if( tight > 0 ){ m_passestight = true;}
  else{ m_passestight = false; }
 }


void Particle::SetAuthor( const int author ){
  m_author = author;
}

void Particle::SetNumberOfBLayerHits( const int hits ){
  m_num_hits_blayer = hits;
}

void Particle::SetNumberOfPixelHits( const int hits ){
  m_num_hits_pixel = hits;
}

void Particle::SetNumberOfSCTHits( const int hits ){
  m_num_hits_sct = hits;
}

void Particle::SetNumberOfTRTHits( const int hits ){
  m_num_hits_trt = hits;
}

void Particle::SetExpectBLayerHit( const int expect_hit ) {
  m_expect_b_layer_hit = expect_hit;
}

void Particle::SetNumberOfSCTDeadSensors( const int dead_sensors){ 
  m_num_sct_dead_sensors = dead_sensors;
}

void Particle::SetNumberOfSCTHoles( const int holes ){
  m_num_sct_holes = holes;
}

void Particle::SetNumberOfPixelDeadSensors( const int dead_sensors){ 
  m_num_pixel_dead_sensors = dead_sensors;
}

void Particle::SetNumberOfPixelHoles( const int holes ){
  m_num_pixel_holes = holes;
}

void Particle::SetNumberOfTRTOutliers( const int outliers ){
  m_num_trt_outliers = outliers;
}

void Particle::SetIDCharge( const float charge ){
  m_id_charge = charge;
}

void Particle::SetIDUnSmearedPt( const float pt ){
  m_id_unsmeared_pt = pt;
}

void Particle::SetIDD0( const float d0 ){
  m_id_d0 = d0;
}

void Particle::SetIDZ0( const float z0 ){
  m_id_z0 =z0;
}

void Particle::SetIDD0Error( const float d0_err ){
  m_id_d0_err = d0_err;
}

void Particle::SetIDZ0Error( const float z0_err ){
  m_id_z0_err = z0_err;
}


void Particle::SetIDVector3( const TVector3 vec ){
  m_id_vector = vec;
}

void Particle::SetCBQoverP( const float QoverP ){
 m_cb_q_over_p = QoverP;
}

void Particle::SetCBTheta( const float theta ){
 m_cb_theta = theta;
}

void Particle::SetCBPhi( const float phi ){
 m_cb_phi = phi;
}

void Particle::SetCBD0( const float d0 ){
  m_cb_d0 = d0;
}

void Particle::SetCBZ0( const float z0 ){
  m_cb_z0 =z0;
}

void Particle::SetMSQoverP( const float QoverP ){
 m_ms_q_over_p = QoverP;
}

void Particle::SetMSTheta( const float theta ){
 m_ms_theta = theta;
}

void Particle::SetMSPhi( const float phi ){
 m_ms_phi = phi;
}

void Particle::SetMSD0( const float d0 ){
  m_ms_d0 = d0;
}

void Particle::SetMSZ0( const float z0 ){
  m_ms_z0 =z0;
}

void Particle::SetMEQoverP( const float QoverP ){
 m_me_q_over_p = QoverP;
}

void Particle::SetMETheta( const float theta ){
 m_me_theta = theta;
}

void Particle::SetMEPhi( const float phi ){
 m_me_phi = phi;
}

void Particle::SetMED0( const float d0 ){
  m_me_d0 = d0;
}

void Particle::SetMEZ0( const float z0 ){
  m_me_z0 =z0;
}

void Particle::SetIEQoverP( const float QoverP ){
 m_ie_q_over_p = QoverP;
}

void Particle::SetIETheta( const float theta ){
 m_ie_theta = theta;
}

void Particle::SetIEPhi( const float phi ){
 m_ie_phi = phi;
}

void Particle::SetIED0( const float d0 ){
  m_ie_d0 = d0;
}

void Particle::SetIEZ0( const float z0 ){
  m_ie_z0 =z0;
}



//-- MC Truth Setter Functions:

void Particle::SetMCTruthBarCode( const int barcode ){
  m_mc_truth_barcode = barcode;
}

void Particle::SetMCTruthMotherBarCode( const int motherbarcode ){
  m_mc_truth_mother_barcode = motherbarcode;
}

void Particle::SetMCTruthMatched( const int truth_matched ){
  m_mc_truth_matched = truth_matched;
}

void Particle::SetMCTruthE( const float truth_e ){
  m_mc_truth_e = truth_e;
}

void Particle::SetMCTruthDeltaR( const float truth_delta_r ){
  m_mc_truth_delta_r = truth_delta_r;
}

void Particle::SetMCTruthEta( const float truth_eta ){
  m_mc_truth_eta = truth_eta;
}

void Particle::SetMCTruthPhi( const float truth_phi ){
  m_mc_truth_phi = truth_phi;
}

void Particle::SetMCTruthPt( const float truth_pt ){
  m_mc_truth_pt = truth_pt;
}


void Particle::SetMCTruthStatus( const int truth_status){
  m_mc_truth_status = truth_status;
}

void Particle::SetMCTruthType( const int truth_type ){
  m_mc_truth_type = truth_type;
}

void Particle::SetMCTruthMotherType( const int truth_mother_type ){
  m_mc_truth_mother_type = truth_mother_type;
}

void Particle::SetSmearedChargeFlip( const double flip )
{
  m_smeared_charge_flip = flip;
}


void Particle::SetUnSmearedPt( const float pt ){
  m_unsmeared_pt = pt;
}


void Particle::SetMCRecSmearedPt( const float pt ){
  m_mc_rec_smeared_pt = pt;
}

void Particle::SetMCRecIESmearedPt( const float pt ){
  m_mc_rec_ie_smeared_pt = pt;
}


void Particle::SetMCRecMESmearedPt( const float pt ){
  m_mc_rec_me_smeared_pt = pt;
}


void Particle::Print(){
  std::cout<< this->Px() <<"\t"<< this->Py() << "\t"  << this->Pz() << "\t" << this->E() << "\t" << this->Charge() <<std::endl;
}


SimpleParticle Particle::operator+(const Particle& b){

  TLorentzVector combinedVec = this->LorentzVector() + b.LorentzVector();
  
  float combinedCharge = this->Charge() + b.Charge();

  SimpleParticle combined(combinedVec, combinedCharge);

  return combined;

}
