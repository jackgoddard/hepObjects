#define MCMUONTRUTH_CXX

#include "MCMuonTruth.h"

#include <iostream>
#include <vector>

using std::cout;
using std::endl;


MCMuonTruth::MCMuonTruth(){
  m_n       = -1;
  m_pt      = 0.0;
  m_m       = 0.0;
  m_eta     = 0.0;
  m_phi     = 0.0;
  m_charge  = 0;
  m_origin  = 0;
  m_barcode = 0;
  m_pdgid   = 0;
  m_type    = 0;
}

MCMuonTruth::~MCMuonTruth(){}



int   MCMuonTruth::N()       const { return m_n;       }
float MCMuonTruth::Pt()      const { return m_pt;      }
float MCMuonTruth::M()       const { return m_m;       }
float MCMuonTruth::Eta()     const { return m_eta;     }
float MCMuonTruth::Phi()     const { return m_phi;     }
float MCMuonTruth::Charge()  const { return m_charge;  }
int   MCMuonTruth::Origin()  const { return m_origin;  }
int   MCMuonTruth::BarCode() const { return m_barcode; }
int   MCMuonTruth::PDGID()   const { return m_pdgid;   }
int   MCMuonTruth::Type()    const { return m_type;    }

TLorentzVector MCMuonTruth::LorentzVector() const { return m_lorentz_vector;}

//-- Setter Functions

void MCMuonTruth::SetN( const int n ){             m_n       = n;       } 
void MCMuonTruth::SetPt( const float pt ) {        m_pt      = pt;      } 
void MCMuonTruth::SetM( const float m ){           m_m       = m;       }
void MCMuonTruth::SetEta( const float eta ){       m_eta     = eta;     }
void MCMuonTruth::SetPhi( const float phi ){       m_phi     = phi;     }
void MCMuonTruth::SetCharge( const float charge ){ m_charge  = charge;  }
void MCMuonTruth::SetOrigin( const int origin ){   m_origin  = origin;  }
void MCMuonTruth::SetBarCode( const int barcode ){ m_barcode = barcode; }
void MCMuonTruth::SetPDGID( const int pdgid ){     m_pdgid   = pdgid;   }
void MCMuonTruth::SetType( const int type ){       m_type    = type;    }

void MCMuonTruth::SetLorentzVector( const TLorentzVector vector ){ m_lorentz_vector = vector;}


