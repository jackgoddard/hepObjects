#define MCBLOCK_CXX

#include "MCBlock.h"

#include <iostream>
#include <vector>

using std::cout;
using std::endl;


MCBlock::MCBlock(){
  m_n  = -1;
  m_pt = 0.0;
  m_m  = 0.0;
  m_eta= 0.0;
  m_phi= 0.0;
  m_charge = 0;
  m_status =0;
  m_barcode = 0;
  m_pdgid   = 0;
}

MCBlock::~MCBlock(){}



int   MCBlock::N()       const { return m_n;       }
float MCBlock::Pt()      const { return m_pt;      }
float MCBlock::M()       const { return m_m;       }
float MCBlock::Eta()     const { return m_eta;     }
float MCBlock::Phi()     const { return m_phi;     }
float MCBlock::Charge()  const { return m_charge;  }
int   MCBlock::Status()  const { return m_status;  }
int   MCBlock::BarCode() const { return m_barcode; }
int   MCBlock::PDGID()   const { return m_pdgid;   }

std::vector<int> MCBlock::ParentIndex() const { return m_parent_index; }
std::vector<int> MCBlock::ChildIndex()  const { return m_child_index;  }

std::vector<int> MCBlock::Parents()  const { return m_parents;  }
std::vector<int> MCBlock::Children() const { return m_children; }

TLorentzVector MCBlock::LorentzVector() const {return m_lorentz_vector;}

//-- Setter Functions

void MCBlock::SetN( const int n ){ m_n = n; }
void MCBlock::SetPt( const float pt ) { m_pt = pt; }
void MCBlock::SetM( const float m ){ m_m = m; }
void MCBlock::SetEta( const float eta ){ m_eta = eta;}
void MCBlock::SetPhi( const float phi ){ m_phi = phi; }
void MCBlock::SetCharge( const float charge ){ m_charge = charge; }
void MCBlock::SetStatus( const int status ){ m_status = status; }
void MCBlock::SetBarCode( const int barcode ){ m_barcode = barcode; }
void MCBlock::SetPDGID( const int pdgid ){ m_pdgid = pdgid; }

void MCBlock::SetParents( const std::vector<int> parents ){   m_parents = parents;   }
void MCBlock::SetParentIndex( const std::vector<int> index ){ m_parent_index = index;}
void MCBlock::SetChildren( const std::vector<int> children){  m_children = children; }
void MCBlock::SetChildIndex( const std::vector<int> index){   m_child_index = index; }

void MCBlock::SetLorentzVector( const TLorentzVector vector ){ m_lorentz_vector = vector;}
