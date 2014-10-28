#include "TruthParticle.h"



TruthParticle::TruthParticle(){
  m_pdgid = 0;
}

TruthParticle::~TruthParticle() {};

TruthParticle::TruthParticle(const TLorentzVector& particle, const float& charge, const int& pdgid) : SimpleParticle( particle, charge){
  m_pdgid = pdgid;
}


int TruthParticle::PDGID() const {return m_pdgid; }
