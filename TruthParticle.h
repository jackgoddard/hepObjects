
#ifndef TRUTH_PARTICLE_H
#define TRUTH_PARTICLE_H

#include "SimpleParticle.h"

class TruthParticle : public SimpleParticle{

 public:
  TruthParticle();
  TruthParticle( const TLorentzVector& particle, const float& charge, const int& pdgid );
  ~TruthParticle();

  int PDGID() const;
  
 private:
  int m_pdgid;




};

#endif
