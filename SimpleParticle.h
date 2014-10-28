//  Simple Particle Class
//  
//  Jack Goddard, February 2011
//

#include <TLorentzVector.h>

#ifndef SIMPLE_PARTICLE_H
#define SIMPLE_PARTICLE_H


class SimpleParticle{

 public:
  
  
  SimpleParticle();
  SimpleParticle( const TLorentzVector&, const float&);
  SimpleParticle( const float&, const float&, const float&, const float&, const float& );
  ~SimpleParticle();
  
  TLorentzVector LorentzVector() const;
  float Charge() const;
  float Px()     const;   
  float Py()     const;
  float Pz()     const;
  float E()      const;
  float Et()     const;
  
  float Phi()       const;
  float Eta()       const;
  float Pt()        const;
  float Rapidity()  const;
  float M()         const;
  float Mass()      const;
  

  void SetLorentzVector(const TLorentzVector&);
  void SetCharge(const float&);
  void SetPtEtaPhiMass(const float&, const float&, const float&, const float&);

  void Print();
  
  SimpleParticle operator+(const SimpleParticle&);    
  
  
 private:

  float          m_charge;   
  TLorentzVector m_particle;

  
  

};


#endif

