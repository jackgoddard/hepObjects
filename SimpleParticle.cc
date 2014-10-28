//
//  Simple Particle Class
//
//
//  Jack Goddard, February 2012
//


#include <iostream>
#include "SimpleParticle.h"

using namespace std;


//-- Constructors:

SimpleParticle::SimpleParticle(){
  m_particle = TLorentzVector( 0, 0, 0, 0);
  m_charge   = 100000;
}

SimpleParticle::SimpleParticle( const TLorentzVector& particle, const float& charge){
  m_particle = particle;
  m_charge   = charge;
}


SimpleParticle::SimpleParticle( const float& part_px, const float& part_py, const float& part_pz, const float& part_E, const float& charge){
  m_particle = TLorentzVector( part_px, part_py, part_pz, part_E );
  m_charge   = charge;
}


//-- Destructor

SimpleParticle::~SimpleParticle() {}


//-- Getter Functions

TLorentzVector SimpleParticle::LorentzVector() const {return m_particle; }

float SimpleParticle::Pt()        const { return m_particle.Pt();       }
float SimpleParticle::Eta()       const { return m_particle.Eta();      }
float SimpleParticle::Phi()       const { return m_particle.Phi();      }
float SimpleParticle::Rapidity()  const { return m_particle.Rapidity(); }
float SimpleParticle::M()         const { return m_particle.M();        } 
float SimpleParticle::Mass()      const { return m_particle.M();        }

float SimpleParticle::Px() const { return m_particle.Px(); }
float SimpleParticle::Py() const { return m_particle.Py(); }
float SimpleParticle::Pz() const { return m_particle.Pz(); }
float SimpleParticle::E()  const { return m_particle.E();  }
float SimpleParticle::Et() const { return m_particle.Et(); }

float SimpleParticle::Charge() const { return m_charge; }



//--- Setter Functions  

void SimpleParticle::SetCharge(const float& charge){
  m_charge = charge;
}

void SimpleParticle::SetLorentzVector(const TLorentzVector& particleVec){
  m_particle = particleVec;
}

void SimpleParticle::SetPtEtaPhiMass(const float& pt, const float& eta, const float& phi, const float& m){
  m_particle.SetPtEtaPhiM( pt, eta, phi, m );
}


void SimpleParticle::Print(){
  std::cout<< this->Px() <<"\t"<< this->Py() << "\t"  << this->Pz() << "\t" << this->E() << "\t" << this->Charge() <<std::endl;
}


SimpleParticle SimpleParticle::operator+(const SimpleParticle& b){

  TLorentzVector combinedVec = this->m_particle + b.LorentzVector();
  
  float combinedCharge = this->m_charge + b.Charge();

  SimpleParticle combined(combinedVec, combinedCharge);

  return combined;

}
