#define TRK_CXX

#include "Trk.h"


//-- Constructor
Trk::Trk(){
  m_trk_d0       = 1;
  m_trk_z0       = 1;
  m_trk_theta    = 1;
  m_trk_phi      = 1;
  m_trk_q_over_p = 1;
  m_trk_pt       = -1;
  m_trk_eta      = 1;
   
}


//-- Destructor

Trk::~Trk(){};


//-- Getter Functions
float Trk::D0()     const { return m_trk_d0;       }
float Trk::Z0()     const { return m_trk_z0;       }
float Trk::Theta()  const { return m_trk_theta;    }
float Trk::Phi()    const { return m_trk_phi;      }
float Trk::QoverP() const { return m_trk_q_over_p; }
float Trk::Pt()     const { return m_trk_pt;       }
float Trk::Eta()    const { return m_trk_eta;      }


//-- Setter Functions

void Trk::SetD0( const float& d0 ){
  m_trk_d0 = d0;
}

void Trk::SetZ0( const float& z0 ){
  m_trk_z0 = z0 ;
}

void Trk::SetTheta( const float& theta ){
  m_trk_theta = theta;
}

void Trk::SetPhi( const float& phi ){
  m_trk_phi = phi;
}

void Trk::SetQoverP( const float& qoverp ){
  m_trk_q_over_p = qoverp;
}

void Trk::SetPt( const float& pt ){
  m_trk_pt = pt;
}

void Trk::SetEta( const float& eta ){
  m_trk_eta = eta;
}
