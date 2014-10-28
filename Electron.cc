//
//  Electron Class
// 
//  Inherits from the SimpleParticle Class
//
//  Jack Goddard
//
//  Created: 16/07/2012
//


#include <iostream>
#include "Electron.h"
#include "SimpleParticle.h"

using namespace std;


Electron::Electron(){
  SimpleParticle();
  m_d0 = 0.0;
  m_z0 = 0.0;
};

Electron::Electron( const TLorentzVector& particle, const float& charge){
  SimpleParticle( particle, charge );
  m_d0 = 0.0;
  m_z0 = 0.0;
};

Electron::~Electron(){};

float Electron::D0()     const { return m_d0; }
float Electron::Z0()     const { return m_z0; }

void Electron::SetD0( const float& d0 ){ m_d0 = d0; }
void Electron::SetZ0( const float& z0 ){ m_z0 = z0; }


