//
//  Electron Class
// 
//  Inherits from the SimpleParticle Class
//
//  Jack Goddard
//
//  Created: 16/07/2012
//

#ifndef ELECTRON_H
#define ELECTRON_H

#include "SimpleParticle.h"

class Electron : public SimpleParticle{
  
 public:
  Electron();
  Electron( const TLorentzVector& particle, const float& charge );
  ~Electron();
  
  float D0() const;
  float Z0() const;


  void SetD0( const float& d0 );
  void SetZ0( const float& z0 );
  



  
  

 private:
  double m_d0;
  double m_z0;

};
#endif
