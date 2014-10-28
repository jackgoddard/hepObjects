#ifndef MCMUONTRUTH_H
#define MCMUONTRUTH_H

#include <vector>
#include <TLorentzVector.h>

class MCMuonTruth{
 public:
  MCMuonTruth();
  ~MCMuonTruth();

//-- Accessor Functions
  int   N()       const;
  float Pt()      const;
  float M()       const;
  float Eta()     const; 
  float Phi()     const;
  float Charge()  const;
  int   Origin()  const;
  int   BarCode() const;
  int   PDGID()   const;
  int   Type()    const;

  std::vector<int>  Parents()  const;
  std::vector<int>  Children() const;

  TLorentzVector LorentzVector() const;

//-- Setter Functions

  void SetN(       const int   );
  void SetPt(      const float );
  void SetM(       const float );
  void SetEta(     const float );
  void SetPhi(     const float );
  void SetCharge(  const float );
  void SetOrigin(  const int   );
  void SetBarCode( const int   );
  void SetPDGID(   const int   );
  void SetType(    const int   );

  void SetLorentzVector( const TLorentzVector );

 private:


  int   m_n;
  float m_pt;
  float m_m;
  float m_eta;
  float m_phi;
  float m_charge;
  int   m_origin;
  int   m_barcode;
  int   m_pdgid; 
  int   m_type;

  TLorentzVector m_lorentz_vector;

};

#endif
