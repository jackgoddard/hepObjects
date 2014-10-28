#ifndef MCBLOCK_H
#define MCBLOCK_H

#include <vector>
#include <TLorentzVector.h>

class MCBlock{
 public:
  MCBlock();
  ~MCBlock();

//-- Accessor Functions
  int   N()       const;
  float Pt()      const;
  float M()       const;
  float Eta()     const; 
  float Phi()     const;
  float Charge()  const;
  int   Status()  const;
  int   BarCode() const;
  int   PDGID()   const;

  int   Type() const;

  std::vector<int>  Parents()    const;
  std::vector<int>  Children()   const;

  std::vector<int> ParentIndex() const;
  std::vector<int> ChildIndex()  const;

  TLorentzVector LorentzVector() const;

//-- Setter Functions

  void SetN(       const int   );
  void SetPt(      const float );
  void SetM(       const float );
  void SetEta(     const float );
  void SetPhi(     const float );
  void SetCharge(  const float );
  void SetStatus(  const int   );
  void SetBarCode( const int   );
  void SetPDGID(   const int   );

  void SetParents(     const std::vector<int> );
  void SetChildren(    const std::vector<int> );
  void SetParentIndex( const std::vector<int> );
  void SetChildIndex(  const std::vector<int> );

  void SetLorentzVector( const TLorentzVector );

 private:


  int   m_n;
  float m_pt;
  float m_m;
  float m_eta;
  float m_phi;
  float m_charge;
  int   m_status;
  int   m_barcode;
  int   m_pdgid; 


  std::vector<int> m_parent_index;
  std::vector<int> m_parents;
  
  std::vector<int> m_child_index;
  std::vector<int> m_children;

  TLorentzVector m_lorentz_vector;

};

#endif
