////
//
//  Vertex Class Header File for the Low Mass Drell-Yan Analysis using D3PDs
//  
//
////

#ifndef VERTEX_H
#define VERTEX_H


class Vertex{
 public:

//-- Constructors
  Vertex();
  Vertex(const float&, const float&, const float&);

//-- Destructor
  ~Vertex();


//-- Getter Functions

  float X() const;
  float Y() const;
  float Z() const;

  float XError() const;
  float YError() const;
  float ZError() const;

  float Chi2() const;
  int   NDof() const;

  int   NumberOfTrks() const;
  float SumPt() const;
  
  int   Type() const;

  std::vector<int>  VertexTrkIndex() const;
  

//-- Setter Functions

  void SetX(const float&);
  void SetY(const float&);
  void SetZ(const float&);

  void SetXError(const float&);
  void SetYError(const float&);
  void SetZError(const float&);


  void SetChi2(const float&);
  void SetNDof(const int& );

  void SetNumberOfTrks(const int&);
  void SetSumPt(const float&);

  void SetType(const int&);
  
  void SetVertexTrkIndex(const std::vector<int>&);

 private:

  float m_x;
  float m_y;
  float m_z;

  float m_x_error;
  float m_y_error;
  float m_z_error;

  float m_chi2;
  int   m_dof;

  float m_sum_pt;

  int   m_numberoftracks;
  int   m_vtx_type;

  std::vector<int> m_vtx_trk_index;

};

#endif
