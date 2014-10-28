#define Vertex_cxx

////
//
//  Vertex Class for the Low Mass Drell-Yan Analysis using D3PDs
//
//
////

#include <iostream>
#include <vector>
#include "Vertex.h"
//#include "Particle.h"



//--- Constructors

Vertex::Vertex(){

  m_x = 0.0;
  m_y = 0.0;
  m_z = 0.0;

}

Vertex::Vertex( const float& pos_x, const float& pos_y, const float& pos_z ){

  m_x = pos_x;
  m_y = pos_y;
  m_z = pos_z;

}


//--- Destructor

Vertex::~Vertex(){}


//--- Getter Functions

float Vertex::X() const { return m_x; }
float Vertex::Y() const { return m_y; }
float Vertex::Z() const { return m_z; }

float Vertex::XError() const { return m_x_error; }
float Vertex::YError() const { return m_y_error; }
float Vertex::ZError() const { return m_z_error; }

float Vertex::Chi2()  const { return m_chi2; }
int   Vertex::NDof()  const { return m_dof;  }

int   Vertex::NumberOfTrks() const { return m_numberoftracks; }

float Vertex::SumPt() const { return m_sum_pt;   }

int   Vertex::Type()  const { return m_vtx_type; }

std::vector<int> Vertex::VertexTrkIndex() const { return m_vtx_trk_index; }


//--- Setter Functions


void Vertex::SetX(const float& pos_x ){
  m_x = pos_x;
}

void Vertex::SetY(const float& pos_y ){
  m_y = pos_y;
}

void Vertex::SetZ(const float& pos_z ){
  m_z = pos_z;
}


void Vertex::SetXError(const float& pos_x_error ){
  m_x_error = pos_x_error;
}

void Vertex::SetYError(const float& pos_y_error ){
  m_y_error = pos_y_error;
}

void Vertex::SetZError(const float& pos_z_error ){
  m_z_error = pos_z_error;
}


void Vertex::SetChi2(const float& chi2){
  m_chi2 = chi2;
}

void Vertex::SetNDof(const int& dof ){
  m_dof = dof;
}


void Vertex::SetNumberOfTrks( const int& numtracks ){
  m_numberoftracks = numtracks;
}


void Vertex::SetSumPt( const float& sum_pt ){
  m_sum_pt = sum_pt;
}


void Vertex::SetType( const int& type ){
  m_vtx_type = type;
}


void Vertex::SetVertexTrkIndex( const std::vector<int>& index ){
  m_vtx_trk_index = index;
}
