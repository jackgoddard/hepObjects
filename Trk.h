#ifndef TRK_H
#define TRK_H


class Trk{

 public:

//-- Constructors
  Trk();

//-- Destructor
  ~Trk();


//-- Getter Functions
  float D0()     const;
  float Z0()     const;
  float Theta()  const;
  float Phi()    const;
  float QoverP() const;
  float Pt()     const;
  float Eta()    const;



//-- Setter Functions

  void SetD0(     const float& );
  void SetZ0(     const float& );
  void SetTheta(  const float& );
  void SetPhi(    const float& );
  void SetQoverP( const float& );
  void SetPt(     const float& );
  void SetEta(    const float& );


 private:

  float m_trk_d0;
  float m_trk_z0;
  float m_trk_theta;
  float m_trk_phi;
  float m_trk_q_over_p;
  float m_trk_pt;
  float m_trk_eta;



};
#endif
