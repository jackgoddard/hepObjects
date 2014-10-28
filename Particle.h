//  Particle Class
//  
//  Jack Goddard, February 2011
//

#include <TROOT.h>
#include <iostream>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "SimpleParticle.h"

#ifndef PARTICLE_H
#define PARTICLE_H


class Particle{
   public:
	Particle();
	Particle( const TLorentzVector particle, const float charge );
	Particle( const float part_px, const float part_py, const float part_pz, const float part_E, const float charge );
	~Particle();
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
	int   Author()    const;
	float SmearedPt() const;
	

	float PtRatio20() const;
	float PtRatio30() const;
	float PtRatio40() const;
	float EtRatio20() const; 
	float EtRatio30() const; 
	float EtRatio40() const; 
	
	float PtCone20() const;
	float PtCone30() const;
	float PtCone40() const;
	float EtCone20() const;
	float EtCone30() const;
	float EtCone40() const;
	float EtCore()   const;
	

	float D0()         const;
	float Z0()         const;
	float D0Error()    const;
	float Z0Error()    const;
	float CovD0Z0()    const;
	float TrkFitChi2() const;

	bool PassesLoose()  const;
	bool PassesMedium() const;
	bool PassesTight()  const;
	bool IsCombined()   const;
	bool PassRel17MCPTrackQuality() const;

	int   MCTruthBarCode()       const;
	int   MCTruthMatched()       const;
	float MCTruthE()             const;
	float MCTruthDeltaR()        const;
	float MCTruthEta()           const;
	float MCTruthPhi()           const;
	float MCTruthPt()            const;
	int   MCTruthStatus()        const;
	int   MCTruthType()          const;
	int   MCTruthMotherType()    const;
	int   MCTruthMotherBarCode() const;


	int NumberOfBLayerHits() const;
	int NumberOfPixelHits()  const;
	int NumberOfSCTHits()    const;
	int NumberOfTRTHits()    const;

	int ExpectBLayerHit()          const;
	int NumberOfSCTDeadSensors()   const;
	int NumberOfPixelDeadSensors() const;
	int NumberOfSCTHoles()         const;
	int NumberOfPixelHoles()       const;
	int NumberOfTRTOutliers()      const;

//---- Component Track Varibles
	float IDQoverP()     const;
	float IDTheta()      const;
	float IDPhi()        const;
	float IDD0()         const;
	float IDZ0()         const;
	float IDD0Error()    const;
	float IDZ0Error()    const;
	float IDD0Significance() const;
	float IDZ0Significance() const;
	float IDCharge() const;
	float IDUnSmearedPt() const;
	TVector3 IDVector3() const;

	float CBQoverP()     const;
	float CBTheta()      const;
	float CBPhi()        const;
	float CBD0()         const;
	float CBZ0()         const;
	TVector3 CBVector3() const;

	float MSQoverP()     const;
	float MSTheta()      const;
	float MSPhi()        const;
	float MSD0()         const;
	float MSZ0()         const;
	TVector3 MSVector3() const;

	float MEQoverP()     const;
	float METheta()      const;
	float MEPhi()        const;
	float MED0()         const;
	float MEZ0()         const;
	TVector3 MEVector3() const;

	float IEQoverP()     const;
	float IETheta()      const;
	float IEPhi()        const;
	float IED0()         const; 
	float IEZ0()         const;
	TVector3 IEVector3() const;

	double SmearedChargeFlip()  const;
	float  UnSmearedPt()        const;
	float  MCRecSmearedPt()   const;
	float  MCRecIESmearedPt() const;
	float  MCRecMESmearedPt() const;

//----- Setter Functions

	void SetLorentzVector(const TLorentzVector);
	void SetCharge(const float);
	void SetPtEtaPhiMass(const float, const float, const float, const float);
	void SetPtCone20(const float);
	void SetPtCone30(const float);
	void SetPtCone40(const float);
	void SetEtCone20(const float);
	void SetEtCone30(const float);
	void SetEtCone40(const float);
	void SetEtCore(const float);
	void SetD0(const float);
	void SetZ0(const float);
	void SetD0Error(const float);
	void SetZ0Error(const float);
	void SetCovD0Z0(const float);
	void SetTrkFitChi2( const float);

	void SetPassesLoose(const int);
	void SetPassesMedium(const int);
	void SetPassesTight(const int);
	void SetIsCombined( const int );
	
	void SetAuthor(const int);

	void  SetMCTruthBarCode(const int);
	void  SetMCTruthMatched(const int);
	void  SetMCTruthE(const float);
	void  SetMCTruthDeltaR(const float);
	void  SetMCTruthEta(const float);
	void  SetMCTruthPhi(const float);
	void  SetMCTruthPt(const float);
	void  SetMCTruthStatus(const int);
	void  SetMCTruthType(const int);
	void  SetMCTruthMotherType(const int);
	void  SetMCTruthMotherBarCode(const int);

//----- Refitted Track Variables
	void SetNumberOfBLayerHits(const int);
	void SetNumberOfPixelHits(const int);
	void SetNumberOfSCTHits(const int);
	void SetNumberOfTRTHits(const int);

        void SetExpectBLayerHit(const int);
        void SetNumberOfSCTDeadSensors(const int);
        void SetNumberOfPixelDeadSensors(const int);
        void SetNumberOfSCTHoles( const int);
	void SetNumberOfPixelHoles( const int);
        void SetNumberOfTRTOutliers( const int);

//---- Component Track Varibles
//	void SetIDQoverP( const float);
//	void SetIDTheta( const float );
//	void SetIDPhi( const float );
	void SetIDD0( const float );
	void SetIDZ0( const float );
	void SetIDD0Error( const float );
	void SetIDZ0Error( const float );
	void SetIDVector3( const TVector3 );
	void SetIDCharge( const float );
	void SetIDUnSmearedPt( const float );

	void SetCBQoverP( const float);
	void SetCBTheta( const float );
	void SetCBPhi( const float );
	void SetCBD0( const float );
	void SetCBZ0( const float );

	void SetMSQoverP( const float);
	void SetMSTheta( const float );
	void SetMSPhi( const float );
	void SetMSD0( const float );
	void SetMSZ0( const float );

	void SetMEQoverP( const float);
	void SetMETheta( const float );
	void SetMEPhi( const float );
	void SetMED0( const float );
	void SetMEZ0( const float );

	void SetIEQoverP( const float);
	void SetIETheta( const float );
	void SetIEPhi( const float );
	void SetIED0( const float );
	void SetIEZ0( const float );

	void SetSmearedChargeFlip( const double );
	void SetUnSmearedPt( const float );
	void SetMCRecSmearedPt( const float );
	void SetMCRecIESmearedPt( const float );
	void SetMCRecMESmearedPt( const float );

	void Print();

	SimpleParticle operator+(const Particle&);	

   private:

	TLorentzVector m_particle;
	float m_charge;	
	float m_ptcone20;
	float m_ptcone30;
	float m_ptcone40;	
	float m_etcone20;
	float m_etcone30;
	float m_etcone40;
	float m_et_core;

	float m_d0;
	float m_z0;
	float m_d0_err;
	float m_z0_err;
	float m_cov_d0z0;
	float m_trk_fit_chi2;

	bool m_passesloose;
	bool m_passesmedium;
	bool m_passestight;
	bool m_is_combined;	

	int m_author;

	int m_num_hits_blayer;
	int m_num_hits_pixel;
	int m_num_hits_sct;
	int m_num_hits_trt;

	int m_expect_b_layer_hit;
	int m_num_pixel_dead_sensors;
	int m_num_sct_dead_sensors;
	int m_num_pixel_holes;
	int m_num_sct_holes;
	int m_num_trt_outliers;


	float m_id_q_over_p;
	float m_id_theta;
	float m_id_phi;
	float m_id_d0;
	float m_id_z0;
	float m_id_d0_err;
	float m_id_z0_err;
	TVector3 m_id_vector;
	float    m_id_charge;
	float    m_id_unsmeared_pt;
	int      m_id_smeared_charge_flip;
	
	float m_cb_q_over_p;
	float m_cb_theta;
	float m_cb_phi;
	float m_cb_d0;
	float m_cb_z0;

	float m_ms_q_over_p;
	float m_ms_theta;
	float m_ms_phi;
	float m_ms_d0;
	float m_ms_z0;


	float m_me_q_over_p;
	float m_me_theta;
	float m_me_phi;
	float m_me_d0;
	float m_me_z0;


	float m_ie_q_over_p;
	float m_ie_theta;
	float m_ie_phi;
	float m_ie_d0;
	float m_ie_z0;

	
	float  m_mc_rec_smeared_pt;
	float  m_mc_rec_ie_smeared_pt;
	float  m_mc_rec_me_smeared_pt;
	double m_smeared_charge_flip;
	float  m_unsmeared_pt;

//----- Monte Carlo Truth Variables:

	int    m_mc_truth_barcode;
	int    m_mc_truth_matched;
	float  m_mc_truth_e;
	float  m_mc_truth_delta_r;
	float  m_mc_truth_eta;
	float  m_mc_truth_phi;
	float  m_mc_truth_pt;
	int    m_mc_truth_status;
	int    m_mc_truth_type;
       	int    m_mc_truth_mother_type;
	int    m_mc_truth_mother_barcode;

};


#endif
