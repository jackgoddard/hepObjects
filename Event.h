//
//  Event Class
//
//  This is the class in which all the hard work is done, objects of other 
//  classes are contained here and pre-selections, selections and calculations
//  are done by the member functions of this class.
//

#include "SimpleParticle.h"
#include "Particle.h"
#include "Vertex.h"
#include "Trk.h"
#include "MCBlock.h"
#include "MCMuonTruth.h"
#include "TH1D.h"
#include "../PhysicsAnalysis.h"


#ifndef EVENT_H
#define EVENT_H

class Event{
 public:
  Event();
  Event(int muons_n, int vpx_n, bool trigger);
  ~Event();

  bool IsVerbose() const;
 
//-- Getter Functions

  int    KFactorMode() const;
  int    NumberOfMuons() const;
  int    NumberOfTruthMuons() const;
  int    NumberOfMCBlockParticles() const;
  int    NumberOfPrimaryVertices()  const;
  int    NumberOfTrks() const;
  bool   TriggerFired() const;
  int    EventNumber()  const;
  int    RunNumber()         const;
  int    MCSimPeriod()       const;
  float  AverageIntPerXing() const;
  unsigned int LArError()    const;
  bool   WantPtSmearedMuons() const;

  Particle       Muon1() const;
  Particle       Muon2() const;
  SimpleParticle Propagator() const;
  Particle       Muon() const;
  Particle       AntiMuon() const;
  

  Particle       MCMuon1() const;
  Particle       MCMuon2() const;
  SimpleParticle MCPropagator() const;

  Vertex ZerothVertex() const;

  std::vector<Vertex> Vertices() const;
  std::vector<Particle> Muons() const; 
  std::vector<Particle> MuonCandidates() const;
  std::vector<Trk> Trks() const;

  std::vector<MCMuonTruth> MCMuonsTruth() const;
  std::vector<MCBlock> MCBlockParticles() const;

  double MuonDeltaZ0()  const;
  double MuonDeltaPhi() const;
  double MuonDeltaEta() const;
  double MuonDeltaR()   const;
  double MuonIDDeltaR() const;

  float MCPdfScale() const;
  float MCPdf1()     const;
  float MCPdf2()     const;
  float MCPdfX1()    const;
  float MCPdfX2()    const;
  int   MCPdfId1()   const;
  int   MCPdfId2()   const;
  
  double KFactorWeight()                const;
  double LumiWeight()                   const;
  double VertexWeight()                 const;
  double MCEventWeight()                const;
  double MuonMCEfficiencyWeight()       const;
  double EventWeight()                  const;
  double PileupWeight()                 const;
  double kGridfn( double, double )      const;
  double kGridfn( double, double, int ) const;
  double TriggerEfficiencyWeight()      const;
  double IsolationEfficiencyWeight()    const;
  double VertexPositionWeight()         const;
  double PtReweight()                   const;
  
  double Muon1RecoSf()      const;
  double Muon2RecoSf()      const;
  double Muon1RecoSfError() const;
  double Muon2RecoSfError() const;  

  vector<double> PrimaryVertexDistributionWeights() const;

  std::map<TString, bool> SystematicVariations() const;

  void PrintWeights()   const;
  void PrintEventInfo() const;
  void SetIsVerbose( const bool );
  
  void DefineSelection();
  bool MakeZSelection();
  bool MakeDYSelection();

  //  void TwoHighestPtMuons( std::vector<Particle>, Particle&, Particle&);
  void TwoHighestPtMuons( std::vector<Particle>);
  void TwoHighestPtMuonsRegardlessOfCharge( const std::vector<Particle>);
  bool EventHasPositiveAndNegativeMuons( const std::vector<Particle> );

  bool BothMuonsPassEtRatioCut( const int, const double);
  bool BothMuonsPassPtRatioCut( const int, const double);
  bool BothMuonsFailEtRatioCut( const int, const double);
  bool BothMuonsFailPtRatioCut( const int, const double);

  bool BothMuonsPassEtaCut( const double );
  bool BothMuonsPassPtCut(  const double );
  bool BothMuonsPassPtCut(  const double, const double );
  bool BothMuonsPassD0Cut(  const double );
  bool BothMuonsPassZ0Cut(  const double );
  bool BothMuonsPassD0SignificanceCut(  const double );
  bool PropagatorPassesMassCut( const double, const double);
  bool PropagatorPassesRapidityCut( const double );
  bool BothMuonsPassPtSumCut( const int, const double );
  bool BothMuonsPassEtSumCut( const int, const double );
  bool BothMuonsFailPtSumCut( const int, const double );
  bool BothMuonsFailEtSumCut( const int, const double );
  bool MuonsPassDeltaZ0Cut(  const double );
  bool BothMuonsPassMSPtCut( const double );
  bool BothMuonsPassMSIDPtDiffCut( const double );
  bool BothMuonsPassZ0ZVtxDiffCut( const double );
  bool MuonPassesRel16MCPTrackQualityCuts( const Particle );
  bool MuonPassesRel17MCPTrackQualityCuts( const Particle );
  bool MuonInLArHole();

  bool BothMuonsPassRel16MCPTrackQualityCuts();
  bool BothMuonsPassRel17MCPTrackQualityCuts();
  void UseMCPtSmearedMuons();


  void CheckVertexTrackIndexRelation();

  void SetNumberOfMuons(const int);

  void SetNumberOfPrimaryVertices( const int );
  void SetNumberOfTrks( const int );
  void SetIfTriggerFired( const bool );
  void SetEventNumber( const int );
  void SetRunNumber( const int );
  void SetMCSimPeriod( const int );
  void SetAverageIntPerXing( const float );

  void SetMuon1(      const Particle       );
  void SetMuon2(      const Particle       );
  void SetPropagator( const SimpleParticle );
  void SetMuon(       const Particle       );
  void SetAntiMuon(   const Particle       );

  void SetMuonAndAntiMuon();

  // void SetMCMuon1(Particle);
  //void SetMCMuon2(Particle);
  // void SetMCPropagator(Particle);

  void SetMCPdfScale(const float);
  void SetMCPdf1(    const float);
  void SetMCPdf2(    const float);
  void SetMCPdfX1(   const float);
  void SetMCPdfX2(   const float);
  void SetMCPdfId1(  const int  );
  void SetMCPdfId2(  const int  );

  
  void SetVertices( const std::vector<Vertex>& );
  void SetTrks(     const std::vector<Trk>    );
  void SetZerothVertex();

  void CorrectForOverLappingCones();
  void CheckSystematicVariations();
  void SetSystematicVariations( const std::map< TString, bool>& );
 
  void SetKFactorWeight( const int );
  void SetLumiWeight(const double, const double);
  void SetLumiWeight(const double, const double, const double);
  void SetVertexWeight( const TH1D* );
  void SetMCEventWeight( const double );
  void SetMuonMCEfficiencyWeight(Analysis::AnalysisMuonEfficiencyScaleFactors* );
  void SetPileupWeight(const double);
  void SetTriggerEfficiencyWeight( const double sf ); 
  void SetTriggerEfficiencyWeight( const std::map<TString, TH2F*>); 
  void SetIsolationEfficiencyWeight( int sampleID, const std::map<TString, TH1D*>);
  void SetVertexPositionWeight( const float& );
  void SetPtReweight( const double, map<TString,TH1D*> );
  
  double CalculatePtReweight( const double );

  void SetPrimaryVertexDistributionWeights( const TH1D*, const TH1D*, TH1D* );

  void ApplyMuonMCEfficiencyScaleFactors();

  void SetLArError( const unsigned int );

  void SetKFactorMode( const int );

  void SetMuons( const std::vector<Particle>& );

  void SetDYMuonCandidates( const std::vector<Particle>& );
  void SetDYMuonCandidatesRegardlessOfCharge( const std::vector<Particle>& );
  
  void PreSelectionD0Cut(std::vector<Particle>&, double);
  void PreSelectionMuonEtRatioCut(std::vector<Particle>&, double);
  void PreSelectionMuonPtRatioCut(std::vector<Particle>&, double);
  
  void SetMakeZSelection(bool);
  void SetMakeDYSelection(bool);
  
  bool PassesMuonPreSelection( const int, const double, const double);
  bool PassesPVWithAssociatedTracksCut( const std::vector<Vertex>, int, int);
  bool BothMuonsFromSameVertex( const double );

  
  bool PassesNonIsolatedSelection();
  bool PassesEventSelection();
  bool PassesPreSelection();
  bool PassesSelectionWithoutIso();
  bool PassesIsolationSelection();
  void CreateEventSelectionCutFlow(TH1D*,  const double );
  void CreatePreSelectionCutFlow(TH1D*,  const double );
  void CreateSelectionCutFlow(TH1D*,     const double );
  void CreateFullSelectionCutFlow(TH1D*,     const double );
  void CreateZEventSelectionCutFlow(TH1D*, const double );
  void CreateZPreSelectionCutFlow(TH1D*, const double );
  void CreateZSelectionCutFlow(TH1D*,    const double );
  void CreateCutList( TProfile* );
  void SetMCPropagatorCandidates(const std::vector<SimpleParticle>);

  void SetMCMuonsTruth( const std::vector<MCMuonTruth>& );
  void SetMCBlockParticles( const std::vector<MCBlock>& );

//-- Variable Declarations: -----------------------------------------------
//-- MC/Weightings
  bool  IsMCEvent;
  bool  UseKFactorWeighting;
  bool  UseVertexWeighting;
  bool  UseLumiWeighting;
  bool  UseMCEventWeighting;
  bool  UseMuonEfficiencyWeighting;
  bool  UsePileupWeighting;
  bool  AnalyseEntireMCBlock; 
  bool  UsePtSmearedMuonsInAnalysis;
  bool  UseTriggerEfficiencyWeighting;
  bool  UseIsolationEfficiencyWeighting;
  bool  UseVertexPositionWeighting;
  bool  UsePtReweighting;

  double DataLuminosity;
  double MCCrossSection;
  int    MCSampleID;


//-- PreSelection
  int    PreSelection_n_vertices;
  int    PreSelection_n_tracks_to_vertex;
  int    PreSelection_n_leptons;
  bool   PreSelection_reject_high_lumi_mc;
  double PreSelection_lepton_min_pt;
  double PreSelection_lepton_max_z0; 
  double PreSelection_lepton_max_d0; 
  double PreSelection_lepton_max_eta; 
  double PreSelection_lepton_max_deltaz0; 
  bool   PreSelection_muons_from_same_vertex;
  double PreSelection_vxp_chi2_cut;
  double PreSelection_max_vtx_z;
  bool   PreSelection_require_trigger;
  double PreSelection_min_MS_pt;
  double PreSelection_min_MS_ID_pt_diff;
  double PreSelection_min_z0_zvtx_diff;
  bool   PreSelection_reject_lar_errors;
  bool   PreSelection_reject_lar_hole;


  //-- Selection
  bool   Selection_require_trigger;
  int    Selection_n_leptons;
  double Selection_muon1_min_pt;
  double Selection_muon2_min_pt;
  double Selection_lepton_max_z0;
  double Selection_lepton_max_d0;
  double Selection_lepton_max_d0_sig;
  double Selection_lepton_max_eta;
  double Selection_min_mass;
  double Selection_max_mass;
  double Selection_max_prop_rap;
  bool   Selection_iso_et_sum;
  bool   Selection_iso_pt_sum;
  bool   Selection_iso_et_ratio;
  bool   Selection_iso_pt_ratio;
  int    Selection_et_cone;
  double Selection_et_cone_ratio_max;
  int    Selection_pt_cone;
  double Selection_pt_cone_ratio_max;
  double Selection_et_sum_max;
  double Selection_pt_sum_max;
  bool   Selection_require_opposite_charge_leptons;
  bool   Selection_require_MCPRel16_track_quality_cuts;
  bool   Selection_require_MCPRel17_track_quality_cuts;
 
  


//-------------------------------------------------------------------------
 private:

  int   m_kfactor_mode;
  
  int   m_vxp_n;
  int   m_muon_n;
  bool  m_trigger;
  int   m_event_number;
  int   m_run_number;
  int   m_mc_sim_period;
  int   m_trk_n;
  float m_av_int_per_xing;

  unsigned int m_lar_error;

  bool m_make_Z_selection;
  bool m_make_DY_selection;
  bool m_is_verbose;


  Particle       m_muon1;
  Particle       m_muon2;
  SimpleParticle m_propagator;
  Particle       m_muon;
  Particle       m_antimuon;
  Vertex         m_zerothvertex;

  std::vector<Particle> m_muons;

  std::vector<Particle> m_muon_candidates;
  std::vector<Particle> m_muon_candidates_regardless_of_charge;

  std::vector<Particle> m_muon_preselected_candidates;

  std::vector<Vertex>   m_vertices;

  std::vector<Trk>      m_trks;
 
  std::vector<double>   m_primary_vertex_distribution_weights;
 
  double m_k_factor_weight;
  double m_lumi_weight;
  double m_vertex_weight;
  double m_mcevent_weight;
  double m_mc_muon_eff_weight;
  double m_pileup_weight;
  double m_trigger_efficiency_weight;
  double m_isolation_efficiency_weight;
  double m_vertex_position_weight;
  double m_mc_pt_reweight;
  
  
  double m_muon1_reco_sf;
  double m_muon2_reco_sf;
  double m_muon1_reco_sf_error;
  double m_muon2_reco_sf_error;


//--- PDF Info
  float m_mc_pdf1;
  float m_mc_pdf2;
  float m_mc_pdf_scale;
  float m_mc_pdf_x1;
  float m_mc_pdf_x2;
  int   m_mc_pdf_id1;
  int   m_mc_pdf_id2;


  std::vector<MCBlock>     m_mc_block_particles;
  std::vector<MCMuonTruth> m_mc_muons_truth;


  std::map< TString, bool > m_syst_variations;


};

#endif
