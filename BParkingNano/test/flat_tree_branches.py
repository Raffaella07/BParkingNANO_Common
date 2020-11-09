from collections import OrderedDict

branches = OrderedDict()


branches['run'           ]     = lambda x : x.run
branches['lumi'          ]     = lambda x : x.lumi
branches['event'         ]     = lambda x : x.event

branches['b_pt'          ]     = lambda x : x.b.pt()
branches['b_eta'         ]     = lambda x : x.b.eta()
branches['b_phi'         ]     = lambda x : x.b.phi()
branches['b_mass'        ]     = lambda x : x.b.mass()
branches['b_sv_chi2'     ]     = lambda x : x.b.sv_chi2()

branches['hnl_pt'        ]     = lambda x : x.b.hnl_fit_pt()
branches['hnl_eta'       ]     = lambda x : x.b.hnl_fit_eta()
branches['hnl_phi'       ]     = lambda x : x.b.hnl_fit_phi()
branches['hnl_mass'      ]     = lambda x : x.b.hnl_fit_mass()
branches['hnl_charge'    ]     = lambda x : x.b.hnl_charge()
branches['hnl_cos2D'     ]     = lambda x : x.b.hnl_fit_cos2D()

branches['mu_pt'         ]     = lambda x : x.b.mu().pt()
branches['mu_eta'        ]     = lambda x : x.b.mu().eta()
branches['mu_phi'        ]     = lambda x : x.b.mu().phi()
branches['mu_mass'       ]     = lambda x : x.b.mu().mass()
branches['mu_isSoft'     ]     = lambda x : x.b.sel_mu_isSoft()
branches['mu_isTight'    ]     = lambda x : x.b.sel_mu_isTight()
branches['mu_isMedium'   ]     = lambda x : x.b.sel_mu_isMedium()
branches['mu_isLoose'    ]     = lambda x : x.b.sel_mu_isLoose()
branches['mu_ip3d'       ]     = lambda x : x.b.sel_mu_ip3d()
branches['mu_sip3d'      ]     = lambda x : x.b.sel_mu_sip3d()
branches['mu_dxy'        ]     = lambda x : x.b.sel_mu_dxy()
branches['mu_dz'         ]     = lambda x : x.b.sel_mu_dz()

branches['pi_pt'         ]     = lambda x : x.b.pi().pt()
branches['pi_eta'        ]     = lambda x : x.b.pi().eta()
branches['pi_phi'        ]     = lambda x : x.b.pi().phi()
branches['pi_mass'       ]     = lambda x : x.b.pi().mass()
branches['pi_dz'         ]     = lambda x : x.b.pi_dz()
branches['pi_dxy'        ]     = lambda x : x.b.pi_dxy()
branches['pi_dzS'        ]     = lambda x : x.b.pi_dzS()
branches['pi_dxyS'       ]     = lambda x : x.b.pi_dxyS()
branches['pi_DCASig'     ]     = lambda x : x.b.pi_DCASig()

branches['trg_mu_pt'     ]     = lambda x : x.b.trg_mu().pt()
branches['trg_mu_eta'    ]     = lambda x : x.b.trg_mu().eta()
branches['trg_mu_phi'    ]     = lambda x : x.b.trg_mu().phi()
branches['trg_mu_mass'   ]     = lambda x : x.b.trg_mu().mass()
branches['trg_mu_ip3d'   ]     = lambda x : x.b.trg_mu_ip3d()
branches['trg_mu_sip3d'  ]     = lambda x : x.b.trg_mu_sip3d()
branches['trg_mu_dxy'    ]     = lambda x : x.b.trg_mu_dxy()
branches['trg_mu_dz'     ]     = lambda x : x.b.trg_mu_dz()

# difference of vertex positions of the two muons
branches['muons_vxdiff'  ]     = lambda x : x.b.mu_vxdiff()
branches['muons_vydiff'  ]     = lambda x : x.b.mu_vydiff()
branches['muons_vzdiff'  ]     = lambda x : x.b.mu_vzdiff()
branches['muons_Lxy'     ]     = lambda x : x.b.mu_Lxy()
branches['muons_Lxyz'    ]     = lambda x : x.b.mu_Lxyz()

# difference of vertex positions of the pion and trigger muon
branches['pi_trgmu_vzdiff']    = lambda x : x.b.pi_mu_vzdiff()

# dR
branches['dr_mu_pi'      ]     = lambda x : x.b.dr_mu_pi()
branches['dr_trgmu_hnl'  ]     = lambda x : x.b.dr_trgmu_hnl()




