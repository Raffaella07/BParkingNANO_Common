from collections import OrderedDict

branches = OrderedDict()


branches['run'  ] = lambda x : x.run
branches['lumi' ] = lambda x : x.lumi
branches['event'] = lambda x : x.event


branches['b_pt'  ] = lambda x : x.b.pt()
branches['b_eta' ] = lambda x : x.b.eta()
branches['b_phi' ] = lambda x : x.b.phi()
branches['b_mass'] = lambda x : x.b.mass()

branches['hnl_pt'  ] = lambda x : x.b.hnl_fit_pt()
branches['hnl_eta' ] = lambda x : x.b.hnl_fit_eta()
branches['hnl_phi' ] = lambda x : x.b.hnl_fit_phi()
branches['hnl_mass'] = lambda x : x.b.hnl_fit_mass()

branches['mu_pt'  ] = lambda x : x.b.mu().pt()
branches['mu_eta' ] = lambda x : x.b.mu().eta()
branches['mu_phi' ] = lambda x : x.b.mu().phi()
branches['mu_mass'] = lambda x : x.b.mu().mass()

branches['pi_pt'  ] = lambda x : x.b.pi().pt()
branches['pi_eta' ] = lambda x : x.b.pi().eta()
branches['pi_phi' ] = lambda x : x.b.pi().phi()
branches['pi_mass'] = lambda x : x.b.pi().mass()

branches['trg_mu_pt'  ] = lambda x : x.b.trg_mu().pt()
branches['trg_mu_eta' ] = lambda x : x.b.trg_mu().eta()
branches['trg_mu_phi' ] = lambda x : x.b.trg_mu().phi()
branches['trg_mu_mass'] = lambda x : x.b.trg_mu().mass()
