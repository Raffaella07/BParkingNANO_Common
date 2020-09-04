# you need to install these two packages
# python -m pip install uproot4  --user
# python -m pip install awkward1 --user


# stuff that I've read while developing this
# https://medium.com/@mgarod/dynamically-add-a-method-to-a-class-in-python-c49204b85bd6



import hashlib
import awkward1 as ak
import numpy as np
import uproot
import uproot4
import ROOT

base_branches = ['pt', 'eta', 'phi', 'mass', 'pdgId']
extra_branches = ['charge', 'dxy', 'dz']
gen_branches = ['status', 'statusFlags', 'genPartIdxMother']
bcand_branches = ['pt', 'eta', 'phi', 'mass', 'sv_chi2', 'sv_prob', 'sv_lxy', 'sv_lxye', 'sv_lxy_sig', 'sv_x', 'sv_y', 'sv_z', 'sv_xe', 'sv_ye', 'sv_ze', 'hnl_cos2D', 'hnl_fit_cos2D', 'hnl_fit_mass', 'hnl_fit_masse', 'hnl_fit_pt', 'hnl_fit_eta', 'hnl_fit_phi', 'hnl_charge']

class Particle(object):

    def __new__(self, *args, **kwargs):
        if len(args)>=4: 
            idx = args[3]
        elif 'idx' in kwargs.keys():
            idx = kwargs['idx']
        else:
            return None        
        if idx<0:
            return None           
        return super(Particle, self).__new__(self)

    def __init__(self, tree, collection, ievent, idx, branches=base_branches):
        
        self.collection = collection    
        self.ievent     = ievent    
        self.idx        = idx       
        
        self.__set_values(tree, branches)
    
    def __set_values(self, tree, branches):
        for ibranch in branches:
            # prepend double underscore to make these attributes 'invisible', sort of private
            setattr(self, '__'+ibranch, tree['_'.join([self.collection, ibranch])].array()[self.ievent][self.idx])

            # function factory trick explained here
            # https://stackoverflow.com/questions/3431676/creating-functions-in-a-loop 
            def function_factory(ibranch=ibranch):
                def f():
                    return getattr(self, '__'+ibranch)
                f.__name__ = ibranch
                return f
                
            setattr(self, ibranch, function_factory(ibranch))
            
                
    def p4(self):
        p4 = ROOT.Math.LorentzVector('<ROOT::Math::PtEtaPhiM4D<double> >')(self.pt(), self.eta(), self.phi(), self.mass())
        return p4

    def __str__(self):
        toprint = 'particle    pt %.2f  eta %.2f  phi %.2f  charge %s  mass %.2f  pdgid %d' %(self.pt(), self.eta(), self.phi(), str(getattr(self, '__charge', 'NaN')), self.mass(), self.pdgId())
        return toprint

class RecoParticle(Particle):

    def __init__(self, tree, collection, ievent, idx, branches=base_branches+extra_branches, gen=None):
        super(RecoParticle, self).__init__(tree, collection, ievent, idx, branches)

        self.__gen_particle = gen
    
    def genp(self):
        return self.__gen_particle
    
    def set_genp(self, genp):
        self.__gen_particle = genp
    
    def __str__(self):
        toprint = super(RecoParticle, self).__str__()
        if self.genp():
            toprint += '\n'
            toprint += '\tgen match ' + self.genp().__str__()
        return toprint
    
class GenParticle(Particle):

    def __init__(self, tree, collection, ievent, idx, branches=base_branches+gen_branches):
        super(GenParticle, self).__init__(tree, collection, ievent, idx, branches)
        
        if self.genPartIdxMother()>=0:
            self.__mother = GenParticle(tree, collection, ievent, self.genPartIdxMother(), branches=branches)
        else:
            self.__mother = None
            
    def mother(self):
        return self.__mother

    def __str__(self):
        toprint = super(GenParticle, self).__str__()
        if self.mother():
            toprint += '\n'
            toprint += '\tmother   ' + self.mother().__str__()
        return toprint

class BCandidate(Particle):

    def __init__(self, tree, collection, ievent, idx, trg_mu, sel_mu, pi, branches=bcand_branches):
        super(BCandidate, self).__init__(tree, collection, ievent, idx, branches)
        
        self.__trg_mu = trg_mu
        self.__sel_mu = sel_mu
        self.__pi     = pi   
        self.__pdgId  = 0 

    def pdgId(self):
        return self.__pdgId
        
    def trg_mu(self):
        return self.__trg_mu
    
    def mu(self):
        return self.__sel_mu    

    def pi(self):
        return self.__pi    

    def checkFirstAncestor(self, pp, target_pdgid):
        if pp is None:
            return False
        if abs(pp.pdgId()) == target_pdgid :
            return True
        if hasattr(pp, 'mother') and pp.mother() is not None:
            mum = pp.mother()
            if self.checkFirstAncestor(mum, target_pdgid):
                return True
        return False

    def isMatched(self, b_pdgid=521, hnl_pdgid=9900015):
        if any([getattr(pp, 'genp', lambda : None)() is None for pp in [self.trg_mu(), self.mu(), self.pi()]]):
            return False
        
        
        fully_matched = self.checkFirstAncestor(self.trg_mu().genp(), b_pdgid  ) + \
                        self.checkFirstAncestor(self.mu().genp()    , b_pdgid  ) + \
                        self.checkFirstAncestor(self.pi().genp()    , b_pdgid  ) + \
                        self.checkFirstAncestor(self.mu().genp()    , hnl_pdgid) + \
                        self.checkFirstAncestor(self.pi().genp()    , hnl_pdgid)
        
        return fully_matched==5
        
    def __str__(self):
        toprint = super(BCandidate, self).__str__()
        toprint += '\n'
        toprint += '\thnl fitted mass %.2f  pt %.2f  eta %.2f  phi %.2f  charge %d' %(self.hnl_fit_mass(), self.hnl_fit_pt(), self.hnl_fit_eta(), self.hnl_fit_phi(), self.hnl_charge())
        toprint += '\n'
        toprint += '\t    cosine %.4f  vtx prob %.4f  Lxy %.2f  Lxy sig %.2f \tisMatched? %r' %(self.hnl_fit_cos2D(), self.sv_prob(), self.sv_lxy(), self.sv_lxy_sig(), self.isMatched())
        return toprint

##########################################################################################
##########################################################################################
           
if __name__ == '__main__': 

    fin = uproot4.open("BParkNANO_mc_10215.root")
    tree = fin["Events"]
    # tree.show()

    for iev in range(10):
        for icand in range(tree['nb'].array()[iev]):
            
            print '='*50
            print 'event %d, candidate %d' %(iev, icand)
                        
            pi_idx      = tree['b_pi_idx'].array()[iev][icand]
            pi          = RecoParticle(tree, 'ProbeTracks', iev, pi_idx)
            pi_genp_idx = tree['ProbeTracks_genPartIdx'].array()[iev][pi_idx]
            pi_genp     = GenParticle(tree, 'GenPart', iev, pi_genp_idx)
            pi.set_genp(pi_genp)

            trg_mu_idx      = tree['b_trg_mu_idx'].array()[iev][icand]
            trg_mu          = RecoParticle(tree, 'Muon', iev, trg_mu_idx)
            trg_mu_genp_idx = tree['Muon_genPartIdx'].array()[iev][trg_mu_idx]
            trg_mu_genp     = GenParticle(tree, 'GenPart', iev, trg_mu_genp_idx)
            trg_mu.set_genp(trg_mu_genp)

            sel_mu_idx      = tree['b_sel_mu_idx'].array()[iev][icand]
            sel_mu          = RecoParticle(tree, 'Muon', iev, sel_mu_idx)
            sel_mu_genp_idx = tree['Muon_genPartIdx'].array()[iev][sel_mu_idx]
            sel_mu_genp     = GenParticle(tree, 'GenPart', iev, sel_mu_genp_idx)
            sel_mu.set_genp(sel_mu_genp) 

            candidate = BCandidate(tree, 'b', iev, icand, trg_mu, sel_mu, pi)
                        
            print  'pion\n', pi, '\n'
            print  'trigger muon\n', trg_mu, '\n'
            print  'selected muon\n', sel_mu, '\n'
            print  'B candidate muon\n', candidate, '\n'

#             import ipdb ; ipdb.set_trace()        
#             raw_input("Press any key to continue")
