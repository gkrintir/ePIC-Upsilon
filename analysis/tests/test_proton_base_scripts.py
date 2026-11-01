"""Synthetic selection regressions; no ROOT files or PyROOT required."""
import ast
import math
import os
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'workflow'))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata


def selected_functions(path, names, extra=None):
    """Execute original function bodies without importing plotting/I/O dependencies."""
    tree = ast.parse((ROOT / path).read_text())
    body = [node for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in names]
    ns = dict(np=np, os=os, proton_tag_mask=proton_tag_mask,
              proton_branches=proton_branches, proton_selection_metadata=proton_selection_metadata)
    if extra:
        ns.update(extra)
    exec(compile(ast.Module(body=body, type_ignores=[]), str(ROOT/path), 'exec'), ns)
    return ns


def sample():
    n=6
    return dict(weight=np.ones(n), weight_visible=np.ones(n), weight_proj=np.ones(n),
                weight_proj_visible=np.ones(n), true_mass=np.full(n,9.4603),
                true_Q2=np.full(n,.001), true_W=np.full(n,50.), reco_W=np.full(n,50.),
                reco_mass=np.full(n,9.4603), sigma_epz=np.full(n,18.),
                true_t=np.full(n,.3), reco_t=np.full(n,.3), reco_t_eX=np.full(n,.3),
                reco_t_BABE=np.zeros(n), reco_t_eXBABE=np.full(n,.3), reco_t_proj=np.full(n,.3),
                true_pT2=np.full(n,.3), reco_pT2=np.full(n,.3), reco_pT=np.full(n,np.sqrt(.3)),
                reco_topology=np.array([1,2,3,4,0,3]),
                reco_tagged_proton_p=np.array([260.,250.,-999.,-999.,-999.,-999.]))


class SelectedMasks(Exception):
    def __init__(self, truth, reco, weight):
        self.truth, self.reco, self.weight = truth, reco, weight


class ProtonBaseScripts(unittest.TestCase):
    def setUp(self):
        self.env = patch.dict(os.environ, {'YPSILON_PROTON_SELECTION':'required'})
        self.env.start()
        self.addCleanup(self.env.stop)

    def test_report_histogram_truth_unchanged_and_tag_not_t_proxy(self):
        ns=selected_functions('scripts/validation/reconstruct_report_histograms_from_tuple.py', ['selection_mask'])
        a=sample()
        a['sigma_epz'][:]=36.
        self.assertEqual(ns['selection_mask'](a,'truth').sum(),6)
        self.assertEqual(ns['selection_mask'](a,'all').sum(),6)
        self.assertEqual(ns['selection_mask'](a,'reco').sum(),2)
        # Tagged |t|=0 remains tagged; no arbitrary t>0 proxy.
        self.assertEqual(ns['selection_mask'](a,'tagged').sum(),2)
        with patch.dict(os.environ, {'YPSILON_PROTON_SELECTION':'inclusive'}):
            self.assertEqual(ns['selection_mask'](a,'reco').sum(),6)
            self.assertEqual(ns['selection_mask'](a,'tagged').sum(),0)

    def test_benchmark_denominator_and_response_efficiency(self):
        names=['make_reco_t','hist_weighted','response_matrix','purity_stability','iterative_bayes','weighted_quantile','analyze_method']
        extra=dict(UPSILON={'1S':dict(mass=9.4603,sigma=.053)},
                   METHODS={'eX':dict(branch='reco_t_eX',label='eX')},
                   branching_fractions_from_args=lambda _: {'1S':.0234},
                   branching_fraction_metadata=lambda _: {'branching_fraction_convention':'test'})
        ns=selected_functions('scripts/validation/benchmark_t_estimators.py',names,extra)
        args=types.SimpleNamespace(state='1S', t_edges=np.array([0.,1.]),sigma_epz_halfwidth=None,
             sigma_epz_halfwidth_fraction=.2,reco_mass_min=8.,reco_mass_max=12.,seed=123,
             train_fraction=.5,n_iterations=4)
        a=sample()
        result=ns['analyze_method'](a,np.ones(6,dtype=bool),'test',9.,'eX',args,20.,100.)[0]
        self.assertEqual(result['sum_truth_visible'],6.)
        self.assertEqual(result['sum_reco_visible'],2.)
        self.assertAlmostEqual(result['average_efficiency'],1/3)
        with patch.dict(os.environ, {'YPSILON_PROTON_SELECTION':'inclusive'}):
            result=ns['analyze_method'](a,np.ones(6,dtype=bool),'test',9.,'eX',args,20.,100.)[0]
        self.assertEqual(result['sum_truth_visible'],6.)
        self.assertEqual(result['sum_reco_visible'],6.)

    def differential(self, slope=False):
        def capture(true,reco,weight,truth_mask,reco_mask,bins,args):
            raise SelectedMasks(truth_mask,reco_mask,weight)
        extra=dict(STATES={'1S':dict(mass=9.4603,sigma=.053)},
             branching_fractions_from_args=lambda _: {'1S':.0234},
             get_visible_xsec_for_state=lambda *_:1., binned_components=capture,
             VARIABLES={'pT2':dict(true_expr='true_pT2',reco_expr='reco_pT2')})
        ns=selected_functions('scripts/cross_sections/extract_differential_gamma_p_bin_by_bin.py',
                              ['extract_for_state','extract_b_slope_for_state','parse_edges'],extra)
        args=types.SimpleNamespace(state_normalize=False,truth_mass_window=.2,tuple_lumi_invfb=10.,
              sigma_epz_center=None,sigma_epz_halfwidth=None,electron_energy=9.,mass_min=None,
              mass_max=None,reco_mass_window_sigma=3.,w_min=20.,w_max=100.,q2_cut=1.,
              output_convention='inclusive_upsilon',b_slope_zeus_mass_window=False,
              b_slope_variable='pT2',b_slope_bins='0,1')
        a=sample()
        with self.assertRaises(SelectedMasks) as caught:
            if slope:
                ns['extract_b_slope_for_state'](a,args,'1S',.1,10000.)
            else:
                ns['extract_for_state'](a,args,'1S',{'pT2':np.array([0.,1.])},.1,10000.)
        check=caught.exception
        self.assertEqual(check.truth.sum(),6)
        self.assertEqual(check.reco.sum(),2)
        self.assertEqual(check.weight.sum(),6.)

    def test_differential_truth_and_weight_unchanged(self):
        self.differential()

    def test_slope_truth_and_weight_unchanged(self):
        self.differential(slope=True)

    def test_w_comparison_preserves_truth_rate_and_reduces_reco_support(self):
        ns=selected_functions('scripts/validation/compare_w_samples.py',
             ['compute_curve','histogram','state_mask','compare_curves'],
             dict(math=math,integrated_flux=lambda *_:.1,
                  sample_branching_fraction=lambda *_:.0234))
        cfg=dict(mass_gev=9.4603,true_mass_halfwidth_gev=.2,reco_mass_halfwidth_gev=.159,
                 lumi_invfb=10.,min_reconstruction_ratio=0.)
        spec=dict(electron_energy=9.,proton_energy=275.,key='test',
                  normalization_mode='all_generated',visible_xsec_fb=1.)
        a=sample()
        item=dict(data=a,weights=a['weight'],settings=cfg,spec=spec)
        tagged=ns['compute_curve'](item,[20.,100.],1.,1.)
        with patch.dict(os.environ, {'YPSILON_PROTON_SELECTION':'inclusive'}):
            inclusive=ns['compute_curve'](item,[20.,100.],1.,1.)
        self.assertEqual(tagged['rows'][0]['true_rows'],6)
        self.assertEqual(tagged['rows'][0]['reco_rows'],2)
        self.assertEqual(inclusive['rows'][0]['reco_rows'],6)
        self.assertEqual(tagged['rows'][0]['sigma_gamma_p_truth_pb'],
                         inclusive['rows'][0]['sigma_gamma_p_truth_pb'])
        with self.assertRaisesRegex(ValueError,'different proton selections'):
            ns['compare_curves'](tagged,inclusive)

    def test_benchmark_loader_validates_before_nan_drop(self):
        a=sample()
        a['reco_topology'][5]=1
        a['reco_tagged_proton_p'][5]=np.nan
        a['reco_mass'][5]=np.nan
        class Branch:
            def __init__(self,name): self.name=name
            def array(self,library): return a[self.name]
        class Tree:
            def __getitem__(self,name): return Branch(name)
        class Source:
            def __enter__(self): return self
            def __exit__(self,*_): pass
            def __getitem__(self,name): return Tree()
        extra=dict(UPSILON={'1S':dict(mass=9.4603,sigma=.053)},
                   branching_fractions_from_args=lambda _: {'1S':.0234},
                   branching_fraction_metadata=lambda _: {},
                   uproot=types.SimpleNamespace(open=lambda _:Source()))
        ns=selected_functions('scripts/validation/benchmark_t_estimators.py',['load_sample'],extra)
        args=types.SimpleNamespace(state='1S',tree='UpsilonTree',truth_mass_window=.2,
                                  q2_cut=1.,state_norm_mode='raw',lumi_invfb=10.)
        with self.assertRaisesRegex(ValueError,'schema mismatch'):
            ns['load_sample']('dummy.root','test',9.,1.,args)

if __name__=='__main__':
    unittest.main()
