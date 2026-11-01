"""Regression tests of actual selection functions without a ROOT installation.

The analysis functions are compiled from their source AST, avoiding unrelated
plotting/ROOT imports. This tests their masks and numerical response accounting;
it is not a detector-data validation.
"""
import ast
from contextlib import contextmanager
import os
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'workflow'))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata


@contextmanager
def policy(value):
    previous = os.environ.get('YPSILON_PROTON_SELECTION')
    os.environ['YPSILON_PROTON_SELECTION'] = value
    try:
        yield
    finally:
        if previous is None:
            os.environ.pop('YPSILON_PROTON_SELECTION', None)
        else:
            os.environ['YPSILON_PROTON_SELECTION'] = previous


def functions(relative, names, extra=None):
    path = ROOT / relative
    tree = ast.parse(path.read_text(), filename=str(path))
    selected = [node for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in names]
    if {node.name for node in selected} != set(names):
        raise AssertionError('Missing tested functions in ' + str(path))
    module = ast.Module(body=selected, type_ignores=[])
    namespace = dict(np=np, pd=pd, Path=Path, MASSES={'1S':9.46030}, proton_tag_mask=proton_tag_mask,
                     proton_branches=proton_branches, proton_selection_metadata=proton_selection_metadata)
    namespace.update(extra or {})
    exec(compile(module, str(path), 'exec'), namespace)
    return namespace


def sample():
    n=6
    values={
        'weight':np.arange(1,n+1,dtype=float), 'weight_proj':np.arange(1,n+1,dtype=float),
        'true_mass':np.full(n,9.46030), 'true_Q2':np.array([.001]*5+[2.]),
        'true_W':np.array([22.,22.,22.,35.,35.,35.]),
        'true_y':np.zeros(n), 'true_pT':np.full(n,.4), 'true_t':np.full(n,.16),
        'true_ePlus_eta':np.zeros(n), 'true_decay_eMinus_eta':np.zeros(n),
        'reco_mass':np.full(n,9.46030), 'reco_W':np.array([22.,22.,22.,35.,35.,35.]),
        'reco_y':np.zeros(n), 'reco_pT':np.full(n,.4), 'sigma_epz':np.full(n,18.),
        'reco_topology':np.array([1,2,0,3,4,1]),
        'reco_tagged_proton_p':np.array([200.,180.,-999.,-999.,-999.,150.]),
    }
    return values


class FakeBranch:
    def __init__(self, values): self.values=values
    def array(self, library):
        if library!='np': raise AssertionError('unexpected library')
        return self.values.copy()


class FakeTree:
    def __init__(self, data): self.data=data
    def keys(self): return self.data.keys()
    def __getitem__(self, key): return FakeBranch(self.data[key])


class FakeSource:
    def __init__(self, data): self.tree=FakeTree(data)
    def __enter__(self): return self
    def __exit__(self, *args): pass
    def __contains__(self,key): return key=='UpsilonTree'
    def __getitem__(self,key):
        if key!='UpsilonTree': raise KeyError(key)
        return self.tree


class TestProtonAllstates(unittest.TestCase):
    def setUp(self):
        self.mode=policy('required');self.mode.__enter__()
        self.data=sample()
        self.args=SimpleNamespace(true_mass_window=.2,mass_window_sigma=3.,sigma_epz_center=18.,
            sigma_epz_halfwidth=1.8,w_min=20.,w_max=40.,w_edges=np.array([20.,30.,40.]),q2_cut=1.)
        self.props={'mass':9.46030,'sigma':.053}

    def tearDown(self): self.mode.__exit__(None,None,None)

    def test_performance_denominator_preserves_all_eligible_truth(self):
        f=functions('scripts/validation/plot_latest_sample_performance.py',['selection_masks'])['selection_masks']
        truth,reco=f(self.data,'1S',9.,1.,20.,40.)
        self.assertEqual(truth.sum(),5)
        self.assertEqual(reco.sum(),2)
        self.assertEqual(self.data['weight'][truth].sum(),15.)
        self.assertEqual(self.data['weight'][reco].sum(),3.)
        with policy('inclusive'):
            old_truth,old_reco=f(self.data,'1S',9.,1.,20.,40.)
        np.testing.assert_array_equal(truth,old_truth)
        self.assertEqual(old_reco.sum(),5)

    def test_w_extraction_and_response_select_same_tagged_category(self):
        extraction=functions('scripts/cross_sections/extract_gamma_p_cross_sections_bin_by_bin.py',['make_mask_for_state'])['make_mask_for_state']
        unfold=functions('scripts/unfolding/unfold_gamma_p_cross_section_vs_w.py',['make_selection_arrays'])['make_selection_arrays']
        truth=extraction(self.data,self.props,self.args,True,reco=False)
        reco=extraction(self.data,self.props,self.args,True,reco=True)
        _,response_truth,response_reco=unfold(pd.DataFrame(self.data),self.props,self.args,True)
        np.testing.assert_array_equal(truth,response_truth)
        np.testing.assert_array_equal(reco,response_reco)
        self.assertEqual(truth.sum(),5)
        self.assertEqual(reco.sum(),2)
        self.data['reco_topology'][:]=3;self.data['reco_tagged_proton_p'][:]=-999.
        self.assertEqual(extraction(self.data,self.props,self.args,True,reco=False).sum(),5)
        self.assertEqual(extraction(self.data,self.props,self.args,True,reco=True).sum(),0)

    def test_t_response_keeps_untagged_events_in_efficiency_denominator(self):
        ns=functions('scripts/unfolding/unfold_t_spectrum.py',['hist_weighted','response_matrix'])
        x=np.array([.1,.1,.3]);weights=np.array([2.,3.,5.])
        tagged=proton_tag_mask({'reco_topology':np.array([1,3,2]),'reco_tagged_proton_p':np.array([200.,-999.,180.])})
        H,_,truth,_,response,eff=ns['response_matrix'](x,x,weights,np.array([0.,.2,.4]),tagged)
        np.testing.assert_array_equal(truth,[5.,5.])
        np.testing.assert_array_equal(H,[[2.,0.],[0.,5.]])
        np.testing.assert_allclose(eff,[.4,1.])
        np.testing.assert_allclose(response.sum(axis=1),eff)

    def test_mass_and_kinematics_loaders_preserve_rows_and_source_weights(self):
        import_source=SimpleNamespace(open=lambda _:FakeSource(self.data))
        mass_path='scripts/reporting/render_note_mass_efficiency.py'
        module=ast.parse((ROOT/mass_path).read_text())
        required=next(ast.literal_eval(node.value) for node in module.body if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='REQUIRED' for t in node.targets))
        loader=functions(mass_path,['load_sample'],{'uproot':import_source,'REQUIRED':required})['load_sample']
        loaded=loader('sample.root','UpsilonTree','9x275',9.,20.,40.,1.)
        self.assertEqual(len(loaded['data']['weight']),6)
        self.assertEqual(loaded['truth'].sum(),5)
        self.assertEqual(loaded['accepted'].sum(),5)
        self.assertEqual(loaded['reco'].sum(),2)
        self.assertEqual(loaded['data']['weight'].sum(),21.)
        self.assertEqual(loaded['audit']['proton_selection']['mode'],'required')
        kloader=functions('scripts/reporting/render_note_kinematics.py',['load_sample'],{'uproot':import_source})['load_sample']
        ka=SimpleNamespace(tree='UpsilonTree',state='1S',q2_cut=1.,stat_mode='mc')
        loaded=kloader('sample.root','9x275',9.,20.,40.,ka)
        self.assertEqual(loaded['truth'].sum(),5)
        self.assertEqual(loaded['projected'].sum(),2)
        self.assertEqual(loaded['selected'].sum(),2)
        self.assertEqual(loaded['data']['weight'].sum(),21.)

    def test_missing_and_inconsistent_tags_fail_closed_but_inclusive_is_explicit(self):
        f=functions('scripts/validation/plot_latest_sample_performance.py',['selection_masks'])['selection_masks']
        missing={k:v for k,v in self.data.items() if k not in proton_branches()}
        with self.assertRaises(ValueError): f(missing,'1S',9.,1.,20.,40.)
        with policy('inclusive'):
            truth,reco=f(missing,'1S',9.,1.,20.,40.)
            self.assertEqual(truth.sum(),5);self.assertEqual(reco.sum(),5)
        self.data['reco_tagged_proton_p'][3]=100.
        with self.assertRaises(ValueError): f(self.data,'1S',9.,1.,20.,40.)


if __name__=='__main__': unittest.main()
