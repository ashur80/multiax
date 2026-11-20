"""Loading of the transcribed test data (Appendix A) and parameters."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parent.parent / 'data'

LOADCASES = ['ax', 'tors', 'd0', 'd90', 'd45']
LABELS = {'ax': 'Axial', 'tors': 'Torsion', 'd0': 'komb., δ = 0°', 'd90': 'komb., δ = 90°',
          'd45': 'komb., δ = 45°', 'all': 'Alle Versuchsergebnisse'}
COLORS = {'ax': 'black', 'tors': 'tab:blue', 'd0': 'tab:green', 'd90': 'tab:red',
          'd45': 'tab:purple', 'all': 'tab:orange'}
MARKERS = {'ax': 'o', 'tors': 's', 'd0': 'D', 'd90': '^', 'd45': 'v'}


def load_params():
    with open(DATA / 'parameters.json', encoding='utf-8') as f:
        return json.load(f)


def load_tests():
    df = pd.read_csv(DATA / 'tests_appendix_a.csv', comment='#')
    df['F_a'] = df['F_a'].fillna(0.0)
    df['M_ta'] = df['M_ta'].fillna(0.0)
    df['delta'] = df['delta'].fillna(0.0)
    df['remark'] = df['remark'].fillna('')
    df['runout'] = df['remark'].isin(['DL', 'DL_crack'])
    df['raised'] = df['remark'].eq('hDL')
    return df


def include_in_statistics(df, runout_policy='exclude_DL'):
    """Mask of tests that enter the boxplots N_calc/N_exp.

    PRELIMINARY (author will confirm): 'exclude_DL' -> real runouts (DL) are only plotted,
    raised runouts (hDL) count with N_B at the raised level. Alternatives: 'exclude_all', 'include_all'.
    """
    if runout_policy == 'exclude_DL':
        return ~df['runout'].to_numpy()
    if runout_policy == 'exclude_all':
        return ~(df['runout'] | df['raised']).to_numpy()
    if runout_policy == 'include_all':
        return np.ones(len(df), bool)
    raise ValueError(runout_policy)
