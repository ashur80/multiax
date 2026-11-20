"""Build the .ipynb notebooks from the percent-format sources in notebooks/src/.

Usage (from the repository root):
    python notebooks/build_notebooks.py            # build only
    python notebooks/build_notebooks.py --execute  # build and execute (outputs stored in the .ipynb)

Cell markers: '# %%' starts a code cell, '# %% [markdown]' a markdown cell whose lines start with '# '.
"""
import sys
from pathlib import Path
import nbformat
from nbformat.v4 import new_notebook, new_code_cell, new_markdown_cell

HERE = Path(__file__).resolve().parent


def parse(src: str):
    cells, kind, buf = [], None, []

    def flush():
        if kind is None:
            return
        text = '\n'.join(buf).strip('\n')
        if kind == 'md':
            text = '\n'.join(l[2:] if l.startswith('# ') else l.lstrip('#') for l in text.splitlines())
            cells.append(new_markdown_cell(text))
        elif text.strip():
            cells.append(new_code_cell(text))

    for line in src.splitlines():
        if line.startswith('# %%'):
            flush()
            kind, buf = ('md' if '[markdown]' in line else 'code'), []
        else:
            buf.append(line)
    flush()
    return cells


def build(execute=False):
    for src in sorted((HERE / 'src').glob('*.py')):
        nb = new_notebook(cells=parse(src.read_text(encoding='utf-8')))
        nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
        out = HERE / (src.stem + '.ipynb')
        if execute:
            from nbconvert.preprocessors import ExecutePreprocessor
            ExecutePreprocessor(timeout=3600, kernel_name='python3').preprocess(nb, {'metadata': {'path': str(HERE)}})
        nbformat.write(nb, out)
        print('wrote', out.name)


if __name__ == '__main__':
    build(execute='--execute' in sys.argv)
