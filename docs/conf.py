import os

project = 'USTB'
copyright = '2025, University of Oslo'
author = 'USTB Contributors'

extensions = [
    'sphinx.ext.autodoc',
    'sphinxcontrib.matlab',
]

this_dir = os.path.dirname(os.path.abspath(__file__))
matlab_src_dir = os.path.abspath(os.path.join(this_dir, '..'))
primary_domain = 'mat'

exclude_patterns = ['_build']

html_theme = 'sphinx_rtd_theme'
html_static_path = []

matlab_short_links = True
matlab_keep_package_prefix = False

rst_prolog = """
.. |rarr| unicode:: U+2192 .. right arrow
"""
