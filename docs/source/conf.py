# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'BSMArt'
copyright = '2026, Mark Goodsell'
author = 'Mark Goodsell'
release = '2.0'

import sys
import os
sys.path.insert(0, os.path.abspath('../../src'))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [  'myst_parser',
                'sphinx.ext.autodoc',
                'sphinx.ext.napoleon',  # Optional: Allows Google/NumPy style docstrings
                'sphinx.ext.viewcode',
                'sphinx.ext.autosummary']

autosummary_generate = True
autodoc_mock_imports = ["torch","cmaes","deap","evosax","imblearn","jax","skimage","pyod","yaml","Higgs","anyBSM","flavio","wilson","pandas","sklearn","mpi4py","psutil","tqdm","rich","numpy","scipy","matplotlib","emcee","wget","requests","seaborn","corner","six","vegas"]

templates_path = ['_templates']
exclude_patterns = []

html_show_sourcelink = False  # dont show link to the documentation source

source_suffix = {
    '.md': 'markdown',
    '.rst': 'restructuredtext',
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'classic'
html_theme = 'sphinx_rtd_theme'

#html_theme_options = {
#    "light_logo": "logo-light.png",
#    "dark_logo": "logo-dark.png",
#}

def process_docstring(app, what, name, obj, options, lines):
    """
    Hook to append formatted __meta__ information to module docstrings.
    """
    if what == "module" and hasattr(obj, "__meta__"):
        meta = obj.__meta__
        
        # Add a blank line separator
        lines.append("")
        lines.append("Information")
        lines.append("-----------")
        lines.append("")

        if "name" in meta:
            lines.append(f"**BSMArt Name:** {meta['name']}")
            lines.append("")

        if "requires" in meta:
            lines.append("**Requires:**")
            for req in meta["requires"]:
                lines.append(f"   * {req}")
            lines.append("")

        if "settings" in meta:
            lines.append("**Settings:**")
            lines.append("")
            # Using a definition list for settings
            for key, desc in meta["settings"].items():
                #lines.append(f"    * {key}")
                
                if isinstance(desc, dict):
                    lines.append(f"{key}")
                    for k,v in desc.items():
                        lines.append(f"   * **{k}**: {v}")
                else:
                    lines.append(f"   * **{key}**: {desc}")


            lines.append("")

def setup(app):
    app.connect("autodoc-process-docstring", process_docstring)
