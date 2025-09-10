import os, sys
sys.path.append(os.path.abspath("../"))

language = "en"

extensions = [
    "sphinx.ext.autodoc",    
    "sphinx.ext.napoleon", 
    "sphinx.ext.viewcode",  
]
# myst_enable_extensions = ["colon_fence"]
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '.history']
html_baseurl = "/guanaco-viz/"

master_doc = "index"

source_suffix = ".rst"

html_theme = "sphinx_rtd_theme"

html_static_path = ["_static"]
html_theme_options = {
    "collapse_navigation": False,
    "sticky_navigation": True,
}

project = 'GUANACO'
author = 'Xiyuan Zhang, Marian Kuddus'
