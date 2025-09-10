import os, sys
sys.path.append(os.path.abspath("../"))  # 可选

language = "en"

# 扩展
extensions = [
    # "myst_parser",            # 支持 Markdown (commented out since we're using RST)
    "sphinx.ext.autodoc",     # 可选：从 Python docstring 生成文档
    "sphinx.ext.napoleon",    # 可选：支持 Google/NumPy 风格注释
    "sphinx.ext.viewcode",    # 可选：显示源码链接
]
# myst_enable_extensions = ["colon_fence"]
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '.history']
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
