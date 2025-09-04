from setuptools import setup, find_packages
import os

# Read README if it exists
readme_path = os.path.join(os.path.dirname(__file__), "README.md")
with open(readme_path, "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name="guanaco-viz",
    version="0.1.0",
    author="Systems Immunometabolism Lab",
    author_email="your.email@example.com",
    description="GUANACO: Interactive visualization tool for single-cell data and genome browser",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Systems-Immunometabolism-Lab/GUANACO_updated",
    
    # Only include your guanaco package and its subpackages
    packages=find_packages(include=["guanaco", "guanaco.*"]),
    include_package_data=True,

    # Include non-code assets
    package_data={
        "guanaco": ["assets/*", "assets/**/*"],
        "guanaco.pages.matrix": ["cvd_color.json"],
    },

    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],

    python_requires=">=3.10",

    install_requires=[
        "dash==2.18.0",
        "dash-bootstrap-components==1.3.0",
        "dash-bio==1.0.2",
        "dash-draggable",
        "dash-ag-grid",
        "plotly>=5.19.0,<6.0.0",
        "pandas>=2.0.0",
        "scipy>=1.9.0",
        "statsmodels==0.14.4",
        "anndata",
        "muon==0.1.7",
        "logomaker==0.8",
        "pyfaidx==0.7.2",
        "pyjaspar==3.0.0",
        "aiobotocore==2.15.2",
        "boto3==1.35.36",
        "botocore==1.35.36",
        "s3transfer==0.10.3",
        "urllib3==1.26.14",
        "tqdm==4.64.1",
        "python-dotenv==1.0.0",
        "click",
        "gunicorn==21.2.0",
    ],

    entry_points={
        'console_scripts': [
            'guanaco=guanaco.cli:main',   # fix here!
        ],
    }

)
