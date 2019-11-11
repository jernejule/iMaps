import os.path
import setuptools

# Get the long description from README.
with open('README.rst', 'r') as fh:
    long_description = fh.read()

# Get package metadata from '__about__.py' file.
about = {}
base_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(base_dir, 'imaps', '__about__.py'), 'r') as fh:
    exec(fh.read(), about)

setuptools.setup(
    name=about['__title__'],
    # use_scm_version=True,
    version=about['__version__'],
    description=about['__summary__'],
    long_description=long_description,
    long_description_content_type='text/x-rst',
    author=about['__author__'],
    author_email=about['__email__'],
    url=about['__url__'],
    license=about['__license__'],
    # Exclude tests from built/installed package.
    packages=setuptools.find_packages(
        exclude=['*.tests', '*.tests.*']
    ),
    python_requires='>=3.6, <3.8',
    install_requires=[
        'matplotlib',
        'numpy',
        'pandas==0.24.2',
        'pybedtools',
        'seaborn',
        'sklearn',
        'resdk>=11.0.1',
   ],
    extras_require={
        'docs': ['sphinx_rtd_theme'],
        'package': ['twine', 'wheel'],
        'test': [
            'check-manifest',
            'isort',
            'mock==1.3.0',
            'pycodestyle~=2.5.0',
            'pydocstyle~=3.0.0',
            'pylint~=2.3.1',
            'pytest-cov',
            'setuptools_scm',
            'twine',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='imaps CLIP iCLIP bioinformatics protein',
)
