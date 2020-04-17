from setuptools import setup, find_packages

setup(
    name="phd3",
    version="1.0.0",
    author=["Matthew R. Hennefarth", "David J. Reilley"],
    package_data={'': ['*.json', '*.j2', '*.config']},
    packages=find_packages(),
    install_requires=[ 'jinja2', 'sklearn', 'propka', 'numpy'], 
    entry_points = {
        'console_scripts' : [
            'setupphd.py=phd3.bin.setupphd:main',
            'runphd.py=phd3.bin.runphd:main',
            'submitphd.py=phd3.bin.submitphd:main',
            'submitdmd.py=phd3.bin.submitdmd:main',
            'setupturbomole.py=phd3.bin.setupturbomole:main',
            'runturbomole.py=phd3.bin.runturbomole:main',
            'submitturbomole.py=phd3.bin.submitturbomole:main',
            'tfe.py=phd3.bin.tfe:main',
            'relabelpdb.py=phd3.bin.relabelpdb:main',
            'setupdmd.py=phd3.bin.setupdmd:main',
            'rundmd.py=phd3.bin.rundmd:main',
            'm2p=phd3.bin.movietopdb:main'
            ]
    },
)
