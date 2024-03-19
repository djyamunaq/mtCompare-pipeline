from setuptools import setup
setup(
    name = 'mtCompare-pipeline',
    version = '0.1.0',
    packages = ['sgetpipe'],
    entry_points = {
        'console_scripts': [
            'mtcomppipe = mtcomppipe.__main__:main'
        ]
    })