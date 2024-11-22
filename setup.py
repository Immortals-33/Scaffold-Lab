from setuptools import setup

setup(
    name="scaffold_lab",
    packages=[
        'data',
        'analysis',
        'openfold'
    ],
    package_dir={
        'data': './data',
        'analysis': './analysis',
        'openfold': './openfold',
    },
)
