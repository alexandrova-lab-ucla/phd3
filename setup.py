from setuptools import setup, find_packages

setup(
    name="phd3",
    version="1.0.0",
    author="Matthew R. Hennefarth",
    packages=find_packages(),
    install_requires=['turbopy', 'jinja2', 'dmdpy']
)