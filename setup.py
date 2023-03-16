from setuptools import setup

setup(
    name="nccr_rnaseq",
    version="0.0.1",
    author="Anna Sintsova",
    author_email="ansintsova@ethz.ch",
    description=("RNAseq Pipeline"),
    license="LICENSE",
    keywords="nccrRna",
    url="https://github.com/MicrobiologyETHZ/NCCR_genomicsPipeline",
    install_requires=[
        'click', 'pyaml'],
    packages=['nccr_rnaseq'],
    entry_points={
        'console_scripts': ['nccrRna=nccr_rnaseq.main:main'],
    }
)

