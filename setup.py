from setuptools import setup

setup(
    name="workflow",
    version="0.0.1",
    author="Anna Sintsova",
    author_email="ansintsova@ethz.ch",
    description=("RNAseq Pipeline"),
    license="LICENSE",
    keywords="nccrRna",
    url="https://github.com/MicrobiologyETHZ/NCCR_genomicsPipeline",
    install_requires=[
        'click', 'pyaml'],
    packages=['workflow'],
    entry_points={
        'console_scripts': ['rnapipe=workflow.main:main'],
    }
)

