# coding: utf-8
import os
from setuptools import setup

def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'figaro', 'VERSION'))
    return versionFile.readline().strip()

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name="figaro",
    version=version(),
    description="FIGARO - An efficient and objective tool for optimizing microbiome rRNA gene trimming parameters",
    long_description=readme(),
    long_description_content_type='text/markdown',
    license="GPLv3",
    url="https://github.com/Zymo-Research/figaro",
    packages=["figaro", "figaro.defaults"],
    package_data={'figaro': ['VERSION']},
    scripts=["bin/figaro"],
    include_package_data=True,
    install_requires=list(req.strip() for req in open("requirements.txt")),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific Engineering :: Bio/Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.7"
    ],
    zip_safe=False
)
