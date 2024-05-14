import setuptools

with open("README.md", "r", encoding="utf-8") as fhand:
    long_description = fhand.read()

setuptools.setup(
    name="TLseq",
    version="0.0.1",
    author="OU",
    author_email="ou@gmail.com",
    description=("A demo packafe for addition"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ouslalu/TL_SEQ_APP",
    project_urls={
        "Bug Tracker": "https://github.com/ouslalu/TL_SEQ_APP/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "pandas",
        "numpy",
        "matplotlib",
        "requests",
        "Biopython",
        "glob",
        "Seaborn",
        "pysam",
        "itertools",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "TLseq = fragment_preprocessing.cli:main",
        ]
    },
)
