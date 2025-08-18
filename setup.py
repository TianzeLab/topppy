import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="topppy",
    version="0.1.0",
    author="Cao Tianze",
    author_email="hnrcao@qq.com",
    description="Python implementation of ToppGene API wrapper",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tianzelab/topppy",
    project_urls={
        "Bug Tracker": "https://github.com/tianzelab/topppy/issues",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=['topppy'],
    python_requires=">=3.12",
    install_requires=['pandas>=2.2.0', 'requests>=2.32.0'],
    license="MIT",
    keywords=['scToppR', 'ToppGene', 'bioinformatics','pandas'],
    package_data={
        "": ["*.csv","*.txt"]
    }
)
