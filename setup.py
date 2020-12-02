import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="geostar-pkg-simold", # Replace with your own username
    version="0.0.1",
    author="Simon Oldfield",
    author_email="simold@dtu.dk",
    description="A package of geospatial and related tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/soldfield/geostar",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering",
    ],
    python_requires='>=3.6',
)