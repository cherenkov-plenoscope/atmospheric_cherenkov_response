import setuptools

with open("README.rst", "r", encoding="utf-8") as f:
    long_description = f.read()

setuptools.setup(
    name="atmospheric_cherenkov_response_sebastian-achim-mueller",
    version="0.0.2",
    description="Estimate response-functions for the atmospheric Cherenkov-method",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/cherenkov-plenoscope/atmospheric_cherenkov_response.git",
    author="Sebastian Achim Mueller, Werner Hofmann",
    author_email="sebastian-achim.mueller@mpi-hd.mpg.de",
    packages=["atmospheric_cherenkov_response"],
    install_requires=["homogeneous_transformation_sebastian-achim-mueller"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
    ],
)
