import os
import setuptools

with open("README.rst", "r", encoding="utf-8") as f:
    long_description = f.read()


with open(os.path.join("atmospheric_cherenkov_response", "version.py")) as f:
    txt = f.read()
    last_line = txt.splitlines()[-1]
    version_string = last_line.split()[-1]
    version = version_string.strip("\"'")


setuptools.setup(
    name="atmospheric_cherenkov_response_cherenkov-plenoscope-project",
    version=version,
    description="Estimate response-functions for the atmospheric Cherenkov-method",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/cherenkov-plenoscope/atmospheric_cherenkov_response",
    author="Sebastian Achim Mueller, Werner Hofmann",
    author_email="sebastian-achim.mueller@mpi-hd.mpg.de",
    packages=[
        "atmospheric_cherenkov_response",
        "atmospheric_cherenkov_response.instrument",
        "atmospheric_cherenkov_response.instrument.toy",
        "atmospheric_cherenkov_response.grid",
        "atmospheric_cherenkov_response.production",
        "atmospheric_cherenkov_response.analysis",
    ],
    package_data={"atmospheric_cherenkov_response": []},
    install_requires=[
        "homogeneous_transformation",
        "json_utils_sebastian-achim-mueller",
        "rename_after_writing",
        "photon_spectra_cherenkov-plenoscope-project",
        "thin_lens",
        "optic_object_wavefronts",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
    ],
)
