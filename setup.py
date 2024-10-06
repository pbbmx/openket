from setuptools import setup, find_packages

setup(
    name="openket",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[],  # Aquí puedes agregar dependencias si tienes
    author="ultrxvioletx",
    description="Free software for manipulating quantum objects in Dirac notation",
    url="https://github.com/pbbmx/openket.git",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)