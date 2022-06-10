from setuptools import setup


setup(
    name="TFBMiner",
    version="1.2.0",
    description="A data acquisition and analysis pipeline for the rapid identification of putative transcription factor-based biosensors.",
    author="Tariq Joosab",
    project_urls={"Source": "https://github.com/tariqjoosab/tfb-miner"},
    python_requires="==3.10.1",
    platforms="OS independent",
    license="MIT License",
    zip_safe=False,
    entry_points={
        "console_scripts": [
            "TFBMiner = TFBMiner.__main__:main",
        ]
    },
    install_requires=[
        line.rstrip() for line in open("requirements.txt", "rt")
    ],
    packages=["TFBMiner"],
    package_dir={"TFBMiner": "TFBMiner"},
    include_package_data=True
    )