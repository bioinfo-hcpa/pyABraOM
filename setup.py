from setuptools import setup
setup(
  name= "PyABraOM",
  packages = ['pyabraom'],
  version = '0.2.0',
  description= "A Python API to communicate with ABraOM database",
  author ="Paola Carneiro and Felipe Colombelli",
  author_email= "bioinfo@hcpa.edu.br",
  url="https://github.com/bioinfo-hcpa/pyABraOM",
  license='GNU General Public License v3.0',
  include_package_data=True,
  install_requires=["requests>=2.23.0",
                    "pandas>=1.1.5",
                    "times",
                    "tqdm"
                    ],
  py_modules=["Search_gene",
              "Search_Region",
              "Variant_ID", 
              'Searches'],
    classifiers=[
    'Development Status :: 3 - Alpha',      
    'Intended Audience :: Developers',      
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 3.8',
  ],
)