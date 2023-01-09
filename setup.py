from setuptools import setup, find_packages

setup(name='flukit',
      version='0.0.2',
      description='the influenza surveillance toolkit',
      author='Ammar Aziz',
      author_email='ammar.aziz@mh.org.au',
      license='GPL3',
      package_data={'flukit': ['config/*', 'config/nextclade_datasets/**/*']},
      install_requires=['biopython==1.80', 'pandas', 'numpy', 'typer', 'rich', 'importlib_resources'],
      python_requires='>=3.9.15',
      include_package_data=True,
      entry_points={"console_scripts": ["flukit = flukit.flukit:app"]}
      )
