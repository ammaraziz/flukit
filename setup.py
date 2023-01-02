from setuptools import setup, find_packages

setup(name='flukit',
      version='0.01',
      description='the influenza surveillance toolkit',
      author='Ammar Aziz',
      author_email='ammar.aziz@mh.org.au',
      license='GPL3',
      packages=['flukit', 'flukit.utils', 'flukit.config'],
      package_data={'flukit': ['config/*']},
      install_requires=['biopython', 'pandas', 'numpy', 'typer', 'rich', 'importlib_resources'],
      include_package_data=True,
      entry_points={"console_scripts": ["flukit = flukit.flukit:app"]}
      )
