#!/usr/bin/env python

from setuptools import setup

setup(name='SurveySim',
      version='1.0',
      install_requires=['numpy>=1.0','scipy','matplotlib>=1.4','seaborn','astropy','pyfits'],
      py_modules=['SurveySimObsFile','SurveySimGUI'],
      packages=['SurveySim'],
      author='Noah Kurinsky',
      author_email='kurinsky@stanford.edu',
      url='http://cosmos2.phy.tufts.edu/~asajina/SurveySim.html',
      scripts=['SurveySimFilters','SurveySimGUI','SurveySimObsFile','SurveySimGUI.py','SurveySimObsFile.py'])
