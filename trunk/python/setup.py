#!/usr/bin/env python

from distutils.core import setup

setup(name='SurveySim',
      version='1.0',
      py_modules=['ModelFile','ObsFile','OutputFile','SurveySim','filters'],
      author='Noah Kurinsky',
      author_email='kurinsky@stanford.edu',
      url='http://cosmos2.phy.tufts.edu/~asajina/SurveySim.html',
      scripts=['SurveySimFilters','SurveySimGUI','SurveySimObsFile'])
