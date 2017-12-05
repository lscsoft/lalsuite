# Copyright (C) 2013 Lindy Blackburn, Reed Essick and Ruslan Vaulin
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation,.

## \defgroup laldetchar_py_idq iDQ pipeline
## \ingroup laldetchar_python
"""
 Python modules for iDQ pipeline. The pipeline is designed for low-latency detection of transient artifacts in GW data. 

event: module for constructing and manipulating glitch events using data (e.g. triggers) from auxiliary and GW channels
idq_gdb_utils: module with utility functions for generating iDQ input to GraceDB.
idq:  main classes and function used to construct the workflow for iDQ pipeline
idq_tables: definitions of tables that store various data characterizing glitch events
idq_summary_plots: functions for making variety of plots for summarizing and interpreting  iDQ pipeline output.
ovl: classes and functions for Ordered-Veto-List (OVL) algorithm used by iDQ pipeline.
pdf_estimation: functions for estimation of probability density function(s) of glitch-rank and computation of corresponding likelihood ratio.
svmkit: utility functions for Support Vector Machine (SVM) - one of the algorithms used by the iDQ pipeline.  
auxmvc: classes for defining jobs and nodes for machine learning algorithms (MLAs)
auxmvc_utils: functions for handling auxmvc feature vectors used in MLAs 
"""

# \author Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>) and Ruslan Vaulin (<ruslan.vaulin@ligo.org)
# \heading{Synopsis}
# ~~~
# from laldetchar.idq import idq
# ~~~
# \heading{Example}
# \code
# from laldetchar.idq import idq
# \endcode

from laldetchar import git_version as version
__author__ = "Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>) and Ruslan Vaulin (<ruslan.vaulin@ligo.org)"
__version__ = version.id
__date__ = version.date
