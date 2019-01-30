#/usr/bin/env python
# Coder: -.- Dr. Stefano Pirro (aka wynstep)

# This script performs a filtering stage according to the quality of variant calls

# importing fixed vars and methods
from vars import * # import file with vars
from functions import * # import file with functions

# import libraries
import os, sys, optparse

# Loading arguments necessary for the analysis
parser = optparse.OptionParser(description='Arguments for filtering the plasma ctDNA mutation file')
parser.add_option('-f','--vcf', dest="vcf", help='vcf file to parse')
parser.add_option('-o','--output', dest="output", help='filtered VCF file name')
args, remaining = parser.parse_args()

# loading VCF file into a
