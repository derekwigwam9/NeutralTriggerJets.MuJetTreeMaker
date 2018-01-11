#!/bin/bash
# 'MakeMuDstJetTree.sh'
# Derek Anderson
# 11.13.2017
#
# Use this to run 'MakeMuDstJetTree.C'
# in batch mode.

root -b -q MakeMuDstJetTree.C\("true"\)
