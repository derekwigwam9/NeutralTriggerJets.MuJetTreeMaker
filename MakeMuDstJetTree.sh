#!/bin/bash
# 'MakeMuDstJetTree.sh'
# Derek Anderson
# 11.13.2017
#
# Use this to run 'MakeMuDstJetTree.C'
# in batch mode.

input="\"../../MuDstMatching/output/merged/pt4ff.matchWithMc.root\""
output="\"pp200r9pt4ff.et920vz55had.r03rm1chrg.root\""

root -b -q MakeMuDstJetTree.C\($input, $output, "true"\)
