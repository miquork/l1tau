# l1tau
Package for determining HB trigger inefficiencies, particularly DoubleIsoTau from 1-leg parameterization.

How to use:
1) root -l -b -q mk_CondFormats.C     [compile code]
2) root -l -b -q mk_L1Tau.C           [process .root files, then plot]

Can also re-run plotting part from #2 with
2a) root -l -b -q drawL1Tau.C\(\"IsoTau34\"\)
2b) root -l -b -q drawL1TauSelection.C\(\"23C\"\)

New efficiency maps go into 'rootfiles' subdirectory. Reference maps for Mjj prefiring are also accessed from there. Plots go in 'pdf' subdirectory.