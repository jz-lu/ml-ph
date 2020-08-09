# Workflow for a semiautomated analysis of electronic and phononic structures

**Remark**: everything is modular, so if you don't want a particular property, don't run that script, with the exception of relaxation, which needs to happen no matter what. If your input is already relaxed, you may of course skip the relaxation. You must run postProcess_relaxation.py (hence the prefaced underscore).

* Run a VASP relaxation on the original input files. Keep the default methods as outlined in Daniel's suggesions (these are in the README on the parent directory ph-practice).
* Run `_postProcess_relaxation.py` upon job completion to organize the files and make the proper directories. Specify in whatever order `'eledos'`, `'phdos'`, `'eleband'`, `'phband'` in the command line to specify which folders should be made.
* To get electronic DOS, run `get_eleDOS_analysis.py`. The output files will be placed in a "results" folder inside the DOS folder (created in the postProcess).
* To get a band structure, run `get_eleBand_analysis.py`. The output files will placed in a "results" folder inside the band folder.
* To get the relevant phonon properties, run `get_ph_analysis.py`. In the command line arguments, in whatever order, specify `'dos'` and/or `'band'` to tell the program what to calculate. As with electronic structures, the results will be placed in the appropriate subfolders.