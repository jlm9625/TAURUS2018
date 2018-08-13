# TAURUS2018
Code and python notebooks for Jordyn Mascarenas-Wells' Taurus summer project 
Main goal is to produce a model that can predict ages/masses/properties of stars using 
Gaia parallaxes and various photometry.

Current status:
- Model includes binaries, age, mass, and uses Gaia G/BP/RP photometry. No [Fe/H] or extinction yet.
- Model spans 0.1-1 Msun and 1-3000 Myr.
- Test run on Taurus window without metallicity or IR mags was mostly successful in bounds of model.
- Most of the code is in python notebooks, model file generation code is in a separate script.

TODO:
- ACR+JMW: Think about planning regular skype meetings.
- ACR: Make a reliable code for cross matching Gaia catalog entries to 2MASS to get IR photometry.
- ACR: Generalize the model generation code and document it so that it's useful for other people
in the future.
- JMW: Use dummy 2MASS numbers with giant errors to incorporate the metallicity axis into the 
models. 
- JMW: Once metallicity is done, think about extinction. 
- JMW+ACR?: Start converting the jupyter note code that is basically finished into stand-alone python 
scripts that can be easily imported and used on e.g. the TACC super computers.
