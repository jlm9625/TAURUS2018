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


General Coding Guidelines for a Good Time:
- Name every file logically so that its easy to understand what its for. Also write comments at 
start of each script/notebook that describes what its for and maybe how to use it.
-Always add and commit at the end of working. And write a good commit message. Remember to push after commit.
-How to git:
Start with getting changes:
git pull 
git status (look for the new and modified stuff)
git add filename (modified and new files)
can also do 
git add *.py or git add *.ipynb which will add all python scripts/notebooks that are new/mod

Once everything is good to go do:
git commit -m 'A useful message on what you did, but not too long'
Send all to the cloud master repository
git push 
 
Login to github and check for pull requests.