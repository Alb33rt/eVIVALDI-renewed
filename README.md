# eVIVALDI-renewed
Renewed version of the bacteriophage spatial interaction model for Python 3

## Translation Check
Under every file that is a Python 2.7 file, translations have been made.
To find where the current translation progress is, `Ctrl+f` and search for `CURRENT PROGRESS`

Send a push request and I will double check.

## Combine with original eVIVALDI
The project is a translated version of eVIVALDI, but under the Python 3 language, along with certain modules that are only available in this language.
You have to download the original eVIVALDI, and change the simulation program only. [Here is the link to eVIVALDI](https://gitlab.pasteur.fr/jsousa/eVIVALDI/-/tree/master/Src_eVIVALDI_9.2)

### Installing Packages
Make a new `virtualenv` named eVIVALDI, and install the following packages.
+ Pillow
+ Seaborn
+ MatplotLib

Then, navigate to the repository named `Src_eVIVALDI_9.2`, and __delete everything inside__

    git clone <thisrepo>
    > The folder will now be filled with new Python 3 Scripts.
 
That updates your eVIVALDI to be compatible with Python 3.
