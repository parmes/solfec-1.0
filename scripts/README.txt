This directory contains Python scripts, modules and packages to extend Solfec.

Some of these are intended to be imported by Solfec input scripts. To tell python how to find these
scripts you have a couple of options:

1. Add the absolute path to the ../scripts directory (including '/scripts' iself) to an environment variable PYTHONPATH.

2. a. Create a text file with a name ending '.pth'.
   b. Add the above absolute path to this file
   c. Save the file to one of Solfec's library directories, e.g:
      
      Windows: C:\Python26\Lib\site-packages
      Linux: /usr/lib64/python2.6/site-packages

On Linux, if you don't have write access to the above "global" site library directory you can add it to a user-library instead.  To find where this is run:
  python -c 'import site;print site.USER_SITE'
Note the path it returns may not exist yet.