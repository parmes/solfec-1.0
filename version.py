# convert git shortehand hash into version string
import subprocess

git = ["git",  "log",  "--pretty=format:'%h'", "-n", "1"]
p = subprocess.Popen(git, stdout=subprocess.PIPE)
commit = p.communicate()[0]
git = ["git", "show", "-s", "--format=%cd", "--date=short", commit.replace("'",'')]
p = subprocess.Popen(git, stdout=subprocess.PIPE)
date = p.communicate()[0].replace('\n','')
with open('version.h','w') as f:
    f.write('#define VERSION_HASH %s\n' % commit.replace("'",'"'))
    f.write('#define VERSION_DATE "%s"' % date)
