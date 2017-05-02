""" File: parray.py """
# lldb array printing script available at:
# http://stackoverflow.com/questions/7062173/view-array-in-lldb-equivalent-of-gdbs-operator-in-xcode-4-1
# usage:
# (lldb) command script import /path/to/parray.py
import lldb
import shlex

@lldb.command('parray', 'command script add -f parray.parray parray')
def parray(debugger, command, result, dict):

    target = debugger.GetSelectedTarget()
    process = target.GetProcess()
    thread = process.GetSelectedThread()
    frame = thread.GetSelectedFrame()

    args = shlex.split(command)
    if len(args) == 2:
        count = int(args[1])
        indices = range(count)
    elif len(args) == 3:
        first = int(args[1])
        count = int(args[2])
        indices = range(first, first + count)
    else:
        print 'Usage: parray ARRAY [FIRST] COUNT'
        return

    for i in indices:
        print frame.GetValueForVariablePath(args[0] + "[" + str(i) + "]")
