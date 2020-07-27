# Problem: want to execute and manipulate linux commands and bash script directly from python.

# Solution 1: using os, which doesn't have output control....
print('START Method 1: os basic')
import os
os.system('echo "Hello there."')
os.system('echo "General Kenobi!"')
# os.system('ls -l') # Notice if we wanted to store any of these directory strings, os is not helpful
print('END Method 1\n')

# Solution 2: OS Pipeline opens a pipeline to the command line to pipe commands (sequential execution with outputs as next input).
# We can read output using .read() or .readlines() but can't run anything.
# This method is outdated by subprocesses since python 3.
# (We're still just using the os import)
print('START Method 2: os pipeline')
pipe = os.popen('echo No, I am your father.')
output = pipe.read()
output = output.strip() # Gets rid of the newline \n char at the end of the stream
print(output)
print('END Method 2\n')

# Solution 3: python subprocesses, most general method and the one we will use. Cannot pipe commands automatically, but we don't need it.
# Instead we can set stdout=subprocess.PIPE and then as input next do stdin = p1.stdout where p1 is name of subprocess object 1.
# There's a bit more to it to get it to work, see this link https://janakiev.com/blog/python-shell-commands/.
import subprocess
print('START Method 3: subprocess')
# subprocess.Popen is a class, so we call its constructor
# If we want to run a linux/shell command, we need to set "shell=True" in the constructor (defaults to False)
print('Example 1')
process1 = subprocess.Popen(['echo', 'More output'],
                     stdout=subprocess.PIPE, 
                     stderr=subprocess.PIPE)
stdout, stderr = process1.communicate() # Reads the process output
print(stdout, stderr)

print('Example 2')
process2 = subprocess.Popen(['ping', '-c 4', 'python.org'], 
                           stdout=subprocess.PIPE,
                           universal_newlines=True) # This converts default byte literal to string literal

# communicarte() waits until the end. We can ping more frequently by looping through.
# while True:
#     output = process2.stdout.readline()
#     print(output.strip())
#     # Do something else
#     return_code = process2.poll()
#     if return_code is not None:
#         print('RETURN CODE', return_code)
#         # Process has finished, read rest of the output 
#         for output in process2.stdout.readlines():
#             print(output.strip())
#         break
# If we run multple piped processes, we need to close all but the last ones
# More info on that: https://stackoverflow.com/questions/23074705/usage-of-stdout-close-in-pythons-subprocess-module-when-piping

# Finally, we can run a batch script using subprocess.run(). It's a simplification of subprocess.Popen.
print('Example 3')
run = subprocess.call('./basic_practice.sh')

print('END Method 3\n')