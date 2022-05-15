import subprocess

cdout = subprocess.run('cd ~/codes/ml-ph/testfiles/Gr6', shell=True, capture_output=True, universal_newlines=True)
print(cdout.stderr)

lsout = subprocess.run(['ls'], capture_output=True, universal_newlines=True)

print(lsout.stdout)
