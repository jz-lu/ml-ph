import subprocess

test = subprocess.run(['mv', 'test', 'line_generator'], capture_output=True, universal_newlines=True)
print(test.stderr == '')