# import subprocess
# import sys

# numPoscars = 4

# for i in range(1, numPoscars + 1):
#     pre = '' # Preface for POSCAR name, is it 00 or 0 or empty string?
#     if i <= 9:
#         pre = '00'
#     elif (10 <= i <= 99):
#         pre = '0'
#     elif 100 <= i <= 999:
#         pre = ''
#     else:
#         sys.exit("Too many POSCAR-displacement files to handle.") # This literally should never happen
    
#     subprocess.run(['mkdir', '%s/disp-%s'%('./', pre+str(i))])


# print('.' + ('/disp-%s'%(str(10))))

print("{%s}"%('hi'))