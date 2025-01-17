'''
This file performs the automatic testing for all the cases in .........
'''

import numpy as np
import glob
import os 
import re
import pytest

def readstp(fn):
  f = open(fn, 'r')
  for i in f:
    exec(i)
  f.close()
  globals().update(locals())

def readin(fn,nin):
  f = open(fn, 'r')
  data = f.readlines()
  return data[nin]

def modin(fn,nin,strIn):
  f = open(fn, 'r')
  data = f.readlines()
  f.close()
  data[nin] = strIn
  with open(fn, 'w') as file:
    file.writelines( data )
  return

def chkRes():
  return
#############################################################
# define execution variables
wfold = os.getcwd()
mainFold = wfold[:[_.start() for _ in re.finditer('/',wfold)][-1]]
#compiler = 'generic-laptop'
compiler = 'generic-nvf'
#compiler = 'generic-gpu'
tfold = 'utils/testing/templateTest/'
#
doDBG = ['0']#,'1']
#doDBG = ['1']#,'1']
#doDBG = ['1']#,'1']
#doCuda= ['1']
doCuda= ['0']
doProc= ['2 2\n']#,'2 2\n']
#doProc= ['1 1\n']#,'2 2\n']
#doProc= ['1 2\n']#,'2 2\n']
#
cdir = os.getcwd()
apps = [ name for name in os.listdir() if os.path.isdir(os.path.join(name))]

os.system('rm '+wfold+'/test.out')
fileOut = open('test.out','w')
for dbgFlag in doDBG:
  for cudaFlag in doCuda:
    for app in apps:
        os.chdir(wfold+'/'+app)
        testcases = [ name for name in os.listdir() if os.path.isdir(os.path.join(name))]
        print(testcases)
        os.chdir(mainFold+'/src')
        #
        os.system('make clean APP='+app+'>/dev/null;'+
                  #'make APP='+app+' ARCH='+compiler+ ' Do_DBG='+dbgFlag+' -j8>/dev/null')
                  'make APP='+app+' ARCH='+compiler+ ' DO_DBG='+dbgFlag+' -j8>/dev/null')
        #
        if 'flutas' in os.listdir():
          print('app '+app+', debug='+dbgFlag+' : compiled-------> PASSED')
          fileOut.write('app '+app+', debug='+dbgFlag+' : compiled -------> PASSED\n')
          for case in testcases:
            #print(case)
            #print(testcases)
            os.chdir(wfold+'/'+app+'/'+case)
            if len(glob.glob('test*'))>0:
              os.system('cp '+mainFold+'/src/flutas .')
              for procN  in doProc:
                tmp = procN.split(' ')
                nmpi = int(tmp[0])*int(tmp[1])
                modin('dns.in',-2,procN)
                os.system('mpirun.openmpi -np '+str(nmpi)+' flutas>/dev/null')
                testout = pytest.main()
                testout = str(testout)[str(testout).find('.')+1:]
                fileOut.write('test '+case+' np='+str(nmpi)+'-------> '+testout+'\n') 
                #os.system('rm -rf data __pycache__ flutas;'+
                os.system('rm -rf data && rm -rf .pytest_cache && rm -rf __pycache__ && rm -rf flutas*;'+
                          'git checkout *.in')
        else:
          print('test '+app+', debug='+dbgFlag+' : compilation failed \n')
          exit()
fileOut.close()
