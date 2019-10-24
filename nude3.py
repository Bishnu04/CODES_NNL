#!/apps/python/3.4.3/bin/python3

# Toshiyuki Gogami
# Nov 2, 2018
# this code helps to run the nude3.cc code which is ued for the data rreduction purpose. To run this code just do pytho3 nude3.py.
# This code need a list of root files that you like to run by nude3.cc code



import sys
import time, os.path
from subprocess import call
#import concurrent.futures
from logging import StreamHandler, Formatter, INFO, getLogger
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures.process import ProcessPoolExecutor 
import numpy as np

trigflag = 5 # (coin)
#trigflag = 4 # (RHRS)
#trigflag = 1 # (LHRS)
nworkers=4 # no of computer for work

#runfile = "second_part_h2.dat" # 2nd half of H2 runs
#runfile = "runlist_h2_1.dat" # file that contains input root file information
#runfile = "H22_T.dat" # for test purpose
#runfile = "run_list_h2_1.dat" # first half of H2 dataemacs 
#runfile = "T_1st.dat" # for test purpose
#runfile = "T_2nd.dat" # T data run # 111369 to 111832 that is 2nd group of T data
runfile = "T_2nd_skipped.dat"
#runfile = "test_list.dat" # for test purpose

thisfile = "nude3.py"

def nude_start(command):
    time.sleep(1.0)
    call(command,shell=True)

def main():
    comlist = []
    #inputfile = open("h2_2.dat","r")
    inputfile = open(runfile,"r")
    lines = inputfile.readlines()
    for line in lines:
        data = line.split()
        #com = "./nude " + data[0]+ " " +data[1]
        com = "root -l -q \"nude3.cc(" + data[0]+ "," +data[1] + "," + str(trigflag)
        com2 = com + ")\";"
        #call(com,shell=True)
        comlist.append(com2)
    with ProcessPoolExecutor(max_workers=nworkers) as executor:
        executor.map(nude_start,comlist)


stime = time.time()
main()
print("\n Jobs were done in %.0f sec \n" % float(time.time()-stime))
