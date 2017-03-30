# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 11:53:49 2017

@author: johan
"""
import subprocess, os, glob#, signal

def make():
    subprocess.check_call('make',cwd='..')

def first_run(folder): #f=folder
    files = glob.glob(folder+'logfile*')
    for file in files:
        os.remove(file)
    files = glob.glob(folder+'*.dat')
    for file in files:
        os.remove(file)
    #subprocess.call(['rm '+folder+'logfile* '+folder+'*.dat'], shell=True)
    subprocess.call(['cp', '../pimd',folder])
    subprocess.call(['cp', '../configuration.cfg',folder])
    cmd = folder+'pimd'
    try:
        #proc = subprocess.Popen("exec "+ex_path, cwd=f, shell=True).wait()
        proc = subprocess.Popen(cmd, cwd=folder)
        #proc = subprocess.Popen(cmd, cwd=f, shell=False, preexec_fn=os.setsid)
        proc.wait()
        #os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
    except:
        #os.kill(proc.pid, signal.SIGINT)
        #os.killpg(os.getpid(), signal.SIGTERM)
        #os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        proc.kill()
        #print(proc.pid)"""

def continue_run(folder,num_runs):
    for i in range(num_runs):
        try:
            cmd = [folder+'pimd', '-c']
            proc = subprocess.Popen(cmd, cwd=folder)
            #proc = subprocess.Popen(cmd, cwd=folder,shell=False, preexec_fn=os.setsid)
            proc.wait()
        except:
            proc.kill()
            #os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        
    #for i in range(4):
        #subprocess.Popen([ex_path,'-c'], cwd=f, shell=True).wait()
#popen = subprocess.Popen(ex_path, stdout=subprocess.PIPE)
#popen.wait()
#output = popen.stdout.read()
#print(output)

if __name__=="__main__":
    make()
    folders=['.','../run1/','../run2/','../run3/','../run4/','../run5/','../run6/']
    f = folders[3]
    first_run(f)
    continue_run(f,9)
    print(f)