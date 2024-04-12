import argparse
import os
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--input','-i',help='xyz format point cloud input')
parser.add_argument('--stepA','-A',default=0.12)
parser.add_argument('--stepB','-B',default=0.0)
parser.add_argument('--stepC','-C',default=0.0)
parser.add_argument('--iternum','-n',default=3)
parser.add_argument('--Ni','-Ni',default=25)
parser.add_argument('--Noise','-Noise',default=0)
args = parser.parse_args()
in_filename = args.input.replace('\\','/')
path=os.getcwd()
os.chdir(path+'/build')
Point=np.loadtxt(in_filename)
N=Point.shape[0]
if os.path.exists('Intermediate'+'/pointonly'+'_PCAestimate.xyz'):
    os.remove('Intermediate'+'/pointonly'+'_PCAestimate.xyz')  
np.savetxt('Intermediate'+'/pointonly.xyz',Point,fmt='%f',delimiter=' ')
if (N>500000000 & N <500000000)|N<2000:
    runPCA_EXE = 'PCA.exe'
    recon_cmdPCA = f"{runPCA_EXE} {'Intermediate'+'/pointonly.xyz'}"
    os.system(recon_cmdPCA)
else:
    runPCA_EXE = 'PCANP.exe'
    recon_cmdPCA = f"{runPCA_EXE} {'Intermediate'+'/pointonly.xyz'}"
    os.system(recon_cmdPCA)
    os.rename('Intermediate'+'/pointonly'+'_NP_PCAestimate.xyz','Intermediate'+'/pointonly'+'_PCAestimate.xyz')
PTS_EXE='xyz2pts.exe'
if os.path.exists('Intermediate'+'/pts/pointonly_PCAestimate.pts'):
    os.remove('Intermediate'+'/pts/pointonly_PCAestimate.pts')  
recon_cmdPTS = f"{PTS_EXE}  {'Intermediate'+'/pointonly_PCAestimate.xyz'}"
os.system(recon_cmdPTS)

in_filename='Intermediate'+'/pts/pointonly'+'_PCAestimate.pts'

run2_EXE = path+r'\build\wave2.exe'
stepA=args.stepA
if(args.stepB==0):
    stepB=args.stepA
else:
    stepB=args.stepB
if(args.stepC==0):
    stepC=args.stepA
else:
    stepC=args.stepC

iternum=args.iternum
Ni=args.Ni
recon_cmd3 = f"{run2_EXE} pts {in_filename}  stepA {stepA} stepB {stepB} stepC {stepC} iternum {iternum} N {Ni}"
os.system(recon_cmd3)
obj='BobjToObj.exe output00.bobj'
os.system(obj)

if os.path.exists(path+r'\build\result\reconstuction.obj'):
    os.remove(path+r'\build\result\reconstuction.obj')   
os.rename(path+r'\build\output00.bobj.obj',path+r'\build\result\reconstuction.obj')

if os.path.exists(path+r'\build\result\orientation.xyz'):
    os.remove(path+r'\build\result\orientation.xyz')  

os.rename(path+r'\build\output2.xyz',path+r'\build\result\orientation.xyz')



