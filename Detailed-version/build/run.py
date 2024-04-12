from time import time
import argparse
import os
import numpy as np
import trimesh
import os
from scipy.spatial import cKDTree
from pytorch3d.io import load_obj
from pytorch3d.structures import Meshes
from pytorch3d.ops import sample_points_from_meshes
from pytorch3d.loss import chamfer_distance
import subprocess
import psutil


def compute_shortest_distances(point_cloud_A, point_cloud_B):
    kdtree_B = cKDTree(point_cloud_B) 

    distances, _ = kdtree_B.query(point_cloud_A)  
    return distances

def isexist(name, path=None):
    '''
    '''
    if path is None:
        path = os.getcwd()
    if os.path.exists(path + '/' + name):
        print("Under the path: " + path + '\n' + name + " is exist")
        return True
    else:
        if (os.path.exists(path)):
            print("Under the path: " + path + '\n' + name + " is not exist")
        else:
            print("This path could not be found: " + path + '\n')
        return False
import numpy as np
from scipy.spatial import KDTree

def find_nearest_points_and_compute_inner_product(points1, normals1, points2, normals2):

    kdtree = KDTree(points2)
    _, indices = kdtree.query(points1)
    inner_products = np.sum(normals1 * normals2[indices], axis=1)

    return inner_products


parser = argparse.ArgumentParser()
parser.add_argument('--input','-i',help='ply format mesh')
parser.add_argument('--Number','-PN',default='20000')
parser.add_argument('--stepA','-A',default=0.16)
parser.add_argument('--stepB','-B',default=0.0)
parser.add_argument('--stepC','-C',default=0.0)
parser.add_argument('--iternum','-n',default=2)
parser.add_argument('--Ni','-Ni',default=25)
parser.add_argument('--Noise','-Noise',default=0)
args = parser.parse_args()
in_filename = args.input.replace('\\','/')
currentpath = os.getcwd().replace("\\",'/') 


Intermediate= os.path.join(currentpath,'Intermediate').replace("\\",'/')

isExists = os.path.exists(Intermediate) 

if not isExists: 
    os.mkdir(Intermediate) 

mesh = trimesh.load_mesh(in_filename)

N=int(args.Number)
points,faces=trimesh.sample.sample_surface_even(mesh,N)
normal=mesh.face_normals[faces]
mu, sigma = 0, 1 
s = np.random.normal(mu, sigma, points.shape) *args.Noise
PointCloud=np.array(points) + s
np.savetxt(Intermediate+'/pointonly'+'.xyz',  np.c_[PointCloud],fmt='%f',delimiter=' ')
np.savetxt(Intermediate+'/point'+'.xyz', np.c_[PointCloud,normal],fmt='%f',delimiter=' ')
sample=Intermediate+'/pointonly'+'.xyz'
if os.path.exists(Intermediate+'/pointonly'+'_PCAestimate.xyz'):
    os.remove(Intermediate+'/pointonly'+'_PCAestimate.xyz')   
if (N>500000000 & N <500000000)|N<8000:
    runPCA_EXE = 'PCA.exe'
    recon_cmdPCA = f"{runPCA_EXE} {sample}"
    os.system(recon_cmdPCA)
else:
    runPCA_EXE = 'PCANP.exe'
    recon_cmdPCA = f"{runPCA_EXE} {sample}"
    os.system(recon_cmdPCA)
    os.rename(Intermediate+'/pointonly'+'_NP_PCAestimate.xyz',Intermediate+'/pointonly'+'_PCAestimate.xyz')
PTS_EXE='xyz2pts.exe'
recon_cmdPTS2 = f"{PTS_EXE}  {Intermediate+'/pointonly'+'_PCAestimate.xyz'}"
os.system(recon_cmdPTS2)
    
in_filename2 =Intermediate+'/pts/'+'pointonly'+'_PCAestimate.pts'
run2_EXE = 'wave2.exe'
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
recon_cmd3 = f"{run2_EXE} pts {in_filename2} stepA {stepA} stepB {stepB} stepC {stepC} iternum {iternum} N {Ni}"
os.system(recon_cmd3)
print(recon_cmd3)

obj='BobjToObj.exe output00.bobj'
os.system(obj)

if os.path.exists(currentpath+r'\result\reconstuction.obj'):
    os.remove(currentpath+r'\result\reconstuction.obj')   
os.rename('output00.bobj.obj',currentpath+r'\result\reconstuction.obj')

if os.path.exists(currentpath+r'\result\orientation.xyz'):
    os.remove(currentpath+r'\result\orientation.xyz')  

os.rename('output2.xyz',currentpath+r'\result\orientation.xyz')

