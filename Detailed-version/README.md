Windows & Linux: The code is tested by Visual Studio Code
The provision of detailed article source code is beneficial for: 
(1) Further research can extract the implementation steps of the core idea of the article, which helps readers understand the article
(2) For ease of execution, some adjustable hyperparameters are not used as input parameters during execution, but can be adjusted in the source code (main.cpp, iterator_daub4.cpp), especially in terms of the depth of octree and parameter update methods. After adjustment, the .exe  can be run again to obtain better results.

[Remind]:
Build \ CMakeCache. txt must be deleted before compilation

[Dataset]:
[mesh]:\data\mesh
[point clouds(.xyz)]:\data\sparse


[Path]:
--cd "current directory"
--cd build

[Build_.exe_file]:

-If you want to run in the Linux Platform:

--cmake
--make

-If you want to run in the Windows Platform: (Cmd/PowerShell)

--cmake -G "MinGW Makefiles" ..
--mingw32-make

Run the generated program:
(1) Input file type is .xyz/.txt(without normal) :[suggest]

(Through the runiWSR.py file)

-- cd..
(Return to：Detailed-version)

[Standard_input_format]: 
-- python runiWSR.py '-i' '.xyz/.txt format point cloud input' '-A' 'step1 of the method(abnormal value point)' '-B' 'step2 of the method(abnormal value point)'(Default:equal to A)
'-C' 'step3 of the method(abnormal gradient point)'(Default:equal to A) '-n' 'number of iterations'(Default:3) '-Ni' 'The number of update rounds included in each iteration'(
Default:25) '-Noise' 'The level of Gaussian noise'

Example data: \data\sparse(4-10K)
Recommended parameters(-A -n -Ni)[may not be the best parameters]
[4-10K]
-A 0.06 -n 2 -Ni 25 e.g. for 
-A 0.08 -n 3 -Ni 25 e.g. for 
-A 0.08 -n 2 -Ni 25 e.g. for 
-A 0.16 -n 3 -Ni 25 e.g. for kitten,BS_4000_tonus, BS_4000_holes,hammer,heels,tennyson,johnwesley

[>=40K]
-A 0.16 -n 3 -Ni 25 ; -A 0.24 -n 4 -Ni 20 ;  -A 0.32 -n 5 -Ni 20 

[Example] python runiWSR.py -i 'D:/Surface_Reconstruction/Detailed-versiondata/sparse/BS_4000_holes.xyz'(need replace) -A 0.06 -n 2 -Ni 25 

[Algorithm output in this mode(A)](reconstruction) reconstrcution.obj {in build/result}
[Algorithm output in this mode(B)](oriented_point_cloud) orientation.xyz {in build/result}


(2) Input file type is mesh(.ply)[Need to sample and obtain point clouds on the mesh]:[suggest]
[method1]:Use the following code for sampling, and then process the saved. xyz using (1)
<!-- points,faces=trimesh.sample.sample_surface_even(mesh,N)
normal=mesh.face_normals[faces]
PointCloud=np.array(points) -->

[method2]:Utilize the Python files that we have already integrated.
--cd build

-- python runiWSR.py '-i' '.ply format mesh' '-PN' 'The number of points want to be sampled' '-A' 'step1 of the method(abnormal value point)' '-B' 'step2 of the method(abnormal value point)'(Default:equal to A) '-C' 'step3 of the method(abnormal gradient point)'(Default:equal to A) '-n' 'number of iterations'(Default:3) '-Ni' 'The number of update rounds included in each iteration'(Default:25) '-Noise' 'The level of Gaussian noise'

[Example] python run.py -i run.py -i "D:/Surface_Reconstruction/dataset/datasets/famous_dense/03_meshes/column.ply"(need replace) -A 0.06 -n 2 -Ni 25 


(3) Input file type is .pts (Standard input format for wavelet surface reconstruction(with inconsistent oriented normal)):

--cd build
(Enter build)

--wave2.exe pts {in_filename}  stepA {stepA} stepB {stepB} stepC {stepC} iternum {iternum} N {Ni}
(runiwsr)

[Example] wave2.exe pts D:\Surface_Reconstruction\WaveletPipeRecon_Linux64\prepare\Detailed-version\build\Intermediate\pts\pointonly_PCAestimate.pts stepA 0.16 stepB 0.16 stepC 0.16 iternum 3 N 20

--BobjToObj.exe output00.bobj
(Convert output format)
[Algorithm output in this mode(A)](reconstruction) output00.bobj.obj {in build}
[Algorithm output in this mode(B)](oriented_point_cloud) output2.xyz {in build}

(4) Input file type is .xyz(with inconsistent oriented normal):

convert the pre packaged xyz2pts.exe to pts, and follow the same steps as (3)

--cd build
--xyz2pts.exe example_PCA.xyz(need replace)

[Algorithm output in this mode(A)](reconstruction) output00.bobj.obj {in build}
[Algorithm output in this mode(B)](oriented_point_cloud) output2.xyz {in build}

