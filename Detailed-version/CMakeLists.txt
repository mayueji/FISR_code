PROJECT(wave2) # project name
cmake_minimum_required(VERSION 3.1)

ADD_EXECUTABLE(
   wave2
average.cpp  coeffs_daub4.cpp  connectneighbors.cpp  eval_haar.cpp  main.cpp           matrix.cpp   octtree.cpp  plane.cpp  readpts.cpp      streamingvisitor.cpp
blur.cpp     coeffs_haar.cpp   eval_daub4.cpp        global.cpp     marchingcubes.cpp  octnode.cpp  parse.cpp    prune.cpp  reconstruct.cpp  traverse.cpp normal.cpp
iterator.cpp iterator_daub4.cpp)
find_package(CGAL REQUIRED)


TARGET_LINK_LIBRARIES(
   wave2                           # to wave link
   m
   CGAL
)

###################

