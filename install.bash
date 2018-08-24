#!/bin/bash

NIMBLE_DIR=`pwd`

echo "Enter the build directory (e.g, ../build), or press enter to create it here: "
read ans
if [ "$ans" = "" ]
then
  BUILD_DIR="build"
else
  BUILD_DIR=$ans
fi

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
pwd

LOOK_FOR_X_OUTPUT=""

#Finds the non-symbolic directory a symbolic link points back to
deref() {
  local file="$@"
  local newfile=`readlink "$file"`
  while [ "$newfile" != "" ]
  do
    file=$newfile
    newfile=`readlink "$file"`
  done
  echo $file
}

look_for_x() {
  local path=$1
  local file=$2
  if [ -d "$path" ]
  then
    if [ -x "$path/$file" ]
    then
      LOOK_FOR_X_OUTPUT="$path/$file"
      return 0
    else
      return 1
    fi
  else
    if [ -x "$path" ]
    then
      LOOK_FOR_X_OUTPUT="$path"
      return 0
    else
      return 1
    fi
  fi
}

LOOK_FOR_D_OUTPUT=""
look_for_d() {
  local path=$1
  if [ -d "$path" ]
  then  
    LOOK_FOR_D_OUTPUT=$path
    return 0
  else
    return 1
  fi
}

MPICC=`which mpicc`
while ! look_for_x $MPICC mpicc
do
  echo "mpicc not found."
  echo "Enter the path of the mpicc compiler (for C files): "
  read ans
  MPICC=$ans
done
MPICC=$LOOK_FOR_X_OUTPUT

MPICXX=`dirname $MPICC`
while ! look_for_x $MPICXX mpicxx
do
  echo "mpic++ not found."
  echo "Enter the path of the mpic++ compiler (for C++ files): "
  read ans
  MPICC=$ans
done
MPICXX=$LOOK_FOR_X_OUTPUT

echo "Using $MPICC for c files"
echo "Using $MPICXX for c++ files"

#try finding the bin folder of seacas
SEACAS_DIR=""
if which exodiff
then
  SEACAS_DIR=`which exodiff`
  SEACAS_DIR=`deref $SEACAS_DIR`
  SEACAS_DIR=`dirname "$SEACAS_DIR"`
  SEACAS_DIR=`dirname "$SEACAS_DIR"`
fi

while ! look_for_d "$SEACAS_DIR/lib" || ! look_for_d "$SEACAS_DIR/include"
do 
  echo "Couldn't find seacas installation location. "
  echo "Enter the path where the seacas/include and seacas/lib directories are located: "
  read ans
  SEACAS_DIR=$ans
done
echo "Found seacas."

print_cmake_command() {
rm -f CMakeCache.txt
echo "#!/bin/sh"
echo
echo "cmake \\"
echo "-D CMAKE_BUILD_TYPE:STRING=Release \\"
echo "-D CMAKE_C_COMPILER:STRING=\"$MPICC\" \\"
echo "-D CMAKE_CXX_COMPILER:STRING=\"$MPICXX\" \\"
echo "-D CMAKE_CXX_FLAGS:STRING=\"-O2 -std=c++11\" \\"
echo "-D EXODUS_INCLUDE_DIR:PATH=\"$SEACAS_DIR/include\" \\"
echo "-D EXODUS_LIB_DIR:PATH=\"$SEACAS_DIR/lib\" \\"
echo "-D USE_PURE_MPI:BOOL=ON \\"
echo "\"$NIMBLE_DIR\""
}

echo "Creating cmake command..."
print_cmake_command > "cmake_command.sh"
sh cmake_command.sh
make -j 4
