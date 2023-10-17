#!/bin/sh
#
# installation of FRET_Embedlab
#
#
README=README.md
ompvar=0
while test $# -gt 0; do
        case "$1" in
                -h|--help)
                        head -n46 $README 
                        exit 0
                        ;;
                -b|--build)     #build directory 
                        shift
                        export build=$1
                        shift
                        ;;
                -fc|--fortrancompiler) # compilatore
                        shift
                        export FORTRANCOMPILER=$1
                        shift
                        ;;
                -omp|--omp)   #open mp 
                        shift
                        ompvar=1
                        ;;
                *)
                        break
                        ;;
        esac
                
done

if [ -z ${build+x} ]; then 
    buildir="build"
else 
    buildir=$build
fi

if [ ! -d $buildir ] ; then
   mkdir $buildir
   cd $buildir
else
   echo "$buildir folder already exists"
   echo "change its namewith -b option"
   echo "Stop"
   exit 0
fi

if [ -z ${FORTRANCOMPILER+x} ]; then 
    export FC=$(which gfortran)
else 
    export FC=$FORTRANCOMPILER
fi

if [ $ompvar -eq 0 ] ; then 
   cmake ..
else
   cmake .. -DENABLE_OMP=ON
fi
cd ..

echo ""
echo " To complete the compilation:"
echo "$ cd $buildir"
echo "$ make"
echo ""
echo " To test the compilation:"
echo "$ cd $buildir"
echo "$ ctest"
echo ""
