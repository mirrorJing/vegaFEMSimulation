#################################################################
# Scons script for vegaFEMSimulation project
# @author: Fei Zhu, 01/25/2014
# Usage: enter root directory of the project in terminal and
#        enter "scons"                
#################################################################

######
import fnmatch
import os
from os.path import basename
import platform
from glob import glob


#OS TYPE
os_name=platform.system()
os_architecture=platform.architecture()[0]

#BUILD TYPE
build_type='release'
#build_type='debug'

#BUILD MSVC PROJECTS FOR WINDOWS
build_msvc=True
#build_msvc=False

#GL LIBRARY
gl_include_path='GL/'
gl_library_path='GL/lib/'

#VegaFEM LIBRARY	
vega_include_path='VegaFEM/include/'
vega_library_path='VegaFEM/lib/'

#Pthread LIBRARY
pthread_include_path='pthread/include/'
pthread_library_path='pthread/lib/'

#NLopt LIBRARY
#nlopt_include_path='nlopt/include/'
#nlopt_library_path='nlopt/lib/'

#Opt++ LIBRARY
#optpp_include_path='optpp/include/'
#optpp_library_path='optpp/lib/'

#ALGLIB LIBRARY
#alglib_include_path='alglib/'

#PATHS MADE PLATFORM SPECIFIC
if os_name=='Linux':
    vega_library_path=vega_library_path+'Linux/'
    #nlopt_library_path=nlopt_library_path+'Linux/'
    #optpp_library_path=optpp_library_path+'Linux/'
    gl_library_path=gl_library_path+'Linux/'
    pthread_library_path=pthread_library_path+'Linux/'
elif os_name=='Darwin':
    vega_library_path=vega_library_path+'Apple/'
    #nlopt_library_path=nlopt_library_path+'Apple/'
    #optpp_library_path=optpp_library_path+'Apple/'
    gl_library_path=gl_library_path+'Apple/'
elif os_name=='Windows':
    vega_library_path=vega_library_path+'Windows/'
    #nlopt_library_path=nlopt_library_path+'Windows/'
    #optpp_library_path=optpp_library_path+'Windows/'
    gl_library_path=gl_library_path+'Windows/'
    pthread_library_path=pthread_library_path+'Windows/'

if os_architecture=='32bit':
    vega_library_path=vega_library_path+'X86/'
    #nlopt_library_path=nlopt_library_path+'X86/'
    #optpp_library_path=optpp_library_path+'X86/'
else:
    vega_library_path=vega_library_path+'X64/'
    #nlopt_library_path=nlopt_library_path+'X64/'
    #optpp_library_path=optpp_library_path+'X64/'

#vegaFEMSimulation
source_filename=Glob('*.cpp',True,False,True)
#source_filename.append(Glob(alglib_include_path+'*.cpp',True,False,True))

target_filename='vegaFEMSimulation'+build_type
if build_type=='release':
  	  target_filename='vegaFEMSimulation'
	#target_filename='secondVegaFEMSimulation'

#LIB FILES
lib_files=[]
lib_files.append(Split('sceneObject integrator elasticForceModel forceModel loadList insertRows lighting performanceCounter configFile volumetricMesh getopts camera graph isotropicHyperelasticFEM stvk corotationalLinearFEM polarDecomposition massSpringSystem objMesh sparseSolver matrix sparseMatrix minivector'))

#lib_files.append('nlopt')
#lib_files.append('opt')
#lib_files.append('newmat')
if os_architecture=='32bit':
    lib_files.append('glui32')
else:
    lib_files.append('glui64')

if os_name=='Linux':
   lib_files.append('glut')
   lib_files.append('GLU')
   lib_files.append('GL')	
   lib_files.append('pthread')
elif os_name=='Darwin':
   lib_files.append('pthread')
elif os_name=='Windows':
   lib_files.append('glut32')
   lib_files.append('glu32')
   lib_files.append('opengl32')
   lib_files.append('pthreadVC2')
   lib_files.append('pthreadVCE2')
   lib_files.append('pthreadVSE2')

#COMPILER
compiler=''
if os_name in ('Linux','Darwin') or (os_name=='Windows' and build_msvc==False):
   compiler='g++'
else:
   compiler='msvc'

#SRC FILES
src_root_path='./'
ignored_src_path=[]
src_files=[]
inc_files=[]
for dir,_,_ in os.walk(src_root_path):
    if dir not in ignored_src_path:
       src_files.extend(glob(os.path.join(dir,'*.cpp')))
       inc_files.extend(glob(os.path.join(dir,'*.h')))

#COMPILER OPTIONS
CC='g++'
CXX='g++'
tools=['gcc', 'g++', 'gnulink']
CPPPATH=[gl_include_path,vega_include_path]
LIBPATH=[gl_library_path,vega_library_path,pthread_library_path]
RPATH=[gl_library_path,pthread_library_path]
LIBS=lib_files
ENV={'PATH':os.environ['PATH']}
if build_type=='release':
    CCFLAGS=['-O3','-fno-strict-aliasing','-std=gnu++0x','-DNDEBUG','-DHAVE_STD','-DHAVE_NAMESPACES']
else:
    CCFLAGS=['-std=gnu++0x','-fno-strict-aliasing','-g','-DHAVE_STD','-DHAVE_NAMESPACES']

env=Environment(CC=CC,CXX=CXX,tools=tools,CCFLAGS=CCFLAGS,CPPPATH=CPPPATH,LIBPATH=LIBPATH,RPATH=RPATH,LIBS=LIBS,ENV=ENV)
#MAC OS: USE OPENGL GLUT FRAMEWORKS SHIPPED WITH THE SYSTEM
if os_name=='Darwin':
    env.Append(LINKFLAGS=['-framework','OpenGL'])
    env.Append(LINKFLAGS=['-framework','GLUT'])

#BUILD
target=env.Program(target_filename,source_filename)

#WINDOWS WORKAROUND: COPY DLLS TO EXECUTIVE DIRECTORY
if os_name=='Windows':
    for rpath in RPATH:
        for dll_name in os.listdir(rpath):
            if dll_name.endswith('.dll'):
                Command(dll_name, rpath+dll_name, Copy("$TARGET", "$SOURCE"))



