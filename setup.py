from setuptools.command.build_ext import build_ext
from setuptools import Extension, setup, find_packages

import os, sys
sys.platform
import platform


# if sys.platform=='linux':
setup(
    	name="depth",
   	version="1.0.0",
    	author="Pavlo mozharovskyi",
    	author_email="pavlo.mozharovskyi@telecom-paris.fr",
   	description="The package provides many procedures for calculating the depth of points in an empirical distribution for many notions of data depth",
   	# long_description=long_description,
   	# long_description_content_type="text/markdown",
   	url="https://github.com/pypa",
    	packages=find_packages(),
    	ext_modules=[
        
   	 ],
    	include_package_data=True,
  	package_data={"UNIX-POSIX": ['ddalpha.so','depth_wrapper.so']},
   	zip_safe=False,
	)

# if sys.platform=='nt' and platform.architecture()[0] == "32bit":
# 	setup(
#     	name="depth",
#     	version="1.0.0",
#     	author="Pavlo mozharovskyi",
#     	author_email="pavlo.mozharovskyi@telecom-paris.fr",
#    	description="The package provides many procedures for calculating the depth of points in an empirical distribution for many notions of data depth",
#    	# long_description=long_description,
#    	# long_description_content_type="text/markdown",
#    	url="https://github.com/pypa",
#     	packages=find_packages(),
#     	ext_modules=[
        
#    	 ],
#     	include_package_data=True,
#   	  package_data={"win32": 		['ddalpha.dll','depth_wrapper.so','libatomic-1.dll','libgcc_s_seh-1.dll','libgfortran-3.dll','libgmp-10.dll','libgmpxx-4.dll','libgomp-1.dll','libgomp-plugin-host_nonshm-1.dll','libquadmath-0.dll','libssp-0.dll','libstdc++-6.dll','libvtv-0.dll','libvtv_stubs-0.dll','libwinpthread-1.dll']},
#    	 zip_safe=False,
# 	)



# if sys.platform=='nt' and platform.architecture()[0] == "64bit":
# 	setup(
#     	name="depth",
#     	version="1.0.0",
#     	author="Pavlo mozharovskyi",
#     	author_email="pavlo.mozharovskyi@telecom-paris.fr",
#    	description="The package provides many procedures for calculating the depth of points in an empirical distribution for many notions of data depth",
#    	# long_description=long_description,
#    	# long_description_content_type="text/markdown",
#    	url="https://github.com/pypa",
#     	packages=find_packages(),
#     	ext_modules=[
        
#    	 ],
#     	include_package_data=True,
#   	  package_data={"win64": 		['ddalpha.dll','depth_wrapper.so','libatomic-1.dll','libgcc_s_seh-1.dll','libgfortran-3.dll','libgmp-10.dll','libgmpxx-4.dll','libgomp-1.dll','libgomp-plugin-host_nonshm-1.dll','libquadmath-0.dll','libssp-0.dll','libstdc++-6.dll','libvtv-0.dll','libvtv_stubs-0.dll','libwinpthread-1.dll']},
#    	 zip_safe=False,
# 	)



	
