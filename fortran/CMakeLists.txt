cmake_minimum_required(VERSION 2.8)
project(apass Fortran C)

# force a non-default underscoring method for easier linking with C code 
set(CMAKE_Fortran_FLAGS "-fno-underscoring")

add_executable(filstat filstat.f)
target_link_libraries(filstat aah)

add_executable(hepp hepp.f)

add_executable(hphot hphot.f)
#doread.f dofind.f doparm.f dophot.f doinit.f
#      dosky.f writehead.f mmm.f)

add_executable(filepp filepp.f)
target_link_libraries(filstat aah)

#
# FILTRAN
#
# Filtran requires several additional libraries for which there are not
# CMake find scripts. Rather than write find scripts, we'll just link
# to them directly from here

find_package(X11 REQUIRED)

# Output from SuperMongo compilation
link_directories("/home/kloppenb/local/lib")
set(PLOTSUB_LIBRARIES "libplotsub.a" CACHE FILEPATH "libplotsub.a")
set(DEVICES_LIBRARIES "libdevices.a" CACHE FILEPATH "libdevices.a")
set(UTILS_LIBRARIES  "libutils.a"   CACHE FILEPATH "libutils.a")

# Output from SLAlib (Positional Astronomy Library)
link_directories("/home/kloppenb/workspace/slalib/build/src")
set(SLA_LIBRARIES  "sla_f.a"   CACHE FILEPATH "libsla_f.a")

add_executable(filtran filtran.f)
target_link_libraries(filtran aah ${PLOTSUB_LIBRARIES} ${DEVICES_LIBRARIES}
  ${UTILS_LIBRARIES} ${SLA_LIBRARIES} ${X11_LIBRARIES})

add_subdirectory("aahlib")

# original Makefiles:
#gfortran -o filstat filstat.f aahlib/liba.a
#gfortran -o hepp hepp.f
#gfortran -o hphot  hphot.f
#gfortran -o hphot  hphot.f doread.f dofind.f doparm.f dophot.f doinit.f dosky.f \
#      writehead.f mmm.f

#gfortran -o filepp filepp.f /2/home/welch/hphot_and_hepp/sources/aahlib/liba.a

#gfortran -o filtran filtran.f /raid/aavsonet/sources/aahlib/liba.a \
#-lplotsub -ldevices -lutils -lsla_f -L/usr/local/lib \
#-L/usr/X11R6/lib -lX11
