cmake_minimum_required(VERSION 2.8)
project(aahlib Fortran)

file(GLOB SOURCE "*.f")
# exclude the "solve3.f" file
list(REMOVE_ITEM SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/solve3.f)

add_library(aah ${SOURCE})

# Original Makefile was as follows:
#g77  -c -fno-automatic bell.f
#g77  -c -fno-automatic disply.f
#g77  -c -fno-automatic extrd.f
#g77  -c -fno-automatic getfilt.f
#g77  -c -fno-automatic helcor.f
#g77  -c -fno-automatic jday.f
#g77  -c -fno-automatic lin.f
#g77  -c -fno-automatic nblank.f
#g77  -c -fno-automatic openit.f
#g77  -c -fno-automatic prompt.f
#g77  -c -fno-automatic rdnames.f
#g77  -c -fno-automatic solve2.f
#g77  -c -fno-automatic ucase.f
#g77  -c -fno-automatic wrtnames.f
