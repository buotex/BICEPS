  project(DIRECTAG)
  include_directories(${Boost_INCLUDE_DIRS}) 
  add_subdirectory(ext/src/expat-2.0.1)
  add_subdirectory(pwiz-src/pwiz)
  #add_subdirectory(pwiz-src/libraries)
  add_subdirectory(freicore)
  add_subdirectory(directag)  





