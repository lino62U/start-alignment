cmake_minimum_required(VERSION 3.14)
project(AlineadorSecuencias CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

# Incluir directorios
include_directories(include)

# Buscar todos los archivos fuente excepto main
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.cpp)
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

# Crear la biblioteca con toda la lógica
add_library(alineador_lib ${SOURCES})

# Ejecutable principal
add_executable(alineador src/main.cpp)
target_link_libraries(alineador alineador_lib)

# -------------------------------
# GoogleTest (opcional)
# -------------------------------
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/heads/main.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)

# Ejecutables de prueba (si hay)
file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS tests/*.cpp)

if(TEST_SOURCES)
  add_executable(test_alineador ${TEST_SOURCES})
  target_link_libraries(test_alineador gtest_main alineador_lib)
  include(GoogleTest)
  gtest_discover_tests(test_alineador)
endif()
