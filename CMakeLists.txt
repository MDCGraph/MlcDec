cmake_minimum_required(VERSION 3.20)
project(MlcDec)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPENMP_EXE_LINKER_FLAGS}")

endif ()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DMY_DEBUG")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

add_executable(MlcDec main.cpp Graphs/Graph.cpp Graphs/Graph.h Header.h Graphs/MultilayerGraph.cpp Graphs/MultilayerGraph.h Utils/ArrayUtils.h Core/MLCore.cpp Core/MLCore.h Utils/Timer.h Core/CoreIndex.h Core/MLCTree.cpp Core/MLCTree.h Core/MLCTreeBuilder.cpp Core/MLCTreeBuilder.h Structures/DataBuf.h Structures/IntLinearHeap.h Core/KC.cpp Core/KC.h Core/MLCHashTable.cpp Core/MLCHashTable.h Structures/Trie.h Core/MLCState.h CmdParser.h Params.h Core/CoreCube.cpp Core/CoreCube.h Utils/RandVecGenerator.h Utils/MemoryUtils.h Core/MLCDfs.h Testing/WDSHelper.h Core/Cube.h DS/Metric.h Testing/MLCSearchHelper.h Testing/MLCDecHelper.h Utils/StringUtils.h Core/MLCTBuilder.h Graphs/GOrder.h Testing/GraphHelper.h Core/MLCTreeBuilder_simp.cpp Core/MLCTreeBuilder_simp.h Core/PCoreIndex.h Utils/Constructs.h Core/PMLCTreeBuilder.h Core/PMLCTreeBuilder.cpp Core/MLCGSmall.h Core/MLCore_simp.cpp Core/MLCore_simp.h Core/MLCTreeBuildOption.h DS/PWeightedDS.cpp DS/PWeightedDS.h DS/WDSPeelingState.h DS/DSMLCTree.h DS/DSMLCTree.cpp DS/DSMLCTreeBuilder.h DS/DSMLCTreeBuilder.cpp)


if (USE_JEMALLOC)
    message("USE JEMALLOC")
    add_compile_definitions("USE_JEMALLOC")
    target_link_libraries(MLC_tree libjemalloc.a dl)
else ()
    message("USE DEFAULT MEMORY ALLOCATOR")
endif ()




# FOR GPROF
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
#set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -pg")
#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

# FOR MEMTRACK
#configure_file(
#        "${PROJECT_SOURCE_DIR}/config.h.in"
#        "${PROJECT_SOURCE_DIR}/config.h"
#)
#
#option(USE_MEMTRACK "Use memtrack library" OFF)
#
#if (USE_MEMTRACK)
#    include_directories("${PROJECT_SOURCE_DIR}/memtrack")
#    add_subdirectory(memtrack)
#    Set(EXTRA_LIBS ${EXTRA_LIBS} memtrack)
#    message("memtrack on")
#else ()
#    message("memtrack off")
#endif (USE_MEMTRACK)