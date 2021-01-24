set(SRC
    data.cpp
    callbacks/lazyCBBenders.cpp
    callbacks/lazyCBBendersMC.cpp
    callbacks/heuristicCBIncumbentUpdate.cpp
    callbacks/usercutCBfUB.cpp
    callbacks/usercutCBfUBx.cpp
    callbacks/usercutCBfUBhpp.cpp
    callbacks/usercutCBBendersMC.cpp
    callbacks/branchCB.cpp
    callbacks/nodeCB.cpp
    models/hpp.cpp
    models/master.cpp
    models/follower.cpp
    models/followerMC.cpp
    models/followerX.cpp
    models/leaderfollower.cpp
    )

set(INC_DIR 
    .
    ${DEPEND_DIR}/include/coin-or
    models
    callbacks)
set(LIBS 
	${DEPEND_DIR}/lib/libCoinUtils.a 
	${ZLIB_LIBRARIES}
	${BZIP2_LIBRARIES}
)

if(CPLEXLIB)
    set(INC_DIR 
        ${INC_DIR} 
        ${CPLEX_INC_DIR})
    set(LIBS 
        ${LIBS}
        ${CPLEX_LIB_DIR}/libilocplex.a
        ${CPLEX_LIB_DIR}/libcplex.a
        ${CPLEX_LIB_CONCERT_DIR}/libconcert.a
        -lm -lpthread -ldl)
    add_definitions(-DDSP_HAS_CPX)
endif(CPLEXLIB)

if(CMAKE_BUILD_TYPE EQUAL "DEBUG")
    add_compile_options(-Wall)
else()
    add_compile_options(-Wno-absolute-value)
endif()

if(CODE_COVERAGE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -coverage")
endif(CODE_COVERAGE)
message(STATUS "Code coverage: ${CODE_COVERAGE}")

# if(APPLE)
    set(LIBS ${LIBS})
        # ${GFORTRANLIB}
        # ${BLAS_LIBRARIES}
        # ${LAPACK_LIBRARIES}
        # -ldl)
# else(APPLE)
#     set(LIBS ${LIBS}
#         ${GFORTRANLIB}
#         ${BLAS_LIBRARIES}
#         ${LAPACK_LIBRARIES}
#         -lquadmath
#         -ldl)
# endif(APPLE)

include_directories(${INC_DIR})

# Build Object
add_library(Obj OBJECT ${SRC})
set_property(TARGET Obj PROPERTY POSITION_INDEPENDENT_CODE ON)
# add_dependencies(Obj Alps Smi Clp)

# Build libbilevel library
# add_library(bilevel SHARED $<TARGET_OBJECTS:Obj>)

# Build runbilevel executable
add_executable(runbilevel bilevel.cpp $<TARGET_OBJECTS:Obj>)

if (APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,-undefined,error")
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-undefined,error")
else(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,--no-undefined")
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
endif(APPLE)
# message(STATUS ${LIBS})
# target_link_libraries(bilevel ${LIBS})
target_link_libraries(runbilevel ${LIBS})
