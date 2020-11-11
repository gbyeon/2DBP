##################################################
# Build and install external packages
##################################################
include(ExternalProject)

set(EXTRA_DIR      "${CMAKE_SOURCE_DIR}/extra")
set(COINUTILS_DIR  "${EXTRA_DIR}/CoinUtils")
set(COIN_INC_DIR   "${DEPEND_DIR}/include/coin")
set(COINUTILS_LIB  "-L${DEPEND_DIR}/lib -lCoinUtils")

set(EXTRA_CXX_FLAGS "-fPIC -Wno-sign-conversion -Wno-float-conversion -Wno-shorten-64-to-32 -Wno-unused-variable -Wno-constant-logical-operand -Wno-dynamic-class-memaccess -Wno-unused-private-field")

##################################################
# Compile CoinUtils
##################################################

ExternalProject_Add(CoinUtils
	SOURCE_DIR ${COINUTILS_DIR}
	CONFIGURE_COMMAND ${COINUTILS_DIR}/configure 
		--prefix=${DEPEND_DIR}
		--enable-static 
		--disable-shared
		ADD_CXXFLAGS=${EXTRA_CXX_FLAGS}
	BUILD_COMMAND make -j1
	INSTALL_COMMAND make -j1 install)
