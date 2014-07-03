add_executable(mi Mutual_information.cpp Mutual_information.hpp MI_main.cpp  digamma.cpp  digamma.hpp  MI_depend.cpp  MI_depend.hpp  )

find_library(LPTHREAD pthread)

target_link_libraries(mi ${LPTHREAD} m)

install(TARGETS mi DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
