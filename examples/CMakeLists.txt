add_executable(matrix01
    01-matrix.cpp
)

add_executable(vector01
    01-vector.cpp
)

add_executable(so301
    01-so3.cpp
)

add_executable(se301
    01-se3.cpp
)

add_executable(plot05
    05-plot.cpp
)

add_executable(rviz05
    05-rviz.cpp
)

target_link_libraries(rviz05
    raylib
)


add_executable(perturbations10
    10-perturbations.cpp
)

add_executable(localization15
    15-localization.cpp
)
