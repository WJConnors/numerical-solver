add_executable(finite-differences
    main.cpp
)

target_link_libraries(finite-differences PRIVATE
	numerical-common
)
