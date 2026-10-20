add_executable(multigrid
	main.cpp
)

target_link_libraries(multigrid PRIVATE
	numerical-common
)
