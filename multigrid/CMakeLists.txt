add_executable(multigrid
	main.cpp
)

target_link_libraries(poisson PRIVATE
	numerical-common
)
